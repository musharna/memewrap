"""FIMO: scan sequences for occurrences of known motifs.

## Prior art, stated honestly

Unlike STREME and TOMTOM, FIMO *is* already available to Python:
`pymemesuite` (0.1.0a4) provides real Cython bindings, no subprocess and no
temporary files. If you can install it, prefer it.

This module exists for the case that package does not cover: driving the
installed `fimo` CLI, which is what a conda/bioconda MEME environment gives you
without a compiler, and doing it across processes because `fimo --text` on a
whole plant promoter set against a full JASPAR database is minutes of CPU that
parallelises perfectly over sequence chunks.

## The header bug in the source, fixed here

The source merged chunk outputs with

    for i, o in enumerate(outs):
        for j, line in enumerate(fh):
            if j == 0 and i > 0:   # keep only the first chunk's header
                continue
            out.write(line)

which assumes chunk 0 produced a header. `fimo --text` writes a header only
when it writes output at all, so an EMPTY chunk 0 -- routine, since chunks are
dealt round-robin and a scan can legitimately find nothing in one -- writes zero
lines. The merge then reaches chunk 1, strips its line 0 as a "duplicate
header", and the merged TSV begins with a data row.

Downstream that file is read with `usecols=["motif_id", "sequence_name",
"score"]`, so pandas takes the first DATA row as the column names and either
raises a confusing KeyError or, worse, silently loses one hit and mislabels
every column. Nothing in the pipeline says "the header was dropped".

The fix is to keep the header from the first chunk that HAS one, and to skip a
leading header line in every chunk regardless of index.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from .tools import find_tool

__all__ = ["FimoError", "run_fimo", "run_fimo_parallel", "build_feature_matrix"]


class FimoError(RuntimeError):
    """FIMO exited non-zero."""


def _read_fasta_records(path: Path) -> list[tuple[str, str]]:
    """(header, sequence) pairs. Deliberately dependency-free.

    Biopython is a heavy dependency to add for splitting a FASTA, and the source
    only used SeqIO.parse/write to round-trip the records unchanged.
    """
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header, seq = line[1:], []
            elif line:
                seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def split_fasta(
    fasta: str | os.PathLike, outdir: str | os.PathLike, n_chunks: int
) -> list[Path]:
    """Deal records round-robin into at most `n_chunks` files.

    Never returns an empty chunk. The source computed `recs[i::n_chunks]` for
    every i in `range(n_chunks)`, so asking for more chunks than there are
    sequences produced empty FASTAs, and `fimo` on an empty input is a
    subprocess launched to do nothing -- or, depending on version, an error that
    reads as a real scan failure.

    The `min(n_chunks, len(records))` bound is what prevents that, and it is
    sufficient on its own: for every i below it, `records[i::n_chunks]` contains
    at least `records[i]`. An additional `if not part: continue` guard was here
    and is gone -- mutation testing showed it unreachable, and an unreachable
    guard reads as protection that is actually being provided by something else.
    """
    if n_chunks < 1:
        raise ValueError(f"n_chunks must be >= 1, got {n_chunks}")
    records = _read_fasta_records(Path(fasta))
    if not records:
        raise ValueError(f"no FASTA records in {fasta}")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for i in range(min(n_chunks, len(records))):
        part = records[i::n_chunks]
        p = outdir / f"chunk_{i}.fa"
        with open(p, "w") as fh:
            for h, s in part:
                fh.write(f">{h}\n{s}\n")
        paths.append(p)
    return paths


def run_fimo(
    fasta: str | os.PathLike,
    motifs: str | os.PathLike,
    out_tsv: str | os.PathLike,
    *,
    thresh: str | float = "1e-4",
    meme_bin: str | os.PathLike | None = None,
) -> Path:
    """Scan one FASTA with `fimo --text`, writing the TSV to `out_tsv`."""
    fasta, motifs, out_tsv = Path(fasta), Path(motifs), Path(out_tsv)
    if not fasta.is_file():
        raise FileNotFoundError(f"FASTA does not exist: {fasta}")
    if not motifs.is_file():
        raise FileNotFoundError(f"motif file does not exist: {motifs}")
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(find_tool("fimo", meme_bin)),
        "--text",
        "--thresh",
        str(thresh),
        str(motifs),
        str(fasta),
    ]
    with open(out_tsv, "w") as fh:
        proc = subprocess.run(
            cmd, stdout=fh, stderr=subprocess.PIPE, text=True, check=False
        )
    if proc.returncode != 0:
        raise FimoError(
            f"fimo exited {proc.returncode}\n"
            f"  command: {' '.join(cmd)}\n"
            f"  stderr: {proc.stderr.strip()[-2000:]}"
        )
    return out_tsv


def merge_fimo_tsvs(parts: list[Path], merged: str | os.PathLike) -> Path:
    """Concatenate FIMO TSVs, keeping exactly one header.

    See the module docstring: the header is taken from the first part that HAS
    one, and a leading header line is stripped from every part -- not from
    "every part except the first", which is what breaks when part 0 is empty.
    """
    merged = Path(merged)
    header: str | None = None
    body: list[str] = []
    for p in parts:
        with open(p) as fh:
            lines = fh.readlines()
        if not lines:
            continue  # an empty chunk contributes nothing, INCLUDING no header
        first = lines[0]
        # FIMO's --text header begins with the motif_id column name.
        if first.startswith("motif_id") or first.startswith("#"):
            if header is None:
                header = first
            body.extend(lines[1:])
        else:
            body.extend(lines)
    with open(merged, "w") as out:
        if header is not None:
            out.write(header)
        out.writelines(body)
    return merged


def run_fimo_parallel(
    fasta: str | os.PathLike,
    motifs: str | os.PathLike,
    outdir: str | os.PathLike,
    *,
    n_chunks: int = 4,
    thresh: str | float = "1e-4",
    meme_bin: str | os.PathLike | None = None,
) -> Path:
    """Scan `fasta` with FIMO across `n_chunks` processes; return the merged TSV."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    chunks = split_fasta(fasta, outdir, n_chunks)

    fimo_exe = str(find_tool("fimo", meme_bin))
    procs = []
    for i, ch in enumerate(chunks):
        o = outdir / f"fimo_{i}.txt"

        # launched concurrently and each child writes to its own stdout for as
        # long as it runs, so the file must stay open past the end of this
        # iteration. It is closed in the collection loop below, after wait().
        fh = open(o, "w")  # noqa: SIM115
        p = subprocess.Popen(
            [fimo_exe, "--text", "--thresh", str(thresh), str(motifs), str(ch)],
            stdout=fh,
            stderr=subprocess.PIPE,
            text=True,
        )
        procs.append((p, fh, o))

    failures = []
    outs: list[Path] = []
    for p, fh, o in procs:
        _, err = p.communicate()
        fh.close()
        if p.returncode != 0:
            # Collect every failure before raising. The source raised on the
            # first, leaving the remaining children running and unwaited.
            failures.append(
                f"{o.name} (rc={p.returncode}): {(err or '').strip()[-400:]}"
            )
        else:
            outs.append(o)
    if failures:
        raise FimoError("FIMO chunk(s) failed:\n  " + "\n  ".join(failures))

    return merge_fimo_tsvs(outs, outdir / "fimo.tsv")


def build_feature_matrix(
    fimo_tsv: str | os.PathLike, gene_ids: list, motif_ids: list
) -> np.ndarray:
    """(n_genes, n_motifs) matrix of the MAX FIMO score per pair; 0 where no hit.

    Rows and columns follow `gene_ids` / `motif_ids` exactly, so the caller's
    label vectors stay aligned. A gene or motif absent from the scan is an
    all-zero row/column rather than a missing one -- dropping it would silently
    shift every downstream index.

    Note that 0 means "no hit at the scan threshold", not "score zero". FIMO
    scores can legitimately be negative, so a real weak hit is DISTINGUISHABLE
    from an absent one only if you scanned at a permissive threshold. This is
    the source's convention and it is fine for the sparse-regression use it was
    written for; it is worth knowing before feeding the matrix to something that
    treats 0 as a magnitude.
    """
    gene_ids, motif_ids = list(gene_ids), list(motif_ids)
    idx = {g: i for i, g in enumerate(gene_ids)}
    mdx = {m: j for j, m in enumerate(motif_ids)}
    X = np.zeros((len(gene_ids), len(motif_ids)))

    df = pd.read_csv(fimo_tsv, sep="\t", usecols=["motif_id", "sequence_name", "score"])
    df = df[df["sequence_name"].isin(idx) & df["motif_id"].isin(mdx)]
    if df.empty:
        return X
    grp = df.groupby(["sequence_name", "motif_id"], sort=False)["score"].max()
    for (g, m), s in grp.items():
        X[idx[g], mdx[m]] = s
    return X
