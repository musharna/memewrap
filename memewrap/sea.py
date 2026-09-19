"""SEA: enrichment of KNOWN motifs in primary sequences relative to controls.

STREME asks "what motif is enriched here?"; SEA asks "which of THESE motifs are
enriched here?" (Bailey & Grant 2021, bioRxiv, doi:10.1101/2021.08.23.457422).
It is the usual next step after a TOMTOM match, and the successor to AME for the
two-set question -- see `memewrap.ame` for when AME is still the one to use.

## A motif missing from the table was TESTED, and lost

SEA reports only motifs passing `--thresh`, and its default is E-value <= 10.
So the row count of `sea.tsv` is the number of motifs that passed, not the
number scored; a motif absent from the table was tested and not enriched, and a
header-only table is a legitimate result, not a failed run. Measured against
SEA 5.5.9 on a two-motif database with one motif planted: both rows come back
at the default threshold (the decoy at p=0.99, E=1.98), and only the planted one
at `thresh=0.05`. `read_sea` therefore returns an empty frame WITH typed
columns for a header-only file, so `df[df.QVALUE < 0.05]` works either way.

`thresh` is an E-value unless `threshold_on` says otherwise. SEA takes that as
two separate flags (`--qvalue`, `--pvalue`) which it lets you pass together;
one keyword with three legal values cannot express the contradiction.

## Why `--text` is the default here when it is not for FIMO

`run_fimo` defaults to `--oc` because `fimo --text` skips the q-value
computation. SEA does not: `sea --text` writes the same table, QVALUE column
populated, to stdout (verified against 5.4.1 and 5.5.9). What `--text` drops is
`sites.tsv`, `sequences.tsv` and the HTML report.

And `--oc` is what breaks. The bioconda build of MEME 5.5.9 this was developed
against (`pl5321he99cc7f_1`, reproduced on a fresh CI install) ships a `sea_template.html` with no data section;
`sea --oc` computes everything, fails rendering the HTML, exits 1 with

    FATAL: Template does not contain data section.

and leaves a header-only `sea.tsv` behind -- which, per the section above, is
indistinguishable from "nothing was enriched" to anything that reads the file
without checking the exit status. The same build of `ame` fails the same way.
`text=False` is available for installs where the template is intact, and the
non-zero exit is raised as `SeaError` with that stderr, never swallowed.

## Underflow

PVALUE/EVALUE/QVALUE are printed in decimal and parse to 0.0 below ~1e-308.
SEA supplies LOG_PVALUE/LOG_EVALUE/LOG_QVALUE (natural log) for that reason;
rank on those, not on the underflowed value.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pandas as pd

from .tools import find_tool

__all__ = ["SeaError", "SEA_DTYPES", "run_sea", "read_sea"]


class SeaError(RuntimeError):
    """SEA exited non-zero, or exited 0 without a table. Carries its stderr."""


# Column names are SEA's own, as TOMTOM's and FIMO's are kept elsewhere in this
# package. Counts are integers; everything statistical is float64.
SEA_DTYPES: dict[str, str] = {
    "RANK": "int64",
    "DB": "string",
    "ID": "string",
    "ALT_ID": "string",
    "CONSENSUS": "string",
    "TP": "int64",
    "TP%": "float64",
    "FP": "int64",
    "FP%": "float64",
    "ENR_RATIO": "float64",
    "SCORE_THR": "float64",
    "PVALUE": "float64",
    "LOG_PVALUE": "float64",
    "EVALUE": "float64",
    "LOG_EVALUE": "float64",
    "QVALUE": "float64",
    "LOG_QVALUE": "float64",
}

_THRESHOLD_FLAGS = {"evalue": None, "qvalue": "--qvalue", "pvalue": "--pvalue"}


def _motif_files(motifs: str | os.PathLike | list) -> list[Path]:
    files = [
        Path(m)
        for m in ([motifs] if isinstance(motifs, (str, os.PathLike)) else motifs)
    ]
    if not files:
        raise ValueError("no motif files given")
    for m in files:
        if not m.is_file():
            raise FileNotFoundError(f"motif file does not exist: {m}")
    return files


def run_sea(
    primary: str | os.PathLike,
    motifs: str | os.PathLike | list,
    outdir: str | os.PathLike,
    *,
    control: str | os.PathLike | None = None,
    thresh: float | None = None,
    threshold_on: str = "evalue",
    order: int | None = None,
    seed: int | None = None,
    text: bool = True,
    extra: list[str] | None = None,
    meme_bin: str | os.PathLike | None = None,
) -> Path:
    """Run SEA; return the path to `sea.tsv` in `outdir`. Parse it with `read_sea`.

    Args:
        primary: FASTA of sequences to test for enrichment IN.
        motifs: one MEME-format motif file, or a list of them.
        outdir: output directory, created if absent.
        control: FASTA of background sequences. None lets SEA build controls by
            shuffling `primary` (preserving `order`-mers) -- pass `seed` to make
            that reproducible across SEA versions' defaults.
        thresh: reporting threshold; None is SEA's default (E-value <= 10).
            Motifs failing it are OMITTED from the table -- see module docstring.
        threshold_on: what `thresh` applies to: "evalue", "qvalue" or "pvalue".
        order: Markov order of the background model and of the shuffle.
        seed: random seed for the shuffle.
        text: True (default) runs `sea --text` and writes its stdout to
            `outdir/sea.tsv`. False runs `sea --oc outdir`, which also writes
            sites.tsv, sequences.tsv and sea.html -- and fails on MEME builds
            with a broken HTML template; see the module docstring.
        extra: additional raw flags appended verbatim.

    Raises:
        ValueError: unknown `threshold_on`, `thresh` out of range, or no motif
            files.
        FileNotFoundError: an input file is missing; the message says which.
        MemeToolNotFound: no `sea` executable.
        SeaError: SEA exited non-zero, or produced no table.
    """
    if threshold_on not in _THRESHOLD_FLAGS:
        raise ValueError(
            f"threshold_on must be one of {sorted(_THRESHOLD_FLAGS)}, "
            f"got {threshold_on!r}"
        )
    if thresh is not None and thresh <= 0:
        raise ValueError(f"thresh must be positive, got {thresh}")
    if thresh is not None and threshold_on != "evalue" and thresh > 1:
        # SEA's own rule. It reports it on stderr and exits non-zero, but only
        # after loading every motif file.
        raise ValueError(
            f"thresh must be in (0, 1] when threshold_on={threshold_on!r}, "
            f"got {thresh}; only an E-value threshold may exceed 1"
        )
    primary = Path(primary)
    if not primary.is_file():
        raise FileNotFoundError(f"primary FASTA does not exist: {primary}")
    if control is not None and not Path(control).is_file():
        raise FileNotFoundError(f"control FASTA does not exist: {control}")
    motif_files = _motif_files(motifs)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tsv = outdir / "sea.tsv"
    sea_exe = str(find_tool("sea", meme_bin))
    # A table from an earlier run in this directory must not outlive a failed
    # one: the caller would read old results as this run's.
    tsv.unlink(missing_ok=True)

    cmd = [sea_exe]
    cmd += ["--text"] if text else ["--oc", str(outdir)]
    cmd += ["--p", str(primary)]
    if control is not None:
        cmd += ["--n", str(control)]
    for m in motif_files:
        cmd += ["--m", str(m)]
    if thresh is not None:
        cmd += ["--thresh", str(thresh)]
    flag = _THRESHOLD_FLAGS[threshold_on]
    if flag is not None:
        cmd.append(flag)
    if order is not None:
        cmd += ["--order", str(order)]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    if extra:
        cmd += list(extra)

    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise SeaError(
            f"sea exited {proc.returncode}\n"
            f"  command: {' '.join(cmd)}\n"
            f"  stderr: {proc.stderr.strip()[-2000:]}"
        )
    if text:
        # Written only after a zero exit, so a failed run cannot leave behind a
        # table that reads as "nothing enriched".
        tsv.write_text(proc.stdout)
    if not tsv.is_file() or not tsv.read_text().startswith("RANK\t"):
        raise SeaError(
            f"sea exited 0 but produced no results table at {tsv}\n"
            f"  command: {' '.join(cmd)}\n"
            f"  stderr: {proc.stderr.strip()[-2000:]}"
        )
    return tsv


def read_enrichment_tsv(
    tsv: str | os.PathLike,
    dtypes: dict[str, str],
    tool: str,
    optional: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Read a SEA/AME results table with declared dtypes; shared with `ame`.

    The trailing provenance block (`# SEA (Simple Enrichment Analysis): ...`)
    and the blank line before it are dropped. A column this package declares
    and the file lacks is an error naming it, not a KeyError three calls later:
    it means a MEME release changed the format and the dtypes need revisiting.
    `optional` types columns that only some runs emit (AME's depend on its
    `--method`); anything else in the file is kept as pandas infers it.
    """
    # dtype= at read time, not astype afterwards: inference would already have
    # turned a motif ID of "007" into 7, and has nothing to infer from on a
    # header-only table (every column comes back `object`).
    df = pd.read_csv(tsv, sep="\t", comment="#", dtype={**(optional or {}), **dtypes})
    missing = [c for c in dtypes if c not in df.columns]
    if missing:
        raise ValueError(
            f"{tsv} is not a {tool} results table: missing column(s) {missing}; "
            f"found {list(df.columns)}"
        )
    return df


def read_sea(tsv: str | os.PathLike) -> pd.DataFrame:
    """Read a `sea.tsv` into a DataFrame typed per `SEA_DTYPES`.

    One row per motif that PASSED the run's reporting threshold, best first.
    An empty frame means no motif passed, and still has every column.
    """
    return read_enrichment_tsv(tsv, SEA_DTYPES, "SEA")
