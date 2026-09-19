"""AME: enrichment of known motifs, by a choice of statistical tests.

AME (McLeay & Bailey 2010, BMC Bioinformatics 11:165,
doi:10.1186/1471-2105-11-165) predates SEA and answers the same two-set question
with Fisher's exact test by default. Prefer `memewrap.sea` for primary-vs-control
enrichment; AME remains the tool when the sequences carry a continuous score
(`--method ranksum|pearson|spearman`) rather than a set membership.

## Leaving out the control does not mean "shuffle"

`sea` with no control shuffles the primary sequences. `ame` with no `--control`
does something else entirely: it treats the FASTA input ORDER as a ranking and
maximises over partitions of it. Measured against AME 5.5.9 on 200 sequences,
90% carrying a planted motif, in no meaningful order:

    ame --control control.fa  ...   PLANTED1  p = 4.37e-97
    ame --control --shuffle-- ...   PLANTED1  p = 1.98e-208
    ame                       ...   PLANTED1  p = 6.99e-2

The third run exits 0 and writes a well-formed table. So `control` is a REQUIRED
keyword here with no default: a path, `AME_SHUFFLE`, or an explicit `None` for
the partition-maximisation mode when the input really is ranked.

## Reporting threshold, `--text`, underflow

As for SEA (see `memewrap.sea`): motifs failing the E-value report threshold
(default 10) are OMITTED from the table; `--text` is the default because
`ame --oc` fails on MEME builds with a broken HTML template and the table is
identical either way; and p-values below ~1e-308 parse to 0.0. AME has no log
columns to fall back on, and reports no q-value: its multiple-testing columns
are `adj_p-value` (corrected for the partitions tried per motif) and `E-value`
(`adj_p-value` times the number of motifs).
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pandas as pd

from .sea import _motif_files, read_enrichment_tsv
from .tools import find_tool

__all__ = ["AmeError", "AME_DTYPES", "AME_SHUFFLE", "run_ame", "read_ame"]

# AME's own keyword for "use shuffled copies of the input as the control".
AME_SHUFFLE = "--shuffle--"


class AmeError(RuntimeError):
    """AME exited non-zero, or exited 0 without a table. Carries its stderr."""


# The columns every `--method` emits. The rest depend on the method.
AME_DTYPES: dict[str, str] = {
    "rank": "int64",
    "motif_DB": "string",
    "motif_ID": "string",
    "motif_alt_ID": "string",
    "consensus": "string",
    "p-value": "float64",
    "adj_p-value": "float64",
    "E-value": "float64",
    "tests": "int64",
}

# Emitted by the default method (fisher) with a control set.
_AME_METHOD_DTYPES: dict[str, str] = {
    "pos": "int64",
    "neg": "int64",
    "TP": "int64",
    "%TP": "float64",
    "FP": "int64",
    "%FP": "float64",
}


def run_ame(
    sequences: str | os.PathLike,
    motifs: str | os.PathLike | list,
    outdir: str | os.PathLike,
    *,
    control: str | os.PathLike | None,
    method: str | None = None,
    scoring: str | None = None,
    evalue_report_threshold: float | None = None,
    seed: int | None = None,
    text: bool = True,
    extra: list[str] | None = None,
    meme_bin: str | os.PathLike | None = None,
) -> Path:
    """Run AME; return the path to `ame.tsv` in `outdir`. Parse it with `read_ame`.

    Args:
        sequences: FASTA of sequences to test for enrichment IN.
        motifs: one MEME-format motif file, or a list of them.
        outdir: output directory, created if absent.
        control: REQUIRED. A control FASTA; `AME_SHUFFLE` to test against
            shuffled copies of `sequences`; or None for AME's partition-
            maximisation mode, which ranks by input order -- see module docstring.
        method: statistical test (`fisher`, `ranksum`, `pearson`, ...).
        scoring: sequence scoring (`avg`, `max`, `sum`, `totalhits`).
        evalue_report_threshold: None is AME's default (10). Motifs above it are
            OMITTED from the table.
        seed: random seed for the shuffle.
        text: True (default) runs `ame --text` and writes its stdout to
            `outdir/ame.tsv`. False runs `ame --oc outdir` (adds sequences.tsv
            and ame.html).
        extra: additional raw flags appended verbatim.

    Raises:
        FileNotFoundError: an input file is missing; the message says which.
        MemeToolNotFound: no `ame` executable.
        AmeError: AME exited non-zero, or produced no table.
    """
    sequences = Path(sequences)
    if not sequences.is_file():
        raise FileNotFoundError(f"sequence FASTA does not exist: {sequences}")
    shuffled = isinstance(control, str) and control == AME_SHUFFLE
    if control is not None and not shuffled and not Path(control).is_file():
        raise FileNotFoundError(f"control FASTA does not exist: {control}")
    motif_files = _motif_files(motifs)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tsv = outdir / "ame.tsv"
    ame_exe = str(find_tool("ame", meme_bin))
    tsv.unlink(missing_ok=True)  # as in run_sea: no stale table after a failure

    cmd = [ame_exe]
    cmd += ["--text"] if text else ["--oc", str(outdir)]
    if control is not None:
        cmd += ["--control", str(control)]
    if method is not None:
        cmd += ["--method", method]
    if scoring is not None:
        cmd += ["--scoring", scoring]
    if evalue_report_threshold is not None:
        cmd += ["--evalue-report-threshold", str(evalue_report_threshold)]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    if extra:
        cmd += list(extra)
    cmd += [str(sequences), *[str(m) for m in motif_files]]

    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise AmeError(
            f"ame exited {proc.returncode}\n"
            f"  command: {' '.join(cmd)}\n"
            f"  stderr: {proc.stderr.strip()[-2000:]}"
        )
    if text:
        tsv.write_text(proc.stdout)
    if not tsv.is_file() or not tsv.read_text().startswith("rank\t"):
        raise AmeError(
            f"ame exited 0 but produced no results table at {tsv}\n"
            f"  command: {' '.join(cmd)}\n"
            f"  stderr: {proc.stderr.strip()[-2000:]}"
        )
    return tsv


def read_ame(tsv: str | os.PathLike) -> pd.DataFrame:
    """Read an `ame.tsv` into a DataFrame; shared columns typed per `AME_DTYPES`.

    One row per motif that PASSED the E-value report threshold, best first. An
    empty frame means none did, and still has the file's columns.
    """
    return read_enrichment_tsv(tsv, AME_DTYPES, "AME", optional=_AME_METHOD_DTYPES)
