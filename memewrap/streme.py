"""STREME: discriminative motif discovery, primary vs control.

## Why this wrapper exists

STREME replaced DREME as the MEME suite's discriminative discovery tool. The
Python ecosystem has not followed: `gimmemotifs` 0.18.4 ships wrappers for 22
discovery tools including `meme.py`, `memew.py` and `dreme.py` -- but no
`streme.py`. `pymemesuite` 0.1.0a4 binds FIMO and the MEME motif parsers and
does not mention STREME at all. So the deprecated predecessor is wrapped
everywhere and its successor is wrapped nowhere.

## The inert-parameter trap this refuses

Two call sites in the source repo differed, and one looked much more careful:

    # scripts/denovo_motif_discovery.py
    --order 2 --thresh 0.05 --nmotifs N --minw .. --maxw ..
    # scripts/cg_grammar_sharing.py
    --nmotifs 15

Measured against the installed STREME 5.5.9 help text, the "careful" one is
passing two flags that are already the defaults (`--order` defaults to 2 for
DNA; `--thresh` defaults to 0.05) -- so the two are behaviourally identical
apart from the width bounds. The apparent divergence was cosmetic.

The real problem is in that same help text:

    --nmotifs <nmotifs>  stop if <nmotifs> motifs have been output;
                         OVERRIDES --thresh if > 0

Passing both, as the source does, makes `--thresh` INERT. The command reads as
"significant motifs, up to N of them" and actually means "exactly N motifs,
significant or not" -- STREME will happily emit its Nth motif at p=0.9. Nothing
in the output announces which rule applied.

So `run_streme` takes them as mutually exclusive and raises. You choose a
count-based run or a significance-based run, and the flag you passed is the one
that acts.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

from .tools import find_tool

__all__ = ["StremeError", "run_streme"]


class StremeError(RuntimeError):
    """STREME exited non-zero. Carries its stderr, which is where it says why."""


def run_streme(
    primary: str | os.PathLike,
    control: str | os.PathLike,
    outdir: str | os.PathLike,
    *,
    minw: int = 8,
    maxw: int = 15,
    nmotifs: int | None = None,
    thresh: float | None = None,
    order: int | None = None,
    dna: bool = True,
    extra: list[str] | None = None,
    meme_bin: str | os.PathLike | None = None,
) -> Path:
    """Run STREME on `primary` against `control`; return the streme.txt path.

    Exactly one of `nmotifs` / `thresh` may be set -- see the module docstring;
    STREME silently ignores `thresh` when `nmotifs > 0`. Passing neither uses
    STREME's own default stopping rule (significance at 0.05).

    Args:
        primary: FASTA of sequences to find enriched motifs IN.
        control: FASTA of background sequences to discriminate AGAINST.
        outdir: output directory (`--oc`, created if absent).
        minw/maxw: motif width bounds.
        nmotifs: stop after this many motifs. Mutually exclusive with `thresh`.
        thresh: p-value significance threshold. Mutually exclusive with `nmotifs`.
        order: background Markov order. STREME's default is 2 for DNA.
        extra: additional raw flags appended verbatim.

    Raises:
        ValueError: both `nmotifs` and `thresh` given, or a width bound is invalid.
        StremeError: STREME exited non-zero.
    """
    if nmotifs is not None and thresh is not None:
        raise ValueError(
            "nmotifs and thresh are mutually exclusive: STREME's --nmotifs "
            "overrides --thresh when > 0, so passing both makes thresh inert "
            "and the run stops at a fixed count regardless of significance. "
            "Pass one."
        )
    if nmotifs is not None and nmotifs <= 0:
        raise ValueError(f"nmotifs must be positive, got {nmotifs}")
    # STREME's documented floor. Catching it here turns a subprocess failure
    # several seconds in into an immediate, readable error.
    if minw < 3:
        raise ValueError(f"minw must be >= 3 (STREME's floor), got {minw}")
    if maxw < minw:
        raise ValueError(f"maxw ({maxw}) is below minw ({minw})")

    primary, control = Path(primary), Path(control)
    for label, p in (("primary", primary), ("control", control)):
        if not p.is_file():
            raise FileNotFoundError(f"{label} FASTA does not exist: {p}")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(find_tool("streme", meme_bin)),
        "--p",
        str(primary),
        "--n",
        str(control),
        "--oc",
        str(outdir),
        "--minw",
        str(minw),
        "--maxw",
        str(maxw),
    ]
    if dna:
        cmd.append("--dna")
    if order is not None:
        cmd += ["--order", str(order)]
    if nmotifs is not None:
        cmd += ["--nmotifs", str(nmotifs)]
    if thresh is not None:
        cmd += ["--thresh", str(thresh)]
    if extra:
        cmd += list(extra)

    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        # The source used check=True. CalledProcessError's message is just the
        # command and the exit code -- STREME's actual complaint ("sequences too
        # short", "alphabet mismatch") sits in .stderr and never gets printed,
        # so the failure reads as unexplained.
        raise StremeError(
            f"streme exited {proc.returncode}\n"
            f"  command: {' '.join(cmd)}\n"
            f"  stderr: {proc.stderr.strip()[-2000:]}"
        )

    out = outdir / "streme.txt"
    if not out.is_file():
        # Exit 0 with no output file. Postcondition, not paranoia: the caller's
        # next step parses this path, and a missing file there surfaces as a
        # confusing parse error rather than "STREME produced nothing".
        raise StremeError(
            f"streme exited 0 but wrote no {out.name} in {outdir} "
            f"(contents: {sorted(p.name for p in outdir.iterdir())})"
        )
    return out
