"""TOMTOM: compare discovered motifs against a reference set.

Not wrapped by `gimmemotifs` (which ships its own motif-similarity code) nor by
`pymemesuite`. The count form below -- "how many of my query motifs match
anything in this target set" -- is the shape a permutation test needs, because
it collapses a whole comparison to one integer that can be recomputed against
each null set.

## The output directory collision

The source's comment records a bug already fixed once:

    # unique out dir per (query, target) so the observed comparison's match
    # detail is not overwritten by subsequent null comparisons

TOMTOM writes to a directory, not a file. A permutation test runs one observed
comparison and then N null comparisons with the SAME query, so a naming scheme
keyed on the query alone has every null overwrite the observed run's detail.
The count returned is still correct -- it is read before the next call -- so
nothing fails; you simply find the wrong `tomtom.tsv` on disk afterwards when
you go to inspect which motifs matched. The default here derives the directory
from both query and target.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pandas as pd

from .tools import find_tool

__all__ = ["TomtomError", "run_tomtom", "match_count"]


class TomtomError(RuntimeError):
    """TOMTOM exited non-zero. Carries its stderr."""


def run_tomtom(
    query: str | os.PathLike,
    targets: str | os.PathLike | list,
    outdir: str | os.PathLike,
    *,
    thresh: float = 0.5,
    dist: str = "ed",
    min_overlap: int = 4,
    no_ssc: bool = False,
    meme_bin: str | os.PathLike | None = None,
) -> Path:
    """Run TOMTOM; return the path to tomtom.tsv.

    Args:
        query: MEME-format motif file to search WITH.
        targets: one path or a list of paths to search AGAINST.
        outdir: output directory (`-oc`).
        thresh: TOMTOM's own reporting threshold (q-value unless `-evalue`).
        dist: column comparison function (`ed`, `pearson`, `sandelin`, ...).
        min_overlap: minimum overlapping columns for a comparison.
        no_ssc: pass `-no-ssc` (disable small-sample correction).
    """
    query = Path(query)
    if not query.is_file():
        raise FileNotFoundError(f"query motif file does not exist: {query}")
    target_list = [
        Path(t)
        for t in ([targets] if isinstance(targets, (str, os.PathLike)) else targets)
    ]
    if not target_list:
        raise ValueError("no target motif files given")
    for t in target_list:
        if not t.is_file():
            raise FileNotFoundError(f"target motif file does not exist: {t}")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [str(find_tool("tomtom", meme_bin)), "-oc", str(outdir)]
    if no_ssc:
        cmd.append("-no-ssc")
    cmd += [
        "-thresh",
        str(thresh),
        "-dist",
        dist,
        "-min-overlap",
        str(min_overlap),
        str(query),
        *[str(t) for t in target_list],
    ]

    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise TomtomError(
            f"tomtom exited {proc.returncode}\n"
            f"  command: {' '.join(cmd)}\n"
            f"  stderr: {proc.stderr.strip()[-2000:]}"
        )

    tsv = outdir / "tomtom.tsv"
    if not tsv.is_file():
        raise TomtomError(f"tomtom exited 0 but wrote no tomtom.tsv in {outdir}")
    return tsv


def read_tomtom(tsv: str | os.PathLike) -> pd.DataFrame:
    """Read a tomtom.tsv into a DataFrame, without its trailing comment block.

    TOMTOM appends free-text provenance after the data, separated by a blank
    line -- verified against real 5.5.9 output, which ends:

        PLANTED1\\tPLANTED1\\t0\\t3.05e-05\\t...\\t+\\n
        \\n
        # Tomtom (Motif Comparison Tool): Version 5.5.9 ...
        # The format of this file is described at ...

    `comment='#'` drops the hash lines and pandas' default `skip_blank_lines`
    drops the separator, so this is a plain read. The source additionally called
    `dropna(subset=["Query_ID"])` to clear a phantom NaN row; mutation testing
    showed that removing it changes nothing on real output, so it is not carried
    over. `test_the_trailing_comment_block_is_dropped` still asserts the
    phantom row is absent -- the behaviour is what matters, and pandas' handling
    of it is a dependency detail that could change.
    """
    return pd.read_csv(tsv, sep="\t", comment="#")


def match_count(
    query: str | os.PathLike,
    target: str | os.PathLike,
    *,
    q_thresh: float = 0.1,
    tomtom_thresh: float = 1.0,
    outdir: str | os.PathLike | None = None,
    meme_bin: str | os.PathLike | None = None,
) -> int:
    """Number of DISTINCT query motifs with a match at q < `q_thresh`.

    Counts query motifs, not matches: a motif matching six targets contributes
    one. That is the quantity a permutation test over target sets compares,
    since the alternative (raw match rows) scales with how large the target set
    happens to be rather than with how much of the query is explained.

    ## Why TOMTOM is run permissively and the filtering happens here

    The source passed `q_thresh` as TOMTOM's own `-thresh` AND re-filtered the
    resulting q-values in pandas at the same value. The second filter could
    therefore never remove a row -- mutation testing confirmed it: deleting it
    changed no result. Two filters at one threshold are one filter and one piece
    of code that looks like it is doing something.

    Worse, it makes `q_thresh` decide what lands on disk. A permutation test
    that later wants a stricter threshold has to re-run every comparison,
    because the detail file only ever held what the first threshold admitted.

    So TOMTOM runs at `tomtom_thresh` (permissive by default) and this function
    applies `q_thresh`. The written `tomtom.tsv` is then the complete
    comparison, re-filterable without re-running, and the pandas filter is the
    one that decides the count.

    `outdir` defaults to a directory beside the query, named for the target, so
    repeated calls with one query and many targets do not overwrite each other.
    """
    if tomtom_thresh < q_thresh:
        raise ValueError(
            f"tomtom_thresh ({tomtom_thresh}) is below q_thresh ({q_thresh}): "
            f"TOMTOM would discard matches this function is meant to count, and "
            f"the count would silently be a floor rather than the answer"
        )
    query, target = Path(query), Path(target)
    if outdir is None:
        # Name from the target's PARENT directory when the file itself is a
        # generic `streme.txt` -- which it always is, so keying on the filename
        # alone would collide for every target.
        tag = target.parent.name or target.stem
        outdir = query.parent / f"tomtom_vs_{tag}"

    tsv = run_tomtom(
        query, target, outdir, thresh=tomtom_thresh, no_ssc=True, meme_bin=meme_bin
    )
    df = read_tomtom(tsv)
    if df.empty:
        return 0
    return int(df.loc[df["q-value"] < q_thresh, "Query_ID"].nunique())
