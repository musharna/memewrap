"""Guards for the TOMTOM wrapper and the match-count statistic."""

from __future__ import annotations

import pytest

from conftest import requires_meme
from memewrap.tomtom import match_count, read_tomtom, run_tomtom

TT_HEADER = (
    "Query_ID\tTarget_ID\tOptimal_offset\tp-value\tE-value\tq-value\t"
    "Overlap\tQuery_consensus\tTarget_consensus\tOrientation\n"
)


def _tt_row(q, t, q_value):
    return f"{q}\t{t}\t0\t1e-06\t1e-04\t{q_value}\t8\tTGTCTCTC\tTGTCTCTC\t+\n"


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------


def test_the_trailing_comment_block_is_dropped(tmp_path):
    """TOMTOM appends free-text provenance after a blank line. Read naively, the
    blank line becomes a row of NaNs and the file has one phantom motif."""
    p = tmp_path / "tomtom.tsv"
    p.write_text(
        TT_HEADER
        + _tt_row("M1", "T1", 0.01)
        + "\n"
        + "# Tomtom (Motif Comparison Tool): Version 5.5.9\n"
        + "# Command line: tomtom -oc out query.meme target.meme\n"
    )
    df = read_tomtom(p)
    assert len(df) == 1
    assert df.iloc[0]["Query_ID"] == "M1"


def test_a_header_only_file_reads_as_no_matches(tmp_path):
    p = tmp_path / "tomtom.tsv"
    p.write_text(TT_HEADER)
    assert read_tomtom(p).empty


# --------------------------------------------------------------------------
# Real execution
# --------------------------------------------------------------------------


@requires_meme
@pytest.mark.meme
def test_a_motif_matches_ITSELF(tmp_path, motif_db):
    """The positive control for the whole comparison path: a motif compared
    against itself must match. If this fails nothing else about match_count
    means anything, because zero is also what a broken invocation returns."""
    assert match_count(motif_db, motif_db, q_thresh=0.5, outdir=tmp_path / "self") == 1


@requires_meme
@pytest.mark.meme
def test_an_unrelated_motif_does_NOT_match(tmp_path, motif_db):
    """The negative half. Without it, match_count could return 1 unconditionally
    -- or the q-value filter could be inert -- and the self-match test above
    would still pass."""
    other = tmp_path / "other.meme"
    rows = []
    for base in "AAAACCCC":
        rows.append(" ".join("0.997" if b == base else "0.001" for b in "ACGT"))
    other.write_text(
        "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
        "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"
        "MOTIF OTHER1 other\n"
        "letter-probability matrix: alength= 4 w= 8 nsites= 100 E= 1e-30\n"
        + "\n".join(rows)
        + "\n"
    )
    n = match_count(motif_db, other, q_thresh=0.001, outdir=tmp_path / "cross")
    assert n == 0, f"an unrelated motif matched at q<0.001 ({n} matches)"


@requires_meme
@pytest.mark.meme
def test_the_default_outdir_does_not_collide_across_targets(tmp_path, motif_db):
    """A permutation test runs one query against many targets. Keying the output
    directory on the query alone -- which the source's comment records as an
    already-fixed bug -- has every null run overwrite the observed run's detail.
    The COUNT stays right, so nothing fails; the evidence on disk is just wrong.
    """
    a = tmp_path / "targ_a"
    b = tmp_path / "targ_b"
    a.mkdir()
    b.mkdir()
    (a / "streme.txt").write_bytes(motif_db.read_bytes())
    (b / "streme.txt").write_bytes(motif_db.read_bytes())
    query = tmp_path / "query.meme"
    query.write_bytes(motif_db.read_bytes())

    match_count(query, a / "streme.txt", q_thresh=0.5)
    match_count(query, b / "streme.txt", q_thresh=0.5)

    assert (tmp_path / "tomtom_vs_targ_a" / "tomtom.tsv").is_file()
    assert (tmp_path / "tomtom_vs_targ_b" / "tomtom.tsv").is_file()


@requires_meme
@pytest.mark.meme
def test_the_q_threshold_is_applied_HERE_not_only_by_tomtom(tmp_path, motif_db):
    """The filter must be able to remove something.

    The source passed q_thresh as TOMTOM's own -thresh and then re-filtered at
    the same value, so the second filter could never drop a row -- mutation
    testing confirmed deleting it changed no result. This asserts the
    python-side filter is load-bearing: the self-match sits at q ~ 6.1e-05,
    TOMTOM reports it under the permissive default, and a stricter q_thresh
    must still exclude it.
    """
    permissive = match_count(
        motif_db, motif_db, q_thresh=0.5, outdir=tmp_path / "loose"
    )
    assert permissive == 1, "precondition: the match is reported at all"

    strict = match_count(motif_db, motif_db, q_thresh=1e-9, outdir=tmp_path / "strict")
    assert strict == 0, (
        "a q_thresh far below the observed q-value did not exclude the match; "
        "the python-side filter is inert"
    )


def test_a_tomtom_threshold_below_the_q_threshold_is_refused(tmp_path, motif_db):
    """Otherwise TOMTOM discards rows this function is supposed to count, and the
    returned number is a floor rather than the answer -- with nothing to
    distinguish 'no matches' from 'matches thrown away upstream'."""
    with pytest.raises(ValueError, match="below q_thresh"):
        match_count(motif_db, motif_db, q_thresh=0.5, tomtom_thresh=0.1)


@requires_meme
@pytest.mark.meme
def test_the_written_tsv_keeps_matches_beyond_the_q_threshold(tmp_path, motif_db):
    """The point of running TOMTOM permissively: the detail file stays complete,
    so a stricter threshold can be applied later without re-running."""
    match_count(motif_db, motif_db, q_thresh=1e-9, outdir=tmp_path / "kept")
    df = read_tomtom(tmp_path / "kept" / "tomtom.tsv")
    assert len(df) == 1, (
        "the match was excluded from the file as well as the count; "
        "re-thresholding would require re-running every comparison"
    )


@requires_meme
@pytest.mark.meme
def test_run_tomtom_returns_a_parseable_tsv(tmp_path, motif_db):
    tsv = run_tomtom(motif_db, motif_db, tmp_path / "out", thresh=0.5)
    assert tsv.is_file() and tsv.name == "tomtom.tsv"
    df = read_tomtom(tsv)
    assert "Query_ID" in df.columns and "q-value" in df.columns


def test_a_missing_query_file_is_named(tmp_path, motif_db):
    with pytest.raises(FileNotFoundError, match="query"):
        run_tomtom(tmp_path / "absent.meme", motif_db, tmp_path / "out")


def test_a_missing_target_file_is_named(tmp_path, motif_db):
    with pytest.raises(FileNotFoundError, match="target"):
        run_tomtom(motif_db, tmp_path / "absent.meme", tmp_path / "out")


def test_an_empty_target_list_is_refused(tmp_path, motif_db):
    with pytest.raises(ValueError, match="no target"):
        run_tomtom(motif_db, [], tmp_path / "out")
