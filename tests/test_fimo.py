"""Guards for the FIMO wrapper, chiefly the header bug the source shipped.

`test_the_SOURCE_merge_loses_the_header_when_chunk_zero_is_empty` reproduces the
original algorithm verbatim and asserts it FAILS. That is the "seen it fail"
control made permanent: without it, the fixed merge would pass its own test on
day one and nobody could tell whether the fixture could ever have caught the
bug. With it, a regression to the old loop turns two tests red -- one for the
wrong output, one because the deliberately-broken control started passing.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from conftest import PLANTED_MOTIF, requires_meme
from memewrap.fimo import (
    build_feature_matrix,
    merge_fimo_tsvs,
    run_fimo_parallel,
    split_fasta,
)

HEADER = "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tmatched_sequence\n"


def _row(motif, seq, score):
    return f"{motif}\tplanted\t{seq}\t1\t8\t+\t{score}\t1e-05\tTGTCTCTC\n"


# --------------------------------------------------------------------------
# The merge -- pure logic, no binary
# --------------------------------------------------------------------------


def _source_merge(parts, merged):
    """The ORIGINAL algorithm, reproduced exactly for use as a failing control.

    for i, o in enumerate(outs):
        for j, line in enumerate(fh):
            if j == 0 and i > 0:   # keep only the first chunk's header
                continue
            out.write(line)
    """
    with open(merged, "w") as out:
        for i, o in enumerate(parts):
            with open(o) as fh:
                for j, line in enumerate(fh):
                    if j == 0 and i > 0:
                        continue
                    out.write(line)
    return merged


def test_a_normal_merge_keeps_exactly_one_header(tmp_path):
    a = tmp_path / "a.txt"
    b = tmp_path / "b.txt"
    a.write_text(HEADER + _row("M1", "seq0", 10.0))
    b.write_text(HEADER + _row("M1", "seq1", 11.0))
    out = merge_fimo_tsvs([a, b], tmp_path / "m.tsv")
    lines = out.read_text().splitlines()
    assert lines[0].startswith("motif_id")
    assert sum(ln.startswith("motif_id") for ln in lines) == 1
    assert len(lines) == 3


def test_an_EMPTY_first_chunk_does_not_cost_the_header(tmp_path):
    """THE bug. FIMO writes no header when a chunk yields no hits, and chunks are
    dealt round-robin, so an empty chunk 0 is routine rather than exotic."""
    a = tmp_path / "a.txt"
    b = tmp_path / "b.txt"
    a.write_text("")  # scanned, found nothing, wrote nothing
    b.write_text(HEADER + _row("M1", "seq1", 11.0))
    out = merge_fimo_tsvs([a, b], tmp_path / "m.tsv")
    lines = out.read_text().splitlines()
    assert lines[0].startswith("motif_id"), (
        f"merged file starts with {lines[0]!r} -- the header was dropped"
    )
    assert len(lines) == 2, "the data row was lost along with the header"


def test_the_SOURCE_merge_loses_the_header_when_chunk_zero_is_empty(tmp_path):
    """The broken-state control: proof this fixture can detect the bug.

    The consequence is asserted the way a caller meets it -- pandas silently
    promotes the first DATA row to column names, so the scan appears to have
    found one hit fewer and every column is mislabelled. No exception is raised.
    """
    a = tmp_path / "a.txt"
    b = tmp_path / "b.txt"
    a.write_text("")
    b.write_text(HEADER + _row("M1", "seq1", 11.0) + _row("M1", "seq2", 12.0))
    broken = _source_merge([a, b], tmp_path / "broken.tsv")

    first = broken.read_text().splitlines()[0]
    assert not first.startswith("motif_id"), (
        "the source algorithm kept the header -- this control no longer "
        "reproduces the bug, so the fixture proves nothing"
    )

    df = pd.read_csv(broken, sep="\t")
    assert "motif_id" not in df.columns, "expected the mislabelled-columns failure"
    assert len(df) == 1, (
        "the promoted data row is consumed as the header: 2 hits read back as 1"
    )


def test_a_chunk_without_a_header_is_still_included(tmp_path):
    """Not every part is guaranteed to carry a header line; body-only parts must
    contribute their rows rather than lose their first one."""
    a = tmp_path / "a.txt"
    b = tmp_path / "b.txt"
    a.write_text(HEADER + _row("M1", "seq0", 10.0))
    b.write_text(_row("M1", "seq1", 11.0))
    out = merge_fimo_tsvs([a, b], tmp_path / "m.tsv")
    assert len(out.read_text().splitlines()) == 3


def test_merging_nothing_yields_an_empty_file(tmp_path):
    out = merge_fimo_tsvs([], tmp_path / "m.tsv")
    assert out.read_text() == ""


# --------------------------------------------------------------------------
# Chunking
# --------------------------------------------------------------------------


def test_splitting_preserves_every_record(tmp_path, small_fasta):
    parts = split_fasta(small_fasta, tmp_path / "chunks", 4)
    total = sum(p.read_text().count(">") for p in parts)
    assert total == small_fasta.read_text().count(">") == 12


def test_more_chunks_than_sequences_yields_no_EMPTY_chunks(tmp_path):
    """`recs[i::n]` produces empty lists once n exceeds the record count, and
    launching FIMO on an empty FASTA is a subprocess started to do nothing."""
    from conftest import make_corpus, write_fasta

    fa = write_fasta(tmp_path / "three.fa", make_corpus(3, 60, seed=7))
    parts = split_fasta(fa, tmp_path / "chunks", 10)
    assert len(parts) == 3
    assert all(p.read_text().count(">") > 0 for p in parts)


def test_an_empty_fasta_is_refused(tmp_path):
    fa = tmp_path / "empty.fa"
    fa.write_text("")
    with pytest.raises(ValueError, match="no FASTA records"):
        split_fasta(fa, tmp_path / "chunks", 2)


# --------------------------------------------------------------------------
# build_feature_matrix
# --------------------------------------------------------------------------


def test_the_matrix_takes_the_MAX_score_per_pair(tmp_path):
    """A motif hitting one gene several times contributes its best hit, not its
    last. A fixture with one hit per pair could not tell the two apart."""
    tsv = tmp_path / "f.tsv"
    tsv.write_text(
        HEADER
        + _row("M1", "g1", 3.0)
        + _row("M1", "g1", 9.0)
        + _row("M1", "g1", 5.0)
        + _row("M2", "g2", 2.0)
    )
    X = build_feature_matrix(tsv, ["g1", "g2"], ["M1", "M2"])
    assert X[0, 0] == 9.0
    assert X[1, 1] == 2.0


def test_rows_and_columns_follow_the_CALLERS_order(tmp_path):
    """Reordering to whatever pandas grouped by would silently misalign the
    caller's label vector against the matrix."""
    tsv = tmp_path / "f.tsv"
    tsv.write_text(HEADER + _row("M2", "g2", 7.0) + _row("M1", "g1", 4.0))
    X = build_feature_matrix(tsv, ["g1", "g2"], ["M1", "M2"])
    assert X[0, 0] == 4.0 and X[1, 1] == 7.0
    assert X[0, 1] == 0.0 and X[1, 0] == 0.0


def test_genes_and_motifs_with_no_hits_become_zero_rows_not_missing_ones(tmp_path):
    tsv = tmp_path / "f.tsv"
    tsv.write_text(HEADER + _row("M1", "g1", 4.0))
    X = build_feature_matrix(tsv, ["g1", "ghost"], ["M1", "unseen"])
    assert X.shape == (2, 2)
    assert np.array_equal(X, np.array([[4.0, 0.0], [0.0, 0.0]]))


def test_ids_absent_from_the_caller_lists_are_ignored(tmp_path):
    """A scan covers the whole FASTA; the caller may model a subset."""
    tsv = tmp_path / "f.tsv"
    tsv.write_text(HEADER + _row("M1", "g1", 4.0) + _row("M1", "not_modelled", 99.0))
    X = build_feature_matrix(tsv, ["g1"], ["M1"])
    assert X.shape == (1, 1) and X[0, 0] == 4.0


def test_a_scan_with_no_hits_gives_a_zero_matrix_of_the_right_shape(tmp_path):
    tsv = tmp_path / "f.tsv"
    tsv.write_text(HEADER)
    X = build_feature_matrix(tsv, ["g1", "g2"], ["M1"])
    assert X.shape == (2, 1) and not X.any()


# --------------------------------------------------------------------------
# Real execution
# --------------------------------------------------------------------------


@requires_meme
@pytest.mark.meme
def test_fimo_finds_the_planted_motif_across_chunks(tmp_path, small_fasta, motif_db):
    """Real FIMO, really chunked, really merged.

    Every sequence carries the planted motif, so a correct run reports a hit for
    all 12 -- which also proves the merge did not drop a chunk. The unit tests
    above use hand-written TSVs; only this one proves the file FIMO actually
    emits still parses with the same header assumptions.
    """
    merged = run_fimo_parallel(
        small_fasta, motif_db, tmp_path / "fimo", n_chunks=4, thresh="1e-3"
    )
    df = pd.read_csv(merged, sep="\t")
    assert "motif_id" in df.columns, f"header missing; got {list(df.columns)[:3]}"
    assert "sequence_name" in df.columns

    found = set(df["sequence_name"].astype(str))
    expected = {f"seq{i}" for i in range(12)}
    assert found == expected, f"missing sequences: {sorted(expected - found)}"

    assert df["matched_sequence"].str.contains(PLANTED_MOTIF).any()


@requires_meme
@pytest.mark.meme
def test_the_real_scan_feeds_build_feature_matrix(tmp_path, small_fasta, motif_db):
    """Closes the loop: the matrix is built from a real FIMO TSV, not a fixture."""
    merged = run_fimo_parallel(
        small_fasta, motif_db, tmp_path / "fimo", n_chunks=3, thresh="1e-3"
    )
    genes = [f"seq{i}" for i in range(12)]
    X = build_feature_matrix(merged, genes, ["PLANTED1"])
    assert X.shape == (12, 1)
    assert (X[:, 0] > 0).all(), (
        f"{int((X[:, 0] == 0).sum())} of 12 sequences scored 0 despite carrying "
        f"the motif"
    )
