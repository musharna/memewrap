"""Guards for the SEA wrapper.

The real-execution tests share one corpus with a known answer: PLANTED1 sits in
90% of the primary sequences and none of the controls, DECOY1 is planted
nowhere. Every enrichment assertion has both halves -- the planted motif IS
significant and the decoy IS NOT -- because a reader that mislabelled a column,
or a wrapper that lost the control set, tends to break only one of them.
"""

from __future__ import annotations

import warnings
from contextlib import suppress

import pytest

from conftest import requires_enrichment, write_fasta
from memewrap import MemeToolNotFound, SeaError, read_sea, run_sea
from memewrap.sea import SEA_DTYPES

# Verbatim from `sea --text` 5.5.9 on the fixture corpus, provenance block and
# the blank line before it included.
REAL_SEA_TSV = (
    "RANK\tDB\tID\tALT_ID\tCONSENSUS\tTP\tTP%\tFP\tFP%\tENR_RATIO\tSCORE_THR\t"
    "PVALUE\tLOG_PVALUE\tEVALUE\tLOG_EVALUE\tQVALUE\tLOG_QVALUE\n"
    "1\tdb.meme\tPLANTED1\talt_PLANTED1\tTGTCTCTC\t156\t86.67\t0\t0.00\t157\t16\t"
    "1.08e-76\t-174.92\t2.16e-76\t-174.23\t2.16e-76\t-174.23\n"
    "2\tdb.meme\tDECOY1\talt_DECOY1\tGATTACAG\t12\t6.67\t24\t13.33\t0.52\t-1e-06\t"
    "9.89e-1\t-0.01\t1.98e0\t0.68\t9.89e-1\t-0.01\n"
    "\n"
    "# SEA (Simple Enrichment Analysis): Version 5.5.9 compiled on May 27 2026 at 12:41:20\n"
    "# The format of this file is described at https://meme-suite.org/meme/doc/sea-output-format.html.\n"
    "# sea --text --p p.fa --n c.fa --m db.meme\n"
)
HEADER = REAL_SEA_TSV.splitlines(keepends=True)[0]


# --------------------------------------------------------------------------
# Parsing -- no binary needed
# --------------------------------------------------------------------------


def test_the_statistics_land_in_their_own_columns_with_their_own_dtypes(tmp_path):
    p = tmp_path / "sea.tsv"
    p.write_text(REAL_SEA_TSV)
    df = read_sea(p)

    assert list(df["ID"]) == ["PLANTED1", "DECOY1"], (
        "the trailing comment block or its blank separator became a row"
    )
    assert {c: str(t) for c, t in df.dtypes.items()} == SEA_DTYPES

    # The three differ on the decoy row (p=0.989, E=1.98, q=0.989 -- and E != p),
    # so reading any one from a neighbouring column is caught here.
    decoy = df.iloc[1]
    assert decoy["PVALUE"] == pytest.approx(0.989)
    assert decoy["EVALUE"] == pytest.approx(1.98)
    assert decoy["QVALUE"] == pytest.approx(0.989)
    assert decoy["LOG_EVALUE"] == pytest.approx(0.68)
    # abs=0: pytest.approx's default absolute tolerance (1e-12) makes
    # 1.08e-76 "equal" 2.16e-76, and a PVALUE/QVALUE swap survived that.
    assert df.iloc[0]["PVALUE"] == pytest.approx(1.08e-76, rel=1e-9, abs=0)
    assert df.iloc[0]["QVALUE"] == pytest.approx(2.16e-76, rel=1e-9, abs=0)
    # ENR_RATIO and SCORE_THR print as bare integers on this row; they are
    # floats, and inference from this row alone would say int64.
    assert df.iloc[0]["ENR_RATIO"] == 157.0 and df["SCORE_THR"].dtype == "float64"
    assert df.iloc[0]["TP"] == 156 and df["TP"].dtype == "int64"


def test_a_header_only_table_is_an_empty_frame_with_typed_columns(tmp_path):
    """No motif passing the threshold is a result. It has to survive the filter
    a caller will apply next, which an all-`object` empty frame does not
    promise."""
    p = tmp_path / "sea.tsv"
    p.write_text(HEADER + "\n# SEA (Simple Enrichment Analysis): Version 5.5.9\n")
    df = read_sea(p)
    assert df.empty
    assert {c: str(t) for c, t in df.dtypes.items()} == SEA_DTYPES
    assert df[df["QVALUE"] < 0.05].empty


def test_a_numeric_looking_motif_id_stays_a_string(tmp_path):
    p = tmp_path / "sea.tsv"
    row = REAL_SEA_TSV.splitlines(keepends=True)[1].replace("PLANTED1", "007", 1)
    assert "\t007\t" in row  # the fixture edit took
    p.write_text(HEADER + row)
    assert list(read_sea(p)["ID"]) == ["007"]


def test_a_table_missing_a_declared_column_is_refused_by_name(tmp_path):
    good = tmp_path / "good.tsv"
    good.write_text(REAL_SEA_TSV)
    assert len(read_sea(good)) == 2  # positive control: the intact table reads

    bad = tmp_path / "bad.tsv"
    bad.write_text(REAL_SEA_TSV.replace("\tQVALUE\t", "\tQ_VALUE\t", 1))
    with pytest.raises(ValueError, match="QVALUE") as e:
        read_sea(bad)
    assert "SEA" in str(e.value)


# --------------------------------------------------------------------------
# Argument validation -- no binary needed
# --------------------------------------------------------------------------


def _tiny(tmp_path):
    p = write_fasta(tmp_path / "p.fa", [("a", "ACGT" * 20)])
    c = write_fasta(tmp_path / "c.fa", [("b", "TGCA" * 20)])
    return p, c


def test_an_unknown_threshold_statistic_is_refused(tmp_path, motif_db):
    p, c = _tiny(tmp_path)
    with pytest.raises(ValueError, match="threshold_on"):
        run_sea(p, motif_db, tmp_path / "o", control=c, threshold_on="fdr")
    # Positive control: each legal value gets past validation. Launching SEA
    # (or not finding it) is the outcome checked for; a ValueError is NOT
    # suppressed, so refusing a legal value fails here.
    for legal in ("evalue", "qvalue", "pvalue"):
        with suppress(SeaError, MemeToolNotFound):
            run_sea(p, motif_db, tmp_path / "o", control=c, threshold_on=legal)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"thresh": 0},
        {"thresh": -1.0},
        {"thresh": 1.5, "threshold_on": "qvalue"},
        {"thresh": 1.5, "threshold_on": "pvalue"},
    ],
)
def test_an_out_of_range_threshold_is_refused_before_launching(
    tmp_path, motif_db, kwargs
):
    p, c = _tiny(tmp_path)
    with pytest.raises(ValueError, match="thresh"):
        run_sea(p, motif_db, tmp_path / "o", control=c, **kwargs)
    # Positive control: 1.5 is a legal E-value threshold.
    with suppress(SeaError, MemeToolNotFound):
        run_sea(p, motif_db, tmp_path / "o", control=c, thresh=1.5)


def test_each_missing_input_is_named(tmp_path, motif_db):
    p, c = _tiny(tmp_path)
    absent = tmp_path / "absent"
    with pytest.raises(FileNotFoundError, match="primary"):
        run_sea(absent, motif_db, tmp_path / "o", control=c)
    with pytest.raises(FileNotFoundError, match="control"):
        run_sea(p, motif_db, tmp_path / "o", control=absent)
    with pytest.raises(FileNotFoundError, match="motif"):
        run_sea(p, [motif_db, absent], tmp_path / "o", control=c)
    with pytest.raises(ValueError, match="no motif files"):
        run_sea(p, [], tmp_path / "o", control=c)


# --------------------------------------------------------------------------
# Real execution
# --------------------------------------------------------------------------


@requires_enrichment
@pytest.mark.meme
def test_sea_finds_the_planted_motif_enriched_and_the_decoy_not(
    tmp_path, primary_fasta, control_fasta, enrichment_db
):
    # thresh=100 (an E-value) so that the decoy is REPORTED: at a strict
    # threshold its absence would be read here as "not significant" whether or
    # not it had been tested at all.
    tsv = run_sea(
        primary_fasta,
        enrichment_db,
        tmp_path / "sea",
        control=control_fasta,
        thresh=100,
    )
    assert tsv == tmp_path / "sea" / "sea.tsv"
    df = read_sea(tsv).set_index("ID")
    assert sorted(df.index) == ["DECOY1", "PLANTED1"]

    planted, decoy = df.loc["PLANTED1"], df.loc["DECOY1"]
    assert planted["RANK"] == 1
    assert planted["PVALUE"] < 1e-20
    assert planted["EVALUE"] < 1e-20
    assert planted["QVALUE"] < 1e-20
    assert planted["TP%"] > 70 and planted["FP%"] < 5

    assert decoy["PVALUE"] > 0.05
    assert decoy["QVALUE"] > 0.05
    assert decoy["EVALUE"] > 0.05
    # E-value = p-value x number of motifs (2). Distinguishes the two columns on
    # real output, not only on the fixture above.
    assert decoy["EVALUE"] == pytest.approx(2 * decoy["PVALUE"], rel=0.02)


@requires_enrichment
@pytest.mark.meme
def test_the_control_set_reaches_sea(
    tmp_path, primary_fasta, control_fasta, planted_control_fasta, enrichment_db
):
    """Against a control carrying the same motif, it is not enriched.

    If `control` were dropped SEA would shuffle the primary sequences instead,
    and the motif would come back at p ~ 1e-79 in BOTH runs.
    """
    same = run_sea(
        primary_fasta,
        enrichment_db,
        tmp_path / "same",
        control=planted_control_fasta,
        thresh=100,
    )
    p_same = read_sea(same).set_index("ID").loc["PLANTED1", "PVALUE"]
    assert p_same > 1e-3, f"motif present in both sets came back enriched: {p_same}"

    clean = run_sea(
        primary_fasta, enrichment_db, tmp_path / "clean", control=control_fasta
    )
    assert read_sea(clean).set_index("ID").loc["PLANTED1", "PVALUE"] < 1e-20


@requires_enrichment
@pytest.mark.meme
def test_sea_without_a_control_shuffles_the_primary_under_the_given_seed(
    tmp_path, primary_fasta, enrichment_db
):
    def run(name, seed):
        tsv = run_sea(
            primary_fasta, enrichment_db, tmp_path / name, seed=seed, thresh=100
        )
        return read_sea(tsv).set_index("ID")

    df = run("a", 1)
    assert df.loc["PLANTED1", "QVALUE"] < 1e-20
    assert df.loc["DECOY1", "QVALUE"] > 0.05
    # The shuffle differs by seed and the p-value with it (3.69e-75 vs 4.85e-81
    # on 5.5.9), so a dropped `seed` shows as two seeds agreeing.
    p1, p1_again, p2 = (
        df.loc["PLANTED1", "PVALUE"],
        run("b", 1).loc["PLANTED1", "PVALUE"],
        run("c", 2).loc["PLANTED1", "PVALUE"],
    )
    assert p1 == p1_again
    assert p1 != p2


@requires_enrichment
@pytest.mark.meme
def test_the_threshold_and_its_statistic_reach_sea(
    tmp_path, primary_fasta, control_fasta, enrichment_db
):
    """The decoy has p ~ 0.99 and E ~ 1.98, so a threshold of 1.0 admits it as a
    p-value and omits it as an E-value. One value, two outcomes: that separates
    "thresh was passed" from "threshold_on was passed"."""

    def ids(name, **kw):
        tsv = run_sea(
            primary_fasta, enrichment_db, tmp_path / name, control=control_fasta, **kw
        )
        return sorted(read_sea(tsv)["ID"])

    assert ids("default") == ["DECOY1", "PLANTED1"]  # SEA's default: E <= 10
    assert ids("as_evalue", thresh=1.0) == ["PLANTED1"]
    assert ids("as_pvalue", thresh=1.0, threshold_on="pvalue") == ["DECOY1", "PLANTED1"]
    assert ids("none_pass", thresh=1e-300) == []


@requires_enrichment
@pytest.mark.meme
def test_a_sea_failure_carries_its_stderr_and_leaves_no_table(
    tmp_path, primary_fasta, control_fasta, enrichment_db
):
    out = tmp_path / "out"
    good = run_sea(primary_fasta, enrichment_db, out, control=control_fasta)
    assert good.is_file()  # positive control, and the stale table for below

    # A FASTA is not a motif file.
    with pytest.raises(SeaError) as e:
        run_sea(primary_fasta, primary_fasta, out, control=control_fasta)
    msg = str(e.value)
    assert "sea exited" in msg and "command:" in msg
    assert "motif" in msg.split("stderr:")[1].lower(), "SEA's own complaint is absent"
    assert not good.exists(), (
        "the earlier run's sea.tsv survived a failed run in the same directory"
    )


# SEA's message when its HTML template has no data section, as the bioconda
# 5.5.9 build `pl5321he99cc7f_1` ships it.
BROKEN_TEMPLATE = "Template does not contain data section"


@requires_enrichment
@pytest.mark.meme
def test_oc_mode_agrees_with_text_mode_or_fails_out_loud(
    tmp_path, primary_fasta, control_fasta, enrichment_db
):
    """`text=False` on an intact install yields the same table plus the side
    outputs. On an install with the broken template `sea --oc` exits 1 and
    leaves a header-only sea.tsv; what must hold there is that the caller gets
    SeaError with SEA's message, not that file."""
    text = read_sea(
        run_sea(primary_fasta, enrichment_db, tmp_path / "t", control=control_fasta)
    )
    try:
        tsv = run_sea(
            primary_fasta,
            enrichment_db,
            tmp_path / "oc",
            control=control_fasta,
            text=False,
        )
    except SeaError as e:
        assert BROKEN_TEMPLATE in str(e), f"sea --oc failed for another reason: {e}"
        # Said out loud so the run's summary shows WHICH branch was exercised.
        warnings.warn(
            "this MEME install's sea --oc is broken (HTML template); the "
            "text=False success path was NOT exercised",
            stacklevel=1,
        )
        return
    oc = read_sea(tsv)
    assert list(oc["ID"]) == list(text["ID"]) == ["PLANTED1", "DECOY1"]
    assert list(oc["QVALUE"]) == list(text["QVALUE"])
    # sequences.tsv, not sites.tsv: SEA 5.4.1 does not write the latter.
    assert (tmp_path / "oc" / "sequences.tsv").is_file()


def test_a_missing_sea_binary_is_a_setup_error_and_writes_nothing(
    tmp_path, monkeypatch, motif_db
):
    p, c = _tiny(tmp_path)
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.delenv("MEME_BIN", raising=False)
    with pytest.raises(MemeToolNotFound, match="'sea'"):
        run_sea(p, motif_db, tmp_path / "o", control=c)
    assert not (tmp_path / "o" / "sea.tsv").exists()
