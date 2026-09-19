"""Guards for the AME wrapper. Same corpus and same two-sided assertions as
`test_sea.py`: PLANTED1 is enriched, DECOY1 is not."""

from __future__ import annotations

import pytest

from conftest import requires_enrichment, write_fasta
from memewrap import AME_SHUFFLE, AmeError, MemeToolNotFound, read_ame, run_ame
from memewrap.ame import AME_DTYPES

# Verbatim from `ame --text --control c.fa p.fa db.meme` 5.5.9 (method=fisher).
REAL_AME_TSV = (
    "rank\tmotif_DB\tmotif_ID\tmotif_alt_ID\tconsensus\tp-value\tadj_p-value\t"
    "E-value\ttests\tFASTA_max\tpos\tneg\tPWM_min\tTP\t%TP\tFP\t%FP\n"
    "1\tdb.meme\tPLANTED1\talt_PLANTED1\tTGTCTCTC\t4.37e-97\t8.08e-95\t1.62e-94\t"
    "185\t200\t200\t200\t222\t185\t92.50\t0\t0.00\n"
    "2\tdb.meme\tDECOY1\talt_DECOY1\tGATTACAG\t1.00e0\t1.00e0\t2.00e0\t"
    "1\t200\t200\t200\t224\t0\t0.00\t1\t0.50\n"
    "\n"
    "# AME (Analysis of Motif Enrichment): Version 5.5.9 compiled on May 27 2026 at 12:40:54\n"
    "# The format of this file is described at https://meme-suite.org/meme/doc/ame-output-format.html.\n"
    "# ame --text --control c.fa p.fa db.meme\n"
)
HEADER = REAL_AME_TSV.splitlines(keepends=True)[0]


# --------------------------------------------------------------------------
# Parsing and validation -- no binary needed
# --------------------------------------------------------------------------


def test_the_three_significance_columns_are_read_from_their_own_columns(tmp_path):
    p = tmp_path / "ame.tsv"
    p.write_text(REAL_AME_TSV)
    df = read_ame(p)
    assert list(df["motif_ID"]) == ["PLANTED1", "DECOY1"]
    assert {c: str(df[c].dtype) for c in AME_DTYPES} == AME_DTYPES

    # All three differ on the planted row, so a shifted column is caught. abs=0
    # because pytest.approx's default ABSOLUTE tolerance is 1e-12, under which
    # every p-value on this row equals every other.
    planted = df.iloc[0]
    assert planted["p-value"] == pytest.approx(4.37e-97, rel=1e-9, abs=0)
    assert planted["adj_p-value"] == pytest.approx(8.08e-95, rel=1e-9, abs=0)
    assert planted["E-value"] == pytest.approx(1.62e-94, rel=1e-9, abs=0)
    assert planted["tests"] == 185
    assert df["%TP"].dtype == "float64" and df["TP"].dtype == "int64"


def test_a_header_only_table_is_an_empty_typed_frame(tmp_path):
    p = tmp_path / "ame.tsv"
    p.write_text(HEADER + "\n# AME (Analysis of Motif Enrichment): Version 5.5.9\n")
    df = read_ame(p)
    assert df.empty
    assert {c: str(df[c].dtype) for c in AME_DTYPES} == AME_DTYPES
    assert df[df["adj_p-value"] < 0.05].empty


def test_a_sea_table_is_not_an_ame_table(tmp_path):
    good = tmp_path / "ame.tsv"
    good.write_text(REAL_AME_TSV)
    assert len(read_ame(good)) == 2  # positive control

    bad = tmp_path / "sea.tsv"
    bad.write_text("RANK\tDB\tID\tPVALUE\n1\tdb\tM1\t0.5\n")
    with pytest.raises(ValueError, match="adj_p-value"):
        read_ame(bad)


def test_control_has_no_default(tmp_path, motif_db):
    """AME without --control ranks by input order -- see the module docstring.
    That mode is reachable with an explicit None, never by omission."""
    s = write_fasta(tmp_path / "s.fa", [("a", "ACGT" * 20)])
    with pytest.raises(TypeError, match="control"):
        run_ame(s, motif_db, tmp_path / "o")  # type: ignore[call-arg]


def test_each_missing_input_is_named(tmp_path, motif_db):
    s = write_fasta(tmp_path / "s.fa", [("a", "ACGT" * 20)])
    absent = tmp_path / "absent"
    with pytest.raises(FileNotFoundError, match="sequence"):
        run_ame(absent, motif_db, tmp_path / "o", control=s)
    with pytest.raises(FileNotFoundError, match="control"):
        run_ame(s, motif_db, tmp_path / "o", control=absent)
    with pytest.raises(FileNotFoundError, match="motif"):
        run_ame(s, absent, tmp_path / "o", control=s)


# --------------------------------------------------------------------------
# Real execution
# --------------------------------------------------------------------------


def _pvalues(tsv):
    return read_ame(tsv).set_index("motif_ID")


@requires_enrichment
@pytest.mark.meme
def test_ame_finds_the_planted_motif_enriched_and_the_decoy_not(
    tmp_path, primary_fasta, control_fasta, enrichment_db
):
    tsv = run_ame(
        primary_fasta,
        enrichment_db,
        tmp_path / "ame",
        control=control_fasta,
        evalue_report_threshold=100,  # so the decoy is reported, not omitted
    )
    assert tsv == tmp_path / "ame" / "ame.tsv"
    df = _pvalues(tsv)
    assert sorted(df.index) == ["DECOY1", "PLANTED1"]

    planted, decoy = df.loc["PLANTED1"], df.loc["DECOY1"]
    assert planted["rank"] == 1
    assert planted["p-value"] < 1e-20
    assert planted["adj_p-value"] < 1e-20
    assert planted["E-value"] < 1e-20
    # adj_p = p x tests, E = adj_p x motifs(2): the three columns are distinct
    # on real output and each is where the reader says it is.
    assert planted["p-value"] < planted["adj_p-value"] < planted["E-value"]
    assert planted["E-value"] == pytest.approx(
        2 * planted["adj_p-value"], rel=0.02, abs=0
    )

    assert decoy["p-value"] > 0.05
    assert decoy["adj_p-value"] > 0.05
    assert decoy["E-value"] > 0.05


@requires_enrichment
@pytest.mark.meme
def test_the_control_set_reaches_ame_and_omitting_it_is_a_different_test(
    tmp_path, primary_fasta, control_fasta, planted_control_fasta, enrichment_db
):
    def adj_p(name, control):
        tsv = run_ame(
            primary_fasta,
            enrichment_db,
            tmp_path / name,
            control=control,
            evalue_report_threshold=100,
            seed=3,
        )
        return _pvalues(tsv).loc["PLANTED1", "adj_p-value"]

    assert adj_p("clean", control_fasta) < 1e-20
    assert adj_p("shuffled", AME_SHUFFLE) < 1e-20
    # A different seed is a different shuffle. Any two seeds can tie (the
    # p-value depends only on how many shuffled sequences happen to match), so
    # this asks for some disagreement among four rather than between two.
    by_seed = set()
    for seed in (1, 2, 3, 4):
        tsv = run_ame(
            primary_fasta,
            enrichment_db,
            tmp_path / f"seed{seed}",
            control=AME_SHUFFLE,
            seed=seed,
        )
        by_seed.add(_pvalues(tsv).loc["PLANTED1", "p-value"])
    assert len(by_seed) > 1, f"four seeds, one p-value: {by_seed}"
    # The motif is in the control too: not enriched.
    assert adj_p("same", planted_control_fasta) > 0.05
    # The trap the required keyword exists for. A motif in 90% of the sequences,
    # and AME's no-control mode does not find it, exits 0 and says nothing.
    assert adj_p("input_order", None) > 0.05


@requires_enrichment
@pytest.mark.meme
def test_the_report_threshold_reaches_ame(
    tmp_path, primary_fasta, control_fasta, enrichment_db
):
    def ids(name, **kw):
        tsv = run_ame(
            primary_fasta, enrichment_db, tmp_path / name, control=control_fasta, **kw
        )
        return sorted(read_ame(tsv)["motif_ID"])

    assert ids("default") == ["DECOY1", "PLANTED1"]  # decoy E = 2, default <= 10
    assert ids("strict", evalue_report_threshold=1.0) == ["PLANTED1"]


@requires_enrichment
@pytest.mark.meme
def test_the_method_reaches_ame(tmp_path, primary_fasta, control_fasta, enrichment_db):
    """Each method writes its own columns; `U` exists only under ranksum."""
    fisher = read_ame(
        run_ame(primary_fasta, enrichment_db, tmp_path / "f", control=control_fasta)
    )
    ranksum = read_ame(
        run_ame(
            primary_fasta,
            enrichment_db,
            tmp_path / "r",
            control=control_fasta,
            method="ranksum",
        )
    )
    assert "U" not in fisher.columns and "TP" in fisher.columns
    assert "U" in ranksum.columns and "TP" not in ranksum.columns


@requires_enrichment
@pytest.mark.meme
def test_an_ame_failure_carries_its_stderr_and_leaves_no_table(
    tmp_path, primary_fasta, control_fasta, enrichment_db
):
    out = tmp_path / "out"
    good = run_ame(primary_fasta, enrichment_db, out, control=control_fasta)
    assert good.is_file()  # positive control, and the stale table for below

    with pytest.raises(AmeError) as e:
        run_ame(primary_fasta, primary_fasta, out, control=control_fasta)
    msg = str(e.value)
    assert "ame exited" in msg and "command:" in msg
    assert "motif" in msg.split("stderr:")[1].lower(), "AME's own complaint is absent"
    assert not good.exists(), "the earlier run's ame.tsv survived a failed run"


def test_a_missing_ame_binary_is_a_setup_error(tmp_path, monkeypatch, motif_db):
    s = write_fasta(tmp_path / "s.fa", [("a", "ACGT" * 20)])
    monkeypatch.setenv("PATH", str(tmp_path))
    monkeypatch.delenv("MEME_BIN", raising=False)
    with pytest.raises(MemeToolNotFound, match="'ame'"):
        run_ame(s, motif_db, tmp_path / "o", control=s)
