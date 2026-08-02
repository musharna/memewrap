"""Guards for the STREME wrapper.

The headline guard is `test_nmotifs_and_thresh_cannot_both_be_given`. The source
passed both at one call site; STREME's own help states that `--nmotifs`
overrides `--thresh` when positive, so that call was requesting a significance
filter and receiving a fixed count. Nothing reports this -- the run succeeds and
the output looks the same either way.
"""

from __future__ import annotations

from contextlib import suppress

import pytest

from conftest import PLANTED_MOTIF, make_corpus, requires_meme, write_fasta
from memewrap import MemeToolNotFound, StremeError, run_streme


# --------------------------------------------------------------------------
# Argument validation -- no binary needed
# --------------------------------------------------------------------------


def test_nmotifs_and_thresh_cannot_both_be_given(tmp_path):
    p = write_fasta(tmp_path / "p.fa", [("a", "ACGT" * 20)])
    c = write_fasta(tmp_path / "c.fa", [("b", "TGCA" * 20)])
    with pytest.raises(ValueError, match="mutually exclusive"):
        run_streme(p, c, tmp_path / "out", nmotifs=5, thresh=0.05)


def test_either_one_alone_is_accepted(tmp_path):
    """Positive control for the test above: a function that rejected BOTH
    arguments unconditionally would satisfy it."""
    p = write_fasta(tmp_path / "p.fa", [("a", "ACGT" * 20)])
    c = write_fasta(tmp_path / "c.fa", [("b", "TGCA" * 20)])
    for kwargs in ({"nmotifs": 5}, {"thresh": 0.05}):
        # Getting as far as launching STREME -- which then objects to this
        # two-sequence fixture -- is the outcome being checked for: it means
        # validation let the call through. A ValueError is NOT suppressed, so a
        # regression that refused either argument alone fails the test here.
        with suppress(StremeError, MemeToolNotFound):
            run_streme(p, c, tmp_path / "out", **kwargs)


def test_a_missing_primary_fasta_is_named(tmp_path):
    c = write_fasta(tmp_path / "c.fa", [("b", "TGCA" * 20)])
    with pytest.raises(FileNotFoundError, match="primary"):
        run_streme(tmp_path / "absent.fa", c, tmp_path / "out")


def test_a_missing_control_fasta_is_named(tmp_path):
    p = write_fasta(tmp_path / "p.fa", [("a", "ACGT" * 20)])
    with pytest.raises(FileNotFoundError, match="control"):
        run_streme(p, tmp_path / "absent.fa", tmp_path / "out")


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"minw": 2}, "minw"),
        ({"minw": 10, "maxw": 8}, "maxw"),
        ({"nmotifs": 0}, "nmotifs"),
        ({"nmotifs": -1}, "nmotifs"),
    ],
)
def test_invalid_widths_and_counts_are_refused_before_launching(
    tmp_path, kwargs, match
):
    p = write_fasta(tmp_path / "p.fa", [("a", "ACGT" * 20)])
    c = write_fasta(tmp_path / "c.fa", [("b", "TGCA" * 20)])
    with pytest.raises(ValueError, match=match):
        run_streme(p, c, tmp_path / "out", **kwargs)


# --------------------------------------------------------------------------
# Real execution
# --------------------------------------------------------------------------


@requires_meme
@pytest.mark.meme
def test_streme_recovers_a_planted_motif(tmp_path, primary_fasta, control_fasta):
    """The end-to-end positive control for this whole package.

    Everything else asserts that a wrapper refuses bad input or parses an
    output. This asserts that the flags reach STREME and that STREME does its
    job on data whose answer is known: a motif planted in 90% of the primary
    sequences and none of the controls comes back.
    """
    out = run_streme(
        primary_fasta, control_fasta, tmp_path / "streme", nmotifs=1, minw=6, maxw=10
    )
    assert out.is_file() and out.name == "streme.txt"
    motifs = [ln for ln in out.read_text().splitlines() if ln.startswith("MOTIF")]
    assert motifs, "STREME reported no motif at all"
    assert PLANTED_MOTIF in motifs[0], (
        f"the planted motif {PLANTED_MOTIF} was not recovered; got {motifs[0]!r}"
    )


@requires_meme
@pytest.mark.meme
def test_nmotifs_bounds_the_number_reported(tmp_path, primary_fasta, control_fasta):
    """Confirms the flag is actually reaching the binary.

    Without this, `test_streme_recovers_a_planted_motif` would pass even if
    every keyword were dropped on the floor -- STREME's defaults would still
    find the motif.
    """
    out = run_streme(
        primary_fasta, control_fasta, tmp_path / "s1", nmotifs=1, minw=6, maxw=10
    )
    n = sum(ln.startswith("MOTIF") for ln in out.read_text().splitlines())
    assert n == 1, f"asked for 1 motif, got {n}"


@requires_meme
@pytest.mark.meme
def test_a_streme_failure_carries_its_stderr(tmp_path):
    """The source used check=True, whose message is the command and exit code
    only -- STREME's actual complaint never reaches the reader."""
    bad = write_fasta(tmp_path / "bad.fa", [("a", "not-a-sequence!!")])
    ctrl = write_fasta(tmp_path / "c.fa", make_corpus(5, 50, seed=9))
    with pytest.raises(StremeError) as e:
        run_streme(bad, ctrl, tmp_path / "out", nmotifs=1, minw=6, maxw=8)
    msg = str(e.value)
    assert "streme exited" in msg
    assert "stderr:" in msg and "command:" in msg
