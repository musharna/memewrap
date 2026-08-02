"""Guards for MEME-database validation.

The case this is really for: a JASPAR download that 404s, or returns an HTML
error page, and lands on disk looking like a file. `verify_meme_db` is what
stands between that and a scan that finds nothing and reads as a result.
"""

from __future__ import annotations

import pytest

from memewrap import verify_meme_db

VALID = (
    "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
    "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"
    "MOTIF M1 first\n"
    "letter-probability matrix: alength= 4 w= 2 nsites= 20 E= 1e-10\n"
    "0.9 0.1 0.0 0.0\n0.0 0.0 0.9 0.1\n\n"
    "MOTIF M2 second\n"
    "letter-probability matrix: alength= 4 w= 2 nsites= 20 E= 1e-10\n"
    "0.1 0.9 0.0 0.0\n0.0 0.0 0.1 0.9\n"
)


def test_a_valid_database_is_counted(tmp_path):
    p = tmp_path / "db.meme"
    p.write_text(VALID)
    db = verify_meme_db(p)
    assert db.ok and db.exists
    assert db.n_motifs == 2
    assert db.meme_version == "4"


def test_a_missing_file_is_not_ok(tmp_path):
    db = verify_meme_db(tmp_path / "absent.meme")
    assert not db.exists and not db.ok and db.n_motifs == 0
    with pytest.raises(FileNotFoundError):
        db.require()


def test_an_HTML_error_page_is_rejected(tmp_path):
    """The shape of a failed JASPAR download."""
    p = tmp_path / "db.meme"
    p.write_text("<html><head><title>404 Not Found</title></head></html>")
    db = verify_meme_db(p)
    assert not db.ok
    with pytest.raises(ValueError, match="not a MEME-format"):
        db.require()


def test_a_file_with_MOTIF_lines_but_no_MEME_HEADER_is_rejected(tmp_path):
    """The source counted `^MOTIF\\s` alone, so any file mentioning motifs at the
    start of a line -- a TOMTOM report, a notes file -- validated as a database.
    """
    p = tmp_path / "notes.txt"
    p.write_text("MOTIF M1 something\nMOTIF M2 another\n")
    db = verify_meme_db(p)
    assert db.n_motifs == 2, "precondition: the source's count is satisfied"
    assert not db.ok, "a file with no MEME version header is not a database"


def test_a_TRUNCATED_database_passes_ok_but_fails_expect_min(tmp_path):
    """`ok` cannot express 'I expected JASPAR plants CORE'. A truncated download
    parses cleanly, finds fewer hits, and looks like a biological result."""
    p = tmp_path / "db.meme"
    p.write_text(VALID)
    db = verify_meme_db(p)
    assert db.ok
    db.require(expect_min=2)  # positive control: the real count is accepted
    with pytest.raises(ValueError, match="expected at least 500"):
        db.require(expect_min=500)


def test_require_returns_self_so_it_can_be_chained(tmp_path):
    p = tmp_path / "db.meme"
    p.write_text(VALID)
    assert verify_meme_db(p).require().n_motifs == 2
