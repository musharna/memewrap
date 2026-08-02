"""Guards for tool discovery.

The predicate matters more than it looks. The source used `os.path.exists`,
which is true for things that cannot be executed, so `verify_tools` reported
success and the failure surfaced later inside a subprocess. Two of the tests
below are exactly those cases -- a directory named `streme`, and a
non-executable file named `streme` -- and both pass under the source predicate.
"""

from __future__ import annotations

import os
import stat

import pytest

from memewrap.tools import (
    MemeToolNotFound,
    find_tool,
    require_tools,
    verify_tools,
)


def _make_exe(path, body="#!/bin/sh\nexit 0\n"):
    path.write_text(body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def test_an_explicit_directory_is_searched_first(tmp_path):
    exe = _make_exe(tmp_path / "streme")
    assert find_tool("streme", meme_bin=tmp_path) == exe


def test_the_MEME_BIN_environment_variable_is_honoured(tmp_path, monkeypatch):
    exe = _make_exe(tmp_path / "fimo")
    monkeypatch.setenv("MEME_BIN", str(tmp_path))
    assert find_tool("fimo") == exe


def test_an_explicit_argument_beats_the_environment(tmp_path, monkeypatch):
    """Otherwise a caller who passes meme_bin= is silently overridden by a stale
    variable in their shell, and gets a different binary than they asked for."""
    env_dir = tmp_path / "env"
    arg_dir = tmp_path / "arg"
    env_dir.mkdir()
    arg_dir.mkdir()
    _make_exe(env_dir / "tomtom")
    wanted = _make_exe(arg_dir / "tomtom")
    monkeypatch.setenv("MEME_BIN", str(env_dir))
    assert find_tool("tomtom", meme_bin=arg_dir) == wanted


def test_PATH_is_the_fallback(tmp_path, monkeypatch):
    _make_exe(tmp_path / "streme")
    monkeypatch.delenv("MEME_BIN", raising=False)
    monkeypatch.setenv("PATH", str(tmp_path))
    assert find_tool("streme").name == "streme"


def test_a_DIRECTORY_named_like_the_tool_is_not_accepted(tmp_path, monkeypatch):
    """Passes under the source's os.path.exists check.

    It also passes a naive os.access(X_OK) check on its own, because the execute
    bit on a directory means 'traversable' -- which is why the implementation
    requires is_file() as well.
    """
    (tmp_path / "streme").mkdir()
    monkeypatch.delenv("MEME_BIN", raising=False)
    monkeypatch.setenv("PATH", "/nonexistent-for-this-test")
    assert os.path.exists(tmp_path / "streme"), "precondition: the source check passes"
    with pytest.raises(MemeToolNotFound):
        find_tool("streme", meme_bin=tmp_path)


def test_a_NON_EXECUTABLE_file_is_not_accepted(tmp_path, monkeypatch):
    """The shape of a failed download: an HTML error page saved as `streme`."""
    p = tmp_path / "streme"
    p.write_text("<html>404</html>")
    p.chmod(0o644)
    monkeypatch.delenv("MEME_BIN", raising=False)
    monkeypatch.setenv("PATH", "/nonexistent-for-this-test")
    assert os.path.exists(p), "precondition: the source check passes"
    with pytest.raises(MemeToolNotFound):
        find_tool("streme", meme_bin=tmp_path)


def test_the_error_names_where_it_looked(tmp_path, monkeypatch):
    """A 'not found' with no search path sends the reader to the wrong machine."""
    monkeypatch.delenv("MEME_BIN", raising=False)
    monkeypatch.setenv("PATH", "/nonexistent-for-this-test")
    with pytest.raises(MemeToolNotFound) as e:
        find_tool("streme", meme_bin=tmp_path)
    assert str(tmp_path) in str(e.value)
    assert "PATH" in str(e.value)
    assert "MEME_BIN" in str(e.value), "the error should say how to fix it"


def test_verify_tools_reports_None_rather_than_raising(tmp_path, monkeypatch):
    _make_exe(tmp_path / "fimo")
    monkeypatch.delenv("MEME_BIN", raising=False)
    monkeypatch.setenv("PATH", "/nonexistent-for-this-test")
    got = verify_tools(("streme", "fimo", "tomtom"), meme_bin=tmp_path)
    assert got["fimo"] is not None
    assert got["streme"] is None and got["tomtom"] is None


def test_require_tools_names_EVERY_missing_tool(tmp_path, monkeypatch):
    """Reporting one at a time means one conda install per re-run."""
    _make_exe(tmp_path / "fimo")
    monkeypatch.delenv("MEME_BIN", raising=False)
    monkeypatch.setenv("PATH", "/nonexistent-for-this-test")
    with pytest.raises(MemeToolNotFound) as e:
        require_tools(("streme", "fimo", "tomtom"), meme_bin=tmp_path)
    msg = str(e.value)
    assert "streme" in msg and "tomtom" in msg
    assert "'fimo'" not in msg, "fimo was present; it should not be listed missing"


def test_require_tools_returns_the_paths_when_all_present(tmp_path):
    """Positive control: the raising tests above are all satisfied by a function
    that raises unconditionally."""
    for t in ("streme", "fimo", "tomtom"):
        _make_exe(tmp_path / t)
    got = require_tools(meme_bin=tmp_path)
    assert sorted(got) == ["fimo", "streme", "tomtom"]
    assert all(p.is_file() for p in got.values())
