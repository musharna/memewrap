"""Locating the MEME suite binaries, and refusing to guess.

The source versions of these wrappers hardcoded

    MEME_BIN = os.path.expanduser("~/miniconda3/envs/meme-suite/bin")

at module scope in four separate scripts. That is not portable to any other
machine, and it is not overridable without editing the file.

Resolution order here: an explicit `meme_bin` argument, then the `MEME_BIN`
environment variable, then `PATH`. A caller who installed the suite through
conda, bioconda, a module system or a container all work without edits.

## The check is EXECUTABILITY, not existence

`verify_tools` in the source was

    {t: os.path.exists(os.path.join(MEME_BIN, t)) for t in (...)}

`os.path.exists` is true for a directory, for a dangling config file, and for a
non-executable text file that happens to be named `streme`. Any of those report
the tool as present and then fail at the point of use -- inside a subprocess,
several minutes into a pipeline, as an OSError rather than a setup error. The
whole point of a verify step is to fail at setup, so the predicate has to be the
one the subprocess will actually apply: a regular file with the execute bit.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

__all__ = [
    "MemeToolNotFound",
    "DEFAULT_TOOLS",
    "meme_bin_dir",
    "find_tool",
    "verify_tools",
    "require_tools",
]

# The three the wrappers in this package drive. `verify_tools` takes any names.
DEFAULT_TOOLS = ("streme", "fimo", "tomtom")


class MemeToolNotFound(RuntimeError):
    """A required MEME-suite executable could not be located."""


def meme_bin_dir(meme_bin: str | os.PathLike | None = None) -> Path | None:
    """The directory to search first, or None to fall back to PATH."""
    if meme_bin is not None:
        return Path(os.path.expanduser(str(meme_bin)))
    env = os.environ.get("MEME_BIN")
    if env:
        return Path(os.path.expanduser(env))
    return None


def _is_executable(p: Path) -> bool:
    # is_file() follows symlinks and excludes directories; os.access adds the bit.
    # Both are needed: a directory named `streme` passes os.access(X_OK) because
    # the execute bit on a directory means "traversable".
    return p.is_file() and os.access(p, os.X_OK)


def find_tool(name: str, meme_bin: str | os.PathLike | None = None) -> Path:
    """Absolute path to a MEME-suite executable.

    Raises MemeToolNotFound with the places searched -- a "not found" that does
    not say where it looked sends the reader to the wrong machine.
    """
    searched: list[str] = []

    d = meme_bin_dir(meme_bin)
    if d is not None:
        candidate = d / name
        if _is_executable(candidate):
            return candidate
        searched.append(str(d))

    on_path = shutil.which(name)
    if on_path:
        return Path(on_path)
    searched.append("PATH")

    raise MemeToolNotFound(
        f"MEME-suite tool {name!r} not found (searched: {', '.join(searched)}). "
        f"Point at the install with the MEME_BIN environment variable or the "
        f"meme_bin= argument, e.g. MEME_BIN=~/miniconda3/envs/meme-suite/bin"
    )


def verify_tools(
    names: tuple[str, ...] = DEFAULT_TOOLS,
    meme_bin: str | os.PathLike | None = None,
) -> dict[str, Path | None]:
    """Map each tool name to its resolved path, or None if unavailable.

    Does not raise -- this is the "tell me what I have" call. Use
    `require_tools` when a missing tool should stop the run.
    """
    out: dict[str, Path | None] = {}
    for n in names:
        try:
            out[n] = find_tool(n, meme_bin)
        except MemeToolNotFound:
            out[n] = None
    return out


def require_tools(
    names: tuple[str, ...] = DEFAULT_TOOLS,
    meme_bin: str | os.PathLike | None = None,
) -> dict[str, Path]:
    """Like `verify_tools` but raises if any are missing, naming ALL of them.

    Reporting every missing tool at once matters: resolving them one error at a
    time means one conda install, one re-run, one new error, repeatedly.
    """
    found = verify_tools(names, meme_bin)
    missing = [n for n, p in found.items() if p is None]
    if missing:
        where = meme_bin_dir(meme_bin)
        raise MemeToolNotFound(
            f"missing MEME-suite tools: {missing} "
            f"(searched {where if where else 'PATH'}{'' if where is None else ' then PATH'})"
        )
    return {n: p for n, p in found.items() if p is not None}
