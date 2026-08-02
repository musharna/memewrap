"""Validating a MEME-format motif database before a pipeline depends on it.

The source's check was:

    n = len(re.findall(r"^MOTIF\\s", fh.read(), flags=re.M))
    return {"n_motifs": n, "ok": n >= 1}

That is the right idea -- it catches the case this is actually guarding, where
a JASPAR download 404s or returns an HTML error page and lands on disk as a
plausible-looking file. One motif is enough to prove it parsed.

Two things are added here. `ok` alone cannot express "I expected the JASPAR
plants CORE set and got a file with one motif in it", so `expect_min` makes the
count assertable. And the MEME version header is checked, because a file can
contain `MOTIF` lines and still not be a MEME database -- a TOMTOM output TSV
mentioning motif names would count several.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

__all__ = ["MotifDb", "verify_meme_db"]

_MOTIF_RE = re.compile(r"^MOTIF\s+\S", re.M)
_VERSION_RE = re.compile(r"^MEME version\s+(\S+)", re.M)


@dataclass(frozen=True)
class MotifDb:
    path: Path
    n_motifs: int
    meme_version: str | None
    exists: bool

    @property
    def ok(self) -> bool:
        """Parsed as a MEME database with at least one motif."""
        return self.exists and self.n_motifs >= 1 and self.meme_version is not None

    def require(self, expect_min: int = 1) -> MotifDb:
        """Raise unless the database holds at least `expect_min` motifs.

        `expect_min` is the difference between "a file that parses" and "the
        database I meant". A truncated download of JASPAR plants CORE parses
        fine and yields three motifs; a scan against it returns almost nothing
        and looks like a biological result.
        """
        if not self.exists:
            raise FileNotFoundError(f"motif database does not exist: {self.path}")
        if self.meme_version is None:
            raise ValueError(
                f"{self.path} has no 'MEME version' header -- it is not a MEME-format "
                f"database. A failed download often lands here as an HTML error page."
            )
        if self.n_motifs < expect_min:
            raise ValueError(
                f"{self.path} has {self.n_motifs} motifs, expected at least "
                f"{expect_min}. A truncated download parses cleanly and simply "
                f"finds fewer hits, which reads as a biological result."
            )
        return self


def verify_meme_db(path: str | os.PathLike) -> MotifDb:
    """Inspect a MEME-format motif database. Does not raise; see `MotifDb.require`."""
    path = Path(path)
    if not path.is_file():
        return MotifDb(path=path, n_motifs=0, meme_version=None, exists=False)
    text = path.read_text(errors="replace")
    version = _VERSION_RE.search(text)
    return MotifDb(
        path=path,
        n_motifs=len(_MOTIF_RE.findall(text)),
        meme_version=version.group(1) if version else None,
        exists=True,
    )
