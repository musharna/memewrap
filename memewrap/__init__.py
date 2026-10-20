"""memewrap -- Python wrappers for the MEME-suite command-line tools.

Extracted from duplicated wrapper code in two unrelated plant-genomics projects,
where the same STREME, TOMTOM and FIMO invocations had been copied across four
scripts with a hardcoded conda path at module scope in each.

Scope, checked on PyPI and GitHub (2026-10-19; the README has the table):

- **memewrap runs the MEME binaries, so results are MEME's.** STREME, TOMTOM and
  FIMO now have Python reimplementations (`pystreme`; `memelite` for TOMTOM and
  FIMO; `lightmotif` for FIMO-style scanning) that are faster but are
  independent code whose results can differ from MEME's. `gimmemotifs` 0.18.4
  wraps `dreme.py`, the tool STREME replaced, and has no TOMTOM wrapper.
- **For FIMO, prefer `pymemesuite` when you can install it.** It is Cython
  bindings to FIMO's code. `memewrap.fimo` drives the CLI, for conda
  environments without a compiler, and parallelises over sequence chunks.
- **SEA and AME** (0.3.0) are wrapped for the step after discovery: given known
  motifs, which are enriched in these sequences relative to those. No other
  Python package found covers them.

Each wrapper fixes something measured in the source it came from -- see the
individual module docstrings.
"""

from .ame import AME_SHUFFLE, AmeError, read_ame, run_ame
from .fimo import FimoError, build_feature_matrix, run_fimo, run_fimo_parallel
from .motifdb import MotifDb, verify_meme_db
from .sea import SeaError, read_sea, run_sea
from .streme import StremeError, run_streme
from .tomtom import TomtomError, match_count, read_tomtom, run_tomtom
from .tools import (
    DEFAULT_TOOLS,
    ENRICHMENT_TOOLS,
    MemeToolNotFound,
    find_tool,
    require_tools,
    verify_tools,
)

__version__ = "0.3.0"

__all__ = [
    "__version__",
    # tools
    "MemeToolNotFound",
    "DEFAULT_TOOLS",
    "ENRICHMENT_TOOLS",
    "find_tool",
    "verify_tools",
    "require_tools",
    # streme
    "StremeError",
    "run_streme",
    # tomtom
    "TomtomError",
    "run_tomtom",
    "read_tomtom",
    "match_count",
    # fimo
    "FimoError",
    "run_fimo",
    "run_fimo_parallel",
    "build_feature_matrix",
    # sea
    "SeaError",
    "run_sea",
    "read_sea",
    # ame
    "AmeError",
    "AME_SHUFFLE",
    "run_ame",
    "read_ame",
    # motifdb
    "MotifDb",
    "verify_meme_db",
]
