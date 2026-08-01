"""memewrap -- Python wrappers for the MEME-suite command-line tools.

Extracted from duplicated wrapper code in two active research repos
(`arf_promoter_analysis` and `phelipanche-fm`), where the same STREME, TOMTOM
and FIMO invocations had been copied across four scripts with a hardcoded conda
path at module scope in each.

Scope, decided against a live PyPI check rather than from memory:

- **STREME and TOMTOM are the gap.** `gimmemotifs` 0.18.4 wraps 22 discovery
  tools -- including `dreme.py`, the tool STREME deprecated and replaced -- but
  ships no `streme.py` and no TOMTOM wrapper. `pymemesuite` 0.1.0a4 does not
  mention either.
- **FIMO is only partly a gap.** `pymemesuite` provides real Cython bindings for
  it; prefer those when you can install them. `memewrap.fimo` drives the CLI, for
  conda environments without a compiler, and parallelises over sequence chunks.

Each wrapper fixes something measured in the source it came from -- see the
individual module docstrings.
"""

from .fimo import FimoError, build_feature_matrix, run_fimo, run_fimo_parallel
from .motifdb import MotifDb, verify_meme_db
from .streme import StremeError, run_streme
from .tomtom import TomtomError, match_count, read_tomtom, run_tomtom
from .tools import (
    DEFAULT_TOOLS,
    MemeToolNotFound,
    find_tool,
    require_tools,
    verify_tools,
)

__version__ = "0.1.0"

__all__ = [
    "__version__",
    # tools
    "MemeToolNotFound",
    "DEFAULT_TOOLS",
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
    # motifdb
    "MotifDb",
    "verify_meme_db",
]
