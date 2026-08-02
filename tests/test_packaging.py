"""The version in the metadata must equal the version in the module.

seqbench shipped 0.2.0 while reporting 0.1.0 because its only version test
asserted the string was non-empty. Asserting equality against the INSTALLED
distribution metadata is the check that would have caught it: pyproject reads
the version dynamically from `memewrap.__version__`, so the two agreeing is
exactly what proves the dynamic read is wired up.
"""

from __future__ import annotations

import re
from importlib.metadata import PackageNotFoundError, version

import pytest

import memewrap


def test_the_module_declares_a_pep440_version():
    assert re.fullmatch(r"\d+\.\d+\.\d+([abc]\d+|\.dev\d+)?", memewrap.__version__), (
        f"not a PEP 440 release version: {memewrap.__version__!r}"
    )


def test_the_installed_metadata_agrees_with_the_module():
    try:
        installed = version("memewrap")
    except PackageNotFoundError:  # pragma: no cover - not installed
        pytest.skip("memewrap is not installed; run `pip install -e .`")
    assert installed == memewrap.__version__, (
        f"distribution metadata says {installed}, module says "
        f"{memewrap.__version__} -- the dynamic version read is broken, and a "
        f"release would ship metadata that disagrees with the code"
    )


def test_the_public_api_is_importable_from_the_top_level():
    """__all__ that names something absent fails only at `from memewrap import X`
    in a user's script, not in any test that imports submodules directly."""
    missing = [n for n in memewrap.__all__ if not hasattr(memewrap, n)]
    assert not missing, f"__all__ names absent attributes: {missing}"
