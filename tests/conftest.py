"""Fixtures: a real planted-motif corpus, and the MEME-availability gate.

## Why the sequences are generated rather than bundled

A wrapper for a discovery tool has to be tested on input where the answer is
KNOWN, otherwise "STREME returned something" is the only available assertion and
it passes whether or not the flags reached the binary. Here a fixed motif is
planted into a defined fraction of the primary sequences and into none of the
controls, so a discovery run has a right answer to be checked against.

## Why the MEME tests are not skipped by default

These wrappers exist to drive external binaries. A suite that mocks
`subprocess.run` proves the mock agrees with the wrapper, which is the one
relationship that cannot break in production. So the real-execution tests run
whenever the suite is installed, and skip -- individually, with a reason -- only
when it is not.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from memewrap.tools import ENRICHMENT_TOOLS, verify_tools

# A motif with no strong self-similarity, so a TOMTOM self-match is meaningful.
PLANTED_MOTIF = "TGTCTCTC"

# The availability gates are fixtures, not `skipif(verify_tools() ...)` markers.
# A marker condition runs package code while conftest is being IMPORTED; any
# exception there (mutmut's forced-fail probe, or a mutant in verify_tools /
# find_tool) becomes a conftest ImportError, pytest exits 4 (usage error), and
# mutmut aborts the whole run instead of scoring one test as failed. Inside a
# fixture the same exception fails the test that needed it.


@pytest.fixture(scope="session")
def _meme_suite() -> None:
    missing = [name for name, path in verify_tools().items() if path is None]
    if missing:
        pytest.skip(
            f"MEME suite not installed (missing: {missing}). "
            f"conda create -n meme-suite -c bioconda -c conda-forge meme, "
            f"then set MEME_BIN to its bin/ directory."
        )


# Separate gate for the SEA/AME tests: `sea` only exists from MEME 5.4.0, so an
# older install should skip those and still run everything above.
@pytest.fixture(scope="session")
def _meme_enrichment() -> None:
    missing = [
        name for name, path in verify_tools(ENRICHMENT_TOOLS).items() if path is None
    ]
    if missing:
        pytest.skip(f"MEME suite enrichment tools not installed (missing: {missing})")


requires_meme = pytest.mark.usefixtures("_meme_suite")
requires_enrichment = pytest.mark.usefixtures("_meme_enrichment")

# Planted nowhere. The negative half of every enrichment assertion: a wrapper
# that reported everything as significant would pass a planted-motif-only test.
DECOY_MOTIF = "GATTACAG"


def _random_dna(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")
    return path


def make_corpus(
    n: int, length: int, seed: int, motif: str | None = None, frac: float = 1.0
) -> list[tuple[str, str]]:
    """`n` sequences; if `motif` is given, plant it in `frac` of them."""
    rng = random.Random(seed)
    out: list[tuple[str, str]] = []
    for i in range(n):
        seq = _random_dna(length, rng)
        if motif and rng.random() < frac:
            pos = rng.randint(0, length - len(motif))
            seq = seq[:pos] + motif + seq[pos + len(motif) :]
        out.append((f"seq{i}", seq))
    return out


@pytest.fixture(scope="session")
def primary_fasta(tmp_path_factory) -> Path:
    """200 sequences, 90% carrying PLANTED_MOTIF."""
    d = tmp_path_factory.mktemp("corpus")
    return write_fasta(
        d / "primary.fa", make_corpus(200, 150, seed=1, motif=PLANTED_MOTIF, frac=0.9)
    )


@pytest.fixture(scope="session")
def control_fasta(tmp_path_factory) -> Path:
    """200 sequences from the same composition, motif absent."""
    d = tmp_path_factory.mktemp("corpus")
    return write_fasta(d / "control.fa", make_corpus(200, 150, seed=2))


@pytest.fixture(scope="session")
def small_fasta(tmp_path_factory) -> Path:
    """12 sequences -- enough to exercise chunking, fast enough to scan."""
    d = tmp_path_factory.mktemp("small")
    return write_fasta(
        d / "small.fa", make_corpus(12, 120, seed=3, motif=PLANTED_MOTIF, frac=1.0)
    )


@pytest.fixture(scope="session")
def motif_db(tmp_path_factory) -> Path:
    """A minimal hand-written MEME database containing PLANTED_MOTIF.

    Written literally rather than produced by STREME so the FIMO and TOMTOM
    tests do not depend on the STREME wrapper working -- a shared upstream
    fixture would let one broken wrapper mark three suites green or red
    together.
    """
    d = tmp_path_factory.mktemp("db")
    p = d / "planted.meme"
    rows = []
    for base in PLANTED_MOTIF:
        rows.append(" ".join("0.997" if b == base else "0.001" for b in "ACGT"))
    p.write_text(
        "MEME version 4\n\n"
        "ALPHABET= ACGT\n\n"
        "strands: + -\n\n"
        "Background letter frequencies\n"
        "A 0.25 C 0.25 G 0.25 T 0.25\n\n"
        "MOTIF PLANTED1 planted\n"
        f"letter-probability matrix: alength= 4 w= {len(PLANTED_MOTIF)} "
        "nsites= 100 E= 1e-30\n" + "\n".join(rows) + "\n"
    )
    return p


def _motif_block(name: str, alt: str, consensus: str) -> str:
    rows = [
        " ".join("0.997" if b == base else "0.001" for b in "ACGT")
        for base in consensus
    ]
    return (
        f"MOTIF {name} {alt}\n"
        f"letter-probability matrix: alength= 4 w= {len(consensus)} "
        "nsites= 100 E= 1e-30\n" + "\n".join(rows) + "\n\n"
    )


@pytest.fixture(scope="session")
def enrichment_db(tmp_path_factory) -> Path:
    """PLANTED1 (in 90% of `primary_fasta`) and DECOY1 (planted nowhere)."""
    p = tmp_path_factory.mktemp("db") / "planted_and_decoy.meme"
    p.write_text(
        "MEME version 4\n\n"
        "ALPHABET= ACGT\n\n"
        "strands: + -\n\n"
        "Background letter frequencies\n"
        "A 0.25 C 0.25 G 0.25 T 0.25\n\n"
        + _motif_block("PLANTED1", "planted", PLANTED_MOTIF)
        + _motif_block("DECOY1", "decoy", DECOY_MOTIF)
    )
    return p


@pytest.fixture(scope="session")
def planted_control_fasta(tmp_path_factory) -> Path:
    """A control set that ALSO carries PLANTED_MOTIF in 90% of sequences.

    Against this control the planted motif is not enriched. It is what makes the
    control argument observable: if a wrapper dropped it, SEA/AME would fall
    back to a shuffled or absent control and report the motif as enriched.
    """
    d = tmp_path_factory.mktemp("corpus")
    return write_fasta(
        d / "control_planted.fa",
        make_corpus(200, 150, seed=5, motif=PLANTED_MOTIF, frac=0.9),
    )
