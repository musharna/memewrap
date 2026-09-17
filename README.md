# memewrap

Python wrappers for the MEME-suite command-line tools: **STREME**, **TOMTOM** and **FIMO**.

Extracted from wrapper code that had been copied across four scripts in two active
research repos, each with `MEME_BIN = os.path.expanduser("~/miniconda3/envs/meme-suite/bin")`
hardcoded at module scope.

## Why this exists (and where it doesn't)

Checked against PyPI rather than assumed:

|                                                                | STREME                                                                                                         | TOMTOM                                | FIMO                       |
| -------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- | ------------------------------------- | -------------------------- |
| [`gimmemotifs`](https://pypi.org/project/gimmemotifs/) 0.18.4  | ✗ — wraps 22 tools including **`dreme.py`**, the tool STREME deprecated and replaced, but ships no `streme.py` | ✗ — has its own motif-similarity code | ✗                          |
| [`pymemesuite`](https://pypi.org/project/pymemesuite/) 0.1.0a4 | ✗                                                                                                              | ✗                                     | ✓ **real Cython bindings** |
| `memewrap`                                                     | ✓                                                                                                              | ✓                                     | ✓ (CLI, chunked)           |

**STREME and TOMTOM are the actual gap.** The ecosystem wrapped the previous
generation of discriminative motif discovery and never followed it forward.

**FIMO is only partly a gap — and if you can install `pymemesuite`, prefer it.**
It binds FIMO natively: no subprocess, no temp files. `memewrap.fimo` drives the
CLI instead, which is what a bioconda MEME environment gives you without a
compiler, and parallelises across sequence chunks.

## Install

```bash
pip install memewrap
```

`pip` gives you the wrappers and **none of the binaries** — the MEME suite is a C
toolchain, not a Python package:

```bash
conda create -n meme-suite -c bioconda -c conda-forge meme
export MEME_BIN=~/miniconda3/envs/meme-suite/bin   # or just put them on PATH
```

Resolution order is `meme_bin=` argument → `$MEME_BIN` → `PATH`. Every wrapper
raises `MemeToolNotFound` naming where it looked, so a missing suite fails at
setup rather than several minutes into a scan.

## Use

```python
from memewrap import run_streme, match_count, run_fimo_parallel, build_feature_matrix

# Discriminative discovery: what is enriched in `primary` relative to `control`?
motifs = run_streme(
    "primary.fa", "control.fa", "out/streme", nmotifs=5, minw=6, maxw=12
)

# How many of my motifs match a reference set? (the shape a permutation test needs)
n = match_count(motifs, "jaspar_plants.meme", q_thresh=0.05)

# Scan, then build a (gene x motif) design matrix of max FIMO scores.
tsv = run_fimo_parallel("promoters.fa", "jaspar_plants.meme", "out/fimo", n_chunks=8)
X = build_feature_matrix(tsv, gene_ids, motif_ids)

# Need q-values? Scan in one process. FIMO computes q-values from the whole
# run's test count, so the chunked scan above CANNOT have them -- its q-value
# column is empty, and run_fimo_parallel warns about that on every call.
tsv = run_fimo("promoters.fa", "jaspar_plants.meme", "out/fimo.tsv", thresh="1e-4")
```

## What each wrapper fixes

These are measured against the source they came from, not stylistic preferences.

**`run_streme` refuses `nmotifs` and `thresh` together.** STREME's own help says
`--nmotifs` _overrides_ `--thresh` when positive. One source call site passed
both, so it read as "significant motifs, up to N" and meant "exactly N motifs,
significant or not" — STREME will emit its Nth motif at p=0.9. Nothing in the
output says which rule applied.

_Related, and worth stating because it looked like a bug and wasn't:_ the two
source call sites appeared to diverge, one passing `--order 2 --thresh 0.05` and
one passing neither. Against STREME 5.5.9 those are both the defaults, so the
divergence was cosmetic.

**`run_fimo` computes q-values; `run_fimo_parallel` says out loud that it can't.**
Until 0.2.0 both paths hardcoded `fimo --text`, which streams hits but skips the
q-value computation — the column is emitted with every cell empty, which
`pd.read_csv` reads as NaN and nothing downstream flags. Chunking is the reason:
FIMO's q-values depend on the whole run's test count, so no per-chunk scan can
produce them. `run_fimo` now runs FIMO normally (`--oc`) and copies its
`fimo.tsv`, q-values included; `text=True` restores the old streaming mode.
`run_fimo_parallel` keeps `--text` (that is what makes it parallel) and warns.

**`run_fimo_parallel` no longer loses the header.** The source kept the header
from chunk 0 and stripped line 0 of every later chunk. `fimo --text` writes a
header only when it writes output, so an empty chunk 0 — routine, since chunks
are dealt round-robin — meant chunk 1's header got stripped as a duplicate.
Downstream, `pd.read_csv` promotes the first data row to column names: one hit
silently lost, every column mislabelled, no exception. The header is now taken
from the first chunk that has one.

**`match_count` applies its threshold where it can be observed.** The source
passed `q_thresh` as TOMTOM's `-thresh` _and_ re-filtered at the same value, so
the second filter could never drop a row. TOMTOM now runs permissively and the
filtering happens here — which also leaves the written `tomtom.tsv` complete, so
a stricter threshold doesn't require re-running every comparison.

**`verify_tools` checks executability, not existence.** `os.path.exists` is true
for a directory named `streme` and for a non-executable HTML error page saved
under that name. Both reported the tool present and then failed inside a
subprocess.

**`verify_meme_db` checks the MEME version header, and takes `expect_min`.**
Counting `^MOTIF` alone validates any file that mentions motifs at line start.
And `ok` cannot express "I expected JASPAR plants CORE and got nine motifs" — a
truncated download parses cleanly, finds fewer hits, and reads as a result.

## Tests

```bash
pip install -e '.[dev]'
pytest
```

57 tests. The MEME-dependent ones **run by default** and skip individually, with
a reason, only when the binaries are absent — a wrapper verified against a mocked
`subprocess.run` proves the mock matches the wrapper, which is the one
relationship that cannot break in production.

The end-to-end control plants a known motif in 90% of the primary sequences and
none of the controls, and asserts STREME recovers it; FIMO is then really run,
really chunked and really merged over sequences whose hits are known.

`tests/test_fimo.py` reproduces the source's merge loop verbatim and asserts it
**fails** — the "never trust a test you haven't seen fail" control made
permanent, so a regression turns two tests red rather than one.

A 15-mutant pass kills 15/15, including every source behaviour listed above. Two
guards were deleted rather than kept after mutation showed them unreachable: an
empty-chunk check already guaranteed by a `min()` bound, and a `dropna` already
handled by pandas' comment and blank-line defaults.

## Licence

MIT.
