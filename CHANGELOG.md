# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions: SemVer.

## [0.3.0] - 2026-09-18

### Added

- **`run_sea` / `read_sea`** (`memewrap.sea`): SEA enrichment of known motifs in primary
  vs control sequences. Returns the `sea.tsv` path; `read_sea` parses it with declared
  dtypes (`SEA_DTYPES`: PVALUE, EVALUE, QVALUE and their logs as float64, counts as
  int64, IDs as strings), drops the trailing `#` provenance block, and returns a typed
  empty frame for a header-only table. `thresh` + `threshold_on="evalue"|"qvalue"|"pvalue"`
  replace SEA's two combinable flags; a q/p threshold above 1 is refused before launch.
- **`run_ame` / `read_ame`** (`memewrap.ame`): the same for AME (`p-value`,
  `adj_p-value`, `E-value`; AME reports no q-value). `control` is a required keyword --
  a FASTA, `AME_SHUFFLE`, or an explicit `None` -- because `ame` without `--control`
  ranks by input order instead of shuffling and exits 0 (measured on 5.5.9: p = 0.07
  for a motif planted in 90% of sequences, against 4.4e-97 with a control file).
- Both default to `--text`, writing stdout to `<outdir>/sea.tsv` / `ame.tsv` only after
  a zero exit and removing any earlier table first. `text=False` runs `--oc`. Observed
  on bioconda MEME 5.5.9 build `pl5321he99cc7f_1`, locally and on a fresh CI install: `sea --oc` and `ame --oc`
  exit 1 with "Template does not contain data section" after writing a header-only TSV.
- `ENRICHMENT_TOOLS = ("sea", "ame")`. `DEFAULT_TOOLS` is unchanged, so `require_tools()`
  behaves as before on MEME installs older than 5.4.0, which have no `sea`.
- Real-execution tests against a two-motif database (planted + decoy): the planted motif
  is significant and the decoy is not; a control set carrying the motif makes `control`
  observable; thresholds, seed and method are each shown to reach the binary. Run against
  MEME 5.5.9 and 5.4.1. 23 mutants run, all killed after one test fix (see README).
- `packaging/bioconda/meta.yaml`: DRAFT recipe, not submitted, not linted, not built.

### Changed

- CI's tool-visibility step now also requires `sea` and `ame`.

## [0.2.0] - 2026-09-16

### Fixed

- **`run_fimo` now produces q-values.** Both FIMO paths hardcoded `--text`, which
  streams hits but skips FIMO's q-value computation; the `q-value` column was emitted
  with every cell empty, and `pd.read_csv` turned that into NaN with no signal. Verified
  against FIMO 5.5.9: `--text` → empty column; `--oc` → populated. Found by an
  independent review panel on 2026-09-15.

### Changed

- `run_fimo(...)` runs FIMO in normal `--oc` mode by default and copies its `fimo.tsv`
  (trailing `#` provenance comments stripped) to `out_tsv`. New keyword `text=True`
  restores the previous streaming behaviour; new keyword `max_stored_scores`
  (default 100 000, FIMO's own default) is passed through because normal mode drops the
  weakest hits beyond it.
- `run_fimo_parallel` still uses `--text` — that is what makes it parallel — and now
  emits a `UserWarning` (`NO_QVALUES_IN_TEXT_MODE`) on every call stating that its
  q-value column is empty by construction.

### Added

- Real-execution tests: default `run_fimo` yields non-empty q-values in [0, 1] and
  `text=True` yields none, in one test; `run_fimo_parallel` warns. Both were run against
  0.1.0 first and failed for the stated reason.

## [0.1.0] - 2026-09-02

Initial release: STREME, TOMTOM, FIMO wrappers; header-loss fix in the chunk merge.
