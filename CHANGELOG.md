# Changelog

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions: SemVer.

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
