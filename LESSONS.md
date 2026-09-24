# Lessons

- 2026-09-24: package code called at conftest import time (a `skipif(verify_tools())` gate) turned mutmut's forced-fail probe into a conftest ImportError, pytest exit 4, and aborted every nightly mutation run since 09-18 with 0 mutants checked. Gates that run package code are fixtures. Caught by: the nightly mutation gate's 0-checked/aborted check.
