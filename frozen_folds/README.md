# Frozen five-fold partitions

`GroupKFold(n_splits=5)` breaks ties among equal-sized Bemis–Murcko scaffold groups using
`numpy.argsort`, an unstable sort whose tie order is not guaranteed identical across numpy
builds. Diffing the actual per-acid fold assignment between a Windows run and a macOS run
(same seven pinned library versions, different OS/Python patch version) found 20 of 71 acids
placed in a different one of the five folds, even though both machines agreed on the fold
sizes (25/12/12/11/11) — this is the confirmed cause of the cross-platform disagreement in the
five-fold-derived numbers (`REPRODUCE_acylation.py`, `build_si_figures.py`).

Rather than recompute the split with a different (order-dependent) algorithm — which would
change roughly two dozen already-reported numbers across the manuscript, SI and figures — the
partition actually used for every reported five-fold result is saved here as data and loaded
by both scripts instead of being recomputed by `GroupKFold` at run time. This is a saved
analysis split, not a pre-registration: the numbers it reproduces are the same numbers already
checked across several rounds of review before this file existed.

- `folds_71acids.csv` — the 71-acid split used wherever the fold membership only depends on
  which acids have a value for the response (the section 8 structure-prediction target, four of
  the five Stage 1 contest responses, and the Supplementary Figure S3 target). Verified
  identical across all of those call sites before this file was created.
- `folds_C3_C41_70acids.csv` — the separate 70-acid split for the Stage 1 contest's
  `C3_C41` (BTFFH − EEDQ) response, which drops one acid lacking a shared amine at both
  conditions and therefore groups differently.

Each loader asserts the frozen file's acid set exactly matches the acids being split before
using it, so a change to the dataset that changes which acids have a value would fail loudly
rather than silently reusing a stale partition.
