# Code and data for *Mapping the Substrate Dependence of Amide-Coupling Activators from a High-Throughput Reaction Screen*

Kuangbiao Liao, Guangzhou National Laboratory, Guangzhou 510005, China

Everything needed to reproduce every figure, table and numerical value in the manuscript and
Supporting Information. The dataset is bundled complete - no external download required.

## Contents

| file | what it does |
|---|---|
| `REPRODUCE_acylation.py` | recomputes every number in the manuscript and SI from the dataset files and reports each as pass or fail |
| `check_documents.py` | checks that tables transcribed into the manuscript and SI agree with the Source Data behind each figure |
| `build_figures.py` | builds main-text Figures 1-4 and their Source Data |
| `build_si_figures.py` | builds Supplementary Figures S1, S3, S4 and their Source Data |
| `build_null_figure.py` | builds Supplementary Figure S2 and its Source Data |
| `acyl_descriptors.py` | the frozen descriptor definitions used by the figure builders |
| `make_dataset.py` | extracted `data/DATASET_acylation_wells.csv` and `data/DATASET_acylation_conditions.xlsx` from the in-house corpus; included for provenance, not needed to reproduce the manuscript |
| `MANUSCRIPT_acylation.md`, `SI_acylation.md` | the text `check_documents.py` checks tables against |
| `source_data/` | one CSV per figure; rows correspond to the plotted elements |
| `frozen_folds/` | the five-fold scaffold partition used for every five-fold result, saved as data so it loads identically on any platform (see `frozen_folds/README.md`) |
| `preregistration/` | the frozen files behind SI Note 2's registration record (see `preregistration/README.md`) |
| `data/` | the two dataset files, complete, plus checksums |
| `LICENSE` | MIT |

## Setup

```
python -m venv env
env/Scripts/activate          # Windows;  source env/bin/activate on Unix
python -m pip install -r requirements.txt
```

## Reproducing

```
python REPRODUCE_acylation.py             # all checks
python REPRODUCE_acylation.py --quick     # skips the permutation checks
python build_figures.py                   # regenerates main-text Figures 1-4 + Source Data
python build_si_figures.py                # regenerates Supplementary Figures S1, S3, S4
python build_null_figure.py               # regenerates Supplementary Figure S2
python check_documents.py                 # figures/ and source_data/ are bundled, so this also
                                           # works before the three lines above are run
```

`REPRODUCE_acylation.py` is deliberately independent: it imports nothing from `build_figures.py`
or `acyl_descriptors.py`, reads no Source Data CSV, and reimplements the acid classifier by a
different algorithm. It does read one frozen input beyond the two dataset files -- the five-fold
scaffold partition in `frozen_folds/`, explained below -- which is checked-in data, not output from
any of this repository's own scripts. A number that agrees only with its own generator has not been
reproduced. Every selection in it is guarded, so a check that matches zero rows aborts rather than
passing silently.

**A cross-platform reproducibility issue, root cause confirmed and fixed by freezing the split.**
On at least one run outside the pinned Windows environment (macOS, same seven pinned library
versions, different OS/Python patch version), 19 of 278 checks disagreed by more than the seeded
permutation noise accounts for, including the two-state descriptor's five-fold result (a one-column
regression, expected R² 0.194, reported as 0.189 on the other machine) - so the disagreement was
not confined to the near-singular 18-descriptor-block and ECFP4-fingerprint fits. Diffing the full
per-acid fold assignment between the two machines confirmed the cause: 20 of the 71 acids landed in
a different one of the five folds, even though both machines agreed on the fold sizes
(25/12/12/11/11). scikit-learn's `GroupKFold` breaks ties among equal-sized scaffold groups using
`numpy.argsort`, an unstable sort whose tie order is not guaranteed identical across numpy builds.

The fix does not change the algorithm or any reported number: the exact five-fold partition already
used for every reported five-fold result (the same one checked across several rounds of review) is
saved as data in `frozen_folds/` and loaded by `REPRODUCE_acylation.py` and `build_si_figures.py`
instead of being recomputed by `GroupKFold` at run time. Loading a CSV and indexing an array
involves no sorting and no floating-point step, so the mechanism removes the platform-dependence by
construction. **This has been independently confirmed on macOS: with the frozen partition loaded
from the same files, all 278 of 278 checks pass, and the 19 disagreements are gone.** "Pass" means
each check's own stated tolerance is met, as it always has - this is not a claim that every
intermediate floating-point value is bit-identical across platforms, only that no reported number
differs by more than its documented tolerance. The loader validates the frozen file on every load
(no duplicate or out-of-range rows, and - where scaffold labels are available - no scaffold split
across two folds).

Permutation and bootstrap procedures are seeded, so repeated runs give identical values.

## Verified environment

- python 3.9.10
- pandas 2.3.3
- numpy 2.0.2
- scipy 1.13.1
- scikit-learn 1.6.1
- matplotlib 3.9.4
- rdkit 2025.09.2
- openpyxl 3.1.5

The permutation p-values were additionally confirmed to be stable across reruns; all randomised procedures are seeded.
