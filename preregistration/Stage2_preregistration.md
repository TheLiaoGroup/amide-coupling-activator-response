# Pre-registration — Stage 2 (out-of-sample prediction from structure)

Written and frozen **before** any Stage 2 quantity was computed. Companion to
`PREREGISTRATION_acylation_stage1.md` (frozen `708e522f…`, Amendment 1 `5debf88d…`).

## Decisions carried in from Stage 1

Settled by the project owner after the Stage 1 report, and fixed here so Stage 2 cannot drift:

1. **The descriptor is two-state** — heteroaromatic ring system versus not. Justified without
   selecting a response, because the binary form beat or tied the three-state form on all five
   label-free responses tried.
2. **The three-way breakdown stays visible** in figures and tables as descriptive, explicitly
   labelled as not established, with its p-values (0.056–0.076) shown.
3. **The electronic co-ordinate is reported beside it**, with r = 0.704 stated. The paper does
   not claim the ring-system assignment is uniquely best.
4. The unresolved Stage 1 gate, the five-response table and the control failure on
   `pc1_resid_loading` are disclosed in the SI.

## Target

**Continuous target `y_familycontrast`** — each acid's mean within-condition z-score over the
quinoline-type conditions (EEDQ C41, IIDQ C95) minus its mean over the eight
fluoroformamidinium/carbodiimide conditions. The two groups are *a priori reagent families*, so
the target is built from chemical prior knowledge and yields, never from an acid class label.

**Binary target** — the sign of that contrast. Positive means the acid prefers the quinoline end.

`y_endcontrast` (the same construction with CDMT/DEPBT/PyOxim in place of the quinolines) is
reported alongside as a robustness check, not as the headline.

## Models, all with structure as the only input

| model | features |
|---|---|
| M1 | the two-state descriptor alone — the paper's claim |
| M2 | the frozen Stage 1 descriptor block D1–D18 |
| M3 | ECFP4 (radius 2, 2048 bits) with ridge, α fixed at 1.0, chosen before running and not tuned |
| M0 | baseline: predict the training mean (continuous) or the majority class (binary) |

## Protocol

- Unit: the acid, n = 71. Splitting: `GroupKFold`, 5 folds, Bemis–Murcko scaffold groups.
  No scaffold is split. Random splitting is forbidden.
- No feature selection, no hyper-parameter tuning inside or outside the folds. α = 1.0 is fixed.
- Scores: out-of-fold R² for the continuous target; out-of-fold ROC-AUC and accuracy for the
  binary target; both against M0.
- Significance: label permutation, 5,000 draws, respecting the fold structure.
- Leave-one-scaffold-out is reported in addition to 5-fold, because the largest scaffold group
  holds 25 of 71 acids and 5-fold therefore hides how much of the score rides on that one group.

## What will be reported

Every number, including bad ones. The briefing's instruction is explicit — report honestly, and
if prediction is weak, say it is weak. Specifically:

- If M1 does not beat M0, that is the result and it is stated in the abstract, not buried.
- If M3 (a generic fingerprint that knows nothing about heteroaromaticity) beats M1, that is
  reported as evidence the named descriptor is not capturing all of the structural signal.
- The per-fold spread is reported, not only the pooled score, so a single lucky fold cannot be
  presented as predictive skill.
- Class balance for the binary target is reported before AUC, since a 71-acid set with an uneven
  split makes accuracy uninformative on its own.
