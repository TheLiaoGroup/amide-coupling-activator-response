# Pre-registration — Stage 3 (the amine axis)

Written and frozen **before** any Stage 3 quantity was computed. Companion to the Stage 1 and
Stage 2 registrations (`708e522f…` / Amendment 1 `5debf88d…`, and `9110c3b6…`).

## The question

The paper claims the activator preference belongs to the **acid**. If the effect is in fact
carried by the amine — if it appears only with some amine types and vanishes with others — then
the central proposition is wrong as stated and must be rewritten. Stage 3 is the test that can
falsify it.

## Amine classification — frozen before computing

All 82 amines classified from `sub_1_smiles` alone, by the nitrogen that forms the amide. Where a
molecule carries more than one candidate nitrogen, the least hindered primary/secondary amine
nitrogen not in an amide, carbamate or sulfonamide is taken; the rule is applied mechanically and
the count of ambiguous cases is reported.

| level | definition |
|---|---|
| **AR-C** | aniline — N attached directly to a carbocyclic aromatic ring |
| **AR-H** | heteroaryl amine — N attached directly to an aromatic ring bearing a heteroatom |
| **AL-1** | aliphatic primary amine — N with two H, attached to sp3 carbon |
| **AL-2** | aliphatic secondary amine — N with one H, two heavy neighbours, not aromatic-attached |

A secondary companion split, reported alongside: **primary vs secondary** by hydrogen count
alone, which is the cruder but less arguable division.

## Tests

Response is the same per-acid-per-amine quantity used in Stage 2, rebuilt at cell level: the
within-condition z-score of a substrate pair's yield, contrasted between the quinoline-type
conditions and the eight fluoroformamidinium/carbodiimide conditions.

- **T1 — three-way.** Acid state (2 levels, per the Stage 1 decision) × amine class (4 levels),
  on the family contrast. Reports the acid-state effect **within each amine class separately**.
- **T2 — homogeneity.** Does the acid-state effect differ across amine classes? Permutation test
  on the dispersion of the four within-class acid-state effects, shuffling amine class labels
  across amines, 20,000 draws.
- **T3 — sign consistency.** The fraction of individual amines for which the acid-state effect
  has the predicted sign, with a binomial test against 0.5.

## Units, and the rule that decides them

The acid is the independent unit for any claim about acids. For T3 the **amine** is the unit,
because the claim there is about amines. No test uses the substrate pair as the unit — 830-style
inflated n is exactly the failure this project already recorded.

## Decision rule, fixed in advance

- **Proposition survives** if the acid-state effect has the predicted sign in all four amine
  classes and T2 finds no significant heterogeneity (p ≥ 0.05).
- **Proposition survives with qualification** if the sign holds in all four but T2 shows
  significant heterogeneity — the effect is general in direction, modulated in size, and the
  modulation is reported in the main text.
- **Proposition must be rewritten** if the acid-state effect reverses sign in any amine class
  containing at least eight amines. In that case the paper reports an acid × amine × activator
  three-way effect and drops the claim that preference is a property of the acid alone.
