# Pre-registration — Stage 1 (descriptor selection) and Stage 1b (cross-family support)

Written and frozen **before** any Stage 1 or Stage 1b quantity was computed. The SHA-256 of this
file is recorded in `WORKING_NOTES_acylation.md` at the moment of freezing; any later edit
invalidates the registration and must be recorded as an amendment with its own hash.

Every descriptor is a function of `sub_2_smiles` alone — the acid. No yield, no condition, and no
class label enters any descriptor.

---

## A. The circularity problem, and the response variable that avoids it

`z_gap` is *defined* as mean z(heteroaromatic) − mean z(non-aromatic). Scoring descriptors on
their ability to explain `z_gap` would hand the incumbent the contest by construction: the
incumbent is the definition of the axis. The Stage 1 gate would then be vacuous — the failure
mode this project has already recorded twice.

So Stage 1 scores descriptors against a response in which **no class label appears**.

**Primary response — `pc1_loading`.** Build the 46 × 71 matrix of within-condition acid z-scores
(each condition's acids standardised across acids, exactly as in §3.3). Double-centre it by
subtracting condition means and acid means, leaving pure substrate × condition interaction.
Take the SVD. Each acid's loading on the first right singular vector is its position on the
data's own leading interaction mode. This is computed without reference to any acid class.

Sign convention: the singular vector is signed so that its correlation with the conditions'
mean yield ordering is not used; instead the sign is fixed so that the loading of the acid
`O=C(O)c1ccco1` (furan-2-carboxylic acid, an arbitrary but fixed anchor chosen before seeing
results) is positive. Sign has no effect on R² or on |correlation|.

**Secondary response — `y_C3_C41`.** Each acid's amine-matched C3 − C41 difference. Label-free,
but it uses two conditions that were themselves selected as the axis extremes, so it is reported
as secondary and never as the basis for the gate (§6 rule 3).

**Diagnostic, not a scoring target.** The correlation between `pc1_loading` and `z_gap`-derived
quantities is reported for interpretation only.

---

## B. Candidate descriptors — frozen list

All eighteen are computed for all 71 acids. Attachment atom = the non-carbonyl neighbour of the
carboxyl carbon. Attached ring system = the fused aromatic ring system containing the attachment
atom, when the attachment atom is aromatic.

**Incumbent**

| id | name | type | definition |
|---|---|---|---|
| D1 | `ring_system_state` | 3-level | A heteroaromatic / B carbocyclic aryl / C non-aromatic, on the fused ring system |

**Collapses of the incumbent — do fewer states do as well?**

| id | name | type | definition |
|---|---|---|---|
| D2 | `aromatic_attachment` | binary | (A or B) vs C |
| D3 | `hetero_vs_rest` | binary | A vs (B or C) |

**Locality variants — does it matter where the heteroatom is?**

| id | name | type | definition |
|---|---|---|---|
| D4 | `attached_ring_hetero` | binary | heteroatom in the single ring bearing the attachment, not the fused system |
| D5 | `molecule_has_arom_hetero` | binary | any aromatic heteroatom anywhere in the molecule |

**Heteroatom identity and count**

| id | name | type | definition |
|---|---|---|---|
| D6 | `hetero_N` / `hetero_O` / `hetero_S` | 3 binaries | element present in attached ring system (entered jointly) |
| D7 | `n_ring_heteroatoms` | integer | count of heteroatoms in attached aromatic ring system |

**Ring system identity**

| id | name | type | definition |
|---|---|---|---|
| D8 | `ring_system_class` | categorical | benzene, naphthalene, pyridine, pyrimidine, pyrazine, quinoline, quinoxaline, thiophene, furan, benzothiophene, benzofuran, imidazole, pyrazole, indole, aliphatic |

**Position — the scale §3.5 rejected, re-tested here rather than assumed dead**

| id | name | type | definition |
|---|---|---|---|
| D9 | `dist_to_ring_hetero` | ordinal | shortest bond path, carboxyl C to nearest aromatic ring heteroatom; no heteroatom coded as its own level |
| D10 | `hetero_adjacent` | binary | a ring heteroatom is bonded directly to the attachment atom (the d2 set) |

**Electronic — computed acidity proxies**

| id | name | type | definition |
|---|---|---|---|
| D11 | `q_hydroxyl_O` | continuous | Gasteiger partial charge on the carboxyl hydroxyl oxygen |
| D12 | `q_carbonyl_C` | continuous | Gasteiger partial charge on the carboxyl carbon |

**Steric — Taft/Charton-like**

| id | name | type | definition |
|---|---|---|---|
| D13 | `alpha_heavy_degree` | integer | heavy-atom degree of the attachment atom |
| D14 | `n_ortho_substituents` | integer | non-hydrogen substituents on atoms adjacent to the attachment atom |
| D15 | `n_rotatable_bonds` | integer | RDKit `NumRotatableBonds` |

**Controls — pre-declared to be expected null. If any of these wins, the analysis is broken, not interesting.**

| id | name | type |
|---|---|---|
| D16 | `MolWt` | continuous |
| D17 | `TPSA` | continuous |
| D18 | `MolLogP` | continuous |

---

## C. Scoring protocol — frozen

- **Unit**: the acid. n = 71.
- **Groups**: Bemis–Murcko scaffold of the acid. Cross-validation is `GroupKFold` on these
  groups, so no scaffold is split across a fold boundary. Random splitting is forbidden.
- **Folds**: 5. If any scaffold group is larger than one fold's worth, GroupKFold's own
  balancing applies; the realised fold sizes are reported.
- **Model**: ordinary least squares on the descriptor's encoding (one-hot with dropped first
  level for categoricals, raw value for continuous, standardised). One descriptor per model —
  this is a comparison of descriptors, not a search for the best combination.
- **Score**: out-of-fold R², pooled across folds against the response's overall mean. Reported
  alongside in-fold R² so overfitting is visible.
- **Significance**: permutation of the descriptor's values across acids **within the CV
  structure**, 10,000 draws, giving a p for each descriptor's out-of-fold R².
- **Every score is reported for every descriptor**, including ones that do badly. No descriptor
  is dropped from the table for performing poorly.

### Stage 1 gate — decision rule fixed in advance

- **Incumbent holds** if D1's out-of-fold R² on the primary response is greater than or equal to
  every alternative's, within one standard error of the fold-to-fold spread.
- **Incumbent is replaced** if any alternative beats D1 by more than one standard error on the
  primary response. The replacement becomes the paper's descriptor and the central proposition
  is rewritten around it.
- **Controls D16–D18 must be non-significant.** If a control is significant, the CV structure is
  reported as broken and Stage 1 restarts rather than proceeding.

---

## D. Stage 1b — does the positive end survive without the quinolines?

The positive end of the axis is EEDQ (C41) and IIDQ (C95), one condition each, no within-reagent
replication, both quinoline-type. Behind them sit three conditions from three unrelated families:
CDMT (C64, triazine) +0.58, DEPBT (C42, phosphonate benzotriazinone) +0.46, PyOxim (C37,
oxime-uronium) +0.30.

**Both quinoline conditions are deleted from the array. The axis is rebuilt on the remaining 44
conditions.** All four tests below run on that reduced array.

- **Test 1 — individual.** For each of C64, C42, C37: `z_gap`, and a p from permuting the acid
  class labels (20,000 draws, acid is the unit).
- **Test 2 — end against end.** Per-acid contrast between the three non-quinoline positive-end
  conditions and the eight fluoroformamidinium/carbodiimide conditions. Acid is the unit;
  significance by permutation of the acid class label.
- **Test 3 — dispersion.** The §3.3 interaction test rerun on the 44 conditions: SD of `z_gap`
  against the permutation null.
- **Test 4 — the yield confound, which must be controlled.** `z_gap` correlates with condition
  mean yield at r = 0.566 across the 46. The two quinoline conditions are *low*-yielding (0.138,
  0.139) and sit at the top, so they are the evidence that the axis is not merely a yield axis.
  Deleting them removes exactly that evidence: CDMT, DEPBT and PyOxim are all high-yielding
  (0.380, 0.346, 0.463). So Test 4 recomputes r(`z_gap`, mean yield) on the 44, regresses `z_gap`
  on condition mean yield, and asks whether the three families stay positive **in the residual**.

### Stage 1b branch — decision rule fixed in advance

**Cross-family support holds**, and goes in the main text, only if all three of:

1. all three of CDMT, DEPBT, PyOxim have positive `z_gap` on the reduced array;
2. Test 2 is significant at p < 0.05 with the acid as the unit;
3. Test 4 leaves the three positive in the residual after condition mean yield is regressed out.

**Otherwise the claim narrows.** If (1) or (2) fails, the axis is a quinoline phenomenon: the
title, abstract and central proposition are restricted to quinoline-type activators and may not
say "activator class" without qualification. If only (3) fails, the positive end is reported as
confounded with condition productivity, the cross-family claim is dropped from the main text, and
the limitation is stated in the main text rather than the SI.

The branch is resolved before a word of the manuscript is written.

---

# AMENDMENT 1 — logged 2026-08-31, after Stage 1 primary scoring, before the amended response was computed

The pre-registered primary response `pc1_loading` is label-free, as intended, but a diagnostic run
after scoring shows it does **not** isolate the construct Stage 1 is about. PC1 of the
double-centred matrix correlates with condition mean yield at **r = +0.879**, and EEDQ and IIDQ —
the two conditions that define the top of the preference axis — rank only 17th and 16th of 46 on
it. PC1 is the array's productivity mode: it separates acids that do relatively better in
high-yielding conditions. That is a real mode and worth reporting, but it is not the preference
mode, so a descriptor contest run on it answers a different question from the one asked.

This is a miss in the original registration, recorded rather than quietly repaired. The
pre-registered result stands and is reported in full; the amendment adds one response, it does
not replace or delete anything.

**Amended primary response — `pc1_resid_loading`.** Before the SVD, project the productivity
direction out of the double-centred matrix: let `v` be the centred vector of condition mean
yields over the 46 conditions, and replace `D` with `D − v (vᵀD)/(vᵀv)`. Take the SVD of the
residual and use each acid's loading on the first right singular vector, signed by the same
furan-2-carboxylic-acid anchor. Still label-free — condition mean yield is a property of the
condition, not of any acid class.

Scoring protocol, groups, folds, model, permutation count and the control requirement are
unchanged. The Stage 1 gate is now read on `pc1_resid_loading`, with `pc1_loading` and
`y_C3_C41` both reported beside it. If the three responses disagree about which descriptor wins,
that disagreement is reported as the result — no response is selected after the fact for giving
the preferred answer.

**Additional pre-specified check.** D5 differs from D3 by exactly one molecule. Any descriptor
whose advantage does not survive deletion of the single molecule that distinguishes it from a
simpler descriptor is reported as one-molecule fragility and is not eligible to replace the
incumbent.
