"""Empirical background distribution for the activator-family effect (SFigure 2, Supporting Info).

This robustness check moved from the main text to the SI in the no-new-experiments revision: the
raw paired difference (main-text Figure 2) and the matched additive comparisons (Figure 3) now
carry the primary evidence, and this background/null comparison supports them from SI Note 10.

Scope: the self-contained acylation dataset only (`data/DATASET_acylation_wells.csv` and
`data/DATASET_acylation_conditions.xlsx`, produced from the in-house corpus by `make_dataset.py`).
The sealed rows are deliberately NOT merged - every other figure in this paper is analysis-only,
and a background figure drawn over a wider substrate set than the rest could not be explained to
a referee.

Orientation: the target group (quinoline minus haloformamidinium/carbodiimide) is directional and
is NOT symmetrised - its sign carries chemical meaning. The background pairs are unordered, so the
sign of their difference is an artefact of enumeration order; they ARE symmetrised, by entering
both (x, y) and (y, x). The symmetrised background therefore has mean 0 and 50% negative BY
CONSTRUCTION, and neither is reported as a finding. The only reportable background quantity is
the fraction that reaches the target group's mean.
"""
import os, sys, itertools
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

_here = os.path.dirname(os.path.abspath(__file__))
# the two self-contained files live in ./data beside this script; HTE_DATA overrides that
# location for a data/ folder placed somewhere else after unpacking the archive
DATA = os.environ.get("HTE_DATA") or os.path.join(_here, "data")
OUT = os.path.dirname(os.path.abspath(__file__))
FIG, SRC = os.path.join(OUT, "figures"), os.path.join(OUT, "source_data")
MIN_ACIDS = 5

plt.rcParams.update({"font.family": "Times New Roman", "font.size": 8, "axes.linewidth": 0.7,
                     "xtick.major.width": 0.7, "ytick.major.width": 0.7, "mathtext.fontset": "stix",
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 400})

# ---- explicit family membership. Not a regex: every member is written out and printed. ----
QUINOLINE = {"EEDQ", "IIDQ"}
HALOFORMAMIDINIUM = {"BTFFH", "TFFH", "TCFH"}
CARBODIIMIDE = {"DCC", "EDC-HCl", "DIC"}

COOH = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
def fused(mol, start):
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()
             if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    s, fr, seen = set(), [r for r in rings if start in r], []
    while fr:
        r = fr.pop()
        if r in seen: continue
        seen.append(r); s |= r
        for o in rings:
            if o not in seen and o & r: fr.append(o)
    return s
def acid_state(smi):
    m = Chem.MolFromSmiles(smi); h = m.GetSubstructMatches(COOH)[0]
    c = m.GetAtomWithIdx(h[0])
    a = [n for n in c.GetNeighbors() if n.GetIdx() not in h[1:]][0]
    if not a.GetIsAromatic(): return "C"
    sysm = fused(m, a.GetIdx())
    return "A" if any(m.GetAtomWithIdx(i).GetSymbol() != "C" for i in sysm) else "B"

print("loading (analysis rows only; sealed deliberately not merged) ...")
ac = pd.read_csv(os.path.join(DATA, "DATASET_acylation_wells.csv"), low_memory=False)
ac = ac.dropna(subset=["yield"]).copy()
cond = pd.read_excel(os.path.join(DATA, "DATASET_acylation_conditions.xlsx"))
cA = cond.set_index("condition_id")
conds = sorted(set(ac.condition_id.unique()) & set(cA.index))
assert len(conds) > 0, "no conditions matched"
act = cA["activator_name"].to_dict()
print(f"  {len(conds)} conditions defined and used -> {len(conds)*(len(conds)-1)//2} unordered pairs")

# ---- print family members for line-by-line checking, and verify equal loading ----
print("\nfamily membership, printed for checking:")
fam_of = {}
for label, members in [("quinoline", QUINOLINE), ("haloformamidinium", HALOFORMAMIDINIUM),
                       ("carbodiimide", CARBODIIMIDE)]:
    picked = [c for c in conds if act[c] in members]
    eqs = sorted({float(cA.loc[c, "activator_equiv"]) for c in picked})
    print(f"  {label:20s} {len(picked):2d} conditions, activator equivalents {eqs}")
    for c in picked:
        fam_of[c] = label
        print(f"      {c:5s} {act[c]:9s} equiv {cA.loc[c, 'activator_equiv']}")
    assert eqs == [1.5], f"{label} is not at a single activator loading: {eqs}"
unassigned = sorted({act[c] for c in conds} - QUINOLINE - HALOFORMAMIDINIUM - CARBODIIMIDE)
print(f"  activators in no target family ({len(unassigned)}): {', '.join(unassigned)}")

acids = sorted(ac.sub_2_smiles.unique())
STATE = {s: acid_state(s) for s in acids}
print(f"\nacids: {len(acids)}  heteroaromatic {sum(v=='A' for v in STATE.values())}  "
      f"non-aromatic {sum(v=='C' for v in STATE.values())}")

# ---- subpair x condition matrix, and a subpair -> acid grouping ----
piv = ac.pivot_table(index="subpair", columns="condition_id", values="yield", aggfunc="mean")
piv = piv.reindex(columns=conds)
M = piv.to_numpy(dtype=float)
sp_acid = ac.drop_duplicates("subpair").set_index("subpair").sub_2_smiles
ai = {a: i for i, a in enumerate(acids)}
rows = np.array([ai[sp_acid[s]] for s in piv.index])
G = np.zeros((len(acids), len(piv)), dtype=float)
G[rows, np.arange(len(piv))] = 1.0
isA = np.array([STATE[a] == "A" for a in acids])
isC = np.array([STATE[a] == "C" for a in acids])
assert isA.sum() >= MIN_ACIDS and isC.sum() >= MIN_ACIDS

recs = []
for i, j in itertools.combinations(range(len(conds)), 2):
    d = M[:, i] - M[:, j]
    ok = ~np.isnan(d)
    if not ok.any():
        continue
    dv = np.where(ok, d, 0.0)
    cnt = G @ ok.astype(float)
    tot = G @ dv
    have = cnt > 0
    if (have & isA).sum() < MIN_ACIDS or (have & isC).sum() < MIN_ACIDS:
        continue
    mean_acid = np.divide(tot, cnt, out=np.zeros_like(tot), where=have)
    diff = mean_acid[have & isA].mean() - mean_acid[have & isC].mean()
    recs.append({"cond_x": conds[i], "cond_y": conds[j],
                 "activator_x": act[conds[i]], "activator_y": act[conds[j]],
                 "family_x": fam_of.get(conds[i], "other"),
                 "family_y": fam_of.get(conds[j], "other"),
                 "n_subpairs": int(ok.sum()),
                 "n_acids_hetero": int((have & isA).sum()),
                 "n_acids_nonarom": int((have & isC).sum()),
                 "diff_enumeration_order": float(diff)})
P = pd.DataFrame(recs)
assert len(P) > 0, "no estimable condition pairs"
assert not P.diff_enumeration_order.isna().any(), "NaN in the pair differences"
print(f"\nestimable condition pairs: {len(P)}")

# ---- target vs background ----
tgt_mask = (((P.family_x == "quinoline") & (P.family_y.isin(["haloformamidinium", "carbodiimide"]))) |
            ((P.family_y == "quinoline") & (P.family_x.isin(["haloformamidinium", "carbodiimide"]))))
P["group"] = np.where(tgt_mask, "target", "background")
# orient the target group as quinoline minus the other arm; this direction is chemical, not clerical
sign = np.where(P.family_x == "quinoline", 1.0, -1.0)
P["diff_oriented"] = np.where(tgt_mask, P.diff_enumeration_order * sign, np.nan)
T = P[P.group == "target"]
B = P[P.group == "background"]
print(f"  target pairs (quinoline vs haloformamidinium/carbodiimide): {len(T)}")
print(f"  background pairs: {len(B)}")
tgt_mean = float(T.diff_oriented.mean())
tgt_pos = int((T.diff_oriented > 0).sum())
print(f"  target mean (quinoline minus other) = {tgt_mean:+.4f}")
print(f"  target positive: {tgt_pos}/{len(T)} = {100*tgt_pos/len(T):.1f}%")

# ---- background, symmetrised ----
bg_fwd = B.diff_enumeration_order.to_numpy()
bg_sym = np.concatenate([bg_fwd, -bg_fwd])
neg_fwd = float((bg_fwd < 0).mean())
neg_rev = float(((-bg_fwd) < 0).mean())
print(f"\n  background negative fraction, enumeration order : {neg_fwd:.4f}")
print(f"  background negative fraction, reversed order    : {neg_rev:.4f}")
print(f"  the two sum to {neg_fwd + neg_rev:.6f} - the sign is an artefact of ordering,")
print(f"  so neither the background mean ({bg_sym.mean():+.2e}) nor its 50.0% sign split is a finding.")
# each of the len(bg_fwd) original pairs contributes at most one of its two mirrored entries to
# a one-sided count above a positive threshold, so counting against the true number of distinct
# pairs (not the doubled symmetrised set) is the correct, non-diluted reaching-fraction
n_reach = int((np.abs(bg_fwd) >= tgt_mean).sum())
reach = n_reach / len(bg_fwd)
print(f"\n  REPORTABLE: background pairs reaching the target mean of {tgt_mean:+.4f}:")
print(f"    {n_reach} of {len(bg_fwd)} distinct background pairs = {100*reach:.2f}%")

P["is_target"] = P.group == "target"
P.to_csv(os.path.join(SRC, "SFigure2_source_data.csv"), index=False)

# ---- figure ----
fig, ax = plt.subplots(figsize=(7.2, 3.2))
bins = np.linspace(-0.32, 0.32, 81)
ax.hist(bg_sym, bins=bins, color="#BFBFBF", edgecolor="white", linewidth=0.2,
        label=f"background condition pairs, symmetrised (n = {len(bg_sym):,})", density=True)
ax.hist(T.diff_oriented, bins=bins, color="#053061", alpha=0.85, edgecolor="white",
        linewidth=0.2, label=f"quinoline $-$ haloformamidinium/carbodiimide (n = {len(T)})",
        density=True)
ax.axvline(0, color="k", lw=0.6, ls=":")
ax.axvline(tgt_mean, color="#053061", lw=1.4)
ax.annotate(f"target mean {tgt_mean:+.3f}\n{100*reach:.1f}% of background reaches it",
            xy=(tgt_mean, ax.get_ylim()[1] * 0.72), xytext=(tgt_mean + 0.055, ax.get_ylim()[1] * 0.80),
            fontsize=7, color="#053061",
            arrowprops=dict(arrowstyle="->", lw=0.7, color="#053061"))
ax.set_xlabel("mean yield difference, heteroaromatic acids $-$ non-aromatic acids")
ax.set_ylabel("density")
ax.set_title("Background distribution over all condition pairs", loc="left",
             fontsize=9, fontweight="bold")
ax.legend(fontsize=6.5, frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "SFigure2.png"), bbox_inches="tight")
plt.close(fig)
print(f"\nwrote {os.path.join(FIG, 'SFigure2.png')}")
print(f"wrote {os.path.join(SRC, 'SFigure2_source_data.csv')}  ({len(P)} rows, one per condition pair)")
