"""Build the Supplementary figures (activator axis, quadrangle interactions, full descriptor/ECFP4
comparison) and their Source Data. These three analyses moved out of the main text in the
no-new-experiments revision: they remain as robustness/detail evidence, cited from SI Notes, not
as primary evidence for the headline claim (the raw paired difference in Figure 2 and the matched
additive comparisons in Figure 3 carry that).

The background/null distribution (Figure 5 of the previous version) is built separately by
`build_null_figure.py`, which now writes SFigure2 instead of a main-text figure.

Reads `data/DATASET_acylation_wells.csv` and `data/DATASET_acylation_conditions.xlsx`.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog("rdApp.*")

_here = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("HTE_DATA") or os.path.join(_here, "data")
OUT = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(OUT, "figures")
SRC = os.path.join(OUT, "source_data")
os.makedirs(FIG, exist_ok=True); os.makedirs(SRC, exist_ok=True)
SEED = 20260901
FOLDS = os.path.join(_here, "frozen_folds")


def load_folds(acid_order, filename, groups=None):
    """Load a frozen five-fold partition (see frozen_folds/README.md) instead of recomputing it
    live with GroupKFold, whose tie-breaking among equal-sized scaffold groups was found to
    differ across numpy builds/platforms even at the same pinned library version (confirmed fixed
    by this mechanism: 278/278 on both the original Windows machine and an independent macOS run).
    Validates the frozen file before using it: no duplicate or malformed rows, fold indices in
    range, acid set exactly matching the current run, and - when `groups` (e.g. Murcko scaffold
    labels, aligned with `acid_order`) is supplied - no scaffold split across two folds."""
    frozen = {}
    with open(os.path.join(FOLDS, filename), encoding="utf-8") as f:
        next(f)
        for lineno, line in enumerate(f, start=2):
            a, k = line.rstrip("\n").rsplit(",", 1)
            assert a not in frozen, f"{filename}:{lineno}: duplicate acid {a!r}"
            k = int(k)
            assert 0 <= k < 5, f"{filename}:{lineno}: fold index {k} out of range for acid {a!r}"
            frozen[a] = k
    acid_order = list(acid_order)
    assert set(frozen) == set(acid_order), f"{filename}: frozen acid set does not match the current run"
    fold_of = np.array([frozen[a] for a in acid_order])
    if groups is not None:
        groups = np.asarray(groups)
        for gg in set(groups):
            assert len(set(fold_of[groups == gg])) == 1, f"{filename}: scaffold {gg!r} split across folds"
    idx = np.arange(len(acid_order))
    return [(idx[fold_of != k], idx[fold_of == k]) for k in range(5)]


plt.rcParams.update({"font.family": "Times New Roman", "font.size": 8, "axes.linewidth": 0.7,
                     "xtick.major.width": 0.7, "ytick.major.width": 0.7, "mathtext.fontset": "stix",
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 400})

FAMILY = {}
for a in ["TFFH", "BTFFH", "TCFH"]: FAMILY[a] = "Haloformamidinium"
for a in ["DCC", "DIC", "EDC-HCl"]: FAMILY[a] = "Carbodiimide"
for a in ["DCMT", "DMTMM", "CDMT"]: FAMILY[a] = "Triazine"
for a in ["EEDQ", "IIDQ"]: FAMILY[a] = "Dihydroquinoline"
for a in ["PyBOP", "PyAOP", "PyBrOP", "PyClOP", "PyClocK", "BOP", "BOPCl", "DPP-Cl",
          "DEPBT", "FDPP", "(EtO)2P(O)CN", "PyOxim"]: FAMILY[a] = "Phosphorus"
for a in ["HATU", "HBTU", "HCTU", "TCTU", "TBTU", "TDBTU", "TNTU", "TPTU", "TSTU",
          "HDMC", "CITU", "PyCIU", "HBPYU", "TOTU", "COMU", "CIP"]: FAMILY[a] = "Aminium / uronium"
for a in ["CDI", "IBCF"]: FAMILY[a] = "Other"
COL = {"Haloformamidinium": "#B2182B", "Carbodiimide": "#D6604D", "Triazine": "#4393C3",
       "Dihydroquinoline": "#053061", "Phosphorus": "#92C5DE",
       "Aminium / uronium": "#BFBFBF", "Other": "#7F7F7F"}
FAMILY_ORDER = ["Haloformamidinium", "Carbodiimide", "Triazine", "Dihydroquinoline",
                "Phosphorus", "Aminium / uronium", "Other"]
STATE_COL = {"A": "#762A83", "B": "#9A9A9A", "C": "#1B7837"}

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
    mol = Chem.MolFromSmiles(smi); m = mol.GetSubstructMatches(COOH)[0]
    c = mol.GetAtomWithIdx(m[0])
    a = [n for n in c.GetNeighbors() if n.GetIdx() not in m[1:]][0]
    if not a.GetIsAromatic(): return "C"
    sysm = fused(mol, a.GetIdx())
    return "A" if any(mol.GetAtomWithIdx(i).GetSymbol() != "C" for i in sysm) else "B"

print("loading data ...")
ac = pd.read_csv(os.path.join(DATA, "DATASET_acylation_wells.csv"), low_memory=False).dropna(subset=["yield"])
cond = pd.read_excel(os.path.join(DATA, "DATASET_acylation_conditions.xlsx"))
cA = cond.set_index("condition_id")
acid_smiles = sorted(ac.sub_2_smiles.unique())
STATE = {s: acid_state(s) for s in acid_smiles}
ac["state"] = ac.sub_2_smiles.map(STATE)

keep = []
for c in sorted(ac.condition_id.unique()):
    g = ac[ac.condition_id == c]
    if pd.notna(cA.loc[c, "additive_CAS"]): continue
    if not (cA.loc[c, "solvent_name"] == "DMF" and cA.loc[c, "temp_C"] == 25.0
            and cA.loc[c, "time_h"] == 2.0): continue
    if len(g) <= 300: continue
    keep.append(c)
rows = []
for c in keep:
    g = ac[ac.condition_id == c]
    am = g.groupby("sub_2_smiles")["yield"].mean()
    s = am.index.map(STATE)
    z = (am - am.mean()) / am.std(ddof=1)
    act = cA.loc[c, "activator_name"]
    base = cA.loc[c, "base_name"] if pd.notna(cA.loc[c, "base_name"]) else "none"
    rows.append({"condition_id": c, "activator": act, "base": base,
                 "reagent_family": FAMILY.get(act, "Other"), "n_wells": len(g),
                 "mean_yield": g["yield"].mean(),
                 "z_hetero": z[s == "A"].mean(), "z_carbocyclic": z[s == "B"].mean(),
                 "z_nonaromatic": z[s == "C"].mean(),
                 "z_gap": z[s == "A"].mean() - z[s == "C"].mean()})
AX = pd.DataFrame(rows).sort_values("z_gap").reset_index(drop=True)
sl, ic = stats.linregress(AX.mean_yield, AX.z_gap)[:2]
AX["z_gap_residual"] = AX.z_gap - (ic + sl * AX.mean_yield)
print("  46-condition axis built:", len(AX), "conditions")

# ================================================================== SFIGURE 1: the activator axis
fig = plt.figure(figsize=(7.2, 7.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.55, 1.0], hspace=0.52, wspace=0.28)
axA = fig.add_subplot(gs[0, :]); axB = fig.add_subplot(gs[1, 0]); axC = fig.add_subplot(gs[1, 1])
o = AX.sort_values("z_gap")
axA.bar(range(len(o)), o.z_gap, color=[COL[f] for f in o.reagent_family], width=0.8,
        edgecolor="white", linewidth=0.3)
axA.axhline(0, color="k", lw=0.7)
axA.set_xticks(range(len(o))); axA.set_xticklabels([f"{a}" for a in o.activator], rotation=90, fontsize=5.6)
axA.set_ylabel("z(heteroaromatic) $-$ z(non-aromatic)")
axA.set_xlim(-0.8, len(o) - 0.2)
axA.set_title("a   How the 46 single-activator, additive-free conditions were chosen",
              loc="left", fontsize=9, fontweight="bold")
axA.legend(handles=[Patch(facecolor=COL[f], label=f) for f in FAMILY_ORDER],
           fontsize=6, frameon=False, ncol=4, loc="upper left", handlelength=1.1)
axB.scatter(AX.mean_yield, AX.z_gap, s=26, c=[COL[f] for f in AX.reagent_family],
            edgecolor="k", linewidth=0.3, zorder=3)
xs = np.linspace(AX.mean_yield.min(), AX.mean_yield.max(), 50)
axB.plot(xs, ic + sl * xs, "k--", lw=0.8, zorder=2)
for cid, dx, dy in [("C41", 6, 4), ("C95", 6, -10)]:
    r = AX[AX.condition_id == cid].iloc[0]
    axB.annotate(r.activator, (r.mean_yield, r.z_gap), textcoords="offset points",
                 xytext=(dx, dy), fontsize=6.5, fontweight="bold", color=COL["Dihydroquinoline"])
axB.axhline(0, color="k", lw=0.5)
axB.set_xlabel("condition mean yield"); axB.set_ylabel("z(hetero) $-$ z(non-arom.)")
axB.set_title("b   Not a productivity axis", loc="left", fontsize=9, fontweight="bold")
axB.text(0.97, 0.05, f"r = {stats.pearsonr(AX.mean_yield, AX.z_gap)[0]:.3f}", transform=axB.transAxes,
         ha="right", fontsize=7)
o2 = AX.sort_values("z_gap_residual")
axC.bar(range(len(o2)), o2.z_gap_residual, color=[COL[f] for f in o2.reagent_family],
        width=0.8, edgecolor="white", linewidth=0.3)
axC.axhline(0, color="k", lw=0.7)
axC.set_xticks([]); axC.set_xlabel("46 conditions, re-ranked")
axC.set_ylabel("residual after removing yield")
axC.set_title("c   Productivity regressed out", loc="left", fontsize=9, fontweight="bold")
fig.savefig(os.path.join(FIG, "SFigure1.png"), bbox_inches="tight"); plt.close(fig)
AX.to_csv(os.path.join(SRC, "SFigure1_source_data.csv"), index=False)
print("  SFigure 1 written;", len(AX), "conditions")

# ================================================================== SFIGURE 4: quadrangles
QUADS = ["C30,C71|C31,C70", "C53,C57|C68,C70", "C31,C32|C70,C73", "C56,C71|C57,C70",
         "C27,C31|C68,C70", "C51,C54|C55,C58", "C57,C58|C70,C73", "C68,C69|C70,C73"]
def comps(cid):
    r = cA.loc[cid]; out = set()
    for role in ("activator", "additive", "base"):
        c_ = r[f"{role}_name"]
        if pd.isna(c_): continue
        e = r[f"{role}_equiv"]
        out.add((c_, round(float(e), 4) if pd.notna(e) else None))
    return out
SHORT = {"2,4,5-TRICHLOROPHENOL SODIUM SALT": "Cl3-PhONa", "pentafluorophenol": "PFP"}
def factor(x, y):
    a = comps(x) - comps(y); b = comps(y) - comps(x)
    f = lambda st: SHORT.get(next(iter(st))[0], next(iter(st))[0]) if st else "none"
    return f(a) + "/" + f(b)
def quad_label(cells):
    l, r = cells.split("|"); a, b = l.split(","); c, d = r.split(",")
    return factor(a, b) + "  x  " + factor(a, c)
P = ac.pivot_table(index=["sub_2_smiles", "subpair"], columns="condition_id", values="yield", aggfunc="mean")
rng = np.random.default_rng(SEED)
qr = []
for cells in QUADS:
    l, r = cells.split("|"); a, b = l.split(","); c, dd_ = r.split(",")
    sub = P[[a, b, c, dd_]].dropna(); v = sub.values
    dd = (v[:, 0] - v[:, 1]) - (v[:, 2] - v[:, 3])
    acid = sub.index.get_level_values(0).values
    per = pd.Series(dd).groupby(acid).mean()
    st = np.array([1 if STATE[x] == "A" else 0 for x in per.index])
    ua = per.index.values; idx = {u: np.where(acid == u)[0] for u in ua}
    cb = np.array([dd[np.concatenate([idx[u] for u in rng.choice(ua, len(ua), True)])].mean() for _ in range(4000)])
    qr.append({"quadrangle": cells, "n_substrate_pairs": len(sub), "n_acids": len(ua),
               "double_difference": dd.mean(), "ci_lo": np.percentile(cb, 2.5), "ci_hi": np.percentile(cb, 97.5),
               "dd_heteroaromatic": per[st == 1].mean(), "dd_other": per[st == 0].mean(), "label": quad_label(cells)})
Q = pd.DataFrame(qr).sort_values("double_difference").reset_index(drop=True)
fig, (c1, c2) = plt.subplots(1, 2, figsize=(7.6, 3.1), gridspec_kw={"width_ratios": [1.35, 1]})
yy = np.arange(len(Q))
c1.barh(yy, Q.double_difference, color=["#B2182B" if v < 0 else "#2166AC" for v in Q.double_difference],
        xerr=[Q.double_difference - Q.ci_lo, Q.ci_hi - Q.double_difference],
        error_kw={"lw": 0.7}, height=0.65, edgecolor="white", linewidth=0.4)
c1.axvline(0, color="k", lw=0.7); c1.set_yticks(yy); c1.set_yticklabels(Q.label, fontsize=6)
c1.set_xlabel("double difference (yield)")
c1.set_title("a   Components are not simply additive", loc="left", fontsize=9, fontweight="bold")
c2.scatter(Q.dd_other, Q.dd_heteroaromatic, s=34, c="#4393C3", edgecolor="k", linewidth=0.4, zorder=3)
lim = [-0.22, 0.22]
c2.plot(lim, lim, "k--", lw=0.7); c2.axhline(0, color="k", lw=0.5); c2.axvline(0, color="k", lw=0.5)
c2.set_xlim(lim); c2.set_ylim(lim)
c2.set_xlabel("double difference, other acids"); c2.set_ylabel("double difference,\nheteroaromatic acids")
c2.set_title("b   Modulated by state, never\n      opposite in sign", loc="left", fontsize=9, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "SFigure4.png"), bbox_inches="tight"); plt.close(fig)
Q.to_csv(os.path.join(SRC, "SFigure4_source_data.csv"), index=False)
print("  SFigure 4 written;", len(Q), "quadrangles")

# ================================================================== SFIGURE 3: full descriptor block, ECFP4
QUIN = [c for c in keep if cA.loc[c, "activator_name"] in ("EEDQ", "IIDQ")]
FCD = [c for c in keep if cA.loc[c, "activator_name"] in ("TFFH", "BTFFH", "TCFH", "DCC", "DIC", "EDC-HCl")]
ZM = pd.DataFrame(index=acid_smiles, columns=keep, dtype=float)
for c in keep:
    am = ac[ac.condition_id == c].groupby("sub_2_smiles")["yield"].mean()
    ZM.loc[am.index, c] = (am - am.mean()) / am.std(ddof=1)
q = ZM[QUIN].mean(axis=1); f = ZM[FCD].mean(axis=1)
T = pd.DataFrame({"acid_smiles": acid_smiles})
T["acid_state"] = T.acid_smiles.map(STATE)
T["hetero"] = (T.acid_state == "A").astype(int)
T["family_contrast"] = T.acid_smiles.map(q - f)
T["murcko"] = [MurckoScaffold.MurckoScaffoldSmiles(smi) or "ACYCLIC" for smi in T.acid_smiles]
T = T.dropna(subset=["family_contrast"]).reset_index(drop=True)
y = T.family_contrast.values; g = T.murcko.values
X1 = T[["hetero"]].values.astype(float)
def oof_pred(X, yv, splits):
    p = np.empty(len(yv))
    for tr, te in splits: p[te] = Ridge(alpha=1.0, solver="cholesky").fit(X[tr], yv[tr]).predict(X[te])
    return p
sp5 = load_folds(T.acid_smiles, "folds_71acids.csv", groups=T.murcko)
spL = list(LeaveOneGroupOut().split(X1, y, g))
p5 = oof_pred(X1, y, sp5)
r2_5 = 1 - ((y - p5) ** 2).sum() / ((y - y.mean()) ** 2).sum()
r2_L = 1 - ((y - oof_pred(X1, y, spL)) ** 2).sum() / ((y - y.mean()) ** 2).sum()
yb = (y > 0).astype(int)
def oof_auc(X, splits):
    p = np.empty(len(yb))
    for tr, te in splits:
        p[te] = LogisticRegression(max_iter=5000, solver="lbfgs").fit(X[tr], yb[tr]).predict_proba(X[te])[:, 1]
    return roc_auc_score(yb, p), p
auc5, prob5 = oof_auc(X1, sp5); aucL, _ = oof_auc(X1, spL)
sys.path.insert(0, OUT)
import acyl_descriptors as AD
from rdkit.Chem import rdFingerprintGenerator
_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
XF = np.array([np.array(_gen.GetFingerprint(Chem.MolFromSmiles(s))) for s in T.acid_smiles])
DT = AD.descriptor_table(list(T.acid_smiles))
_numcols = [c for c in DT.columns if c.startswith("D") and c[1].isdigit() and DT[c].dtype != object]
_catcols = [c for c in DT.columns if c.startswith("D") and c[1].isdigit() and DT[c].dtype == object and c != "D8_ring_system_frag"]
XD = np.hstack([DT[_numcols].astype(float).values] +
               [pd.get_dummies(DT[c].astype(str), drop_first=True).astype(float).values for c in _catcols])
def scores(X):
    r2 = 1 - ((y - oof_pred(X, y, sp5)) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    a5, _ = oof_auc(X, sp5); aL, _ = oof_auc(X, spL)
    return r2, a5, aL
r2D, aucD5, aucDL = scores(XD)
r2F, aucF5, aucFL = scores(XF)
mods = pd.DataFrame({"model": ["baseline", "two-state\ndescriptor", "structural\nblock (18)", "ECFP4"],
                     "R2_5fold": [0.0, r2_5, r2D, r2F], "AUC_5fold": [0.5, auc5, aucD5, aucF5],
                     "AUC_LOSO": [0.5, aucL, aucDL, aucFL]})
fig, (b1, b2) = plt.subplots(1, 2, figsize=(6.0, 3.0))
xx = np.arange(4)
b1.bar(xx - 0.2, mods.AUC_5fold, 0.38, color="#2166AC", label="5-fold", edgecolor="white", lw=0.4)
b1.bar(xx + 0.2, mods.AUC_LOSO, 0.38, color="#92C5DE", label="leave-one-scaffold-out", edgecolor="white", lw=0.4)
b1.axhline(0.5, color="k", lw=0.8, ls="--")
b1.set_xticks(xx); b1.set_xticklabels(mods.model, fontsize=6.0, rotation=30, ha="right")
b1.set_ylim(0.35, 0.8); b1.set_ylabel("out-of-sample AUC (quinoline vs eight)", fontsize=8.5)
b1.legend(fontsize=6, frameon=False, loc="upper left")
b1.set_title("a   Classification, this target", loc="left", fontsize=9, fontweight="bold")
b2.bar(xx, mods.R2_5fold, color=["#BFBFBF", "#4393C3", "#7F7F7F", "#B2182B"], edgecolor="white", lw=0.4)
b2.axhline(0, color="k", lw=0.7)
b2.set_xticks(xx); b2.set_xticklabels(mods.model, fontsize=6.0, rotation=30, ha="right")
b2.set_ylabel("out-of-fold $R^2$, 5-fold")
b2.set_title("b   Regression, this target", loc="left", fontsize=9, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "SFigure3.png"), bbox_inches="tight"); plt.close(fig)
mods.to_csv(os.path.join(SRC, "SFigure3_source_data.csv"), index=False)
print("  SFigure 3 written; R2 5-fold %.4f  LOSO %.4f  AUC %.4f / %.4f" % (r2_5, r2_L, auc5, aucL))

print("\nfigures ->", FIG)
print("source  ->", SRC)
