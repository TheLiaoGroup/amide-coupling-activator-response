"""Build the four main figures and their Source Data from the self-contained acylation dataset.

Single output root: OUT below. The manuscript package reads from exactly this directory; there is
no second location. Run this, then `build_si_figures.py`, then `REPRODUCE_acylation.py` to verify
the text.

Reads `data/DATASET_acylation_wells.csv` and `data/DATASET_acylation_conditions.xlsx` (produced
from the in-house corpus by `make_dataset.py`), not the multi-transformation corpus.

Figure 1 - data coverage and response overview.
Figure 2 - the raw per-acid paired yield difference, EEDQ/IIDQ against BTFFH (the primary evidence).
Figure 3 - additive effects at matched activator+base, with the acid-class interaction that does
           not survive correction.
Figure 4 - a fair comparison of the two-state descriptor against the carboxyl-carbon charge, same
           split, same target, plus a baseline; the raw pairwise target the charge descriptor wins.
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import LeaveOneGroupOut
RDLogger.DisableLog("rdApp.*")

_here = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("HTE_DATA") or os.path.join(_here, "data")
OUT = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(OUT, "figures")
SRC = os.path.join(OUT, "source_data")
os.makedirs(FIG, exist_ok=True); os.makedirs(SRC, exist_ok=True)
SEED = 20260901

plt.rcParams.update({"font.family": "Times New Roman", "font.size": 8, "axes.linewidth": 0.7,
                     "xtick.major.width": 0.7, "ytick.major.width": 0.7, "mathtext.fontset": "stix",
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 400})

STATE_COL = {"A": "#762A83", "B": "#9A9A9A", "C": "#1B7837"}
LBL = {"A": "Heteroaromatic", "B": "Carbocyclic aryl", "C": "Non-aromatic"}

def flush_yaxis_to_title(ax, fig, renderer=None):
    """A loc="left" title sits at the axes' own left spine (x=0 in axes-fraction coordinates),
    but y-tick/y-axis-label text is drawn further left than that spine - so the visible plot
    (tick numbers included) was overhanging past its own title's start. Move the title left to
    start flush with that content instead of shrinking the axes to meet the title: for a panel
    whose tick labels are wide (e.g. this figure's long reagent-pair names), forcing the plot
    itself over to meet the title would crush the plotted area down to a sliver."""
    if renderer is None:
        renderer = fig.canvas.get_renderer()
    left_title = getattr(ax, "_left_title", None)
    if left_title is None or not left_title.get_text():
        return
    box = ax.get_position()
    fig_w = fig.bbox.width
    spine_x0_disp = fig.transFigure.transform((box.x0, box.y0))[0]
    ybb = ax.yaxis.get_tightbbox(renderer)
    if ybb is None or ybb.x0 >= spine_x0_disp:
        return
    box_w_disp = (box.x1 - box.x0) * fig_w
    left_title.set_x((ybb.x0 - spine_x0_disp) / box_w_disp)

def flush_title_to_title(ax_move, ax_target, fig, renderer=None):
    """Move ax_move's loc="left" title to start at the same x as ax_target's, without touching
    either axes' position - so a full-width panel's own plot content (e.g. an aspect="equal"
    heatmap matplotlib centres in its row) can stay centred while its title still lines up with
    the panel below it."""
    if renderer is None:
        renderer = fig.canvas.get_renderer()
    t_move, t_target = ax_move._left_title, ax_target._left_title
    target_x0_disp = t_target.get_window_extent(renderer=renderer).x0
    box = ax_move.get_position()
    fig_w = fig.bbox.width
    spine_x0_disp = fig.transFigure.transform((box.x0, box.y0))[0]
    box_w_disp = (box.x1 - box.x0) * fig_w
    t_move.set_x((target_x0_disp - spine_x0_disp) / box_w_disp)

# ------------------------------------------------------------------ data and shared functions
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
def q_carbonyl(smi):
    m = Chem.MolFromSmiles(smi); hit = m.GetSubstructMatches(COOH)[0]
    mh = Chem.AddHs(m); AllChem.ComputeGasteigerCharges(mh)
    return float(mh.GetAtomWithIdx(hit[0]).GetDoubleProp("_GasteigerCharge"))
def matched(ac, cx, cy):
    """acid-mean yield at cx minus at cy, over amines matched within each acid"""
    x, y = ac[ac.condition_id == cx], ac[ac.condition_id == cy]; o = {}
    for a in sorted(set(x.sub_2_smiles) & set(y.sub_2_smiles)):
        xa, ya = x[x.sub_2_smiles == a], y[y.sub_2_smiles == a]
        com = set(xa.sub_1_smiles) & set(ya.sub_1_smiles)
        if com:
            o[a] = (xa[xa.sub_1_smiles.isin(com)].groupby("sub_1_smiles")["yield"].mean().mean()
                    - ya[ya.sub_1_smiles.isin(com)].groupby("sub_1_smiles")["yield"].mean().mean())
    return pd.Series(o)

print("loading data ...")
ac = pd.read_csv(os.path.join(DATA, "DATASET_acylation_wells.csv"), low_memory=False).dropna(subset=["yield"])
cond = pd.read_excel(os.path.join(DATA, "DATASET_acylation_conditions.xlsx"))
cA = cond.set_index("condition_id")
acid_smiles = sorted(ac.sub_2_smiles.unique())
amine_smiles = sorted(ac.sub_1_smiles.unique())
STATE = {s: acid_state(s) for s in acid_smiles}
ac["state"] = ac.sub_2_smiles.map(STATE)

# ================================================================== FIGURE 1: data coverage
# figure width is 7.6 in on all four main figures (matching Figure 3, the widest), and every
# equivalent text role uses the same explicit fontsize, so that scaling all four to the same
# display width for the manuscript leaves their text the same size
FIGW = 7.6
fig = plt.figure(figsize=(FIGW, 8.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1.9, 1.0], hspace=0.4, wspace=0.32)
axA = fig.add_subplot(gs[0, :]); axB = fig.add_subplot(gs[1, 0]); axC = fig.add_subplot(gs[1, 1])

acid_i = {a: i for i, a in enumerate(acid_smiles)}
amine_i = {a: i for i, a in enumerate(amine_smiles)}
cov = np.zeros((len(acid_smiles), len(amine_smiles)), dtype=bool)
for a, b in ac[["sub_2_smiles", "sub_1_smiles"]].drop_duplicates().itertuples(index=False):
    cov[acid_i[a], amine_i[b]] = True
# order rows by acid state so the coverage pattern is readable, columns left as encountered
order = sorted(range(len(acid_smiles)), key=lambda i: (STATE[acid_smiles[i]], i))
covo = cov[order, :]
axA.imshow(covo, aspect="equal", cmap=matplotlib.colors.ListedColormap(["#F0F0F0", "#2166AC"]),
           interpolation="none")
axA.set_xlabel(f"{len(amine_smiles)} amines"); axA.set_ylabel(f"{len(acid_smiles)} acids\n(rows grouped by state)")
n_pairs = int(cov.sum()); n_full = len(acid_smiles) * len(amine_smiles)
axA.set_title(f"a   578 measured acid-amine pairs of {n_full} possible "
              f"({100*n_pairs/n_full:.1f}% coverage)", loc="left", fontsize=9, fontweight="bold")
# state boundaries
bounds = np.cumsum([sum(1 for i in order if STATE[acid_smiles[i]] == s) for s in "ABC"])
for b in bounds[:-1]:
    axA.axhline(b - 0.5, color="white", lw=1.2)
axA.set_yticks([])

fam_counts = pd.Series({s: sum(1 for x in acid_smiles if STATE[x] == s) for s in "ABC"})
axB.bar(range(3), [fam_counts[s] for s in "ABC"], color=[STATE_COL[s] for s in "ABC"],
        edgecolor="white", linewidth=0.4)
axB.set_xticks(range(3)); axB.set_xticklabels([LBL[s] for s in "ABC"], fontsize=7, rotation=15)
axB.set_ylabel("acids"); axB.set_title("b   Acid classes", loc="left", fontsize=9, fontweight="bold")
for i, s in enumerate("ABC"):
    axB.text(i, fam_counts[s] + 0.5, str(fam_counts[s]), ha="center", fontsize=7)

axC.hist(ac["yield"], bins=60, color="#4393C3", edgecolor="none")
axC.set_xlabel("yield"); axC.set_ylabel("wells")
axC.set_title("c   53,233 calibrated yields", loc="left", fontsize=9, fontweight="bold")
axC.text(0.97, 0.92, f"median = {ac['yield'].median():.3f}", transform=axC.transAxes, ha="right", fontsize=6.5)

# panel a's imshow is aspect="equal" (square cells), so matplotlib centres its narrower-than-cell
# plot within the full-width row - left as is, this is the right look for the top full-width panel.
# panels b and c each get their own tick labels flush with their own title first; panel a's title
# then moves (its plot does not) to line up with panel b's, now-final, title position.
fig.canvas.draw()
flush_yaxis_to_title(axB, fig)
flush_yaxis_to_title(axC, fig)
fig.canvas.draw()
flush_title_to_title(axA, axB, fig)

fig.savefig(os.path.join(FIG, "Figure1.png"), bbox_inches="tight"); plt.close(fig)
COV = pd.DataFrame({"acid_smiles": [a for a, _ in ac[["sub_2_smiles","sub_1_smiles"]].drop_duplicates().itertuples(index=False)],
                     "amine_smiles": [b for _, b in ac[["sub_2_smiles","sub_1_smiles"]].drop_duplicates().itertuples(index=False)]})
COV["acid_state"] = COV.acid_smiles.map(STATE)
COV.to_csv(os.path.join(SRC, "Figure1_source_data.csv"), index=False)
print(f"  Figure 1 written; {n_pairs} pairs, {len(acid_smiles)} acids, {len(amine_smiles)} amines")

# ================================================================== FIGURE 2: raw paired differences
d_eedq = matched(ac, "C41", "C3")   # EEDQ minus BTFFH
d_iidq = matched(ac, "C95", "C3")   # IIDQ minus BTFFH
F2 = pd.DataFrame({"acid_smiles": sorted(set(d_eedq.index) | set(d_iidq.index))})
F2["acid_state"] = F2.acid_smiles.map(STATE)
F2["EEDQ_minus_BTFFH"] = F2.acid_smiles.map(d_eedq)
F2["IIDQ_minus_BTFFH"] = F2.acid_smiles.map(d_iidq)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(FIGW, 3.4), sharey=True)
for ax, col, lab, quin in [(a1, "EEDQ_minus_BTFFH", "EEDQ $-$ BTFFH", "EEDQ"),
                           (a2, "IIDQ_minus_BTFFH", "IIDQ $-$ BTFFH", "IIDQ")]:
    sub = F2.dropna(subset=[col]).sort_values(col).reset_index(drop=True)
    ax.bar(range(len(sub)), sub[col], color=[STATE_COL[s] for s in sub.acid_state],
           width=0.85, edgecolor="none")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks([]); ax.set_xlabel(f"{len(sub)} acids, sorted")
    ax.set_title(lab, loc="left", fontsize=9, fontweight="bold")
    frac_pos = float((sub[col] > 0).mean())
    ax.text(0.03, 0.94, f"{100*frac_pos:.0f}% of acids favour {quin} (a dihydroquinoline)",
            transform=ax.transAxes, fontsize=6.5, va="top")
a1.set_ylabel("per-acid yield difference (amine-matched)")
a1.set_title("a   " + a1.get_title(loc="left"), loc="left", fontsize=9, fontweight="bold")
a2.set_title("b   " + a2.get_title(loc="left"), loc="left", fontsize=9, fontweight="bold")
a1.legend(handles=[Patch(facecolor=STATE_COL[s], label=LBL[s]) for s in "ABC"],
          fontsize=6.5, frameon=False, loc="lower right")
fig.tight_layout()
fig.canvas.draw()
flush_yaxis_to_title(a1, fig)
flush_yaxis_to_title(a2, fig)
fig.savefig(os.path.join(FIG, "Figure2.png"), bbox_inches="tight"); plt.close(fig)
F2.to_csv(os.path.join(SRC, "Figure2_source_data.csv"), index=False)
print(f"  Figure 2 written; {len(F2)} acids")

# ================================================================== FIGURE 3: additive effects
def recipe_key(row):
    return (row.activator_CAS, round(float(row.activator_equiv), 3),
            row.base_CAS if pd.notna(row.base_CAS) else None,
            round(float(row.base_equiv), 3) if pd.notna(row.base_equiv) else None)
cA2 = cA.copy(); cA2["_recipe"] = [recipe_key(r) for _, r in cA2.iterrows()]
add_pairs = []
for _, grp in cA2.groupby("_recipe"):
    no, wi = grp[grp.additive_CAS.isna()], grp[grp.additive_CAS.notna()]
    if len(no) and len(wi):
        for nid in no.index:
            for wid in wi.index:
                add_pairs.append((nid, wid))

rng3 = np.random.default_rng(SEED)
rng3p = np.random.default_rng(20260912)  # matches REPRODUCE_acylation.py section 16 exactly,
# so the exported p/q agree with the text to the same Monte Carlo draw, not just in conclusion
add_rows = []
for nid, wid in add_pairs:
    d = -matched(ac, nid, wid)  # with-additive minus no-additive
    st = np.array([STATE[a] for a in d.index]); isA = st == "A"
    boot = np.array([rng3.choice(d.values, len(d), replace=True).mean() for _ in range(4000)])
    row = {"no_additive": nid, "with_additive": wid, "activator": cA.loc[nid, "activator_name"],
           "additive": cA.loc[wid, "additive_name"], "n_acids": len(d),
           "mean_diff_pp": 100 * float(d.mean()),
           "ci_lo_pp": 100 * float(np.percentile(boot, 2.5)), "ci_hi_pp": 100 * float(np.percentile(boot, 97.5))}
    if isA.sum() >= 3 and (~isA).sum() >= 3:
        obs = float(d.values[isA].mean() - d.values[~isA].mean())
        row["interaction_pp"] = 100 * obs
        k = int(isA.sum())
        null = np.array([(lambda q: q[:k].mean() - q[k:].mean())(rng3p.permutation(d.values))
                          for _ in range(20000)])
        row["p"] = (np.sum(np.abs(null) >= abs(obs)) + 1) / 20001
    add_rows.append(row)
ADD = pd.DataFrame(add_rows).sort_values("mean_diff_pp").reset_index(drop=True)
ADD["label"] = ADD.activator + " + " + ADD.additive.fillna("none")
HEADLINE = {("C77", "C78"), ("C87", "C88"), ("C47", "C60"), ("C20", "C19")}
ADD["headline"] = [(n, w) in HEADLINE for n, w in zip(ADD.no_additive, ADD.with_additive)]
# Benjamini-Hochberg q-values across the 21 interaction p-values, exported for auditability
# (the plot's "0 of 21, q<0.05" title text is independently checked against this by
# check_documents.py, rather than being a plain string no one can look behind)
_pmask = ADD.p.notna()
_pvals = ADD.loc[_pmask, "p"].values
_order = np.argsort(_pvals)
_m = len(_pvals)
_q_sorted = np.empty(_m)
_ps = _pvals[_order]
_q_sorted[-1] = _ps[-1]
for _i in range(_m - 2, -1, -1):
    _q_sorted[_i] = min(_q_sorted[_i + 1], _ps[_i] * _m / (_i + 1))
_q = np.empty(_m); _q[_order] = _q_sorted
ADD.loc[_pmask, "q_bh"] = _q
n_sig = int((ADD.q_bh < 0.05).sum())
print(f"  additive x acid-class interactions surviving BH correction at q<0.05: {n_sig} of {_m}")

fig, (c1, c2) = plt.subplots(1, 2, figsize=(FIGW, 4.4), gridspec_kw={"width_ratios": [1.4, 1]})
yy = np.arange(len(ADD))
colors = ["#B2182B" if h else "#9DB8D2" for h in ADD.headline]
c1.barh(yy, ADD.mean_diff_pp, xerr=[ADD.mean_diff_pp - ADD.ci_lo_pp, ADD.ci_hi_pp - ADD.mean_diff_pp],
        color=colors, error_kw={"lw": 0.6, "ecolor": "0.3"}, height=0.7, edgecolor="white", linewidth=0.3)
for y, v, hi, lo in zip(yy, ADD.mean_diff_pp, ADD.ci_hi_pp, ADD.ci_lo_pp):
    if v >= 0:
        c1.text(hi + 0.8, y, f"{v:+.1f}", ha="left", va="center", fontsize=6.2)
    elif v < -3:
        # placed just inside the bar's own tip (white-on-red), not to its left, because this one
        # strongly negative bar (HATU + HOAT) reaches close enough to the axis edge that an
        # outside label there collides with that row's own category name
        c1.text(v + 0.5, y, f"{v:+.1f}", ha="left", va="center", fontsize=6.2, color="white")
    else:
        c1.text(lo - 0.8, y, f"{v:+.1f}", ha="right", va="center", fontsize=6.2)
c1.axvline(0, color="k", lw=0.7)
c1.set_yticks(yy); c1.set_yticklabels(ADD.label, fontsize=6.2)
c1.set_xlabel("additive effect, acid-equal-weighted (percentage points of yield)")
c1.set_title("a   21 matched additive comparisons", loc="left", fontsize=9, fontweight="bold")
c1.legend(handles=[Patch(facecolor="#B2182B", label="named in text"),
                    Patch(facecolor="#9DB8D2", label="other 17")],
          fontsize=6.5, frameon=False, loc="lower right")

interax = ADD.dropna(subset=["interaction_pp"]).copy()
c2.scatter(interax.mean_diff_pp, interax.interaction_pp, s=30,
           color=["#B2182B" if h else "#4393C3" for h in interax.headline], edgecolor="k", linewidth=0.3)
c2.axhline(0, color="k", lw=0.6, ls=":")
c2.set_xlabel("overall additive effect (pp)"); c2.set_ylabel("heteroaromatic $-$ other acids (pp)")
c2.set_title("b   No comparison survives correction", loc="left", fontsize=9, fontweight="bold")
c2.text(0.97, 0.04, f"{n_sig} of {_m} at $q$ < 0.05\n(Benjamini–Hochberg)", transform=c2.transAxes,
        ha="right", va="bottom", fontsize=6.5)
fig.tight_layout()
fig.canvas.draw()
flush_yaxis_to_title(c1, fig)
flush_yaxis_to_title(c2, fig)
fig.savefig(os.path.join(FIG, "Figure3.png"), bbox_inches="tight"); plt.close(fig)
ADD.to_csv(os.path.join(SRC, "Figure3_source_data.csv"), index=False)
print(f"  Figure 3 written; {len(ADD)} additive comparisons")

# ================================================================== FIGURE 4: fair descriptor comparison
targets = {"EEDQ - BTFFH": d_eedq, "IIDQ - BTFFH": d_iidq}
rows4, scatter_rows = [], []
for label, tgt in targets.items():
    tgt = tgt.dropna()
    acidsT = list(tgt.index); yv = tgt.values
    st = np.array([1.0 if STATE[a] == "A" else 0.0 for a in acidsT])
    chg = np.array([q_carbonyl(a) for a in acidsT])
    murcko = np.array([MurckoScaffold.MurckoScaffoldSmiles(a) or "ACYCLIC" for a in acidsT])
    splits = list(LeaveOneGroupOut().split(st.reshape(-1, 1), yv, murcko))
    def ols_oof(*xcols):
        Xd = np.hstack([np.ones((len(yv), 1))] + [c.reshape(-1, 1) for c in xcols])
        p = np.empty(len(yv))
        for tr, te in splits:
            b = np.linalg.pinv(Xd[tr]) @ yv[tr]; p[te] = Xd[te] @ b
        return p, 1 - ((yv - p) ** 2).sum() / ((yv - yv.mean()) ** 2).sum()
    def baseline_oof():
        p = np.empty(len(yv))
        for tr, te in splits: p[te] = yv[tr].mean()
        return 1 - ((yv - p) ** 2).sum() / ((yv - yv.mean()) ** 2).sum()
    p_s, r2s = ols_oof(st)
    p_c, r2c = ols_oof(chg)
    _, r2j = ols_oof(st, chg)  # does the structural descriptor add anything once charge is in
    # the model, under the same leave-one-scaffold-out protocol as every other number here?
    r2b = baseline_oof()
    rows4.append({"target": label, "model": "baseline (train-fold mean)", "R2": r2b})
    rows4.append({"target": label, "model": "two-state descriptor", "R2": r2s})
    rows4.append({"target": label, "model": "carboxyl-carbon charge", "R2": r2c})
    rows4.append({"target": label, "model": "structural + charge (joint)", "R2": r2j})
    for a, obs, pred in zip(acidsT, yv, p_c):
        scatter_rows.append({"target": label, "acid_smiles": a, "acid_state": STATE[a],
                              "observed": obs, "oof_predicted_charge": pred})

MODS = pd.DataFrame(rows4)
SCAT = pd.DataFrame(scatter_rows)
fig, (d1, d2) = plt.subplots(1, 2, figsize=(FIGW, 3.2), gridspec_kw={"width_ratios": [1.15, 1]})
targs = list(targets.keys())
xw = 0.19
_model_order = ["baseline (train-fold mean)", "two-state descriptor", "carboxyl-carbon charge",
                 "structural + charge (joint)"]
_model_colors = ["#BFBFBF", "#4393C3", "#B2182B", "#35978F"]
for j, m in enumerate(_model_order):
    xs = np.arange(2) + (j - 1.5) * xw
    vals = [MODS[(MODS.target == t) & (MODS.model == m)].R2.iloc[0] for t in targs]
    d1.bar(xs, vals, xw, color=_model_colors[j], edgecolor="white", linewidth=0.4,
           label=m.split(" (")[0])
    for x, v in zip(xs, vals):
        d1.text(x, v + (0.008 if v >= 0 else -0.008), f"{v:.3f}", ha="center",
                 va="bottom" if v >= 0 else "top", fontsize=6.5)
d1.axhline(0, color="k", lw=0.6)
d1.set_xticks(range(2)); d1.set_xticklabels(targs)
d1.set_ylabel("out-of-fold $R^2$ (leave-one-scaffold-out)")
d1.set_ylim(bottom=-0.075, top=0.40)
d1.legend(fontsize=6, frameon=False, loc="upper left")
d1.set_title("a   Same split, same target", loc="left", fontsize=9, fontweight="bold")

s0 = SCAT[SCAT.target == "EEDQ - BTFFH"]
d2.scatter(s0.observed, s0.oof_predicted_charge, s=24,
           c=[STATE_COL[s] for s in s0.acid_state], edgecolor="k", linewidth=0.25)
lim = [s0.observed.min() - 0.03, s0.observed.max() + 0.03]
d2.plot(lim, lim, "k--", lw=0.7); d2.set_xlim(lim)
d2.set_xlabel("observed"); d2.set_ylabel("out-of-fold prediction (charge)")
r2_disp = MODS[(MODS.target == "EEDQ - BTFFH") & (MODS.model == "carboxyl-carbon charge")].R2.iloc[0]
d2.text(0.04, 0.92, f"$R^2$ = {r2_disp:.3f}", transform=d2.transAxes, fontsize=6.5, va="top")
d2.text(0.55, 0.94, "EEDQ $-$ BTFFH", transform=d2.transAxes, fontsize=6.5, ha="center", va="top")
d2.legend(handles=[plt.Line2D([0], [0], marker="o", linestyle="", markerfacecolor=STATE_COL[s],
                               markeredgecolor="k", markeredgewidth=0.25, markersize=6, label=LBL[s])
                    for s in "ABC"], fontsize=6.5, frameon=False, loc="lower right")
d2.set_title("b   Carboxyl-charge, held-out scaffolds", loc="left", fontsize=9, fontweight="bold")
fig.tight_layout()
fig.canvas.draw()
flush_yaxis_to_title(d1, fig)
flush_yaxis_to_title(d2, fig)
fig.savefig(os.path.join(FIG, "Figure4.png"), bbox_inches="tight"); plt.close(fig)
MODS.to_csv(os.path.join(SRC, "Figure4_source_data.csv"), index=False)
SCAT.to_csv(os.path.join(SRC, "Figure4_scatter_data.csv"), index=False)
print("  Figure 4 written;", dict(zip(MODS.model.unique(), [MODS[MODS.model==m].R2.round(3).tolist() for m in MODS.model.unique()])))

print("\nfigures ->", FIG)
print("source  ->", SRC)
