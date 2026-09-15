"""Recompute every number in MANUSCRIPT_acylation.md and SI_acylation.md from DATASET_acylation_*.

This script is an INDEPENDENT reimplementation. It imports nothing from `build_figures.py` or
`acyl_descriptors.py`, and reads no Source Data CSV - only the two self-contained acylation dataset
files (`data/DATASET_acylation_wells.csv` and `data/DATASET_acylation_conditions.xlsx`, produced
from the released corpus by `make_dataset.py`) and the frozen five-fold scaffold partition
(`frozen_folds/`, checked-in data, not output from any script in this repository - see its
README.md for why the split is frozen rather than recomputed). A number that agrees with its own
generator is not reproduced.

Every check asserts that it matched a non-empty set of rows before comparing values, so a check
that silently selects nothing fails instead of passing.

    <python> REPRODUCE_acylation.py            all checks
    <python> REPRODUCE_acylation.py --quick    skip the permutation checks

Verified under Python 3.9.10 with pandas 2.3.3, numpy 2.0.2, scipy 1.13.1, scikit-learn 1.6.1,
rdkit 2025.09.2.
"""
import os, sys, hashlib
import numpy as np, pandas as pd
from scipy import stats
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.metrics import roc_auc_score
RDLogger.DisableLog("rdApp.*")

_here = os.path.dirname(os.path.abspath(__file__))
# the two self-contained files live in ./data beside this script; HTE_DATA overrides that
# location for a data/ folder placed somewhere else after unpacking the archive
DATA = os.environ.get("HTE_DATA") or os.path.join(_here, "data")
QUICK = "--quick" in sys.argv
SEED = 20260901
RESULTS = []
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


def check(name, got, want, tol=None, lo=None, hi=None):
    """Record one check. Either an exact/tolerance match, or a band."""
    if lo is not None or hi is not None:
        ok = (lo is None or got >= lo) and (hi is None or got <= hi)
        exp = f"[{lo}, {hi}]"
    elif tol is None:
        ok = got == want
        exp = str(want)
    else:
        ok = abs(got - want) <= tol
        exp = f"{want} +/- {tol}"
    RESULTS.append((ok, name, got, exp))
    print(("  PASS  " if ok else "  FAIL  ") + f"{name:<62s} got {got!s:<24s} expected {exp}")
    return ok


def guard(n, what):
    """Rule: a check that matched nothing must fail loudly, not pass quietly."""
    if n <= 0:
        RESULTS.append((False, f"EMPTY SELECTION: {what}", 0, "> 0"))
        print(f"  FAIL  EMPTY SELECTION: {what}")
        raise SystemExit("a selection matched zero rows; aborting rather than reporting success")
    return n


# ---------------------------------------------------------------- load
print("\n== loading released data ==")
raw = pd.read_csv(os.path.join(DATA, "DATASET_acylation_wells.csv"), low_memory=False)
cond = pd.read_excel(os.path.join(DATA, "DATASET_acylation_conditions.xlsx"))
guard(len(raw), "Acylation|1 rows")
ac = raw.dropna(subset=["yield"]).copy()

print("\n== 1. the array ==")
check("wells in Acylation|1 (before yield filter)", len(raw), 53241)
check("wells with a missing yield", int(raw["yield"].isna().sum()), 8)
check("substrate pairs", raw.subpair.nunique(), 578)
check("carboxylic acids", raw.sub_2_smiles.nunique(), 71)
check("amines", raw.sub_1_smiles.nunique(), 82)
check("conditions", raw.condition_id.nunique(), 94)
check("plates", raw.plate.nunique(), 603)
check("all conditions DMF / 25 C / 2 h", bool((cond.solvent_name == "DMF").all()
      and (cond.temp_C == 25.0).all() and (cond.time_h == 2.0).all()), True)

# side identity - which partner is the acid
u = raw.drop_duplicates("subpair")
guard(len(u), "unique substrate pairs for the side check")
check("fraction of pairs with carboxylic acid on sub_2",
      round(float(u.sub_2_FG.fillna("").str.lower().str.contains("carboxylic acid").mean()), 3), 0.965)
check("fraction of pairs with carboxylic acid on sub_1",
      round(float(u.sub_1_FG.fillna("").str.lower().str.contains("carboxylic acid").mean()), 3), 0.000)

cA = cond.set_index("condition_id")
guard(len(cA), "Acylation conditions in the condition table")
act_of = cA["activator_name"].to_dict()
check("distinct activators across the 94 conditions",
      len({act_of[c] for c in raw.condition_id.unique()}), 41)

# ---------------------------------------------------------------- acid state
print("\n== 2. acid state assignment ==")
_cooh = Chem.MolFromSmarts("C(=O)[OX2H1]")


def state_of(smi):
    """heteroaromatic / carbocyclic aryl / non-aromatic, by fused aromatic ring system.
    Independent implementation: ring systems are built by union-find over shared atoms."""
    m = Chem.MolFromSmiles(smi)
    hit = m.GetSubstructMatches(_cooh)
    if not hit:
        return "?"
    cidx, oidx, ohidx = hit[0]
    nb = [a for a in m.GetAtomWithIdx(cidx).GetNeighbors() if a.GetIdx() not in (oidx, ohidx)]
    if not nb:
        return "C"
    at = nb[0]
    if not at.GetIsAromatic():
        return "C"
    arom = [set(r) for r in m.GetRingInfo().AtomRings()
            if all(m.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    parent = list(range(len(arom)))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    for i in range(len(arom)):
        for j in range(i + 1, len(arom)):
            if arom[i] & arom[j]:
                parent[find(i)] = find(j)
    grp = [i for i in range(len(arom)) if at.GetIdx() in arom[i]]
    if not grp:
        return "C"
    root = find(grp[0])
    sysm = set().union(*[arom[i] for i in range(len(arom)) if find(i) == root])
    return "A" if any(m.GetAtomWithIdx(i).GetSymbol() != "C" for i in sysm) else "B"


acids = sorted(raw.sub_2_smiles.unique())
guard(len(acids), "unique acids")
STATE = {s: state_of(s) for s in acids}
vc = pd.Series(list(STATE.values())).value_counts()
check("heteroaromatic acids", int(vc.get("A", 0)), 21)
check("carbocyclic aryl acids", int(vc.get("B", 0)), 18)
check("non-aromatic acids", int(vc.get("C", 0)), 32)
ac["state"] = ac.sub_2_smiles.map(STATE)

# ---------------------------------------------------------------- the axis
print("\n== 3. the activator axis ==")
keep = []
for c in sorted(ac.condition_id.unique()):
    g = ac[ac.condition_id == c]
    if pd.notna(cA.loc[c, "additive_CAS"]):
        continue
    if len(g) <= 300:
        continue
    keep.append(c)
guard(len(keep), "conditions passing the axis filter")
check("conditions on the axis", len(keep), 46)

rows = []
for c in keep:
    g = ac[ac.condition_id == c]
    am = g.groupby("sub_2_smiles")["yield"].mean()
    st = am.index.map(STATE)
    z = (am - am.mean()) / am.std(ddof=1)
    rows.append({"cid": c, "activator": act_of[c], "n": len(g), "meanY": g["yield"].mean(),
                 "zg": z[st == "A"].mean() - z[st == "C"].mean(),
                 "rawg": am[st == "A"].mean() - am[st == "C"].mean()})
AX = pd.DataFrame(rows).set_index("cid")
check("smallest well count among the 46", int(AX.n.min()), 539)
check("largest well count among the 46", int(AX.n.max()), 577)
r_raw, p_raw = stats.pearsonr(AX.rawg, AX.meanY)
check("r(raw class gap, condition mean yield)", round(float(r_raw), 3), 0.430, tol=0.002)
check("p for that correlation", round(float(p_raw), 4), 0.0029, tol=0.0005)
check("r(standardised axis, raw class gap)",
      round(float(stats.pearsonr(AX.zg, AX.rawg)[0]), 3), 0.896, tol=0.002)
check("r(axis, condition mean yield)",
      round(float(stats.pearsonr(AX.zg, AX.meanY)[0]), 3), 0.566, tol=0.002)
check("axis minimum", round(float(AX.zg.min()), 2), -0.91, tol=0.005)
check("axis maximum", round(float(AX.zg.max()), 2), 0.77, tol=0.005)

FCD = ("TFFH", "BTFFH", "TCFH", "DCC", "DIC", "EDC-HCl")
fcd = AX[AX.activator.isin(FCD)]
qn = AX[AX.activator.isin(("EEDQ", "IIDQ"))]
guard(len(fcd), "fluoroformamidinium/carbodiimide conditions")
guard(len(qn), "quinoline conditions")
check("fluoroformamidinium + carbodiimide conditions", len(fcd), 8)
check("their maximum (whole set below -0.40)", round(float(fcd.zg.max()), 2), -0.40, tol=0.01)
check("their minimum", round(float(fcd.zg.min()), 2), -0.91, tol=0.01)
check("quinoline conditions", len(qn), 2)
check("quinoline minimum (whole set above +0.74)", round(float(qn.zg.min()), 2), 0.74, tol=0.01)
check("no overlap between the two families", bool(fcd.zg.max() < qn.zg.min()), True)
# self-referential claim in the text: "nothing else in the array lies above the quinolines"
check("conditions above the lower quinoline", int((AX.zg > qn.zg.min()).sum()), 1)
# self-referential claim: the largest scaffold group holds 25 of the 71 acids
_sc = pd.Series([MurckoScaffold.MurckoScaffoldSmiles(a) or "ACYCLIC" for a in acids])
check("Bemis-Murcko scaffold groups", int(_sc.nunique()), 32)
check("largest scaffold group", int(_sc.value_counts().max()), 25)
for a, want in [("BTFFH", [-0.83, -0.68]), ("TFFH", [-0.91, -0.63]),
                ("DCMT", [-0.64, -0.58]), ("EDC-HCl", [-0.51, -0.46])]:
    got = sorted(round(float(v), 2) for v in AX[AX.activator == a].zg)
    guard(len(got), f"{a} conditions")
    check(f"{a} both conditions", got, sorted(want))
for a, want in [("EEDQ", 0.77), ("IIDQ", 0.74), ("CDMT", 0.58), ("DEPBT", 0.46), ("PyOxim", 0.30)]:
    v = AX[AX.activator == a].zg
    guard(len(v), f"{a} condition")
    check(f"{a} axis coordinate", round(float(v.iloc[0]), 2), want, tol=0.006)
for a, want in [("EEDQ", 0.139), ("IIDQ", 0.138), ("HATU", 0.49), ("HCTU", 0.48)]:
    v = AX[AX.activator == a].meanY
    guard(len(v), f"{a} mean yield")
    check(f"{a} mean yield", round(float(v.max()), 3 if want < 0.2 else 2), want, tol=0.006)
for a, want in [("DMTMM", 0.008), ("IBCF", 0.012)]:
    check(f"{a} mean yield", round(float(AX[AX.activator == a].meanY.iloc[0]), 3), want, tol=0.001)
dc = sorted(round(float(v), 3) for v in AX[AX.activator == "DCMT"].meanY)
check("DCMT mean yields", dc, [0.009, 0.014])

# residualised axis
sl, ic = stats.linregress(AX.meanY, AX.zg)[:2]
AX["res"] = AX.zg - (ic + sl * AX.meanY)
order = list(AX.sort_values("res", ascending=False).index)
for a, rank, val in [("EEDQ", 1, 1.03), ("IIDQ", 2, 1.01), ("CDMT", 4, 0.52), ("DEPBT", 5, 0.44)]:
    cid = AX[AX.activator == a].index[0]
    check(f"{a} rank on the residualised axis", order.index(cid) + 1, rank)
    check(f"{a} residual value", round(float(AX.loc[cid, "res"]), 2), val, tol=0.006)

# ---------------------------------------------------------------- permutations
zmat = {}
for c in keep:
    am = ac[ac.condition_id == c].groupby("sub_2_smiles")["yield"].mean()
    zmat[c] = (am - am.mean()) / am.std(ddof=1)
lab = np.array([STATE[a] for a in acids])
Z = np.full((len(keep), len(acids)), np.nan)
ai = {a: i for i, a in enumerate(acids)}
for k, c in enumerate(keep):
    Z[k, [ai[a] for a in zmat[c].index]] = zmat[c].values


def gaps(mat, labels, rows_idx):
    A = np.nanmean(np.where(labels == "A", mat[rows_idx], np.nan), axis=1)
    C = np.nanmean(np.where(labels == "C", mat[rows_idx], np.nan), axis=1)
    return A - C


if not QUICK:
    print("\n== 4. interaction significance (permutation) ==")
    rng = np.random.default_rng(SEED)
    allrows = list(range(len(keep)))
    obs = float(np.nanstd(gaps(Z, lab, allrows), ddof=1))
    check("SD of the axis across the 46 conditions", round(obs, 3), 0.396, tol=0.002)
    null = np.array([np.nanstd(gaps(Z, rng.permutation(lab), allrows), ddof=1) for _ in range(20000)])
    check("permutation null median", round(float(np.median(null)), 3), 0.202, tol=0.006)
    p = (np.sum(null >= obs) + 1) / 20001
    check("p for the interaction", float(p), None, hi=0.001)

    print("\n== 5. yield floor sensitivity ==")
    _FLOOR_WANT = {0.01: (0.398, 0.195, 0.551), 0.02: (0.401, 0.188, 0.543),
                   0.05: (0.381, 0.179, 0.450), 0.10: (0.399, 0.156, 0.440)}
    for floor in [0.01, 0.02, 0.05, 0.10]:
        sel = AX[AX.meanY >= floor]
        idx = [keep.index(c) for c in sel.index]
        guard(len(idx), f"conditions above a {floor} yield floor")
        o = float(np.nanstd(gaps(Z, lab, idx), ddof=1))
        nl = np.array([np.nanstd(gaps(Z, rng.permutation(lab), idx), ddof=1) for _ in range(10000)])
        check(f"floor {floor:.2f}: p", float((np.sum(nl >= o) + 1) / 10001), None, hi=0.0005)
        s2 = sel.sort_values("zg", ascending=False)
        top2 = [act_of[c] for c in s2.index[:2]]
        check(f"floor {floor:.2f}: top two activators", sorted(top2), ["EEDQ", "IIDQ"])
        _w_disp, _w_null, _w_r = _FLOOR_WANT[floor]
        check(f"floor {floor:.2f}: axis dispersion", round(o, 3), _w_disp, tol=0.002)
        check(f"floor {floor:.2f}: permuted null median", round(float(np.median(nl)), 3), _w_null, tol=0.006)
        _r_floor = float(stats.pearsonr(sel.zg, sel.meanY)[0])
        check(f"floor {floor:.2f}: r(axis, yield)", round(_r_floor, 3), _w_r, tol=0.006)

    print("\n== 6. cross-family support (quinolines deleted) ==")
    red = AX.drop(index=[AX[AX.activator == a].index[0] for a in ("EEDQ", "IIDQ")])
    check("conditions on the reduced array", len(red), 44)
    check("r(axis, yield) on the reduced array",
          round(float(stats.pearsonr(red.zg, red.meanY)[0]), 3), 0.716, tol=0.003)
    idx = [keep.index(c) for c in red.index]
    o = float(np.nanstd(gaps(Z, lab, idx), ddof=1))
    check("SD on the reduced array", round(o, 3), 0.353, tol=0.003)
    nl = np.array([np.nanstd(gaps(Z, rng.permutation(lab), idx), ddof=1) for _ in range(20000)])
    check("null median on the reduced array", round(float(np.median(nl)), 3), 0.199, tol=0.006)
    check("p on the reduced array", round(float((np.sum(nl >= o) + 1) / 20001), 3), 0.002, tol=0.001)
    for a, want, band in [("CDMT", 0.039, 0.02), ("DEPBT", 0.105, 0.03), ("PyOxim", 0.290, 0.04)]:
        cid = AX[AX.activator == a].index[0]
        k = keep.index(cid)
        ob = float(gaps(Z, lab, [k])[0])
        nl2 = np.array([gaps(Z, rng.permutation(lab), [k])[0] for _ in range(20000)])
        pv = (np.sum(np.abs(nl2) >= abs(ob)) + 1) / 20001
        check(f"{a} alone: p", round(float(pv), 3), want, tol=band)
    trio = [keep.index(AX[AX.activator == a].index[0]) for a in ("CDMT", "DEPBT", "PyOxim")]
    eight = [keep.index(c) for c in fcd.index]
    contrast = np.nanmean(Z[trio], axis=0) - np.nanmean(Z[eight], axis=0)
    for s, want in [("A", 0.66), ("B", -0.07), ("C", -0.39)]:
        check(f"trio-minus-eight contrast, {s}",
              round(float(np.nanmean(contrast[lab == s])), 2), want, tol=0.006)
    obs2 = float(np.nanmean(contrast[lab == "A"]) - np.nanmean(contrast[lab == "C"]))
    check("trio-minus-eight, heteroaromatic minus non-aromatic", round(obs2, 2), 1.05, tol=0.006)
    nl3 = np.empty(20000)
    for i in range(20000):
        pp = rng.permutation(lab)
        nl3[i] = np.nanmean(contrast[pp == "A"]) - np.nanmean(contrast[pp == "C"])
    check("its p", float((np.sum(np.abs(nl3) >= abs(obs2)) + 1) / 20001), None, hi=0.001)
    bs = np.array([np.nanmean(rng.choice(contrast[lab == "A"], (lab == "A").sum(), True))
                   - np.nanmean(rng.choice(contrast[lab == "C"], (lab == "C").sum(), True))
                   for _ in range(20000)])
    check("its CI low", round(float(np.percentile(bs, 2.5)), 2), 0.63, tol=0.03)
    check("its CI high", round(float(np.percentile(bs, 97.5)), 2), 1.45, tol=0.03)
    slr, icr = stats.linregress(red.meanY, red.zg)[:2]
    red = red.assign(res2=red.zg - (icr + slr * red.meanY))
    ro = list(red.sort_values("res2", ascending=False).index)
    for a, rank in [("CDMT", 2), ("DEPBT", 3), ("PyOxim", 10)]:
        check(f"{a} rank on the residualised reduced array",
              ro.index(AX[AX.activator == a].index[0]) + 1, rank)

# ---------------------------------------------------------------- Figure 3
print("\n== 7. the three single-variable conditions ==")
sets = {c: ac[ac.condition_id == c] for c in ("C3", "C41", "C95")}
for c, g in sets.items():
    guard(len(g), f"wells in {c}")
fl = []
for a in sorted(set.intersection(*[set(g.sub_2_smiles) for g in sets.values()])):
    per = {c: g[g.sub_2_smiles == a] for c, g in sets.items()}
    com = set.intersection(*[set(p.sub_1_smiles) for p in per.values()])
    if not com:
        continue
    r = {"state": STATE[a]}
    for c, p in per.items():
        r[c] = p[p.sub_1_smiles.isin(com)].groupby("sub_1_smiles")["yield"].mean().mean()
    fl.append(r)
FL = pd.DataFrame(fl)
guard(len(FL), "acids in the three-way matched set")
mm = FL.groupby("state")[["C3", "C41", "C95"]].mean()
for s, want in [("A", [0.089, 0.189, 0.170]), ("B", [0.114, 0.100, 0.101]),
                ("C", [0.183, 0.081, 0.091])]:
    check(f"absolute yields, state {s}",
          [round(float(mm.loc[s, c]), 3) for c in ("C3", "C41", "C95")], want)
for s, want in [("A", 21), ("B", 18), ("C", 30)]:
    check(f"acids in the matched set, state {s}", int((FL.state == s).sum()), want)

# the text's class-mean sentence ("EEDQ raises the heteroaromatic mean from X to Y") uses the
# two-way BTFFH/EEDQ matched population (n=70, same as Figure 2), not the three-way set above
# (BTFFH+EEDQ+IIDQ simultaneously, n=69) - an earlier draft silently mixed the two, giving a wrong
# non-aromatic n (30 instead of 31) and wrong percentages. Recompute the two-way population
# directly and pin it separately so the two matched sets cannot be confused again.
def matched_levels(cx, cy):
    x, y = sets[cx], sets[cy]
    ox, oy = {}, {}
    for a in sorted(set(x.sub_2_smiles) & set(y.sub_2_smiles)):
        xa, ya = x[x.sub_2_smiles == a], y[y.sub_2_smiles == a]
        com = set(xa.sub_1_smiles) & set(ya.sub_1_smiles)
        if com:
            ox[a] = xa[xa.sub_1_smiles.isin(com)].groupby("sub_1_smiles")["yield"].mean().mean()
            oy[a] = ya[ya.sub_1_smiles.isin(com)].groupby("sub_1_smiles")["yield"].mean().mean()
    return pd.Series(ox), pd.Series(oy)
_eedq2, _btffh2 = matched_levels("C41", "C3")
_df2 = pd.DataFrame({"EEDQ": _eedq2, "BTFFH": _btffh2}); _df2["state"] = _df2.index.map(STATE)
for s, want_n, want_btffh, want_eedq in [("A", 21, 0.088, 0.191), ("C", 31, 0.187, 0.094)]:
    sub = _df2[_df2.state == s]
    check(f"two-way BTFFH/EEDQ matched set, state {s}: n", int(len(sub)), want_n)
    check(f"two-way BTFFH/EEDQ matched set, state {s}: BTFFH mean", round(float(sub.BTFFH.mean()), 3), want_btffh, tol=0.001)
    check(f"two-way BTFFH/EEDQ matched set, state {s}: EEDQ mean", round(float(sub.EEDQ.mean()), 3), want_eedq, tol=0.001)


def paired(cA_, cB_):
    a, b = sets[cA_], sets[cB_]
    out = []
    for acid in sorted(set(a.sub_2_smiles) & set(b.sub_2_smiles)):
        x, y = a[a.sub_2_smiles == acid], b[b.sub_2_smiles == acid]
        com = sorted(set(x.sub_1_smiles) & set(y.sub_1_smiles))
        if not com:
            continue
        out.append({"state": STATE[acid],
                    "d": x[x.sub_1_smiles.isin(com)].groupby("sub_1_smiles")["yield"].mean().mean()
                       - y[y.sub_1_smiles.isin(com)].groupby("sub_1_smiles")["yield"].mean().mean()})
    return pd.DataFrame(out)

# an earlier draft said "most heteroaromatic acids do not [favour EEDQ/IIDQ]" - the opposite of
# what the per-acid win/loss counts show. Pin the counts directly, overall and restricted to the
# heteroaromatic class, so this direction cannot silently flip again.
for cB_, label, want_win, want_lose, want_tie, want_het_win, want_het_lose in [
        ("C41", "EEDQ", 26, 39, 5, 14, 5), ("C95", "IIDQ", 27, 40, 2, 14, 6)]:
    D = paired(cB_, "C3")  # D.d = quinoline minus BTFFH
    win = int((D.d > 0).sum()); lose = int((D.d < 0).sum()); tie = int((D.d == 0).sum())
    check(f"{label} vs BTFFH win/lose/tie", [win, lose, tie], [want_win, want_lose, want_tie])
    het = D[D.state == "A"]
    hwin = int((het.d > 0).sum()); hlose = int((het.d < 0).sum())
    check(f"{label} vs BTFFH, heteroaromatic acids only: win/lose", [hwin, hlose], [want_het_win, want_het_lose])

rng2 = np.random.default_rng(SEED)
for cB_, want, lo_, hi_ in [("C41", -0.196, -0.272, -0.121), ("C95", -0.172, -0.244, -0.103)]:
    D = paired("C3", cB_)
    guard(len(D), f"paired acids for C3 vs {cB_}")
    x = D[D.state == "A"].d.values
    y = D[D.state == "C"].d.values
    obs3 = float(x.mean() - y.mean())
    check(f"BTFFH minus {cB_}: heteroaromatic minus non-aromatic", round(obs3, 3), want, tol=0.002)
    if not QUICK:
        bs = np.array([rng2.choice(x, len(x), True).mean() - rng2.choice(y, len(y), True).mean()
                       for _ in range(20000)])
        check(f"  CI low for C3-{cB_}", round(float(np.percentile(bs, 2.5)), 2),
              round(lo_, 2), tol=0.02)
        check(f"  CI high for C3-{cB_}", round(float(np.percentile(bs, 97.5)), 2),
              round(hi_, 2), tol=0.02)
        pool = np.concatenate([x, y]); k = len(x)
        nl = np.array([(lambda q: q[:k].mean() - q[k:].mean())(rng2.permutation(pool))
                       for _ in range(20000)])
        check(f"  p for C3-{cB_}", float((np.sum(np.abs(nl) >= abs(obs3)) + 1) / 20001),
              None, hi=0.001)

# the carbocyclic-vs-non-aromatic p range quoted in the text
if not QUICK:
    ps = []
    FL2 = FL.copy()
    FL2["d41"] = FL2.C3 - FL2.C41
    FL2["d95"] = FL2.C3 - FL2.C95
    for col in ("d41", "d95"):
        x = FL2[FL2.state == "B"][col].values
        y = FL2[FL2.state == "C"][col].values
        obs4 = x.mean() - y.mean()
        pool = np.concatenate([x, y]); k = len(x)
        nl = np.array([(lambda q: q[:k].mean() - q[k:].mean())(rng2.permutation(pool))
                       for _ in range(20000)])
        ps.append((np.sum(np.abs(nl) >= abs(obs4)) + 1) / 20001)
    for col in ("d41", "d95"):
        D = paired("C3", "C41" if col == "d41" else "C95")
        x = D[D.state == "B"].d.values
        y = D[D.state == "C"].d.values
        obs4 = x.mean() - y.mean()
        pool = np.concatenate([x, y]); k = len(x)
        nl = np.array([(lambda q: q[:k].mean() - q[k:].mean())(rng2.permutation(pool))
                       for _ in range(20000)])
        ps.append((np.sum(np.abs(nl) >= abs(obs4)) + 1) / 20001)
    check("carbocyclic vs non-aromatic: smallest p over matching schemes",
          round(float(min(ps)), 3), None, lo=0.020, hi=0.040)
    check("carbocyclic vs non-aromatic: largest p over matching schemes",
          round(float(max(ps)), 3), None, lo=0.045, hi=0.085)

# ---------------------------------------------------------------- prediction
print("\n== 8. prediction from structure ==")
QUIN = [c for c in keep if act_of[c] in ("EEDQ", "IIDQ")]
EIGHT = [c for c in keep if act_of[c] in FCD]
guard(len(QUIN), "quinoline conditions for the target")
guard(len(EIGHT), "fluoroformamidinium/carbodiimide conditions for the target")
ZM = pd.DataFrame(index=acids, columns=keep, dtype=float)
for c in keep:
    ZM.loc[zmat[c].index, c] = zmat[c].values
target = (ZM[QUIN].mean(axis=1) - ZM[EIGHT].mean(axis=1)).dropna()
guard(len(target), "acids with a target value")
T = pd.DataFrame({"acid": target.index, "y": target.values})
T["hetero"] = [1.0 if STATE[a] == "A" else 0.0 for a in T.acid]
T["murcko"] = [MurckoScaffold.MurckoScaffoldSmiles(a) or "ACYCLIC" for a in T.acid]
X = T[["hetero"]].values
y = T.y.values
g = T.murcko.values


def oof_r2(Xm, yv, gv, splits):
    # solver fixed to the exact closed-form Cholesky solution (not "auto") so the near-singular
    # fits on the 18-descriptor block and the 2,048-bit fingerprint do not depend on which solver
    # a given scikit-learn build happens to pick
    p = np.empty(len(yv))
    for tr, te in splits:
        p[te] = Ridge(alpha=1.0, solver="cholesky").fit(Xm[tr], yv[tr]).predict(Xm[te])
    return 1 - ((yv - p) ** 2).sum() / ((yv - yv.mean()) ** 2).sum(), p


def oof_auc(Xm, yb, gv, splits):
    p = np.empty(len(yb))
    for tr, te in splits:
        p[te] = LogisticRegression(max_iter=5000, C=1.0, solver="lbfgs").fit(Xm[tr], yb[tr]).predict_proba(Xm[te])[:, 1]
    return roc_auc_score(yb, p), p


s5 = load_folds(T.acid, "folds_71acids.csv", groups=T.murcko)
sL = list(LeaveOneGroupOut().split(X, y, g))

# fold membership diagnostic: prints a fingerprint of exactly which acid lands in which of the
# five folds, so a run on a different machine can be compared directly rather than only compared
# on the downstream R2/AUC. A genuine partition difference would show up here; a pure numerical
# (solver/BLAS) difference would not.
_fold_of = {}
for _k, (_tr, _te) in enumerate(s5):
    for _i in _te:
        _fold_of[T.acid.iloc[_i]] = _k
_fp_input = "|".join(f"{a}:{_fold_of[a]}" for a in sorted(_fold_of))
_fp = hashlib.sha256(_fp_input.encode()).hexdigest()[:16]
print(f"  five-fold sizes: {[len(te) for _, te in s5]}  fold-membership fingerprint: {_fp}")
if os.environ.get("DUMP_FOLDS"):
    for _a in sorted(_fold_of):
        print(f"    FOLD_DUMP section8 {_a}\t{_fold_of[_a]}")
check("five-fold count", len(s5), 5)
check("every scaffold group stays inside one fold",
      all(len({_fold_of[T.acid.iloc[_i]] for _i in range(len(T)) if g[_i] == gg}) == 1
          for gg in set(g)), True)
r5, _ = oof_r2(X, y, g, s5)
rL, _ = oof_r2(X, y, g, sL)
yb = (y > 0).astype(int)
a5, prob = oof_auc(X, yb, g, s5)
aL, _ = oof_auc(X, yb, g, sL)
check("out-of-fold R2, five-fold", round(float(r5), 3), 0.194, tol=0.002)
check("out-of-fold R2, leave-one-scaffold-out", round(float(rL), 3), 0.179, tol=0.002)
check("AUC, five-fold", round(float(a5), 3), 0.635, tol=0.002)
check("AUC, leave-one-scaffold-out", round(float(aL), 3), 0.542, tol=0.003)
check("accuracy at 0.5", round(float((((prob > 0.5).astype(int)) == yb).mean()), 3), 0.704, tol=0.002)
check("majority-class baseline accuracy", round(float(max(yb.mean(), 1 - yb.mean())), 3), 0.606, tol=0.002)
tp = int((((prob > 0.5).astype(int)) == 1)[yb == 1].sum())
tn = int((((prob > 0.5).astype(int)) == 0)[yb == 0].sum())
check("acids not preferring the quinoline end, correctly called", [tn, int((yb == 0).sum())], [36, 43])
check("acids preferring the quinoline end, correctly called", [tp, int((yb == 1).sum())], [14, 28])

# richer models
from rdkit.Chem import rdFingerprintGenerator, Descriptors, rdMolDescriptors, Crippen
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
XF = np.array([np.array(gen.GetFingerprint(Chem.MolFromSmiles(a))) for a in T.acid])
rF, _ = oof_r2(XF, y, g, s5)
check("ECFP4 ridge out-of-fold R2", round(float(rF), 3), -0.137, tol=0.003)

if not QUICK:
    rng3 = np.random.default_rng(SEED)
    nl = np.array([oof_r2(X[rng3.permutation(len(T))], y, g, s5)[0] for _ in range(5000)])
    check("permutation p for the R2", float((np.sum(nl >= r5) + 1) / 5001), None, hi=0.002)

# electronic co-ordinate
print("\n== 9. the electronic co-ordinate ==")


def q_carbonyl(smi):
    m = Chem.MolFromSmiles(smi)
    hit = m.GetSubstructMatches(_cooh)[0]
    mh = Chem.AddHs(m)
    AllChem.ComputeGasteigerCharges(mh)
    return float(mh.GetAtomWithIdx(hit[0]).GetDoubleProp("_GasteigerCharge"))


Q = pd.DataFrame({"acid": acids})
Q["q"] = [q_carbonyl(a) for a in Q.acid]
Q["state"] = [STATE[a] for a in Q.acid]
Q["het"] = (Q.state == "A").astype(float)
guard(len(Q), "acids with a Gasteiger charge")
check("r(carboxyl-carbon charge, heteroaromatic state)",
      round(float(stats.pearsonr(Q.q, Q.het)[0]), 3), 0.704, tol=0.002)
gm = Q.groupby("state").q.mean()
check("charge class means", [round(float(gm[s]), 3) for s in "ABC"], [0.350, 0.336, 0.313])

# ---------------------------------------------------------------- amines
print("\n== 10. the amine axis ==")
_amide = Chem.MolFromSmarts("[NX3][CX3]=[OX1]")
_sulf = Chem.MolFromSmarts("[NX3]S(=O)=O")


def amine_class(smi):
    m = Chem.MolFromSmiles(smi)
    bad = {x[0] for patt in (_amide, _sulf) for x in m.GetSubstructMatches(patt)}
    cands = [a for a in m.GetAtoms()
             if a.GetSymbol() == "N" and a.GetIdx() not in bad and not a.GetIsAromatic()
             and a.GetTotalNumHs() > 0]
    if not cands:
        return "SULF"
    a = sorted(cands, key=lambda x: (len(x.GetNeighbors()), x.GetIdx()))[0]
    ar = [n for n in a.GetNeighbors() if n.GetIsAromatic()]
    if ar:
        for r in m.GetRingInfo().AtomRings():
            if ar[0].GetIdx() in r and all(m.GetAtomWithIdx(i).GetIsAromatic() for i in r):
                if any(m.GetAtomWithIdx(i).GetSymbol() != "C" for i in r):
                    return "AR-H"
        return "AR-C"
    return "AL-1" if a.GetTotalNumHs() >= 2 else "AL-2"


amines = sorted(raw.sub_1_smiles.unique())
ACL = {s: amine_class(s) for s in amines}
cnt = pd.Series(list(ACL.values())).value_counts()
for k, v in [("AR-C", 34), ("AL-1", 18), ("AL-2", 12), ("AR-H", 11), ("SULF", 7)]:
    check(f"amines in class {k}", int(cnt.get(k, 0)), v)

sub = ac[ac.condition_id.isin(keep)].copy()
sub["z"] = sub.groupby("condition_id")["yield"].transform(lambda v: (v - v.mean()) / v.std(ddof=1))
qq = sub[sub.condition_id.isin(QUIN)].groupby(["sub_2_smiles", "sub_1_smiles"])["z"].mean()
ff = sub[sub.condition_id.isin(EIGHT)].groupby(["sub_2_smiles", "sub_1_smiles"])["z"].mean()
cell = (qq - ff).dropna().rename("c").reset_index()
guard(len(cell), "cells with both arms")
cell["het"] = cell.sub_2_smiles.map(lambda s: 1 if STATE[s] == "A" else 0)
cell["cl"] = cell.sub_1_smiles.map(ACL)
rng4 = np.random.default_rng(SEED)
_CI_WANT = {"AR-C": (0.536, 1.185), "AR-H": (0.348, 1.370), "AL-1": (0.671, 1.813),
            "AL-2": (0.047, 1.438), "SULF": (-0.303, 0.549)}
for cl, want in [("AR-C", 0.85), ("AR-H", 0.85), ("AL-1", 1.24), ("AL-2", 0.71), ("SULF", 0.12)]:
    s = cell[cell.cl == cl]
    guard(len(s), f"cells for amine class {cl}")
    per = s.groupby(["sub_2_smiles", "het"])["c"].mean().reset_index()
    het1 = per[per.het == 1].c.values
    het0 = per[per.het == 0].c.values
    e = het1.mean() - het0.mean()
    check(f"acid-state contrast within {cl}", round(float(e), 2), want, tol=0.006)
    # percentile bootstrap over acids, same convention as the trio-vs-eight contrast above
    # (resample each het/other group independently, 20000 draws)
    _bs = np.array([rng4.choice(het1, len(het1), True).mean() - rng4.choice(het0, len(het0), True).mean()
                     for _ in range(20000)])
    _w_lo, _w_hi = _CI_WANT[cl]
    check(f"{cl}: 95% CI low", round(float(np.percentile(_bs, 2.5)), 3), _w_lo, tol=0.02)
    check(f"{cl}: 95% CI high", round(float(np.percentile(_bs, 97.5)), 3), _w_hi, tol=0.02)
signs = []
for am, s in cell.groupby("sub_1_smiles"):
    x = s[s.het == 1].c
    y = s[s.het == 0].c
    if len(x) >= 3 and len(y) >= 3:
        signs.append(x.mean() - y.mean())
guard(len(signs), "amines with at least three acids of each class")
signs = np.array(signs)
check("amines testable individually", len(signs), 30)
check("amines with the predicted sign", int((signs > 0).sum()), 27)
check("binomial p", float(stats.binomtest(int((signs > 0).sum()), len(signs), 0.5).pvalue),
      8.4e-6, tol=0.05e-6)

# SI Note 8: homogeneity of the acid-state contrast across amine classes. Test statistic is the
# SD of the per-class contrasts; the null permutes amine class labels among the amines eligible
# for each class set, holding class sizes fixed.
if not QUICK:
    _FOUR = ["AR-C", "AR-H", "AL-1", "AL-2"]
    _FIVE = _FOUR + ["SULF"]

    def _class_matrix(amine_list):
        sub = cell[cell.sub_1_smiles.isin(amine_list)]
        M = sub.pivot_table(index="sub_2_smiles", columns="sub_1_smiles", values="c")
        M = M.reindex(columns=amine_list)
        het_mask = np.array([STATE[a] == "A" for a in M.index])
        return M.values, het_mask

    def _homogeneity_p(amine_list, cls_list, seed, n=20000):
        vals_mat, het_mask = _class_matrix(amine_list)
        labels0 = np.array([ACL[a] for a in amine_list])

        def contrasts(labels):
            out = []
            for cl in cls_list:
                colmask = labels == cl
                with np.errstate(invalid="ignore"):
                    per_acid = np.nanmean(vals_mat[:, colmask], axis=1)
                out.append(np.nanmean(per_acid[het_mask]) - np.nanmean(per_acid[~het_mask]))
            return np.array(out)

        obs = contrasts(labels0)
        obs_sd = obs.std(ddof=1)
        rng8 = np.random.default_rng(seed)
        null = np.array([contrasts(rng8.permutation(labels0)).std(ddof=1) for _ in range(n)])
        return obs, float((np.sum(null >= obs_sd) + 1) / (n + 1))

    _am4 = [a for a, c in ACL.items() if c in _FOUR and a in set(cell.sub_1_smiles)]
    _am5 = [a for a, c in ACL.items() if c in _FIVE and a in set(cell.sub_1_smiles)]
    _obs4, _p4 = _homogeneity_p(_am4, _FOUR, 20260901)
    _obs5, _p5 = _homogeneity_p(_am5, _FIVE, 20260902)
    check("amine-class contrasts, four registered classes",
          [round(float(v), 2) for v in _obs4], [0.85, 0.85, 1.24, 0.71])
    check("amine-class homogeneity permutation p, four classes", round(_p4, 3), 0.609, tol=0.01)
    check("amine-class contrasts, five classes with sulfonamides",
          [round(float(v), 2) for v in _obs5], [0.85, 0.85, 1.24, 0.71, 0.12])
    check("amine-class homogeneity permutation p, five classes", round(_p5, 3), 0.145, tol=0.01)

# ---------------------------------------------------------------- quadrangles
print("\n== 11. quadrangles ==")
import itertools


def canon(cid):
    r = cA.loc[cid]
    out = set()
    for role in ("activator", "additive", "base"):
        c_ = r[f"{role}_CAS"]
        if pd.isna(c_):
            continue
        e = r[f"{role}_equiv"]
        out.add(("r", c_, round(float(e), 4) if pd.notna(e) else None))
    out.add(("s", r["solvent_CAS"], None))
    out.add(("t", None, float(r["temp_C"])))
    out.add(("h", None, float(r["time_h"])))
    return frozenset(out)


cids = sorted(ac.condition_id.unique())
CN = {c: canon(c) for c in cids}
edge = {}
for a, b in itertools.combinations(cids, 2):
    d1, d2 = CN[a] - CN[b], CN[b] - CN[a]
    if len(d1) == 1 and len(d2) == 1:
        edge[(a, b)] = (next(iter(d1)), next(iter(d2)))
        edge[(b, a)] = (next(iter(d2)), next(iter(d1)))
quads, seen = [], set()
for (a, b), L1 in edge.items():
    if a > b:
        continue
    for (c, d), L1b in edge.items():
        if c > d or (c, d) <= (a, b) or L1b != L1 or len({a, b, c, d}) != 4:
            continue
        if (a, c) not in edge or (b, d) not in edge or edge[(a, c)] != edge[(b, d)]:
            continue
        fs = frozenset((a, b, c, d))
        if fs in seen:
            continue
        seen.add(fs)
        quads.append((a, b, c, d))
guard(len(quads), "strict quadrangles")
P = ac.pivot_table(index=["sub_2_smiles", "subpair"], columns="condition_id",
                   values="yield", aggfunc="mean")
res = []
for a, b, c, d in quads:
    s = P[[a, b, c, d]].dropna()
    if len(s) == 0:
        continue
    v = s.values
    dd = (v[:, 0] - v[:, 1]) - (v[:, 2] - v[:, 3])
    res.append({"q": f"{a},{b}|{c},{d}", "n": len(s), "dd": dd.mean(),
                "acids": s.index.get_level_values(0).values, "vals": dd})
check("strict quadrangles within the screen", len(res), 77)
check("smallest shared-pair count", int(min(r["n"] for r in res)), 507, tol=0)
big = [r for r in res if abs(r["dd"]) >= 0.10]
check("quadrangles reaching |double difference| 0.10", len(big), 8)
check("largest three |double difference|",
      [float(v) for v in sorted(round(abs(r["dd"]), 3) for r in big)[-3:]], [0.113, 0.115, 0.116])
check("signs are mixed", bool(any(r["dd"] > 0 for r in big) and any(r["dd"] < 0 for r in big)), True)
# self-referential claim: "all eight lie inside the carbodiimide sub-block"
_carb = []
for r in big:
    a_, b_ = r["q"].split("|")[0].split(",")
    c_, d_ = r["q"].split("|")[1].split(",")
    acts = {act_of[x] for x in (a_, b_, c_, d_)}
    _carb.append(acts <= {"DCC", "DIC", "EDC-HCl"})
check("all eight strongest quadrangles are carbodiimide-only", [sum(_carb), len(big)], [8, 8])

if not QUICK:
    rng5 = np.random.default_rng(SEED)
    ndiff, nflip, shifts = 0, 0, []
    for r in big:
        per = pd.Series(r["vals"]).groupby(r["acids"]).mean()
        st = np.array([1 if STATE[x] == "A" else 0 for x in per.index])
        e = per[st == 1].mean() - per[st == 0].mean()
        shifts.append(abs(e))
        nl = np.array([(lambda q: per[q == 1].mean() - per[q == 0].mean())(rng5.permutation(st))
                       for _ in range(20000)])
        if (np.sum(np.abs(nl) >= abs(e)) + 1) / 20001 < 0.05:
            ndiff += 1
        if np.sign(per[st == 1].mean()) != np.sign(per[st == 0].mean()):
            nflip += 1
    check("quadrangles differing by acid state at p < 0.05", ndiff, 6)
    check("quadrangles changing sign between acid states", nflip, 0)
    check("mean |double difference|", round(float(np.mean([abs(r["dd"]) for r in big])), 3), 0.112, tol=0.002)
    check("mean shift attributable to acid state", round(float(np.mean(shifts)), 3), 0.066, tol=0.002)

# ------------------------------------------------- background over condition pairs
print()
print("== 12. background distribution over condition pairs ==")
QUINOLINE = {"EEDQ", "IIDQ"}
HALOF = {"BTFFH", "TFFH", "TCFH"}
CARBO = {"DCC", "EDC-HCl", "DIC"}
all_conds = sorted(set(ac.condition_id.unique()) & set(cA.index))
guard(len(all_conds), "conditions defined and used")
check("conditions entering the pair enumeration", len(all_conds), 94)

fam94 = {}
for lab, mem in [("quinoline", QUINOLINE), ("haloformamidinium", HALOF), ("carbodiimide", CARBO)]:
    picked = [c for c in all_conds if act_of[c] in mem]
    guard(len(picked), f"conditions in family {lab}")
    eqs = sorted({float(cA.loc[c, "activator_equiv"]) for c in picked})
    check(f"activator equivalents in {lab}", eqs, [1.5])
    for c in picked:
        fam94[c] = lab
check("quinoline conditions", sum(1 for v in fam94.values() if v == "quinoline"), 2)
check("haloformamidinium conditions", sum(1 for v in fam94.values() if v == "haloformamidinium"), 6)
check("carbodiimide conditions", sum(1 for v in fam94.values() if v == "carbodiimide"), 42)

_piv = ac.pivot_table(index="subpair", columns="condition_id", values="yield", aggfunc="mean")
_piv = _piv.reindex(columns=all_conds)
_M = _piv.to_numpy(dtype=float)
_spacid = ac.drop_duplicates("subpair").set_index("subpair").sub_2_smiles
_ai = {a: i for i, a in enumerate(acids)}
_rows = np.array([_ai[_spacid[x]] for x in _piv.index])
_G = np.zeros((len(acids), len(_piv)))
_G[_rows, np.arange(len(_piv))] = 1.0
_isA = np.array([STATE[a] == "A" for a in acids])
_isC = np.array([STATE[a] == "C" for a in acids])
import itertools as _it
_diffs, _grp, _sign = [], [], []
for i, j in _it.combinations(range(len(all_conds)), 2):
    d = _M[:, i] - _M[:, j]
    ok = ~np.isnan(d)
    if not ok.any():
        continue
    cnt = _G @ ok.astype(float)
    tot = _G @ np.where(ok, d, 0.0)
    have = cnt > 0
    if (have & _isA).sum() < 5 or (have & _isC).sum() < 5:
        continue
    ma = np.divide(tot, cnt, out=np.zeros_like(tot), where=have)
    _diffs.append(ma[have & _isA].mean() - ma[have & _isC].mean())
    fx, fy = fam94.get(all_conds[i], "other"), fam94.get(all_conds[j], "other")
    tgt = ((fx == "quinoline" and fy in ("haloformamidinium", "carbodiimide")) or
           (fy == "quinoline" and fx in ("haloformamidinium", "carbodiimide")))
    _grp.append("target" if tgt else "background")
    _sign.append(1.0 if fx == "quinoline" else -1.0)
_diffs = np.array(_diffs); _grp = np.array(_grp); _sign = np.array(_sign)
guard(len(_diffs), "estimable condition pairs")
check("no NaN among the pair differences", bool(np.isfinite(_diffs).all()), True)
check("estimable condition pairs", len(_diffs), 4371)
_t = _diffs[_grp == "target"] * _sign[_grp == "target"]
_b = _diffs[_grp == "background"]
guard(len(_t), "target condition pairs")
guard(len(_b), "background condition pairs")
check("target condition pairs", len(_t), 96)
check("background condition pairs", len(_b), 4275)
check("target mean, quinoline minus other", round(float(_t.mean()), 4), 0.1114, tol=0.0002)
check("target pairs that are positive", [int((_t > 0).sum()), len(_t)], [95, 96])

# the sign of an unordered background pair is an enumeration artefact: the negative fraction in
# one order and in the reverse order must sum to exactly one
_neg_fwd = float((_b < 0).mean())
_neg_rev = float(((-_b) < 0).mean())
check("background negative fraction, enumeration order", round(_neg_fwd, 4), 0.5460, tol=0.0002)
check("background negative fraction, reversed order", round(_neg_rev, 4), 0.4540, tol=0.0002)
check("the two orders sum to one", round(_neg_fwd + _neg_rev, 9), 1.0)
check("no background difference is exactly zero", int((_b == 0).sum()), 0)

_bs = np.concatenate([_b, -_b])
check("symmetrised background size", len(_bs), 8550)
check("symmetrised background mean is zero by construction",
      bool(abs(_bs.mean()) < 1e-12), True)
check("symmetrised background is 50% negative by construction",
      round(float((_bs < 0).mean()), 6), 0.5)
# each of the 4,275 original pairs contributes at most one of its two mirrored entries to a
# one-sided count above a positive threshold, so the count out of 8,550 symmetrised entries and
# the count of distinct ORIGINAL pairs whose |difference| reaches that threshold are the same
# number; 4,275 real pairs is the correct denominator, not the doubled 8,550
check("symmetrised background reaching the target mean",
      [int((_bs >= _t.mean()).sum()), len(_bs)], [293, 8550])
_n_reach = int((np.abs(_b) >= _t.mean()).sum())
check("distinct background pairs reaching the target mean in absolute value",
      [_n_reach, len(_b)], [293, 4275])
_reach = _n_reach / len(_b)
check("that fraction, per cent", round(100 * _reach, 2), 6.85, tol=0.01)

print()
print("== 13. the two disclosures added to Limitations ==")

# (a) the 96 target pairs share only two quinoline conditions
_qc = [c for c in all_conds if act_of[c] in QUINOLINE]
_fc = [c for c in all_conds if act_of[c] in (HALOF | CARBO)]
guard(len(_qc), "quinoline conditions")
guard(len(_fc), "comparison conditions")
check("comparison conditions on the other arm", len(_fc), 48)
_ci = {c: i for i, c in enumerate(all_conds)}
for _q, _wantmean, _wantpos in [("EEDQ", 0.1229, 48), ("IIDQ", 0.1000, 47)]:
    cid = [c for c in _qc if act_of[c] == _q][0]
    vals = []
    for f in _fc:
        d = _M[:, _ci[cid]] - _M[:, _ci[f]]
        ok = ~np.isnan(d)
        cnt = _G @ ok.astype(float)
        tot = _G @ np.where(ok, d, 0.0)
        have = cnt > 0
        ma = np.divide(tot, cnt, out=np.zeros_like(tot), where=have)
        vals.append(ma[have & _isA].mean() - ma[have & _isC].mean())
    vals = np.array(vals)
    guard(len(vals), f"{_q} target pairs")
    check(f"{_q} target pairs", len(vals), 48)
    check(f"{_q} mean over its pairs", round(float(vals.mean()), 4), _wantmean, tol=0.0002)
    check(f"{_q} positive", [int((vals > 0).sum()), len(vals)], [_wantpos, 48])

# (b) the effect sits in the low-yielding part of the condition space
_chal = set()
for a in acids:
    if STATE[a] != "A":
        continue
    m = Chem.MolFromSmiles(a)
    h = m.GetSubstructMatches(_cooh)[0]
    at = [n for n in m.GetAtomWithIdx(h[0]).GetNeighbors() if n.GetIdx() not in h[1:]][0]
    arom = [set(r) for r in m.GetRingInfo().AtomRings()
            if all(m.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    sysm = set()
    frontier = [r for r in arom if at.GetIdx() in r]
    seen = []
    while frontier:
        r = frontier.pop()
        if r in seen:
            continue
        seen.append(r); sysm |= r
        for o in arom:
            if o not in seen and o & r:
                frontier.append(o)
    els = {m.GetAtomWithIdx(i).GetSymbol() for i in sysm} - {"C"}
    if els & {"O", "S"}:
        _chal.add(a)
guard(len(_chal), "chalcogen-heteroaromatic acids")
check("acids whose ring system carries O or S", len(_chal), 6)

_per = ac.groupby(["condition_id", "sub_2_smiles"])["yield"].mean().unstack()
_Acols = [a for a in _per.columns if STATE[a] == "A"]
_Ccols = [a for a in _per.columns if STATE[a] == "C"]
_Bcols = [a for a in _per.columns if STATE[a] == "B"]
_Hcols = [a for a in _per.columns if a in _chal]
_Ncols = [a for a in _Acols if a not in _chal]
for lab, cols, want_eedq, want_c3, want_best, want_rank in [
        ("heteroaromatic", _Acols, 0.187, 0.090, 0.487, 19),
        ("O/S subset", _Hcols, 0.288, 0.076, 0.637, 19)]:
    m = _per[cols].mean(axis=1)
    guard(len(m), f"condition means for {lab}")
    m46 = m.loc[keep].sort_values(ascending=False)
    check(f"{lab}: mean under EEDQ", round(float(m.loc["C41"]), 3), want_eedq, tol=0.0006)
    check(f"{lab}: mean under BTFFH C3", round(float(m.loc["C3"]), 3), want_c3, tol=0.0006)
    check(f"{lab}: best over the 46 axis conditions", round(float(m46.iloc[0]), 3), want_best, tol=0.0006)
    check(f"{lab}: EEDQ rank of 46", list(m46.index).index("C41") + 1, want_rank)

# the O/S acids are the easiest class under the most productive quarter of conditions
_ovr = ac.groupby("condition_id")["yield"].mean()
_topq = list(_ovr.sort_values(ascending=False).index[:len(all_conds) // 4])
guard(len(_topq), "top-quartile conditions")
check("conditions in the top quartile", len(_topq), 23)
_first = 0
for c in _topq:
    vals = {"chal": _per.loc[c, _Hcols].mean(), "N": _per.loc[c, _Ncols].mean(),
            "B": _per.loc[c, _Bcols].mean(), "C": _per.loc[c, _Ccols].mean()}
    if max(vals, key=vals.get) == "chal":
        _first += 1
check("top-quartile conditions where the O/S acids are the highest class",
      [_first, len(_topq)], [len(_topq), len(_topq)])

print()
print("== 14. Table 1: the descriptor contest, one candidate set across five responses ==")
t1_COOH = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")


def t1_ringsys(mol, idx):
    R = [set(r) for r in mol.GetRingInfo().AtomRings()
         if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    s, fr, seen = set(), [r for r in R if idx in r], []
    while fr:
        r = fr.pop()
        if r in seen: continue
        seen.append(r); s |= r
        for o in R:
            if o not in seen and o & r: fr.append(o)
    return s, seen

def t1_desc(smi):
    m = Chem.MolFromSmiles(smi); h = m.GetSubstructMatches(t1_COOH)[0]
    c, oc, oh = h
    at = [n for n in m.GetAtomWithIdx(c).GetNeighbors() if n.GetIdx() not in (oc, oh)][0]
    d = {}
    arom = at.GetIsAromatic()
    if arom:
        sysm, rings = t1_ringsys(m, at.GetIdx())
        het = sorted({m.GetAtomWithIdx(i).GetSymbol() for i in sysm} - {"C"})
        own = [r for r in rings if at.GetIdx() in r][0]
        ownhet = {m.GetAtomWithIdx(i).GetSymbol() for i in own} - {"C"}
        d["D1"] = "A" if het else "B"
        d["D4"] = int(bool(ownhet))
        d["D6N"], d["D6O"], d["D6S"] = int("N" in het), int("O" in het), int("S" in het)
        d["D7"] = len([i for i in sysm if m.GetAtomWithIdx(i).GetSymbol() != "C"])
        d["D8"] = Chem.MolFragmentToSmiles(m, atomsToUse=sorted(sysm), canonical=True)
        d["D10"] = int(any(n.GetSymbol() != "C" and n.GetIsAromatic() for n in at.GetNeighbors()))
        d["D14"] = sum(1 for n in at.GetNeighbors() if n.GetIdx() in sysm
                       for q in n.GetNeighbors() if q.GetIdx() not in sysm and q.GetIdx() != c)
    else:
        d["D1"] = "C"; d["D4"] = 0; d["D6N"] = d["D6O"] = d["D6S"] = 0; d["D7"] = 0
        d["D8"] = "aliphatic"; d["D10"] = 0
        d["D14"] = sum(1 for n in at.GetNeighbors() if n.GetIdx() != c
                       for q in n.GetNeighbors() if q.GetIdx() not in (at.GetIdx(), c))
    d["D2"] = int(arom); d["D3"] = int(d["D1"] == "A")
    d["D5"] = int(any(a.GetIsAromatic() and a.GetSymbol() != "C" for a in m.GetAtoms()))
    dm = Chem.GetDistanceMatrix(m)
    tg = [a.GetIdx() for a in m.GetAtoms() if a.GetIsAromatic() and a.GetSymbol() != "C" and a.IsInRing()]
    d["D9"] = int(min(dm[c][t] for t in tg)) if tg else -1
    mh = Chem.AddHs(m); AllChem.ComputeGasteigerCharges(mh)
    d["D11"] = float(mh.GetAtomWithIdx(oh).GetDoubleProp("_GasteigerCharge"))
    d["D12"] = float(mh.GetAtomWithIdx(c).GetDoubleProp("_GasteigerCharge"))
    d["D13"] = at.GetDegree()
    d["D15"] = rdMolDescriptors.CalcNumRotatableBonds(m)
    d["D16"] = Descriptors.MolWt(m); d["D17"] = rdMolDescriptors.CalcTPSA(m)
    d["D18"] = Crippen.MolLogP(m)
    d["murcko"] = MurckoScaffold.MurckoScaffoldSmiles(mol=m) or "ACYCLIC"
    return d


D = pd.DataFrame([t1_desc(a) for a in acids], index=acids)
guard(len(D), 'descriptor table')
actf = act_of


# 'keep' is the 46-condition axis, already built above
ZM = pd.DataFrame(index=acids, columns=keep, dtype=float)
for c in keep:
    am = ac[ac.condition_id == c].groupby("sub_2_smiles")["yield"].mean()
    ZM.loc[am.index, c] = (am - am.mean()) / am.std(ddof=1)
my = ac[ac.condition_id.isin(keep)].groupby("condition_id")["yield"].mean().loc[keep].values

Zf = ZM.T.to_numpy(dtype=float)                      # conditions x acids
Dc = Zf - np.nanmean(Zf, axis=1, keepdims=True)
Dc = Dc - np.nanmean(Dc, axis=0, keepdims=True)
Df = np.where(np.isnan(Dc), 0.0, Dc)
anchor = acids.index("O=C(O)c1ccco1")
U, S, Vt = np.linalg.svd(Df, full_matrices=False)
if Vt[0][anchor] < 0: Vt, U = -Vt, -U
R = {}
R["pc1"] = pd.Series(Vt[0], index=acids)

v = (my - my.mean()).reshape(-1, 1)
Dr = Df - v @ (v.T @ Df) / (v.T @ v)
U2, S2, Vt2 = np.linalg.svd(Dr, full_matrices=False)
if Vt2[0][anchor] < 0: Vt2 = -Vt2
R["pc1_resid"] = pd.Series(Vt2[0], index=acids)
def matched(cx, cy):
    x, y = ac[ac.condition_id == cx], ac[ac.condition_id == cy]; o = {}
    for a in sorted(set(x.sub_2_smiles) & set(y.sub_2_smiles)):
        xa, ya = x[x.sub_2_smiles == a], y[y.sub_2_smiles == a]
        com = set(xa.sub_1_smiles) & set(ya.sub_1_smiles)
        if com:
            o[a] = (xa[xa.sub_1_smiles.isin(com)].groupby("sub_1_smiles")["yield"].mean().mean()
                    - ya[ya.sub_1_smiles.isin(com)].groupby("sub_1_smiles")["yield"].mean().mean())
    return pd.Series(o)
R["C3_C41"] = matched("C3", "C41")
QU = [c for c in keep if actf[c] in ("EEDQ", "IIDQ")]
FC = [c for c in keep if actf[c] in ("TFFH", "BTFFH", "TCFH", "DCC", "DIC", "EDC-HCl")]
TRIO = [c for c in keep if actf[c] in ("CDMT", "DEPBT", "PyOxim")]
R["quin_vs_eight"] = ZM[QU].mean(axis=1) - ZM[FC].mean(axis=1)
R["trio_vs_eight"] = ZM[TRIO].mean(axis=1) - ZM[FC].mean(axis=1)

SPEC = [("D1", "cat"), ("D2", "num"), ("D3", "num"), ("D4", "num"), ("D5", "num"),
        (["D6N", "D6O", "D6S"], "multi"), ("D7", "num"), ("D8", "cat"), ("D9", "cat"),
        ("D10", "num"), ("D11", "num"), ("D12", "num"), ("D13", "num"), ("D14", "num"),
        ("D15", "num"), ("D16", "num"), ("D17", "num"), ("D18", "num")]
def design(col, kind, fr):
    if kind == "num": X = fr[[col]].astype(float).values
    elif kind == "multi": X = fr[col].astype(float).values
    else: X = pd.get_dummies(fr[col].astype(str), drop_first=True).astype(float).values
    return np.hstack([np.ones((len(fr), 1)), X])
def oof(X, y, splits):
    p = np.empty(len(y)); f = []
    for tr, te in splits:
        b = np.linalg.pinv(X[tr]) @ y[tr]; p[te] = X[te] @ b
        f.append(1 - ((y[te] - p[te]) ** 2).sum() / ((y[te] - y[tr].mean()) ** 2).sum())
    return 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum(), np.array(f)
def response_folds(name, acid_order, groups):
    fname = "folds_C3_C41_70acids.csv" if name == "C3_C41" else "folds_71acids.csv"
    return load_folds(acid_order, fname, groups=groups)
print("\ncandidate encodings scored per response: %d" % len(SPEC))
print("%-16s %8s %6s   %s" % ("response", "D1 R2", "rank", "controls (MolWt/TPSA/logP)"))
for name, y in R.items():
    fr_idx = D.join(y.rename("y")).dropna(subset=["y"])
    fr = fr_idx.reset_index(drop=True)
    yy = fr.y.values.astype(float)
    splits = response_folds(name, fr_idx.index, fr.murcko)
    sc = {}
    for col, kind in SPEC:
        key = col if isinstance(col, str) else "D6"
        sc[key] = oof(design(col, kind, fr), yy, splits)[0]
    order = sorted(sc, key=lambda k: -sc[k])
    ctrl = max(sc["D16"], sc["D17"], sc["D18"])
    print("%-16s %+8.4f %4d/%d   best control %+.4f  %s" % (
        name, sc["D1"], order.index("D1") + 1, len(sc), ctrl,
        "NULL" if ctrl < sc["D1"] else "CONTROL BEATS D1"))
    # D12 (carboxyl-carbon charge) is printed on every response, not only where D1 leads: on
    # two of the five (C3_C41, quin_vs_eight) it ranks first outright, ahead of D1 - a fact the
    # manuscript's Table 1 does not surface because it only reports D1's own rank
    print("%-16s %+8.4f %4d/%d   (D12, for comparison)" % ("", sc["D12"], order.index("D12") + 1, len(sc)))
    _, f1 = oof(design("D1", "cat", fr), yy, splits)
    _, f3 = oof(design("D3", "num", fr), yy, splits)
    d = f3 - f1; se = d.std(ddof=1) / np.sqrt(len(d))
    print("%-16s   D3 %+0.4f rank %d/%d ctrl %+.4f %s  paired t(D3-D1) = %+.2f" % (
        "", sc["D3"], order.index("D3") + 1, len(sc), ctrl, "NULL" if ctrl < sc["D3"] else "BEATS D3",
        d.mean() / se if se > 0 else 0))

guard(len(R), "responses built")
check("r(leading mode condition scores, condition mean yield)",
      round(float(np.corrcoef(U[:, 0], my)[0, 1]), 3), 0.879, tol=0.002)
WANT = {"pc1": (0.057, 5), "pc1_resid": (-0.109, 15), "C3_C41": (0.140, 7),
        "quin_vs_eight": (0.156, 7), "trio_vs_eight": (0.268, 1)}
CTRL_NULL = {"pc1": True, "pc1_resid": False, "C3_C41": True,
             "quin_vs_eight": True, "trio_vs_eight": True}
for name, y in R.items():
    fr_full = D.join(y.rename("y")).dropna(subset=["y"])
    fr = fr_full.reset_index(drop=True)
    guard(len(fr), "acids with a value for response " + name)
    yy = fr.y.values.astype(float)
    splits = response_folds(name, fr_full.index, fr.murcko)
    check(name + ": no NaN in the response", bool(np.isfinite(yy).all()), True)
    sc = {}
    for col, kind in SPEC:
        key = col if isinstance(col, str) else "D6"
        sc[key] = oof(design(col, kind, fr), yy, splits)[0]
    check(name + ": candidate encodings scored", len(sc), 18)
    order = sorted(sc, key=lambda k: -sc[k])
    w_r2, w_rank = WANT[name]
    check(name + ": D1 out-of-fold R2", round(float(sc["D1"]), 3), w_r2, tol=0.0015)
    check(name + ": D1 rank of 18", order.index("D1") + 1, w_rank)
    ctrl = max(sc["D16"], sc["D17"], sc["D18"])
    check(name + ": controls stay below D1", bool(ctrl < sc["D1"]), CTRL_NULL[name])
    check(name + ": D3 vs D1", bool(sc["D3"] >= sc["D1"]),
          False if name == "trio_vs_eight" else True)
    # D3 (the two-state form actually used in the main text) gets its own R2/rank/control check,
    # not just a comparison against D1 - this is what the SI's full five-response table reports
    W_D3 = {"pc1": (0.076, 4, True), "pc1_resid": (0.015, 4, False), "C3_C41": (0.178, 5, True),
            "quin_vs_eight": (0.190, 5, True), "trio_vs_eight": (0.265, 3, True)}
    _w_r2_3, _w_rank_3, _w_null_3 = W_D3[name]
    check(name + ": D3 out-of-fold R2", round(float(sc["D3"]), 3), _w_r2_3, tol=0.0015)
    check(name + ": D3 rank of 18", order.index("D3") + 1, _w_rank_3)
    check(name + ": controls stay below D3", bool(ctrl < sc["D3"]), _w_null_3)
    if name == "trio_vs_eight":
        check(name + ": D1 leads D3 by less than 0.005",
              round(float(sc["D1"] - sc["D3"]), 4), 0.0029, tol=0.0015)
    if name in ("pc1", "quin_vs_eight"):
        _, f1 = oof(design("D1", "cat", fr), yy, splits)
        _, f3 = oof(design("D3", "num", fr), yy, splits)
        dd = f3 - f1; se = dd.std(ddof=1) / np.sqrt(len(dd))
        check(name + ": paired t(D3-D1)", round(float(dd.mean() / se), 2),
              6.14 if name == "pc1" else 2.08, tol=0.02)
    # SI Note 1's D9/D10 rank claim, checked against the same one-hot ("cat") encoding SPEC
    # actually uses for both -- not the ordinal/numeric encoding D9's own type label in the
    # candidate table might suggest
    W_D9D10 = {"pc1": (6, 2), "pc1_resid": (18, 2), "C3_C41": (12, 6),
               "quin_vs_eight": (11, 4), "trio_vs_eight": (6, 2)}
    _w_d9, _w_d10 = W_D9D10[name]
    check(name + ": D9 rank of 18", order.index("D9") + 1, _w_d9)
    check(name + ": D10 rank of 18", order.index("D10") + 1, _w_d10)
    # SI Note 2's "control that fired": TPSA (D17) on the productivity-orthogonalised response
    if name == "pc1_resid":
        check(name + ": D17 (TPSA) out-of-fold R2", round(float(sc["D17"]), 3), 0.075, tol=0.0015)
        check(name + ": D17 rank of 18", order.index("D17") + 1, 1)
        if not QUICK:
            _rng17 = np.random.default_rng(20260903)
            _X17 = design("D17", "num", fr)
            _null17 = np.array([oof(_X17[_rng17.permutation(len(yy))], yy, splits)[0]
                                 for _ in range(20000)])
            _p17 = float((np.sum(_null17 >= sc["D17"]) + 1) / 20001)
            check(name + ": D17 permutation p", round(_p17, 4), 0.0015, tol=0.0006)

print()
print("== 14b. SI Note 5: the label-free mode is preprocessing-dependent ==")
# Eight combinations of three binary preprocessing choices (within-condition z-score, complete-case
# vs mean-filled acids, whether the condition-mean-yield direction is projected out) applied to the
# 46-condition activator axis. Sign is oriented so the leading mode's acid-loading correlates
# non-negatively with the heteroaromatic indicator -- SVD has no fixed sign otherwise, and a
# different, arbitrary convention (e.g. an anchor molecule) reproduces the identical finding as its
# mirror image (rank k of 46 becomes rank 47-k, r flips sign) without changing anything real.
_rawM = pd.DataFrame(index=acids, columns=keep, dtype=float)
for _c in keep:
    _rawM[_c] = ac[ac.condition_id == _c].groupby("sub_2_smiles")["yield"].mean()
_complete_acids = list(_rawM.index[_rawM.isna().sum(axis=1) == 0])
guard(len(_complete_acids), "complete-case acids for SI Note 5")
check("complete-case acids (no missing condition)", len(_complete_acids), 65)
_EEDQ_C = [c for c in keep if actf[c] == "EEDQ"][0]
_IIDQ_C = [c for c in keep if actf[c] == "IIDQ"][0]
_HET5 = {a for a in acids if STATE[a] == "A"}


def _note5_score(M):
    Zf = M.T.to_numpy(dtype=float)
    Dc = Zf - Zf.mean(axis=1, keepdims=True)
    Dc = Dc - Dc.mean(axis=0, keepdims=True)
    U5, S5, Vt5 = np.linalg.svd(Dc, full_matrices=False)
    het_vals = np.array([1 if a in _HET5 else 0 for a in M.index])
    if np.corrcoef(Vt5[0], het_vals)[0, 1] < 0:
        Vt5, U5 = -Vt5, -U5
    var_share = (S5[0] ** 2) / (S5 ** 2).sum()
    cond_scores = pd.Series(U5[:, 0], index=keep).sort_values(ascending=False)
    eedq_rank = list(cond_scores.index).index(_EEDQ_C) + 1
    iidq_rank = list(cond_scores.index).index(_IIDQ_C) + 1
    r = np.corrcoef(Vt5[0], het_vals)[0, 1]
    return len(M.index), round(float(var_share), 3), eedq_rank, iidq_rank, round(float(r), 3)


def _note5_project(M):
    v = (my - my.mean()).reshape(-1, 1)
    Zf = M.T.to_numpy(dtype=float)
    Dc = Zf - np.nanmean(Zf, axis=1, keepdims=True)
    Dc = Dc - np.nanmean(Dc, axis=0, keepdims=True)
    Dc = Dc - v @ (v.T @ np.where(np.isnan(Dc), 0, Dc)) / (v.T @ v)
    return pd.DataFrame(Dc.T, index=M.index, columns=M.columns)


_NOTE5_WANT = {
    (True, "complete", True): (65, 0.169, 2, 1, 0.227),
    (True, "complete", False): (65, 0.244, 5, 2, 0.376),
    (True, "mean", True): (71, 0.171, 1, 5, 0.194),
    (True, "mean", False): (71, 0.223, 14, 22, 0.339),
    (False, "complete", True): (65, 0.221, 8, 9, 0.008),
    (False, "complete", False): (65, 0.496, 23, 25, 0.149),
    (False, "mean", True): (71, 0.242, 39, 37, 0.065),
    (False, "mean", False): (71, 0.523, 23, 27, 0.136),
}
for _zscore in (True, False):
    _base = _rawM.copy()
    if _zscore:
        _base = (_base - _base.mean(axis=0, skipna=True)) / _base.std(axis=0, ddof=1, skipna=True)
    for _fill in ("complete", "mean"):
        if _fill == "complete":
            _M = _base.loc[_complete_acids].copy()
        else:
            _M = _base.apply(lambda col: col.fillna(col.mean(skipna=True)), axis=0)
        for _project in (True, False):
            _Mp = _note5_project(_M) if _project else _M
            _n, _v, _e, _i, _r = _note5_score(_Mp)
            _wn, _wv, _we, _wi, _wr = _NOTE5_WANT[(_zscore, _fill, _project)]
            _label = f"z={_zscore} fill={_fill} proj={_project}"
            check(f"Note 5 [{_label}]: acids", _n, _wn)
            check(f"Note 5 [{_label}]: PC1 variance share", round(_v, 3), _wv, tol=0.002)
            check(f"Note 5 [{_label}]: EEDQ rank of 46", _e, _we)
            check(f"Note 5 [{_label}]: IIDQ rank of 46", _i, _wi)
            check(f"Note 5 [{_label}]: r(loading, heteroaryl)", round(_r, 3), _wr, tol=0.002)

print()
print("== 15. the heteroaromatic split, and the boundary among productive reagents ==")

# which heteroatom does the ring system carry?
_os_acids, _n_acids = set(), set()
for a in acids:
    if STATE[a] != "A":
        continue
    m = Chem.MolFromSmiles(a)
    h = m.GetSubstructMatches(_cooh)[0]
    at = [x for x in m.GetAtomWithIdx(h[0]).GetNeighbors() if x.GetIdx() not in h[1:]][0]
    arom = [set(r) for r in m.GetRingInfo().AtomRings()
            if all(m.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    sysm, front, seen = set(), [r for r in arom if at.GetIdx() in r], []
    while front:
        r = front.pop()
        if r in seen:
            continue
        seen.append(r); sysm |= r
        for o in arom:
            if o not in seen and o & r:
                front.append(o)
    els = {m.GetAtomWithIdx(i).GetSymbol() for i in sysm} - {"C"}
    (_os_acids if els & {"O", "S"} else _n_acids).add(a)
guard(len(_os_acids), "O/S heteroaromatic acids")
guard(len(_n_acids), "N-only heteroaromatic acids")
check("heteroaromatic acids with O or S in the ring system", len(_os_acids), 6)
check("heteroaromatic acids with N only", len(_n_acids), 15)

# Table 2, absolute yields on the three-way matched set
_sub = FL.copy()
_sub["acid"] = [r["acid"] for r in fl] if False else None
_rows = []
for a in sorted(set.intersection(*[set(g.sub_2_smiles) for g in sets.values()])):
    per = {c: g[g.sub_2_smiles == a] for c, g in sets.items()}
    com = set.intersection(*[set(x.sub_1_smiles) for x in per.values()])
    if not com:
        continue
    r = {"acid": a}
    for c, x in per.items():
        r[c] = x[x.sub_1_smiles.isin(com)].groupby("sub_1_smiles")["yield"].mean().mean()
    _rows.append(r)
_T2 = pd.DataFrame(_rows)
guard(len(_T2), "acids in the three-way matched set for Table 2")
for lab, members, want in [("O/S", _os_acids, [0.070, 0.301, 0.241]),
                           ("N only", _n_acids, [0.096, 0.145, 0.142])]:
    g = _T2[_T2.acid.isin(members)]
    guard(len(g), f"matched acids in subset {lab}")
    check(f"Table 2 {lab}: n", len(g), 6 if lab == "O/S" else 15)
    check(f"Table 2 {lab}: BTFFH / EEDQ / IIDQ",
          [round(float(g[c].mean()), 3) for c in ("C3", "C41", "C95")], want)

# position on the ordering, each subset against the non-aromatic acids
_QU = [c for c in keep if act_of[c] in ("EEDQ", "IIDQ")]
_FC = [c for c in keep if act_of[c] in FCD]
guard(len(_QU), "quinoline conditions"); guard(len(_FC), "comparison conditions")
_ZM = pd.DataFrame(index=acids, columns=keep, dtype=float)
for c in keep:
    am = ac[ac.condition_id == c].groupby("sub_2_smiles")["yield"].mean()
    _ZM.loc[am.index, c] = (am - am.mean()) / am.std(ddof=1)
_contrast = _ZM[_QU].mean(axis=1) - _ZM[_FC].mean(axis=1)
_base = _contrast[[a for a in acids if STATE[a] == "C"]].values
_rng15 = np.random.default_rng(20260910)
for lab, members, want_d, want_p in [("O/S", _os_acids, 1.95, 0.0001), ("N only", _n_acids, 1.05, 0.0022)]:
    x = _contrast[[a for a in acids if a in members]].values
    obs = x.mean() - _base.mean()
    check(f"position vs non-aromatic, {lab}", round(float(obs), 2), want_d, tol=0.006)
    if not QUICK:
        pool = np.concatenate([x, _base]); k = len(x)
        nl = np.array([(lambda q: q[:k].mean() - q[k:].mean())(_rng15.permutation(pool))
                       for _ in range(20000)])
        pv = (np.sum(np.abs(nl) >= abs(obs)) + 1) / 20001
        check(f"  its p, {lab}", round(float(pv), 4), want_p, tol=0.0025)

# the boundary: clean single-variable pairs among the most productive conditions
print()
_TOP = {"HCTU": "C23", "HATU": "C20", "PyOxim": "C37", "HDMC": "C7", "TCTU": "C8"}
def _paired(cx, cy):
    x, y = ac[ac.condition_id == cx], ac[ac.condition_id == cy]
    out = {}
    for a in sorted(set(x.sub_2_smiles) & set(y.sub_2_smiles)):
        xa, ya = x[x.sub_2_smiles == a], y[y.sub_2_smiles == a]
        com = set(xa.sub_1_smiles) & set(ya.sub_1_smiles)
        if len(com) < 3:
            continue
        out[a] = (xa[xa.sub_1_smiles.isin(com)].groupby("sub_1_smiles")["yield"].mean().mean()
                  - ya[ya.sub_1_smiles.isin(com)].groupby("sub_1_smiles")["yield"].mean().mean())
    return pd.Series(out)
_WANT = {("HCTU", "HATU"): (0.023, 0.48), ("HCTU", "PyOxim"): (-0.021, 0.41),
         ("HCTU", "HDMC"): (0.001, 0.97), ("HCTU", "TCTU"): (0.007, 0.75),
         ("HATU", "PyOxim"): (-0.043, 0.051)}
_flat, _marginal = 0, 0
for (u, v), (wd, wp) in _WANT.items():
    d = _paired(_TOP[u], _TOP[v])
    guard(len(d), f"acids paired for {u} vs {v}")
    hx = d[[a for a in d.index if STATE[a] == "A"]].values
    cx = d[[a for a in d.index if STATE[a] == "C"]].values
    obs = hx.mean() - cx.mean()
    check(f"{u} - {v}: class difference", round(float(obs), 3), wd, tol=0.0015)
    if not QUICK:
        pool = np.concatenate([hx, cx]); k = len(hx)
        nl = np.array([(lambda q: q[:k].mean() - q[k:].mean())(_rng15.permutation(pool))
                       for _ in range(20000)])
        pv = (np.sum(np.abs(nl) >= abs(obs)) + 1) / 20001
        check(f"  its p", round(float(pv), 2), wp, tol=0.03)
        if (u, v) == ("HATU", "PyOxim"):
            # this one sits on the threshold; assert the band, not the side it falls on
            check("  HATU-PyOxim p sits at the conventional threshold",
                  float(pv), None, lo=0.03, hi=0.08)
        elif pv >= 0.35:
            _flat += 1
if not QUICK:
    check("productive-end comparisons with no class difference (p >= 0.35)", _flat, 4)

# ---------------------------------------------------------------- additive effects
print()
print("== 16. additive effects: matched no-additive vs with-additive condition pairs ==")


def _recipe_key(row):
    return (row.activator_CAS, round(float(row.activator_equiv), 3),
            row.base_CAS if pd.notna(row.base_CAS) else None,
            round(float(row.base_equiv), 3) if pd.notna(row.base_equiv) else None)


cA["_recipe"] = [_recipe_key(r) for _, r in cA.iterrows()]
_add_pairs = []
for _key, _grp in cA.groupby("_recipe"):
    _no = _grp[_grp.additive_CAS.isna()]
    _wi = _grp[_grp.additive_CAS.notna()]
    if len(_no) == 0 or len(_wi) == 0:
        continue
    for _nid in _no.index:
        for _wid in _wi.index:
            _add_pairs.append((_nid, _wid))
guard(len(_add_pairs), "additive-present-vs-absent condition pairs")
check("additive-present-vs-absent condition pairs sharing an activator+base", len(_add_pairs), 21)

_rng16 = np.random.default_rng(20260912)
_add_rows = []
for _nid, _wid in _add_pairs:
    _d = matched(_nid, _wid)  # matched(cx,cy) = cx minus cy; here no-additive minus with-additive
    _d = -_d  # report as with-additive minus no-additive, i.e. the additive's own effect
    _st = np.array([STATE[a] for a in _d.index])
    _isA = _st == "A"
    _row = {"no_add": _nid, "with_add": _wid, "activator": cA.loc[_nid, "activator_name"],
            "additive": cA.loc[_wid, "additive_name"], "n_acids": len(_d),
            "mean_diff_pp": round(100 * float(_d.mean()), 2)}
    if _isA.sum() >= 3 and (~_isA).sum() >= 3:
        _obs = float(_d.values[_isA].mean() - _d.values[~_isA].mean())
        _pool = _d.values.copy(); _k = int(_isA.sum())
        _nl = np.array([(lambda q: q[:_k].mean() - q[_k:].mean())(_rng16.permutation(_pool))
                        for _ in range(20000)])
        _row["interaction_pp"] = round(100 * _obs, 2)
        _row["p"] = (np.sum(np.abs(_nl) >= abs(_obs)) + 1) / 20001
    _add_rows.append(_row)
ADD = pd.DataFrame(_add_rows)
guard(len(ADD), "additive comparisons scored")

for _nid, _wid, _want in [("C77", "C78", 32.74), ("C87", "C88", 23.12),
                          ("C47", "C60", 30.63), ("C20", "C19", -8.93)]:
    _v = float(ADD.loc[(ADD.no_add == _nid) & (ADD.with_add == _wid), "mean_diff_pp"].iloc[0])
    check(f"additive effect, {_nid} -> {_wid}", round(_v, 2), _want, tol=0.02)

# does the additive effect's sign/size track activator chemistry? carbodiimides (DCC, EDC-HCl)
# are textbook-known to need a nucleophile trap (HOAt/HOBt/Oxyma-type/active-ester-forming
# additives) to suppress O-acylisourea rearrangement and racemisation [9-11,22-25]; preformed
# aminium/uronium/phosphonium activators already generate a reactive species without one
_carbodiimide = {"DCC", "EDC-HCl"}
_non_ae_additive = {"DMAP"}  # not active-ester-forming; a nucleophilic acyl-transfer catalyst instead
ADD["_is_cbdi"] = ADD.activator.isin(_carbodiimide)
ADD["_ae_additive"] = ~ADD.additive.isin(_non_ae_additive)
_grp_cbdi_ae = ADD.loc[ADD._is_cbdi & ADD._ae_additive, "mean_diff_pp"]
_grp_cbdi_other = ADD.loc[ADD._is_cbdi & ~ADD._ae_additive, "mean_diff_pp"]
_grp_other = ADD.loc[~ADD._is_cbdi, "mean_diff_pp"]
guard(len(_grp_cbdi_ae), "carbodiimide + active-ester-additive comparisons")
guard(len(_grp_other), "non-carbodiimide-activator comparisons")
check("carbodiimide + active-ester additive: count", len(_grp_cbdi_ae), 12)
check("carbodiimide + active-ester additive: all positive", bool((_grp_cbdi_ae > 0).all()), True)
check("carbodiimide + active-ester additive: smallest effect (pp)", round(float(_grp_cbdi_ae.min()), 1), 3.7, tol=0.1)
check("non-carbodiimide activator: count", len(_grp_other), 8)
check("non-carbodiimide activator: largest effect (pp)", round(float(_grp_other.max()), 1), 0.3, tol=0.1)
check("the two groups do not overlap", bool(_grp_cbdi_ae.min() > _grp_other.max()), True)
check("carbodiimide + DMAP (not active-ester-forming) breaks the pattern",
      round(float(_grp_cbdi_other.iloc[0]), 1), 0.2, tol=0.1)
_u16, _pmw16 = stats.mannwhitneyu(_grp_cbdi_ae, _grp_other, alternative="greater")
check("Mann-Whitney p (carbodiimide+active-ester-additive > non-carbodiimide), one-sided",
      bool(_pmw16 < 0.0001), True)
print(f"  carbodiimide+active-ester additive: n={len(_grp_cbdi_ae)}, {_grp_cbdi_ae.min():.1f} to {_grp_cbdi_ae.max():.1f} pp, mean {_grp_cbdi_ae.mean():.1f}")
print(f"  non-carbodiimide activator: n={len(_grp_other)}, {_grp_other.min():.1f} to {_grp_other.max():.1f} pp, mean {_grp_other.mean():.1f}")
print(f"  Mann-Whitney p = {_pmw16:.2e}")

if not QUICK:
    _pvals = ADD["p"].dropna().sort_values().values
    _m16 = len(_pvals)
    guard(_m16, "additive comparisons with an interaction test")
    check("additive comparisons with an interaction test", _m16, 21)
    _q16 = np.empty(_m16)
    _q16[-1] = _pvals[-1]
    for _i in range(_m16 - 2, -1, -1):
        _q16[_i] = min(_q16[_i + 1], _pvals[_i] * _m16 / (_i + 1))
    check("additive x acid-class interactions surviving BH correction at q<0.05",
          int((_q16 < 0.05).sum()), 0)
    # this exact value (permutation Monte Carlo, so given a wider tolerance than an exact-formula
    # statistic) is quoted in the main text, so it is pinned here rather than only printed
    _smallest_row = ADD.loc[ADD.p == ADD.p.min()].iloc[0]
    check("smallest interaction p", round(float(_pvals[0]), 4), 0.006, tol=0.004)
    check("smallest interaction p belongs to", [_smallest_row.activator, _smallest_row.additive],
          ["HATU", "HOAT"])
    print(f"  smallest interaction p = {_pvals[0]:.4f} (q = {_q16[0]:.3f})")

# ---------------------------------------------------------------- fair descriptor comparison
print()
print("== 17. two-state descriptor against carboxyl-charge, same split, same target ==")
_targets17 = {"EEDQ - BTFFH": (matched("C3", "C41") * -1, 0.164, 0.279, 0.255, 0.15),
              "IIDQ - BTFFH": (matched("C3", "C95") * -1, 0.139, 0.264, 0.236, 0.11)}
for _label, (_tgt, _want_s, _want_c, _want_j, _want_p) in _targets17.items():
    guard(len(_tgt), f"acids for the {_label} target")
    _acidsT = list(_tgt.index)
    _yv = _tgt.values
    _st17 = np.array([1.0 if STATE[a] == "A" else 0.0 for a in _acidsT])
    _chg17 = np.array([q_carbonyl(a) for a in _acidsT])
    _mk17 = np.array([MurckoScaffold.MurckoScaffoldSmiles(a) or "ACYCLIC" for a in _acidsT])
    _splits17 = list(LeaveOneGroupOut().split(_st17.reshape(-1, 1), _yv, _mk17))

    def _ols_oof(*xcols):
        Xd = np.hstack([np.ones((len(_yv), 1))] + [c.reshape(-1, 1) for c in xcols])
        p = np.empty(len(_yv))
        for tr, te in _splits17:
            b = np.linalg.pinv(Xd[tr]) @ _yv[tr]
            p[te] = Xd[te] @ b
        return 1 - ((_yv - p) ** 2).sum() / ((_yv - _yv.mean()) ** 2).sum()

    def _baseline_oof():
        p = np.empty(len(_yv))
        for tr, te in _splits17:
            p[te] = _yv[tr].mean()
        return 1 - ((_yv - p) ** 2).sum() / ((_yv - _yv.mean()) ** 2).sum()

    _r2b = _baseline_oof()
    _r2s = _ols_oof(_st17)
    _r2c = _ols_oof(_chg17)
    _r2j = _ols_oof(_st17, _chg17)
    print(f"  {_label} (n={len(_tgt)}): baseline R2={_r2b:+.3f}  two-state R2={_r2s:+.3f}  "
          f"charge R2={_r2c:+.3f}  joint R2={_r2j:+.3f}")
    # a train-fold-mean baseline is expected to land slightly negative under leave-one-group-out
    # with this few acids - that is ordinary cross-validation noise, not a bug
    check(f"{_label}: baseline (train-fold mean) R2, near zero", round(float(_r2b), 3),
          None, lo=-0.08, hi=0.02)
    check(f"{_label}: two-state R2, leave-one-scaffold-out", round(float(_r2s), 3), _want_s, tol=0.01)
    check(f"{_label}: carboxyl-charge R2, leave-one-scaffold-out", round(float(_r2c), 3), _want_c, tol=0.01)
    check(f"{_label}: charge outperforms two-state on this target", bool(_r2c > _r2s), True)
    # does the structural descriptor add anything once charge is already in the model, under the
    # same leave-one-scaffold-out protocol used throughout - or does adding a correlated, weaker
    # predictor just add variance? pinned because the manuscript now reads the answer off this
    check(f"{_label}: joint (structural + charge) R2", round(float(_r2j), 3), _want_j, tol=0.01)
    check(f"{_label}: joint model does not beat charge alone", bool(_r2j <= _r2c), True)

    def _residualize(_x, _on):
        _Xd = np.hstack([np.ones((len(_on), 1)), _on.reshape(-1, 1)])
        _b = np.linalg.pinv(_Xd) @ _x
        return _x - _Xd @ _b
    _partial_r = float(np.corrcoef(_residualize(_st17, _chg17), _residualize(_yv, _chg17))[0, 1])
    check(f"{_label}: partial r(structural, target | charge)", round(_partial_r, 2), _want_p, tol=0.02)

# ---------------------------------------------------------------- summary
print("\n" + "=" * 78)
npass = sum(1 for ok, *_ in RESULTS if ok)
print(f"{npass} of {len(RESULTS)} checks reproduced from DATASET_*.csv")
if npass != len(RESULTS):
    print("\nFAILURES:")
    for ok, name, got, exp in RESULTS:
        if not ok:
            print(f"  {name}: got {got}, expected {exp}")
    sys.exit(1)
print("all checks passed")
