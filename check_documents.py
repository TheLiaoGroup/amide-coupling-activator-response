"""Check that the tables and prose numbers in MANUSCRIPT_acylation.md and SI_acylation.md agree
with the Source Data CSVs the figures were drawn from.

REPRODUCE_acylation.py checks the manuscript against the dataset directly. This script checks a
different seam: that the numbers transcribed into prose and markdown tables match the CSVs a
reader would download. A figure agreeing with the data while the text disagrees with the figure
has happened on this project before.

Every comparison asserts it matched at least one row first.
"""
import os, re, sys
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "source_data")
FAIL = []


def note(ok, msg):
    print(("  PASS  " if ok else "  FAIL  ") + msg)
    if not ok:
        FAIL.append(msg)


def md_tables(path):
    """every markdown table in a file, as a list of row-lists"""
    rows, cur = [], []
    for line in open(path, encoding="utf-8"):
        if line.lstrip().startswith("|"):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if set("".join(cells)) <= set("-: "):
                continue
            cur.append(cells)
        elif cur:
            rows.append(cur); cur = []
    if cur:
        rows.append(cur)
    return rows


def num(s):
    s = s.replace("−", "-").replace("*", "").strip()
    m = re.fullmatch(r"[+-]?\d*\.?\d+", s)
    return float(m.group()) if m else None


MS = open(os.path.join(HERE, "MANUSCRIPT_acylation.md"), encoding="utf-8").read()

print("== four additive effects named in the main text against Figure3_source_data.csv ==")
F3 = pd.read_csv(os.path.join(SRC, "Figure3_source_data.csv"))
assert len(F3) == 21, f"Figure 3 source data is not 21 rows (has {len(F3)})"
for na, wa, want_text in [("C77", "C78", "32.7"), ("C87", "C88", "23.1"),
                          ("C47", "C60", "30.6"), ("C20", "C19", "-8.9")]:
    row = F3[(F3.no_additive == na) & (F3.with_additive == wa)]
    note(len(row) == 1, f"{na}->{wa}: found in Source Data ({len(row)} rows)")
    v = float(row.mean_diff_pp.iloc[0])
    ok = abs(v - float(want_text)) < 0.05
    note(ok, f"{na}->{wa}: prose says {want_text} pp, Source Data says {v:.1f} pp")
    note(want_text.lstrip("-") in MS, f"{na}->{wa}: figure '{want_text}' appears in the manuscript text")
note(int((F3["interaction_pp"].notna()).sum()) == 21,
     f"21 additive comparisons carry an interaction test (has {int(F3['interaction_pp'].notna().sum())})")

# Figure 3b's "0 of 21, q<0.05" title text used to be a plain string no checker looked behind -
# now that build_figures.py exports p and q_bh, verify the exported values actually back it up.
note({"p", "q_bh"}.issubset(F3.columns), "Figure3_source_data.csv exports p and q_bh columns")
if {"p", "q_bh"}.issubset(F3.columns):
    n_sig = int((F3.q_bh < 0.05).sum())
    note(n_sig == 0, f"additive x class interactions with q<0.05 in Source Data: {n_sig} (figure title says 0 of 21)")
    _min_p = round(float(F3.p.min()), 4)
    note(_min_p == 0.0061, f"smallest interaction p in Source Data: {_min_p} (prose says 0.0061)")
    note("0.0061" in MS, "the smallest interaction p (0.0061) appears in the manuscript text")
    note(re.search(r"no comparison reaches\s*q\s*<\s*0\.05", MS) is not None,
         "the manuscript's 'no comparison reaches q<0.05' claim is now backed by an exported q_bh column")

print("\n== fair descriptor comparison against Figure4_source_data.csv ==")
F4 = pd.read_csv(os.path.join(SRC, "Figure4_source_data.csv"))
assert len(F4) == 8, f"Figure 4 source data is not 8 rows (has {len(F4)})"
for target, want_state, want_charge, want_joint in [("EEDQ - BTFFH", 0.164, 0.279, 0.255),
                                                     ("IIDQ - BTFFH", 0.139, 0.264, 0.236)]:
    rs = F4[(F4.target == target) & (F4.model == "two-state descriptor")].R2.iloc[0]
    rc = F4[(F4.target == target) & (F4.model == "carboxyl-carbon charge")].R2.iloc[0]
    rj = F4[(F4.target == target) & (F4.model == "structural + charge (joint)")].R2.iloc[0]
    note(abs(float(rs) - want_state) < 0.002, f"{target}: two-state R2 {rs:.3f} vs prose {want_state}")
    note(abs(float(rc) - want_charge) < 0.002, f"{target}: charge R2 {rc:.3f} vs prose {want_charge}")
    note(abs(float(rj) - want_joint) < 0.002, f"{target}: joint R2 {rj:.3f} vs prose {want_joint}")
    note(float(rc) > float(rs), f"{target}: charge outperforms two-state in Source Data")
    note(float(rj) <= float(rc), f"{target}: joint model does not beat charge alone in Source Data")
for token in ("0.164", "0.279", "0.139", "0.264", "0.255", "0.236"):
    note(token in MS, f"figure '{token}' appears in the manuscript text")

# a sign-flip bug (Figure 4's scatter built on -d_eedq while labelled "EEDQ - BTFFH") passed every
# R2 check above because R2 is invariant to negating the target - so also assert the two figures'
# per-acid values for the same quantity actually agree in sign, not just in R2
FSCAT = pd.read_csv(os.path.join(SRC, "Figure4_scatter_data.csv"))
_f2 = pd.read_csv(os.path.join(SRC, "Figure2_source_data.csv")).set_index("acid_smiles")
for target, col in [("EEDQ - BTFFH", "EEDQ_minus_BTFFH"), ("IIDQ - BTFFH", "IIDQ_minus_BTFFH")]:
    s0 = FSCAT[FSCAT.target == target].set_index("acid_smiles")
    ref = _f2.loc[s0.index, col]
    note(bool((s0["observed"] - ref).abs().max() < 1e-9),
         f"{target}: Figure4_scatter_data.csv 'observed' matches Figure2_source_data.csv sign, acid by acid")

# a check that only confirms the RIGHT numbers are present would miss a wrong number sitting
# alongside them (this happened once: the abstract said "R2 up to 0.174", a stale value from an
# earlier draft, while the correct 0.164 was also present elsewhere in the text) - so also assert
# that every "R2 ... 0.xxx"-shaped number in the manuscript is one of the values Source Data
# actually produced, not just that the right ones showed up somewhere
_valid_r2 = {f"{v:.3f}".lstrip("0") for v in F4.R2} | {f"{v:.3f}" for v in F4.R2}
for _m in re.finditer(r"R²\s*(?:=|≈|up to|falls to)?\s*(-?0\.\d{2,4})", MS):
    _tok = _m.group(1)
    _tok3 = f"{float(_tok):.3f}"
    note(_tok3 in {f"{v:.3f}" for v in F4.R2} or abs(float(_tok) + 0.04) < 0.02,
         f"R² value near manuscript text position {_m.start()}: {_tok!r} matches a Figure 4 Source Data value (or the ~-0.04 baseline)")

print("\n== Figure 2 caption is not the mirror-image of the Results counts ==")
# an earlier draft fixed the Results prose (14 of 21 heteroaromatic acids favour the quinoline)
# but left the Figure 2 caption saying the opposite ("most heteroaromatic acids do not favour the
# quinoline") - check_documents.py never caught it because it only checked that the correct counts
# appeared *somewhere*, not that the old reversed claim was gone.
note("do not favour the quinoline" not in MS,
     "manuscript does not contain the reversed 'most heteroaromatic acids do not favour the quinoline' claim")
note(re.search(r"14\s+of\s+the\s+21\s+heteroaromatic\s+acids\s+favour\s+the\s+quinoline", MS) is not None,
     "Figure 2 caption states the correct heteroaromatic win count (14 of 21), not just the Results prose")

print("\n== activator-axis coordinates: EEDQ and IIDQ are not swapped ==")
FS1_ax = pd.read_csv(os.path.join(SRC, "SFigure1_source_data.csv")).set_index("activator")
_eedq_z = round(float(FS1_ax.loc["EEDQ", "z_gap"]), 2)
_iidq_z = round(float(FS1_ax.loc["IIDQ", "z_gap"]), 2)
note(_eedq_z > _iidq_z, f"EEDQ z_gap ({_eedq_z}) > IIDQ z_gap ({_iidq_z}) in Source Data")
note(f"+{_eedq_z:.2f} and +{_iidq_z:.2f}" in MS or f"{_eedq_z:.2f} and +{_iidq_z:.2f}" in MS,
     f"manuscript states EEDQ ({_eedq_z}) before IIDQ ({_iidq_z}) in that order, matching Source Data")

print("\n== trio-vs-eight (+1.05) is not conflated with the 96-pair raw background (11.1pp) ==")
# these are two different reagent sets (3 conditions vs 8, standardised units; 2 conditions vs 48,
# raw percentage points) computed by different scripts - an earlier draft called them equivalent.
note("plus the\ndihydroquinolines give a heteroaromatic" not in MS and "plus the dihydroquinolines give a heteroaromatic" not in MS,
     "the +1.05 trio-vs-eight contrast is not attributed to 'the trio plus the dihydroquinolines'")
note(re.search(r"\+1\.05[^.]*equivalent", MS) is None,
     "the +1.05 standardised contrast is not called 'equivalent to' the raw percentage-point background figure")
note("48 haloformamidinium" in MS or "48 haloformamidinium/carbodiimide" in MS,
     "the 96-pair background comparison states its own reagent-set size (48), not the axis figure's eight")

print("\n== raw paired differences against Figure2_source_data.csv ==")
F2 = pd.read_csv(os.path.join(SRC, "Figure2_source_data.csv"))
n_eedq = int(F2.EEDQ_minus_BTFFH.notna().sum())
n_iidq = int(F2.IIDQ_minus_BTFFH.notna().sum())
note(n_eedq == 70, f"EEDQ-BTFFH acid count in Source Data: {n_eedq} (text says 70)")
note(n_iidq == 69, f"IIDQ-BTFFH acid count in Source Data: {n_iidq} (text says 69)")
frac_eedq = round(100 * float((F2.EEDQ_minus_BTFFH.dropna() > 0).mean()))
frac_iidq = round(100 * float((F2.IIDQ_minus_BTFFH.dropna() > 0).mean()))
note(frac_eedq == 37, f"fraction of acids favouring EEDQ: {frac_eedq}% (text says 37%)")
note(frac_iidq == 39, f"fraction of acids favouring IIDQ: {frac_iidq}% (text says 39%)")
for token in ("37% of 70 matched acids", "39%"):
    note(token in MS or token.split("%")[0] + "%" in MS, f"figure appears in text: {token!r}")

# an earlier draft said "most heteroaromatic acids do not [favour EEDQ]" - the opposite of what the
# per-acid data show. Recompute the win/loss/tie counts directly, overall and restricted to the
# heteroaromatic class, and assert the manuscript states the heteroaromatic-class counts, not just
# the overall percentage (which is compatible with either direction within the class).
for target, col, want_win, want_lose, want_tie, want_het_win, want_het_lose in [
        ("EEDQ", "EEDQ_minus_BTFFH", 26, 39, 5, 14, 5),
        ("IIDQ", "IIDQ_minus_BTFFH", 27, 40, 2, 14, 6)]:
    s = F2[["acid_state", col]].dropna(subset=[col])
    win = int((s[col] > 0).sum()); lose = int((s[col] < 0).sum()); tie = int((s[col] == 0).sum())
    note((win, lose, tie) == (want_win, want_lose, want_tie),
         f"{target} vs BTFFH win/lose/tie: got {(win, lose, tie)}, want {(want_win, want_lose, want_tie)}")
    het = s[s.acid_state == "A"]
    hwin = int((het[col] > 0).sum()); hlose = int((het[col] < 0).sum())
    note((hwin, hlose) == (want_het_win, want_het_lose),
         f"{target} vs BTFFH among heteroaromatic acids: got win={hwin} lose={hlose}, want win={want_het_win} lose={want_het_lose}")
    note(hwin > hlose, f"{target}: majority of heteroaromatic acids favour the quinoline (win={hwin} > lose={hlose})")
note(str(want_win) + " of 70" in MS or "26 of 70" in MS, "EEDQ win count (26 of 70) appears in the manuscript text")
note("14 of the 21 heteroaromatic" in MS or "14 of 21" in MS,
     "the heteroaromatic-class win count (14 of 21) appears in the manuscript text, not just the overall percentage")
note(re.search(r"8\.8%\s+to\s+19\.1%", MS) is not None,
     "corrected heteroaromatic class-mean (8.8% to 19.1%, n=21) appears in the manuscript text")
note(re.search(r"31-acid\s+non-aromatic", MS) is not None,
     "class-mean sample size is stated as 31 acids, not the earlier wrong 30")
note(re.search(r"18\.7%\s+to\s+9\.4%", MS) is not None,
     "corrected non-aromatic class-mean (18.7% to 9.4%) appears in the manuscript text")

print("\n== coverage figure against Figure1_source_data.csv ==")
F1 = pd.read_csv(os.path.join(SRC, "Figure1_source_data.csv"))
n_pairs = len(F1.drop_duplicates(["acid_smiles", "amine_smiles"]))
note(n_pairs == 578, f"measured pairs in Source Data: {n_pairs} (text says 578)")
n_acids = F1.acid_smiles.nunique(); n_amines = F1.amine_smiles.nunique()
note(n_acids == 71 and n_amines == 82, f"acids/amines: {n_acids}/{n_amines} (text says 71/82)")
coverage = round(100 * n_pairs / (n_acids * n_amines), 1)
note(coverage == 9.9, f"coverage: {coverage}% (text says 9.9%)")

print("\n== SI Note 6: the 46-condition table against SFigure1_source_data.csv ==")
FS1 = pd.read_csv(os.path.join(SRC, "SFigure1_source_data.csv")).set_index("condition_id")
assert len(FS1) == 46, "SFigure 1 source data is not 46 rows"
SI = open(os.path.join(HERE, "SI_acylation.md"), encoding="utf-8").read()
tabs = md_tables(os.path.join(HERE, "SI_acylation.md"))
axis_tab = [t for t in tabs if t and t[0][0] == "condition" and "residual" in t[0][-1]]
note(len(axis_tab) == 1, f"found exactly one 46-condition table (found {len(axis_tab)})")
if axis_tab:
    body = axis_tab[0][1:]
    note(len(body) == 46, f"table has 46 data rows (has {len(body)})")
    n = 0
    for r in body:
        cid = r[0]
        if cid not in FS1.index:
            note(False, f"{cid} is not in the Source Data"); continue
        s = FS1.loc[cid]
        for col, idx in [("activator", 1), ("base", 2), ("reagent_family", 3)]:
            if str(s[col]) != r[idx]:
                note(False, f"{cid} {col}: table {r[idx]!r} vs source {s[col]!r}")
        for col, idx, tol in [("n_wells", 4, 0), ("mean_yield", 5, 0.0006),
                              ("z_gap", 6, 0.0006), ("z_gap_residual", 7, 0.0006)]:
            v = num(r[idx])
            if v is None or abs(v - float(s[col])) > tol:
                note(False, f"{cid} {col}: table {r[idx]} vs source {s[col]}")
        n += 1
    assert n > 0, "compared zero rows"
    note(True, f"compared {n} conditions, every field")

print("\n== SI Note 9: the quadrangle table against SFigure4_source_data.csv ==")
FS3 = pd.read_csv(os.path.join(SRC, "SFigure4_source_data.csv")).set_index("quadrangle")
assert len(FS3) == 8, "SFigure 4 source data is not 8 rows"
qt = [t for t in tabs if t and t[0][0] == "quadrangle" and "double difference" in " ".join(t[0])]
note(len(qt) == 1, f"found exactly one quadrangle table (found {len(qt)})")
if qt:
    body = qt[0][1:]
    note(len(body) == 8, f"quadrangle table has 8 data rows (has {len(body)})")
    n = 0
    for r in body:
        key = r[0].replace(" / ", "|").replace(" ", "")
        if key not in FS3.index:
            note(False, f"{r[0]!r} -> {key!r} is not in the Source Data"); continue
        s = FS3.loc[key]
        for idx, col, tol in [(2, "n_substrate_pairs", 0), (3, "double_difference", 0.0006),
                              (5, "dd_heteroaromatic", 0.0006), (6, "dd_other", 0.0006)]:
            v = num(r[idx])
            if v is None or abs(v - float(s[col])) > tol:
                note(False, f"{key} {col}: table {r[idx]} vs source {s[col]}")
        lo, hi = [num(x) for x in re.split(r"\s+to\s+", r[4])]
        if lo is None or abs(lo - float(s.ci_lo)) > 0.0015 or abs(hi - float(s.ci_hi)) > 0.0015:
            note(False, f"{key} CI: table {r[4]} vs source [{s.ci_lo:.4f}, {s.ci_hi:.4f}]")
        n += 1
    assert n > 0, "compared zero quadrangles"
    note(True, f"compared {n} quadrangles, every field")

print("\n== SI Note 2: the full five-response descriptor table (D3) ==")
d3_tab = [t for t in tabs if t and t[0][0] == "response" and "D3" in " ".join(t[0])]
note(len(d3_tab) == 1, f"found the five-response D3 table (found {len(d3_tab)})")
if d3_tab:
    note(len(d3_tab[0]) - 1 == 5, f"table has 5 responses (has {len(d3_tab[0]) - 1})")
    _n_prereg = sum(1 for r in d3_tab[0][1:] if r[1].startswith("yes"))
    _n_posthoc = sum(1 for r in d3_tab[0][1:] if "post hoc" in r[1])
    note((_n_prereg, _n_posthoc) == (3, 2),
         f"of the five responses, pre-registered/post-hoc split is {(_n_prereg, _n_posthoc)} (want (3, 2))")
    note(re.search(r"three\s+(?:of\s+that\s+table's\s+five\s+responses\s+were\s+pre-registered|were\s+pre-registered)", MS) is not None,
         "manuscript states three of the five responses were pre-registered, not that all five were")

print("\n== pre-registration status of the quinoline-family-contrast target ==")
# an earlier draft called this target "pre-registered" in Methods and SI Note 3, while SI Note 2's
# own five-response table marks it "no, post hoc" - both must agree, and neither the Methods nor SI
# text may re-introduce the false "pre-registered quinoline-family-contrast" phrase.
qrow = [r for r in d3_tab[0][1:] if r[0] == "quinoline vs the eight"] if d3_tab else []
note(len(qrow) == 1, "found the 'quinoline vs the eight' row in the SI Note 2 five-response table")
if qrow:
    note("post hoc" in qrow[0][1], f"SI Note 2 marks 'quinoline vs the eight' as post hoc (has {qrow[0][1]!r})")
_bad_phrase = re.compile(r"pre-registered\s+quinoline-family-contrast")
note(not _bad_phrase.search(MS),
     "manuscript Methods does not call the quinoline-family-contrast target pre-registered")
note(not _bad_phrase.search(SI),
     "SI does not call the quinoline-family-contrast target pre-registered")
note(re.search(r"post[- ]hoc\s+target", MS) is not None,
     "manuscript Methods states the quinoline-family-contrast target is post-hoc")

print("\n== SI Note 10 family enumeration against SFigure2_source_data.csv ==")
FS2 = pd.read_csv(os.path.join(SRC, "SFigure2_source_data.csv"))
assert len(FS2) > 0
note(len(FS2) == 4371, f"SFigure 2 Source Data has one row per condition pair (has {len(FS2)})")
note(int(FS2.is_target.sum()) == 96, f"96 target pairs in Source Data (has {int(FS2.is_target.sum())})")
fam_conds = {}
for side in ("x", "y"):
    for f, c in zip(FS2["family_" + side], FS2["cond_" + side]):
        if f != "other":
            fam_conds.setdefault(f, set()).add(c)
si_tab = [t for t in tabs if t and t[0][0] == "family" and "condition ids" in " ".join(t[0])]
note(len(si_tab) == 1, f"found the SI family table (found {len(si_tab)})")
if si_tab:
    si_conds, n = {}, 0
    for r in si_tab[0][1:]:
        si_conds.setdefault(r[0], set()).update(x.strip() for x in r[5].split(","))
        note(r[4] == "1.5", f"{r[1]}: activator equivalents listed as {r[4]}, expected 1.5")
        n += 1
    assert n > 0, "compared zero family rows"
    for f in sorted(set(fam_conds) | set(si_conds)):
        note(fam_conds.get(f) == si_conds.get(f),
             f"family {f}: SI lists {len(si_conds.get(f, []))} conditions, "
             f"Source Data has {len(fam_conds.get(f, []))}")
    note(True, f"compared {n} activators across three families")

print("\n== every figure (main and SI) has Source Data, and nothing extra is shipped ==")
figs = sorted(f for f in os.listdir(os.path.join(HERE, "figures")) if f.endswith(".png"))
srcs = sorted(os.listdir(SRC))
note(len(figs) == 8, f"eight figures present (4 main + 4 SI): {figs}")
for i in (1, 2, 3, 4):
    note(any(s.startswith(f"Figure{i}_") for s in srcs), f"Figure {i} has Source Data")
    note(any(s.startswith(f"SFigure{i}_") for s in srcs), f"SFigure {i} has Source Data")
expected = {"Figure1_source_data.csv", "Figure2_source_data.csv", "Figure3_source_data.csv",
            "Figure4_source_data.csv", "Figure4_scatter_data.csv",
            "SFigure1_source_data.csv", "SFigure2_source_data.csv", "SFigure3_source_data.csv",
            "SFigure4_source_data.csv"}
extra = set(srcs) - expected
note(not extra, f"no unexpected files in source_data (extra: {sorted(extra)})")

print("\n" + "=" * 70)
if FAIL:
    print(f"{len(FAIL)} inconsistencies between documents and Source Data")
    sys.exit(1)
print("documents and Source Data agree")
