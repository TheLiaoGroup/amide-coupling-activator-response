"""Extract the acylation-only, self-contained dataset from the in-house corpus.

Single source of truth for `data/DATASET_acylation_wells.csv` and
`data/DATASET_acylation_conditions.xlsx`. Every other script in this project (REPRODUCE_acylation.py,
build_figures.py, build_null_figure.py) reads only these two files - never the multi-transformation
corpus - so the code/data package is self-contained and does not require a separate ~58 MB download.

Two tables, two readers:
  - `DATASET_acylation_wells.csv`      - one row per well. The reaction condition is given only as
    `condition_id`; look up what that condition contains in the other file.
  - `DATASET_acylation_conditions.xlsx` - one row per condition_id (94 rows): activator, additive and
    base by CAS/name/equivalents, plus solvent, temperature and time (constant across all 94 - DMF,
    25 degrees C, 2 h - stated once here rather than repeated 53,241 times in the wells table).

Columns dropped from the corpus's own `DATASET_analysis.csv`/`DATASET_conditions.csv` (`screen`,
`rxn_class`, `transformation`, the per-well `reagent1_id..reagent6_id`, `solvent_id`, `temp_id`,
`time_id`) are either constant for this screen or fully determined by `condition_id`; nothing about
a well is lost, since every dropped value is recoverable from the two tables here.

    <python> make_dataset.py

Reads the three released files from $HTE_DATA (falling back to the path this analysis was
originally run from). Run once; commit the two output files to `data/`.
"""
import os
import pandas as pd

_here = os.path.dirname(os.path.abspath(__file__))
_orig = os.path.join("D:" + os.sep, "Google Drive", "ManuscriptForSubmission", "HTE data",
                     "SUBMISSION", "06_Data_and_Code")
DATA_IN = os.environ.get("HTE_DATA") or _orig
DATA_OUT = os.path.join(_here, "data")
os.makedirs(DATA_OUT, exist_ok=True)

print("== reading the released corpus ==")
an = pd.read_csv(os.path.join(DATA_IN, "DATASET_analysis.csv"), low_memory=False)
cond = pd.read_csv(os.path.join(DATA_IN, "DATASET_conditions.csv"), low_memory=False)
cmpd = pd.read_csv(os.path.join(DATA_IN, "DATASET_compounds.csv"), low_memory=False)
cas_of = dict(zip(cmpd["Compound ID"], cmpd["cas"]))

raw = an[an.screen == "Acylation|1"].copy()
assert len(raw) == 53241, f"expected 53241 acylation wells, got {len(raw)}"

# ---------------------------------------------------------------- wells table
WELLS_COLS = ["plate", "condition_id", "subpair", "sub_1_smiles", "sub_2_smiles", "product_smiles",
              "sub_1_FG", "sub_2_FG", "yield", "yield_clf", "check_passed", "allocation"]
wells = raw[WELLS_COLS].copy()
wells_path = os.path.join(DATA_OUT, "DATASET_acylation_wells.csv")
wells.to_csv(wells_path, index=False)
print(f"wrote {wells_path}: {len(wells)} rows, {len(wells.columns)} columns, "
      f"{wells.condition_id.nunique()} distinct conditions")

# ---------------------------------------------------------------- conditions table
cA = cond[(cond.reaction_type == "acylation") & (cond.condition_id.isin(raw.condition_id.unique()))]
assert set(cA.condition_id) == set(raw.condition_id.unique()), "wells and conditions tables disagree"
assert len(cA) == 94, f"expected 94 acylation conditions, got {len(cA)}"
# constant across every acylation condition - checked, not assumed
assert (cA["Standardized Name"] == "DMF").all()
assert (cA.solvent_id == "CMP0552").all()
assert (cA.temp_id == 25.0).all() and (cA.time_id == 2.0).all()
assert (cA.Reagent4_id.isna() & cA.Reagent5_id.isna() & cA.Reagent6_id.isna()).all(), \
    "a 4th/5th/6th reagent slot is used - the activator/additive/base role mapping below is incomplete"

rows = []
for _, r in cA.iterrows():
    rows.append({
        "condition_id": r.condition_id,
        "activator_CAS": r.Reagent1_HTE_CAS, "activator_name": r["Reagent1_Standardized Name"],
        "activator_equiv": r.Reagent1_HTE_Equiv,
        "additive_CAS": r.Reagent2_HTE_CAS, "additive_name": r["Reagent2_Standardized Name"],
        "additive_equiv": r.Reagent2_HTE_Equiv,
        "base_CAS": r.Reagent3_HTE_CAS, "base_name": r["Reagent3_Standardized Name"],
        "base_equiv": r.Reagent3_HTE_Equiv,
        "solvent_CAS": cas_of.get(r.solvent_id), "solvent_name": r["Standardized Name"],
        "temp_C": r.temp_id, "time_h": r.time_id,
        "note_en": r.note_en if pd.notna(r.note_en) else "",
    })
conditions = pd.DataFrame(rows).sort_values("condition_id",
    key=lambda s: s.str.extract(r"(\d+)")[0].astype(int)).reset_index(drop=True)

conditions_path = os.path.join(DATA_OUT, "DATASET_acylation_conditions.xlsx")
with pd.ExcelWriter(conditions_path, engine="openpyxl") as xl:
    conditions.to_excel(xl, sheet_name="conditions", index=False)
    ws = xl.sheets["conditions"]
    ws.freeze_panes = "B2"
    ws.auto_filter.ref = ws.dimensions
    for i, col in enumerate(conditions.columns, start=1):
        width = max(len(col), conditions[col].astype(str).map(len).max()) + 2
        ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = min(width, 40)
print(f"wrote {conditions_path}: {len(conditions)} rows, {len(conditions.columns)} columns")

# ---------------------------------------------------------------- read back and print, per CLAUDE.md
print("\n== read back ==")
print(pd.read_csv(wells_path, nrows=3).to_string())
print()
print(pd.read_excel(conditions_path).head(3).to_string())
