# Data files

Two files, both complete and self-contained - no external download is needed to reproduce this
manuscript. They are the acylation subset only of a larger in-house dataset; the extraction script
(`make_dataset.py`) is included for provenance.

| file | rows | role |
|---|---|---|
| `DATASET_acylation_wells.csv` | 53,241 | one row per well: plate, `condition_id`, substrate pair, SMILES, functional groups, yield, QC flag, train/test allocation |
| `DATASET_acylation_conditions.xlsx` | 94 | one row per `condition_id`: activator / additive / base by CAS, name and equivalents, solvent, temperature, time |

The reaction condition in the wells table is given only as `condition_id`; look up what that
condition contains in the conditions file. Solvent (DMF), temperature (25 degrees C) and time (2 h)
are constant across all 94 conditions and are stated once there rather than repeated 53,241 times.

