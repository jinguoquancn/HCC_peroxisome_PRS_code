# -*- coding: utf-8 -*-
"""Build ICGC-LIRI-JP cohort: expression (tumor) + donor survival clinical."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os
import numpy as np
import pandas as pd

RAW = str(PROJECT_ROOT / "data" / "raw")
OUT = str(PROJECT_ROOT / "data" / "processed")

expr = pd.read_csv(os.path.join(RAW, "ICGC_LIRI_expression.csv"), index_col=0)
ann = pd.read_csv(os.path.join(RAW, "ICGC_LIRI_sample_annotation.csv"))
clin = pd.read_csv(os.path.join(RAW, "ICGC_LIRI_clinical.csv"))

# primary tumor specimens
tum = ann[ann["specimen_type"].str.contains("Primary tumour", na=False)]
print("primary tumour specimens:", len(tum), flush=True)
tum_ids = tum["icgc_specimen_id"].astype(str)
expr_t = expr[[c for c in expr.columns if c in set(tum_ids)]] if "gene_id" not in expr.columns else expr[[c for c in expr.columns if c != "gene_id"]][[c for c in expr.columns[1:] if c in set(tum_ids)]]
# expr columns are icgc_specimen_id
spec = [c for c in expr.columns if c in set(tum_ids)]
expr_t = expr[spec]
print("expr_t:", expr_t.shape, flush=True)

# map specimen -> donor
donor_map = tum.set_index("icgc_specimen_id")["icgc_donor_id"].to_dict()
donors = [donor_map.get(s, np.nan) for s in spec]
clin.set_index("icgc_donor_id", inplace=True)

rows = []
for s, d in zip(spec, donors):
    if pd.isna(d):
        continue
    c = clin.loc[d]
    rows.append({
        "specimen": s, "donor": d,
        "os_time_days": float(c["donor_survival_time"]) if not pd.isna(c["donor_survival_time"]) else np.nan,
        "os_event": 1 if str(c["donor_vital_status"]).lower() == "deceased" else 0,
        "age": float(c["donor_age_at_diagnosis"]) if not pd.isna(c["donor_age_at_diagnosis"]) else np.nan,
        "sex": c["donor_sex"],
        "stage": c["donor_tumour_stage_at_diagnosis"],
    })
cl = pd.DataFrame(rows)
cl["os_months"] = cl["os_time_days"] / 30.44
print("ICGC survival table:", cl.shape, "events:", cl["os_event"].sum(), flush=True)
cl.to_csv(os.path.join(OUT, "ICGC_LIRI_clinical_built.csv"), index=False)
expr_t.to_csv(os.path.join(OUT, "ICGC_LIRI_expression_tumor.csv"))
print("ICGC DONE", flush=True)
