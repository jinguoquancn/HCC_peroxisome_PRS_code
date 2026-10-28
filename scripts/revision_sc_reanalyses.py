# -*- coding: utf-8 -*-
"""PLOS ONE major-revision single-cell reanalyses (PONE-D-26-43453).

Addresses:
  Reviewer #1.1 - batch effects across the five scRNA-seq datasets
  Reviewer #2.7 - patient-level dependence of the malignant vs non-malignant
                  hepatocyte comparison (cells are not independent replicates)

Reads ONLY .obs (backed mode) from the annotated h5ad files, so it does not
recompute expression. Writes NEW files only:
  data/processed/revision_sc_patient_level.csv
  data/processed/revision_sc_batch.csv
  data/processed/revision_sc_variance.json
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, re, json, sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data", "processed")
SC = os.path.join(DATA, "sc")

DSETS = ["GSE149614", "GSE151530", "GSE125449", "GSE146115", "GSE98638"]
MAL, HEP = "Malignant", "Hepatocyte"


def patient_from_index(ds, idx):
    """Derive a patient/sample identifier when obs carries no sample column."""
    if ds == "GSE146115":               # HCC1_COL1_ROW1 -> HCC1
        return pd.Series([str(x).split("_")[0] for x in idx], index=idx)
    if ds == "GSE98638":               # PTH13 / TTH22 -> PTH / TTH
        return pd.Series([re.match(r"([A-Za-z]+)", str(x)).group(1) for x in idx], index=idx)
    return pd.Series([str(x) for x in idx], index=idx)


def load_obs(ds):
    import anndata as ad
    a = ad.read_h5ad(os.path.join(SC, f"{ds}_annot.h5ad"), backed="r")
    cols = [c for c in ["celltype", "pxs", "pxs_core", "sample", "Sample", "S_ID", "Type"]
            if c in a.obs.columns]
    o = a.obs[cols].copy()
    idx = o.index.astype(str)
    del a
    o["dataset"] = ds
    if "sample" in o.columns:
        o["patient"] = o["sample"].astype(str)
    elif "Sample" in o.columns:
        o["patient"] = o["Sample"].astype(str)
    else:
        o["patient"] = patient_from_index(ds, idx).values
    o["celltype"] = o["celltype"].astype(str)
    return o


def cohens_d(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan
    sp = np.sqrt(((nx - 1) * x.var(ddof=1) + (ny - 1) * y.var(ddof=1)) / (nx + ny - 2))
    return float((x.mean() - y.mean()) / sp) if sp > 0 else np.nan


def main():
    # ------------------------------------------------------------------ load
    frames = []
    for ds in DSETS:
        o = load_obs(ds)
        print(f"[load] {ds}: {len(o)} cells, {o.patient.nunique()} patients/samples",
              flush=True)
        frames.append(o)
    obs = pd.concat(frames, ignore_index=True)
    obs["pxs"] = pd.to_numeric(obs["pxs"], errors="coerce")
    obs = obs.dropna(subset=["pxs", "celltype"])

    # ---------------------------------------------------------- 1. batch table
    # (a) cell-level malignant vs hepatocyte, per dataset (as originally done)
    batch_rows = []
    for ds in DSETS:
        s = obs[obs.dataset == ds]
        m = s.loc[s.celltype == MAL, "pxs"].values
        h = s.loc[s.celltype == HEP, "pxs"].values
        if len(m) < 3 or len(h) < 3:
            batch_rows.append({"dataset": ds, "n_malignant": len(m), "n_hepatocyte": len(h),
                               "mean_malignant": np.nan, "mean_hepatocyte": np.nan,
                               "diff": np.nan, "cohens_d": np.nan, "mwu_p": np.nan,
                               "direction": "insufficient"})
            continue
        p = stats.mannwhitneyu(m, h).pvalue
        batch_rows.append({
            "dataset": ds, "n_malignant": len(m), "n_hepatocyte": len(h),
            "mean_malignant": float(np.mean(m)), "mean_hepatocyte": float(np.mean(h)),
            "diff": float(np.mean(m) - np.mean(h)), "cohens_d": cohens_d(m, h),
            "mwu_p": float(p),
            "direction": "malignant_lower" if np.mean(m) < np.mean(h) else "malignant_higher"})
    batch = pd.DataFrame(batch_rows)

    # (b) dataset-wide scale: is PXS comparable across datasets? (no normalisation)
    scale = obs.groupby("dataset")["pxs"].agg(["count", "mean", "std", "median"]).reset_index()
    scale.columns = ["dataset", "n_cells", "pxs_mean", "pxs_sd", "pxs_median"]
    batch = batch.merge(scale, on="dataset", how="left")
    batch.to_csv(os.path.join(DATA, "revision_sc_batch.csv"), index=False)

    # ------------------------------------------------- 2. patient-level analysis
    # within-dataset z-scoring removes any dataset-level offset/scale difference
    obs["pxs_z_within"] = obs.groupby("dataset")["pxs"].transform(
        lambda v: (v - v.mean()) / (v.std() + 1e-12))

    def patient_means(score_col):
        g = (obs[obs.celltype.isin([MAL, HEP])]
             .groupby(["dataset", "patient", "celltype"])[score_col]
             .agg(["mean", "count"]).reset_index())
        w = g.pivot(index=["dataset", "patient"], columns="celltype", values="mean")
        n = g.pivot(index=["dataset", "patient"], columns="celltype", values="count")
        w.columns = [f"mean_{c}" for c in w.columns]
        n.columns = [f"n_{c}" for c in n.columns]
        return w.join(n).reset_index()

    pat = patient_means("pxs")
    pat.to_csv(os.path.join(DATA, "revision_sc_patient_level.csv"), index=False)

    # paired test across patients (within dataset) for both raw and z-scored
    paired_rows = []
    tabs = {}
    for tag, col in [("raw", "pxs"), ("z", "pxs_z_within")]:
        t = patient_means(col)
        t = t.dropna(subset=[f"mean_{MAL}", f"mean_{HEP}"])
        t["diff"] = t[f"mean_{MAL}"] - t[f"mean_{HEP}"]
        tabs[tag] = t
        # per dataset paired test
        for ds in DSETS:
            s = t[t.dataset == ds]
            if len(s) >= 3:
                w = stats.wilcoxon(s[f"mean_{MAL}"], s[f"mean_{HEP}"])
                pv = float(w.pvalue)
            else:
                pv = np.nan
            paired_rows.append({
                "score": "PXS" if tag == "raw" else "PXS (z within dataset)",
                "dataset": ds,
                "n_patients_paired": int(len(s)),
                "mean_diff_malignant_minus_hepatocyte": float(s["diff"].mean()) if len(s) else np.nan,
                "n_patients_malignant_lower": int((s["diff"] < 0).sum()) if len(s) else 0,
                "n_patients_malignant_higher": int((s["diff"] > 0).sum()) if len(s) else 0,
                "wilcoxon_p": pv})
        # pooled across datasets
        if tag == "raw":
            pooled = t
            w = stats.wilcoxon(pooled[f"mean_{MAL}"], pooled[f"mean_{HEP}"])
            pooled_p = float(w.pvalue)
            pooled_n = int(len(pooled))
            pooled_diff = float(pooled["diff"].mean())
            pooled_lower = int((pooled["diff"] < 0).sum())
            pooled_higher = int((pooled["diff"] > 0).sum())
    paired = pd.DataFrame(paired_rows)
    paired.to_csv(os.path.join(DATA, "revision_sc_patient_tests.csv"), index=False)

    # ------------------------------------------------- 3. variance decomposition
    sub = obs[obs.celltype.str.strip() != ""].copy()
    y = sub["pxs"].values
    sst = float(np.sum((y - y.mean()) ** 2))

    def eta2(cols):
        X = pd.get_dummies(sub[cols].astype(str), drop_first=False).astype(float).values
        X = np.column_stack([np.ones(len(X)), X])
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        ssr = float(np.sum(resid ** 2))
        return float((sst - ssr) / sst) if sst > 0 else np.nan

    var = {
        "n_cells": int(len(sub)),
        "n_datasets": int(sub.dataset.nunique()),
        "n_patients": int(sub.groupby("dataset").patient.nunique().sum()),
        "R2_dataset_only": round(eta2(["dataset"]), 4),
        "R2_celltype_only": round(eta2(["celltype"]), 4),
        "R2_dataset_plus_celltype": round(eta2(["dataset", "celltype"]), 4),
        "R2_dataset_plus_celltype_plus_patient": round(
            eta2(["dataset", "celltype", "patient"]), 4),
        "pooled_malignant_mean_PXS": round(float(
            obs.loc[obs.celltype == MAL, "pxs"].mean()), 4),
        "pooled_hepatocyte_mean_PXS": round(float(
            obs.loc[obs.celltype == HEP, "pxs"].mean()), 4),
        "patient_level_pooled": {
            "n_patients_paired": pooled_n,
            "mean_diff_malignant_minus_hepatocyte": round(pooled_diff, 4),
            "n_patients_malignant_lower": pooled_lower,
            "n_patients_malignant_higher": pooled_higher,
            "wilcoxon_p": pooled_p},
        "n_patients_by_dataset": obs.groupby("dataset").patient.nunique().to_dict(),
    }
    json.dump(var, open(os.path.join(DATA, "revision_sc_variance.json"), "w"), indent=1)

    # ------------------------------------------------------------- reporting
    print("\n===== PER-DATASET (cell-level) malignant vs hepatocyte PXS =====")
    print(batch[["dataset", "n_malignant", "n_hepatocyte", "mean_malignant",
                 "mean_hepatocyte", "diff", "cohens_d", "mwu_p", "direction",
                 "pxs_mean", "pxs_sd"]].to_string(index=False))
    print("\n===== PATIENT-LEVEL paired (within-patient) =====")
    print(paired.to_string(index=False))
    print("\n===== VARIANCE DECOMPOSITION =====")
    print(json.dumps({k: v for k, v in var.items() if k != "n_patients_by_dataset"}, indent=1))
    print("\nn patients by dataset:", var["n_patients_by_dataset"])
    print("\nrevision_sc_* saved")


if __name__ == "__main__":
    main()
