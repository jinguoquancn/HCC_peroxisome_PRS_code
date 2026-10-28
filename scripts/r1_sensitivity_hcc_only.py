# -*- coding: utf-8 -*-
"""Reviewer 1, point 1: HCC-only sensitivity analysis for GSE151530.

GSE151530 (Zhang et al., Gastroenterology 2021) deposits 46 samples: 32 HCC
labelled H08-H77 and 14 intrahepatic cholangiocarcinoma labelled C25-C76.
The primary analysis used the series as deposited; this script repeats the
malignant-versus-non-malignant comparison restricted to the HCC samples and
writes the comparison table.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

SC = str(PROJECT_ROOT / "data" / "processed" / "sc")
OUT = str(PROJECT_ROOT / "data" / "processed")


def patient_level(df, label):
    """df: columns patient, mean_hep, mean_mal, n_hep, n_mal."""
    d = df[(df["n_hep"] > 0) & (df["n_mal"] > 0)].copy()
    diff = d["mean_mal"] - d["mean_hep"]
    n_low = int((diff < 0).sum())
    n_high = int((diff > 0).sum())
    if len(d) >= 5:
        p = float(stats.wilcoxon(diff).pvalue)
    else:
        p = float("nan")
    return {"subset": label, "n_patients_paired": int(len(d)),
            "mean_diff_malignant_minus_hepatocyte": float(diff.mean()) if len(d) else float("nan"),
            "n_patients_malignant_lower": n_low,
            "n_patients_malignant_higher": n_high,
            "wilcoxon_p": p}


def main():
    import anndata as ad
    a = ad.read_h5ad(os.path.join(SC, "GSE151530_annot.h5ad"), backed="r")
    obs = a.obs[["S_ID", "Sample", "celltype", "pxs"]].copy()
    obs["histology"] = np.where(obs["Sample"].astype(str).str.startswith("C"),
                                "iCCA", "HCC")
    print(obs.groupby("histology")["Sample"].nunique().to_string())
    print(obs["celltype"].value_counts().to_string())

    rows, cell_rows = [], []
    for hist, sub in obs.groupby("histology"):
        mal = sub[sub["celltype"] == "Malignant"]["pxs"].astype(float)
        hep = sub[sub["celltype"] == "Hepatocyte"]["pxs"].astype(float)
        if len(mal) > 2 and len(hep) > 2:
            d = (mal.mean() - hep.mean()) / np.sqrt(
                ((len(mal) - 1) * mal.var() + (len(hep) - 1) * hep.var())
                / (len(mal) + len(hep) - 2))
            p = stats.mannwhitneyu(mal, hep).pvalue
        else:
            d, p = np.nan, np.nan
        cell_rows.append({"subset": hist, "n_malignant": int(len(mal)),
                          "n_hepatocyte": int(len(hep)),
                          "mean_PXS_malignant": float(mal.mean()) if len(mal) else np.nan,
                          "mean_PXS_hepatocyte": float(hep.mean()) if len(hep) else np.nan,
                          "cohens_d": float(d), "mwu_p_cell_level": float(p)})

        g = sub.groupby(["Sample", "celltype"])["pxs"].agg(["mean", "size"]).reset_index()
        wide = g.pivot(index="Sample", columns="celltype", values="mean")
        cnt = g.pivot(index="Sample", columns="celltype", values="size")
        t = pd.DataFrame({
            "patient": wide.index,
            "mean_hep": wide.get("Hepatocyte"),
            "mean_mal": wide.get("Malignant"),
            "n_hep": cnt.get("Hepatocyte"),
            "n_mal": cnt.get("Malignant")}).fillna({"n_hep": 0, "n_mal": 0})
        rows.append(patient_level(t, hist))

    ct = pd.DataFrame(cell_rows)
    pt = pd.DataFrame(rows)
    ct.to_csv(os.path.join(OUT, "revision_sc_hcc_only_cells.csv"), index=False)
    pt.to_csv(os.path.join(OUT, "revision_sc_hcc_only_patient.csv"), index=False)
    print(ct.to_string(index=False))
    print(pt.to_string(index=False))
    json.dump({"cells": ct.to_dict("records"), "patients": pt.to_dict("records")},
              open(os.path.join(OUT, "revision_sc_hcc_only.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
