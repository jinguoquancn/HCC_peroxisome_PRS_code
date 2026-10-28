# -*- coding: utf-8 -*-
"""Reviewer 1, point 2: multivariable model including tumour stage.

Fits three Cox models on the TCGA-LIHC training cohort:
  M1  age + sex                         (demographic baseline)
  M2  age + sex + AJCC stage            (clinical baseline)
  M3  age + sex + AJCC stage + PRS      (added value of the PRS)
and reports Harrell C-index, the likelihood-ratio test for adding the PRS,
scaled Schoenfeld residual tests, and an alternative model that uses
Child-Pugh instead of stage (sensitivity).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import json
import os
import sys

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data", "processed")
CLIN = os.path.join(
    ROOT, "data", "raw", "extracted_tcga",
    "gdac.broadinstitute.org_LIHC.Merge_Clinical.Level_1.2016012800.0.0",
    "LIHC.merged_only_clinical_clin_format.txt")

import core  # noqa: E402
from analysis_tcga import load  # noqa: E402

STAGE_MAP = {"stage i": 1, "stage ii": 2, "stage iii": 3, "stage iiia": 3,
             "stage iiib": 3, "stage iiic": 3, "stage iv": 4, "stage iva": 4,
             "stage ivb": 4}


def load_stage():
    d = pd.read_csv(CLIN, sep="\t", header=None, low_memory=False)
    attr = d[0].astype(str).str.strip()
    body = d.drop(columns=[0])
    bc = body[attr == "patient.bcr_patient_barcode"].iloc[0].values
    st = body[attr == "patient.stage_event.pathologic_stage"].iloc[0].values
    s = pd.Series(st, index=[str(b).lower().strip() for b in bc])
    s = s.map(lambda x: STAGE_MAP.get(str(x).lower().strip(), np.nan))
    return s[~s.index.duplicated()]


def fit(df, cols, label):
    d = df[["T", "E"] + cols].dropna()
    cph = CoxPHFitter(penalizer=1e-6)
    cph.fit(d, duration_col="T", event_col="E")
    return cph, d, label

def main():
    E, cl, _ = load()
    coef = pd.read_csv(os.path.join(DATA, "TCGA_PRS_coefficients.csv"), index_col=0).iloc[:, 0]
    prs = core.risk_score(E, coef)

    stage = load_stage()
    df = pd.DataFrame({
        "T": cl["os_years"], "E": cl["os_event"].astype(int),
        "age_sd": (cl["age"] - cl["age"].mean()) / cl["age"].std(),
        "male": (cl["gender"].astype(str).str.lower() == "male").astype(float),
        "stage": cl.index.to_series().str.lower().str[:12].map(stage),
        "prs_sd": (prs - prs.mean()) / prs.std(),
    })
    df["stage"] = pd.to_numeric(df["stage"], errors="coerce")
    print("stage available:", df["stage"].notna().sum(), "of", len(df))
    print(df["stage"].value_counts().sort_index().to_string())

    out = {}
    models = {
        "M1_demographic": ["age_sd", "male"],
        "M2_clinical_stage": ["age_sd", "male", "stage"],
        "M3_clinical_stage_PRS": ["age_sd", "male", "stage", "prs_sd"],
    }
    summ = {}
    for name, cols in models.items():
        cph, d, _ = fit(df, cols, name)
        c = float(cph.concordance_index_)
        rows = []
        for term in cols:
            rows.append({"model": name, "term": term,
                         "HR": float(np.exp(cph.params_[term])),
                         "CI_low": float(np.exp(cph.confidence_intervals_.loc[term].iloc[0])),
                         "CI_high": float(np.exp(cph.confidence_intervals_.loc[term].iloc[1])),
                         "p": float(cph.summary.loc[term, "p"])})
        summ[name] = {"n": int(len(d)), "events": int(d["E"].sum()),
                      "c_index": c, "terms": rows}
        print(f"\n{name}: n={len(d)} events={int(d['E'].sum())} C={c:.3f}")
        print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    # likelihood-ratio test for adding PRS to the clinical model
    m2 = CoxPHFitter(penalizer=1e-6).fit(
        df[["T", "E"] + models["M2_clinical_stage"]].dropna(), "T", "E")
    d3 = df[["T", "E"] + models["M3_clinical_stage_PRS"]].dropna()
    m3 = CoxPHFitter(penalizer=1e-6).fit(d3, "T", "E")
    lrt = 2 * (m3.log_likelihood_ - m2.log_likelihood_)
    from scipy.stats import chi2
    p_lrt = float(chi2.sf(lrt, 1))
    print(f"\nLRT add PRS: chi2={lrt:.2f} p={p_lrt:.3g}")

    # proportional hazards on the full model
    ph = proportional_hazard_test(m3, d3, time_transform="rank")
    ph_rows = [{"term": t, "p": float(ph.summary.loc[t, "p"])} for t in ph.summary.index]
    print(pd.DataFrame(ph_rows).to_string(index=False))

    # sensitivity: stage dichotomised as advanced (III/IV) versus early (I/II)
    df["stage_adv"] = (df["stage"] >= 3).astype(float).where(df["stage"].notna())
    m4 = CoxPHFitter(penalizer=1e-6).fit(
        df[["T", "E", "age_sd", "male", "stage_adv", "prs_sd"]].dropna(), "T", "E")
    alt = {t: {"HR": float(np.exp(m4.params_[t])),
               "CI_low": float(np.exp(m4.confidence_intervals_.loc[t].iloc[0])),
               "CI_high": float(np.exp(m4.confidence_intervals_.loc[t].iloc[1])),
               "p": float(m4.summary.loc[t, "p"])}
           for t in m4.params_.index}
    print("\nAdvanced-stage sensitivity:", alt, "C=", round(m4.concordance_index_, 3))

    # sensitivity: Child-Pugh instead of stage (guard: NaN must stay NaN)
    cp = cl["child_pugh"].str.extract(r"^([abc])", expand=False)
    df["cp"] = cp.map({"a": 1, "b": 2, "c": 3})
    cs = CoxPHFitter(penalizer=1e-6).fit(
        df[["T", "E", "age_sd", "male", "cp", "prs_sd"]].dropna(), "T", "E")
    cp_summary = {t: {"HR": float(np.exp(cs.params_[t])),
                      "CI_low": float(np.exp(cs.confidence_intervals_.loc[t].iloc[0])),
                      "CI_high": float(np.exp(cs.confidence_intervals_.loc[t].iloc[1])),
                      "p": float(cs.summary.loc[t, "p"])} for t in cs.params_.index}
    print("\nChild-Pugh sensitivity:", cp_summary, "C=", round(cs.concordance_index_, 3))

    json.dump({"stage_models": summ,
               "lrt_add_PRS": {"chi2": float(lrt), "df": 1, "p": p_lrt},
               "ph_test_rank": ph_rows,
               "advanced_stage_sensitivity": {"n": int(m4._n_examples),
                                              "c_index": float(m4.concordance_index_),
                                              "terms": alt},
               "child_pugh_sensitivity": {"n": int(cs._n_examples),
                                          "c_index": float(cs.concordance_index_),
                                          "terms": cp_summary},
               "stage_counts": {int(k): int(v) for k, v in
                                df["stage"].value_counts().sort_index().items()},
               "n_stage_missing": int(df["stage"].isna().sum())},
              open(os.path.join(DATA, "revision_multivariable_stage.json"), "w"), indent=1)
    rows = [r for m in summ.values() for r in m["terms"]]
    pd.DataFrame(rows).to_csv(os.path.join(DATA, "revision_multivariable_stage.csv"), index=False)
    print("\nwrote revision_multivariable_stage.json / .csv")


if __name__ == "__main__":
    main()
