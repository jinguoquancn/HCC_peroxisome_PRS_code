# -*- coding: utf-8 -*-
"""External validation of the PRS across GEO cohorts + ICGC-LIRI-JP.
Computes risk scores, KM split, HR, C-index, time-dependent ROC per cohort."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, sys, json, warnings
import numpy as np
import pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")

sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
import core

DATA = str(PROJECT_ROOT / "data" / "processed")
rng = np.random.default_rng(20260829)


def load_cohort_expr(name):
    df = pd.read_csv(os.path.join(DATA, f"{name}_expression.csv"), index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    return np.log2(df + 1)


def load_cohort_clin(name):
    f = os.path.join(DATA, f"{name}_clinical.csv")
    if not os.path.exists(f):
        return None
    df = pd.read_csv(f)
    return df


def map_survival(df, cohort):
    """Normalize survival columns -> time (years), event."""
    t, e = None, None
    if cohort == "GSE14520":
        t = df["os_months"] / 12
        e = df["os_status"]
    elif cohort == "GSE76427":
        t = df["os_months"] / 12
        e = df["os_event"]
    elif cohort == "GSE10186":
        t = df["os_months"] / 12
        e = df["os_event"]
    elif cohort == "GSE116174":
        if "os_months" in df.columns:
            t = pd.to_numeric(df["os_months"], errors="coerce") / 12
            e = pd.to_numeric(df["os_event"], errors="coerce")
    return t, e


def main():
    prs = pd.read_csv(os.path.join(DATA, "TCGA_PRS_coefficients.csv"), index_col=0)["0"]
    prs = prs.dropna()
    print("PRS genes:", len(prs), flush=True)

    cohorts = ["GSE14520", "GSE76427", "GSE116174", "GSE10186",
               "GSE36376", "GSE25097", "GSE54236", "GSE87630"]
    results = []
    for name in cohorts:
        ep = os.path.join(DATA, f"{name}_expression.csv")
        if not os.path.exists(ep):
            print("[skip]", name, flush=True)
            continue
        expr = load_cohort_expr(name)
        # tumor samples only: GSE14520/76427/10186 etc. use all (matrices are tumor)
        risk = core.risk_score(expr, prs)
        clin = load_cohort_clin(name)
        if clin is None or name not in ("GSE14520", "GSE76427", "GSE10186", "GSE116174"):
            # store risk scores for expression-level analyses
            results.append({"cohort": name, "n": len(risk), "risk_available": True,
                            "survival": False})
            risk.to_csv(os.path.join(DATA, f"{name}_risk.csv"))
            print(f"== {name}: risk scores for {len(risk)} samples (no survival)", flush=True)
            continue
        t, e = map_survival(clin, name)
        if t is None:
            results.append({"cohort": name, "n": len(risk), "risk_available": True, "survival": False})
            continue
        # align clinical to expression samples
        if "sid" in clin.columns:
            key_col = "sid"
        elif "sample" in clin.columns:
            key_col = "sample"
        elif "Affy_GSM" in clin.columns:
            key_col = "Affy_GSM"
        else:
            key_col = clin.columns[0]
        clin["sid"] = clin[key_col].astype(str)
        if name == "GSE14520":
            # strip leading 'GSM'? already there; expression cols are GSM ids
            pass
        risk = risk.reindex(clin["sid"].values)
        d = pd.DataFrame({"risk": risk.values, "t": np.asarray(t, float), "e": np.asarray(e, float)})
        d = d.dropna()
        if name == "GSE14520":
            d = d[d["t"] > 0]
        if len(d) < 30:
            results.append({"cohort": name, "n": len(d), "survival": True, "km": False})
            continue
        med = np.median(d["risk"])
        grp = (d["risk"] > med).astype(int)
        hr, lo, hi, p, n = core.cox_hr(d["t"], d["e"], grp)
        c, cn = core.cindex(d["t"], d["e"], d["risk"].values)
        auc1, auc3, auc5 = (core.td_auc(d["t"], d["e"], d["risk"].values, th) for th in (1, 3, 5))
        results.append({"cohort": name, "n": len(d), "survival": True, "km": True,
                        "HR": hr, "HR_low": lo, "HR_high": hi, "p_logrank": p,
                        "C_index": c, "AUC_1y": auc1, "AUC_3y": auc3, "AUC_5y": auc5,
                        "median_risk": med})
        d.to_csv(os.path.join(DATA, f"{name}_validation_surv.csv"))
        print(f"== {name}: n={len(d)} HR={hr:.2f} [{lo:.2f}-{hi:.2f}] p={p:.4f} C={c:.3f} AUC1/3/5={auc1:.2f}/{auc3:.2f}/{auc5:.2f}", flush=True)

    # ICGC
    ep = os.path.join(DATA, "ICGC_LIRI_expression_tumor.csv")
    cp = os.path.join(DATA, "ICGC_LIRI_clinical_built.csv")
    if os.path.exists(ep) and os.path.exists(cp):
        expr = pd.read_csv(ep, index_col=0)
        clin = pd.read_csv(cp)
        risk = core.risk_score(expr, prs)
        d = pd.DataFrame({"specimen": clin["specimen"], "t": clin["os_time_days"] / 365.25,
                          "e": clin["os_event"], "age": clin["age"], "stage": clin["stage"]})
        risk = risk.reindex(d["specimen"].values)
        d["risk"] = risk.values
        d = d.dropna(subset=["risk", "t"])
        if len(d) > 30:
            med = np.median(d["risk"])
            grp = (d["risk"] > med).astype(int)
            hr, lo, hi, p, n = core.cox_hr(d["t"], d["e"], grp)
            c, cn = core.cindex(d["t"], d["e"], d["risk"].values)
            auc1, auc3, auc5 = (core.td_auc(d["t"], d["e"], d["risk"].values, th) for th in (1, 3, 5))
            results.append({"cohort": "ICGC-LIRI-JP", "n": len(d), "survival": True, "km": True,
                            "HR": hr, "HR_low": lo, "HR_high": hi, "p_logrank": p,
                            "C_index": c, "AUC_1y": auc1, "AUC_3y": auc3, "AUC_5y": auc5,
                            "median_risk": med})
            d.to_csv(os.path.join(DATA, "ICGC_LIRI_validation_surv.csv"))
            print(f"== ICGC-LIRI-JP: n={len(d)} HR={hr:.2f} p={p:.4f} C={c:.3f}", flush=True)

    res = pd.DataFrame(results)
    res.to_csv(os.path.join(DATA, "external_validation_summary.csv"), index=False)
    print(res.to_string(), flush=True)
    print("VALIDATION DONE", flush=True)


if __name__ == "__main__":
    main()
