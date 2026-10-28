# -*- coding: utf-8 -*-
"""PLOS ONE major-revision bulk reanalyses, part 2 (PONE-D-26-43453).

Addresses:
  Reviewer #2.4 - cohort-specific performance with CIs, event counts, follow-up,
                  meta-analysis heterogeneity
  Reviewer #2.5 - proportional hazards assumption; documented missing data / events
  Reviewer #2.3 - confounding by tumor purity / cellular composition
  Reviewer #2.1 - exact N and exclusions per analysis

Writes NEW files only:
  data/processed/revision_external_rigor.csv
  data/processed/revision_meta_heterogeneity.json
  data/processed/revision_ph_test.csv
  data/processed/revision_purity_confounding.csv
  data/processed/revision_n_table.csv
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, sys, json

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data", "processed")

SEED = 42
COHORTS = ["GSE14520", "GSE76427", "GSE116174", "GSE10186", "ICGC-LIRI-JP"]


def ci_bootstrap(t, e, risk, n_boot=1000, seed=SEED):
    """Bootstrap CI for the C-index (positional arrays -> no index reuse trap)."""
    from lifelines.utils import concordance_index
    t = np.asarray(t, float)
    e = np.asarray(e, float)
    r = np.asarray(risk, float)
    ok = np.isfinite(t) & np.isfinite(e) & np.isfinite(r)
    t, e, r = t[ok], e[ok], r[ok]
    # NOTE: lifelines expects predicted *survival* times (higher = longer survival),
    # so a risk score (higher = shorter survival) must be negated.
    base = concordance_index(t, -r, e)
    rng = np.random.default_rng(seed)
    n = len(t)
    bs = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        tt, ee, rr = t[i], e[i], r[i]
        if len(np.unique(ee)) < 2:
            continue
        try:
            bs.append(concordance_index(tt, -rr, ee))
        except Exception:
            continue
    bs = np.asarray(bs, float)
    bs = bs[np.isfinite(bs)]
    # degenerate-check (see lesson: silent 0.5 / flat distribution)
    if len(bs) < 50 or np.ptp(bs) < 0.02:
        return base, np.nan, np.nan, len(bs)
    return base, float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), len(bs)


def median_followup(t, e):
    """Reverse Kaplan-Meier median follow-up."""
    from lifelines import KaplanMeierFitter
    t = np.asarray(t, float)
    e = np.asarray(e, float)
    ok = np.isfinite(t) & np.isfinite(e)
    kmf = KaplanMeierFitter().fit(t[ok], event_observed=1 - e[ok])
    m = kmf.median_survival_time_
    return float(m) if np.isfinite(m) else np.nan


def dl_meta(hrs, los, his):
    """DerSimonian-Laird random-effects pooling of log HRs."""
    y = np.log(np.asarray(hrs, float))
    se = (np.log(np.asarray(his, float)) - np.log(np.asarray(los, float))) / (2 * 1.96)
    w = 1 / se ** 2
    y_fixed = np.sum(w * y) / np.sum(w)
    Q = float(np.sum(w * (y - y_fixed) ** 2))
    k = len(y)
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    wr = 1 / (se ** 2 + tau2)
    y_re = np.sum(wr * y) / np.sum(wr)
    se_re = float(np.sqrt(1 / np.sum(wr)))
    z = y_re / se_re
    return {"k": k, "pooled_HR": float(np.exp(y_re)),
            "CI_low": float(np.exp(y_re - 1.96 * se_re)),
            "CI_high": float(np.exp(y_re + 1.96 * se_re)),
            "p": float(2 * stats.norm.sf(abs(z))),
            "Q": round(Q, 3), "I2_percent": round(I2, 1),
            "tau2": round(tau2, 5)}


def external_rigor():
    rows = []
    for c in COHORTS:
        f = os.path.join(DATA, f"{c}_validation_surv.csv")
        if c == "ICGC-LIRI-JP":
            f = os.path.join(DATA, "ICGC_LIRI_validation_surv.csv")
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f)
        d = d.rename(columns={d.columns[0]: "id"})
        t, e, r = d["t"].values, d["e"].values, d["risk"].values
        ok = np.isfinite(t) & np.isfinite(e) & np.isfinite(r) & (t > 0)
        t, e, r = t[ok], e[ok], r[ok]
        ev = int(np.sum(e))
        # continuous (per SD of z-standardised risk ~ per unit)
        from lifelines import CoxPHFitter
        df = pd.DataFrame({"T": t, "E": e, "R": r}).reset_index(drop=True)
        cph = CoxPHFitter(penalizer=1e-6).fit(df, "T", "E")
        hr = float(np.exp(cph.params_["R"]))
        lo = float(np.exp(cph.confidence_intervals_.loc["R"].iloc[0]))
        hi = float(np.exp(cph.confidence_intervals_.loc["R"].iloc[1]))
        # median split
        med = np.median(r)
        g = (r > med).astype(int)
        from lifelines.statistics import logrank_test
        lr = logrank_test(t[g == 0], t[g == 1], e[g == 0], e[g == 1])
        dkm = pd.DataFrame({"T": t, "E": e, "G": g})
        cph2 = CoxPHFitter(penalizer=1e-6).fit(dkm, "T", "E")
        hr2 = float(np.exp(cph2.params_["G"]))
        lo2 = float(np.exp(cph2.confidence_intervals_.loc["G"].iloc[0]))
        hi2 = float(np.exp(cph2.confidence_intervals_.loc["G"].iloc[1]))
        ci, clo, chi, nb = ci_bootstrap(t, e, r)
        rows.append({"cohort": c, "n": int(len(t)), "n_events": ev,
                     "median_followup_years": round(median_followup(t, e), 2),
                     "HR_per_unit_continuous": round(hr, 3),
                     "HR_cont_CI_low": round(lo, 3), "HR_cont_CI_high": round(hi, 3),
                     "HR_median_split": round(hr2, 3),
                     "HR_split_CI_low": round(lo2, 3), "HR_split_CI_high": round(hi2, 3),
                     "logrank_p_median_split": float(lr.p_value),
                     "C_index": round(ci, 3), "C_index_CI_low": round(clo, 3),
                     "C_index_CI_high": round(chi, 3), "n_boot_ok": nb})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(DATA, "revision_external_rigor.csv"), index=False)
    meta = dl_meta(df["HR_median_split"], df["HR_split_CI_low"], df["HR_split_CI_high"])
    json.dump(meta, open(os.path.join(DATA, "revision_meta_heterogeneity.json"), "w"), indent=1)
    print("===== EXTERNAL VALIDATION RIGOR =====")
    print(df.to_string(index=False))
    print("\n===== META (DL random effects, median split) =====")
    print(json.dumps(meta, indent=1))
    return df, meta


def ph_test():
    from lifelines import CoxPHFitter
    from lifelines.statistics import proportional_hazard_test
    d = pd.read_csv(os.path.join(DATA, "TCGA_LIHC_clinical_with_PRS.csv"), index_col=0)
    d = d.dropna(subset=["os_years", "os_event"])
    d["gender_f"] = (d["gender"].astype(str).str.lower() == "female").astype(float)
    d["child_pugh_ord"] = d["child_pugh"].astype(str).str.lower().map({"a": 1.0, "b": 2.0, "c": 3.0})
    d["PRS_sd"] = (d["PRS"] - d["PRS"].mean()) / d["PRS"].std()
    d["age_sd"] = (d["age"] - d["age"].mean()) / d["age"].std()
    cols = ["os_years", "os_event", "PRS_sd", "age_sd", "gender_f", "child_pugh_ord"]
    dd = d[cols].dropna().astype(float).reset_index(drop=True)
    cph = CoxPHFitter(penalizer=1e-6).fit(dd, "os_years", "os_event")
    rows = []
    for transform, label in [("rank", "rank"), ("km", "km")]:
        try:
            res = proportional_hazard_test(cph, dd, time_transform=transform)
            s = res.summary
            for term in s.index:
                rows.append({"time_transform": label, "term": term,
                             "test_statistic": float(s.loc[term, "test_statistic"]),
                             "p": float(s.loc[term, "p"])})
        except Exception as e:
            rows.append({"time_transform": label, "term": f"ERR {e}", "test_statistic": np.nan, "p": np.nan})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(DATA, "revision_ph_test.csv"), index=False)
    print("\n===== PROPORTIONAL HAZARDS (Schoenfeld) =====")
    print(out.to_string(index=False))
    print("n complete-case:", len(dd), " events:", int(dd.os_event.sum()))
    return out


# --------------------------------------------------------------- purity proxy
IMMUNE_MARKERS = ["PTPRC", "CD3D", "CD3E", "CD2", "CD8A", "CD8B", "CD4",
                  "CD19", "MS4A1", "CD79A", "CD79B", "IGHM",
                  "LYZ", "CD14", "FCGR3A", "ITGAX", "CD1C",
                  "NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1"]
STROMAL_MARKERS = ["COL1A1", "COL1A2", "COL3A1", "COL5A2", "DCN", "LUM", "FBLN1",
                   "PDGFRB", "ACTA2", "TAGLN", "PECAM1", "VWF", "KDR", "ENG", "CDH5"]
CELL_SCORE_SETS = {
    "CD8_T": ["CD8A", "CD8B", "GZMK", "PRF1"],
    "NK": ["KLRD1", "KLRF1", "NKG7", "GNLY", "NCAM1"],
    "B": ["CD19", "MS4A1", "CD79A", "CD79B", "IGHM"],
    "Neutrophil": ["FCGR3B", "CSF3R", "S100A8", "S100A9"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "ENG", "CDH5"],
    "SingleGene_HAVCR2": ["HAVCR2"],
    "SingleGene_SIGLEC15": ["SIGLEC15"],
}


def partial_spearman(x, y, z):
    """Spearman partial correlation of x,y controlling z."""
    def rk(v):
        return stats.rankdata(v)
    rx, ry, rz = rk(x), rk(y), rk(z)
    import numpy.linalg as la
    def resid(a, b):
        A = np.column_stack([np.ones(len(b)), b])
        beta, *_ = la.lstsq(A, a, rcond=None)
        return a - A @ beta
    ex, ey = resid(rx, rz), resid(ry, rz)
    r = np.corrcoef(ex, ey)[0, 1]
    n = len(x)
    if n > 3 and abs(r) < 1:
        tstat = r * np.sqrt((n - 3) / (1 - r ** 2))
        p = 2 * stats.t.sf(abs(tstat), n - 3)
    else:
        p = np.nan
    return float(r), float(p)


def purity_confounding():
    cl = pd.read_csv(os.path.join(DATA, "TCGA_LIHC_clinical_with_PRS.csv"), index_col=0)
    ex = pd.read_csv(os.path.join(DATA, "TCGA_LIHC_rsem_normalized_matrix.csv"), index_col=0)
    # sample columns are full barcodes; clinical 'sample' is short
    # clinical 'sample' is lowercase 12-char (tcga-2v-a95s); expr cols are full 28-char barcodes
    cl["short"] = cl["sample"].astype(str).str.lower().str[:12]
    ex.columns = [c[:12].lower() for c in ex.columns]
    ex = ex.loc[:, ~ex.columns.duplicated()]
    cl = cl.loc[~cl["short"].duplicated()]
    common = [c for c in cl["short"] if c in ex.columns]
    ex = ex[common]
    cl = cl.set_index("short").loc[common]
    print(f"[purity] matched samples: {len(common)}")

    def zmean(genes):
        g = [x for x in genes if x in ex.index]
        z = ex.loc[g].apply(lambda v: (v - v.mean()) / (v.std() + 1e-12), axis=1)
        return z.mean(axis=0), len(g)

    imm, n_imm = zmean(IMMUNE_MARKERS)
    str_, n_str = zmean(STROMAL_MARKERS)
    noncancer = (imm.values + str_.values) / 2.0        # higher = less pure tumour
    prs = cl["PRS"].values

    rho, p = stats.spearmanr(prs, noncancer)
    rows = [{"feature": "PURITY_PROXY(non-tumour content)", "n_genes_immune": n_imm,
             "n_genes_stromal": n_str, "spearman_rho_with_PRS": round(float(rho), 3),
             "p_unadjusted": float(p), "spearman_partial_rho_given_purity": np.nan,
             "p_adjusted_for_purity": np.nan}]
    print(f"[purity] immune markers used={n_imm}, stromal={n_str}")
    for name, genes in CELL_SCORE_SETS.items():
        sc, ng = zmean(genes)
        sc = sc.values
        r0, p0 = stats.spearmanr(prs, sc)
        r1, p1 = partial_spearman(prs, sc, noncancer)
        rows.append({"feature": name, "n_genes_immune": ng, "n_genes_stromal": np.nan,
                     "spearman_rho_with_PRS": round(float(r0), 3), "p_unadjusted": float(p0),
                     "spearman_partial_rho_given_purity": round(float(r1), 3),
                     "p_adjusted_for_purity": float(p1)})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(DATA, "revision_purity_confounding.csv"), index=False)
    print("\n===== PURITY / COMPOSITION CONFOUNDING =====")
    print(out.to_string(index=False))
    return out


def n_table():
    rows = [
        {"analysis": "TCGA-LIHC primary tumours", "n": 371, "excluded": 0,
         "note": "all tumor samples after expression/clinical merge"},
        {"analysis": "TCGA-LIHC with OS follow-up", "n": 370, "excluded": 1,
         "note": "one tumour lacks os_time/os_event"},
        {"analysis": "Fig 8A KM (median split)", "n": 370, "excluded": 1,
         "note": "185 low-risk + 185 high-risk"},
        {"analysis": "Multivariable Cox complete-case", "n": 238, "excluded": 133,
         "note": "56 events; 132 missing Child-Pugh + 1 missing OS"},
        {"analysis": "Immune infiltration / checkpoint", "n": 371, "excluded": 0,
         "note": "median split 186 low + 185 high"},
        {"analysis": "Mutation frequency", "n": 371, "excluded": 0,
         "note": "186 low + 185 high"},
        {"analysis": "PXS tumour vs non-tumour", "n": 423, "excluded": 0,
         "note": "371 tumour + 52 adjacent non-tumour"},
        {"analysis": "External GSE14520", "n": 221, "excluded": np.nan, "note": ""},
        {"analysis": "External GSE76427", "n": 115, "excluded": np.nan, "note": ""},
        {"analysis": "External GSE116174", "n": 64, "excluded": np.nan, "note": ""},
        {"analysis": "External GSE10186", "n": 80, "excluded": np.nan, "note": ""},
        {"analysis": "External ICGC-LIRI-JP", "n": 203, "excluded": np.nan, "note": ""},
    ]
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(DATA, "revision_n_table.csv"), index=False)
    print("\n===== N PER ANALYSIS =====")
    print(out.to_string(index=False))
    return out


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "ext"):
        external_rigor()
    if which in ("all", "ph"):
        ph_test()
    if which in ("all", "purity"):
        purity_confounding()
    if which in ("all", "n"):
        n_table()
    print("\nBULK REANALYSES 2 DONE")
