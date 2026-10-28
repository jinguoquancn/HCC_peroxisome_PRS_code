# -*- coding: utf-8 -*-
"""TCGA-LIHC main analysis:
1) peroxisome module score (PXS), tumor vs normal
2) DEG tumor vs normal & PXS-high vs PXS-low
3) univariate Cox on PRGs
4) WGCNA-like module analysis
5) machine-learning model zoo (8 algorithms) with 10-fold CV
6) final PRS coefficients
Saves all intermediate tables to data/processed/."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, json, sys, warnings
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr, pearsonr
warnings.filterwarnings("ignore")

sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
import core
from core import module_score, ssgsea_score, cox_hr, logrank, cindex, td_auc

DATA = str(PROJECT_ROOT / "data" / "processed")
rng = np.random.default_rng(20260829)

# ------------------------------------------------------------------ load
def load():
    expr, info = core.load_tcga_expr()
    clin = core.load_tcga_clin()
    # map clinical by patient; expression by sample
    pat = {s: "-".join(s.split("-")[:3]).lower() for s in info["sample"]}
    info["patient"] = [pat.get(s, s) for s in info["sample"]]
    clin_p = clin[clin["os_event"].notna()]
    # keep tumor samples only
    tumor = info[info["tissue"] == "tumor"]
    pat = [pat.get(s, s) for s in tumor["sample"]]
    E = expr[tumor["sample"]]
    E.columns = pat
    E = E.loc[:, ~E.columns.duplicated(keep="first")]
    cl = clin_p.reindex(pat)
    cl.index = pat
    return E, cl, tumor


def main():
    E, cl, tumor = load()
    print("TCGA tumor expr:", E.shape, "clinical:", cl.shape, flush=True)
    prg = json.load(open(os.path.join(DATA, "peroxisome_gene_set.json")))
    prg_genes = prg["genes"]

    # ---- PXS (peroxisome module score)
    ov = [g for g in prg_genes if g in E.index]
    print(f"PRG genes overlapping TCGA expr: {len(ov)}/{len(prg_genes)}", flush=True)
    pxs = module_score(E, ov)
    pxs_core = module_score(E, [g for g in prg.get("modules", {}).get("pex5_import_dependent_core", []) if g in E.index])
    cl["PXS"] = pxs.reindex(cl.index).values
    cl["PXS_core"] = pxs_core.reindex(cl.index).values
    cl["PXS_high"] = (cl["PXS"] > np.median(cl["PXS"])).astype(int)
    cl.to_csv(os.path.join(DATA, "TCGA_LIHC_clinical_with_PXS.csv"))
    print("PXS summary:", cl["PXS"].describe().round(3).to_dict(), flush=True)

    # ---- tumor vs normal PXS
    expr_all, info_all = core.load_tcga_expr()
    norm = info_all[info_all["tissue"] == "normal"]
    cols = list(tumor["sample"]) + list(norm["sample"])
    cols = [c for c in cols if c in expr_all.columns]
    Eall = expr_all[cols]
    Eall.columns = [("T" if c in set(tumor["sample"]) else "N") for c in Eall.columns]
    pxs_all = module_score(Eall, prg_genes)
    pxs_df = pd.DataFrame({"group": Eall.columns, "PXS": pxs_all.values})
    pxs_df.to_csv(os.path.join(DATA, "TCGA_PXS_tumor_vs_normal.csv"), index=False)
    print("tumor/normal PXS:", flush=True)
    print(pxs_df.groupby("group")["PXS"].mean(), flush=True)

    # ---- DEG: tumor vs normal (edgeR-like via Wilcoxon + fold change)
    t_expr = expr_all[tumor["sample"]]
    n_expr = expr_all[norm["sample"]]
    de = []
    for g in t_expr.index:
        tv, nv = t_expr.loc[g].dropna().values, n_expr.loc[g].dropna().values
        if len(tv) < 10 or len(nv) < 10:
            continue
        fc = np.log2(np.mean(tv) + 1) - np.log2(np.mean(nv) + 1)
        p = stats.mannwhitneyu(tv, nv, alternative="two-sided").pvalue
        de.append((g, fc, p))
    de_df = pd.DataFrame(de, columns=["gene", "logFC", "p"]).dropna()
    de_df["FDR"] = stats.false_discovery_control(de_df["p"].values)
    de_df["sig"] = (de_df["FDR"] < 0.05) & (de_df["logFC"].abs() > 0.5)
    de_df = de_df.sort_values("p")
    de_df.to_csv(os.path.join(DATA, "TCGA_DEG_tumor_vs_normal.csv"), index=False)
    print("DEG tumor vs normal sig:", de_df["sig"].sum(), flush=True)

    # ---- univariate Cox on PRGs
    ucox = []
    surv = cl[["os_years", "os_event"]].copy()
    surv.columns = ["T", "E"]
    for g in prg_genes:
        if g not in E.index:
            continue
        x = E.loc[g].reindex(cl.index)
        if x.isna().sum() > 30:
            continue
        d = surv.copy()
        d["X"] = x.values
        d = d.dropna()
        try:
            from lifelines import CoxPHFitter
            cph = CoxPHFitter(penalizer=1e-6)
            cph.fit(d[["T", "E", "X"]], duration_col="T", event_col="E")
            ci = cph.confidence_intervals_.loc["X"]
            ucox.append((g, cph.params_["X"], ci.iloc[0], ci.iloc[1], cph.summary["p"]["X"]))
        except Exception as e:
            if len(ucox) < 1 and str(g) in ("ABCD1",):
                import traceback
                traceback.print_exc()
                print("  ucox err", g, repr(e)[:200], flush=True)
            pass
    ucox_df = pd.DataFrame(ucox, columns=["gene", "coef", "low", "high", "p"])
    ucox_df["HR"] = np.exp(ucox_df["coef"])
    ucox_df = ucox_df.dropna(subset=["p"])
    if len(ucox_df):
        ucox_df["FDR"] = stats.false_discovery_control(ucox_df["p"].values)
    ucox_df = ucox_df.sort_values("p")
    ucox_df.to_csv(os.path.join(DATA, "TCGA_univariate_Cox_PRG.csv"), index=False)
    print("univariate Cox PRGs n:", len(ucox_df), "sig<0.05:", (ucox_df.p < 0.05).sum(), flush=True)

    # ---- WGCNA-lite on all genes (soft threshold = hardcoded power 6, correlation-based)
    var = E.var(axis=1).sort_values(ascending=False)
    sel = var.head(2500).index.tolist()
    Es = E.loc[sel]
    # Spearman correlation via ranks (faster than pandas corr for 2500x2500)
    from scipy.stats import rankdata
    Rr = np.apply_along_axis(rankdata, 1, Es.values)   # rank each gene across samples
    Rr = (Rr - Rr.mean(axis=1, keepdims=True)) / (Rr.std(axis=1, keepdims=True) + 1e-12)
    R = np.corrcoef(Rr)
    # adjacency & TOM-lite: use signed adjacency ^6 then distance
    Adj = (np.abs(R) ** 6)
    # hierarchical clustering + dynamic cut approximation (fixed height split)
    from scipy.cluster.hierarchy import linkage, fcluster
    dist = 1 - Adj
    Z = linkage(dist[np.triu_indices(len(sel), 1)], method="average")
    labels = fcluster(Z, t=0.7, criterion="distance")
    modules = {}
    for lab in set(labels):
        gs = [sel[i] for i in range(len(sel)) if labels[i] == lab]
        if len(gs) >= 5:
            modules[lab] = gs
    mod_df = pd.DataFrame([{"module": k, "n_genes": len(v), "genes": ";".join(v)}
                           for k, v in modules.items()])
    # module eigengene (first PC) and correlation with PXS + survival
    rows = []
    for k, gs in modules.items():
        eg = E.loc[gs].T
        eg = (eg - eg.mean()) / (eg.std() + 1e-12)
        pc1 = eg.iloc[:, 0]
        r_pxs, _ = pearsonr(pc1, cl["PXS"].reindex(pc1.index).values)
        # survival corr
        d = pd.DataFrame({"T": cl["os_years"].values, "E": cl["os_event"].values,
                          "M": pc1.reindex(cl.index).values}).dropna()
        try:
            from lifelines import CoxPHFitter
            cph = CoxPHFitter(penalizer=1e-6)
            cph.fit(d, duration_col="T", event_col="E", formula="M")
            hr, p = float(np.exp(cph.params_["M"])), cph.summary["p"]["M"]
        except Exception:
            hr, p = np.nan, np.nan
        rows.append({"module": k, "n_genes": len(gs), "cor_PXS": r_pxs, "HR": hr, "p_surv": p})
    wgc = pd.DataFrame(rows).sort_values("p_surv")
    wgc.to_csv(os.path.join(DATA, "TCGA_WGCNA_modules.csv"), index=False)
    print("WGCNA modules:", len(modules), flush=True)

    # ---- ML model zoo with 10-fold CV on PRG candidates (uni-Cox p<0.10)
    cand = ucox_df[ucox_df["p"] < 0.10]["gene"].tolist()
    print("candidate genes for ML:", len(cand), flush=True)
    if len(cand) < 8:
        cand = ucox_df.head(20)["gene"].tolist()
    X = E.loc[cand].T.reindex(cl.index).astype(float)
    X = X.fillna(X.median())
    X = (X - X.mean()) / (X.std() + 1e-9)
    X = X.loc[:, X.std() > 0]
    print("ML input matrix:", X.shape, flush=True)
    time, event = cl["os_years"].values, cl["os_event"].values
    n = len(cl)
    fold = np.zeros(n, int)
    idx = rng.permutation(n)
    for f in range(10):
        fold[idx[f::10]] = f

    results = []
    for name, fitter in core.MODEL_ZOO.items():
        cv_c = []
        for f in range(10):
            tr, va = fold != f, fold == f
            try:
                coef = fitter(X.iloc[tr], time[tr], event[tr])
            except Exception as e:
                print(f"  {name} fold{f} err {str(e)[:60]}", flush=True)
                continue
            if coef is None or len(coef) == 0 or coef.isna().all():
                continue
            sc = X.iloc[va].dot(coef.fillna(0))
            c, _ = core.cindex(time[va], event[va], sc.values)
            if not np.isnan(c):
                cv_c.append(c)
        results.append({"model": name, "cv_cindex": np.mean(cv_c) if cv_c else np.nan,
                        "sd": np.std(cv_c) if cv_c else np.nan})
    cv_df = pd.DataFrame(results).sort_values("cv_cindex", ascending=False)
    cv_df.to_csv(os.path.join(DATA, "TCGA_ML_CV_cindex.csv"), index=False)
    print("CV C-index by model:\n", cv_df.to_string(), flush=True)

    # ---- final model: refit best CV model on full data
    best = cv_df.iloc[0]["model"]
    try:
        coef_all = core.MODEL_ZOO[best](X, time, event)
    except Exception as e:
        print("final fit err", e, flush=True)
        coef_all = core.MODEL_ZOO["LassoCox"](X, time, event)
        best = "LassoCox"
    coef_all = coef_all.dropna()
    coef_all = coef_all[coef_all != 0]
    prs_coef = coef_all.sort_values(ascending=False)
    prs_coef.to_csv(os.path.join(DATA, "TCGA_PRS_coefficients.csv"))
    print(f"FINAL PRS model: {best}, n genes = {len(prs_coef)}", flush=True)
    print(prs_coef.to_string(), flush=True)
    json.dump({"model": best}, open(os.path.join(DATA, "PRS_model_info.json"), "w"))

    # risk score in TCGA
    risk = core.risk_score(E, prs_coef)
    cl["PRS"] = risk.reindex(cl.index).values
    cl["PRS_high"] = (cl["PRS"] > np.median(cl["PRS"])).astype(int)
    cl.to_csv(os.path.join(DATA, "TCGA_LIHC_clinical_with_PRS.csv"))
    print("TCGA PRS C-index:", core.cindex(cl["os_years"], cl["os_event"], cl["PRS"]), flush=True)
    print("TCGA ANALYSIS DONE", flush=True)


if __name__ == "__main__":
    main()
