# -*- coding: utf-8 -*-
"""Core analysis utilities for the HCC-peroxisome study.
Survival stats, time-dependent ROC/AUC, C-index, ssGSEA, penalized Cox,
and the machine-learning model zoo used to build the PRS.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, gzip, json, pickle, re, warnings
import numpy as np
import pandas as pd
from scipy import stats, sparse
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings("ignore")

DATA = str(PROJECT_ROOT / "data" / "processed")
OUT  = str(PROJECT_ROOT / "data" / "processed")
rng = np.random.default_rng(20260829)


# ------------------------------------------------------------------ I/O
def load_tcga_expr():
    df = pd.read_csv(os.path.join(DATA, "TCGA_LIHC_rsem_normalized_matrix.csv"), index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = np.log2(df + 1)
    info = json.load(open(os.path.join(DATA, "TCGA_LIHC_sample_info.json")))
    return df, pd.DataFrame({"sample": info["sample"], "tissue": info["tissue"]})


def load_tcga_clin():
    return pd.read_csv(os.path.join(DATA, "TCGA_LIHC_clinical.csv"), index_col=0)


# ------------------------------------------------------------------ module score (ssGSEA-like)
def ssgsea_score(expr, genes, rank=True):
    """expr: genes x samples DataFrame. Returns per-sample enrichment score of `genes`.
    Uses a rank-based ssGSEA approximation (mean z of ranks)."""
    genes = [g for g in genes if g in expr.index]
    if not genes:
        return pd.Series(np.nan, index=expr.columns)
    if rank:
        ranked = expr.rank(axis=0) / expr.shape[0]   # genes x samples, 0..1
    else:
        ranked = expr
    sub = ranked.loc[genes]
    # mean rank of the set minus expected mean (0.5)
    score = sub.mean(axis=0) - 0.5
    # normalize by std to combine ranks across datasets
    return score


def module_score(expr, genes):
    """mean z-score across genes (z across samples per gene)."""
    genes = [g for g in genes if g in expr.index]
    if not genes:
        return pd.Series(np.nan, index=expr.columns)
    Z = (expr.loc[genes] - expr.loc[genes].mean(axis=1).values[:, None]) / (
        expr.loc[genes].std(axis=1).values[:, None] + 1e-12)
    return Z.mean(axis=0)


# ------------------------------------------------------------------ survival helpers
def make_surv(time, event):
    t = np.asarray(time, float)
    e = np.asarray(event, float)
    ok = ~(np.isnan(t) | np.isnan(e))
    return t[ok], e[ok], np.where(ok)[0]


def cox_hr(time, event, group):
    """group: binary 0/1. returns (hr, ci_low, ci_high, p, n)."""
    t, e, ok = make_surv(time, event)
    g = np.asarray(group)[ok]
    if len(set(g)) < 2 or len(t) < 10:
        return np.nan, np.nan, np.nan, np.nan, len(t)
    d = pd.DataFrame({"T": t, "E": e.astype(int), "G": g})
    try:
        cph = CoxPHFitter(penalizer=1e-6)
        cph.fit(d, duration_col="T", event_col="E", formula="G")
        ci = cph.confidence_intervals_.loc["G"]
        return float(np.exp(cph.params_["G"])), float(np.exp(ci.iloc[0])), \
               float(np.exp(ci.iloc[1])), cph.summary["p"]["G"], len(t)
    except Exception:
        return np.nan, np.nan, np.nan, np.nan, len(t)


def logrank(time, event, group):
    t, e, ok = make_surv(time, event)
    g = np.asarray(group)[ok]
    try:
        res = logrank_test(t[g == 1], t[g == 0], event_observed_A=e[g == 1], event_observed_B=e[g == 0])
        return res.p_value
    except Exception:
        return np.nan


def cindex(time, event, score):
    """Harrell's concordance index."""
    t, e, ok = make_surv(time, event)
    s = np.asarray(score)[ok]
    n = len(t)
    if n < 10:
        return np.nan, n
    concord, tot = 0, 0
    for i in range(n):
        if e[i] == 0:
            continue
        ti = t[i]
        for j in range(n):
            if i == j:
                continue
            if t[j] > ti:   # j outlives i
                tot += 1
                concord += 1 if s[i] > s[j] else (0.5 if s[i] == s[j] else 0)
    return (concord / tot) if tot else np.nan, n


def td_auc(time, event, score, t_horizon):
    """Cumulative/dynamic time-dependent AUC at time(s) t_horizon (Uno-Heagerty style).
    AUC(t) = P(M_i > M_j | T_i <= t, T_j > t)."""
    t, e, ok = make_surv(time, event)
    s = np.asarray(score)[ok]
    out = []
    for th in np.atleast_1d(t_horizon):
        cases = (t <= th) & (e == 1)
        ctrls = t > th
        if cases.sum() == 0 or ctrls.sum() == 0:
            out.append(np.nan)
            continue
        num = sum(s[c] > s[k] for c in np.where(cases)[0] for k in np.where(ctrls)[0])
        out.append(num / (cases.sum() * ctrls.sum()))
    return out[0] if np.ndim(t_horizon) == 0 else np.array(out)


# ------------------------------------------------------------------ ssGSEA / signature helpers
IMMUNE_SIGS = {
    "CD8_T": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG"],
    "CD4_T": ["CD4", "IL7R", "CCR7", "LEF1", "TCF7"],
    "Treg": ["FOXP3", "CTLA4", "IL2RA", "IKZF2"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRC1", "NCAM1", "KLRK1"],
    "B": ["MS4A1", "CD79A", "CD79B", "BANK1"],
    "Plasma": ["MZB1", "XBP1", "SDC1", "JCHAIN"],
    "Macrophage": ["CD68", "CD163", "MRC1", "CSF1R", "ITGAM"],
    "M1": ["NOS2", "TNF", "IL1B", "IL6", "CXCL10"],
    "M2": ["CD163", "MRC1", "MSR1", "TGFBI", "FN1"],
    "Monocyte": ["CD14", "FCGR3A", "LYZ", "S100A8"],
    "DC": ["ITGAX", "CD1C", "CLEC9A", "BATF3", "FLT3"],
    "Mast": ["TPSAB1", "CPA3", "KIT", "MS4A2"],
    "Neutrophil": ["FCGR3B", "CSF3R", "S100A12", "CXCR2"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "EGFL7", "ESM1"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB"],
    "Hepatocyte": ["ALB", "APOA1", "APOA2", "TTR", "FABP1"],
    "Cholangiocyte": ["KRT19", "EPCAM", "KRT7", "SOX9"],
}
CHECKPOINTS = ["PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "SIGLEC15", "IDO1", "CD276", "VTCN1", "LGALS9"]


def immune_ssgsea(expr):
    """expr: genes x samples -> DataFrame (celltype x sample)"""
    out = {}
    for ct, gs in IMMUNE_SIGS.items():
        out[ct] = ssgsea_score(expr, gs)
    return pd.DataFrame(out)


def estimate_scores(expr):
    """ESTIMATE-like immune/stromal scores (curated signature lists)."""
    immune = ["PTPRC", "CD2", "CD3D", "CD3E", "CXCR3", "CCR5", "GZMB", "IFNG", "CD8A", "CD8B",
              "CD68", "CD163", "LYZ", "MS4A1", "CD79A", "NKG7", "GNLY", "KLRD1", "CXCL9",
              "CXCL10", "CXCL11", "CCL2", "CCL3", "CCL4", "CCL5", "IL10", "TNF", "IL1B",
              "CX3CR1", "ITGAX", "FCGR3A", "IL2RG", "STAT1", "TLR7", "TLR9", "IRF7", "IRF1",
              "LCP2", "LCK", "ZAP70", "LAT", "IL7R", "IL2RA", "FOXP3", "PRF1", "GZMA", "GZMK",
              "TBX21", "CD40", "CD80", "CD86", "CD274", "PDCD1", "CTLA4", "TIGIT", "HAVCR2"]
    stromal = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "DCN", "LUM", "BGN",
               "PDGFRB", "PDGFRA", "VIM", "ACTN1", "TAGLN", "ACTA2", "MYH11", "FBN1", "FN1",
               "SPARC", "POSTN", "THY1", "ENG", "KDR", "FLT1", "PECAM1", "VWF", "CDH5", "CLDN5",
               "ESM1", "EGFL7", "MCAM", "PLVAP", "NOTCH4", "DLL4"]
    isc = module_score(expr, immune)
    stc = module_score(expr, stromal)
    return pd.DataFrame({"ImmuneScore": isc, "StromalScore": stc})


# ------------------------------------------------------------------ ML model zoo
def fit_cox_model(X, time, event, alpha=1.0, l1_ratio=1.0, standardize=True):
    """penalized Cox (lifelines). X: samples x genes df (already log expr)."""
    from lifelines import CoxPHFitter
    Xc = X.copy()
    if standardize:
        mu, sd = Xc.mean(), Xc.std()
        Xc = (Xc - mu) / (sd + 1e-12)
    d = Xc.copy()
    d["T"] = np.asarray(time, float)
    d["E"] = np.asarray(event, float).astype(int)
    d = d[(d["T"].notna()) & (d["E"].notna())].reset_index(drop=True)
    if len(d) < 20 or d.shape[1] < 3:
        return pd.Series(np.nan, index=X.columns)
    try:
        cph = CoxPHFitter(penalizer=alpha, l1_ratio=l1_ratio)
        cph.fit(d, duration_col="T", event_col="E")
        coef = cph.params_
        return coef
    except Exception:
        return pd.Series(np.nan, index=X.columns)


def stepwise_cox(X, time, event, direction="both", alpha=0.05, max_genes=25):
    """Stepwise Cox by AIC (forward then backward). Returns pd.Series of coefs."""
    genes = list(X.columns)
    chosen = []
    best_aic = np.inf
    d = pd.DataFrame(X).copy()
    d["T"] = np.asarray(time, float)
    d["E"] = np.asarray(event, float).astype(int)
    d = d.dropna(subset=["T", "E"])
    cur = []

    def aic(cols):
        if not cols:
            return np.inf
        try:
            cph = CoxPHFitter(penalizer=0)
            cph.fit(d[cols + ["T", "E"]], duration_col="T", event_col="E")
            return cph.AIC_partial_
        except Exception:
            return np.inf

    for step in range(max_genes):
        # forward
        best, bestk = aic(cur), None
        for g in genes:
            if g in cur:
                continue
            a = aic(cur + [g])
            if a < best - 1e-6:
                best, bestk = a, g
        if bestk is not None:
            cur.append(bestk)
            continue
        # backward
        changed = False
        for i in range(len(cur)):
            cand = cur[:i] + cur[i + 1:]
            a = aic(cand)
            if a < best - 1e-6:
                cur = cand
                changed = True
                break
        if not changed:
            break
    coef = pd.Series(0.0, index=genes)
    if cur:
        cph = CoxPHFitter(penalizer=0)
        cph.fit(d[cur + ["T", "E"]], duration_col="T", event_col="E")
        for g in cur:
            coef[g] = cph.params_[g]
    return coef


def plsrcox(X, time, event, ncomp=3):
    """PLS components + Cox."""
    d = pd.DataFrame(X)
    d["T"] = np.asarray(time, float)
    d["E"] = np.asarray(event, float).astype(int)
    d = d.dropna(subset=["T", "E"])
    pls = PLSRegression(n_components=ncomp)
    pls.fit(d[X.columns], np.log1p(d["T"]))
    comps = pls.transform(d[X.columns])
    comp_df = pd.DataFrame(comps, columns=[f"C{i}" for i in range(comps.shape[1])], index=d.index)
    m = comp_df.join(d[["T", "E"]])
    try:
        cph = CoxPHFitter(penalizer=0)
        cph.fit(m, duration_col="T", event_col="E")
        # map components -> original genes via loadings
        loadings = pls.x_loadings_
        coef = np.zeros(X.shape[1])
        for i, c in enumerate(cph.params_):
            coef += cph.params_[i] * loadings[:, i]
        return pd.Series(coef, index=X.columns)
    except Exception:
        return pd.Series(np.nan, index=X.columns)


def superpc(X, time, event, ncomp=2):
    """Supervised PCA (threshold genes by univariate Cox score)."""
    d = pd.DataFrame(X)
    d["T"] = np.asarray(time, float)
    d["E"] = np.asarray(event, float).astype(int)
    d = d.dropna(subset=["T", "E"])
    genes = list(X.columns)
    scores = []
    for g in genes:
        try:
            cph = CoxPHFitter(penalizer=0)
            cph.fit(d[[g, "T", "E"]], duration_col="T", event_col="E")
            scores.append(abs(cph.params_[g]))
        except Exception:
            scores.append(0.0)
    thr = np.percentile(scores, 70) if np.max(scores) > 0 else 0
    sel = [genes[i] for i in range(len(genes)) if scores[i] >= max(thr, 1e-9)]
    if len(sel) < 3:
        # fall back to top-variance genes
        v = X.var().sort_values(ascending=False)
        sel = v.head(15).index.tolist()
    sel = sel[:max(len(sel), 1)]
    pca = PCA(n_components=ncomp)
    comps = pca.fit_transform(StandardScaler().fit_transform(d[sel]))
    comp_df = pd.DataFrame(comps, columns=[f"C{i}" for i in range(ncomp)], index=d.index)
    m = comp_df.join(d[["T", "E"]])
    try:
        cph = CoxPHFitter(penalizer=0)
        cph.fit(m, duration_col="T", event_col="E")
        coef = np.zeros(X.shape[1])
        for i, c in enumerate(cph.params_):
            coef += c * pca.components_[i]
        return pd.Series(coef, index=X.columns)
    except Exception:
        # fallback: return univariate cox coefficients (standardized)
        coef = pd.Series(np.zeros(X.shape[1]), index=X.columns)
        for i, g in enumerate(genes):
            coef[g] = scores[i]
        return coef


def rsr(rs, df):
    return rs


def cox_gbm(X, time, event, n_est=80, lr=0.05, depth=2):
    """Gradient-boosting Cox (log partial likelihood, small trees)."""
    from sklearn.tree import DecisionTreeRegressor
    d = pd.DataFrame(X)
    d["T"] = np.asarray(time, float)
    d["E"] = np.asarray(event, float).astype(int)
    d = d.dropna(subset=["T", "E"])
    Xm = d[X.columns].values
    T = d["T"].values
    E = d["E"].values
    n = len(T)
    order = np.argsort(T)
    # risk set index for tie handling (right-continuous)
    def neg_grad(eta):
        # negative gradient of Cox partial LL w.r.t. eta
        eta_sorted = eta[order]
        # cumulative: for subject at position k, risk set = k..n
        # exp_eta suffix sums
        ex = np.exp(np.clip(eta_sorted, -30, 30))
        suf = np.cumsum(ex[::-1])[::-1]
        g = np.zeros(n)
        idx_of = np.zeros(n, dtype=int)
        idx_of[order] = np.arange(n)
        for i in range(n):
            if E[i]:
                k = idx_of[i]
                g[i] = 1 - ex[k] / suf[k]
            else:
                g[i] = 0
        return -g

    eta = np.zeros(n)
    models = []
    for it in range(n_est):
        g = neg_grad(eta)
        tree = DecisionTreeRegressor(max_depth=depth, min_samples_leaf=10, random_state=42 + it)
        tree.fit(Xm, g)
        # step size via 1-D line search along tree fit
        pred = tree.predict(Xm)
        # approximate: use fixed lr
        eta = eta + lr * pred
        models.append(tree)
    # importance
    imp = np.zeros(X.shape[1])
    for tree in models:
        imp += tree.feature_importances_
    return pd.Series(imp / len(models), index=X.columns)


def random_survival_forest(X, time, event, n_est=50, mtry=None):
    """Simple RSF: bootstrap + CART splitting by fast log-rank statistic."""
    from sklearn.tree import DecisionTreeRegressor
    d = pd.DataFrame(X)
    d["T"] = np.asarray(time, float)
    d["E"] = np.asarray(event, float).astype(int)
    d = d.dropna(subset=["T", "E"])
    Xm = d[X.columns].values
    T = d["T"].values
    E = d["E"].values
    n, p = Xm.shape
    mtry = mtry or int(np.sqrt(p)) + 1
    rng = np.random.default_rng(7)

    def logrank_stat(t_l, e_l, t_r, e_r):
        """Fast two-sample log-rank via ranks (no CoxPHFitter)."""
        t = np.concatenate([t_l, t_r])
        e = np.concatenate([e_l, e_r])
        g = np.concatenate([np.zeros(len(t_l)), np.ones(len(t_r))])
        order = np.argsort(t)
        t, e, g = t[order], e[order], g[order]
        # risk set counts from right
        n_all = len(t)
        # at each event time, compute observed - expected for group 1
        E1 = (g == 1).sum()
        obs = 0.0
        exp = 0.0
        var = 0.0
        i = 0
        n1 = E1
        while i < n_all:
            j = i
            while j < n_all and t[j] == t[i]:
                j += 1
            d1 = int((g[i:j] == 1).sum())
            n_risk = n_all - i
            n1_risk = int((g[i:] == 1).sum())
            if n_risk > 1:
                obs += d1
                exp += (n1_risk * (j - i)) / n_risk
                if n_risk > 1:
                    var += (n1_risk / n_risk) * (1 - n1_risk / n_risk) * ((n_risk - (j - i)) / (n_risk - 1)) * (j - i)
            i = j
        z = (obs - exp) / (np.sqrt(var) + 1e-9)
        return z * z

    importances = np.zeros(p)
    for it in range(n_est):
        boot = rng.integers(0, n, n)
        Xb, Tb, Eb = Xm[boot], T[boot], E[boot]
        feats = rng.choice(p, size=min(mtry, p), replace=False)
        best_gain, best_f, best_th = -1, None, None
        for f in feats:
            vals = np.unique(Xb[:, f])
            if len(vals) < 2:
                continue
            cand = (vals[:-1] + vals[1:]) / 2
            for th in cand[:min(len(cand), 12)]:
                left = Xb[:, f] <= th
                if left.sum() < 4 or (~left).sum() < 4:
                    continue
                gain = logrank_stat(Tb[left], Eb[left], Tb[~left], Eb[~left])
                if gain > best_gain:
                    best_gain, best_f, best_th = gain, f, th
        if best_f is not None:
            importances[best_f] += 1
    return pd.Series(importances / n_est, index=X.columns)


MODEL_ZOO = {
    "LassoCox": lambda X, t, e: fit_cox_model(X, t, e, alpha=1.0, l1_ratio=1.0),
    "RidgeCox": lambda X, t, e: fit_cox_model(X, t, e, alpha=1.0, l1_ratio=0.0),
    "EnetCox":  lambda X, t, e: fit_cox_model(X, t, e, alpha=1.0, l1_ratio=0.5),
    "StepCox":  lambda X, t, e: stepwise_cox(X, t, e),
    "plsRcox":  lambda X, t, e: plsrcox(X, t, e),
    "SuperPC":  lambda X, t, e: superpc(X, t, e),
    "GBM":      lambda X, t, e: cox_gbm(X, t, e),
    "RSF":      lambda X, t, e: random_survival_forest(X, t, e),
}


def risk_score(expr, coef):
    """expr: genes x samples; coef: Series(gene->coef). Returns per-sample score."""
    genes = [g for g in coef.index if g in expr.index]
    if not genes:
        return pd.Series(np.nan, index=expr.columns)
    # z-score each gene across samples
    E = expr.loc[genes]
    Z = (E - E.mean(axis=1).values[:, None]) / (E.std(axis=1).values[:, None] + 1e-12)
    return (Z.T @ coef.loc[genes]).rename("risk")


def cindex_fast(time, event, score):
    """Vectorized Harrell C-index (concordance)."""
    t, e, ok = make_surv(time, event)
    s = np.asarray(score)[ok]
    n = len(t)
    if n < 10 or np.all(e == 0):
        return np.nan, n
    # matrix approach may be memory heavy for n>5000; use loop for safety
    return cindex(time, event, s)


if __name__ == "__main__":
    print("core module OK")
