# -*- coding: utf-8 -*-
"""PLOS ONE major-revision reanalyses (PONE-D-26-43453).

Writes NEW files only (originals untouched):
  data/processed/revision_nestedCV_cindex.csv   - leakage-free CV, SuperPC fixed
  data/processed/revision_multivariable_cox.csv - PRS + age + stage
  data/processed/revision_fdr.csv               - BH-FDR for immune/checkpoint/mutation
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
sys.path.insert(0, HERE)
DATA = os.path.join(ROOT, "data", "processed")
import core  # noqa: E402

SEED = 42


# ---------------------------------------------------------------- SuperPC (fixed)
def superpc_fixed(X, time, event, ncomp=1, frac=0.7):
    """Supervised PCA with convergent Cox fits and guarded PCA."""
    from lifelines import CoxPHFitter
    d = pd.DataFrame(X).copy()
    d["T"] = np.asarray(time, float)
    d["E"] = np.asarray(event, float).astype(int)
    d = d.dropna(subset=["T", "E"])
    genes = list(X.columns)
    scores = []
    for g in genes:
        try:
            cph = CoxPHFitter(penalizer=1e-4)
            cph.fit(d[[g, "T", "E"]], duration_col="T", event_col="E")
            scores.append(abs(float(cph.params_[g])))
        except Exception:
            scores.append(0.0)
    scores = np.asarray(scores, float)
    thr = np.percentile(scores, frac * 100) if np.max(scores) > 0 else 0.0
    sel = [genes[i] for i in range(len(genes)) if scores[i] >= max(thr, 1e-9)]
    if len(sel) < 3:
        v = X.var().sort_values(ascending=False)
        sel = v.head(min(15, X.shape[1])).index.tolist()
    from sklearn.decomposition import PCA
    Z = X[sel].fillna(X[sel].median())
    Z = (Z - Z.mean()) / (Z.std() + 1e-9)
    ncomp = max(1, min(ncomp, Z.shape[1], Z.shape[0] - 1))
    pc = PCA(n_components=ncomp).fit_transform(Z.values)
    # Cox on PC(s)
    dd = pd.DataFrame(pc, columns=[f"PC{i+1}" for i in range(ncomp)])
    dd["T"] = np.asarray(time, float)[: len(dd)]
    dd["E"] = np.asarray(event, float).astype(int)[: len(dd)]
    dd = dd.dropna(subset=["T", "E"])
    cph = CoxPHFitter(penalizer=1e-4).fit(dd, "T", "E")
    return cph.params_


# ---------------------------------------------------------------- data
def load_tcga():
    from analysis_tcga import load  # reuse the original loader
    return load()


def nested_cv():
    """10-fold CV with PRG screening REPEATED inside each training fold."""
    E, cl, _ = load_tcga()
    prg = json.load(open(os.path.join(DATA, "peroxisome_gene_set.json")))["genes"]
    genes = [g for g in prg if g in E.index]
    print(f"PRG overlap: {len(genes)}", flush=True)

    time = cl["os_years"].values.astype(float)
    event = cl["os_event"].values.astype(int)
    n = len(cl)
    rng = np.random.default_rng(SEED)
    fold = np.zeros(n, int)
    idx = rng.permutation(n)
    for f in range(10):
        fold[idx[f::10]] = f

    zoo = dict(core.MODEL_ZOO)
    zoo["SuperPC"] = lambda X, t, e: superpc_fixed(X, t, e)

    rows = []
    for name, fitter in zoo.items():
        cs = []
        for f in range(10):
            tr, va = fold != f, fold == f
            # --- screen INSIDE the training fold only (fixes leakage)
            sub = pd.DataFrame({g: E.loc[g, :].reindex(cl.index).values for g in genes},
                               index=cl.index).iloc[tr]
            sub["T"], sub["E"] = time[tr], event[tr]
            ps = []
            for g in genes:
                try:
                    cph = CoxPHFitter(penalizer=1e-4).fit(
                        sub[[g, "T", "E"]].dropna(), "T", "E")
                    ps.append((g, cph.summary["p"][g]))
                except Exception:
                    ps.append((g, 1.0))
            ps = pd.DataFrame(ps, columns=["gene", "p"]).sort_values("p")
            cand = ps[ps["p"] < 0.10]["gene"].tolist()
            if len(cand) < 8:
                cand = ps.head(20)["gene"].tolist()
            Xtr = E.loc[cand].T.reindex(cl.index[tr]).astype(float)
            Xtr = Xtr.fillna(Xtr.median())
            Xtr = (Xtr - Xtr.mean()) / (Xtr.std() + 1e-9)
            Xva = E.loc[cand].T.reindex(cl.index[va]).astype(float)
            Xva = (Xva - Xtr.mean()) / (Xtr.std() + 1e-9)  # train-fold scaling
            Xva = Xva.fillna(0)
            try:
                coef = fitter(Xtr, time[tr], event[tr])
                if coef is None or len(coef) == 0 or coef.isna().all():
                    continue
                sc = Xva.dot(coef.fillna(0))
                c, _ = core.cindex(time[va], event[va], sc.values)
                if not np.isnan(c):
                    cs.append(c)
            except Exception as e:
                print(f"  {name} fold{f}: {str(e)[:60]}", flush=True)
        rows.append({"model": name, "n_folds": len(cs),
                     "nested_cv_cindex": float(np.mean(cs)) if cs else np.nan,
                     "sd": float(np.std(cs)) if cs else np.nan})
        print(rows[-1], flush=True)
    out = pd.DataFrame(rows).sort_values("nested_cv_cindex", ascending=False)
    out.to_csv(os.path.join(DATA, "revision_nestedCV_cindex.csv"), index=False)
    print("saved revision_nestedCV_cindex.csv", flush=True)


# ---------------------------------------------------------------- multivariable Cox
def multivariable():
    """PRS + age + gender + Child-Pugh, complete case, vs clinical-only baseline.
    NOTE: tissue column in clinical_with_PRS.csv is corrupted (all 'normal');
    all 371 rows are tumor samples (verified previously)."""
    from scipy.stats import chi2 as _chi2
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    d = pd.read_csv(os.path.join(DATA, "TCGA_LIHC_clinical_with_PRS.csv"), index_col=0)
    d = d.dropna(subset=["os_years", "os_event"])
    d["gender_f"] = (d["gender"].astype(str).str.lower() == "female").astype(float)
    d["child_pugh_ord"] = d["child_pugh"].astype(str).str.lower().map(
        {"a": 1.0, "b": 2.0, "c": 3.0})
    d["PRS_sd"] = (d["PRS"] - d["PRS"].mean()) / d["PRS"].std()
    d["age_sd"] = (d["age"] - d["age"].mean()) / d["age"].std()
    cols_full = ["os_years", "os_event", "PRS_sd", "age_sd", "gender_f", "child_pugh_ord"]
    clin_cols = ["age_sd", "gender_f", "child_pugh_ord"]
    d_full = d[cols_full].dropna().astype(float)

    def fit(dd, cols):
        return CoxPHFitter(penalizer=1e-6).fit(dd[cols], "os_years", "os_event")

    cph_full = fit(d_full, cols_full)
    summ = cph_full.summary[["exp(coef)", "exp(coef) lower 95%",
                             "exp(coef) upper 95%", "p"]]
    rows = [{"term": t, "HR": float(summ.loc[t, "exp(coef)"]),
             "CI_low": float(summ.loc[t, "exp(coef) lower 95%"]),
             "CI_high": float(summ.loc[t, "exp(coef) upper 95%"]),
             "p": float(summ.loc[t, "p"])} for t in summ.index]
    d_cli = d_full[["os_years", "os_event"] + clin_cols]
    cph_cli = fit(d_cli, ["os_years", "os_event"] + clin_cols)
    r_full = -cph_full.predict_partial_hazard(d_full[cols_full[2:]])
    r_cli = -cph_cli.predict_partial_hazard(d_cli[clin_cols])
    c_full = concordance_index(d_full["os_years"], r_full, d_full["os_event"])
    c_cli = concordance_index(d_cli["os_years"], r_cli, d_cli["os_event"])
    lr_stat = 2 * (cph_full.log_likelihood_ - cph_cli.log_likelihood_)
    lr_p = float(_chi2.sf(lr_stat, df=1))
    meta = {"n_complete_case": int(len(d_full)),
            "n_events": int(d_full["os_event"].sum()),
            "n_child_pugh_missing": int(d["child_pugh_ord"].isna().sum()),
            "c_index_clinical": round(float(c_cli), 4),
            "c_index_clinical_plus_PRS": round(float(c_full), 4),
            "LRT_add_PRS_stat": round(float(lr_stat), 3),
            "LRT_add_PRS_p": lr_p}
    pd.DataFrame(rows).to_csv(os.path.join(DATA, "revision_multivariable_cox.csv"),
                              index=False)
    json.dump(meta, open(os.path.join(DATA, "revision_multivariable_summary.json"), "w"),
              indent=1)
    print("MULTIVARIABLE:", json.dumps(meta), flush=True)
    print(pd.DataFrame(rows).to_string(), flush=True)


# ---------------------------------------------------------------- FDR
def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    n = len(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def fdr():
    rows = []
    imm = pd.read_csv(os.path.join(DATA, "immune_infiltration_PRS.csv"))
    imm = imm[imm["group"] == "p"].reset_index(drop=True)
    q = bh(imm["mean"].values)
    for s, p, qq in zip(imm["set"], imm["mean"], q):
        rows.append({"family": "immune_ssGSEA", "feature": s, "p": float(p), "BH_FDR": float(qq)})
    cp = pd.read_csv(os.path.join(DATA, "checkpoints_PRS.csv"))
    q = bh(cp["p"].values)
    for g, p, qq in zip(cp["gene"], cp["p"], q):
        rows.append({"family": "checkpoint", "feature": g, "p": float(p), "BH_FDR": float(qq)})
    mut = pd.read_csv(os.path.join(DATA, "mutation_PRS.csv"))
    genes = sorted(mut["gene"].dropna().unique())
    ps = []
    for g in genes:
        t = mut[mut["gene"] == g].set_index("group")
        try:
            tab = [[t.loc["low", "n_mut"], t.loc["low", "n_total"] - t.loc["low", "n_mut"]],
                   [t.loc["high", "n_mut"], t.loc["high", "n_total"] - t.loc["high", "n_mut"]]]
            ps.append(stats.fisher_exact(tab)[1])
        except Exception:
            ps.append(np.nan)
    ok = ~np.isnan(ps)
    q = bh(np.asarray(ps)[ok])
    for g, p, qq in zip(np.asarray(genes)[ok], np.asarray(ps)[ok], q):
        rows.append({"family": "mutation_fisher", "feature": g, "p": float(p), "BH_FDR": float(qq)})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(DATA, "revision_fdr.csv"), index=False)
    print(out[out["p"] < 0.05].to_string(), flush=True)
    print("saved revision_fdr.csv", flush=True)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "cv"):
        nested_cv()
    if which in ("all", "cox"):
        multivariable()
    if which in ("all", "fdr"):
        fdr()
    print("REANALYSES DONE")
