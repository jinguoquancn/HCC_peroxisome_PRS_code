# -*- coding: utf-8 -*-
"""Reviewer 1, point 6: recover the SuperPC entry for Fig 2B.

The original run left SuperPC empty because the reference implementation failed
inside every fold. This script re-runs the identical 10-fold loop with the
convergent SuperPC implementation and writes TCGA_ML_CV_cindex_v2.csv, which
adds the recovered row without touching the original table.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data", "processed")

import core  # noqa: E402
from analysis_tcga import load  # noqa: E402
from revision_reanalyses import superpc_fixed  # noqa: E402


def superpc_score(Xtr, ttr, etr, Xva):
    """Fit supervised PC1 on the training fold and score the held-out fold."""
    from lifelines import CoxPHFitter
    from sklearn.decomposition import PCA

    ok = np.isfinite(np.asarray(ttr, float)) & np.isfinite(np.asarray(etr, float))
    Xtr, ttr, etr = Xtr[ok], np.asarray(ttr, float)[ok], np.asarray(etr, float)[ok]
    genes = list(Xtr.columns)
    scores = []
    d = pd.DataFrame(Xtr.values, columns=genes)
    d["T"], d["E"] = ttr, etr.astype(int)
    for g in genes:
        try:
            cph = CoxPHFitter(penalizer=1e-3)
            cph.fit(d[[g, "T", "E"]], duration_col="T", event_col="E")
            scores.append(abs(float(cph.params_[g])))
        except Exception:
            scores.append(0.0)
    scores = np.asarray(scores, float)
    thr = np.percentile(scores, 70)
    sel = [genes[i] for i in range(len(genes)) if scores[i] >= max(thr, 1e-9)]
    if len(sel) < 3:
        sel = list(pd.Series(Xtr.var()).sort_values(ascending=False).head(15).index)

    mu = Xtr[sel].mean()
    sd = Xtr[sel].std().replace(0, 1e-9)
    Ztr = ((Xtr[sel] - mu) / sd).values
    Zva = ((Xva[sel] - mu) / sd).values
    pca = PCA(n_components=1).fit(Ztr)
    p1tr, p1va = pca.transform(Ztr)[:, 0], pca.transform(Zva)[:, 0]
    dd = pd.DataFrame({"PC1": p1tr, "T": ttr, "E": etr.astype(int)})
    cph = CoxPHFitter(penalizer=1e-3).fit(dd, "T", "E")
    beta = float(cph.params_["PC1"])
    return p1va * beta


def main():
    E, cl, _ = load()
    ucox = pd.read_csv(os.path.join(DATA, "TCGA_univariate_Cox_PRG.csv"))
    cand = ucox[ucox["p"] < 0.10]["gene"].tolist()
    X = E.loc[cand].T.reindex(cl.index).astype(float)
    X = X.fillna(X.median())
    X = (X - X.mean()) / (X.std() + 1e-9)
    X = X.loc[:, X.std() > 0]
    time, event = cl["os_years"].values, cl["os_event"].values
    n = len(cl)

    rng = np.random.default_rng(42)
    fold = np.zeros(n, int)
    idx = rng.permutation(n)
    for f in range(10):
        fold[idx[f::10]] = f

    cv = []
    for f in range(10):
        tr, va = fold != f, fold == f
        try:
            sc = superpc_score(X.iloc[tr], time[tr], event[tr], X.iloc[va])
            c, _ = core.cindex(time[va], event[va], sc)
            if not np.isnan(c):
                cv.append(c)
        except Exception as e:  # pragma: no cover
            print("fold", f, "failed:", str(e)[:70])
    print(f"SuperPC folds={len(cv)} mean={np.mean(cv) if cv else float('nan'):.4f} "
          f"sd={np.std(cv) if cv else float('nan'):.4f}")

    base = pd.read_csv(os.path.join(DATA, "TCGA_ML_CV_cindex.csv"))
    v2 = base.copy()
    v2.loc[v2["model"] == "SuperPC", "cv_cindex"] = np.mean(cv) if cv else np.nan
    v2.loc[v2["model"] == "SuperPC", "sd"] = np.std(cv) if cv else np.nan
    v2 = v2.sort_values("cv_cindex", ascending=False)
    v2.to_csv(os.path.join(DATA, "TCGA_ML_CV_cindex_v2.csv"), index=False)
    print(v2.to_string(index=False))
    json.dump({"superpc_n_folds": len(cv),
               "superpc_cv_cindex": float(np.mean(cv)) if cv else None,
               "superpc_sd": float(np.std(cv)) if cv else None},
              open(os.path.join(DATA, "superpc_cv.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
