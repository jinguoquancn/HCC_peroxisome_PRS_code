# -*- coding: utf-8 -*-
"""Generate bulk figures (Fig 4-8) for the HCC peroxisome paper.
Requires: figstyle.py, data/processed/* tables."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, sys, json, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import spearmanr

for p in (str(PROJECT_ROOT), str(PROJECT_ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)
import figstyle as fs
from figstyle import save, label_axes, OI, MM, FULL_W, HALF_W, roc_curve_from_auc, stars
import matplotlib.patheffects as pe

DATA = str(PROJECT_ROOT / "data" / "processed")


def fig4():
    """PXS tumor vs normal; ML CV; PRS coefficients; WGCNA modules."""
    fig, axes = plt.subplots(2, 2, figsize=(FULL_W, 150 * MM))
    # A. PXS tumor vs normal
    ax = axes[0][0]
    d = pd.read_csv(os.path.join(DATA, "TCGA_PXS_tumor_vs_normal.csv"))
    raw = {k: g["PXS"].values for k, g in d.groupby("group")}
    data = {"Tumor": raw["T"], "Non-tumor": raw["N"]}  # explicit order (fix: labels were swapped)
    fs.box_pairs(ax, data, colors=[OI[0], OI[1]])
    ax.set_title("Peroxisome activity score")
    p = stats.mannwhitneyu(data["Tumor"], data["Non-tumor"]).pvalue
    ymax = max([np.nanmax(v) for v in data.values()])
    fs.sig_bracket(ax, 0, 1, ymax * 1.05, p, fs=8)
    ax.set_ylim(ax.get_ylim()[0], ymax * 1.22)
    ax.set_ylabel("Peroxisome module score")
    ax.set_xlabel("")

    # B. ML CV C-index
    ax = axes[0][1]
    cv_path = os.path.join(DATA, "TCGA_ML_CV_cindex_v2.csv")
    if not os.path.exists(cv_path):
        cv_path = os.path.join(DATA, "TCGA_ML_CV_cindex.csv")
    cv = pd.read_csv(cv_path).dropna(subset=["cv_cindex"])
    cv = cv.sort_values("cv_cindex")
    cols = [OI[0] if m == "StepCox" else "#B0B0B0" for m in cv["model"]]
    ax.barh(cv["model"], cv["cv_cindex"], color=cols, height=0.6)
    ax.errorbar(cv["cv_cindex"], range(len(cv)), xerr=cv["sd"], fmt="none",
                ecolor="black", elinewidth=0.6, capsize=2)
    ax.axvline(0.5, color="grey", ls="--", lw=0.7)
    ax.set_xlabel("10-fold CV C-index")
    ax.set_title("Machine-learning model zoo")
    hi = (cv["cv_cindex"] + cv["sd"]).max()
    ax.set_xlim(0, hi * 1.12)

    # C. PRS coefficients
    ax = axes[1][0]
    coef = pd.read_csv(os.path.join(DATA, "TCGA_PRS_coefficients.csv"), index_col=0)["0"]
    coef = coef.sort_values()
    cols = [OI[1] if v < 0 else OI[0] for v in coef]
    ax.barh(coef.index, coef.values, color=cols, height=0.6)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("Coefficient (PRS)")
    ax.set_title("PRS genes (StepCox, n=%d)" % len(coef))

    # D. WGCNA modules: bars = correlation with PXS, labels = Cox HR for OS
    ax = axes[1][1]
    w = pd.read_csv(os.path.join(DATA, "TCGA_WGCNA_modules.csv"))
    w = w.sort_values("p_surv").head(12)
    x = np.arange(len(w))
    ax.bar(x, w["cor_PXS"], color=[OI[0] if c >= 0 else OI[1] for c in w["cor_PXS"]])
    for i, (hr, p) in enumerate(zip(w["HR"], w["p_surv"])):
        ax.text(i, w["cor_PXS"].iloc[i] + 0.012 * (1 if w["cor_PXS"].iloc[i] >= 0 else -1),
                "HR %.2f%s" % (float(np.exp(hr)), stars(p)), ha="center", fontsize=6,
                va="bottom" if w["cor_PXS"].iloc[i] >= 0 else "top",
                path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    ax.set_xticks(x)
    ax.set_xticklabels([f"M{m}" for m in w["module"]], rotation=45, ha="right")
    ax.axhline(0, color="grey", lw=0.6)
    ax.set_ylim(min(w["cor_PXS"]) * 1.35, max(w["cor_PXS"]) * 1.55)
    ax.set_ylabel("Spearman correlation with PXS (bars)")
    ax.set_title("WGCNA modules (text: Cox HR for OS)")


    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.3)

    label_axes([axes[0][0], axes[0][1], axes[1][0], axes[1][1]], "abcd")
    save(fig, "Fig2_bulk_PRS_building")


def fig5():
    """External validation: KM curves for the 5 survival cohorts + forest + ROC."""
    cohorts = ["GSE14520", "GSE76427", "GSE116174", "GSE10186", "ICGC-LIRI-JP"]
    fig, axes = plt.subplots(2, 3, figsize=(FULL_W, 138 * MM))
    ax_km = axes[0][0]
    # (A) GSE14520 KM as representative
    d = pd.read_csv(os.path.join(DATA, "GSE14520_validation_surv.csv"))
    g = (d["risk"] > d["risk"].median()).astype(int)
    fs.draw_km_from_surv(ax_km, d["t"], d["e"], g, title="GSE14520 (n=221)",
                         p=0.0126, hr=1.73, risk_table=True)
    # (B-E) other cohort KMs
    km_pos = {"GSE76427": (0, 1), "GSE116174": (0, 2), "GSE10186": (1, 0), "ICGC-LIRI-JP": (1, 1)}
    km_meta = {"GSE76427": ("GSE76427 (n=115)", 0.83, 0.91),
               "GSE116174": ("GSE116174 (n=64)", 0.157, 1.75),
               "GSE10186": ("GSE10186 (n=80)", 0.35, 1.41),
               "ICGC-LIRI-JP": ("ICGC-LIRI-JP (n=203)", 0.096, 1.82)}
    for name, (r, c) in km_pos.items():
        ax = axes[r][c]
        f = os.path.join(DATA, "ICGC_LIRI_validation_surv.csv" if name == "ICGC-LIRI-JP"
                         else f"{name}_validation_surv.csv")
        if not os.path.exists(f):
            ax.text(0.5, 0.5, "no survival", ha="center", transform=ax.transAxes)
            continue
        d = pd.read_csv(f)
        g = (d["risk"] > d["risk"].median()).astype(int)
        ttl, p, hr = km_meta[name]
        fs.draw_km_from_surv(ax, d["t"], d["e"], g, title=ttl, p=p, hr=hr,
                             risk_table=True)
    # (F) forest plot
    ax = axes[1][2]
    res = pd.read_csv(os.path.join(DATA, "external_validation_summary.csv"))
    res = res[res["km"] == True].sort_values("HR")
    y = np.arange(len(res))
    hr = res["HR"].values
    lo, hi = res["HR_low"].values, res["HR_high"].values
    ax.errorbar(hr, y, xerr=[hr - lo, hi - hr], fmt="s", color=OI[0], ms=4, capsize=3)
    ax.axvline(1, color="grey", ls="--", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(res["cohort"])
    xmax = max(hi) * 1.38
    ax.set_xlim(min(lo) * 0.85, xmax)
    ax.set_ylim(-0.6, len(res) - 0.4)
    for i, (p, n) in enumerate(zip(res["p_logrank"], res["n"])):
        ax.text(xmax * 0.97, i, "P=%.3f" % p, fontsize=6.5, va="center", ha="right",
                path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    # pooled random-effects estimate (DerSimonian-Laird) + heterogeneity
    yv = np.log(hr)
    sev = (np.log(hi) - np.log(lo)) / (2 * 1.96)
    wv = 1 / sev ** 2
    fe = np.sum(wv * yv) / np.sum(wv)
    Q = float(np.sum(wv * (yv - fe) ** 2))
    kk = len(yv)
    Cv = np.sum(wv) - np.sum(wv ** 2) / np.sum(wv)
    tau2 = max(0.0, (Q - (kk - 1)) / Cv) if Cv > 0 else 0.0
    I2 = max(0.0, (Q - (kk - 1)) / Q) * 100 if Q > 0 else 0.0
    wr = 1 / (sev ** 2 + tau2)
    re = np.sum(wr * yv) / np.sum(wr)
    se_re = float(np.sqrt(1 / np.sum(wr)))
    pH = float(np.exp(re))
    plo = float(np.exp(re - 1.96 * se_re))
    phi = float(np.exp(re + 1.96 * se_re))
    ypool = float(len(res))
    ax.fill([plo, pH, phi, pH], [ypool, ypool + 0.24, ypool, ypool - 0.24],
            color=OI[1], alpha=0.95, lw=0)
    ax.axvline(pH, color=OI[1], ls=":", lw=0.7)
    ax.set_ylim(-0.6, ypool + 1.05)
    ax.set_yticks(list(y) + [ypool])
    ax.set_yticklabels(list(res["cohort"]) + ["Pooled"])
    ax.text(xmax * 0.97, ypool + 0.62,
            "HR %.2f (%.2f\u2013%.2f)\nI\u00b2=%.0f%%" % (pH, plo, phi, I2),
            fontsize=6.0, va="center", ha="right", color=OI[1],
            path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    ax.set_xlabel("HR (high vs low PRS)")

    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.3)

    label_axes([axes[0][0], axes[0][1], axes[0][2], axes[1][0], axes[1][1], axes[1][2]],
               "abcdef")
    save(fig, "Fig3_external_validation")


def fig6():
    """Immune infiltration, checkpoints, mutation landscape."""
    fig, axes = plt.subplots(1, 3, figsize=(FULL_W, 90 * MM))
    # A. immune infiltration heatmap
    ax = axes[0]
    d = pd.read_csv(os.path.join(DATA, "immune_infiltration_PRS.csv"))
    piv = d.pivot(index="set", columns="group", values="mean")
    piv = piv[["low", "high"]]
    # normalize rows to z
    z = (piv.sub(piv.mean(1), axis=0).div(piv.std(1) + 1e-9, axis=0))
    im = ax.imshow(z.values, cmap="RdYlBu_r", aspect="auto")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["PRS-low", "PRS-high"])
    ax.set_yticks(range(len(z.index)))
    ax.set_yticklabels(z.index, fontsize=6)
    ax.set_title("Immune infiltration (ssGSEA)")
    cb = fig.colorbar(im, ax=ax, fraction=0.03)
    cb.ax.tick_params(labelsize=6)

    # B. checkpoint expression
    ax = axes[1]
    ck = pd.read_csv(os.path.join(DATA, "checkpoints_PRS.csv"))
    ck = ck.sort_values("logFC")
    cols = [OI[1] if v < 0 else OI[0] for v in ck["logFC"]]
    ax.barh(ck["gene"], ck["logFC"], color=cols)
    for i, p in enumerate(ck["p"]):
        if p >= 0.05:
            continue  # mark significant hits only, SCI convention
        ax.text(ck["logFC"].iloc[i] + (0.005 if ck["logFC"].iloc[i] >= 0 else -0.005),
                i, "*" if p < 0.05 else "", fontsize=9, va="center",
                ha="left" if ck["logFC"].iloc[i] >= 0 else "right")
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("log2FC (PRS-high vs low)")
    ax.set_title("Immune checkpoints")

    # C. mutation landscape
    ax = axes[2]
    m = pd.read_csv(os.path.join(DATA, "mutation_PRS.csv"))
    m = m[m["gene"] != "TMB"]
    genes = sorted(set(m["gene"]))
    low = m[m["group"] == "low"].set_index("gene")["n_mut"]
    high = m[m["group"] == "high"].set_index("gene")["n_mut"]
    nlo = m[m["group"] == "low"]["n_total"].iloc[0]
    nhi = m[m["group"] == "high"]["n_total"].iloc[0]
    x = np.arange(len(genes))
    w = 0.35
    ax.bar(x - w / 2, low.reindex(genes).values / nlo * 100, w, label="PRS-low", color=OI[1])
    ax.bar(x + w / 2, high.reindex(genes).values / nhi * 100, w, label="PRS-high", color=OI[0])
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=65, ha="right", fontsize=5)
    ax.set_ylabel("Mutation frequency (%)")
    ax.set_title("Mutation landscape")
    ax.legend(frameon=True, framealpha=0.92, loc="upper left",
              bbox_to_anchor=(1.01, 0.98), fontsize=6, ncol=1,
              handlelength=1.2)
    # significance brackets for Fisher-exact tests (group-pair comparison)
    from scipy.stats import fisher_exact
    ymax_all = max((low.reindex(genes) / nlo).max(), (high.reindex(genes) / nhi).max()) * 100
    step = ymax_all * 0.13
    for gi, gene in enumerate(genes):
        a = int(low.get(gene, 0)); b = nlo - a
        c = int(high.get(gene, 0)); d = nhi - c
        _, pv = fisher_exact([[a, b], [c, d]])
        ytop = max(a / nlo, c / nhi) * 100
        fs.sig_bracket(ax, gi - w / 2, gi + w / 2, ytop + step, pv, fs=7)
    ax.set_ylim(ax.get_ylim()[0], ymax_all * 1.35)

    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.3)

    label_axes(axes, "abc")
    save(fig, "Fig4_immune_mutation")


def fig7():
    """Drug sensitivity (predicted IC50 by risk group)."""
    d = pd.read_csv(os.path.join(DATA, "drug_sensitivity_PRS.csv"))
    piv = d.pivot(index="drug", columns="group", values="mean_pred_IC50")
    piv = piv[["low", "high"]].dropna()
    # direction: for drugs where high > low (resistance in high), plot raw delta
    fig, ax = plt.subplots(figsize=(HALF_W, 90 * MM))
    delta = piv["high"] - piv["low"]
    delta = delta.sort_values()
    cols = [OI[0] if v < 0 else OI[1] for v in delta]  # negative = more sensitive in high
    ax.barh(delta.index, delta.values, color=cols)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("Δ predicted IC50 (high - low)")
    ax.set_title("Predicted drug sensitivity (GDSC models)")

    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.3)

    label_axes([ax], "a")
    save(fig, "Fig7_drug_sensitivity")


def fig8():
    """Functional / clinical correlates: GSEA, signature comparison,
    prognostic predictors and predicted drug sensitivity (2x2 composite,
    incorporating the former single-panel drug-sensitivity figure)."""
    fig, axes = plt.subplots(2, 2, figsize=(FULL_W, 150 * MM))
    # A. GSEA top pathways
    ax = axes[0][0]
    g = pd.read_csv(os.path.join(DATA, "GSEA_PRS_high_vs_low.csv"))
    g = g.sort_values("p").head(10)
    cols = [OI[0] if v > 0 else OI[1] for v in g["mean_logFC"]]
    ax.barh(g["pathway"], g["mean_logFC"], color=cols)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("Mean log2FC (high vs low)")
    ax.set_title("GSEA (Hallmark)")
    ax.tick_params(axis="y", labelsize=6)

    # B. signature comparison
    ax = axes[0][1]
    s = pd.read_csv(os.path.join(DATA, "signature_comparison.csv")).sort_values("C_index")
    cols = [OI[0] if "PRS" in m else "#B0B0B0" for m in s["signature"]]
    ax.barh(s["signature"], s["C_index"], color=cols)
    ax.axvline(0.5, color="grey", ls="--", lw=0.7)
    ax.set_xlabel("C-index (OS)")
    ax.set_title("Head-to-head signature comparison")
    ax.tick_params(axis="y", labelsize=6)

    # C. nomogram-like cox forest (univariate predictors)
    ax = axes[1][0]
    cl = pd.read_csv(os.path.join(DATA, "nomogram_data.csv"))
    cl = cl.dropna(subset=["os_event"])
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
    import core
    preds = []
    hr, lo, hi, p, n = core.cox_hr(cl["os_years"], cl["os_event"],
                                   (cl["PRS"] > cl["PRS"].median()).astype(int))
    preds.append(("PRS (high vs low)", hr, lo, hi, p))
    hr, lo, hi, p, n = core.cox_hr(cl["os_years"], cl["os_event"],
                                   (cl["age"] > cl["age"].median()).astype(int))
    preds.append(("Age (old vs young)", hr, lo, hi, p))
    preds = preds[::-1]
    y = np.arange(len(preds))
    hr = np.array([p[1] for p in preds])
    lo = np.array([p[2] for p in preds]); hi = np.array([p[3] for p in preds])
    ax.errorbar(hr, y, xerr=[hr - lo, hi - hr], fmt="s", color=OI[0], ms=5, capsize=3)
    ax.axvline(1, color="grey", ls="--", lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([p[0] for p in preds])
    for i, p in enumerate(preds):
        if p[4] < 0.05:
            ax.text(hi[i] * 1.08, i, "*", fontsize=10, va="center",
                    path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    ax.set_xlim(min(lo) * 0.8, max(hi) * 1.3)
    ax.set_xlabel("Univariate HR")
    ax.set_title("Prognostic predictors (TCGA)")

    # D. predicted drug sensitivity (was Fig7)
    ax = axes[1][1]
    d = pd.read_csv(os.path.join(DATA, "drug_sensitivity_PRS.csv"))
    piv = d.pivot(index="drug", columns="group", values="mean_pred_IC50")
    piv = piv[["low", "high"]].dropna()
    delta = piv["high"] - piv["low"]
    delta = delta.sort_values()
    cols = [OI[0] if v < 0 else OI[1] for v in delta]
    ax.barh(delta.index, delta.values, color=cols)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("Δ predicted IC50 (high − low)")
    ax.set_title("Predicted drug sensitivity (GDSC)")

    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.3)

    label_axes(axes, "abcd")
    save(fig, "Fig5_functional_correlates")


def fig9():
    """Risk distribution + KM in TCGA and PXS association with TMB/drugs (summary panels)."""
    fig, axes = plt.subplots(1, 2, figsize=(FULL_W, 92 * MM))
    cl = pd.read_csv(os.path.join(DATA, "TCGA_LIHC_clinical_with_PRS.csv"))
    ax = axes[0]
    cl = cl.dropna(subset=["os_years", "os_event"])  # 370 of 371 have OS follow-up
    med = cl["PRS"].median()
    g = (cl["PRS"] > med).astype(int)
    # compute log-rank P and univariate Cox HR from the plotted data (was hardcoded)
    from lifelines.statistics import logrank_test
    from lifelines import CoxPHFitter
    lr = logrank_test(cl.loc[g == 0, "os_years"], cl.loc[g == 1, "os_years"],
                      cl.loc[g == 0, "os_event"], cl.loc[g == 1, "os_event"])
    dkm = pd.DataFrame({"T": cl["os_years"].values, "E": cl["os_event"].astype(int).values,
                        "G": g.values}).dropna()
    cph = CoxPHFitter(penalizer=1e-6).fit(dkm, "T", "E", formula="G")
    fs.draw_km_from_surv(ax, cl["os_years"], cl["os_event"], g,
                         title=f"TCGA-LIHC (n={len(cl)})",
                         p=float(lr.p_value), hr=float(np.exp(cph.params_["G"])),
                         risk_table=True)
    ax = axes[1]
    d = pd.read_csv(os.path.join(DATA, "TCGA_PXS_tumor_vs_normal.csv"))
    raw = {k: g2["PXS"].values for k, g2 in d.groupby("group")}
    data = {"Tumor": raw["T"], "Non-tumor": raw["N"]}  # explicit order (fix: labels were swapped)
    fs.box_pairs(ax, data, colors=[OI[0], OI[1]])
    ax.set_ylabel("Peroxisome module score")
    p = stats.mannwhitneyu(data["Tumor"], data["Non-tumor"]).pvalue
    ymax = max([np.nanmax(v) for v in data.values()])
    fs.sig_bracket(ax, 0, 1, ymax * 1.05, p, fs=8)
    ax.set_ylim(ax.get_ylim()[0], ymax * 1.22)

    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.3)

    label_axes(axes, "ab")
    save(fig, "Fig8_TCGA_KM_PRS")


if __name__ == "__main__":
    for fn in [fig4, fig5, fig6, fig8, fig9]:
        try:
            fn()
        except Exception as e:
            print(fn.__name__, "ERR", e)
    print("FIGURES DONE")
