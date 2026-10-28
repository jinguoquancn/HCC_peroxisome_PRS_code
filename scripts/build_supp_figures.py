# -*- coding: utf-8 -*-
"""Generate the S1-S12 Supporting Information figures for the PLOS ONE
manuscript. All panels are drawn from real intermediate results under
data/processed (and the sc/ subfolder); no synthetic values are injected.

Outputs are written to figures/supplementary/ as TIF (300 dpi, LZW) and PNG.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, sys, json, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
warnings.filterwarnings("ignore")

ROOT = str(PROJECT_ROOT)
for p in (ROOT, os.path.join(ROOT, "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)
import figstyle as fs
from figstyle import save, label_axes, OI, MM, FULL_W, HALF_W, stars

DATA = os.path.join(ROOT, "data/processed")
SC = os.path.join(DATA, "sc")
OUTDIR = os.path.join(ROOT, "figures/supplementary")
os.makedirs(OUTDIR, exist_ok=True)
fs.OUT = OUTDIR   # redirect figstyle.save() to the supplementary folder

DSETS = ["GSE149614", "GSE125449", "GSE146115", "GSE151530", "GSE98638"]


def load_obs(ds):
    import scanpy as sc
    f = os.path.join(SC, f"{ds}_annot.h5ad")
    if not os.path.exists(f):
        return None
    a = sc.read_h5ad(f, backed="r")
    ob = a.obs[["n_genes", "pct_mt", "celltype"]].copy()
    ob["dataset"] = ds
    return ob


# ---------------------------------------------------------------- S1 Fig
def s1_qc_composition():
    """Single-cell QC metrics and dataset composition."""
    parts = [load_obs(d) for d in DSETS]
    parts = [p for p in parts if p is not None]
    ob = pd.concat(parts, ignore_index=True)
    fig, axes = plt.subplots(2, 2, figsize=(FULL_W, 140 * MM))

    # A. cells per dataset
    ax = axes[0][0]
    n = ob.groupby("dataset")["n_genes"].count().reindex(DSETS)
    bars = ax.bar(range(len(n)), n.values, color=[OI[i % 8] for i in range(len(n))], width=0.62)
    ax.set_xticks(range(len(n)))
    ax.set_xticklabels(n.index, rotation=40, ha="right", fontsize=6)
    ax.set_ylabel("Cells after QC")
    ax.set_title("Cells per scRNA-seq dataset")
    for i, v in enumerate(n.values):
        ax.text(i, v + v * 0.02, f"{int(v):,}", ha="center", va="bottom", fontsize=6)
    ax.text(0.98, 0.94, f"Total: {len(ob):,}", transform=ax.transAxes,
            ha="right", fontsize=7, fontweight="bold")

    # B. n_genes distribution
    ax = axes[0][1]
    vals = [ob.loc[ob["dataset"] == d, "n_genes"].dropna() for d in DSETS]
    bp = ax.boxplot(vals, positions=range(len(DSETS)), widths=0.55, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=0.8),
                    boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    for i, b in enumerate(bp["boxes"]):
        b.set_facecolor(OI[i % 8]); b.set_alpha(0.65)
    ax.set_xticks(range(len(DSETS)))
    ax.set_xticklabels(DSETS, rotation=40, ha="right", fontsize=6)
    ax.set_ylabel("Genes detected per cell")
    ax.set_title("Library complexity")

    # C. % mito distribution
    ax = axes[1][0]
    vals = [ob.loc[ob["dataset"] == d, "pct_mt"].dropna() for d in DSETS]
    bp = ax.boxplot(vals, positions=range(len(DSETS)), widths=0.55, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=0.8),
                    boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    for i, b in enumerate(bp["boxes"]):
        b.set_facecolor(OI[(i + 2) % 8]); b.set_alpha(0.65)
    ax.set_xticks(range(len(DSETS)))
    ax.set_xticklabels(DSETS, rotation=40, ha="right", fontsize=6)
    ax.set_ylabel("Mitochondrial %")
    ax.set_title("Mitochondrial fraction")

    # D. n_genes vs pct_mt scatter (downsampled)
    ax = axes[1][1]
    for d in DSETS:
        sub = ob.loc[ob["dataset"] == d]
        if len(sub) > 4000:
            sub = sub.sample(4000, random_state=7)
        ax.scatter(sub["n_genes"], sub["pct_mt"], s=1.2, alpha=0.25, label=d)
    ax.set_xlabel("Genes detected per cell")
    ax.set_ylabel("Mitochondrial %")
    ax.set_title("QC trade-off")
    ax.legend(frameon=False, fontsize=5.5, loc="upper right")

    fig.tight_layout(pad=1.0, h_pad=1.6, w_pad=1.3)
    label_axes(axes, "abcd")
    save(fig, "S1Fig_QC_composition")


# ---------------------------------------------------------------- S2 Fig
def s2_ml_cv():
    """Cross-validation performance of the eight ML algorithms and PRS coefficients."""
    cv = pd.read_csv(os.path.join(DATA, "TCGA_ML_CV_cindex.csv"))
    co = pd.read_csv(os.path.join(DATA, "TCGA_PRS_coefficients.csv"))
    co.columns = ["gene", "coef"]
    co = co.sort_values("coef")
    fig, axes = plt.subplots(1, 2, figsize=(FULL_W, 88 * MM),
                             gridspec_kw={"width_ratios": [1.15, 1.0]})

    ax = axes[0]
    ax.bar(range(len(cv)), cv["cv_cindex"], yerr=cv["sd"],
           color=OI[0], edgecolor="white", linewidth=0.5, capsize=2,
           error_kw=dict(elinewidth=0.8, ecolor="#333333"))
    ax.set_xticks(range(len(cv)))
    ax.set_xticklabels(cv["model"], rotation=40, ha="right", fontsize=6)
    ax.set_ylabel("Cross-validated C-index")
    ax.set_title("Eight-algorithm CV performance (10-fold)")
    ax.axhline(0.5, color="grey", lw=0.6, ls="--")
    ax.set_ylim(0.4, 0.85)
    for i, v in enumerate(cv["cv_cindex"]):
        ax.text(i, v + 0.012, f"{v:.3f}", ha="center", fontsize=5.5)

    ax = axes[1]
    cols = [OI[0] if c > 0 else OI[1] for c in co["coef"]]
    ax.barh(co["gene"], co["coef"], color=cols, edgecolor="white", linewidth=0.5)
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_xlabel("LASSO Cox coefficient")
    ax.set_title("PRS gene coefficients (n = 10)")
    for i, v in enumerate(co["coef"]):
        ax.text(v + 0.012 if v >= 0 else v - 0.012, i, f"{v:.3f}",
                va="center", ha="left" if v >= 0 else "right", fontsize=5.5)

    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.6)
    label_axes(axes, "ab")
    save(fig, "S2Fig_ML_CV")


# ---------------------------------------------------------------- S3 Fig
def compute_bulk_pxs(cohort):
    """Mean per-gene z-score of peroxisome genes across samples of a bulk
    cohort, with NaN-safe handling for sparse/missing expression values."""
    expr = pd.read_csv(os.path.join(DATA, cohort), index_col=0)
    genes = json.load(open(os.path.join(DATA, "peroxisome_gene_set.json")))["genes"]
    g = [x for x in genes if x in expr.index]
    E = expr.loc[g]
    E = E.dropna(how="all").fillna(0)              # NaN-safe: drop empty rows, zero-fill missing entries
    mu = E.mean(axis=1)
    sd = E.std(axis=1).replace(0, 1e-9).fillna(1e-9)
    z = E.sub(mu, axis=0).div(sd, axis=0)
    pxs = z.mean(axis=0).fillna(0).values
    return pxs, E.shape[0]


def s3_pxs_bulk():
    """Peroxisome module score distributions across bulk cohorts."""
    cohorts = {
        "TCGA-LIHC": (os.path.join(DATA, "TCGA_PXS_tumor_vs_normal.csv"), "builtin"),
        "GSE14520": ("GSE14520_expression.csv", "calc"),
        "GSE76427": ("GSE76427_expression.csv", "calc"),
        "GSE116174": ("GSE116174_expression.csv", "calc"),
        "GSE10186": ("GSE10186_expression.csv", "calc"),
        "ICGC-LIRI-JP": ("ICGC_LIRI_expression_tumor.csv", "calc"),
    }
    fig, axes = plt.subplots(1, 2, figsize=(FULL_W, 92 * MM),
                             gridspec_kw={"width_ratios": [1.5, 1.0]})
    names, dists = [], []
    for cname, (f, mode) in cohorts.items():
        if mode == "builtin":
            d = pd.read_csv(f)
            v = d["PXS"].dropna().values
        else:
            v, _ = compute_bulk_pxs(f)
        names.append(cname); dists.append(v)

    ax = axes[0]
    bp = ax.boxplot(dists, positions=range(len(names)), widths=0.6, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=0.8),
                    boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    for i, b in enumerate(bp["boxes"]):
        b.set_facecolor(OI[i % 8]); b.set_alpha(0.65)
        vals = dists[i]
        xs = np.random.normal(i, 0.06, min(len(vals), 500))
        ax.scatter(xs, np.random.choice(vals, min(len(vals), 500), replace=False),
                   s=1.5, color="black", alpha=0.18, zorder=3)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=30, ha="right", fontsize=6.5)
    ax.set_ylabel("PXS (mean z-score, 131 PRGs)")
    ax.set_title("Peroxisome activity across bulk cohorts")

    ax = axes[1]
    t = pd.read_csv(os.path.join(DATA, "TCGA_PXS_tumor_vs_normal.csv"))
    gv = [t.loc[t["group"] == g, "PXS"].values for g in ["T", "N"]]
    bp = ax.boxplot(gv, positions=[0, 1], widths=0.55, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=0.8),
                    boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    for i, b in enumerate(bp["boxes"]):
        b.set_facecolor([OI[0], OI[1]][i]); b.set_alpha(0.7)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Tumor", "Non-tumor"])
    ax.set_ylabel("PXS")
    ax.set_title("TCGA: tumor vs non-tumor")
    _, p = stats.mannwhitneyu(gv[0], gv[1])
    ymax = max(np.percentile(gv[0], 75), np.percentile(gv[1], 75))
    fs.sig_bracket(ax, 0, 1, ymax + (ymax - min(np.min(gv[0]), np.min(gv[1]))) * 0.06, p, fs=7)

    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.5)
    label_axes(axes, "ab")
    save(fig, "S3Fig_PXS_bulk")


# ---------------------------------------------------------------- S4 Fig
def s4_external_validation():
    """External validation details: KM curves for the four additional cohorts
    plus a forest plot of the five cohorts with survival data."""
    sumv = pd.read_csv(os.path.join(DATA, "external_validation_summary.csv"))
    has = sumv[sumv["HR"].notna()]
    fig = plt.figure(figsize=(FULL_W, 150 * MM))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1],
                          hspace=0.55, wspace=0.42)
    km_meta = {
        "GSE76427": ("GSE76427_validation_surv.csv", "GSE76427 (n=115)"),
        "GSE116174": ("GSE116174_validation_surv.csv", "GSE116174 (n=64)"),
        "GSE10186": ("GSE10186_validation_surv.csv", "GSE10186 (n=80)"),
        "ICGC-LIRI-JP": ("ICGC_LIRI_validation_surv.csv", "ICGC-LIRI-JP (n=203)"),
    }
    for k, (cname, (f, ttl)) in enumerate(km_meta.items()):
        d = pd.read_csv(os.path.join(DATA, f))
        med = d["risk"].median()
        grp = (d["risk"] > med).astype(int)
        row = has[has["cohort"] == cname]
        hr = row["HR"].values[0] if len(row) else None
        p = row["p_logrank"].values[0] if len(row) else None
        ax = fig.add_subplot(gs[k // 3, k % 3])
        fs.draw_km_from_surv(ax, d["t"], d["e"], grp, title=ttl, p=p, hr=hr)
        ax.set_xlabel("Time (years)")

    ax = fig.add_subplot(gs[1, 2])
    order = has.sort_values("HR")["cohort"]
    hr = has.set_index("cohort").loc[order, "HR"].values
    lo = has.set_index("cohort").loc[order, "HR_low"].values
    hi = has.set_index("cohort").loc[order, "HR_high"].values
    y = np.arange(len(order))
    ax.errorbar(hr, y, xerr=[hr - lo, hi - hr], fmt="D", color=OI[0], ms=4,
                capsize=3, elinewidth=0.9, markeredgecolor="white")
    ax.axvline(1, color="grey", ls="--", lw=0.7)
    ax.set_yticks(y); ax.set_yticklabels(order, fontsize=6.5)
    ax.set_xlabel("Hazard ratio (95% CI)")
    ax.set_title("Forest plot")
    ax.set_xlim(0, max(hi) * 1.15)

    fig.tight_layout(pad=1.0)
    # manual panel letters because axes created via add_subplot + add_gridspec
    letters = ["A", "B", "C", "D", "E"]
    axs = [fig.axes[0], fig.axes[1], fig.axes[2], fig.axes[3], fig.axes[4]]
    label_axes(axs, letters)
    save(fig, "S4Fig_external_validation")


# ---------------------------------------------------------------- S5 Fig
def s5_immune():
    """Immune infiltration ssGSEA scores and checkpoint expression by risk group."""
    imm = pd.read_csv(os.path.join(DATA, "immune_infiltration_PRS.csv"))
    ck = pd.read_csv(os.path.join(DATA, "checkpoints_PRS.csv"))
    piv = imm[imm["group"] != "p"].pivot(index="set", columns="group", values="mean")
    pvals = imm[imm["group"] == "p"].set_index("set")["mean"]
    sets = piv.index.tolist()
    fig, axes = plt.subplots(1, 2, figsize=(FULL_W, 100 * MM),
                             gridspec_kw={"width_ratios": [1.15, 1.0]})

    ax = axes[0]
    x = np.arange(len(sets))
    w = 0.36
    ax.bar(x - w / 2, piv["low"], w, color=OI[0], label="PRS low", alpha=0.85)
    ax.bar(x + w / 2, piv["high"], w, color=OI[1], label="PRS high", alpha=0.85)
    for i, s in enumerate(sets):
        p = pvals.get(s, 1)
        y = max(piv["low"][s], piv["high"][s])
        fs.sig_bracket(ax, i - w / 2, i + w / 2, y + 0.02, p, fs=5.5)
    ax.set_xticks(x)
    ax.set_xticklabels(sets, rotation=45, ha="right", fontsize=5.5)
    ax.set_ylabel("Mean ssGSEA score")
    ax.set_title("Immune infiltration by risk group")
    ax.legend(frameon=False, fontsize=6, loc="upper right")
    ymax = max(piv.max().max(), 0.05)
    ax.set_ylim(-0.35, ymax * 1.45)

    ax = axes[1]
    ck = ck.sort_values("logFC")
    cols = [OI[0] if v > 0 else OI[1] for v in ck["logFC"]]
    ax.barh(ck["gene"], ck["logFC"], color=cols, edgecolor="white", linewidth=0.5)
    for i, (v, p) in enumerate(zip(ck["logFC"], ck["p"])):
        if p < 0.05:
            ax.text(v + 0.01 if v >= 0 else v - 0.01, i, "*", va="center",
                    ha="left" if v >= 0 else "right", fontsize=8, fontweight="bold")
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_xlabel("Mean difference (high - low)")
    ax.set_title("Checkpoint genes")

    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.6)
    label_axes(axes, "ab")
    save(fig, "S5Fig_immune")


# ---------------------------------------------------------------- S6 Fig
def s6_mutations():
    """Mutational landscape and TP53 association by risk group."""
    m = pd.read_csv(os.path.join(DATA, "mutation_PRS.csv"))
    genes = [g for g in m["gene"].unique() if g != "TMB"]
    rows = []
    for g in genes:
        lo = m[(m["gene"] == g) & (m["group"] == "low")]
        hi = m[(m["gene"] == g) & (m["group"] == "high")]
        if len(lo) and len(hi):
            fr_lo = lo["n_mut"].iloc[0] / lo["n_total"].iloc[0]
            fr_hi = hi["n_mut"].iloc[0] / hi["n_total"].iloc[0]
            tt = [[hi["n_mut"].iloc[0], hi["n_total"].iloc[0] - hi["n_mut"].iloc[0]],
                  [lo["n_mut"].iloc[0], lo["n_total"].iloc[0] - lo["n_mut"].iloc[0]]]
            _, p = stats.fisher_exact(tt)
            rows.append({"gene": g, "fr_lo": fr_lo, "fr_hi": fr_hi, "p": p})
    dg = pd.DataFrame(rows).sort_values("fr_hi", ascending=False)
    fig, axes = plt.subplots(1, 2, figsize=(FULL_W, 96 * MM))

    ax = axes[0]
    x = np.arange(len(dg))
    w = 0.36
    ax.bar(x - w / 2, dg["fr_lo"], w, color=OI[0], label="PRS low", alpha=0.85)
    ax.bar(x + w / 2, dg["fr_hi"], w, color=OI[1], label="PRS high", alpha=0.85)
    for i, (_, r) in enumerate(dg.iterrows()):
        y = max(r["fr_lo"], r["fr_hi"])
        if r["p"] < 0.05:
            fs.sig_bracket(ax, i - w / 2, i + w / 2, y + 0.02, r["p"], fs=6)
    ax.set_xticks(x); ax.set_xticklabels(dg["gene"], rotation=45, ha="right")
    ax.set_ylabel("Mutation frequency")
    ax.set_title("Recurrent mutations")
    ax.legend(frameon=False, fontsize=6, loc="upper right")
    ax.set_ylim(0, dg[["fr_lo", "fr_hi"]].max().max() * 1.45)

    ax = axes[1]
    tmb = m[m["gene"] == "TMB"]
    vals = dict(zip(tmb["group"], tmb["n_mut"]))
    ax.bar([0, 1], [vals.get("low", np.nan), vals.get("high", np.nan)],
           color=[OI[0], OI[1]], width=0.5)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["PRS low", "PRS high"])
    ax.set_ylabel("Mean TMB")
    ax.set_title("Tumor mutational burden")

    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.6)
    label_axes(axes, "ab")
    save(fig, "S6Fig_mutations")


# ---------------------------------------------------------------- S7 Fig
def s7_gsea():
    """Hallmark GSEA-lite results (BH-FDR)."""
    g = pd.read_csv(os.path.join(DATA, "GSEA_PRS_high_vs_low.csv"))
    g["nes_dir"] = np.sign(g["mean_logFC"])
    sig = g.sort_values("FDR")
    up = sig[sig["nes_dir"] > 0].head(10)
    dn = sig[sig["nes_dir"] < 0].head(10)
    fig, axes = plt.subplots(1, 2, figsize=(FULL_W, 110 * MM))

    ax = axes[0]
    up = up.iloc[::-1]
    ax.barh(up["pathway"], up["mean_logFC"], color=OI[0], edgecolor="white", linewidth=0.5)
    for i, (_, r) in enumerate(up.iterrows()):
        ax.text(0.0002, i, f"FDR={r['FDR']:.1e}", va="center", fontsize=5)
    ax.set_xlabel("Mean log2FC (high vs low)")
    ax.set_title("Top enriched in PRS-high (n=10)")
    ax.set_xlim(0, up["mean_logFC"].max() * 1.6)

    ax = axes[1]
    dn = dn.iloc[::-1]
    ax.barh(dn["pathway"], dn["mean_logFC"], color=OI[1], edgecolor="white", linewidth=0.5)
    for i, (_, r) in enumerate(dn.iterrows()):
        ax.text(-0.0002, i, f"FDR={r['FDR']:.1e}", va="center", ha="right", fontsize=5)
    ax.set_xlabel("Mean log2FC (high vs low)")
    ax.set_title("Top depleted in PRS-high (n=10)")
    ax.set_xlim(dn["mean_logFC"].min() * 1.6, 0)

    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.6)
    label_axes(axes, "ab")
    save(fig, "S7Fig_GSEA")


# ---------------------------------------------------------------- S8 Fig
def s8_signature_comparison():
    """Head-to-head comparison of the PRS with published signatures."""
    s = pd.read_csv(os.path.join(DATA, "signature_comparison.csv")).sort_values("C_index")
    fig, ax = plt.subplots(figsize=(HALF_W, 100 * MM))
    cols = [OI[0] if "PRS" in x else "#999999" for x in s["signature"]]
    ax.barh(s["signature"], s["C_index"], color=cols, edgecolor="white", linewidth=0.5)
    for i, v in enumerate(s["C_index"]):
        ax.text(v + 0.004, i, f"{v:.3f}", va="center", fontsize=6)
    ax.set_xlabel("C-index (TCGA-LIHC)")
    ax.set_title("PRS vs seven published signatures")
    ax.set_xlim(0, max(s["C_index"].max(), 0.6) * 1.15)
    fig.tight_layout(pad=1.0)
    save(fig, "S8Fig_signature_comparison")


# ---------------------------------------------------------------- S9 Fig
def s9_malignant_axis():
    """Malignant-cell ordering axis and CNV inference (GSE151530)."""
    pt = pd.read_csv(os.path.join(SC, "GSE151530_malignant_pseudotime.csv"))
    cn = pd.read_csv(os.path.join(SC, "GSE151530_cnv.csv"))
    fig, axes = plt.subplots(1, 3, figsize=(FULL_W, 92 * MM),
                             gridspec_kw={"width_ratios": [1.05, 1.0, 0.95]})

    ax = axes[0]
    ax.scatter(pt["pseudotime"], pt["pxs"], s=1, alpha=0.12, color="black", linewidths=0)
    rho, p = stats.spearmanr(pt["pseudotime"], pt["pxs"])
    ax.set_xlabel("Peroxisome-correlated ordering axis")
    ax.set_ylabel("PXS")
    ax.set_title(f"Malignant cells (n={len(pt)}, rho={rho:.2f})")
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))

    ax = axes[1]
    d = pt.sort_values("pseudotime")
    n = max(20, len(d) // 50)
    bins = np.array_split(d, n)
    m = np.array([b["pxs"].mean() for b in bins])
    pv = np.array([b["pseudotime"].mean() for b in bins])
    ax.fill_between(pv, m, color=OI[0], alpha=0.7)
    ax.set_xlabel("Peroxisome-correlated ordering axis")
    ax.set_ylabel("Mean PXS (binned)")
    ax.set_title("PXS along the ordering axis")

    ax = axes[2]
    cn["malig"] = (cn["celltype"] == "Malignant").astype(int)
    m_ = cn[cn["malig"] == 1]["cnv_amp"]
    n_ = cn[cn["malig"] == 0]["cnv_amp"]
    ax.bar([0, 1], [n_.mean(), m_.mean()], color=[OI[1], OI[0]], width=0.5,
           yerr=[[n_.sem(), m_.sem()], [n_.sem(), m_.sem()]], capsize=3,
           error_kw=dict(elinewidth=0.8, ecolor="black"))
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Non-malignant", "Malignant"])
    ax.set_ylabel("CNV score")
    ax.set_title(f"Malignant {m_.mean():.3f} vs {n_.mean():.3f}")
    _, p = stats.mannwhitneyu(n_, m_)
    fs.sig_bracket(ax, 0, 1, max(n_.max(), m_.max()) + 0.02, p, fs=7)

    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.4)
    label_axes(axes, "abc")
    save(fig, "S9Fig_malignant_axis")


# ---------------------------------------------------------------- S10 Fig
def s10_communication():
    """Cell-cell communication network and immune-state scores."""
    comm = pd.read_csv(os.path.join(SC, "GSE149614_communication.csv"))
    imm = pd.read_csv(os.path.join(SC, "immune_states_scores.csv"))
    fig = plt.figure(figsize=(FULL_W, 120 * MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.3)

    ax = fig.add_subplot(gs[0, 0])
    dt = comm[comm["site"] == "Tumor"]
    agg = dt.groupby("pathway")["ligand_expr"].mean().sort_values(ascending=False).head(15)
    ax.barh(agg.index[::-1], agg.values[::-1], color=OI[0], edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Mean ligand expression")
    ax.set_title("Top signaling pathways (tumor)")

    ax = fig.add_subplot(gs[0, 1])
    pv = imm.pivot(index="state", columns="celltype", values="mean_score").fillna(0)
    im = ax.imshow(pv.values, cmap="RdYlBu_r", aspect="auto")
    ax.set_xticks(range(len(pv.columns)))
    ax.set_xticklabels(pv.columns, rotation=45, ha="right", fontsize=5.5)
    ax.set_yticks(range(len(pv.index)))
    ax.set_yticklabels(pv.index, fontsize=6)
    ax.set_title("Immune state scores")
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)

    fig.tight_layout(pad=1.0)
    label_axes([fig.axes[0], fig.axes[1]], ["A", "B"])
    save(fig, "S10Fig_communication")


# ---------------------------------------------------------------- S11 Fig
def s11_tcga():
    """TCGA risk stratification and tumor-non-tumor peroxisome contrast."""
    cl = pd.read_csv(os.path.join(DATA, "TCGA_LIHC_clinical_with_PRS.csv"))
    # NOTE: TCGA_LIHC_clinical_with_PRS.csv has 'tissue' column uniformly set
    # to 'normal' in the cached file; every row is in fact a tumor sample with
    # PRS/PXS/survival. We therefore filter on the survival columns directly.
    tum = cl.dropna(subset=["os_years", "os_event", "PRS_high"]).copy()
    tum = tum[tum["PRS_high"].notna()]
    px = pd.read_csv(os.path.join(DATA, "TCGA_PXS_tumor_vs_normal.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(FULL_W, 92 * MM))

    ax = axes[0]
    g = tum["PRS_high"].astype(int).values
    fs.draw_km_from_surv(ax, tum["os_years"], tum["os_event"], g,
                         title="TCGA-LIHC (n=%d)" % len(tum),
                         xlabel="Time (years)")
    fig.axes[0].texts  # noqa

    ax = axes[1]
    gv = [px.loc[px["group"] == g2, "PXS"].values for g2 in ["T", "N"]]
    bp = ax.boxplot(gv, positions=[0, 1], widths=0.55, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=0.8),
                    boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    for i, b in enumerate(bp["boxes"]):
        b.set_facecolor([OI[0], OI[1]][i]); b.set_alpha(0.7)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Tumor", "Non-tumor"])
    ax.set_ylabel("PXS")
    ax.set_title("Peroxisome activity")
    _, p = stats.mannwhitneyu(gv[0], gv[1])
    fs.sig_bracket(ax, 0, 1, max(np.percentile(gv[0], 75), np.percentile(gv[1], 75)) + 0.1, p, fs=7)

    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.6)
    label_axes(axes, "ab")
    save(fig, "S11Fig_TCGA_risk")


# ---------------------------------------------------------------- S12 Fig
def s12_gdsc():
    """GDSC predicted drug sensitivity by risk group."""
    d = pd.read_csv(os.path.join(DATA, "drug_sensitivity_PRS.csv"))
    drugs = sorted(d["drug"].unique())
    lo = d[d["group"] == "low"].set_index("drug")["mean_pred_IC50"].reindex(drugs)
    hi = d[d["group"] == "high"].set_index("drug")["mean_pred_IC50"].reindex(drugs)
    fig, ax = plt.subplots(figsize=(FULL_W, 110 * MM))
    x = np.arange(len(drugs))
    w = 0.36
    ax.bar(x - w / 2, lo.values, w, color=OI[0], label="PRS low", alpha=0.85)
    ax.bar(x + w / 2, hi.values, w, color=OI[1], label="PRS high", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels([dd.split("_")[0] for dd in drugs], rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("Predicted mean log(IC50)")
    ax.set_title("GDSC ridge-regression predictions (hypothesis-generating)")
    ax.legend(frameon=False, fontsize=6)
    fig.tight_layout(pad=1.0)
    save(fig, "S12Fig_GDSC")


if __name__ == "__main__":
    jobs = [s1_qc_composition, s2_ml_cv, s3_pxs_bulk, s4_external_validation,
            s5_immune, s6_mutations, s7_gsea, s8_signature_comparison,
            s9_malignant_axis, s10_communication, s11_tcga, s12_gdsc]
    for fn in jobs:
        try:
            fn()
        except Exception as e:
            print(f"{fn.__name__} ERR {e}")
    print("SUPP FIGURES DONE")
