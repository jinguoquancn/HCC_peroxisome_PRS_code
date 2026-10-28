# -*- coding: utf-8 -*-
"""Single-cell figure generation (Figures 1-3) for the HCC peroxisome paper.
Reads outputs from sc_analysis.py + annotated h5ads. Produces:
- Fig1: scRNA atlas overview, cell type dot plot, PXS distribution, PXS by dataset/celltype
- Fig2: malignant heterogeneity (PXS heatmap, CNV, pseudotime)
- Fig3: cell-cell communication, immune states"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, sys, json, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

for p in (str(PROJECT_ROOT), str(PROJECT_ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)
import figstyle as fs
from figstyle import save, label_axes, OI, MM, FULL_W, HALF_W, stars

OUT = str(PROJECT_ROOT / "data" / "processed")
SC = str(PROJECT_ROOT / "data" / "processed" / "sc")
DATA = str(PROJECT_ROOT / "data" / "genesets")

DSETS = ["GSE149614", "GSE125449", "GSE146115", "GSE151530", "GSE98638"]


def load_annot(ds):
    import scanpy as sc
    f = os.path.join(SC, f"{ds}_annot.h5ad")
    if not os.path.exists(f):
        return None
    return sc.read_h5ad(f)


def fig1():
    """Overview: cell-type distribution, PXS by dataset/celltype."""
    fig, axes = plt.subplots(2, 2, figsize=(FULL_W, 172 * MM),
                             gridspec_kw={"width_ratios": [1.30, 0.70]})
    # A. Stacked bar of cell types per dataset
    ax = axes[0][0]
    rows = []
    for ds in DSETS:
        a = load_annot(ds)
        if a is None:
            continue
        ct = a.obs["celltype"].astype(str).value_counts()
        for c, n in ct.items():
            rows.append({"dataset": ds, "celltype": c, "n": n})
    df = pd.DataFrame(rows)
    pv = df.pivot(index="dataset", columns="celltype", values="n").fillna(0)
    pv = pv.div(pv.sum(1), axis=0)
    pv.plot(kind="bar", stacked=True, ax=ax, legend=False, colormap="tab20",
            width=0.7)
    ax.set_ylabel("Proportion")
    ax.set_title("Cell type composition across 5 scRNA datasets")
    ax.tick_params(axis="x", rotation=45, labelsize=6)
    # custom legend: single column, OUTSIDE the panel's right side (no overlap)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="center left", bbox_to_anchor=(1.005, 0.5),
              frameon=False, fontsize=5, ncol=1, handlelength=1.0,
              handletextpad=0.4, borderaxespad=0.2)

    # B. PXS by cell type (mean ± SE)
    ax = axes[0][1]
    cdf = pd.read_csv(os.path.join(SC, "pxs_by_celltype.csv"))
    cts = sorted(set(cdf["celltype"]))
    means = cdf.groupby("celltype")["mean_pxs"].mean().reindex(cts)
    se = cdf.groupby("celltype")["mean_pxs"].sem().reindex(cts).fillna(0)
    cols = [OI[0] if c == "Malignant" else OI[1] if c == "Hepatocyte" else "#888888" for c in cts]
    # symmetric error bars (mean ± SE) with the y-axis starting at 0
    ax.barh(cts, means.values, xerr=se.values, color=cols, edgecolor="white",
            linewidth=0.5, error_kw=dict(elinewidth=0.8, capsize=2, ecolor="#444444"))
    ax.set_xlabel("Mean peroxisome score (z)")
    ax.set_title("Peroxisome activity by cell type (meta-analysis)")
    hi = (means + se).max()
    ax.set_xlim(0, hi * 1.15)

    # C. PXS per cell (sample scatter)
    ax = axes[1][0]
    all_rows = []
    for ds in DSETS:
        a = load_annot(ds)
        if a is None or "pxs" not in a.obs.columns:
            continue
        for ct in ["Malignant", "Hepatocyte", "T", "Myeloid", "Endothelial", "Fibroblast"]:
            if ct not in a.obs["celltype"].unique():
                continue
            m = a.obs["celltype"] == ct
            all_rows.append(pd.DataFrame({
                "pxs": a.obs.loc[m, "pxs"].values,
                "ct": ct, "dataset": ds
            }))
    if all_rows:
        dd = pd.concat(all_rows, ignore_index=True)
        order = ["Malignant", "Hepatocyte", "T", "Myeloid", "Endothelial", "Fibroblast"]
        dd["ct"] = pd.Categorical(dd["ct"], order)
        for i, c in enumerate(order):
            sub = dd[dd["ct"] == c]
            x = np.random.normal(i, 0.08, len(sub))
            ax.scatter(x, sub["pxs"], s=2, alpha=0.3, c="black")
        ax.boxplot([dd[dd["ct"] == c]["pxs"] for c in order], positions=range(len(order)),
                   widths=0.6, patch_artist=True,
                   medianprops=dict(color="white", lw=0.6),
                   boxprops=dict(lw=0.4), whiskerprops=dict(lw=0.4), showfliers=False)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=45, ha="right", fontsize=5)
        ax.set_ylabel("PXS (per cell)")
        ax.set_title("PXS per cell across datasets")
        # Malignant vs every other cell type -> staircase bracket + stars/ns
        from scipy.stats import mannwhitneyu
        m_v = dd[dd["ct"] == "Malignant"]["pxs"].dropna()
        gvals = [dd[dd["ct"] == c]["pxs"].dropna() for c in order]
        ymin, ymax = dd["pxs"].min(), dd["pxs"].max()
        step = (ymax - ymin) * 0.10
        for gi in range(1, len(order)):
            if len(gvals[gi]) < 5:
                continue
            _, p = mannwhitneyu(m_v, gvals[gi])
            yb = ymax + step * gi
            fs.sig_bracket(ax, 0, gi, yb, p, fs=6)
        ax.set_ylim(ymin - (ymax - ymin) * 0.05, ymax + step * (len(order) + 1))

    # D. Malignant vs normal in GSE149614 (site)
    ax = axes[1][1]
    f = os.path.join(SC, "GSE149614_pxs_tumor_vs_adjacent.csv")
    if os.path.exists(f):
        df = pd.read_csv(f)
        sites = ["Tumor", "Normal", "PVTT", "Lymph"]
        cts = ["Malignant", "Hepatocyte", "T", "Myeloid", "Endothelial", "Fibroblast"]
        M = df.pivot(index="celltype", columns="site", values="mean_pxs").reindex(cts)
        im = ax.imshow(M[sites].values, cmap="RdYlBu_r", aspect="auto")
        ax.set_xticks(range(len(sites)))
        ax.set_xticklabels(sites)
        ax.set_yticks(range(len(cts)))
        ax.set_yticklabels(cts, fontsize=6)
        ax.set_title("GSE149614: PXS by site × celltype")
        fig.colorbar(im, ax=ax, fraction=0.03)

    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.3)

    label_axes(axes, "abcd")
    save(fig, "Fig1_scRNA_overview")


def fig2():
    """Malignant cell heterogeneity: PXS heatmap in CNV subgroups, pseudotime."""
    fig, axes = plt.subplots(1, 3, figsize=(FULL_W, 90 * MM),
                             gridspec_kw={"width_ratios": [1.35, 1.05, 1.05]})
    # A. CNV amplitude by cell type (GSE151530, malignant-rich)
    ax = axes[0]
    f = os.path.join(SC, "GSE151530_cnv.csv")
    if os.path.exists(f):
        d = pd.read_csv(f)
        cts = ["Malignant", "Hepatocyte", "T", "Myeloid", "Endothelial", "Fibroblast"]
        means = d.groupby("celltype")["cnv_amp"].mean().reindex(cts)
        ax.bar(means.index, means.values,
               color=[OI[0] if c == "Malignant" else OI[1] for c in cts])
        ax.set_xticklabels(cts, rotation=55, ha="right", fontsize=5)
        ax.set_ylabel("Mean CNV amplitude")
        ax.set_title("GSE151530: CNV by cell type")
        # significance vs Malignant for every other cell type (staircase brackets)
        from scipy.stats import mannwhitneyu
        mv = d[d["celltype"] == "Malignant"]["cnv_amp"].dropna()
        gv = [d[d["celltype"] == c]["cnv_amp"].dropna() for c in cts]
        ymin, ymax = d["cnv_amp"].min(), d["cnv_amp"].max()
        step = (ymax - ymin) * 0.12
        for gi, c in enumerate(cts):
            if c == "Malignant" or len(gv[gi]) < 5:
                continue
            _, p = mannwhitneyu(mv, gv[gi])
            yb = ymax + step * gi
            fs.sig_bracket(ax, 0, gi, yb, p, fs=5)
        ax.set_ylim(ymin - (ymax - ymin) * 0.1, ymax + step * (len(cts) + 1))

    # B. Pseudotime scatter (GSE151530 malignant)
    ax = axes[1]
    f = os.path.join(SC, "GSE151530_malignant_pseudotime.csv")
    if os.path.exists(f):
        d = pd.read_csv(f)
        ax.scatter(d["pseudotime"], d["pxs"], s=1, alpha=0.12, c="black", linewidths=0)
        from scipy.stats import spearmanr
        rho, p = spearmanr(d["pseudotime"], d["pxs"])
        ax.set_xlabel("Pseudotime")
        ax.set_ylabel("PXS")
        ax.set_title(f"GSE151530 malignant pseudotime (n={len(d)}, rho={rho:.2f})")
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))

    # C. PXS density along pseudotime
    ax = axes[2]
    if os.path.exists(f):
        d = pd.read_csv(f).sort_values("pseudotime")
        n = max(20, len(d) // 50)
        bins = np.array_split(d, n)
        m = np.array([b["pxs"].mean() for b in bins])
        pt = np.array([b["pseudotime"].mean() for b in bins])
        ax.fill_between(pt, m, color=OI[0], alpha=0.7)
        ax.set_xlabel("Pseudotime")
        ax.set_ylabel("Mean PXS (binned)")
        ax.set_title("PXS evolution along pseudotime (GSE151530)")

    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.3)

    label_axes(axes, "abc")
    save(fig, "Fig6_malignant_heterogeneity")


def fig3():
    """Cell-cell communication + immune states."""
    fig, axes = plt.subplots(2, 2, figsize=(FULL_W, 140 * MM))
    # A. Communication network (sender -> receiver for Malignant->Myeloid in tumor)
    ax = axes[0][0]
    f = os.path.join(SC, "GSE149614_communication.csv")
    if os.path.exists(f):
        d = pd.read_csv(f)
        d_t = d[d["site"] == "Tumor"]
        # top LR pairs Malignant->Myeloid
        m_to_mye = d_t[(d_t["sender"] == "Malignant") & (d_t["receiver"] == "Myeloid")]
        m_to_mye = m_to_mye.sort_values("ligand_expr", ascending=False).head(15)
        ax.barh(m_to_mye["pathway"], m_to_mye["ligand_expr"], color=OI[0])
        ax.set_xlabel("Malignant ligand expr")
        ax.set_title("Malignant → Myeloid signaling (tumor)")

    # B. Myeloid -> T signaling
    ax = axes[0][1]
    if os.path.exists(f):
        d_t = pd.read_csv(f)
        d_t = d_t[d_t["site"] == "Tumor"]
        m_to_t = d_t[(d_t["sender"] == "Myeloid") & (d_t["receiver"] == "T")]
        m_to_t = m_to_t.sort_values("ligand_expr", ascending=False).head(15)
        ax.barh(m_to_t["pathway"], m_to_t["ligand_expr"], color=OI[1])
        ax.set_xlabel("Myeloid ligand expr")
        ax.set_title("Myeloid → T signaling (tumor)")

    # C. Immune states heatmap
    ax = axes[1][0]
    f = os.path.join(SC, "immune_states_scores.csv")
    if os.path.exists(f):
        d = pd.read_csv(f)
        pv = d.pivot(index="state", columns="celltype", values="mean_score").fillna(0)
        im = ax.imshow(pv.values, cmap="RdYlBu_r", aspect="auto")
        ax.set_xticks(range(len(pv.columns)))
        ax.set_xticklabels(pv.columns, rotation=45, ha="right", fontsize=5.5)
        ax.set_yticks(range(len(pv.index)))
        ax.set_yticklabels(pv.index, fontsize=6)
        ax.set_title("Immune state scores")
        fig.colorbar(im, ax=ax, fraction=0.03)
    # D. Mal-T LR pairs
    ax = axes[1][1]
    f_comm = os.path.join(SC, "GSE149614_communication.csv")
    if os.path.exists(f_comm):
        d = pd.read_csv(f_comm)
        d_t = d[d["site"] == "Tumor"]
        m_to_t2 = d_t[(d_t["sender"] == "Malignant") & (d_t["receiver"] == "T")]
        m_to_t2 = m_to_t2.sort_values("ligand_expr", ascending=False).head(15)
        ax.barh(m_to_t2["pathway"], m_to_t2["ligand_expr"], color=OI[0])
        ax.set_xlabel("Malignant ligand expr")
        ax.set_title("Malignant → T signaling (tumor)")

    fig.tight_layout(pad=1.0, h_pad=1.5, w_pad=1.3)

    label_axes(axes, "abcd")
    save(fig, "Fig7_communication_immune")


if __name__ == "__main__":
    for fn in [fig1, fig2, fig3]:
        try:
            fn()
        except Exception as e:
            print(fn.__name__, "ERR", e)
    print("SC FIGURES DONE")
