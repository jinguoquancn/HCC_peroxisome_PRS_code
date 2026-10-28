# -*- coding: utf-8 -*-
"""Single-cell core analysis for the HCC peroxisome paper (Figure 1-3 data).
Loads annotated h5ads, computes peroxisome scores across cell types, CNV-based
malignant cell inference, pseudotime, cell-cell communication, immune states.
Saves all intermediate tables under data/processed/sc/."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, json, sys, warnings
import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import issparse
warnings.filterwarnings("ignore")

import scanpy as sc
sc.settings.verbosity = 0

OUT = str(PROJECT_ROOT / "data" / "processed" / "sc")
GENESETS = str(PROJECT_ROOT / "data" / "genesets")
DATA = str(PROJECT_ROOT / "data" / "processed")

# ------------------------------------------------------------------ helpers
def getX(adata, genes):
    g = [x for x in genes if x in adata.var_names]
    if not g:
        return None
    X = adata[:, g].X
    if issparse(X):
        X = X.toarray()
    return X, g


def score(adata, genes):
    r = getX(adata, genes)
    if r is None:
        return np.zeros(adata.n_obs)
    X, g = r
    # log-normalized already in raw h5ads; use mean expression as module score
    m = X.mean(axis=1)
    return m


def getpos():
    return json.load(open(os.path.join(GENESETS, "gene_pos.json")))


# ------------------------------------------------------------------ load annotated
def load_annot(ds):
    f = os.path.join(OUT, f"{ds}_annot.h5ad")
    if not os.path.exists(f):
        return None
    a = sc.read_h5ad(f)
    return a


def main():
    prg = json.load(open(os.path.join(DATA, "peroxisome_gene_set.json")))
    prg_genes = prg["genes"]
    dsets = ["GSE149614", "GSE125449", "GSE146115", "GSE151530", "GSE98638"]

    # ---- 1. per-cell peroxisome score & per celltype summary
    celltype_rows = []
    for ds in dsets:
        a = load_annot(ds)
        if a is None:
            print("[skip] annot missing", ds, flush=True)
            continue
        print(f"[load] {ds} {a.shape}", flush=True)
        if "pxs" not in a.obs.columns:
            s = score(a, prg_genes)
            a.obs["pxs"] = s
        ct = a.obs["celltype"].astype(str)
        pxs = a.obs["pxs"].values
        for c in sorted(set(ct)):
            m = ct == c
            if m.sum() < 10:
                continue
            celltype_rows.append({"dataset": ds, "celltype": c, "n": int(m.sum()),
                                  "mean_pxs": float(np.nanmean(pxs[m])),
                                  "frac_high": float(np.mean(pxs[m] > np.median(pxs)))})
    cdf = pd.DataFrame(celltype_rows)
    cdf.to_csv(os.path.join(OUT, "pxs_by_celltype.csv"), index=False)
    print("PXS by celltype saved", flush=True)

    # ---- 2. GSE149614 malignant vs non-malignant + tumor vs adjacent PXS
    a = load_annot("GSE149614")
    if a is not None:
        st = a.obs["site"].astype(str).values if "site" in a.obs.columns else ["tumor"] * a.n_obs
        ct = a.obs["celltype"].astype(str).values
        pxs = a.obs["pxs"].values
        rows = []
        for si in ["Tumor", "Normal", "PVTT", "Lymph"]:
            for c in ["Malignant", "Hepatocyte", "T", "Myeloid", "Endothelial", "Fibroblast"]:
                m = (st == si) & (ct == c)
                if m.sum() < 5:
                    continue
                rows.append({"site": si, "celltype": c, "n": int(m.sum()),
                             "mean_pxs": float(np.nanmean(pxs[m]))})
        pd.DataFrame(rows).to_csv(os.path.join(OUT, "GSE149614_pxs_tumor_vs_adjacent.csv"), index=False)
        print("GSE149614 site/celltype PXS saved", flush=True)

    # ---- 3. CNV inference for malignant cells (GSE149614 + GSE151530)
    pos = getpos()
    for ds in ["GSE149614", "GSE151530"]:
        a = load_annot(ds)
        if a is None:
            continue
        ct = a.obs["celltype"].astype(str).values
        mal = ct == "Malignant"
        # reference = non-immune, non-malignant (hepatocytes when possible, else all non-malignant epithelial)
        cts = pd.Series(ct)
        ref = cts.isin(["Hepatocyte", "Fibroblast", "Endothelial"])
        # infer CNV: mean expression of windows of genes on each chromosome
        # build gene->chromosome ordered list
        genes = [str(x) for x in a.var_names]
        ginfo = [(g, pos.get(g)) for g in genes if pos.get(g)]
        # only use protein-coding-ish expressed genes
        ginfo = [(g, p) for g, p in ginfo if p and p[0] in [str(i) for i in range(1, 23)]]
        ginfo.sort(key=lambda x: (int(x[1][0]), x[1][1]))
        chr_genes = {}
        for g, p in ginfo:
            chr_genes.setdefault(p[0], []).append(g)
        # per-cell CNV score = mean of z-scored expr over genes in each chr window (scaled by ref)
        Xall = a.X.toarray() if issparse(a.X) else np.array(a.X)
        # log1p already applied; z-score per gene across cells
        Xz = (Xall - Xall.mean(0)) / (Xall.std(0) + 1e-6)
        cnv_rows = []
        for ch in sorted(chr_genes, key=int):
            gl = [g for g in chr_genes[ch] if g in a.var_names]
            if len(gl) < 20:
                continue
            idx = [list(a.var_names).index(g) for g in gl]
            sub = Xz[:, idx]
            # smooth within chromosome (rolling window 101 genes)
            sm = pd.DataFrame(sub).rolling(101, min_periods=11, center=True).mean().values
            # relative to reference cells
            refv = np.nanmean(sm[ref], axis=0)
            score_c = np.nanmean(sm - refv, axis=1)
            cnv_rows.append((ch, score_c))
        cnv_mat = np.vstack([c for _, c in cnv_rows]).T   # cells x chr
        cnv_amp = cnv_mat.max(axis=1) - cnv_mat.min(axis=1)
        cnv_sum = cnv_mat.sum(axis=1)
        a.obs["cnv_amp"] = cnv_amp
        a.obs["cnv_score"] = cnv_sum
        # comparison: malignant vs non-malignant
        res = pd.DataFrame({"celltype": ct, "cnv_amp": cnv_amp, "cnv_score": cnv_sum})
        res.to_csv(os.path.join(OUT, f"{ds}_cnv.csv"), index=False)
        print(f"{ds} CNV: malignant amp {cnv_amp[mal].mean():.3f} vs non-mal {cnv_amp[~mal].mean():.3f}", flush=True)

    # ---- 4. pseudotime on malignant cells (diffusion/pca along PXS)
    for ds in ["GSE149614", "GSE151530"]:
        a = load_annot(ds)
        if a is None:
            continue
        ct = a.obs["celltype"].astype(str).values
        mal = ct == "Malignant"
        if mal.sum() < 50:
            continue
        am = a[mal].copy()
        sc.pp.normalize_total(am, target_sum=1e4)
        sc.pp.log1p(am)
        sc.pp.highly_variable_genes(am, n_top_genes=1500, flavor="seurat", subset=True)
        sc.pp.scale(am, max_value=10)
        sc.tl.pca(am, svd_solver="arpack")
        # build a 1D pseudotime from PC1..PC3 via diffusion-ish path ordering on PXS
        pxs = am.obs["pxs"].values
        X = am.obsm["X_pca"][:, :3]
        # order by first principal component that correlates with PXS
        corrs = [np.corrcoef(X[:, i], pxs)[0, 1] for i in range(3)]
        axis = int(np.argmax(np.abs(corrs)))
        sgn = np.sign(corrs[axis])
        pseudo = sgn * X[:, axis]
        pseudo = (pseudo - pseudo.min()) / (pseudo.max() - pseudo.min() + 1e-9)
        am.obs["pseudotime"] = pseudo
        # save pseudo + pxs + author celltype for the subtype analysis
        out = pd.DataFrame({"cell": am.obs_names, "pseudotime": pseudo,
                            "pxs": pxs, "celltype": ct[mal],
                            "author": am.obs.get("author_celltype", "").values if "author_celltype" in am.obs else ""})
        out.to_csv(os.path.join(OUT, f"{ds}_malignant_pseudotime.csv"), index=False)
        print(f"{ds} pseudotime done n={mal.sum()}", flush=True)

    # ---- 5. cell-cell communication (curated LR) on GSE149614
    lr = json.load(open(os.path.join(GENESETS, "lr_db.json")))
    a = load_annot("GSE149614")
    if a is not None:
        ct = a.obs["celltype"].astype(str).values
        site = a.obs["site"].astype(str).values if "site" in a.obs.columns else ["tumor"] * a.n_obs
        # cell types to include
        cts = ["Malignant", "Hepatocyte", "T", "NK", "B", "Plasma", "Myeloid", "Mast",
               "Endothelial", "Fibroblast", "Cholangiocyte", "Proliferating"]
        comm_rows = []
        for si in ["Tumor", "Normal", "PVTT", "Lymph"]:
            m = site == si
            if m.sum() < 100:
                continue
            for pair in lr:
                L, R = pair["ligand"], pair["receptor"]
                if L not in a.var_names or R not in a.var_names:
                    continue
                Lx = a[:, L].X.toarray().ravel()[m]
                Rx = a[:, R].X.toarray().ravel()[m]
                ctm = ct[m]
                for sender in cts:
                    for receiver in cts:
                        sl = ctm == sender
                        rr = ctm == receiver
                        if sl.sum() < 5 or rr.sum() < 5:
                            continue
                        # ligand: fraction of sender cells expressing > 0 (z-score > 0)
                        # receptor: fraction of receiver cells expressing > 0
                        ls = float((Lx[sl] > 0).mean())
                        rr_ = float((Rx[rr] > 0).mean())
                        if ls > 0.2 and rr_ > 0.2:
                            comm_rows.append({"site": si, "sender": sender, "receiver": receiver,
                                              "ligand": L, "receptor": R, "pathway": pair["pathway"],
                                              "ligand_expr": float(ls), "receptor_expr": float(rr_)})
        print(f"GSE149614 site={si} candidate pairs: {len(comm_rows)}", flush=True)
        comm = pd.DataFrame(comm_rows)
        comm.to_csv(os.path.join(OUT, "GSE149614_communication.csv"), index=False)
        print(f"GSE149614 communication pairs: {len(comm)}", flush=True)

    # ---- 6. immune cell state scores (GSE149614)
    if a is not None:
        states = {
            "Tex": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "CXCL13"],
            "Treg": ["FOXP3", "IL2RA", "IKZF2", "CTLA4", "BATF"],
            "M1": ["IL1B", "TNF", "CXCL9", "CXCL10", "NOS2"],
            "M2": ["CD163", "MSR1", "MRC1", "TGFBI", "C1QB"],
            "CAF_myofib": ["ACTA2", "TAGLN", "POSTN", "COL1A1"],
            "ECM_CAF": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"],
            "EC_angiogenic": ["ESM1", "ANGPT2", "ACKR1", "VWF"],
            "Neutrophil": ["FCGR3B", "CSF3R", "S100A8", "S100A9"],
            "Cholangiocyte": ["KRT19", "KRT7", "EPCAM"],
        }
        rows = []
        for name, genes in states.items():
            g = [x for x in genes if x in a.var_names]
            if not g:
                continue
            X = a[:, g].X.toarray() if issparse(a[:, g].X) else a[:, g].X
            s = X.mean(1)
            # per celltype
            for c in ["Malignant", "T", "Myeloid", "Endothelial", "Fibroblast", "NK", "B"]:
                m = ct == c
                if m.sum() < 10:
                    continue
                rows.append({"celltype": c, "state": name, "n": int(m.sum()),
                             "mean_score": float(np.nanmean(s[m]))})
        pd.DataFrame(rows).to_csv(os.path.join(OUT, "immune_states_scores.csv"), index=False)
        print("immune states saved", flush=True)

    print("SC ANALYSIS DONE", flush=True)


if __name__ == "__main__":
    main()
