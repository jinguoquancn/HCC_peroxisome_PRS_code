# -*- coding: utf-8 -*-
"""Annotate cell types for each cached single-cell dataset using module scores
(Seurat AddModuleScore-lite via scanpy score_genes) + author labels when present.
Adds per-cell peroxisome module score. Saves annotated h5ad + annotation csv."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sc.settings.verbosity = 0
OUT = str(PROJECT_ROOT / "data" / "processed" / "sc")

MARKERS = {
    "Hepatocyte": ["ALB", "APOA1", "APOA2", "APOC3", "TTR", "FABP1", "SERPINA1", "CYP1A2", "ADH1B"],
    "Malignant": ["AFP", "GPC3", "SPP1", "MDK", "GSTP1"],
    "Cholangiocyte": ["KRT19", "KRT7", "EPCAM", "SOX9", "CFTR"],
    "T": ["CD3D", "CD3E", "CD2", "TRAC", "IL7R", "CD3G"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRC1", "KLRF1", "FGFBP2"],
    "B": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD19"],
    "Plasma": ["MZB1", "IGHG1", "JCHAIN", "SDC1", "XBP1", "TNFRSF17"],
    "Myeloid": ["LYZ", "CD68", "CD14", "CSF1R", "C1QA", "FCGR3A", "ITGAM"],
    "Mast": ["TPSAB1", "CPA3", "KIT", "MS4A2", "TPSB2"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "EGFL7", "ESM1", "RAMP2"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB", "COL3A1"],
    "Proliferating": ["MKI67", "TOP2A", "PCNA", "BIRC5"],
}
MAIN_ORDER = ["Hepatocyte", "Malignant", "Cholangiocyte", "T", "NK", "B", "Plasma",
              "Myeloid", "Mast", "Endothelial", "Fibroblast", "Proliferating"]


def annotate(adata, author_label=None):
    scores = {}
    for ct in MAIN_ORDER:
        g = [x for x in MARKERS[ct] if x in adata.var_names]
        try:
            sc.tl.score_genes(adata, g, score_name=f"score_{ct}", use_raw=False)
        except Exception:
            pass
        if f"score_{ct}" in adata.obs.columns:
            scores[ct] = adata.obs[f"score_{ct}"].values
    S = pd.DataFrame(scores, index=adata.obs_names)
    lab = S.idxmax(axis=1)

    # author labels for quality: map to our schema
    auth_map = {"Hepatocyte": "Hepatocyte", "Malignant cells": "Malignant", "Malignant": "Malignant",
                "T cells": "T", "T": "T", "CD4 T": "T", "CD8 T": "T", "NK": "NK", "B cells": "B",
                "B": "B", "Plasma": "Plasma", "Myeloid": "Myeloid", "Macrophage": "Myeloid",
                "TAMs": "Myeloid", "TAM": "Myeloid", "Mast": "Mast", "TECs": "Endothelial",
                "Endothelial": "Endothelial", "CAFs": "Fibroblast", "Fibroblast": "Fibroblast",
                "CAF": "Fibroblast", "Cholangiocyte": "Cholangiocyte", "Proliferating": "Proliferating"}
    if author_label is not None and author_label in adata.obs.columns:
        auth = adata.obs[author_label].astype(str)
        mapped = auth.map(auth_map).fillna("unknown")
        # use author label where confidently mapped
        conf = mapped != "unknown"
        lab = pd.Series(np.where(conf, mapped, lab), index=adata.obs_names)
        adata.obs["author_celltype"] = auth.values
    adata.obs["celltype"] = lab.values
    return adata


def main():
    prg = json.load(open(str(PROJECT_ROOT / "data" / "processed" / "peroxisome_gene_set.json")))
    prg_genes = prg["genes"]
    core = [x for x in prg.get("modules", {}).get("pex5_import_dependent_core", [])]
    for ds in ["GSE149614", "GSE125449", "GSE146115", "GSE151530", "GSE98638"]:
        cache = os.path.join(OUT, f"{ds}_raw.h5ad")
        if not os.path.exists(cache):
            print("[skip] not processed:", ds, flush=True)
            continue
        out = os.path.join(OUT, f"{ds}_annot.h5ad")
        if os.path.exists(out):
            print("[skip] annotated:", ds, flush=True)
            continue
        adata = sc.read_h5ad(cache)
        auth = "celltype" if ds == "GSE149614" else ("Type" if ds in ("GSE151530", "GSE125449") else None)
        adata = annotate(adata, auth)
        g = [x for x in prg_genes if x in adata.var_names]
        if g:
            sc.tl.score_genes(adata, g, score_name="pxs", use_raw=False)
        cc = [x for x in core if x in adata.var_names]
        if cc:
            sc.tl.score_genes(adata, cc, score_name="pxs_core", use_raw=False)
        adata.write(out)
        print(f"[ok] annotated {ds}: {adata.shape} celltypes=\n"
              f"{adata.obs['celltype'].value_counts().to_string()}", flush=True)
        cols = ["celltype", "pxs", "pxs_core"] + (["author_celltype"] if "author_celltype" in adata.obs else [])
        adata.obs[cols].to_csv(os.path.join(OUT, f"{ds}_annotation.csv"))
    print("SC ANNOTATION DONE", flush=True)


if __name__ == "__main__":
    main()
