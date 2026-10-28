# -*- coding: utf-8 -*-
"""Build GSE149614 h5ad from the preloaded npz + metadata + pipeline."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, sys, time
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

RAW = str(PROJECT_ROOT / "data" / "raw")
SC = str(PROJECT_ROOT / "data" / "processed" / "sc")
os.makedirs(SC, exist_ok=True)

t0 = time.time()
X = sp.load_npz(os.path.join(SC, "GSE149614_counts.npz"))
cells = np.load(os.path.join(SC, "GSE149614_cells.npy"), allow_pickle=True)
genes = np.load(os.path.join(SC, "GSE149614_genes.npy"), allow_pickle=True)
print(f"loaded {X.shape} in {time.time()-t0:.1f}s", flush=True)
adata = sc.AnnData(X.T, obs=pd.DataFrame(index=[str(x) for x in cells]),
                   var=pd.DataFrame(index=[str(x) for x in genes]))
adata.obs_names_make_unique()
adata.var_names_make_unique()
# metadata
meta = pd.read_csv(os.path.join(RAW, "GSE149614_metadata.txt.gz"), sep="\t")
mc = meta.columns.tolist()
print("metadata cols:", mc[:8], "n_rows:", len(meta), flush=True)
key = mc[0]
meta = meta.set_index(key)
meta.index = meta.index.astype(str)
common = adata.obs_names.intersection(meta.index)
print(f"common cells {len(common)}; adata before {adata.shape}", flush=True)
# subset adata first
adata = adata[common].copy()
# build a meta DataFrame indexed by common
meta_sub = meta.loc[common]
print(f"meta_sub {meta_sub.shape}; adata after subset {adata.shape}", flush=True)
for c in mc[1:]:
    adata.obs[c] = meta_sub[c].values
adata.obs["dataset"] = "GSE149614"
adata.obs["platform"] = "10x_Genomics"
print(f"after obs assign {adata.shape}", flush=True)
# ENSG->symbol if needed
if adata.var_names[0].startswith("ENSG"):
    from ensg_map import to_symbols
    new, _ = to_symbols(adata.var_names)
    adata.var_names = new
adata.obs_names_make_unique()
adata.var_names_make_unique()
# pipeline
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=7000)
sc.pp.filter_genes(adata, min_cells=3)
mt = adata.var_names.str.startswith("MT-")
if mt.sum() > 0:
    adata.obs["pct_mt"] = (adata[:, mt].X.sum(1) / adata.X.sum(1)).A.ravel()
    adata = adata[adata.obs["pct_mt"] < 15].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2500, flavor="seurat", subset=False)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata, resolution=1.0, key_added="leiden")
sc.tl.umap(adata)
print(f"pipeline done {adata.shape}", flush=True)
adata.write(os.path.join(SC, "GSE149614_raw.h5ad"))
print("DONE", flush=True)
