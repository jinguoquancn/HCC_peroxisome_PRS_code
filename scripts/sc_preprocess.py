# -*- coding: utf-8 -*-
"""scRNA-seq pre-processing for the 5 HCC single-cell datasets.
Each dataset -> cached .h5ad (log-normalized counts, PCA, UMAP, leiden, annotation)
Output also a pooled 'atlas' for the major figures (downsampled).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, gzip, json, pickle, time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, dpi_save=150)

RAW = str(PROJECT_ROOT / "data" / "raw")
OUT = str(PROJECT_ROOT / "data" / "processed" / "sc")
os.makedirs(OUT, exist_ok=True)

GENE_POS = None  # gene -> (chr, start) loaded on demand


def load_gene_pos():
    global GENE_POS
    if GENE_POS is None:
        p = os.path.join(OUT, "gene_pos.json")
        if os.path.exists(p):
            GENE_POS = json.load(open(p))
        else:
            GENE_POS = {}
    return GENE_POS


def save_gene_pos(d):
    GENE_POS = d
    json.dump(d, open(os.path.join(OUT, "gene_pos.json"), "w"))


# ------------------------------------------------------------------ loaders
_PROBE = ["ALB", "GAPDH", "CD3D", "LYZ", "VWF", "MKI67"]


def orient_cells_genes(df):
    """Ensure df is genes x cells (rows=genes, cols=cells)."""
    in_idx = sum(g in df.index for g in _PROBE)
    in_col = sum(g in df.columns for g in _PROBE)
    if in_col > in_idx:
        df = df.T
    return df


def _finish(adata):
    # drop non-gene rows/columns (e.g. 'Gene/Cell' header) and dedupe
    bad = [g for g in adata.var_names if "/" in str(g) or str(g).lower() in ("gene", "cell", "nan", "")]
    if bad:
        adata = adata[:, ~adata.var_names.isin(bad)].copy()
    # map Ensembl -> symbol if needed
    if adata.n_vars and str(adata.var_names[0]).startswith("ENSG"):
        from ensg_map import to_symbols
        new, _ = to_symbols(adata.var_names)
        adata.var_names = new
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    return adata


def load_gse149614():
    path = os.path.join(RAW, "GSE149614_HCC_count.txt.gz")
    import gzip as _gz
    rng = np.random.default_rng(7)
    # header line = all cell barcodes
    with _gz.open(path, "rt", errors="ignore") as f:
        header = f.readline().rstrip("\n").split("\t")
    cells = [h for h in header]
    # subsample cells (positions 1..n-1) to keep memory/time tractable;
    # data field j <-> header[j], field 0 = gene name (becomes index)
    pos = np.arange(1, len(cells))
    keep = np.sort(rng.choice(pos, min(len(pos), 15000), replace=False))
    cells = [cells[i] for i in keep]
    usecols = [0] + [int(i) for i in keep]
    chunks = []
    genes_all = []
    with pd.read_csv(path, sep="\t", usecols=usecols, converters={0: str},
                     dtype=np.float32, engine="c", chunksize=2000) as reader:
        for chunk in reader:
            chunk.index.name = None
            genes_all.extend(chunk.index.astype(str).tolist())
            chunks.append(sparse.csr_matrix(chunk.values.astype(np.float32)))
    X = sparse.vstack(chunks).tocsr()   # genes x cells
    adata = ad.AnnData(X.T, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes_all))
    print("  GSE149614 sample n_cells:", len(cells), "n_genes:", len(genes_all), flush=True)
    meta = pd.read_csv(os.path.join(RAW, "GSE149614_metadata.txt.gz"), sep="\t", index_col=0)
    meta.index = meta.index.astype(str)
    common = adata.obs_names.intersection(meta.index)
    adata = adata[common].copy()
    for c in meta.columns:
        adata.obs[c] = meta.loc[common, c].values
    adata.obs["dataset"] = "GSE149614"
    adata.obs["platform"] = "10x_Genomics"
    return _finish(adata)


def load_gse125449():
    mtx = os.path.join(RAW, "GSE125449_Set1_matrix.mtx.gz")
    bc = os.path.join(RAW, "GSE125449_Set1_barcodes.tsv.gz")
    g = os.path.join(RAW, "GSE125449_Set1_genes.tsv.gz")
    sm = os.path.join(RAW, "GSE125449_Set1_samples.txt.gz")
    adata = sc.read_mtx(mtx)   # may be cells x genes or genes x cells
    genes = pd.read_csv(g, sep="\t", header=None, index_col=0)
    barcodes = pd.read_csv(bc, sep="\t", header=None, index_col=0)
    genes.index.name = None
    barcodes.index.name = None
    if adata.n_obs == len(genes) and adata.n_vars == len(barcodes):
        adata = adata.T
    adata.var_names = genes.index.astype(str)
    adata.obs_names = barcodes.index.astype(str)
    samples = pd.read_csv(sm, sep="\t")
    if "Cell Barcode" in samples.columns:
        samples = samples.set_index("Cell Barcode")
    sm_map = samples["Sample"].astype(str).to_dict()
    if "Type" in samples.columns:
        ty_map = samples["Type"].astype(str).to_dict()
    adata.obs["sample"] = [sm_map.get(b, "Set1") for b in adata.obs_names]
    if "Type" in samples.columns:
        adata.obs["Type"] = [ty_map.get(b, "unknown") for b in adata.obs_names]
    adata.obs["dataset"] = "GSE125449"
    adata.obs["platform"] = "Microwell-seq"
    return _finish(adata)


def load_gse146115():
    path = os.path.join(RAW, "GSE146115_count_with_ERCC.txt.gz")
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = orient_cells_genes(df)
    df.index.name = None
    df.columns.name = None
    df = df[~df.index.astype(str).str.startswith("ERCC")]
    X = sparse.csr_matrix(df.values.T.astype(np.float32))
    adata = ad.AnnData(X, obs=pd.DataFrame(index=df.columns.astype(str)))
    adata.var_names = df.index.astype(str)
    adata.obs["dataset"] = "GSE146115"
    adata.obs["platform"] = "Smart-seq2"
    return _finish(adata)


def load_gse151530():
    mtx = os.path.join(RAW, "GSE151530_matrix.mtx.gz")
    bc = os.path.join(RAW, "GSE151530_barcodes.tsv.gz")
    g = os.path.join(RAW, "GSE151530_genes.tsv.gz")
    info = os.path.join(RAW, "GSE151530_Info.txt.gz")
    adata = sc.read_mtx(mtx)
    genes = pd.read_csv(g, sep="\t", header=None, index_col=0)
    barcodes = pd.read_csv(bc, sep="\t", header=None, index_col=0)
    genes.index.name = None
    barcodes.index.name = None
    if adata.n_vars == len(barcodes) and adata.n_obs == len(genes):
        adata = adata.T  # genes x cells mtx
    adata.var_names = genes.index.astype(str)
    adata.obs_names = barcodes.index.astype(str)
    info_df = pd.read_csv(info, sep="\t")
    if "Cell" in info_df.columns:
        info_df = info_df.set_index("Cell")
    info_df.index = info_df.index.astype(str)
    common = adata.obs_names.intersection(info_df.index)
    print(f"  GSE151530 info common: {len(common)}/{adata.n_obs}", flush=True)
    adata = adata[common].copy()
    for c in info_df.columns:
        adata.obs[c] = info_df.loc[common, c].values
    adata.obs["dataset"] = "GSE151530"
    adata.obs["platform"] = "Smart-seq2"
    return _finish(adata)


def load_gse98638():
    # prefer count matrix, fall back to TPM
    cands = [("GSE98638_Tcell_count.txt.gz", "counts"), ("GSE98638_Tcell_TPM.txt.gz", "TPM")]
    for fn, kind in cands:
        path = os.path.join(RAW, fn)
        if not (os.path.exists(path) and os.path.getsize(path) > 0):
            continue
        try:
            df = pd.read_csv(path, sep="\t", index_col=0)
            if "symbol" in df.columns:
                # format: geneID | symbol | cells...
                df = df.set_index("symbol").drop(columns=df.columns[0], errors="ignore")
                df.index.name = None
                df = df.apply(pd.to_numeric, errors="coerce")
            df = orient_cells_genes(df)
            df.index.name = None
            df.columns.name = None
            X = sparse.csr_matrix(df.values.T.astype(np.float32))
            adata = ad.AnnData(X, obs=pd.DataFrame(index=df.columns.astype(str)))
            adata.var_names = df.index.astype(str)
            adata.obs["dataset"] = "GSE98638"
            adata.obs["platform"] = "Smart-seq2"
            adata.obs["data_kind"] = kind
            return _finish(adata)
        except Exception as e:
            print("gse98638 read err", fn, e, flush=True)
    raise FileNotFoundError("GSE98638 data missing")


LOADERS = {"GSE149614": load_gse149614, "GSE125449": load_gse125449,
           "GSE146115": load_gse146115, "GSE151530": load_gse151530,
           "GSE98638": load_gse98638}


# ------------------------------------------------------------------ pipeline
def pipeline(adata, min_genes=200, max_genes=7000, mito_frac=0.25):
    adata.var_names_make_unique()
    adata.obs["n_genes"] = (adata.X > 0).sum(1).A1 if sparse.issparse(adata.X) else (adata.X > 0).sum(1)
    adata.obs["total_counts"] = adata.X.sum(1).A1 if sparse.issparse(adata.X) else adata.X.sum(1)
    mt = adata.var_names.str.startswith("MT-")
    if mt.sum() > 10:
        adata.obs["pct_mt"] = (adata[:, mt].X.sum(1).A1 / np.maximum(adata.obs["total_counts"], 1)) * 100
    else:
        adata.obs["pct_mt"] = 0.0
    keep = (adata.obs["n_genes"] >= min_genes) & (adata.obs["n_genes"] <= max_genes)
    if mt.sum() > 10:
        keep &= (adata.obs["pct_mt"] < mito_frac * 100)
    adata = adata[keep].copy()
    # normalize + log1p
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # hvg
    sc.pp.highly_variable_genes(adata, n_top_genes=2500, flavor="seurat", subset=False)
    # PCA on scaled HVG
    hvg = adata.var_names[adata.var["highly_variable"]].tolist()
    Xh = adata[:, hvg].X.toarray() if sparse.issparse(adata[:, hvg].X) else adata[:, hvg].X
    from sklearn.preprocessing import StandardScaler
    Xs = StandardScaler().fit_transform(Xh)
    adata.obsm["X_pca_in"] = Xs  # store scaled HVG matrix for clustering
    return adata


def main(which):
    for ds in which:
        t0 = time.time()
        cache = os.path.join(OUT, f"{ds}_raw.h5ad")
        if os.path.exists(cache):
            print(f"[skip] {ds} cached", flush=True)
            continue
        adata = LOADERS[ds]()
        print(f"[load] {ds} {adata.shape}", flush=True)
        adata = pipeline(adata)
        print(f"[pipeline] {ds} {adata.shape} in {time.time()-t0:.0f}s", flush=True)
        # simple leiden clusters for reference
        Xs = adata.obsm["X_pca_in"]
        from sklearn.decomposition import PCA
        pca = PCA(n_components=30, random_state=42)
        adata.obsm["X_pca"] = pca.fit_transform(Xs)
        adata.uns["pca_variance"] = pca.explained_variance_ratio_.tolist()
        sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15)
        sc.tl.leiden(adata, resolution=1.0, key_added="leiden")
        sc.tl.umap(adata)
        adata.write(cache)
        print(f"[ok] {ds} cached {time.time()-t0:.0f}s", flush=True)
    print("SC PREPROCESS DONE", flush=True)


if __name__ == "__main__":
    import sys
    which = sys.argv[1:] or list(LOADERS)
    main(which)
