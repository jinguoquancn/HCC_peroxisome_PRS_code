# -*- coding: utf-8 -*-
"""Minimal loader for GSE149614: read subsampled count matrix -> sparse npz.
Avoids anndata/scanpy imports to reduce memory during the big parse."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import gzip, os, sys, time
import numpy as np
import pandas as pd
import scipy.sparse as sp

RAW = str(PROJECT_ROOT / "data" / "raw")
OUT = str(PROJECT_ROOT / "data" / "processed" / "sc")
os.makedirs(OUT, exist_ok=True)

N_CELLS = int(sys.argv[1]) if len(sys.argv) > 1 else 12000
path = os.path.join(RAW, "GSE149614_HCC_count.txt.gz")
rng = np.random.default_rng(7)
with gzip.open(path, "rt", errors="ignore") as f:
    header = f.readline().rstrip("\n").split("\t")
pos = np.arange(1, len(header))
keep = np.sort(rng.choice(pos, min(len(pos), N_CELLS), replace=False))
cells = [header[i] for i in keep]
usecols = [0] + [int(i) for i in keep]
t0 = time.time()
chunks = []
genes = []
with pd.read_csv(path, sep="\t", usecols=usecols, converters={0: str},
                 dtype=np.float32, engine="c", chunksize=2000) as reader:
        for chunk in reader:
            chunk.index.name = None
            genes.extend(chunk.index.astype(str).tolist())
            # drop the first data column (gene-name string col added by usecols[0])
            arr = chunk.values[:, 1:].astype(np.float32)
            chunks.append(sp.csr_matrix(arr))
        print("chunk", len(genes), "rows in", round(time.time() - t0, 1), "s", flush=True)
X = sp.vstack(chunks).tocsr()
print("read done", X.shape, "in", round(time.time() - t0, 1), "s", flush=True)
sp.save_npz(os.path.join(OUT, "GSE149614_counts.npz"), X)
np.save(os.path.join(OUT, "GSE149614_cells.npy"), np.array(cells, dtype=object), allow_pickle=True)
np.save(os.path.join(OUT, "GSE149614_genes.npy"), np.array(genes, dtype=object), allow_pickle=True)
print("SAVED", flush=True)
