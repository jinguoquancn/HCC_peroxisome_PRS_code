# -*- coding: utf-8 -*-
"""Download all public datasets used in this study (throttled, resumable)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import sys, os, time, threading
sys.path.insert(0, os.path.dirname(__file__))
from dl import par_dl

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/{p}/{g}/{k}/{f}"

# ---------------- bulk cohorts (series matrix) ----------------
BULK = [
    ("GSE10186",  "GSE10nnn",  "matrix"),
    ("GSE60502",  "GSE60nnn",  "matrix"),
    ("GSE16757",  "GSE16nnn",  "matrix"),
    ("GSE39791",  "GSE39nnn",  "matrix"),
    ("GSE54236",  "GSE54nnn",  "matrix"),
    ("GSE116174", "GSE116nnn", "matrix"),
    ("GSE87630",  "GSE87nnn",  "matrix"),
    ("GSE76427",  "GSE76nnn",  "matrix"),
]

# ---------------- single-cell cohorts ----------------
SC = [
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162616/suppl/GSE162616_gene_count.txt.gz",
     "GSE162616_gene_count.txt.gz"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146115/suppl/GSE146115_HCC1-2-5-9_count_with_ERCC.txt.gz",
     "GSE146115_HCC_count.txt.gz"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/GSE125449_Set1_matrix.mtx.gz",
     "GSE125449_Set1_matrix.mtx.gz"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/GSE125449_Set1_barcodes.tsv.gz",
     "GSE125449_Set1_barcodes.tsv.gz"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/GSE125449_Set1_genes.tsv.gz",
     "GSE125449_Set1_genes.tsv.gz"),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98638/suppl/GSE98638_HCC.TCell.S5063.count.txt.gz",
     "GSE98638_Tcell_count.txt.gz"),
]

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "bulk"):
        for g, p, k in BULK:
            f = f"{g}_series_matrix.txt.gz"
            par_dl(GEO.format(p=p, g=g, k=k, f=f), f, nthreads=6)
    if which in ("all", "sc"):
        for u, out in SC:
            par_dl(u, out, nthreads=6)
