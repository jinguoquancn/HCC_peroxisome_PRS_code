# -*- coding: utf-8 -*-
"""Download every external dataset used in the HCC-peroxisome study."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, sys, time, threading, urllib.request

RAW = str(PROJECT_ROOT / "data" / "raw")
os.makedirs(RAW, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}

# url, outname, nthreads
JOBS = [
    # ---------------- TCGA-LIHC (Broad GDAC Firehose, level 3 / level 1) ----------------
    ("https://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LIHC/20160128/"
     "gdac.broadinstitute.org_LIHC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz",
     "TCGA_LIHC_rsem_normalized.tar.gz", 6),
    ("https://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LIHC/20160128/"
     "gdac.broadinstitute.org_LIHC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz",
     "TCGA_LIHC_rsem_raw.tar.gz", 6),
    ("https://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LIHC/20160128/"
     "gdac.broadinstitute.org_LIHC.Merge_Clinical.Level_1.2016012800.0.0.tar.gz",
     "TCGA_LIHC_clinical.tar.gz", 4),
    ("https://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LIHC/20160128/"
     "gdac.broadinstitute.org_LIHC.Mutation_Packager_Calls.Level_3.2016012800.0.0.tar.gz",
     "TCGA_LIHC_mutation_calls.tar.gz", 4),
    # ---------------- GEO bulk validation cohorts ----------------
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/matrix/GSE14520-GPL3921_series_matrix.txt.gz",
     "GSE14520-GPL3921_series_matrix.txt.gz", 4),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/matrix/GSE14520-GPL571_series_matrix.txt.gz",
     "GSE14520-GPL571_series_matrix.txt.gz", 2),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/suppl/GSE14520_Extra_Supplement.txt.gz",
     "GSE14520_Extra_Supplement.txt.gz", 2),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76427/matrix/GSE76427_series_matrix.txt.gz",
     "GSE76427_series_matrix.txt.gz", 6),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116174/matrix/GSE116174_series_matrix.txt.gz",
     "GSE116174_series_matrix.txt.gz", 6),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116174/suppl/GSE116174_HCC-64-u133_plus_2_clinical_data.xls.gz",
     "GSE116174_clinical.xls.gz", 2),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87630/matrix/GSE87630_series_matrix.txt.gz",
     "GSE87630_series_matrix.txt.gz", 6),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10186/matrix/GSE10186_series_matrix.txt.gz",
     "GSE10186_series_matrix.txt.gz", 2),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16757/matrix/GSE16757_series_matrix.txt.gz",
     "GSE16757_series_matrix.txt.gz", 2),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39791/matrix/GSE39791_series_matrix.txt.gz",
     "GSE39791_series_matrix.txt.gz", 4),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54236/matrix/GSE54236_series_matrix.txt.gz",
     "GSE54236_series_matrix.txt.gz", 6),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60502/matrix/GSE60502_series_matrix.txt.gz",
     "GSE60502_series_matrix.txt.gz", 2),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE36nnn/GSE36376/matrix/GSE36376_series_matrix.txt.gz",
     "GSE36376_series_matrix.txt.gz", 8),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25097/matrix/GSE25097_series_matrix.txt.gz",
     "GSE25097_series_matrix.txt.gz", 8),
    # ---------------- single-cell cohorts ----------------
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz",
     "GSE149614_HCC_count.txt.gz", 8),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.metadata.updated.txt.gz",
     "GSE149614_metadata.txt.gz", 2),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/GSE125449_Set1_matrix.mtx.gz",
     "GSE125449_Set1_matrix.mtx.gz", 4),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/GSE125449_Set1_barcodes.tsv.gz",
     "GSE125449_Set1_barcodes.tsv.gz", 1),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/GSE125449_Set1_genes.tsv.gz",
     "GSE125449_Set1_genes.tsv.gz", 1),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125449/suppl/GSE125449_Set1_samples.txt.gz",
     "GSE125449_Set1_samples.txt.gz", 1),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/GSE151530_matrix.mtx.gz",
     "GSE151530_matrix.mtx.gz", 8),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/GSE151530_barcodes.tsv.gz",
     "GSE151530_barcodes.tsv.gz", 1),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/GSE151530_genes.tsv.gz",
     "GSE151530_genes.tsv.gz", 1),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/GSE151530_Info.txt.gz",
     "GSE151530_Info.txt.gz", 1),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146115/suppl/GSE146115_HCC1-2-5-9_count_with_ERCC.txt.gz",
     "GSE146115_count_with_ERCC.txt.gz", 4),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98638/suppl/GSE98638_HCC.TCell.S5063.count.txt.gz",
     "GSE98638_Tcell_count.txt.gz", 4),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98638/suppl/GSE98638_HCC.TCell.S5063.TPM.txt.gz",
     "GSE98638_Tcell_TPM.txt.gz", 4),
    ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162616/suppl/GSE162616_gene_count.txt.gz",
     "GSE162616_gene_count.txt.gz", 4),
]


def _head(url):
    try:
        req = urllib.request.Request(url, headers=UA, method="HEAD")
        with urllib.request.urlopen(req, timeout=90) as r:
            return int(r.headers["Content-Length"]), r.headers.get("Accept-Ranges")
    except Exception:
        return None, None


def par_dl(url, out, nthreads=4, tries=8):
    out = os.path.join(RAW, out)
    if os.path.exists(out) and os.path.getsize(out) > 0:
        print(f"[skip] {os.path.basename(out)} {os.path.getsize(out)/1e6:.2f} MB", flush=True)
        return
    size, acc = _head(url)
    if size is None:
        # HEAD unsupported -> single threaded plain GET
        for t in range(tries):
            try:
                req = urllib.request.Request(url, headers=UA)
                with urllib.request.urlopen(req, timeout=3600) as r, open(out, "wb") as f:
                    while True:
                        d = r.read(1 << 18)
                        if not d:
                            break
                        f.write(d)
                print(f"[ok  ] {os.path.basename(out)} (plain) {os.path.getsize(out)/1e6:.2f} MB", flush=True)
                return
            except Exception as e:
                print(f"[retry{ t}] {os.path.basename(out)} {e}", flush=True)
                time.sleep(4 + 4 * t)
        return
    if acc != "bytes":
        nthreads = 1
    print(f"[get ] {os.path.basename(out)} {size/1e6:.1f} MB x{nthreads}", flush=True)
    chunk = (size + nthreads - 1) // nthreads
    bounds = [(i * chunk, min(size - 1, (i + 1) * chunk - 1)) for i in range(nthreads)]

    def worker(idx):
        a, b = bounds[idx]
        part = f"{out}.part{idx}"
        have = os.path.getsize(part) if os.path.exists(part) else 0
        for t in range(tries):
            try:
                hdr = dict(UA)
                mode = "wb"
                if have and have < (b - a + 1):
                    hdr["Range"] = f"bytes={a+have}-{b}"
                    mode = "ab"
                else:
                    hdr["Range"] = f"bytes={a}-{b}"
                    have = 0
                req = urllib.request.Request(url, headers=hdr)
                with urllib.request.urlopen(req, timeout=1800) as r, open(part, mode) as f:
                    while True:
                        d = r.read(1 << 16)
                        if not d:
                            break
                        f.write(d)
                if os.path.getsize(part) >= (b - a + 1):
                    return True
                have = os.path.getsize(part)
            except Exception as e:
                time.sleep(3 + 4 * t)
        return False

    t0 = time.time()
    ths = [threading.Thread(target=worker, args=(i,)) for i in range(nthreads)]
    for t in ths:
        t.start()
    for t in ths:
        t.join()
    with open(out, "wb") as f:
        for i in range(nthreads):
            p = f"{out}.part{i}"
            if not os.path.exists(p):
                raise RuntimeError("missing part " + p)
            with open(p, "rb") as g:
                f.write(g.read())
            os.remove(p)
    dt = max(time.time() - t0, 1)
    print(f"[ok  ] {os.path.basename(out)} {os.path.getsize(out)/1e6:.2f} MB in {dt:.0f}s "
          f"({os.path.getsize(out)/1e6/dt*1000:.0f} KB/s)", flush=True)


if __name__ == "__main__":
    order = sys.argv[1] if len(sys.argv) > 1 else "all"
    jobs = JOBS
    if order == "tcga":
        jobs = [j for j in JOBS if j[1].startswith("TCGA")]
    elif order == "bulk":
        jobs = [j for j in JOBS if j[1].startswith("GSE") and "series_matrix" in j[1]]
    elif order == "sc":
        jobs = [j for j in JOBS if not j[1].startswith("TCGA") and "series_matrix" not in j[1]]
    for url, out, nt in jobs:
        try:
            par_dl(url, out, nt)
        except Exception as e:
            print(f"[FAIL] {out}: {e}", flush=True)
    print("ALL DONE", flush=True)
