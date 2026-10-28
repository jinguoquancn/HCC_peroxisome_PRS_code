# -*- coding: utf-8 -*-
"""Single-threaded downloader for Broad GDAC (no ranged requests) + integrity check."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, sys, time, gzip, urllib.request

RAW = str(PROJECT_ROOT / "data" / "raw")
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
B = "https://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LIHC/20160128/"
JOBS = [
    ("gdac.broadinstitute.org_LIHC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz",
     "TCGA_LIHC_rsem_normalized.tar.gz"),
    ("gdac.broadinstitute.org_LIHC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz",
     "TCGA_LIHC_rsem_raw.tar.gz"),
    ("gdac.broadinstitute.org_LIHC.Merge_Clinical.Level_1.2016012800.0.0.tar.gz",
     "TCGA_LIHC_clinical.tar.gz"),
    ("gdac.broadinstitute.org_LIHC.Mutation_Packager_Calls.Level_3.2016012800.0.0.tar.gz",
     "TCGA_LIHC_mutation_calls.tar.gz"),
]


def ok_gz(p):
    try:
        with gzip.open(p, "rb") as f:
            while f.read(1 << 20):
                pass
        return True
    except Exception:
        return False


def dl(url, out):
    out = os.path.join(RAW, out)
    if os.path.exists(out) and os.path.getsize(out) > 1000 and ok_gz(out):
        print("[skip OK]", os.path.basename(out), flush=True)
        return
    if os.path.exists(out):
        os.remove(out)
    for t in range(6):
        try:
            t0 = time.time()
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=3600) as r, open(out, "wb") as f:
                while True:
                    d = r.read(1 << 18)
                    if not d:
                        break
                    f.write(d)
            dt = max(time.time() - t0, 1)
            good = ok_gz(out)
            print(f"[ok] {os.path.basename(out)} {os.path.getsize(out)/1e6:.1f}MB {dt:.0f}s "
                  f"({os.path.getsize(out)/1e6/dt*1000:.0f} KB/s) integrity={'OK' if good else 'BAD'}", flush=True)
            if good:
                return
        except Exception as e:
            print(f"[retry{t}] {os.path.basename(out)} {str(e)[:70]}", flush=True)
            time.sleep(5 + 5 * t)


if __name__ == "__main__":
    for u, o in JOBS:
        dl(B + u, o)
    print("GDAC DONE", flush=True)
