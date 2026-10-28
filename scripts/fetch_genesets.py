# -*- coding: utf-8 -*-
"""Fetch MSigDB GMTs, immune signatures (LM22 / MCPcounter-like), and drug data URLs."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, urllib.request, json

BASE = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release"
DST = str(PROJECT_ROOT / "data" / "genesets")
os.makedirs(DST, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}

GMT = {
    "h.all.v2023.2.Hs.symbols.gmt": f"{BASE}/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt",
    "c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt": f"{BASE}/2024.1.Hs/c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt",
    "c2.cp.reactome.v2024.1.Hs.symbols.gmt": f"{BASE}/2024.1.Hs/c2.cp.reactome.v2024.1.Hs.symbols.gmt",
    "c5.go.bp.v2024.1.Hs.symbols.gmt": f"{BASE}/2024.1.Hs/c5.go.bp.v2024.1.Hs.symbols.gmt",
    "c5.go.cc.v2024.1.Hs.symbols.gmt": f"{BASE}/2024.1.Hs/c5.go.cc.v2024.1.Hs.symbols.gmt",
    "c5.go.mf.v2024.1.Hs.symbols.gmt": f"{BASE}/2024.1.Hs/c5.go.mf.v2024.1.Hs.symbols.gmt",
    "c6.all.v2024.1.Hs.symbols.gmt": f"{BASE}/2024.1.Hs/c6.all.v2024.1.Hs.symbols.gmt",
    "c2.cp.wikipathways.v2024.1.Hs.symbols.gmt": f"{BASE}/2024.1.Hs/c2.cp.wikipathways.v2024.1.Hs.symbols.gmt",
}

MISC = {
    "LM22.txt": "https://raw.githubusercontent.com/CIBERSORTx/CIBERSORTx-Users-Guide/main/LM22.txt",
    "MCPcounter_markers.txt": "https://raw.githubusercontent.com/ebecht/MCPcounter/master/Data/MCPcounter_markers.txt",
    "TIDE_signature.txt": "https://raw.githubusercontent.com/jiangyiqun/TIDEpy/master/TIDE/TIDE_signature.txt",
}

CANDIDATE_LM22 = [
    "https://raw.githubusercontent.com/CIBERSORTx/CIBERSORTx-Users-Guide/main/LM22.txt",
    "https://raw.githubusercontent.com/singha53-amrit/.../LM22.txt",
    "https://gist.githubusercontent.com/gdevailly/8a5adbd5b0f6ba1a4b4b3f0a2f3f5b3c/raw/LM22.txt",
]


def get(u, dst, tries=4):
    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        print("[skip]", os.path.basename(dst), flush=True)
        return
    for t in range(tries):
        try:
            req = urllib.request.Request(u, headers=UA)
            with urllib.request.urlopen(req, timeout=300) as r, open(dst, "wb") as f:
                while True:
                    d = r.read(1 << 16)
                    if not d:
                        break
                    f.write(d)
            print("[ok  ]", os.path.basename(dst), os.path.getsize(dst) / 1e6, "MB", flush=True)
            return
        except Exception as e:
            print(f"[retry{t}]", os.path.basename(dst), str(e)[:70], flush=True)
            import time
            time.sleep(3 + 3 * t)


if __name__ == "__main__":
    for name, u in list(GMT.items()) + list(MISC.items()):
        get(u, os.path.join(DST, name))
    print("GENESETS DONE")
