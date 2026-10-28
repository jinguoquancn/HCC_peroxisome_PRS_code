# -*- coding: utf-8 -*-
"""Download GPL SOFT annotation files and extract probe->gene symbol mapping."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, gzip, re, urllib.request, time

RAW = str(PROJECT_ROOT / "data" / "raw")
ANN = str(PROJECT_ROOT / "data" / "genesets")
os.makedirs(ANN, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}

GPLS = ["GPL3921", "GPL571", "GPL570", "GPL96", "GPL10558", "GPL6947",
        "GPL6480", "GPL5474", "GPL10687"]


def gpl_dir(g):
    num = g[3:]
    n = len(num)
    if n <= 3:
        return "GPLnnn"
    if n == 4:
        return "GPL" + num[:1] + "nnn"
    return "GPL" + num[:2] + "nnn"


def dl(g):
    d = gpl_dir(g)
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{d}/{g}/soft/{g}_family.soft.gz"
    out = os.path.join(ANN, f"{g}_family.soft.gz")
    if os.path.exists(out) and os.path.getsize(out) > 1000:
        return out
    for t in range(6):
        try:
            t0 = time.time()
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=3600) as r, open(out, "wb") as f:
                while True:
                    b = r.read(1 << 18)
                    if not b:
                        break
                    f.write(b)
            print(f"[ok] {g} {os.path.getsize(out)/1e6:.1f}MB in {time.time()-t0:.0f}s", flush=True)
            return out
        except Exception as e:
            print(f"[retry{t}] {g} {str(e)[:60]}", flush=True)
            time.sleep(4 + 4 * t)
    return None


def extract_mapping(g):
    """Parse SOFT family file platform table -> {probe_id: gene_symbol}"""
    path = os.path.join(ANN, f"{g}_family.soft.gz")
    if not os.path.exists(path):
        return {}
    header = None
    probe_idx = sym_idx = None
    mapping = {}
    in_table = False
    with gzip.open(path, "rt", encoding="utf8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                for i, h in enumerate(header):
                    hl = h.strip().lower().replace(" ", "_")
                    if hl == "id":
                        probe_idx = i
                    if hl in ("gene_symbol", "symbol", "gene", "gene_name", "genesymbol"):
                        sym_idx = i
                if probe_idx is None:
                    return {}
                continue
            if probe_idx is not None and sym_idx is not None and len(parts) > max(probe_idx, sym_idx):
                pid = parts[probe_idx].strip()
                sym = parts[sym_idx].strip()
                if sym and sym not in ("---", "NA") and pid:
                    # keep first non-dup (or join)
                    if pid not in mapping:
                        mapping[pid] = sym
    print(f"    {g}: {len(mapping)} probes mapped", flush=True)
    return mapping


if __name__ == "__main__":
    all_map = {}
    for g in GPLS:
        dl(g)
        all_map[g] = extract_mapping(g)
    import pickle
    with open(os.path.join(ANN, "gpl_mapping.pkl"), "wb") as f:
        pickle.dump(all_map, f)
    print("GPL MAPPING DONE", flush=True)
