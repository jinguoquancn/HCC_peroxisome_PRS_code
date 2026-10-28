# -*- coding: utf-8 -*-
"""Download all GPL SOFT files (parallel, resume, integrity-checked) and extract
probe->gene-symbol mappings. Saves data/genesets/gpl_mapping.pkl."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, sys, gzip, pickle, time, urllib.request

sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from dl import par_dl

RAW = str(PROJECT_ROOT / "data" / "raw")
ANN = str(PROJECT_ROOT / "data" / "genesets")
os.makedirs(ANN, exist_ok=True)

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


def ok_gz(p):
    try:
        with gzip.open(p, "rb") as f:
            while f.read(1 << 20):
                pass
        return True
    except Exception:
        return False


def download_all():
    for g in GPLS:
        out = os.path.join(RAW, f"{g}_family.soft.gz")
        if os.path.exists(out) and os.path.getsize(out) > 1000:
            if ok_gz(out):
                print(f"[skip OK] {g}", flush=True)
                continue
            print(f"[bad] re-downloading {g}", flush=True)
            os.remove(out)
        d = gpl_dir(g)
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{d}/{g}/soft/{g}_family.soft.gz"
        try:
            par_dl(url, f"{g}_family.soft.gz", nthreads=8, tries=12)
        except Exception as e:
            print(f"[FAIL] {g}: {e}", flush=True)
        if not ok_gz(out):
            print(f"[CORRUPT] {g} even after download", flush=True)
    print("GPL DL DONE", flush=True)


def extract_mapping(g):
    path = os.path.join(RAW, f"{g}_family.soft.gz")
    if not os.path.exists(path):
        return {}
    header = None
    probe_idx = sym_idx = None
    mapping = {}
    in_table = False
    try:
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
                        mapping.setdefault(pid, sym)
    except Exception as e:
        print(f"  extract err {g}: {e}", flush=True)
        return {}
    print(f"    {g}: {len(mapping)} probes mapped", flush=True)
    return mapping


def main():
    download_all()
    all_map = {}
    for g in GPLS:
        all_map[g] = extract_mapping(g)
    with open(os.path.join(ANN, "gpl_mapping.pkl"), "wb") as f:
        pickle.dump(all_map, f)
    print("GPL MAPPING DONE", flush=True)


if __name__ == "__main__":
    main()
