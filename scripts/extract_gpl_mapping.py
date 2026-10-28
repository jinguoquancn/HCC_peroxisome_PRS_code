# -*- coding: utf-8 -*-
"""Extract probe->gene symbol mappings from downloaded GPL annot/soft files.
Saves data/genesets/gpl_mapping.pkl"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, gzip, pickle

RAW = str(PROJECT_ROOT / "data" / "raw")
ANN = str(PROJECT_ROOT / "data" / "genesets")
GPLS = ["GPL3921", "GPL571", "GPL570", "GPL96", "GPL10558", "GPL6947",
        "GPL6480", "GPL5474", "GPL10687"]


def extract_mapping(path):
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
                if "platform_table_begin" in line:
                    in_table = True
                    continue
                if "platform_table_end" in line:
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
                        if hl in ("gene_symbol", "symbol", "gene", "gene_name", "genesymbol", "genenamesymbol"):
                            sym_idx = i
                    if probe_idx is None:
                        return {}
                    continue
                if probe_idx is not None and sym_idx is not None and len(parts) > max(probe_idx, sym_idx):
                    pid = parts[probe_idx].strip()
                    sym = parts[sym_idx].strip()
                    if sym and sym not in ("---", "NA", "na") and pid:
                        mapping.setdefault(pid, sym)
    except Exception as e:
        print(f"  extract err {os.path.basename(path)}: {e}", flush=True)
        return {}
    return mapping


def main():
    all_map = {}
    for g in GPLS:
        # prefer annot.gz
        cands = [os.path.join(RAW, f"{g}.annot.gz"),
                 os.path.join(RAW, f"{g}_family.soft.gz")]
        mp = {}
        for p in cands:
            mp = extract_mapping(p)
            if mp:
                print(f"  {g}: {len(mp)} probes via {os.path.basename(p)}", flush=True)
                break
        all_map[g] = mp
    with open(os.path.join(ANN, "gpl_mapping.pkl"), "wb") as f:
        pickle.dump(all_map, f)
    print("GPL MAPPING DONE", flush=True)


if __name__ == "__main__":
    main()
