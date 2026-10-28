# -*- coding: utf-8 -*-
"""Parse GEO series-matrix files -> per-cohort gene expression + clinical/survival tables.
Outputs pickled dict + CSV for each cohort in data/processed/."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, re, gzip, pickle, json, sys
import numpy as np

RAW = str(PROJECT_ROOT / "data" / "raw")
OUT = str(PROJECT_ROOT / "data" / "processed")
os.makedirs(OUT, exist_ok=True)

# cohort -> which probes (GPL) matter, and survival fields to keep
COHORTS = [
    "GSE14520-GPL3921", "GSE76427", "GSE116174", "GSE54236",
    "GSE87630", "GSE10186", "GSE36376", "GSE25097", "GSE39791", "GSE60502", "GSE16757",
]


def stripq(s):
    return s.strip().strip('"').strip()


def parse_series_matrix(path):
    """Return dict with keys: samples, characteristics (per sample dict), value matrix rows."""
    with gzip.open(path, "rt", encoding="utf8", errors="ignore") as f:
        lines = f.readlines()

    sample_ids = []
    char_keys = {}   # key -> list per sample
    mat = []         # rows: (id_ref, [values])
    in_table = False

    for ln in lines:
        ln = ln.rstrip("\n")
        if ln.startswith("!Sample_geo_accession"):
            sample_ids = [stripq(x) for x in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics"):
            vals = [stripq(x) for x in ln.split("\t")[1:]]
            char_keys.setdefault("char", [])
            char_keys["char"].append(vals)
        elif ln.startswith("!Sample_source_name"):
            char_keys["source_name"] = [stripq(x) for x in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_title"):
            char_keys["title"] = [stripq(x) for x in ln.split("\t")[1:]]
        elif ln.startswith("!series_matrix_table_begin"):
            in_table = True
        elif in_table:
            if ln.startswith("!series_matrix_table_end"):
                in_table = False
                break
            parts = ln.split("\t")
            probe = stripq(parts[0])
            vals = parts[1:]
            mat.append((probe, vals))

    # assemble per-sample characteristics dictionary
    per_sample = {}
    for i, sid in enumerate(sample_ids):
        per_sample[sid] = {}
        for k, vals in char_keys.items():
            if k == "char":
                # merge all characteristic lines for this sample
                d = {}
                for vl in vals:
                    if i < len(vl):
                        for item in vl[i].split("|"):
                            if ":" in item:
                                kk, vv = item.split(":", 1)
                                d[kk.strip()] = vv.strip()
                per_sample[sid]["char"] = d
            elif i < len(vals):
                per_sample[sid][k] = vals[i]

    # matrix rows -> float (also strip quotes)
    probes = [m[0] for m in mat]
    X = []
    for m in mat:
        row = []
        for v in m[1]:
            vv = stripq(v)
            try:
                row.append(float(vv))
            except ValueError:
                row.append(np.nan)
        X.append(row)
    X = np.array(X, dtype=np.float64)
    return {"samples": sample_ids, "per_sample": per_sample, "probes": probes, "X": X}


def collapse_rows(X, probes, gene_for_probe, method="max"):
    """Map probe-level to gene-level by max-mean."""
    genes = []
    keep = []
    for p in probes:
        g = gene_for_probe(p)
        if g:
            keep.append((p, g))
    gene_list = sorted(set(g for _, g in keep))
    out = np.zeros((len(gene_list), X.shape[1]))
    for gi, g in enumerate(gene_list):
        idxs = [i for i, (p, gg) in enumerate(keep) if gg == g]
        sub = X[idxs, :]
        with np.errstate(invalid="ignore"):
            means = np.nanmean(sub, axis=1)
        best = int(np.nanargmax(means))
        out[gi, :] = sub[best, :]
    return gene_list, out


def main():
    out_all = {}
    for name in COHORTS:
        path = os.path.join(RAW, f"{name}_series_matrix.txt.gz")
        if not os.path.exists(path):
            print("missing", name, flush=True)
            continue
        try:
            d = parse_series_matrix(path)
        except Exception as e:
            print("ERR parse", name, e, flush=True)
            continue
        print(f"== {name}: {len(d['samples'])} samples, {len(d['probes'])} probes, "
              f"{d['X'].shape} mat", flush=True)
        # print characteristics keys for first sample to inspect clinically relevant fields
        s0 = d["samples"][0] if d["samples"] else None
        if s0:
            print("   chars:", json.dumps(d["per_sample"][s0], ensure_ascii=False)[:500], flush=True)
        out_all[name] = d

    with open(os.path.join(OUT, "geo_series_raw.pkl"), "wb") as f:
        pickle.dump(out_all, f)


if __name__ == "__main__":
    main()
