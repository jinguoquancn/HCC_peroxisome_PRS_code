# -*- coding: utf-8 -*-
"""Build per-cohort gene-expression matrices + survival tables for the GEO
external validation cohorts. Requires gpl_mapping.pkl (probe->symbol)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, re, gzip, json, pickle, sys
import numpy as np
import pandas as pd
import openpyxl

sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
RAW = str(PROJECT_ROOT / "data" / "raw")
OUT = str(PROJECT_ROOT / "data" / "processed")
ANN = str(PROJECT_ROOT / "data" / "genesets")


def load_mapping():
    p = os.path.join(ANN, "gpl_mapping.pkl")
    if os.path.exists(p):
        return pickle.load(open(p, "rb"))
    return {}


def expr_to_genes(d, mapping, log2=False):
    """d: parsed series matrix dict. mapping: probe->symbol."""
    probes = d["probes"]
    X = d["X"]
    m = []
    for i, p in enumerate(probes):
        pk = str(p).replace("_PM_", "_")   # Affy PM-only probe naming
        g = mapping.get(pk)
        if g:
            g = str(g).split("///")[0].strip()   # multiple symbols -> first
            if g:
                m.append((pk, g, X[i]))
    gene_list = sorted(set(g for _, g, _ in m))
    out = np.full((len(gene_list), X.shape[1]), np.nan)
    for gix, g in enumerate(gene_list):
        rows = [x for _, gg, x in m if gg == g]
        with np.errstate(invalid="ignore"):
            means = np.nanmean(np.vstack(rows), axis=1)   # per-probe mean
        if not np.any(~np.isnan(means)):
            out[gix] = rows[0]
            continue
        best = int(np.nanargmax(means))
        out[gix] = rows[best]
    df = pd.DataFrame(out, index=gene_list, columns=d["samples"])
    if log2:
        df = np.log2(df + 1)
    return df


def parse_gse14520_clinical():
    """GSE14520 Extra_Supplement clinical (OS/RFS)."""
    f = os.path.join(RAW, "GSE14520_Extra_Supplement.txt.gz")
    df = pd.read_csv(f, sep="\t", compression="gzip", dtype=str)
    out = pd.DataFrame({
        "Affy_GSM": df["Affy_GSM"].astype(str),
        "tissue": df["Tissue Type"],
        "gender": df["Gender"],
        "age": pd.to_numeric(df["Age"], errors="coerce"),
        "TNM": df["TNM staging"],
        "BCLC": df["BCLC staging"],
        "AFP": df.get("AFP (>/<=300ng/ml)"),
        "os_status": pd.to_numeric(df["Survival status"], errors="coerce"),
        "os_months": pd.to_numeric(df["Survival months"], errors="coerce"),
        "rfs_status": pd.to_numeric(df["Recurr status"], errors="coerce"),
        "rfs_months": pd.to_numeric(df["Recurr months"], errors="coerce"),
    })
    return out


def parse_gse116174_clinical():
    f = os.path.join(RAW, "GSE116174_clinical.xls")
    if not os.path.exists(f):
        import gzip as _gz
        raw = _gz.open(f + ".gz", "rb").read()
        open(f, "wb").write(raw)
    xls = pd.ExcelFile(f)
    df = xls.parse("mRNA")
    df = df.rename(columns={"Follow_up time(month)": "os_months",
                            "Event_death": "os_event"})
    return df


def main():
    mapping = load_mapping()
    print("GPL mappings loaded:", {k: len(v) for k, v in mapping.items()}, flush=True)

    raw = pickle.load(open(os.path.join(OUT, "geo_series_raw.pkl"), "rb"))
    new = pickle.load(open(os.path.join(OUT, "geo_new.pkl"), "rb")) if os.path.exists(os.path.join(OUT, "geo_new.pkl")) else {}
    raw.update(new)
    # rename cohort keys with GPL suffix
    for k in list(raw):
        if k.startswith("GSE14520"):
            pass

    cohorts = {
        "GSE14520": "GPL3921",
        "GSE76427": "GPL10558",
        "GSE116174": "GPL570",
        "GSE10186": "GPL5474",
        "GSE36376": "GPL10558",
        "GSE25097": "GPL10687",
        "GSE87630": "GPL6947",
        "GSE54236": "GPL6480",
        "GSE60502": "GPL96",
        "GSE16757": "GPL6102",
        "GSE39791": "GPL10558",
    }
    res = {}
    for cohort, gpl in cohorts.items():
        key = cohort if cohort in raw else (cohort + "-" + gpl)
        if key not in raw:
            print("missing", cohort, flush=True)
            continue
        d = raw[key]
        mp = mapping.get(gpl, {})
        if not mp:
            print(f"  no mapping for {gpl}, skipping {cohort}", flush=True)
            continue
        expr = expr_to_genes(d, mp, log2=False)
        print(f"== {cohort}: expr {expr.shape}", flush=True)
        res[cohort] = {"expr": expr, "meta": d["per_sample"]}
        expr.to_csv(os.path.join(OUT, f"{cohort}_expression.csv"))

    # survival tables
    surv = {}
    # GSE14520
    if os.path.exists(os.path.join(RAW, "GSE14520_Extra_Supplement.txt.gz")):
        cl = parse_gse14520_clinical()
        cl.to_csv(os.path.join(OUT, "GSE14520_clinical.csv"), index=False)
        surv["GSE14520"] = cl
        print("GSE14520 clinical n=", len(cl), flush=True)
    # GSE76427 (from characteristics)
    if "GSE76427" in res:
        ps = res["GSE76427"]["meta"]
        rows = []
        for s, v in ps.items():
            c = v.get("char", {})
            def fnum(x):
                try:
                    return float(x)
                except Exception:
                    return np.nan
            if c.get("tissue", "") == "primary hepatocellular carcinoma tumor":
                rows.append({"sample": s, "os_time": fnum(c.get("duryears_os", np.nan)),
                             "os_event": fnum(c.get("event_os", np.nan)),
                             "rfs_time": fnum(c.get("duryears_rfs", np.nan)),
                             "rfs_event": fnum(c.get("event_rfs", np.nan)),
                             "age": c.get("age (years)", ""), "gender": c.get("gender (1=m, 2=f)", ""),
                             "bclc": c.get("bclc_staging", ""), "tnm": c.get("tnm_staging_clinical", "")})
        sdf = pd.DataFrame(rows)
        sdf["os_months"] = sdf["os_time"] * 12
        sdf["rfs_months"] = sdf["rfs_time"] * 12
        sdf.to_csv(os.path.join(OUT, "GSE76427_clinical.csv"), index=False)
        surv["GSE76427"] = sdf
        print("GSE76427 clinical n=", len(sdf), flush=True)
    # GSE10186
    if "GSE10186" in res:
        ps = res["GSE10186"]["meta"]
        rows = []
        for s, v in ps.items():
            c = v.get("char", {})
            st = str(c.get("survival status (0", "")).strip()
            try:
                ev = int(st.split(":")[-1].strip())
            except Exception:
                ev = 1 if "dead" in st.lower() else 0
            rows.append({"sample": s, "os_days": float(c.get("survival time (days)", np.nan)),
                         "os_event": ev, "hbsag": c.get("hepatitis b virus (hbs antigen)", ""),
                         "hcv": c.get("hepatitis c virus", ""), "mvi": c.get("microvascular invasion", ""),
                         "satellite": c.get("satellite lesions", "")})
        sdf = pd.DataFrame(rows)
        sdf["os_months"] = sdf["os_days"] / 30.44
        sdf.to_csv(os.path.join(OUT, "GSE10186_clinical.csv"), index=False)
        surv["GSE10186"] = sdf
        print("GSE10186 clinical n=", len(sdf), flush=True)
    # GSE116174
    if os.path.exists(os.path.join(RAW, "GSE116174_clinical.xls")):
        df = parse_gse116174_clinical()
        # map sampleID -> GSM via series matrix source_name/title
        d = raw["GSE116174"]
        name2gsm = {}
        for k, v in d["per_sample"].items():
            for key in ("source_name", "title"):
                if v.get(key):
                    name2gsm[str(v[key]).strip()] = k
        df["GSM"] = df["sampleID"].astype(str).map(name2gsm)
        df["sid"] = df["GSM"]
        df.to_csv(os.path.join(OUT, "GSE116174_clinical.csv"), index=False)
        surv["GSE116174"] = df
        print("GSE116174 clinical n=", len(df), "cols:", df.columns.tolist(), flush=True)

    with open(os.path.join(OUT, "geo_cohorts.pkl"), "wb") as f:
        pickle.dump(res, f)
    with open(os.path.join(OUT, "geo_survival.pkl"), "wb") as f:
        pickle.dump(surv, f)
    print("GEO BUILD DONE", flush=True)


if __name__ == "__main__":
    main()
