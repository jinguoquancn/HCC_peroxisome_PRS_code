# -*- coding: utf-8 -*-
"""Process TCGA-LIHC: clinical/survival table, RSEM expression matrix (tumor+normal),
mutation matrix (for TMB & oncoprint). Saves to data/processed/."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, re, gzip, tarfile, pickle, json
import numpy as np
import pandas as pd

RAW = str(PROJECT_ROOT / "data" / "raw")
OUT = str(PROJECT_ROOT / "data" / "processed")
EX = os.path.join(RAW, "extracted_tcga")
os.makedirs(OUT, exist_ok=True)


# ---------------------------------------------------------------- clinical
def parse_clinical():
    base = EX
    fdir = [d for d in os.listdir(base) if "Merge_Clinical" in d][0]
    f = os.path.join(base, fdir, "LIHC.clin.merged.txt")
    df = pd.read_csv(f, sep="\t", header=None, index_col=0)
    # GDAC merged files lack a header; sample ids live in patient.bcr_patient_barcode row
    br = df.index[df.index.astype(str).str.contains("bcr_patient_barcode", na=False)]
    if len(br):
        df.columns = df.loc[br[0]].tolist()
        df = df.drop(index=br[0])
    samples = [str(c) for c in df.columns]
    row = lambda k: df.loc[k].to_dict() if k in df.index else {}
    vd = row("patient.days_to_death")
    vlf = row("patient.days_to_last_followup")
    vs = row("patient.vital_status")
    age = row("patient.age_at_initial_pathologic_diagnosis")
    gender = row("patient.gender")
    race = row("patient.race")
    cpg = row("patient.child_pugh_classification_grade")
    fib = row("patient.fibrosis_ishak_score")

    def num(x):
        try:
            return float(x)
        except Exception:
            return np.nan

    recs = {}
    for s in samples:
        t = s.split("-")
        patient = "-".join(t[:3])
        t_d = num(vd.get(s, np.nan))
        t_l = num(vlf.get(s, np.nan))
        # follow-up records may extend last followup
        fu_days = []
        for c in df.index:
            if re.match(rf"patient\.follow_ups\.follow_up-\d+\.days_to_last_followup", c):
                v = num(df.loc[c].get(s, np.nan))
                if not np.isnan(v):
                    fu_days.append(v)
            if re.match(rf"patient\.follow_ups\.follow_up-\d+\.days_to_death", c):
                v = num(df.loc[c].get(s, np.nan))
                if not np.isnan(v) and np.isnan(t_d):
                    t_d = v
        if not np.isnan(t_l):
            fu_days.append(t_l)
        fu = max(fu_days) if fu_days else np.nan
        alive = str(vs.get(s, "")).lower()
        if alive in ("dead", "deceased"):
            os_time, os_ev = t_d, 1
        elif alive in ("alive", "living"):
            os_time, os_ev = fu, 0
        else:
            os_time, os_ev = (t_d, 1) if not np.isnan(t_d) else (fu, 0) if not np.isnan(fu) else (np.nan, np.nan)
        recs[patient] = {
            "patient": patient, "sample": s, "tissue": "tumor" if len(t) > 3 and t[3][:2] in ("01", "04", "06", "07", "08", "09") else "normal",
            "os_time_days": os_time, "os_event": os_ev,
            "age": num(age.get(s, np.nan)), "gender": gender.get(s, ""), "race": race.get(s, ""),
            "child_pugh": cpg.get(s, ""), "fibrosis": fib.get(s, ""),
        }
    clin = pd.DataFrame(list(recs.values())).set_index("patient")
    # stage / grade from the "merged_only_clinical_clin_format" (has TNM)
    f2 = os.path.join(base, fdir, "LIHC.merged_only_clinical_clin_format.txt")
    if os.path.exists(f2):
        d2 = pd.read_csv(f2, sep="\t", index_col=0)
        for k, col in [("pathologic_stage", "patient.stage_event.tnm_categories.pathologic_categories.pathologic_stage"),
                       ("tumor_grade", "patient.tumor_grade"),
                       ("tumor_status", "patient.tumor_status"),
                       ("days_to_death", "patient.days_to_death"),
                       ("days_to_last_followup", "patient.days_to_last_followup")]:
            if col in d2.index:
                vals = d2.loc[col].to_dict()
                for s, v in vals.items():
                    p = "-".join(s.split("-")[:3])
                    if p in clin.index:
                        clin.loc[p, k] = v
    clin["os_years"] = clin["os_time_days"] / 365.25
    clin.to_csv(os.path.join(OUT, "TCGA_LIHC_clinical.csv"))
    print("clinical:", clin.shape, flush=True)
    return clin


# ---------------------------------------------------------------- expression
def parse_rsem(tag="normalized"):
    fdir = [d for d in os.listdir(EX) if "Merge_rnaseqv2" in d and "normalized" in d][0]
    f = os.path.join(EX, fdir)
    # list actual data file
    files = [x for x in os.listdir(f) if x.endswith(".data.txt")]
    if not files:
        # tar may not be extracted yet
        return None
    data_f = os.path.join(f, files[0])
    # GDAC format: first line has "gene" then sample ids; may have "gene_id" rows
    df = pd.read_csv(data_f, sep="\t", index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    # drop the "gene" placeholder row(s)
    df = df[~df.index.isin(["gene", "scaled_estimate"])]
    print("rsem", tag, df.shape, flush=True)
    return df


# ---------------------------------------------------------------- mutation
def parse_mutation():
    # Mutation_Packager_Calls tar
    tars = [x for x in os.listdir(RAW) if x.startswith("TCGA_LIHC_mutation_calls") and x.endswith(".tar.gz")]
    mut = None
    for ta in tars:
        p = os.path.join(RAW, ta)
        try:
            with tarfile.open(p, "r:gz") as tf:
                members = [m for m in tf.getmembers() if m.name.endswith(".maf.txt")]
                if not members:
                    members = [m for m in tf.getmembers() if m.name.endswith(".txt")]
                parts = []
                for m in members[:400]:
                    fh = tf.extractfile(m)
                    if fh is None:
                        continue
                    df = pd.read_csv(fh, sep="\t", comment="#", nrows=200000, low_memory=False)
                    cols = [c.lower() for c in df.columns]
                    if "hugo_symbol" not in cols:
                        continue
                    # sample from filename: .../TCGA-BC-A10W-01.maf.txt -> patient id (lowercase)
                    samp = os.path.basename(m.name).split(".")[0]
                    base = "-".join(samp.split("-")[:3]).lower()
                    df = df.copy()
                    df.columns = cols
                    df["tumor_sample"] = base
                    parts.append(df)
                if not parts:
                    continue
                allc = set()
                for df in parts:
                    allc.update(df.columns)
                allc = sorted(allc)
                mut = pd.concat([df.reindex(columns=allc) for df in parts], ignore_index=True)
                print("mutation maf rows:", mut.shape, flush=True)
                break
        except Exception as e:
            print("mut err", e, flush=True)
        if mut is not None:
            break
    if mut is None:
        return None
    # build sample-level binary mutation matrix for top genes + TMB (nonsynonymous)
    nonsyn = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
              "In_Frame_Del", "In_Frame_Ins", "Nonstop_Mutation", "Translation_Start_Site",
              "Splice_Site"}
    m = mut[mut["variant_classification"].isin(nonsyn)] if "variant_classification" in mut.columns else mut
    m = m.drop_duplicates(["tumor_sample", "hugo_symbol"])
    sample_genes = m.pivot_table(index="tumor_sample", columns="hugo_symbol", values="variant_classification",
                                 aggfunc=lambda x: "x" if len(x) else None)
    sample_genes = sample_genes.fillna("")
    sample_genes = sample_genes.replace("x", "x")  # keep marker 'x', NaN->""
    # TMB = # non-syn mutations / 30 Mb (WXS ~30Mb)
    tmb = m.groupby("tumor_sample").size() / 30.0
    print("mutation matrix:", sample_genes.shape, flush=True)
    return sample_genes, tmb


if __name__ == "__main__":
    clin = parse_clinical()
    print(clin[["os_time_days", "os_event", "age", "gender"]].head(), flush=True)
    print("n unique patients:", clin.shape[0], flush=True)
    # expression: extract tars if present
    for tag, tar_name in [("normalized", "TCGA_LIHC_rsem_normalized.tar.gz"),
                          ("raw", "TCGA_LIHC_rsem_raw.tar.gz")]:
        p = os.path.join(RAW, tar_name)
        if os.path.exists(p) and not any("Merge_rnaseqv2" in d for d in os.listdir(EX)):
            try:
                with tarfile.open(p) as tf:
                    tf.extractall(EX)
                print("extracted", tar_name, flush=True)
            except Exception as e:
                print("extract fail", tar_name, e, flush=True)
    expr = parse_rsem()
    if expr is not None:
        expr.to_csv(os.path.join(OUT, "TCGA_LIHC_rsem_normalized_matrix.csv"))
    mut = parse_mutation()
    if mut is not None:
        mg, tmb = mut
        mg.to_csv(os.path.join(OUT, "TCGA_LIHC_mutation_matrix.csv"))
        tmb.to_frame("TMB").to_csv(os.path.join(OUT, "TCGA_LIHC_TMB.csv"))
    print("TCGA PROCESSING DONE", flush=True)
