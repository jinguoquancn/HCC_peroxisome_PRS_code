# -*- coding: utf-8 -*-
"""Rebuild Supplementary_Tables.xlsx so that the sheet names and contents
match the S1-S12 Table list declared in the PLOS ONE manuscript
(manuscript/plos_one/final_manuscript.md -> ## Supporting information).

Outputs:
  results/supplementary/Supplementary_Tables_PLOSONE.xlsx
  manuscript/plos_one/Supplementary_Tables.xlsx (copy for submission)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment

ROOT = str(PROJECT_ROOT)
DATA = os.path.join(ROOT, "data/processed")
SC = os.path.join(DATA, "sc")
GENESETS = os.path.join(ROOT, "data/genesets")
OUT_MAIN = os.path.join(ROOT, "results/supplementary/Supplementary_Tables_PLOSONE.xlsx")
OUT_COPY = os.path.join(ROOT, "manuscript/plos_one/Supplementary_Tables.xlsx")
os.makedirs(os.path.dirname(OUT_MAIN), exist_ok=True)
os.makedirs(os.path.dirname(OUT_COPY), exist_ok=True)

IMMUNE_SETS = {
    "CD8_T": ["CD8A", "CD8B", "GZMB", "PRF1", "IFNG", "GZMA"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "T_exhausted": ["PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX"],
    "NK": ["NKG7", "GNLY", "KLRD1", "NCR1"],
    "B": ["MS4A1", "CD79A", "CD19", "CD79B"],
    "Macrophage_M1": ["IL1B", "TNF", "CXCL9", "CXCL10"],
    "Macrophage_M2": ["CD163", "MSR1", "MRC1", "TGFBI"],
    "Neutrophil": ["FCGR3B", "CSF3R", "S100A8", "S100A9"],
    "DC": ["CD1C", "CLEC9A", "ITGAX", "BATF3"],
    "CAF": ["COL1A1", "DCN", "LUM", "PDGFRB"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "KDR"],
    "PD1_axis": ["PDCD1", "CD274", "PDCD1LG2", "CTLA4"],
}

# ---- S1 Table
def t_prs_genes():
    co = pd.read_csv(os.path.join(DATA, "TCGA_PRS_coefficients.csv"))
    co.columns = ["Gene", "Coefficient"]
    co["HR"] = np.exp(co["Coefficient"])
    co["Direction"] = np.where(co["Coefficient"] > 0, "Risk", "Protective")
    co = co.sort_values("Coefficient", ascending=False).reset_index(drop=True)
    co.insert(0, "Rank", co.index + 1)
    return co

# ---- S2 Table
def t_univariate_cox():
    d = pd.read_csv(os.path.join(DATA, "TCGA_univariate_Cox_PRG.csv"))
    d = d.sort_values("p").reset_index(drop=True)
    d["Rank"] = d.index + 1
    return d[["Rank", "gene", "coef", "low", "high", "HR", "p", "FDR"]]

# ---- S3 Table
def t_wgcna():
    d = pd.read_csv(os.path.join(DATA, "TCGA_WGCNA_modules.csv"))
    return d

# ---- S4 Table
def t_external_validation():
    d = pd.read_csv(os.path.join(DATA, "external_validation_summary.csv"))
    d = d[d["HR"].notna()].copy()
    d["Significant"] = d["p_logrank"] < 0.05
    return d[["cohort", "n", "HR", "HR_low", "HR_high", "p_logrank",
              "C_index", "AUC_1y", "AUC_3y", "AUC_5y", "Significant"]]

# ---- S5 Table
def t_immune():
    d = pd.read_csv(os.path.join(DATA, "immune_infiltration_PRS.csv"))
    mean_pv = d[d["group"] != "p"].pivot(index="set", columns="group", values="mean")
    med_pv = d[d["group"] != "p"].pivot(index="set", columns="group", values="median")
    pvals = d[d["group"] == "p"].set_index("set")["mean"]
    out = pd.DataFrame({
        "ssGSEA_set": mean_pv.index,
        "mean_PRS_low": mean_pv["low"].values,
        "mean_PRS_high": mean_pv["high"].values,
        "median_PRS_low": med_pv["low"].values,
        "median_PRS_high": med_pv["high"].values,
        "Mann_Whitney_p": [pvals.get(s, np.nan) for s in mean_pv.index],
    })
    out["Direction_PRS_high_minus_low"] = out["mean_PRS_high"] - out["mean_PRS_low"]
    # BH-FDR across the 12-set family (revision: reviewer #1.2)
    fdr = pd.read_csv(os.path.join(DATA, "revision_fdr.csv"))
    fdr = fdr[fdr["family"] == "immune_ssGSEA"].set_index("feature")["BH_FDR"]
    out["BH_FDR_q"] = [fdr.get(s, np.nan) for s in out["ssGSEA_set"]]
    out["Significant_p<0.05"] = out["Mann_Whitney_p"] < 0.05
    out["Significant_BH_FDR_q<0.05"] = out["BH_FDR_q"] < 0.05
    return out

# ---- S6 Table
def t_mutations():
    d = pd.read_csv(os.path.join(DATA, "mutation_PRS.csv"))
    rows = []
    for g in d["gene"].unique():
        if g == "TMB":
            continue
        lo = d[(d["gene"] == g) & (d["group"] == "low")]
        hi = d[(d["gene"] == g) & (d["group"] == "high")]
        if not len(lo) or not len(hi):
            continue
        n_mut_lo, n_tot_lo = int(lo["n_mut"].iloc[0]), int(lo["n_total"].iloc[0])
        n_mut_hi, n_tot_hi = int(hi["n_mut"].iloc[0]), int(hi["n_total"].iloc[0])
        table = [[n_mut_hi, n_tot_hi - n_mut_hi],
                 [n_mut_lo, n_tot_lo - n_mut_lo]]
        _, p = stats.fisher_exact(table)
        rows.append({
            "Gene": g, "Mut_PRS_low": n_mut_lo, "Total_low": n_tot_lo,
            "Mut_PRS_high": n_mut_hi, "Total_high": n_tot_hi,
            "Frac_PRS_low": n_mut_lo / n_tot_lo,
            "Frac_PRS_high": n_mut_hi / n_tot_hi,
            "Fisher_exact_p": p,
            "Significant_p<0.05": p < 0.05,
        })
    tmb = d[d["gene"] == "TMB"]
    if len(tmb):
        rows.append({
            "Gene": "TMB (mean)", "Mut_PRS_low": float(tmb[tmb["group"] == "low"]["n_mut"].iloc[0]),
            "Total_low": int(tmb[tmb["group"] == "low"]["n_total"].iloc[0]),
            "Mut_PRS_high": float(tmb[tmb["group"] == "high"]["n_mut"].iloc[0]),
            "Total_high": int(tmb[tmb["group"] == "high"]["n_total"].iloc[0]),
            "Frac_PRS_low": np.nan, "Frac_PRS_high": np.nan,
            "Fisher_exact_p": np.nan, "Significant_p<0.05": False,
        })
    outm = pd.DataFrame(rows)
    fdr = pd.read_csv(os.path.join(DATA, "revision_fdr.csv"))
    fdr = fdr[fdr["family"] == "mutation_fisher"].set_index("feature")["BH_FDR"]
    outm["BH_FDR_q"] = [fdr.get(g, np.nan) for g in outm["Gene"]]
    outm["Significant_BH_FDR_q<0.05"] = outm["BH_FDR_q"] < 0.05
    return outm

# ---- S7 Table
def t_gsea():
    return pd.read_csv(os.path.join(DATA, "GSEA_PRS_high_vs_low.csv"))

# ---- S8 Table
def t_signature_comparison():
    d = pd.read_csv(os.path.join(DATA, "signature_comparison.csv"))
    d = d.sort_values("C_index", ascending=False).reset_index(drop=True)
    d["Rank"] = d.index + 1
    return d

# ---- S9 Table
def t_immune_panel():
    rows = []
    for k, v in IMMUNE_SETS.items():
        rows.append({"ssGSEA_panel": k, "n_genes": len(v),
                     "Genes": ", ".join(v), "Source": "Curated from prior literature"})
    return pd.DataFrame(rows)

# ---- S10 Table
def t_lr_pairs():
    # curated database
    cur = pd.DataFrame(LR_PAIRS, columns=["Ligand", "Receptor", "Pathway"])
    cur["Source"] = "build_lr_db.py curated"
    # observed in GSE149614 communication
    obs = pd.read_csv(os.path.join(SC, "GSE149614_communication.csv"))
    obs_pairs = (obs[["pathway", "ligand", "receptor"]]
                 .drop_duplicates()
                 .rename(columns={"pathway": "Pathway",
                                  "ligand": "Ligand",
                                  "receptor": "Receptor"}))
    obs_pairs["Source"] = "GSE149614 (observed)"
    cur2 = cur[["Pathway", "Ligand", "Receptor", "Source"]]
    out = pd.concat([cur2, obs_pairs], ignore_index=True)
    out = out.drop_duplicates(["Pathway", "Ligand", "Receptor"]).reset_index(drop=True)
    out.insert(0, "Index", out.index + 1)
    return out

# ---- S11 Table
def t_drug_sensitivity():
    d = pd.read_csv(os.path.join(DATA, "drug_sensitivity_PRS.csv"))
    pv = d.pivot(index="drug", columns="group", values="mean_pred_IC50")
    pv = pv.rename(columns={"low": "mean_pred_IC50_PRS_low",
                            "high": "mean_pred_IC50_PRS_high"})
    med = d.pivot(index="drug", columns="group", values="median_pred_IC50")
    med = med.rename(columns={"low": "median_pred_IC50_PRS_low",
                              "high": "median_pred_IC50_PRS_high"})
    out = pd.concat([pv, med], axis=1).reset_index()
    out["Drug"] = out["drug"].str.replace(r"_\d+$", "", regex=True)
    out["Note"] = "Hypothesis-generating (GDSC ridge-regression)"
    return out

# ---- S12 Table
def t_meta_analysis():
    d = pd.read_csv(os.path.join(DATA, "external_validation_summary.csv"))
    d = d[d["HR"].notna()].copy()
    y = np.log(d["HR"].values)
    se = (np.log(d["HR_high"].values) - np.log(d["HR_low"].values)) / (2 * 1.96)
    w = 1 / se ** 2
    pooled_fe = np.sum(w * y) / np.sum(w)
    Q = np.sum(w * (y - pooled_fe) ** 2)
    k = len(y)
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0, (Q - (k - 1)) / c)
    w_re = 1 / (se ** 2 + tau2)
    pooled_re = np.sum(w_re * y) / np.sum(w_re)
    se_re = np.sqrt(1 / np.sum(w_re))
    z = pooled_re / se_re
    p_pooled = 2 * (1 - stats.norm.cdf(abs(z)))
    HR = np.exp(pooled_re)
    CI_low = np.exp(pooled_re - 1.96 * se_re)
    CI_high = np.exp(pooled_re + 1.96 * se_re)

    inp = d[["cohort", "n", "HR", "HR_low", "HR_high", "p_logrank"]].copy()
    inp.columns = ["Cohort", "n", "HR", "HR_95CI_low", "HR_95CI_high", "LogRank_p"]
    inp["SE_logHR"] = np.round(se, 4)
    inp["Weight_fixed"] = np.round(w, 3)
    inp["Weight_random"] = np.round(w_re, 3)

    # DerSimonian-Laird random-effects meta-analysis of the five cohort HRs.
    # The pooled estimate is computed from the per-cohort log HRs and their
    # standard errors derived from the reported 95% CIs; Cochran's Q, I2 and
    # tau2 are reported as heterogeneity measures (revision: reviewer #2.4).
    summary = pd.DataFrame({
        "Statistic": [
            "Pooled_HR (DerSimonian-Laird random effects)",
            "Pooled_95CI_low",
            "Pooled_95CI_high",
            "Pooled_p",
            "k_cohorts",
            "tau2",
            "Q_Cochran",
            "I2_percent",
        ],
        "Value": [
            round(HR, 3), round(CI_low, 3), round(CI_high, 3),
            f"{p_pooled:.3f}",
            k, round(tau2, 5), round(Q, 3),
            round(max(0, (Q - (k - 1)) / Q) * 100, 2) if Q > 0 else 0.0,
        ],
    })
    out = pd.concat([summary,
                     pd.DataFrame({"Statistic": ["", "--- input cohorts ---"],
                                   "Value": ["", ""]}),
                     inp.reset_index(drop=True)], ignore_index=True)
    return out


# need LR_PAIRS from build_lr_db (imported)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import build_lr_db  # noqa: E402
LR_PAIRS = build_lr_db.LR_PAIRS


# ---- S13 Table (nested CV + multivariable Cox + proportional hazards)
def t_nested_cv_multivariable():
    rows = []
    nc = pd.read_csv(os.path.join(DATA, "revision_nestedCV_cindex.csv"))
    for _, r in nc.iterrows():
        ci = r["nested_cv_cindex"]
        rows.append(("Nested cross-validation (10-fold, leakage-free)", r["model"],
                     "out-of-fold C-index", None if pd.isna(ci) else round(float(ci), 4)))
        rows.append(("Nested cross-validation (10-fold, leakage-free)", r["model"],
                     "folds used", int(r["n_folds"])))
    mx = pd.read_csv(os.path.join(DATA, "revision_multivariable_cox.csv"))
    for _, r in mx.iterrows():
        rows.append(("Multivariable Cox (complete case)", r["term"], "HR (95% CI)",
                     "%.2f (%.2f-%.2f)" % (r["HR"], r["CI_low"], r["CI_high"])))
        rows.append(("Multivariable Cox (complete case)", r["term"], "p", float(r["p"])))
    meta = json.load(open(os.path.join(DATA, "revision_multivariable_summary.json")))
    for k, lab in [("n_complete_case", "n complete-case"),
                   ("n_events", "events"),
                   ("n_child_pugh_missing", "Child-Pugh missing"),
                   ("c_index_clinical", "C-index clinical only"),
                   ("c_index_clinical_plus_PRS", "C-index clinical + PRS"),
                   ("LRT_add_PRS_stat", "LRT chi-square for adding PRS"),
                   ("LRT_add_PRS_p", "LRT p")]:
        rows.append(("Clinical baseline vs baseline + PRS", lab, "", meta[k]))
    phf = os.path.join(DATA, "revision_ph_test.csv")
    if os.path.exists(phf):
        ph = pd.read_csv(phf)
        for _, r in ph.iterrows():
            rows.append(("Proportional hazards (Schoenfeld, %s time)" % r["time_transform"],
                         r["term"], "p", float(r["p"])))
    return pd.DataFrame(rows, columns=["Block", "Item", "Detail", "Value"])


# ---- S14 Table (single-cell dataset effects + patient-level tests)
def t_sc_patient_level():
    rows = []
    b = pd.read_csv(os.path.join(DATA, "revision_sc_batch.csv"))
    for _, r in b.iterrows():
        rows.append(("Cell level (cells as observations)", r["dataset"],
                     "n malignant / n hepatocyte",
                     "%d / %d" % (r["n_malignant"], r["n_hepatocyte"])))
        rows.append(("Cell level (cells as observations)", r["dataset"],
                     "mean PXS malignant / hepatocyte",
                     "%.4f / %.4f" % (r["mean_malignant"], r["mean_hepatocyte"])))
        rows.append(("Cell level (cells as observations)", r["dataset"], "Cohen's d",
                     None if pd.isna(r["cohens_d"]) else round(float(r["cohens_d"]), 3)))
        rows.append(("Cell level (cells as observations)", r["dataset"],
                     "Mann-Whitney p", None if pd.isna(r["mwu_p"]) else float(r["mwu_p"])))
        rows.append(("Cell level (cells as observations)", r["dataset"],
                     "direction", r["direction"]))
    t = pd.read_csv(os.path.join(DATA, "revision_sc_patient_tests.csv"))
    for _, r in t.iterrows():
        rows.append(("Patient level (paired pseudobulk, %s)" % r["score"], r["dataset"],
                     "n patients paired", int(r["n_patients_paired"])))
        rows.append(("Patient level (paired pseudobulk, %s)" % r["score"], r["dataset"],
                     "mean difference (malignant - hepatocyte)",
                     round(float(r["mean_diff_malignant_minus_hepatocyte"]), 4)))
        rows.append(("Patient level (paired pseudobulk, %s)" % r["score"], r["dataset"],
                     "patients lower / higher",
                     "%d / %d" % (r["n_patients_malignant_lower"],
                                  r["n_patients_malignant_higher"])))
        rows.append(("Patient level (paired pseudobulk, %s)" % r["score"], r["dataset"],
                     "Wilcoxon signed-rank p",
                     None if pd.isna(r["wilcoxon_p"]) else float(r["wilcoxon_p"])))
    var = json.load(open(os.path.join(DATA, "revision_sc_variance.json")))
    for k in ["R2_dataset_only", "R2_celltype_only", "R2_dataset_plus_celltype"]:
        rows.append(("Variance decomposition of PXS", k.replace("_", " "), "R2", var[k]))
    pl = var["patient_level_pooled"]
    rows.append(("Pooled patient-level test (all datasets)", "n patients paired", "",
                 pl["n_patients_paired"]))
    rows.append(("Pooled patient-level test (all datasets)", "patients lower / higher", "",
                 "%d / %d" % (pl["n_patients_malignant_lower"],
                              pl["n_patients_malignant_higher"])))
    rows.append(("Pooled patient-level test (all datasets)", "Wilcoxon p", "", pl["wilcoxon_p"]))
    return pd.DataFrame(rows, columns=["Block", "Item", "Detail", "Value"])


# ---- S15 Table (external validation rigour + heterogeneity)
def t_external_rigor():
    d = pd.read_csv(os.path.join(DATA, "revision_external_rigor.csv"))
    meta = json.load(open(os.path.join(DATA, "revision_meta_heterogeneity.json")))
    d["Het_I2_percent"] = np.nan
    d["Het_tau2"] = np.nan
    d["Het_Q"] = np.nan
    pooled = {c: np.nan for c in d.columns}
    pooled.update({
        "cohort": "POOLED (DerSimonian-Laird random effects)",
        "HR_median_split": meta["pooled_HR"],
        "HR_split_CI_low": meta["CI_low"],
        "HR_split_CI_high": meta["CI_high"],
        "logrank_p_median_split": meta["p"],
        "Het_I2_percent": meta["I2_percent"],
        "Het_tau2": meta["tau2"],
        "Het_Q": meta["Q"],
    })
    return pd.concat([d, pd.DataFrame([pooled])], ignore_index=True)


# ---- S16 Table (tumour-purity proxy)
def t_purity_confounding():
    d = pd.read_csv(os.path.join(DATA, "revision_purity_confounding.csv"))
    return d.rename(columns={
        "feature": "Feature",
        "n_genes_immune": "n_immune_markers",
        "n_genes_stromal": "n_stromal_markers",
        "spearman_rho_with_PRS": "Spearman_rho_with_PRS",
        "p_unadjusted": "p_unadjusted",
        "spearman_partial_rho_given_purity": "Spearman_partial_rho_controlling_purity",
        "p_adjusted_for_purity": "p_adjusted_for_purity",
    })


# ---- S6 Table (checkpoint expression + mutation frequency, long format, with BH-FDR)
def t_checkpoint_mutation():
    rows = []
    fdr = pd.read_csv(os.path.join(DATA, "revision_fdr.csv"))
    qcp = fdr[fdr["family"] == "checkpoint"].set_index("feature")["BH_FDR"]
    cp = pd.read_csv(os.path.join(DATA, "checkpoints_PRS.csv"))
    for _, r in cp.iterrows():
        rows.append(("Checkpoint expression (PRS-high vs PRS-low)", r["gene"],
                     "log2FC %+.2f" % r["logFC"], float(r["p"]),
                     float(qcp.get(r["gene"], np.nan))))
    qm = fdr[fdr["family"] == "mutation_fisher"].set_index("feature")["BH_FDR"]
    mut = pd.read_csv(os.path.join(DATA, "mutation_PRS.csv"))
    for g in mut["gene"].unique():
        if g == "TMB":
            continue
        lo = mut[(mut["gene"] == g) & (mut["group"] == "low")]
        hi = mut[(mut["gene"] == g) & (mut["group"] == "high")]
        if not len(lo) or not len(hi):
            continue
        nl, tl = int(lo["n_mut"].iloc[0]), int(lo["n_total"].iloc[0])
        nh, th = int(hi["n_mut"].iloc[0]), int(hi["n_total"].iloc[0])
        p = stats.fisher_exact([[nh, th - nh], [nl, tl - nl]])[1]
        rows.append(("Mutation frequency (Fisher exact)", g,
                     "PRS-low %d/%d vs PRS-high %d/%d" % (nl, tl, nh, th),
                     float(p), float(qm.get(g, np.nan))))
    return pd.DataFrame(rows, columns=["Block", "Feature", "Effect", "p", "BH_FDR_q"])


TABLES = [
    ("S1 Table. PRS gene list and Cox coefficients.", t_prs_genes),
    ("S2 Table. Univariate Cox screening of PRGs with BH-adjusted q-values.", t_univariate_cox),
    ("S3 Table. WGCNA module definitions and associations with PXS and survival.", t_wgcna),
    ("S4 Table. External validation statistics (HR, 95% CI, p, C-index, AUC) per cohort.", t_external_validation),
    ("S5 Table. Immune-infiltration ssGSEA scores by PRS group: mean per group, direction of effect, Mann-Whitney p and Benjamini-Hochberg q across the 12-set family.", t_immune),
    ("S6 Table. Checkpoint expression and mutation frequencies by PRS group with Benjamini-Hochberg q across the 10-gene checkpoint and 8-gene mutation families.", t_checkpoint_mutation),
    ("S7 Table. Hallmark pathway enrichment statistics (BH-FDR).", t_gsea),
    ("S8 Table. Head-to-head signature comparison C-indices.", t_signature_comparison),
    ("S9 Table. Curated 12-signature ssGSEA panel definitions.", t_immune_panel),
    ("S10 Table. Ligand-receptor pairs retained in cell-cell communication analysis.", t_lr_pairs),
    ("S11 Table. GDSC predicted IC50 by risk group.", t_drug_sensitivity),
    ("S12 Table. DerSimonian-Laird random-effects meta-analysis input and output (cohort HRs, SEs, weights, Q, I2, tau2, pooled estimate).", t_meta_analysis),
    ("S13 Table. Nested cross-validation, multivariable Cox regression, likelihood-ratio test and proportional-hazards tests.", t_nested_cv_multivariable),
    ("S14 Table. Dataset-stratified and patient-level single-cell comparison of malignant versus non-malignant hepatocyte peroxisome scores, with contributing patient numbers, per-cohort effect direction and variance decomposition.", t_sc_patient_level),
    ("S15 Table. External validation rigour: per-cohort C-index with 95% bootstrap CI, event counts, median follow-up, and meta-analysis heterogeneity statistics.", t_external_rigor),
    ("S16 Table. Tumour-purity proxy analysis: PRS-purity correlation and Spearman partial correlations of immune and checkpoint features controlling the purity proxy.", t_purity_confounding),
]


def write_xlsx(path):
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    header_fill = PatternFill("solid", fgColor="DCE6F1")
    title_font = Font(bold=True, size=11, color="1F4E79")
    for i, (title, fn) in enumerate(TABLES, start=1):
        ws = wb.create_sheet(f"S{i}_Table")
        df = fn()
        # title row
        c = ws.cell(row=1, column=1, value=title)
        c.font = title_font
        c.fill = header_fill
        ws.merge_cells(start_row=1, start_column=1,
                       end_row=1, end_column=max(1, df.shape[1]))
        # data
        for r, row in enumerate(df.itertuples(index=False), start=2):
            for c_, v in enumerate(row, start=1):
                cell = ws.cell(row=r, column=c_, value=(None if pd.isna(v) else v))
        # freeze + column widths
        ws.freeze_panes = "A3"
        for col in ws.columns:
            col_letter = col[1].column_letter if len(col) > 1 else col[0].column_letter
            mx = max((len(str(c.value)) if c.value is not None else 0) for c in col)
            ws.column_dimensions[col_letter].width = min(40, max(10, mx + 2))
    wb.save(path)
    return [(f"S{i}", t[0], fn().shape) for i, (t, fn) in enumerate(TABLES, start=1)]


if __name__ == "__main__":
    summary = write_xlsx(OUT_MAIN)
    print("wrote", OUT_MAIN)
    for s in summary:
        print(" ", s)
    # copy to manuscript/plos_one
    import shutil
    shutil.copyfile(OUT_MAIN, OUT_COPY)
    print("copied ->", OUT_COPY)
    print("SUPP TABLES DONE")
