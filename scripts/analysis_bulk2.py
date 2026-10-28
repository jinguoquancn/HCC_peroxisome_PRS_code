# -*- coding: utf-8 -*-
"""Bulk downstream analyses (Figures 6-8 data):
immune deconvolution, checkpoint profile, mutation landscape/TMB, drug
sensitivity (GDSC oncoPredict-style), GSEA, nomogram, pan-cancer PRS,
published-signature comparison."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, sys, json, warnings
import numpy as np
import pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")

sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
import core

DATA = str(PROJECT_ROOT / "data" / "processed")
GENESETS = str(PROJECT_ROOT / "data" / "genesets")

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

CHECKPOINTS = ["CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "SIGLEC15", "IDO1", "TDO2"]


def load_clin_prs():
    cl = pd.read_csv(os.path.join(DATA, "TCGA_LIHC_clinical_with_PRS.csv"))
    if "Unnamed: 0" in cl.columns:
        cl = cl.set_index("Unnamed: 0")
    cl.index = [str(x).lower() for x in cl.index]
    return cl


def load_expr():
    return core.load_tcga_expr()


def main():
    cl = load_clin_prs()
    expr, info = load_expr()
    # tumor samples aligned by patient
    pat = {s: "-".join(s.split("-")[:3]).lower() for s in info["sample"]}
    tum = info[info["tissue"] == "tumor"]
    E = expr[tum["sample"]]
    E.columns = [pat.get(s, s) for s in tum["sample"]]
    E = E.loc[:, ~E.columns.duplicated(keep="first")]
    E = E[cl.index]
    prs = cl["PRS"]
    # risk group
    med = cl["PRS"].median()
    grp = (cl["PRS"] > med).astype(int)
    cl["risk_group"] = grp

    # ---- 1. immune cell infiltration (ssGSEA)
    rows = []
    for name, gs in IMMUNE_SETS.items():
        g = [x for x in gs if x in E.index]
        if not g:
            continue
        sc = core.ssgsea_score(E, g)
        for gk, lab in [(0, "low"), (1, "high")]:
            v = sc[grp == gk]
            rows.append({"set": name, "group": lab, "mean": v.mean(), "median": v.median()})
        p = stats.mannwhitneyu(sc[grp == 0], sc[grp == 1]).pvalue
        rows.append({"set": name, "group": "p", "mean": p})
    imm = pd.DataFrame(rows)
    imm.to_csv(os.path.join(DATA, "immune_infiltration_PRS.csv"), index=False)
    print("immune infiltration done", flush=True)

    # ---- 2. checkpoint expression
    ck = pd.DataFrame({"gene": [], "logFC": [], "p": []})
    chk = []
    for g in CHECKPOINTS:
        if g not in E.index:
            continue
        hi, lo = E.loc[g, grp == 1], E.loc[g, grp == 0]
        p = stats.mannwhitneyu(lo, hi).pvalue
        lfc = hi.mean() - lo.mean()
        chk.append({"gene": g, "logFC": lfc, "p": p,
                    "mean_high": hi.mean(), "mean_low": lo.mean()})
    pd.DataFrame(chk).to_csv(os.path.join(DATA, "checkpoints_PRS.csv"), index=False)
    print("checkpoints done", flush=True)

    # ---- 3. mutation landscape (TCGA mutation calling)
    maf = os.path.join(DATA, "TCGA_LIHC_mutations.csv")
    if os.path.exists(maf):
        mut = pd.read_csv(maf, index_col=0)  # samples x genes binary ('x')
        mut.index = [str(x).lower() for x in mut.index]
        res = []
        for gene in ["TP53", "CTNNB1", "ARID1A", "TERT", "RB1", "AXIN1", "PTEN", "NFE2L2"]:
            if gene not in mut.columns:
                continue
            for gk, lab in [(0, "low"), (1, "high")]:
                pts = [str(p).lower() for p in cl.index[grp == gk]]
                n = int((mut.loc[[p for p in pts if p in mut.index], gene] == "x").sum())
                res.append({"gene": gene, "group": lab, "n_mut": n, "n_total": len(pts)})
        # TMB
        tmb = pd.read_csv(os.path.join(DATA, "TCGA_LIHC_TMB.csv"), index_col=0)
        tmb.index = [str(x).lower() for x in tmb.index]
        tmb = tmb.iloc[:, 0]
        for gk, lab in [(0, "low"), (1, "high")]:
            pts = [str(p).lower() for p in cl.index[grp == gk]]
            v = tmb[[p for p in pts if p in tmb.index]]
            res.append({"gene": "TMB", "group": lab, "n_mut": float(v.mean()) if len(v) else np.nan,
                        "n_total": len(v)})
        pd.DataFrame(res).to_csv(os.path.join(DATA, "mutation_PRS.csv"), index=False)
        print("mutation landscape done", flush=True)

    # ---- 4. drug sensitivity (GDSC ridge regression per drug)
    gdsc_e = os.path.join(str(PROJECT_ROOT / "data" / "raw"), "GDSC_expression.txt")
    gdsc_i = os.path.join(str(PROJECT_ROOT / "data" / "raw"), "GDSC_IC50.txt")
    if os.path.exists(gdsc_e) and os.path.exists(gdsc_i):
        G = pd.read_csv(gdsc_e, sep="\t", index_col=0)
        IC = pd.read_csv(gdsc_i, sep="\t", index_col=0)
        common = G.columns.intersection(IC.index)
        G, IC = G[common], IC.loc[common]
        # drugs of interest
        drugs = [c for c in IC.columns if any(k in c.lower() for k in
                 ["sorafenib", "lenvatinib", "regorafenib", "cabozantinib", "5-fluorouracil",
                  "gemcitabine", "oxaliplatin", "cisplatin", "doxorubicin", "etoposide",
                  "paclitaxel", "docetaxel", "everolimus", "temsirolimus", "erlotinib", "lapatinib"])]
        genes_common = [g for g in G.index if g in E.index]
        # subset to PRS genes + variance genes
        prs_genes = pd.read_csv(os.path.join(DATA, "TCGA_PRS_coefficients.csv"), index_col=0).index.tolist()
        sel_genes = list(dict.fromkeys(prs_genes + [g for g in genes_common if g in E.index][:4000]))
        Ge = G.loc[[g for g in sel_genes if g in G.index]]
        res = []
        for dr in drugs:
            ic = IC[dr].astype(float)
            d = pd.DataFrame({"y": ic})
            d = d[~d["y"].isna()]
            if len(d) < 30:
                continue
            y = d["y"].values
            Xg = Ge[d.index].T.values
            Xg = (Xg - Xg.mean(0)) / (Xg.std(0) + 1e-9)
            Xg = np.nan_to_num(Xg)
            from sklearn.linear_model import Ridge
            rid = Ridge(alpha=1.0)
            rid.fit(Xg, y)
            # predict IC50 for TCGA tumor samples
            Xt = E.loc[[g for g in Ge.index if g in E.index]].T.values
            Xt = (Xt - Xt.mean(0)) / (Xt.std(0) + 1e-9)
            Xt = np.nan_to_num(Xt)
            pred = rid.predict(Xt)
            for gk, lab in [(0, "low"), (1, "high")]:
                m = grp.values == gk
                res.append({"drug": dr, "group": lab,
                            "mean_pred_IC50": pred[m].mean(),
                            "median_pred_IC50": np.median(pred[m])})
        pd.DataFrame(res).to_csv(os.path.join(DATA, "drug_sensitivity_PRS.csv"), index=False)
        print("drug sensitivity done", flush=True)

    # ---- 5. GSEA-lite (high vs low PRS) — top 200 hallmark-ish pathway genesets from MSigDB hallmark subset
    gsea_path = os.path.join(GENESETS, "hallmark_curated.json")
    if os.path.exists(gsea_path):
        hs = json.load(open(gsea_path))
        res = []
        lfc = np.log2(E.loc[:, grp == 1].mean(1) + 1) - np.log2(E.loc[:, grp == 0].mean(1) + 1)
        ranks = lfc.sort_values(ascending=False)
        allg = ranks.index.tolist()
        rank_idx = {g: i for i, g in enumerate(allg)}
        for name, gs in list(hs.items())[:50]:
            g = [x for x in gs if x in rank_idx]
            if len(g) < 15:
                continue
            pos = np.mean([rank_idx[x] for x in g])
            p = stats.mannwhitneyu(lfc[g].values, lfc[~lfc.index.isin(g)].values).pvalue
            res.append({"pathway": name, "n_genes": len(g), "mean_rank": pos,
                        "mean_logFC": lfc[g].mean(), "p": p})
        gsea = pd.DataFrame(res)
        gsea["FDR"] = stats.false_discovery_control(gsea["p"].values)
        gsea.to_csv(os.path.join(DATA, "GSEA_PRS_high_vs_low.csv"), index=False)
        print("GSEA done", flush=True)

    # ---- 6. nomogram cox data (age, stage, PRS)
    nom = cl[["os_years", "os_event", "age", "PRS"]].copy()
    if "TNM" in cl.columns:
        nom["stage"] = cl["TNM"]
    nom.to_csv(os.path.join(DATA, "nomogram_data.csv"), index=False)
    print("nomogram data saved", flush=True)

    # ---- 7. pan-cancer PRS: compute PRS across TCGA pan-cancer (expression matrices optional)
    # ---- 8. published signature comparison
    sig_dir = os.path.join(GENESETS, "published_signatures.json")
    if os.path.exists(sig_dir):
        sigs = json.load(open(sig_dir))
        res = []
        for sname, genes in sigs.items():
            g = [x for x in genes if x in E.index]
            if len(g) < 3:
                continue
            s = core.risk_score(E, pd.Series(1.0, index=g))  # equal-weight score
            c, _ = core.cindex(cl["os_years"], cl["os_event"], s.values)
            c = max(c, 1 - c)   # direction-agnostic discrimination
            res.append({"signature": sname, "n_genes": len(g), "C_index": c})
        ours, _ = core.cindex(cl["os_years"], cl["os_event"], cl["PRS"].values)
        ours = max(ours, 1 - ours)
        res.append({"signature": "PRS (this study)", "n_genes": len(pd.read_csv(os.path.join(DATA, "TCGA_PRS_coefficients.csv"))),
                    "C_index": ours})
        pd.DataFrame(res).sort_values("C_index", ascending=False).to_csv(
            os.path.join(DATA, "signature_comparison.csv"), index=False)
        print("signature comparison done", flush=True)

    print("BULK DOWNSTREAM DONE", flush=True)


if __name__ == "__main__":
    main()
