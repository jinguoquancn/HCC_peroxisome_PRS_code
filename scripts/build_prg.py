# -*- coding: utf-8 -*-
"""Build the peroxisome-related gene set (PRG) used throughout the study.
Combines: UniProt GO annotations, MSigDB (GOCC/KEGG/Reactome peroxisome terms),
and literature-curated peroxisomal genes. Saves curated gene list + functional modules.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, re, json, gzip, urllib.request, time

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "genesets")
OUT = os.path.join(HERE, "..", "data", "processed")
os.makedirs(OUT, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}


def uniprot_genes(go_id, size=600):
    url = (f"https://rest.uniprot.org/uniprotkb/search?query=go:{go_id}%20AND%20"
           f"organism_id:9606%20AND%20reviewed:true&fields=gene_primary&format=tsv&size={size}")
    for t in range(4):
        try:
            r = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120)
            txt = r.read().decode("utf8", "ignore")
            return [l.strip() for l in txt.splitlines()[1:] if l.strip()]
        except Exception:
            time.sleep(2 + 2 * t)
    return []


def parse_gmt(path):
    """return dict term -> [genes]"""
    out = {}
    with open(path, encoding="utf8", errors="ignore") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                out[p[0]] = p[2:]
    return out


def main():
    prg = {}   # gene -> source (set of tags)
    def add(genes, tag):
        for g in genes:
            g = g.upper().strip()
            if g and re.match(r"^[A-Z0-9]+$", g):
                prg.setdefault(g, set()).add(tag)

    # 1) UniProt GO
    for go, tag in [("0005777", "GOCC:peroxisome"), ("0005782", "GOCC:peroxisomal_matrix"),
                    ("0005778", "GOCC:peroxisomal_membrane"), ("0005779", "GOCC:peroxisomal_membrane_legacy"),
                    ("0007031", "GOBP:peroxisome_organization"), ("0016558", "GOBP:protein_import_into_peroxisome_matrix"),
                    ("0042765", "GOBP:GPCR"), ("0006635", "GOBP:fatty_acid_beta_oxidation"),
                    ("0008610", "GOBP:lipid_biosynthetic"), ("0019370", "GOBP:leukotriene")]:
        genes = uniprot_genes(go)
        print(f"UniProt GO:{go} -> {len(genes)} genes", flush=True)
        add(genes, tag)

    # 2) MSigDB GMTs
    cc = parse_gmt(os.path.join(DATA, "c5.go.cc.v2024.1.Hs.symbols.gmt"))
    for term, genes in cc.items():
        t = term.split("%")[0]
        if "PEROXISOM" in t.upper():
            add(genes, "MSigDB:" + t)
    for lib, fname in [("kegg", "c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt"),
                       ("reactome", "c2.cp.reactome.v2024.1.Hs.symbols.gmt")]:
        gmt = parse_gmt(os.path.join(DATA, fname))
        for term, genes in gmt.items():
            t = term.split("%")[0]
            if "PEROXISOM" in t.upper():
                add(genes, f"MSigDB:{t}")

    # 3) literature-curated peroxisomal core (PEX machinery + PTS1/PTS2 cargo + beta-ox + ether lipid + H2O2)
    curated = {
        "PEX machinery": ["PEX1","PEX2","PEX3","PEX5","PEX5L","PEX6","PEX7","PEX10","PEX11A","PEX11B","PEX11G",
                          "PEX12","PEX13","PEX14","PEX16","PEX19","PEX26","PEX28","PEX29","PEX30","PEX31","PEX32",
                          "TPR1","TYSND1","ABCD1","ABCD2","ABCD3","ABCD4","SLC25A17","PXMP2","PXMP4","PEX10","FAR1","FAR2",
                          "TTC1","AIMP1","UBE2D3","PEXD1"],
        "beta-oxidation": ["ACOX1","ACOX2","ACOX3","HSD17B4","SCP2","ACAA1","EHHADH","DECR2","ACBD5","CROT","CRAT",
                           "PECR","ECH1","ECI2","HMGCL","HMGCLL1","ACOT2","ACOT4","ACOT6","ACOT8","ACAA1","MFP2","GPD1"],
        "ether lipid": ["GNPAT","AGPS","FAR1","FAR2","TMEM189","TKT","DOHH","SLC25A17","TMEM14C"],
        "ROS/H2O2": ["CAT","PRDX5","SOD1","SOD2","GPX1","GSTK1","GSTZ1","HAO1","HAO2","DAO","DDO","PIPOX","PAOX",
                     "AOC1","FMO1","MAOA","MAOB","SCOX","XDH"],
        "bile acid/misc": ["AMACR","PHYH","AGXT","IDI1","IDI2","FDPS","FDFT1","NSDHL","SC5D","DHCR24","MVD","MVK",
                           "PMVK","HMGCS1","HMGCR","GRHPR","DEPP1","HSDL2","MTARC2","MFF","GDAP1","MAPKBP1","ACOT1",
                           "ACAD11","DECR1","VWA8","ECH1","NUDT7","NUDT19","NUDT12","CLN3","NPC1","APOA1","APOE"],
    }
    for tag, genes in curated.items():
        add(genes, "curated:" + tag)

    # 4) save
    genes_sorted = sorted(prg)
    print(f"\nTOTAL PRG = {len(genes_sorted)}", flush=True)
    out_json = {
        "genes": genes_sorted,
        "source": {g: sorted(list(s)) for g, s in prg.items()},
    }
    with open(os.path.join(OUT, "peroxisome_gene_set.json"), "w", encoding="utf8") as f:
        json.dump(out_json, f, indent=1)

    # functional submodule assignment (priority order)
    modules = {"peroxisomal_import_machinery": [], "fatty_acid_beta_oxidation": [],
               "ether_lipid_synthesis": [], "ros_detoxification": [], "bile_acid_and_others": []}
    def in_(gene, key):
        return any(k == key or k.startswith(key) for k in prg.get(gene, set()))

    for g in genes_sorted:
        if "curated:PEX machinery" in prg[g] or "GOBP:protein_import_into_peroxisome_matrix" in prg[g]:
            modules["peroxisomal_import_machinery"].append(g)
        elif "curated:beta-oxidation" in prg[g] or "GOBP:fatty_acid_beta_oxidation" in prg[g]:
            modules["fatty_acid_beta_oxidation"].append(g)
        elif "curated:ether lipid" in prg[g]:
            modules["ether_lipid_synthesis"].append(g)
        elif "curated:ROS/H2O2" in prg[g]:
            modules["ros_detoxification"].append(g)
        else:
            modules["bile_acid_and_others"].append(g)

    # PTS1/PTS2-import-dependent core set (PEX5 cargo) for method innovation #3
    pex5_dependent = ["CAT","ACOX1","ACOX2","ACOX3","HSD17B4","SCP2","ACAA1","EHHADH","DECR2","AGPS","GNPAT",
                      "PRDX5","PHYH","AMACR","DAO","DDO","PIPOX","PAOX","HAO1","HAO2","CROT","CRAT","PECR",
                      "ECH1","ECI2","HMGCL","HMGCLL1","ACOT8","IDI1","FDPS","GPD1","PEX5","PEX7","PEX14",
                      "SLC25A17","ABCD1","ABCD2","ABCD3","TYSND1"]
    pex5_dependent = [g for g in pex5_dependent if g in genes_sorted]
    modules["pex5_import_dependent_core"] = pex5_dependent
    out_json["modules"] = {k: sorted(set(v)) for k, v in modules.items()}
    with open(os.path.join(OUT, "peroxisome_gene_set.json"), "w", encoding="utf8") as f:
        json.dump(out_json, f, indent=1)
    for k, v in modules.items():
        print(f"  module {k}: {len(set(v))} genes", flush=True)
    print("PRG SET DONE", flush=True)


if __name__ == "__main__":
    main()
