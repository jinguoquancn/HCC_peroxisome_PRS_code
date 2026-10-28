# -*- coding: utf-8 -*-
"""Fetch gene -> (chromosome, start) positions via mygene (genomic_pos) for
CNV-style analysis on single-cell data. Saves data/genesets/gene_pos.json."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, json, time, urllib.request

ANN = str(PROJECT_ROOT / "data" / "genesets")
os.makedirs(ANN, exist_ok=True)


def mygene_pos(genes):
    genes = list(dict.fromkeys(g for g in genes if g and "/" not in g))
    out = {}
    for st in range(0, len(genes), 800):
        chunk = genes[st:st + 800]
        body = json.dumps({"q": chunk, "scopes": "symbol", "fields": ["symbol", "genomic_pos"],
                           "species": "human"}).encode()
        for t in range(6):
            try:
                req = urllib.request.Request("https://mygene.info/v3/query", data=body,
                                             headers={"Content-Type": "application/json",
                                                      "User-Agent": "Mozilla/5.0"})
                res = json.load(urllib.request.urlopen(req, timeout=180))
                for r in res:
                    q = r.get("query")
                    if "notfound" in r or not r.get("genomic_pos"):
                        continue
                    gp = r["genomic_pos"]
                    if isinstance(gp, dict) and gp.get("chr"):
                        out[q] = (gp["chr"].replace("chr", ""), gp.get("start"))
                break
            except Exception as e:
                print("retry", st, e, flush=True)
                time.sleep(4 + 4 * t)
        print("pos", st + 800, flush=True)
    return out


def main():
    # union of genes from sc datasets + PRG + common genes
    import scanpy as sc
    base = str(PROJECT_ROOT / "data" / "processed" / "sc")
    genes = set()
    for ds in ["GSE125449", "GSE146115", "GSE149614", "GSE151530", "GSE98638"]:
        f = os.path.join(base, f"{ds}_raw.h5ad")
        if os.path.exists(f):
            a = sc.read_h5ad(f)
            genes.update(str(x) for x in a.var_names)
    prg = json.load(open(str(PROJECT_ROOT / "data" / "processed" / "peroxisome_gene_set.json")))
    genes.update(prg["genes"])
    print("total unique genes:", len(genes), flush=True)
    pos = mygene_pos(list(genes))
    json.dump(pos, open(os.path.join(ANN, "gene_pos.json"), "w"))
    print("GENE POS DONE n=", len(pos), flush=True)


if __name__ == "__main__":
    main()
