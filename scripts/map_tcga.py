# -*- coding: utf-8 -*-
"""Map TCGA RSEM Entrez IDs -> gene symbols (mygene.info) and save expression matrix."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import pandas as pd, os, json, urllib.request, urllib.parse, time
BASE = str(PROJECT_ROOT)
EX = os.path.join(BASE, "data/raw/extracted_tcga")
OUT = os.path.join(BASE, "data/processed")
os.makedirs(OUT, exist_ok=True)

# short-path copy of the RSEM matrix (long path breaks Windows open)
F = os.path.join(BASE, "data/raw", "TCGA_LIHC_rsem_normalized.data.txt")
print("reading", F, flush=True)
df = pd.read_csv(F, sep="\t", index_col=0)
df = df[df.index != "gene_id"]   # drop the 'gene_id' placeholder row
df.index = [str(i) for i in df.index]
ent = [i.split("|")[-1] for i in df.index]

sym = {}
for st in range(0, len(ent), 1000):
    chunk = ent[st:st + 1000]
    body = json.dumps({"q": chunk, "scopes": "entrezgene", "fields": "symbol",
                       "species": "human"}).encode()
    for t in range(6):
        try:
            req = urllib.request.Request("https://mygene.info/v3/query", data=body,
                                         headers={"Content-Type": "application/json",
                                                  "User-Agent": "Mozilla/5.0"})
            res = json.load(urllib.request.urlopen(req, timeout=180))
            for r in res:
                if "notfound" not in r and r.get("symbol"):
                    sym[r.get("query")] = r["symbol"]
            break
        except Exception as e:
            print("retry", st, e, flush=True)
            time.sleep(4 + 4 * t)
    print("mapped", st + 1000, flush=True)
print("total mapped", len(sym), flush=True)

df["symbol"] = [sym.get(e, "") for e in ent]
df2 = df[df["symbol"] != ""].drop_duplicates("symbol", keep="first").set_index("symbol")
df2 = df2.apply(pd.to_numeric, errors="coerce")
print("final expr:", df2.shape, flush=True)
df2.to_csv(os.path.join(OUT, "TCGA_LIHC_rsem_normalized_matrix.csv"))

cols = [c for c in df2.columns]
def tissue_of(c):
    try:
        p = c.split("-")
        return "tumor" if len(p) > 3 and p[3][:2] in ("01", "04", "06", "07", "08", "09") else "normal"
    except Exception:
        return "normal"
tiss = [tissue_of(c) for c in cols]
print(pd.Series(tiss).value_counts(), flush=True)
json.dump({"sample": cols, "tissue": tiss}, open(os.path.join(OUT, "TCGA_LIHC_sample_info.json"), "w"))
print("DONE", flush=True)
