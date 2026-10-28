# -*- coding: utf-8 -*-
"""Shared gene-id utilities: Ensembl->symbol via mygene batch POST (q array)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import json, os, time, urllib.request

_cache = {}


def mygene_symbols(ids, scopes="symbol"):
    """Batch-query mygene for symbols. Returns {query_id: symbol}."""
    ids = [str(i) for i in ids]
    out = {}
    for st in range(0, len(ids), 800):
        chunk = ids[st:st + 800]
        body = json.dumps({"q": chunk, "scopes": scopes, "fields": "symbol",
                           "species": "human"}).encode()
        for t in range(6):
            try:
                req = urllib.request.Request("https://mygene.info/v3/query", data=body,
                                             headers={"Content-Type": "application/json",
                                                      "User-Agent": "Mozilla/5.0"})
                res = json.load(urllib.request.urlopen(req, timeout=180))
                for r in res:
                    if "notfound" not in r and r.get("symbol"):
                        out[r.get("query")] = r["symbol"]
                break
            except Exception as e:
                print("  mygene retry", st, e, flush=True)
                time.sleep(4 + 4 * t)
    return out


def to_symbols(var_names):
    """Convert an index of gene ids to symbols when needed.
    Returns (new_index, converted_bool)."""
    v = [str(x) for x in var_names]
    if v and v[0].startswith("ENSG"):
        key = "ensg"
        d = mygene_symbols(v, scopes="ensembl.gene")
        new = [d.get(x, x) for x in v]
        return new, True
    return v, False
