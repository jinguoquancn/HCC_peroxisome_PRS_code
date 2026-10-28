# -*- coding: utf-8 -*-
"""Probe GEO series metadata (title/platform/sample count/supplementary files) and
write a compact report to data/raw/geo_probe.txt"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, re, json, urllib.request, concurrent.futures as cf

UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
RAW = str(PROJECT_ROOT / "data" / "raw")
ACC = ["GSE14520", "GSE36376", "GSE25097", "GSE10186", "GSE16757", "GSE39791",
       "GSE54236", "GSE60502", "GSE76427", "GSE87630", "GSE116174",
       "GSE149614", "GSE151530", "GSE156625", "GSE156337", "GSE125449",
       "GSE146115", "GSE98638", "GSE166635", "GSE162616"]


def get_txt(u, t=45):
    try:
        return urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=t).read().decode("utf8", "ignore")
    except Exception as e:
        return "ERR " + str(e)


def head_size(u, t=30):
    try:
        r = urllib.request.urlopen(urllib.request.Request(u, headers=UA, method="HEAD"), timeout=t)
        return int(r.headers.get("Content-Length", 0))
    except Exception:
        return -1


def probe(acc):
    s = get_txt(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full")
    out = {"acc": acc}
    if s.startswith("ERR"):
        out["error"] = s
        return out
    for key, name in [("!Series_title", "title"), ("!Series_platform_id", "platform"),
                      ("!Series_sample_count", "n_sample"), ("!Series_summary", "summary"),
                      ("!Series_overall_design", "design")]:
        m = re.search(rf"{key} = (.+)", s)
        if m:
            out[name] = m.group(1)[:400]
    # pubmed ids
    out["pmid"] = re.findall(r"!Series_pubmed_id = (\d+)", s)
    # supplementary files
    sup = re.findall(r"!Series_supplementary_file = (.+)", s)
    with cf.ThreadPoolExecutor(8) as ex:
        sizes = list(ex.map(head_size, sup))
    out["supp"] = [(f.split("/")[-1], f, sz) for f, sz in zip(sup, sizes)]
    # series matrix size
    mdir = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/matrix/"
    for fn in [f"{acc}_series_matrix.txt.gz", f"{acc}-GPL570_series_matrix.txt.gz",
               f"{acc}-GPL571_series_matrix.txt.gz", f"{acc}-GPL3921_series_matrix.txt.gz"]:
        sz = head_size(mdir + fn)
        if sz > 0:
            out.setdefault("matrix", []).append((fn, sz))
    return out


if __name__ == "__main__":
    os.makedirs(RAW, exist_ok=True)
    res = []
    with cf.ThreadPoolExecutor(5) as ex:
        for r in ex.map(probe, ACC):
            res.append(r)
            print("###", r["acc"], "|", r.get("title", "")[:120], "| PMID", r.get("pmid"))
            print("   platform:", r.get("platform"), " n:", r.get("n_sample"))
            for fn, sz in r.get("matrix", []):
                print(f"   MATRIX {fn} {sz/1e6:.1f} MB")
            tot = 0
            for fn, url, sz in r.get("supp", [])[:40]:
                if sz > 0:
                    tot += sz
                print(f"   SUPP {fn} {sz/1e6:.1f} MB")
            print(f"   [supp n={len(r.get('supp',[]))} total={tot/1e6:.0f} MB]")
            print()
    json.dump(res, open(os.path.join(RAW, "geo_probe.json"), "w"), indent=1)
