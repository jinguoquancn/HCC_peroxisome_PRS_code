#!/usr/bin/env python3
"""PubMed E-utilities helper: search + fetch MEDLINE records and print formatted citations."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import sys, json, time, urllib.parse, urllib.request, ssl, re

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def get(url):
    for _ in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60, context=CTX) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            time.sleep(2)
    return ""


def esearch(term, retmax=30, mindate=None, maxdate=None, sort="relevance"):
    p = {"db": "pubmed", "term": term, "retmode": "json", "retmax": str(retmax), "sort": sort}
    if mindate:
        p["mindate"], p["maxdate"], p["datetype"] = mindate, maxdate or "2026", "pdat"
    url = BASE + "esearch.fcgi?" + urllib.parse.urlencode(p)
    d = json.loads(get(url))
    return d.get("esearchresult", {}).get("idlist", [])


def efetch(ids):
    if not ids:
        return []
    url = BASE + "efetch.fcgi?" + urllib.parse.urlencode(
        {"db": "pubmed", "id": ",".join(ids), "retmode": "text", "rettype": "medline"})
    txt = get(url)
    recs, cur = [], None
    for line in txt.split("\n"):
        if line.startswith("PMID- "):
            if cur:
                recs.append(cur)
            cur = {"PMID": line[6:].strip(), "AU": [], "AID": []}
            continue
        if cur is None:
            continue
        if line.startswith("      ") or line.startswith("    "):
            if cur.get("_k"):
                cur[cur["_k"]] = (cur.get(cur["_k"], "") + " " + line.strip()).strip()
            continue
        m = re.match(r"^([A-Z0-9]{2,4})\s+-\s(.*)$", line)
        if not m:
            continue
        k, v = m.group(1), m.group(2).strip()
        cur["_k"] = k
        if k in ("AU", "FAU"):
            cur["AU"].append(v)
        elif k in ("AID", "LID"):
            cur["AID"].append(v)
        else:
            cur[k] = v
    if cur:
        recs.append(cur)
    return recs


def fmt(r):
    au = r["AU"][0] if r["AU"] else "?"
    au = re.sub(r"\s+", "", au.split(",")[0])
    year = (r.get("DP", "????") + " ")[:4]
    ti = r.get("TI", "").strip().rstrip(".")
    jt = r.get("JT", r.get("TA", "")).strip()
    vi = r.get("VI", "")
    ip = r.get("IP", "")
    pg = r.get("PG", "")
    vol = vi + (f"({ip})" if ip else "") + (f":{pg}" if pg else "")
    doi = ""
    for a in r.get("AID", []):
        if "doi" in a:
            doi = a.replace(" [doi]", "").strip()
    pid = "PMID: " + r["PMID"] + (f"; {doi}" if doi else "")
    return f'- {au} {year}. {ti}. {jt}. {vol}. {pid}.'


def run(term, retmax=30, mindate=None, maxdate=None, sort="relevance"):
    ids = esearch(term, retmax, mindate, maxdate, sort)
    recs = efetch(ids)
    print(f"### QUERY: {term}  (n={len(recs)})")
    for r in recs:
        print(fmt(r))
    print()
    return recs


if __name__ == "__main__":
    for t in sys.argv[1:]:
        run(t, retmax=int(__import__("os").environ.get("N", "30")))
