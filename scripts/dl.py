# -*- coding: utf-8 -*-
"""Utility helpers: throttled parallel downloader with retries."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import os, time, threading, urllib.request, urllib.error, json, sys

RAW = str(PROJECT_ROOT / "data" / "raw")
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}


def _head(url):
    req = urllib.request.Request(url, headers=UA, method="HEAD")
    with urllib.request.urlopen(req, timeout=90) as r:
        return int(r.headers["Content-Length"]), r.headers.get("Accept-Ranges")


def par_dl(url, out, nthreads=4, tries=6, log=print):
    """Download url -> RAW/out using ranged parallel GET with retries."""
    out = os.path.join(RAW, out)
    if os.path.exists(out) and os.path.getsize(out) > 0:
        log(f"[skip] {os.path.basename(out)} exists ({os.path.getsize(out)/1e6:.2f} MB)")
        return out
    size, acc = _head(url)
    log(f"[get ] {os.path.basename(out)}  {size/1e6:.2f} MB  ranges={acc}")
    if acc != "bytes":
        nthreads = 1
    chunk = (size + nthreads - 1) // nthreads
    bounds = [(i * chunk, min(size - 1, (i + 1) * chunk - 1)) for i in range(nthreads)]

    def worker(idx):
        a, b = bounds[idx]
        part = f"{out}.part{idx}"
        if os.path.exists(part):
            have = os.path.getsize(part)
        else:
            have = 0
        for t in range(tries):
            try:
                hdr = dict(UA)
                mode = "wb"
                if have and have < (b - a + 1):
                    hdr["Range"] = f"bytes={a+have}-{b}"
                    mode = "ab"
                else:
                    hdr["Range"] = f"bytes={a}-{b}"
                    have = 0
                    mode = "wb"
                req = urllib.request.Request(url, headers=hdr)
                with urllib.request.urlopen(req, timeout=1200) as r, open(part, mode) as f:
                    while True:
                        d = r.read(1 << 15)
                        if not d:
                            break
                        f.write(d)
                if os.path.getsize(part) >= (b - a + 1):
                    return True
                have = os.path.getsize(part)
            except Exception as e:
                time.sleep(3 + 4 * t)
        return False

    t0 = time.time()
    ths = [threading.Thread(target=worker, args=(i,)) for i in range(nthreads)]
    for t in ths:
        t.start()
    for t in ths:
        t.join()
    with open(out, "wb") as f:
        for i in range(nthreads):
            p = f"{out}.part{i}"
            if not os.path.exists(p):
                raise RuntimeError(f"missing part {p}")
            with open(p, "rb") as g:
                f.write(g.read())
            os.remove(p)
    dt = time.time() - t0
    log(f"[ok  ] {os.path.basename(out)} {os.path.getsize(out)/1e6:.2f} MB in {dt:.0f} s "
        f"({os.path.getsize(out)/1e6/max(dt,1)*1000:.0f} KB/s)")
    return out


def get_json(url, tries=4):
    for t in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.loads(r.read())
        except Exception:
            time.sleep(3 + 4 * t)
    raise RuntimeError("failed " + url)


def post_json(url, payload, tries=4, timeout=600):
    body = json.dumps(payload).encode()
    for t in range(tries):
        try:
            hdr = dict(UA)
            hdr["Content-Type"] = "application/json"
            req = urllib.request.Request(url, data=body, headers=hdr, method="POST")
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read())
        except Exception as e:
            if t == tries - 1:
                raise
            time.sleep(5 + 5 * t)
