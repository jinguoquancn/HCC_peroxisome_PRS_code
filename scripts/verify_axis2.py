# -*- coding: utf-8 -*-
"""Verify Fig1B asymmetric error bars and Fig4B x-origin."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import sys
for p in (str(PROJECT_ROOT), str(PROJECT_ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import fig_bulk, fig_sc

captured = {}
def fs_(fig, name):
    captured["fig"] = fig

old = fig_sc.save
fig_sc.save = fs_
fig_sc.fig1()
fig_sc.save = old
ax = captured["fig"].axes[1]
print("Fig1B xlim:", round(float(ax.get_xlim()[0]), 3), "-", round(float(ax.get_xlim()[1]), 3))
# bars: asymmetric xerr -> error lines should start at mean (right only)
for c in ax.containers:
    if hasattr(c, "datavalues"):
        vals = getattr(c, "datavalues")
        print("  bar data range:", round(float(vals.min()), 3), "-", round(float(vals.max()), 3))
plt.close(captured["fig"])

captured = {}
old = fig_bulk.save
fig_bulk.save = fs_
fig_bulk.fig4()
fig_bulk.save = old
ax = captured["fig"].axes[1]
print("Fig4B xlim:", round(float(ax.get_xlim()[0]), 3), "-", round(float(ax.get_xlim()[1]), 3),
      "(should start at 0)")
plt.close(captured["fig"])
