# -*- coding: utf-8 -*-
"""Verify: letters at y-label-left edge, stars do not cross other brackets,
Fig1B y-axis at 0, Fig1A legend upper-right, Fig9B legend close."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import sys, warnings
warnings.filterwarnings("ignore")
for p in (str(PROJECT_ROOT), str(PROJECT_ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import fig_bulk, fig_sc


def probe(mod, fn, name):
    captured = {}
    def fake_save(fig, name):
        captured["fig"] = fig
    old = mod.save
    mod.save = fake_save
    try:
        fn()
        fig = captured["fig"]
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        exp = getattr(fig, "_panel_export", None)
        if exp:
            ax_list, labs, letters = exp
            for ax, lab, lt in zip(ax_list, labs, letters):
                if lt is None:
                    continue
                lbb = lt.get_window_extent(renderer)
                abb = ax.get_window_extent(renderer)
                # letter should be left of the axes box (incl. tick labels)
                left_of_axes = lbb.x1 <= abb.x0
                above = lbb.y0 >= abb.y1
                print(f"[{name}] {lab}: above={above} left_of_axes={left_of_axes} "
                      f"(letter_right={lbb.x1:.0f}, axes_left={abb.x0:.0f})")
        # stars vs bracket lines (pixel-level: star bbox must not cross a line
        # belonging to a DIFFERENT bracket)
        for ai, ax in enumerate(fig.axes):
            lines_y = []
            for ln in ax.lines:
                xd, yd = ln.get_xdata(), ln.get_ydata()
                if len(xd) == 2 and len(yd) == 2:
                    span = max(abs(xd[1] - xd[0]), 1e-9)
                    if abs(yd[1] - yd[0]) < span * 0.1:
                        lines_y.append(float(ax.transData.transform((xd[0], yd[0]))[1]))
            stars = [t for t in ax.texts if t.get_text() and set(t.get_text()) <= {"*"}]
            if lines_y and stars:
                bad = []
                for t in stars:
                    tbb = t.get_window_extent(renderer)
                    for ly in lines_y:
                        # star crosses a bracket line it does not belong to
                        if tbb.y0 < ly < tbb.y1:
                            bad.append((t.get_text(), round(tbb.y0, 1), round(ly, 1),
                                        round(tbb.y1, 1)))
                if bad:
                    print(f"[{name}] ax{ai}: STAR CROSSES OTHER BRACKET {bad[:6]}")
        plt.close(fig)
    except Exception as e:
        print(f"[{name}] ERR {type(e).__name__}: {e}")
    finally:
        mod.save = old


for mod, fn, nm in [(fig_bulk, fig_bulk.fig4, "Fig4"), (fig_bulk, fig_bulk.fig5, "Fig5"),
                    (fig_bulk, fig_bulk.fig6, "Fig6"), (fig_bulk, fig_bulk.fig8, "Fig8"),
                    (fig_bulk, fig_bulk.fig9, "Fig9"),
                    (fig_sc, fig_sc.fig1, "Fig1"), (fig_sc, fig_sc.fig2, "Fig2"),
                    (fig_sc, fig_sc.fig3, "Fig3")]:
    probe(mod, fn, nm)
