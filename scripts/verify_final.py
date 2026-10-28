# -*- coding: utf-8 -*-
"""Final verification of the 5 fixes across all main figures."""
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
        # 1) letters above axes, left-shifted, no overlap with title
        if exp:
            ax_list, labs, letters = exp
            for ax, lab, lt in zip(ax_list, labs, letters):
                if lt is None:
                    continue
                lbb = lt.get_window_extent(renderer)
                abb = ax.get_window_extent(renderer)
                lx0 = lbb.x0
                ax0 = abb.x0
                above = lbb.y0 >= abb.y1
                t_overlap = ax.get_title() and lbb.overlaps(ax.title.get_window_extent(renderer))
                print(f"[{name}] {lab}: above={above} left_offset={ax0-lx0:.0f}px title_overlap={t_overlap}")
        # 2) text overlap (70% shrink)
        issues = []
        for ax in fig.axes:
            texts = [t for t in ax.texts if t.get_text().strip()]
            texts += [t for t in ax.get_xticklabels() + ax.get_yticklabels()
                      if t.get_text().strip()]
            if ax.get_title():
                texts.append(ax.title)
            leg = ax.get_legend()
            if leg is not None:
                texts += [t for t in leg.get_texts() if t.get_text().strip()]
            n = len(texts)
            for i in range(n):
                for j in range(i + 1, n):
                    try:
                        if texts[i].get_window_extent(renderer).shrunk(0.7, 0.7).overlaps(
                           texts[j].get_window_extent(renderer).shrunk(0.7, 0.7)):
                            issues.append((texts[i].get_text()[:12], texts[j].get_text()[:12]))
                    except Exception:
                        pass
        print(f"[{name}] text-overlap: {'CLEAN' if not issues else issues[:8]}")
        plt.close(fig)
    except Exception as e:
        print(f"[{name}] ERR {type(e).__name__}: {e}")
    finally:
        mod.save = old


for mod, fn, nm in [(fig_bulk, fig_bulk.fig4, "Fig4"), (fig_bulk, fig_bulk.fig5, "Fig5"),
                    (fig_bulk, fig_bulk.fig6, "Fig6"), (fig_bulk, fig_bulk.fig7, "Fig7"),
                    (fig_bulk, fig_bulk.fig8, "Fig8"), (fig_bulk, fig_bulk.fig9, "Fig9"),
                    (fig_sc, fig_sc.fig1, "Fig1"), (fig_sc, fig_sc.fig2, "Fig2"),
                    (fig_sc, fig_sc.fig3, "Fig3")]:
    probe(mod, fn, nm)
