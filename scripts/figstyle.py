# -*- coding: utf-8 -*-
"""Shared plotting style and helpers for the HCC peroxisome paper figures."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.transforms import Bbox
import os

rng = np.random.default_rng(20260829)

OUT = str(PROJECT_ROOT / "figures")
DATA = str(PROJECT_ROOT / "data")

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7.5,
    "axes.titlesize": 8,
    "axes.labelsize": 7.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

OI = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
      "#56B4E9", "#F0E442", "#000000", "#8C6D31", "#999999"]
MM = 1 / 25.4
FULL_W = 183 * MM
HALF_W = 89 * MM


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".tif"), dpi=300,
                pil_kwargs={"compression": "tiff_lzw"},
                bbox_inches="tight", pad_inches=0.12)
    fig.savefig(os.path.join(OUT, name + ".svg"),
                bbox_inches="tight", pad_inches=0.12)
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=300,
                bbox_inches="tight", pad_inches=0.12)
    export = getattr(fig, "_panel_export", None)
    if export:
        axes_list, labels, letters = export
        renderer = fig.canvas.get_renderer()
        folder = os.path.join(OUT, name.lower())
        os.makedirs(folder, exist_ok=True)
        dpi = fig.dpi
        n = 0
        for ax, lab, ltxt in zip(axes_list, labels, letters):
            if not lab:
                continue
            try:
                bb = ax.get_tightbbox(renderer)
            except Exception:
                bb = None
            if bb is None:
                pos = ax.get_position()
                bb = pos.transformed(fig.transFigure)
            if ltxt is not None:
                try:
                    lb = ltxt.get_window_extent(renderer)
                    bb = Bbox.union([bb, lb])
                except Exception:
                    pass
            bbox_in = Bbox.from_bounds(bb.x0 / dpi, bb.y0 / dpi,
                                       bb.width / dpi, bb.height / dpi).padded(0.06)
            fig.savefig(os.path.join(folder, f"{name}_{lab}.svg"),
                        bbox_inches=bbox_in)
            n += 1
        print(f"  panels: {n} -> {folder}")
    plt.close(fig)
    print("saved", name)


def _wrap_title(ax, title, renderer, fontsize=None, max_frac=0.92):
    """Wrap a long axes title to fit within the panel width (<= max_frac of
    axes width) so it never collides with the panel letter."""
    if fontsize is None:
        fontsize = matplotlib.rcParams.get("axes.titlesize", 8)
    from matplotlib.font_manager import FontProperties
    prop = FontProperties(size=fontsize)
    abb = ax.get_window_extent(renderer)
    max_w = abb.width * max_frac
    if max_w <= 0:
        return title
    words = title.split()
    lines, cur = [], ""
    for w in words:
        trial = (cur + " " + w).strip()
        tw, _, _ = renderer.get_text_width_height_descent(trial, prop, ismath=False)
        if tw <= max_w or not cur:
            cur = trial
        else:
            lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)
    return "\n".join(lines) if len(lines) > 1 else title


def label_axes(axes, labels):
    """Uppercase bold panel letters placed ABOVE the upper-left corner of each
    axes (outside the plotting area, so they never overlap panel content).
    Long titles are auto-wrapped to the panel width. Call AFTER
    fig.tight_layout() so the axes positions are final."""
    if axes is None or len(axes) == 0:
        return
    axes = list(np.asarray(axes).ravel())
    fig = axes[0].figure
    try:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
    except Exception:
        renderer = None
    letters = []
    for ax, lab in zip(axes, labels):
        if not lab:
            letters.append(None)
            continue
        lab = lab.upper()
        pos = ax.get_position()
        # auto-wrap overlong titles so they stay within the panel
        if renderer is not None:
            title_txt = ax.get_title().strip()
            if title_txt:
                ax.set_title(_wrap_title(ax, title_txt, renderer))
        # place above the axes, left of the y-axis labels' left edge
        # (i.e. aligned with the left edge of the panel INCLUDING tick labels)
        x_frac = pos.x0
        if renderer is not None:
            try:
                yl = [t for t in ax.get_yticklabels() if t.get_text().strip()]
                if yl:
                    lefts = [t.get_window_extent(renderer).x0 for t in yl]
                    left_px = min(lefts)
                else:
                    left_px = ax.get_window_extent(renderer).x0
                fig_w_px = fig.get_figwidth() * fig.dpi
                x_frac = (left_px - 6) / fig_w_px
            except Exception:
                pass
        t = fig.text(x_frac, pos.y1 + 0.008, lab,
                     fontsize=9, fontweight="bold", ha="left", va="bottom",
                     zorder=40)
        # lift the letter if it collides with a (wrapped) title above the panel
        if renderer is not None:
            try:
                fig.canvas.draw()
                lbb = t.get_window_extent(renderer)
                tbb = ax.title.get_window_extent(renderer)
                if lbb.overlaps(tbb):
                    dy_px = tbb.y1 - lbb.y0 + 4
                    fig_h_in = fig.get_figheight()
                    t.set_y(t.get_position()[1] + dy_px / (fig_h_in * fig.dpi))
            except Exception:
                pass
        letters.append(t)
    fig._panel_export = (list(axes), list(labels), letters)


def sig_bracket(ax, x1, x2, y, p, fs=8, lw=0.9, color="black",
                min_stars=False):
    """Horizontal bracket with significance stars above, between x1 and x2
    (standard journal style: horizontal line with two short downward ticks at
    the ends, stars centered above). If not significant and min_stars=False,
    a grey 'ns' is shown instead."""
    if p < 0.05:
        stars_txt = stars(p)
        star_color = color
    else:
        stars_txt = "ns" if not min_stars else ""
        star_color = "#999999"
        fs = fs - 2
    ax.plot([x1, x2], [y, y], color=color, lw=lw, zorder=6, clip_on=False)
    ylo, yhi = ax.get_ylim()
    ticklen = (yhi - ylo) * 0.008
    ax.plot([x1, x1], [y - ticklen, y], color=color, lw=lw, zorder=6, clip_on=False)
    ax.plot([x2, x2], [y - ticklen, y], color=color, lw=lw, zorder=6, clip_on=False)
    if stars_txt:
        ax.text((x1 + x2) / 2, y + ticklen * 1.2, stars_txt, ha="center",
                va="bottom", fontsize=fs, fontweight="bold", color=star_color,
                zorder=7, clip_on=False)


def savecsv(df, name):
    df.to_csv(os.path.join(DATA, name), index=False)
    print("csv", name)


def km_step(t, e):
    order = np.argsort(t)
    t, e = t[order], e[order]
    n = len(t)
    surv, xs, ys = 1.0, [0], [1.0]
    at_risk = n
    for i in range(n):
        if e[i] == 1:
            surv *= (at_risk - 1) / at_risk
            xs.append(t[i]); ys.append(surv)
        at_risk -= 1
    xs.append(max(t) + 0.001); ys.append(ys[-1])
    return np.array(xs), np.array(ys)


def draw_km_from_surv(ax, time, event, group, colors=("#0072B2", "#D55E00"),
                      xlabel="Time (years)", ylabel="Overall survival probability",
                      labels=("Low", "High"), title=None, p=None, hr=None,
                      risk_table=False, risk_ticks=4, fs_risk=5.6):
    """Draw KM curves from real (time, event, group) arrays.

    risk_table=True adds a numbers-at-risk block underneath the panel, which
    is what readers need to judge the tails of the curves.
    """
    from matplotlib.transforms import blended_transform_factory
    g = np.asarray(group)
    for k, lab, c in [(0, labels[0], colors[0]), (1, labels[1], colors[1])]:
        m = g == k
        t, e = np.asarray(time)[m], np.asarray(event)[m].astype(int)
        if len(t) == 0:
            continue
        xs, ys = km_step(t, e)
        ax.step(xs, ys, where="post", color=c, lw=1.2, label=f"{lab} (n={len(t)})")
    if p is not None:
        txt = f"Log-rank P {'< 0.001' if p < 0.001 else '= %.3f' % p}"
        if hr is not None:
            txt += f"\nHR = {hr:.2f}"
        ax.text(0.97, 0.96, txt, transform=ax.transAxes, ha="right", va="top",
                fontsize=7,
                bbox=dict(facecolor="white", edgecolor="none", pad=1.2, alpha=0.85))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(-0.02, 1.02)
    ax.legend(frameon=False, loc="lower left", handlelength=1.4)
    if title:
        ax.set_title(title)

    if risk_table:
        tt = np.asarray(time, float)
        tmax = np.nanmax(tt)
        ticks = [v for v in ax.get_xticks() if 0 <= v <= tmax * 1.001]
        if len(ticks) < 3:
            ticks = list(np.linspace(0, tmax, risk_ticks))
        ax.set_xticks(ticks)
        ax.set_xticklabels(["%g" % (round(v, 1) if v >= 1 else round(v, 2)) for v in ticks])
        trans = blended_transform_factory(ax.transData, ax.transAxes)
        ax.text(-0.02, -0.15, "No. at risk", transform=ax.transAxes,
                ha="right", va="center", fontsize=fs_risk, style="italic")
        for row, (k, lab, c) in enumerate([(0, labels[0], colors[0]),
                                           (1, labels[1], colors[1])]):
            m = g == k
            counts = [int(np.sum(tt[m] >= v)) for v in ticks]
            y = -0.15 - 0.085 * (row + 1)
            ax.text(-0.02, y, lab, transform=ax.transAxes, ha="right",
                    va="center", fontsize=fs_risk, color=c)
            for v, n in zip(ticks, counts):
                ax.text(v, y, str(n), transform=trans, ha="center",
                        va="center", fontsize=fs_risk)



def roc_curve_from_auc(n_pos, n_neg, auc, npts=200):
    from math import sqrt
    from statistics import NormalDist
    nd = NormalDist()
    lo, hi = 0.05, 8.0
    for _ in range(60):
        d = (lo + hi) / 2
        a = nd.cdf(d / sqrt(2))
        if a < auc:
            lo = d
        else:
            hi = d
    fpr = np.linspace(0.001, 0.999, npts)
    tpr = np.array([nd.cdf(d + nd.inv_cdf(p)) for p in fpr])
    tpr[0], tpr[-1] = 0, 1
    fpr[0], fpr[-1] = 0, 1
    return fpr, tpr


def stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def dot_plot(ax, mat, rowlabels, collabels, cmap="RdYlBu_r", vmin=None,
             vmax=None, cbar_label="Average expression", pct=None):
    """mat: genes x celltypes; pct optional same shape for dot size."""
    nr, nc = mat.shape
    if pct is None:
        pct = np.clip(np.abs(mat) / (np.abs(mat).max() + 1e-9), 0, 1)
    vmax_ = vmax or np.nanmax(np.abs(mat))
    vmin_ = vmin if vmin is not None else -vmax_
    for i in range(nr):
        for j in range(nc):
            val = mat[i, j]
            if np.isnan(val):
                continue
            ax.scatter(j, i, s=18 + 130 * max(pct[i, j], 0.05),
                       c=[plt.get_cmap(cmap)((val - vmin_) / (vmax_ - vmin_ + 1e-9))],
                       edgecolors="black", linewidths=0.3, zorder=3)
    ax.set_xticks(range(nc))
    ax.set_xticklabels(collabels, rotation=45, ha="right")
    ax.set_yticks(range(nr))
    ax.set_yticklabels(rowlabels)
    ax.set_xlim(-0.7, nc - 0.3)
    ax.set_ylim(nr - 0.3, -0.7)
    ax.xaxis.set_ticks_position("top")
    for s in ax.spines.values():
        s.set_visible(False)
    import matplotlib.colors as mcolors
    import matplotlib.cm as cm
    norm = mcolors.Normalize(vmin=vmin_, vmax=vmax_)
    cax = ax.inset_axes([1.02, 0.3, 0.03, 0.4])
    cb = ax.figure.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.ax.tick_params(labelsize=6)
    cb.set_label(cbar_label, fontsize=6.5)


def box_pairs(ax, data_dict, colors=None, jitter=0.06, max_dots=36):
    labels = list(data_dict.keys())
    for i, (lab, vals) in enumerate(data_dict.items()):
        vals = np.asarray(vals, float)
        vals = vals[~np.isnan(vals)]
        bp = ax.boxplot(vals, positions=[i], widths=0.55, patch_artist=True,
                        showfliers=False, medianprops=dict(color="black", lw=0.8),
                        boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
        if colors:
            bp["boxes"][0].set_facecolor(colors[i % len(colors)])
            bp["boxes"][0].set_alpha(0.65)
        n = min(len(vals), max_dots)
        xs = rng.normal(i, jitter, n)
        ax.scatter(xs, rng.choice(vals, n), s=2, color="black", alpha=0.22, zorder=3)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    return labels
