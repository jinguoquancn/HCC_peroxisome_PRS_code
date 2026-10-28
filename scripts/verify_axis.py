# -*- coding: utf-8 -*-
"""Precise overflow check for Fig1B / Fig4B (barh error bars)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import PROJECT_ROOT
import numpy as np
from PIL import Image


def check(fname, x0, x1, ytop, ybot, name):
    a = np.asarray(Image.open(f"figures/{fname}.png").convert("L"))
    h, w = a.shape
    reg = a[int(h * ytop):int(h * ybot), int(w * x0):int(w * x1)]
    colsum = (reg < 150).sum(axis=0)
    # y-axis = the tallest dark column (usually full height)
    if not colsum.max():
        print(f"{name}: nothing dark")
        return
    axis_x = int(np.argmax(colsum))
    # left of axis: exclude the tick labels band (first ~28px) and check for
    # a horizontal bar that spans wide continuously (error-bar overhang)
    left = reg[:, :max(axis_x - 3, 0)]
    dark = left < 120
    # error bar crossing = dark row whose longest run reaches the very left edge
    # tick labels are short runs, error bars are continuous ~10-20px
    maxrun = 0
    at_left = 0
    for ri, row in enumerate(dark):
        idx = np.where(row)[0]
        if len(idx):
            runs = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
            for r in runs:
                if len(r) > maxrun:
                    maxrun = len(r)
                if r[0] <= 2:  # touches left boundary of region
                    at_left += 1
    status = "OK" if maxrun < 14 and at_left == 0 else "CHECK"
    print(f"{name}: axis_col={axis_x} (region h={reg.shape[0]}), maxrun={maxrun}px, "
          f"touches_left={at_left} -> {status}")


if __name__ == "__main__":
    check("Fig1_scRNA_overview", 0.50, 0.99, 0.05, 0.47, "Fig1B")
    check("Fig4_bulk_PRS_building", 0.50, 0.99, 0.05, 0.47, "Fig4B")
