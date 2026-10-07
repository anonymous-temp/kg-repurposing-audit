#!/usr/bin/env python3
"""FigS_placebo_all: change in per-disease AP relative to masking when a negated set is added on top of masking (Table S22a).
For each negated set: stopped pairs (filled), degree-matched placebo of the same size (open), uniform placebo (cross).
Rows: graph x scorer; columns: outcome."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import supp_common as C
S = C.S
S.setup()

R = C.placebo_rows()
ROWS = [("Hetionet", "graph"), ("Hetionet", "mf"), ("PrimeKG", "graph"), ("PrimeKG", "mf")]
SETS = {g: [p for p in C.POL_ORDER if any(r["graph"] == g and r["set"] == p for r in R)] for g in ("Hetionet", "PrimeKG")}
TITLES = {"random": "Random edge", "compound": "Compound disjoint", "e2": "External approvals"}
OFF = {"real": 0.24, "degree": 0.0, "uniform": -0.24}

fig = plt.figure(figsize=(6.55, 8.0))
nH, nP = len(SETS["Hetionet"]), len(SETS["PrimeKG"])
gs = fig.add_gridspec(4, 3, height_ratios=[nH + 0.8, nH + 0.8, nP + 0.8, nP + 0.8], hspace=0.5, wspace=0.14,
                      left=0.2, right=0.985, top=0.955, bottom=0.075)
letters = iter("ABCDEFGHIJKL")
axes = {}
for i, (g, sc) in enumerate(ROWS):
    sets = SETS[g]
    ys = np.arange(len(sets))[::-1]
    for j, o in enumerate(C.OUTCOME_KEYS):
        ax = fig.add_subplot(gs[i, j]); axes[(i, j)] = ax
        sel = {r["set"]: r for r in R if r["graph"] == g and r["outcome"] == o and r["scorer"] == sc}
        assert set(sel) == set(sets), (g, o, sc, set(sets) - set(sel))
        for y, p in zip(ys, sets):
            r = sel[p]
            for kind in ("real", "degree", "uniform"):
                d, lo, hi = r[kind]
                kw = dict(fmt=S.MK[sc] if kind != "uniform" else "x", ms=3.1 if kind != "uniform" else 3.4, elinewidth=0.7, capsize=1.3, capthick=0.7,
                          color=S.COL[sc] if kind != "uniform" else S.MUTED, mec=S.COL[sc] if kind != "uniform" else S.MUTED,
                          mfc=S.COL[sc] if kind == "real" else "white")
                if kind == "uniform":
                    kw.pop("mfc"); kw["mew"] = 0.9
                ax.errorbar(d, y + OFF[kind], xerr=[[d - lo], [hi - d]], **kw)
        for y in ys[:-1]:
            ax.axhline(y - 0.5, color=S.LIGHT, lw=0.4, zorder=0)
        ax.axvline(0, color=S.MUTED, lw=0.6, ls=(0, (3, 2)), zorder=0)
        lo_all = min([0.0] + [r[k][1] for r in sel.values() for k in ("real", "degree", "uniform")])
        hi_all = max([0.0] + [r[k][2] for r in sel.values() for k in ("real", "degree", "uniform")])
        pad = 0.06 * (hi_all - lo_all)
        ax.set_xlim(lo_all - pad, hi_all + pad)
        ax.set_ylim(-0.65, len(sets) - 0.35)
        ax.set_yticks(ys)
        ax.set_yticklabels([C.PL[p] for p in sets] if j == 0 else [])
        ax.tick_params(axis="y", length=0)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.tick_params(axis="x", labelsize=6)
        if i == 3:
            ax.set_xlabel("Δ per-disease AP vs masking")
        S.panel(ax, next(letters), x=-0.02, y=1.03)
        if i == 0:
            ax.set_title(TITLES[o], fontsize=7, pad=11, fontweight="normal")
fig.canvas.draw()
for i, (g, sc) in enumerate(ROWS):                                    # row labels at the far left
    pos = axes[(i, 0)].get_position()
    fig.text(0.02, (pos.y0 + pos.y1) / 2, f"{g}\n{S.LAB[sc]}", rotation=90, ha="center", va="center", fontsize=7,
             color=S.COL[sc], fontweight="bold")
hand = [Line2D([], [], marker="o", color="k", mfc="k", ls="none", ms=3.4, label="Stopped pairs"),
        Line2D([], [], marker="o", color="k", mfc="white", ls="none", ms=3.4, label="Placebo, same drug and disease counts"),
        Line2D([], [], marker="x", color=S.MUTED, ls="none", ms=3.6, mew=0.9, label="Placebo, uniform")]
fig.legend(handles=hand, loc="lower center", ncol=3, bbox_to_anchor=(0.55, -0.003))
S.save(fig, C.OUT, "FigS_placebo_all")
print("saved FigS_placebo_all")
