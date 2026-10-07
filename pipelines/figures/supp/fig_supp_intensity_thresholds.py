#!/usr/bin/env python3
"""FigS_intensity_thresholds: testing-intensity policies (typed + <=1, <=2, <=4 registered trials), change in per-disease AP
relative to masking with 95% paired disease-cluster bootstrap intervals (Table S22b), all scorers.
Rows: graph (Hetionet A-C, PrimeKG D-F); columns: outcome. Y groups = thresholds, markers = scorers.
Y labels give the number of negated pairs and the percentage of them that are recorded or approved indications."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import supp_common as C
S = C.S
S.setup()

R = C.intensity_rows()
POLS = ["typed_lowint1", "typed_lowint2", "typed_lowint4"]
THR = {"typed_lowint1": "≤1 trial", "typed_lowint2": "≤2 trials", "typed_lowint4": "≤4 trials"}
TITLES = {"random": "random edge", "compound": "compound disjoint", "e2": "external approvals"}
SCORERS = {"Hetionet": ["graph", "hybrid", "mf", "degree"], "PrimeKG": ["graph", "mf", "degree"]}

# number of negated pairs / percentage are properties of the policy and graph; check they are constant over outcome and scorer
nn = {}
for r in R:
    nn.setdefault((r["graph"], r["policy"]), set()).add((r["n"], r["pct"]))
assert all(len(v) == 1 for v in nn.values()), nn
NN = {k: next(iter(v)) for k, v in nn.items()}

fig, axes = plt.subplots(2, 3, figsize=(6.7, 5.1), sharey="row")
letters = iter("ABCDEF")
ys = np.arange(len(POLS))[::-1]
for i, g in enumerate(("Hetionet", "PrimeKG")):
    scs = SCORERS[g]
    offs = np.linspace(0.27, -0.27, len(scs))
    for j, o in enumerate(C.OUTCOME_KEYS):
        ax = axes[i][j]
        for k, sc in enumerate(scs):
            for y, p in zip(ys, POLS):
                row = [r for r in R if r["graph"] == g and r["policy"] == p and r["outcome"] == o and r["scorer"] == sc]
                assert len(row) == 1, (g, p, o, sc, len(row))
                d, lo, hi = row[0]["d_mask"]
                S.errorbar(ax, d, y + offs[k], lo, hi, S.COL[sc], S.MK[sc], ms=3.0)
        for y in ys[:-1]:
            ax.axhline(y - 0.5, color=S.LIGHT, lw=0.4, zorder=0)
        ax.axvline(0, color=S.MUTED, lw=0.6, ls=(0, (3, 2)), zorder=0)
        ax.set_ylim(-0.6, len(POLS) - 0.4)
        ax.set_yticks(ys)
        ax.set_yticklabels([f"{THR[p]}\n{NN[(g, p)][0]:,.0f} ({NN[(g, p)][1]:.0f}%)" for p in POLS])
        ax.tick_params(axis="y", length=0)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        S.panel(ax, next(letters), title=f"{g}, {TITLES[o]}")
        if j == 0:
            ax.set_ylabel("Registered trials per pair;\nnegated pairs (% recorded or approved)", fontsize=6)
        if i == 1:
            ax.set_xlabel("Δ per-disease AP vs masking")
order = ["graph", "hybrid", "mf", "degree"]
handles = [plt.Line2D([], [], marker=S.MK[s], color=S.COL[s], mfc=S.COL[s], ms=3.6, lw=0.8, label=S.LAB[s]) for s in order]
fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.55, 0.0))
fig.tight_layout(rect=(0, 0.04, 1, 1), w_pad=1.0, h_pad=1.6)
S.save(fig, C.OUT, "FigS_intensity_thresholds")
print("saved FigS_intensity_thresholds")
for k, v in NN.items():
    print(k, v)
