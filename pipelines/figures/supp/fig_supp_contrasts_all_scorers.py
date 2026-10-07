#!/usr/bin/env python3
"""FigS_contrasts_all_scorers: change in per-disease AP relative to no write-back (95% paired disease-cluster bootstrap
interval) for every scorer and policy. Hetionet (A-C): Table S5 (flat, typed, typed + scope, mask all), Table S18b (role policies),
Table S22b (typed + <=2 trials). PrimeKG (D-F): Table S16b, Table S22b (typed + <=2 trials)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import supp_common as C
S = C.S
S.setup()

POLS = ["flat_negative", "typed", "typed_scoped", "typed_role", "typed_role_scoped", "typed_lowint2", "mask_all"]
TITLES = {"random": "held-out, random edge", "compound": "held-out, compound disjoint", "e2": "external approvals"}
SC_ORDER = ["graph", "hybrid", "mf", "degree"]                  # top to bottom within each policy row
DATA = [("Hetionet", C.contrasts_hetionet(), ["graph", "hybrid", "mf", "degree"]),
        ("PrimeKG", C.contrasts_primekg(), ["graph", "mf", "degree"])]

fig, axes = plt.subplots(2, 3, figsize=(6.7, 7.2), sharey=True)
letters = iter("ABCDEF")
ys = np.arange(len(POLS))[::-1]
for r, (gname, D, scorers) in enumerate(DATA):
    offs = np.linspace(0.30, -0.30, len(scorers))
    for c, o in enumerate(C.OUTCOME_KEYS):
        ax = axes[r][c]
        for k, sc in enumerate(scorers):
            for y, p in zip(ys, POLS):
                d, lo, hi = D[(o, p, sc)]["perdis"]
                S.errorbar(ax, d, y + offs[k], lo, hi, S.COL[sc], S.MK[sc], ms=2.9)
        for y in ys[:-1]:
            ax.axhline(y - 0.5, color=S.LIGHT, lw=0.4, zorder=0)
        ax.axvline(0, color=S.MUTED, lw=0.6, ls=(0, (3, 2)), zorder=0)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.set_yticks(ys); ax.set_yticklabels([C.PL[p] for p in POLS]); ax.tick_params(axis="y", length=0)
        ax.set_ylim(-0.6, len(POLS) - 0.4)
        S.panel(ax, next(letters), title=f"{gname}, {TITLES[o]}")
        if r == 1:
            ax.set_xlabel("Δ per-disease AP")
handles = [plt.Line2D([], [], marker=S.MK[s], color=S.COL[s], mfc=S.COL[s], ms=3.6, lw=0.8, label=S.LAB[s]) for s in SC_ORDER]
fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.55, 0.0))
fig.tight_layout(rect=(0, 0.035, 1, 1), w_pad=1.0, h_pad=1.6)
S.save(fig, C.OUT, "FigS_contrasts_all_scorers")
print("saved FigS_contrasts_all_scorers")
