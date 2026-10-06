#!/usr/bin/env python
"""Figure 5: change in per-disease AP relative to no write-back for each policy (graph head and label-only MF),
Hetionet (A-C) and PrimeKG (D-F)."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as S
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

HET = sys.argv[1] if len(sys.argv) > 1 else "work/results_role"
PKG = sys.argv[2] if len(sys.argv) > 2 else "work/results_primekg"
OUT = sys.argv[3] if len(sys.argv) > 3 else "work/figures"
S.setup()
POL = ["flat_negative", "typed", "typed_scoped", "typed_role", "typed_role_scoped", "mask_all"]
PL = {"flat_negative": "Flat negative", "typed": "Typed", "typed_scoped": "Typed + scope", "typed_role": "Typed + role",
      "typed_role_scoped": "Typed + role + scope", "mask_all": "Mask all"}
fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.4), sharey=True)
letters = iter("ABCDEF")
for row, (res, gname) in enumerate([(HET, "Hetionet"), (PKG, "PrimeKG")]):
    for col, (key, title) in enumerate([("random", "held-out, random edge"), ("compound", "held-out, compound disjoint"), ("e2", "external approvals")]):
        ax = axes[row][col]; B = json.load(open(f"{res}/boot_{key}.json")); ys = list(range(len(POL)))[::-1]
        for k, m in enumerate(["graph", "mf"]):
            for y, p in zip(ys, POL):
                c = B["contrasts"][f"{p} - no_writeback [{m}]"]["macroAP"]
                S.errorbar(ax, c["diff"], y + (0.14 if k == 0 else -0.14), c["lo"], c["hi"], S.COL[m], S.MK[m],
                           label=S.LAB[m] if (row == 0 and col == 0 and y == ys[0]) else None, ms=3.2)
        ax.axvline(0, color=S.MUTED, lw=0.6, ls=(0, (3, 2)))
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4)); ax.set_yticks(ys); ax.set_yticklabels([PL[p] for p in POL]); ax.tick_params(axis="y", length=0)
        S.panel(ax, next(letters), title=f"{gname}, {title}")
        if row == 1:
            ax.set_xlabel("Δ per-disease AP")
fig.legend(*axes[0][0].get_legend_handles_labels(), loc="lower center", ncol=2, bbox_to_anchor=(0.55, -0.01))
fig.tight_layout(rect=(0, 0.04, 1, 1), w_pad=1.0, h_pad=1.6)
S.save(fig, OUT, "Figure5_policies")
print("saved")
