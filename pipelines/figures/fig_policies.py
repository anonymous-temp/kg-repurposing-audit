#!/usr/bin/env python
"""Figure 4: change in per-disease AP relative to no write-back, for every write-back policy,
in Hetionet and PrimeKG (held-out treatments in two tasks and external approved indications)."""
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

HET = sys.argv[1] if len(sys.argv) > 1 else "work/results_role"
PKG = sys.argv[2] if len(sys.argv) > 2 else "work/results_primekg"
OUT = sys.argv[3] if len(sys.argv) > 3 else "work/figures"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "pdf.fonttype": 42, "axes.linewidth": 0.6})
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
COL = {"graph": "#2a78d6", "mf": "#1baf7a"}; MK = {"graph": "o", "mf": "D"}; LAB = {"graph": "Graph head", "mf": "Label-only MF"}
POL = ["flat_negative", "typed", "typed_scoped", "typed_role", "typed_role_scoped", "mask_all"]
PL = {"flat_negative": "Flat negative", "typed": "Typed", "typed_scoped": "Typed + scope", "typed_role": "Typed + role",
      "typed_role_scoped": "Typed + role + scope", "mask_all": "Mask all"}
fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.6), sharey=True)
letters = iter("abcdef")
for row, (res, gname) in enumerate([(HET, "Hetionet"), (PKG, "PrimeKG")]):
    for col, (key, title) in enumerate([("random", "E1 random edge"), ("compound", "E1 compound disjoint"),
                                        ("e2", "E2 external approvals")]):
        ax = axes[row][col]; B = json.load(open(f"{res}/boot_{key}.json"))
        ys = list(range(len(POL)))[::-1]
        for k, m in enumerate(["graph", "mf"]):
            off = 0.14 if k == 0 else -0.14
            for y, p in zip(ys, POL):
                c = B["contrasts"].get(f"{p} - no_writeback [{m}]")
                if c is None:
                    continue
                c = c["macroAP"]
                ax.plot([c["lo"], c["hi"]], [y + off] * 2, color=COL[m], lw=1.2, solid_capstyle="round")
                ax.plot(c["diff"], y + off, MK[m], ms=4.0, color=COL[m], mec="white", mew=0.8, label=LAB[m] if (y == ys[0] and row == 0 and col == 0) else None)
        ax.axvline(0, color=MUTED, lw=0.7); ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.set_yticks(ys); ax.set_yticklabels([PL[p] for p in POL])
        ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
        ax.set_title(f"{next(letters)}  {gname}, {title}", loc="left", fontsize=7.0, fontweight="bold")
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        ax.tick_params(length=2, colors=INK)
        if row == 1:
            ax.set_xlabel("Δ per-disease AP vs no write-back", color=MUTED)
fig.legend(*axes[0][0].get_legend_handles_labels(), frameon=False, loc="lower center", ncol=2, bbox_to_anchor=(0.55, -0.02), fontsize=6.5)
fig.tight_layout(rect=(0, 0.04, 1, 1), w_pad=0.8, h_pad=1.2)
fig.savefig(f"{OUT}/Figure4_policies.pdf", bbox_inches="tight"); fig.savefig(f"{OUT}/Figure4_policies.png", dpi=600, bbox_inches="tight")
print("saved")
