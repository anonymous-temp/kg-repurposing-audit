#!/usr/bin/env python
"""Figure: what no-write-back scorers rank highly - later-tested or later-approved pairs (both graphs)."""
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = sys.argv[1] if len(sys.argv) > 1 else "work/results_testedness"
OUT = sys.argv[2] if len(sys.argv) > 2 else "work/figures"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "pdf.fonttype": 42, "axes.linewidth": 0.6})
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
COL = {"graph": "#2a78d6", "mf": "#1baf7a", "degree": "#eda100"}       # categorical slots 1, 3, 4 (validated)
MK = {"graph": "o", "mf": "D", "degree": "^"}
LAB = {"graph": "Graph head", "mf": "Label-only MF", "degree": "Degree reference"}
ROWS = [("A_vs_N", "Approved vs never tested"), ("F_vs_N", "Later failure vs never tested"), ("T_vs_N", "Newly tested vs never tested"),
        ("A_vs_T", "Approved vs newly tested"), ("A_vs_F", "Approved vs later failure")]
graphs = [g for g in ["hetionet", "primekg"] if os.path.exists(f"{RES}/{g}.json")]
fig, axes = plt.subplots(1, len(graphs), figsize=(3.4 * len(graphs), 2.7), sharey=True, squeeze=False)
for ax, g in zip(axes[0], graphs):
    R = json.load(open(f"{RES}/{g}.json"))
    n = R["n"]
    ys = list(range(len(ROWS)))[::-1]
    for k, m in enumerate(["graph", "mf", "degree"]):
        off = (1 - k) * 0.22
        for y, (key, _) in zip(ys, ROWS):
            a = R["auroc"][f"{m}|{key}"]
            ax.plot([a["lo"], a["hi"]], [y + off] * 2, color=COL[m], lw=1.2, solid_capstyle="round")
            ax.plot(a["est"], y + off, MK[m], ms=4.2, color=COL[m], mec="white", mew=0.8, label=LAB[m] if y == ys[0] else None)
    ax.axvline(0.5, color=MUTED, lw=0.7)
    ax.axhline(1.5, color=MUTED, lw=0.5, ls=(0, (2, 2)))
    ax.set_yticks(ys); ax.set_yticklabels([r[1] for r in ROWS])
    ax.set_xlim(0.35, 0.85); ax.set_xlabel("AUROC (95% CI)", color=MUTED)
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    title = {"hetionet": "a  Hetionet", "primekg": "b  PrimeKG"}[g]
    ax.set_title(title, loc="left", fontsize=7.5, fontweight="bold", pad=10)
    ax.text(0.0, 1.005, f"n = {n['A']:,} approved, {n['F']:,} failed, {n['T']:,} newly tested", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=6.0, color=MUTED)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    ax.tick_params(length=2, colors=INK)
axes[0][0].legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5 if len(graphs) == 1 else 1.05, -0.24), ncol=3, fontsize=6.5)
fig.subplots_adjust(left=0.24, right=0.99, top=0.84, bottom=0.27, wspace=0.10)
fig.savefig(f"{OUT}/Figure3_tested_vs_approved.pdf", bbox_inches="tight"); fig.savefig(f"{OUT}/Figure3_tested_vs_approved.png", dpi=600, bbox_inches="tight")
print("saved", graphs)
