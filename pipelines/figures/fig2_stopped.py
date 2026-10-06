#!/usr/bin/env python
"""Figure 2: share of stopped drug-disease pairs that are recorded or approved indications."""
import os, sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = sys.argv[1] if len(sys.argv) > 1 else "work/results"
OUT = sys.argv[2] if len(sys.argv) > 2 else "work/figures"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "pdf.fonttype": 42, "axes.linewidth": 0.6})
INK, MUTED, GRID, BLUE = "#0b0b0b", "#52514e", "#e6e5e1", "#2a78d6"
d = pd.read_csv(f"{RES}/D2_stopped_pairs_approval.tsv", sep="\t")
d = d[d.phases == "phase 1-3 only"]
order = ["efficacy", "safety", "design", "operational", "uninformative", "no reason given"]
labels = {"efficacy": "Efficacy (Negative)", "safety": "Safety or side effects", "design": "Study design",
          "operational": "Operational", "uninformative": "Uninformative", "no reason given": "No reason given"}
fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.4), sharey=True)
for ax, (level, title) in zip(axes, [("Hetionet", "a  Hetionet pairs (recorded treatment or approval)"),
                                     ("Open Targets", "b  All Open Targets pairs (approval)")]):
    sub = d[d.level == level].set_index("group").loc[order]
    y = range(len(order))[::-1]
    ax.barh(list(y), sub.share.values * 100, height=0.55, color=BLUE, edgecolor="none")
    for yi, (_, r) in zip(y, sub.iterrows()):
        ax.text(r.share * 100 + 0.8, yi, f"{r.share * 100:.0f}%  ({int(r.pairs_recorded_or_approved):,}/{int(r.pairs):,})",
                va="center", ha="left", fontsize=6.3, color=INK)
    ax.set_yticks(list(y)); ax.set_yticklabels([labels[g] for g in order])
    ax.set_xlim(0, 62 if level == "Hetionet" else 32)
    ax.set_xlabel("Stopped pairs that are recorded or approved indications (%)", color=MUTED)
    ax.set_title(title, loc="left", fontsize=7.5, fontweight="bold")
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    ax.tick_params(length=0, colors=INK)
fig.tight_layout(w_pad=2.5)
os.makedirs(OUT, exist_ok=True)
fig.savefig(f"{OUT}/Figure2_stopped_pairs.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/Figure2_stopped_pairs.png", dpi=600, bbox_inches="tight")
print("saved")
