#!/usr/bin/env python
"""Figure 2: share of stopped drug-disease pairs that are recorded treatments or approved indications,
by failure type (Hetionet, PrimeKG, all Open Targets pairs) and by the role of the drug in scientific stops."""
import json, os, sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = sys.argv[1] if len(sys.argv) > 1 else "work/results"
PKG = sys.argv[2] if len(sys.argv) > 2 else "work/results_primekg"
ROLE = sys.argv[3] if len(sys.argv) > 3 else "work/results_role"
OUT = sys.argv[4] if len(sys.argv) > 4 else "work/figures"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "pdf.fonttype": 42, "axes.linewidth": 0.6})
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
C1, C2 = "#2a78d6", "#eb6834"                       # categorical slots 1-2 (validated)
order = ["efficacy", "safety", "design", "operational", "uninformative", "no reason given"]
labels = {"efficacy": "Efficacy (Negative)", "safety": "Safety or side effects", "design": "Study design",
          "operational": "Operational", "uninformative": "Uninformative", "no reason given": "No reason given"}
d = pd.read_csv(f"{RES}/D2_stopped_pairs_approval.tsv", sep="\t"); d = d[d.phases == "phase 1-3 only"]
p = pd.read_csv(f"{PKG}/D2_stopped_pairs_approval.tsv", sep="\t")
rh = json.load(open(f"{ROLE}/role_descriptive.json"))["wb_scientific"]; rp = json.load(open(f"{PKG}/descriptive.json"))["roles_scientific"]


def style(ax, title, xmax, xlabel):
    ax.set_xlim(0, xmax); ax.set_xlabel(xlabel, color=MUTED)
    ax.set_title(title, loc="left", fontsize=7.5, fontweight="bold")
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    ax.tick_params(length=0, colors=INK)


fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.4))
X = "Stopped pairs that are recorded or approved indications (%)"
for ax, (sub, title, xmax) in zip(axes.flat[:3], [(d[d.level == "Hetionet"], "a  Hetionet pairs", 62),
                                                   (p, "b  PrimeKG pairs", 62),
                                                   (d[d.level == "Open Targets"], "c  All Open Targets pairs (approval only)", 32)]):
    sub = sub.set_index("group").loc[order]; y = list(range(len(order)))[::-1]
    ax.barh(y, sub.share.values * 100, height=0.55, color=C1, edgecolor="none")
    for yi, (_, r) in zip(y, sub.iterrows()):
        ax.text(r.share * 100 + 0.8, yi, f"{r.share * 100:.0f}%  ({int(r.pairs_recorded_or_approved):,}/{int(r.pairs):,})", va="center", ha="left", fontsize=6.2, color=INK)
    ax.set_yticks(y); ax.set_yticklabels([labels[g] for g in order]); style(ax, title, xmax, X if ax is axes.flat[2] else "")
ax = axes.flat[3]
cats = [("any investigational", "Investigational agent"), ("comparator/backbone only", "Comparator or background\ntherapy only")]
other_h = [k for k in rh if k.startswith("other")][0]; other_p = [k for k in rp if k.startswith("other")][0]
cats.append(("__other__", "Other or unresolved\nroles only"))
y = [2, 1, 0]
for k, (key, lab) in enumerate(cats):
    for j, (R, okey, col, name) in enumerate([(rh, other_h, C1, "Hetionet"), (rp, other_p, C2, "PrimeKG")]):
        v = R[okey if key == "__other__" else key]
        yy = y[k] + (0.17 if j == 0 else -0.17)
        ax.barh(yy, v["pct"], height=0.32, color=col, edgecolor="white", linewidth=0.8, label=name if k == 0 else None)
        ax.text(v["pct"] + 0.8, yy, f"{v['pct']:.0f}%  ({v['recorded_or_approved']}/{v['pairs']})", va="center", ha="left", fontsize=6.2, color=INK)
ax.set_yticks(y); ax.set_yticklabels([c[1] for c in cats])
style(ax, "d  Scientific stops before 2015, by drug role", 100, X); ax.set_xticks([0, 20, 40, 60])
ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 1.0), fontsize=6.3, handlelength=1.2)
fig.tight_layout(w_pad=2.0, h_pad=1.6)
os.makedirs(OUT, exist_ok=True)
fig.savefig(f"{OUT}/Figure2_stopped_pairs.pdf", bbox_inches="tight"); fig.savefig(f"{OUT}/Figure2_stopped_pairs.png", dpi=600, bbox_inches="tight")
print("saved")
