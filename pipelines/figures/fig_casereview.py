#!/usr/bin/env python
"""Figure 5: why reviewed indications acquired their cleanest failure record (80 cases, both graphs).

One square per reviewed indication, in rows by primary code; fill gives the graph, a black outline marks the cases
whose failed trial plausibly contradicts the indication (solid) or might (dashed). An example per code is printed
on the right. Input: case_taxonomy_verified.tsv (paper_results/case_review)."""
import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as S
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

TSV = sys.argv[1] if len(sys.argv) > 1 else "paper_results/case_review/case_taxonomy_verified.tsv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "work/figures"
S.setup()
d = pd.read_csv(TSV, sep="\t")
LAB = {"E": "Stop not attributable to the drug", "B2": "New combination or add-on", "A2": "Other stage, line or goal",
       "A1": "Special population", "C": "Active comparator", "D": "Narrow endpoint or subtype", "F": "Other"}
EX = {"E": "Sacubitril–valsartan vs enalapril, stopped early for benefit",
      "B2": "Lapatinib in phase 1 combinations, stopped for tolerability",
      "A2": "Transdermal estradiol in prostate cancer after docetaxel",
      "A1": "Biweekly docetaxel in women aged ≥70 years with breast cancer",
      "C": "First-line trabectedin vs doxorubicin in soft-tissue sarcoma",
      "D": "Fingolimod in primary progressive multiple sclerosis",
      "F": "Olanzapine with fluoxetine in treatment-resistant depression"}
order = d.primary.value_counts().sort_values(ascending=False, kind="stable").index.tolist()
order = sorted(order, key=lambda c: (-int((d.primary == c).sum()), list(LAB).index(c)))
fig, ax = plt.subplots(figsize=(6.9, 2.55))
W, H = 0.82, 0.72
for i, code in enumerate(order):
    y = len(order) - 1 - i
    sub = d[d.primary == code].assign(g=lambda x: x.graph.map({"hetionet": 0, "primekg": 1})).sort_values(["g", "case_id"])
    for j, (_, r) in enumerate(sub.iterrows()):
        col = S.COL["Hetionet"] if r.graph == "hetionet" else S.COL["PrimeKG"]
        ax.add_patch(Rectangle((j, y - H / 2), W, H, facecolor=col, edgecolor="none"))
        if r.contradicts in ("yes", "unclear"):
            ax.add_patch(Rectangle((j - 0.06, y - H / 2 - 0.06), W + 0.12, H + 0.12, facecolor="none", edgecolor=S.INK,
                                   lw=1.1, ls="-" if r.contradicts == "yes" else (0, (1.5, 1.0))))
    ax.text(-0.6, y, f"{code}  {LAB[code]}", ha="right", va="center", fontsize=6.6)
    ax.text(len(sub) + 0.4, y, f"{len(sub)}", ha="left", va="center", fontsize=6.6, fontweight="bold")
    ax.text(26.2, y, EX[code], ha="left", va="center", fontsize=6.2, color=S.MUTED, style="italic")
ax.text(26.2, len(order) - 0.25, "Example", ha="left", va="bottom", fontsize=6.6, color=S.INK)
ax.set_xlim(-0.3, 47); ax.set_ylim(-0.7, len(order) - 0.2)
ax.axis("off")
handles = [Rectangle((0, 0), 1, 1, facecolor=S.COL["Hetionet"]), Rectangle((0, 0), 1, 1, facecolor=S.COL["PrimeKG"]),
           Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=S.INK, lw=1.1),
           Rectangle((0, 0), 1, 1, facecolor="white", edgecolor=S.INK, lw=1.1, ls=(0, (1.5, 1.0)))]
n_yes = int((d.contradicts == "yes").sum()); n_unc = int((d.contradicts == "unclear").sum())
ax.legend(handles, ["Hetionet (40 cases)", "PrimeKG (40 cases)", f"Trial contradicts the indication ({n_yes})", f"Unclear ({n_unc})"],
          loc="upper center", bbox_to_anchor=(0.42, -0.02), ncol=4, handlelength=1.0, handleheight=1.0, columnspacing=1.6)
S.save(fig, OUT, "Figure5_case_review")
print(d.primary.value_counts().to_dict(), d.contradicts.value_counts().to_dict())
