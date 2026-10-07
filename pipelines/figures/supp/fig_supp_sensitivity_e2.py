#!/usr/bin/env python3
"""FigS_sensitivity_e2: Hetionet graph head on external approved indications (E2) and later scientific failures (E3) under the main
analysis (write-back cut-off 2015) and two sensitivity analyses (cut-off 2017; phase 4 trials excluded). Table S12 (supp_S7_S14.md)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import supp_common as C
S = C.S
S.setup()

_, rows = C.table_dicts(C.F_S7S14, "Table S12.")
# (label in the table, label in the legend, marker)
ANALYSES = [("Main analysis (cut-off 2015)", "Main analysis (cut-off 2015)", "o"),
            ("Cut-off 2017", "Write-back cut-off 2017", "v"),
            ("Phase 1–3 trials only", "Phase 4 trials excluded", "P")]
POLS = ["no_writeback", "flat_negative", "typed", "typed_scoped", "mask_all"]
D = {}
for r in rows:
    D[(r["Analysis"], C.policy_key(r["Policy"]))] = {
        "pooled": C.num(r["Pooled AP"]), "perdis": C.num(r["Per-disease AP"]), "pct_sci": C.num(r["Median percentile, scientific stratum"]),
        "pct_other": C.num(r["Median percentile, other-stop stratum"]), "e3": C.num(r["E3 AUROC"])}
assert all((a, p) in D for a, _, _ in ANALYSES for p in POLS), "missing analysis/policy cell"
assert {k[0] for k in D} == {a for a, _, _ in ANALYSES}, {k[0] for k in D}

PANELS = [("perdis", "E2 per-disease AP"), ("pooled", "E2 pooled AP"), ("e3", "E3 AUROC"),
          ("pct_sci", "Median percentile,\nscientific-stop stratum"), ("pct_other", "Median percentile,\nother-stop stratum")]
fig, axes = plt.subplots(2, 3, figsize=(6.7, 4.5), sharey=True)
ys = np.arange(len(POLS))[::-1]
offs = [0.29, 0.0, -0.29]
for k, (key, lab) in enumerate(PANELS):
    ax = axes[k // 3][k % 3]
    for (a, _, mk), off in zip(ANALYSES, offs):
        xs = [D[(a, p)][key] for p in POLS]
        ax.plot(xs, ys + off, mk, color=S.COL["graph"], mfc=S.COL["graph"], mec=S.COL["graph"], ms=3.4 if mk != "P" else 3.8, ls="none")
    for y in ys[:-1]:
        ax.axhline(y - 0.5, color=S.LIGHT, lw=0.4, zorder=0)
    ax.set_ylim(-0.6, len(POLS) - 0.4)
    ax.set_yticks(ys)
    ax.set_yticklabels([C.PL[p] for p in POLS])      # inner panels are hidden by sharey
    ax.tick_params(axis="y", length=0)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.set_xlabel(lab)
    S.panel(ax, "ABCDE"[k], x=-0.02, y=1.03)
axF = axes[1][2]
axF.axis("off")
hand = [Line2D([], [], marker=mk, color=S.COL["graph"], mfc=S.COL["graph"], ls="none", ms=4, label=lab) for _, lab, mk in ANALYSES]
axF.legend(handles=hand, loc="center left", bbox_to_anchor=(0.0, 0.62), title="Hetionet, graph head", title_fontsize=6.5, alignment="left")
fig.tight_layout(w_pad=1.0, h_pad=1.4)
S.save(fig, C.OUT, "FigS_sensitivity_e2")
print("saved FigS_sensitivity_e2")
