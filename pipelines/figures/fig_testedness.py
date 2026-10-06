#!/usr/bin/env python
"""Figure 4: AUROC between groups of unlabelled pairs for scorers fitted without write-back (both graphs)."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as S
import matplotlib.pyplot as plt

RES = sys.argv[1] if len(sys.argv) > 1 else "work/results_testedness"
OUT = sys.argv[2] if len(sys.argv) > 2 else "work/figures"
S.setup()
ROWS = [("A_vs_N", "Approved vs untested"), ("F_vs_N", "Later failure vs untested"), ("T_vs_N", "Newly tested vs untested"),
        ("A_vs_T", "Approved vs newly tested"), ("A_vs_F", "Approved vs later failure")]
fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6), sharey=True)
for ax, g, letter, name in [(axes[0], "hetionet", "A", "Hetionet"), (axes[1], "primekg", "B", "PrimeKG")]:
    R = json.load(open(f"{RES}/{g}.json")); n = R["n"]
    ys = list(range(len(ROWS)))[::-1]
    for k, m in enumerate(["graph", "mf", "degree"]):
        for y, (key, _) in zip(ys, ROWS):
            a = R["auroc"][f"{m}|{key}"]
            S.errorbar(ax, a["est"], y + (1 - k) * 0.22, a["lo"], a["hi"], S.COL[m], S.MK[m], label=S.LAB[m] if y == ys[0] else None)
    ax.axvline(0.5, color=S.MUTED, lw=0.6, ls=(0, (3, 2)))
    ax.axhline(1.5, color=S.LIGHT, lw=0.6)
    ax.set_yticks(ys); ax.set_yticklabels([r[1] for r in ROWS]); ax.tick_params(axis="y", length=0)
    ax.set_xlim(0.35, 0.85); ax.set_xlabel("AUROC (95% CI)")
    S.panel(ax, letter, title=f"{name}: {n['A']:,} approved, {n['F']:,} failed, {n['T']:,} newly tested")
fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", ncol=3, bbox_to_anchor=(0.55, -0.01))
fig.tight_layout(rect=(0, 0.07, 1, 1), w_pad=1.2)
S.save(fig, OUT, "Figure4_tested_vs_approved")
print("saved")
