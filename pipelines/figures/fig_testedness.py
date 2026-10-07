#!/usr/bin/env python
"""Figure 5: AUROC between groups of unlabelled pairs for scorers fitted without write-back (A, B), and per-disease AP on
external approved indications of these scorers and of scorers built only from registry history (C), both graphs."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as S
import matplotlib.pyplot as plt

RES = sys.argv[1] if len(sys.argv) > 1 else "work/results_testedness"
OUT = sys.argv[2] if len(sys.argv) > 2 else "work/figures"
SEL = sys.argv[3] if len(sys.argv) > 3 else "work/results_selection"
S.setup()
ROWS = [("A_vs_N", "Approved vs untested"), ("F_vs_N", "Later failure vs untested"), ("T_vs_N", "Newly tested vs untested"),
        ("A_vs_T", "Approved vs newly tested"), ("A_vs_F", "Approved vs later failure")]
fig = plt.figure(figsize=(7.0, 4.9)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.95], hspace=0.85, wspace=0.08)
axes = [fig.add_subplot(gs[0, 0])]; axes.append(fig.add_subplot(gs[0, 1], sharey=axes[0])); plt.setp(axes[1].get_yticklabels(), visible=False)
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
fig.legend(*axes[0].get_legend_handles_labels(), loc="upper center", ncol=3, bbox_to_anchor=(0.55, 0.505))
axC = fig.add_subplot(gs[1, :])
pos = axC.get_position(); axC.set_position([pos.x0 + 0.18, pos.y0 - 0.02, pos.width * 0.62, pos.height * 0.92])
E = {g: json.load(open(f"{SEL}/e2_scorers_{g}.json")) for g in ("hetionet", "primekg")}
H = {g: json.load(open(f"{SEL}/selection_{g}.json"))["history"] for g in ("hetionet", "primekg")}
CR = [("degree", "Degree reference", "m"), ("mf", "Label-only MF", "m"), ("graph", "Graph head", "m"),
      ("trial_count", "Number of trials before 2015", "h"), ("stopped_trial_count", "Number of stopped trials before 2015", "h"), ("any_stop", "Any stopped trial before 2015", "h")]
ys = list(range(len(CR)))[::-1]
for g, name, off in [("hetionet", "Hetionet", 0.14), ("primekg", "PrimeKG", -0.14)]:
    for y, (k, lab, kind) in zip(ys, CR):
        if kind == "m":
            v, (lo, hi) = E[g][k]["macroAP"], E[g][k]["ci"]
        else:
            v, (lo, hi) = H[g][k]["macroAP"], H[g][k]["macroAP_ci"]
        S.errorbar(axC, v, y + off, lo, hi, S.COL[name], S.MK[name], label=name if y == ys[0] else None, ms=3.2)
axC.axhline(2.5, color=S.LIGHT, lw=0.6)
axC.set_yticks(ys); axC.set_yticklabels([c[1] for c in CR]); axC.tick_params(axis="y", length=0)
axC.set_xlim(0, 0.36); axC.set_xlabel("Per-disease AP for external approved indications (95% CI)")
axC.legend(loc="lower right", fontsize=6.3)
S.panel(axC, "C", title="Scorers fitted without write-back and scorers using registry history only")
S.save(fig, OUT, "Figure5_tested_vs_approved")
print("saved")
