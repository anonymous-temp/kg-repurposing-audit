#!/usr/bin/env python
"""Figure 4: testing intensity and the meaning of a stop (both graphs).

A, B  share of tested pairs that are recorded or approved indications, by number of registered phase 1-3 trials
      started before 2015, with and without an efficacy or safety stop (Wilson 95% intervals)
C     share of indications among written-back pairs with an efficacy or safety stop, by number of trials
D     odds ratio of an efficacy or safety stop for being an indication, crude and adjusted (95% disease-cluster
      bootstrap intervals, 1,000 draws)
E     share of indications among pairs with an efficacy or safety stop, by phase of the stopped trial
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as S
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SEL = sys.argv[1] if len(sys.argv) > 1 else "work/results_selection"
OUT = sys.argv[2] if len(sys.argv) > 2 else "work/figures"
S.setup()
R = {g: json.load(open(f"{SEL}/selection_{g.lower()}.json")) for g in ("Hetionet", "PrimeKG")}
BINS = ["1", "2", "3-4", "5-9", "10-19", "20+"]; BLAB = ["1", "2", "3–4", "5–9", "10–19", "≥20"]

fig = plt.figure(figsize=(7.0, 4.7))
gs = fig.add_gridspec(2, 6, height_ratios=[1, 1], hspace=0.75, wspace=1.6)
axA, axB, axD = fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6])
axC, axE = fig.add_subplot(gs[1, 2:4]), fig.add_subplot(gs[1, 4:6])


def vbar(ax, x, k, ci, color, marker, filled=True):
    mfc = color if filled else "white"
    ax.errorbar(x, 100 * k, yerr=[[100 * (k - ci[0])], [100 * (ci[1] - k)]], fmt=marker, ms=3.4, color=color, mfc=mfc, mec=color,
                elinewidth=0.8, capsize=1.5, capthick=0.8)


for ax, g, letter in [(axA, "Hetionet", "A"), (axB, "PrimeKG", "B")]:
    rows = R[g]["intensity"]["by_trial_count"]
    for s, off, filled in [(0, -0.12, False), (1, 0.12, True)]:
        xs, ys = [], []
        for i, b in enumerate(BINS):
            r = next(x for x in rows if x["bin"] == b and x["scientific_stop"] == s)
            vbar(ax, i + off, r["share"], r["ci"], S.COL[g], S.MK[g], filled); xs.append(i + off); ys.append(100 * r["share"])
        ax.plot(xs, ys, color=S.COL[g], lw=0.7, ls="-" if filled else "--")
    m = R[g]["intensity"]["models"]
    ax.text(0.02, 0.97, f"Crude OR {m['crude']['terms']['sci']['OR']:.2f}\nAdjusted OR {m['adjusted']['terms']['sci']['OR']:.2f}",
            transform=ax.transAxes, va="top", fontsize=6.3, color=S.INK)
    ax.set_xticks(range(len(BINS))); ax.set_xticklabels(BLAB); ax.set_ylim(0, 100)
    ax.set_xlabel("Registered trials before 2015"); ax.set_ylabel("Indications (%)")
    S.panel(ax, letter, title=g)
hand = [Line2D([], [], marker="o", color="k", mfc="k", ls="-", ms=3.2, lw=0.7, label="Efficacy or safety stop"),
        Line2D([], [], marker="o", color="k", mfc="white", ls="--", ms=3.2, lw=0.7, label="No such stop")]
axA.legend(handles=hand, loc="upper left", bbox_to_anchor=(0.0, 0.78), handlelength=1.6, fontsize=6.2)

SPECS = [("crude", "sci", "Crude"), ("adjusted", "sci", "Adjusted for number of trials"),
         ("adjusted_popularity", "sci", "+ trials of the drug and the disease"), ("adjusted_phase12_intensity", "sci", "Trials counted in phase 1–2 only"),
         ("adjusted_oncology", "sci", "Cancers only"), ("adjusted_non_oncology", "sci", "Other diseases only")]
ys = np.arange(len(SPECS))[::-1]
for g, off in [("Hetionet", 0.14), ("PrimeKG", -0.14)]:
    for y, (k, t, _) in zip(ys, SPECS):
        e = R[g]["intensity"]["models"][k]["terms"][t]
        S.errorbar(axC, e["OR"], y + off, e["lo"], e["hi"], S.COL[g], S.MK[g], ms=3.2)
axC.axvline(1, color=S.MUTED, lw=0.6, ls=":")
axC.set_xscale("log"); axC.set_xlim(0.15, 7); axC.set_xticks([0.2, 0.5, 1, 2, 5]); axC.set_xticklabels(["0.2", "0.5", "1", "2", "5"])
axC.set_yticks(ys); axC.set_yticklabels([s[2] for s in SPECS]); axC.tick_params(axis="y", length=0)
axC.set_xlabel("Odds ratio for being an indication (log scale)")
S.panel(axC, "D", x=-1.25, title="Odds ratio of an efficacy or safety stop")

for g, off in [("Hetionet", -0.1), ("PrimeKG", 0.1)]:
    rows = R[g]["purity"]["wb_scientific"]["by_trial_count"]; xs, yv = [], []
    for i, b in enumerate(BINS):
        r = next(x for x in rows if x["bin"] == b)
        vbar(axD, i + off, r["share"], r["ci"], S.COL[g], S.MK[g]); xs.append(i + off); yv.append(100 * r["share"])
    axD.plot(xs, yv, color=S.COL[g], lw=0.7)
axD.set_xticks(range(len(BINS))); axD.set_xticklabels(BLAB); axD.set_ylim(0, 100)
axD.set_xlabel("Registered trials before 2015"); axD.set_ylabel("Indications (%)")
S.panel(axD, "C", title="Written-back pairs, efficacy or safety stop")

PH = ["1", "2", "3", "4"]
for g, off in [("Hetionet", -0.08), ("PrimeKG", 0.08)]:
    rows = R[g]["purity"]["scientific_by_phase"]; xs, yv = [], []
    for i, ph in enumerate(PH):
        r = next(x for x in rows if x["phase"] == ph)
        vbar(axE, i + off, r["share"], r["ci"], S.COL[g], S.MK[g]); xs.append(i + off); yv.append(100 * r["share"])
    axE.plot(xs, yv, color=S.COL[g], lw=0.7)
axE.set_xticks(range(4)); axE.set_xticklabels(PH); axE.set_ylim(0, 100); axE.set_xlim(-0.5, 3.5)
axE.set_xlabel("Phase of the stopped trial"); axE.set_ylabel("Indications (%)")
S.panel(axE, "E", title="Pairs with an efficacy or safety stop")
hand = [Line2D([], [], marker=S.MK[g], color=S.COL[g], ls="-", ms=3.4, lw=0.7, label=g) for g in ("Hetionet", "PrimeKG")]
fig.legend(handles=hand, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.02))
S.save(fig, OUT, "Figure4_testing_intensity")
print("saved")
