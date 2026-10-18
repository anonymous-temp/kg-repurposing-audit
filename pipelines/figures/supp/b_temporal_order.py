#!/usr/bin/env python
"""FigS_temporal_order (release 0.7.0): when stopped trials of indications started, and the odds ratio of a scientific stop by
testing-intensity stratum (Hetionet; paper_results/robustness, written by pipelines/26_robustness_checks.py).
A  start year of each stopped phase 1-3 trial (before 2015) of an indication, relative to the earliest dated European (EMA) or
   Japanese (PMDA) approval for the same drug and disease; trial-pair rows, efficacy or safety stops shaded
B  odds ratio of an efficacy or safety stop for being an indication within six strata of registered trials (Haldane-Anscombe
   correction for display), with the stratified Mantel-Haenszel estimate and the estimate adjusted for log trial count"""
import json, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import b_common as C

S = C.S
S.setup()
ROB = os.environ.get("ROBUSTNESS", "paper_results/robustness")
rows = pd.read_csv(f"{ROB}/temporal_rows_hetionet.tsv", sep="\t")
T = json.load(open(f"{ROB}/temporal_order_hetionet.json"))["indications"]
M = json.load(open(f"{ROB}/mh_or_hetionet.json"))
ADJ = json.load(open(f"{C.SEL}/selection_hetionet.json"))["intensity"]["models"]["adjusted"]["terms"]["sci"]
assert abs(ADJ["OR"] - M["adjusted_or_log_trials"]) < 1e-6
d = rows.start_year - rows.dated_approval
assert len(d) == T["all_stops"]["dated"]["in_or_after"]["n"]
assert int((d >= 0).sum()) == T["all_stops"]["dated"]["in_or_after"]["k"]

fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.6), gridspec_kw=dict(left=0.08, right=0.98, top=0.84, bottom=0.2, wspace=0.42))
ax = axes[0]
bins = np.arange(-20.5, 21.5, 1)
dc = d.clip(-20, 20)
ax.hist(dc, bins=bins, color=S.LIGHT, edgecolor="white", lw=0.3, label="Other stops")
ax.hist(dc[rows.sci.astype(bool)], bins=bins, color=S.COL["Hetionet"], edgecolor="white", lw=0.3, label="Efficacy or safety stops")
ax.axvline(-0.5, ymax=0.74, color=S.MUTED, ls="--", lw=0.7)
a, s_ = T["all_stops"]["dated"]["in_or_after"], T["scientific"]["dated"]["in_or_after"]
ax.text(0.99, 0.97, f"Started in or after the approval year:\nall stops {100 * a['share']:.0f}% ({a['k']:,}/{a['n']:,})\n"
        f"efficacy or safety stops {100 * s_['share']:.0f}% ({s_['k']}/{s_['n']})", transform=ax.transAxes, ha="right", va="top",
        fontsize=6, linespacing=1.3)
ax.set_ylim(0, 1.55 * np.histogram(dc, bins=bins)[0].max())
ax.set_xlim(-20.5, 20.5)
ax.set_xlabel("Trial start year minus year of dated approval (EMA or PMDA)")
ax.set_ylabel("Stopped trial–pair rows")
ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.66), fontsize=6, handlelength=1.2)
S.panel(ax, "A", title="Stopped trials of indications relative to approval")

ax = axes[1]
st = M["strata"]
y = np.arange(len(st))[::-1]
for yi, r in zip(y, st):
    lo, hi = r["ci"]
    ax.errorbar(r["or_haldane"], yi, xerr=[[r["or_haldane"] - lo], [hi - r["or_haldane"]]], fmt="o", ms=3.4,
                color=S.COL["Hetionet"], elinewidth=0.8, capsize=1.6, capthick=0.8)
    n = r["stop_indication"] + r["stop_other"] + r["nostop_indication"] + r["nostop_other"]
    ax.text(1.02, yi, f"{n:,}", transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=6)
yp = -1.3
lo, hi = M["mh_ci"]
ax.errorbar(M["mh_or"], yp, xerr=[[M["mh_or"] - lo], [hi - M["mh_or"]]], fmt="D", ms=3.8, color=S.INK, elinewidth=0.9, capsize=1.8)
ax.errorbar(ADJ["OR"], yp - 0.8, xerr=[[ADJ["OR"] - ADJ["lo"]], [ADJ["hi"] - ADJ["OR"]]], fmt="s", ms=3.4, mfc="white", mec=S.INK,
            color=S.INK, elinewidth=0.9, capsize=1.8)
ax.axvline(1, color=S.MUTED, ls=":", lw=0.7)
ax.set_xscale("log")
ax.set_xlim(0.005, 4)
ax.set_xticks([0.01, 0.1, 0.5, 1, 2])
ax.set_xticklabels(["0.01", "0.1", "0.5", "1", "2"])
ax.set_yticks(list(y) + [yp, yp - 0.8])
lab = {"1": "1 trial", "2": "2 trials", "3-4": "3–4 trials", "5-9": "5–9 trials", "10-19": "10–19 trials", ">=20": "≥20 trials"}
ax.set_yticklabels([lab[r["stratum"]] for r in st] + ["Mantel–Haenszel", "Adjusted for log trials"])
ax.set_ylim(yp - 1.4, len(st) - 0.4)
ax.text(1.02, len(st) - 0.5, "Pairs", transform=ax.get_yaxis_transform(), va="bottom", ha="left", fontsize=6, color=S.MUTED)
ax.set_xlabel("Odds ratio of an efficacy or safety stop (log scale)")
S.panel(ax, "B", title="Within strata of registered trials")
S.save(fig, C.OUT, "FigS_temporal_order")
print("FigS_temporal_order written to", C.OUT)
