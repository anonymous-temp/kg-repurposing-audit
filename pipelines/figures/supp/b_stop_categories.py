#!/usr/bin/env python
"""FigS_stop_categories: Hetionet pairs with a stopped phase 1-3 trial, share that are recorded treatments or
approved indications (Table S6 and the sentence under it, supp_S3_S6.md).
A  by Open Targets stop-reason category (Wilson 95% intervals from the counts of Table S6)
B  by registry status (counts reconstructed from the rounded percentages given in the sentence under Table S6)"""
import re
import matplotlib.pyplot as plt
import b_common as C

S = C.S
S.setup()
PATH = f"{C.MS}/supp_S3_S6.md"
PRETTY = {"Insufficient_Enrollment": "Insufficient enrolment", "Business_Administrative": "Business or administrative",
          "No reason given": "No reason given", "Negative": "Negative (lack of efficacy)", "Logistics_Resources": "Logistics or resources",
          "Study_Design": "Study design", "Invalid_Reason": "Invalid reason", "Safety_Sideeffects": "Safety or side effects",
          "Another_Study": "Another study", "Study_Staff_Moved": "Study staff moved", "Regulatory": "Regulatory",
          "Uncategorised": "Uncategorised", "Covid19": "COVID-19", "No_Context": "No context"}

hdr, rows = C.read_table(PATH, r"Table S6\.")
cat = []
for r in rows:
    name, trials, pairs, k, share = r[0], C.num(r[1]), int(C.num(r[2])), int(C.num(r[3])), C.num(r[4])
    assert abs(100 * k / pairs - share) <= 0.51, (name, k, pairs, share)       # stated share agrees with the counts
    lo, hi = C.wilson(k, pairs)
    cat.append(dict(name=PRETTY[name], trials=int(trials), n=pairs, k=k, p=k / pairs, lo=lo, hi=hi))
cat.sort(key=lambda d: -d["p"])

sent = C.lines_after(PATH, r"Table S6\.", 1)[0]
stat = []
for m in re.finditer(r"(terminated|withdrawn|suspended) ([\d,]+) trials on ([\d,]+) pairs \((\d+)%", sent):
    n = int(C.num(m.group(3))); share = float(m.group(4)); k = round(share / 100 * n)
    lo, hi = C.wilson(k, n)
    stat.append(dict(name=m.group(1).capitalize(), trials=int(C.num(m.group(2))), n=n, k=k, p=k / n, lo=lo, hi=hi))
assert len(stat) == 3
stat.sort(key=lambda d: -d["p"])

fig = plt.figure(figsize=(6.7, 3.9))
gs = fig.add_gridspec(1, 1, left=0.34, right=0.86, top=0.93, bottom=0.12)
axA = fig.add_subplot(gs[0])
XMAX = 90


def dots(ax, data, title, letter, show_trials=False):
    ys = list(range(len(data)))[::-1]
    for y, d in zip(ys, data):
        ax.hlines(y, 0, XMAX, color=S.LIGHT, lw=0.5, zorder=0)
        S.errorbar(ax, 100 * d["p"], y, 100 * d["lo"], 100 * d["hi"], S.COL["Hetionet"], "o", ms=4.0)
        ax.text(1.02, y, f"{100 * d['p']:.0f}%  ({d['k']:,})", transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=6.3)
    ax.text(1.02, len(data) - 0.35, "Share (n recorded\nor approved)", transform=ax.get_yaxis_transform(), va="bottom", ha="left",
            fontsize=6, color=S.MUTED)
    labs = [f"{d['name']}  (n = {d['n']:,})" for d in data]
    ax.set_yticks(ys); ax.set_yticklabels(labs); ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, XMAX); ax.set_ylim(-0.6, len(data) - 0.4)
    ax.set_xticks(range(0, 81, 20))
    ax.set_title(title, loc="left", x=-0.50, fontsize=7)
    ax.spines["left"].set_visible(False)


dots(axA, cat, "Open Targets stop-reason category", "A")
axA.set_xlabel("Recorded treatments or approved indications (%, Wilson 95% interval)")
S.save(fig, C.OUT, "FigS_stop_categories")
print("saved", [(d["name"], d["k"], d["n"]) for d in stat])
