#!/usr/bin/env python
"""Figure: failed pairs in repoDB (output of pipelines/23_repodb_audit.py).

A  stated reasons of the failed trial records, typed by the Open Targets stop-reason classifier
B  share of failed pairs that are indications at three levels of disease matching
C  share of mapped failed pairs that are indications, by number of registered trials before 2017
D  the same share by phase of the failed trial and by stop reason
Usage: fig_repodb.py REPODB_AUDIT_JSON OUT_DIR [NAME]
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as S
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

AUD = sys.argv[1] if len(sys.argv) > 1 else "paper_results/repodb/repodb_audit.json"
OUT = sys.argv[2] if len(sys.argv) > 2 else "work/figures"
NAME = sys.argv[3] if len(sys.argv) > 3 else "Figure_repodb"
S.setup()
R = json.load(open(AUD))
G = "#1baf7a"          # repoDB series (third slot of the house palette)

fig = plt.figure(figsize=(7.0, 4.6))
gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.55, width_ratios=[1.05, 1])
axA, axB, axC, axD = (fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]))

# ---- A: stop reasons of failed records
T = R["records_by_stop_type"]
n = sum(T.values())
rows = [("Efficacy", "efficacy"), ("Safety", "safety"), ("Operational", "operational"), ("Design", "design"),
        ("Uninformative", "uninformative"), ("No reason given", "not_reported"), ("Text not in Open Targets", "unmatched")]
for i, (lab, k) in enumerate(rows):
    v = T.get(k, 0)
    col = G if k in ("efficacy", "safety") else (S.MUTED if k != "unmatched" else S.LIGHT)
    axA.barh(i, 100 * v / n, height=0.62, color=col, edgecolor="none")
    axA.text(100 * v / n + 1.2, i, f"{100 * v / n:.1f}% ({v:,})", va="center", fontsize=6.2, color=S.INK)
axA.set_yticks(range(len(rows))); axA.set_yticklabels([r[0] for r in rows]); axA.invert_yaxis()
axA.set_xlim(0, 75); axA.set_xlabel(f"Failed trial records (%; n = {n:,})")
axA.tick_params(axis="y", length=0)
S.panel(axA, "A", x=-0.42, title="Stated reason for stopping")


def dot(ax, x, s, filled=True, color=G, marker="o"):
    k = s["share"] * 100
    ax.errorbar(x, k, yerr=[[k - 100 * s["ci95"][0]], [100 * s["ci95"][1] - k]], fmt=marker, ms=3.6, color=color,
                mfc=color if filled else "white", mec=color, elinewidth=0.8, capsize=1.6, capthick=0.8)
    return k


# ---- B: share that are indications by level of matching
M = R["mapped_pairs"]
lev = [("Same UMLS\nconcept\n(all pairs)", R["exact_cui_all_pairs"]), ("Same disease\nconcept\n(mapped pairs)", M["same"]),
       ("Same concept or\nstage, line or\nsubtype (mapped)", M["same_or_broader"])]
for i, (lab, s) in enumerate(lev):
    k = dot(axB, i, s)
    axB.text(i, 100 * s["ci95"][1] + 1.0, f"{k:.1f}%\n{s['k']}/{s['n']:,}", va="bottom", ha="center", fontsize=6.0, color=S.INK)
axB.set_xticks(range(3)); axB.set_xticklabels([l for l, _ in lev], fontsize=6.0)
axB.set_xlim(-0.5, 2.5); axB.set_ylim(0, 30); axB.set_ylabel("Failed pairs that are indications (%)")
S.panel(axB, "B", title="Level of disease matching")

# ---- C: by number of registered trials
BINS = ["1", "2", "3-4", "5-9", "10-19", ">=20"]; BLAB = ["1", "2", "3–4", "5–9", "10–19", "≥20"]
for key, filled, off in [("by_trials_same", True, -0.09), ("by_trials", False, 0.09)]:
    xs, ys = [], []
    for i, b in enumerate(BINS):
        ys.append(dot(axC, i + off, R[key][b], filled)); xs.append(i + off)
    axC.plot(xs, ys, color=G, lw=0.7, ls="-" if filled else "--")
axC.set_xticks(range(len(BINS))); axC.set_xticklabels(BLAB); axC.set_ylim(0, 70)
axC.set_xlabel("Registered trials of the pair before 2017"); axC.set_ylabel("Failed pairs that are indications (%)")
hand = [Line2D([], [], marker="o", color=G, mfc=G, ls="-", ms=3.4, lw=0.7, label="Same concept"),
        Line2D([], [], marker="o", color=G, mfc="white", ls="--", ms=3.4, lw=0.7, label="Same concept or stage, line or subtype")]
axC.legend(handles=hand, loc="upper left", fontsize=6.0, handlelength=1.8)
S.panel(axC, "C", title="Testing intensity")

# ---- D: by phase and by stop reason
groups = [("Phase 1", R["by_phase_same"]["Phase 1"], R["by_phase"]["Phase 1"]),
          ("Phase 2", R["by_phase_same"]["Phase 2"], R["by_phase"]["Phase 2"]),
          ("Phase 3", R["by_phase_same"]["Phase 3"], R["by_phase"]["Phase 3"]),
          ("Efficacy\nor safety", R["by_reason_same"]["any_scientific_stop"], R["by_reason"]["any_scientific_stop"]),
          ("Operational\nonly", R["by_reason_same"]["only_operational_stops"], R["by_reason"]["only_operational_stops"])]
xpos = [0, 1, 2, 3.6, 4.9]
for x, (lab, s_same, s_broad) in zip(xpos, groups):
    dot(axD, x - 0.09, s_same, True); dot(axD, x + 0.09, s_broad, False)
axD.axvline(2.8, color=S.LIGHT, lw=0.6)
axD.text(1.0, 43, "Phase of the failed trial", ha="center", va="top", fontsize=6.2, color=S.MUTED)
axD.text(4.25, 43, "Stop reasons of the pair", ha="center", va="top", fontsize=6.2, color=S.MUTED)
axD.set_xticks(xpos); axD.set_xticklabels([g[0] for g in groups], fontsize=6.0)
axD.set_ylim(0, 45); axD.set_xlim(-0.5, 5.4); axD.set_ylabel("Failed pairs that are indications (%)")
S.panel(axD, "D", title="Phase and stop reason")
S.save(fig, OUT, NAME)
print("saved", NAME)
