#!/usr/bin/env python
"""FigS_permutation_null: degree-preserving permutation null for external approved indications among unlabelled stopped pairs
(Table S20e in supp_S20_S21.md, cross-checked against selection_{hetionet,primekg}.json).
The JSON files hold only summaries of the 1,000 permutations (mean, 2.5th/97.5th percentile, ratio interval, one-sided P), not the
per-permutation counts, so no histograms can be drawn.
A, B  observed number of approved indications vs null mean with 95% range (log axis), Hetionet and PrimeKG
C     observed / expected ratio with its 95% range (infinite upper limits drawn as arrows)"""
import json, re, math
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import b_common as C

S = C.S
S.setup()
SEL = C.SEL
G = ("Hetionet", "PrimeKG")
KEY = {"All stopped pairs": "wb_all", "Scientific stop": "wb_scientific", "… same concept": "wb_scientific_scoped",
       "… investigational agent": "wb_scientific_inv", "… both": "wb_scientific_inv_scoped", "Other stops only": "wb_other"}
LABEL = {"All stopped pairs": "All stopped pairs", "Scientific stop": "Efficacy or safety stop", "… same concept": "\u2026 same disease concept",
         "… investigational agent": "\u2026 investigational drug", "… both": "\u2026 both", "Other stops only": "Other stops only"}
J = {g: json.load(open(f"{SEL}/selection_{g.lower()}.json"))["null"] for g in G}

_, rows = C.read_table(f"{C.MS}/supp_S20_S21.md", r"Table S20e\.")
R = {}
order = []
for r in rows:
    s, g = r[0], r[1]
    exp = re.match(r"([\d.]+) \(([\d.]+)–([\d.]+)\)", r[4])
    rat = re.match(r"([\d.]+) \(([\d.]+)–(\S+)\)", r[5])
    d = dict(pairs=int(C.num(r[2])), obs=int(C.num(r[3])), em=float(exp.group(1)), elo=float(exp.group(2)), ehi=float(exp.group(3)),
             rm=float(rat.group(1)), rlo=float(rat.group(2)), rhi=math.inf if rat.group(3) == "∞" else float(rat.group(3)))
    j = J[g][KEY[s]]                                                   # cross-check against the JSON summaries
    assert j["pairs"] == d["pairs"] and j["observed"] == d["obs"], (s, g)
    assert abs(j["null_mean"] - d["em"]) < 0.051 and abs(j["null_lo"] - d["elo"]) <= 0.5 and abs(j["null_hi"] - d["ehi"]) <= 0.5, (s, g)
    assert abs(j["ratio_observed_to_null"] - d["rm"]) < 0.051, (s, g)
    R[(s, g)] = d
    if s not in order:
        order.append(s)
assert len(order) == 6 and len(R) == 12
P = {g: [J[g][KEY[s]]["p_one_sided"] for s in order] for g in G}

fig, axes = plt.subplots(1, 3, figsize=(6.7, 3.15), gridspec_kw=dict(left=0.245, right=0.985, top=0.80, bottom=0.17, wspace=0.50))
ys = list(range(len(order)))[::-1]
XMIN = 0.5
for ax, g, letter in ((axes[0], "Hetionet", "A"), (axes[1], "PrimeKG", "B")):
    for y, s in zip(ys, order):
        d = R[(s, g)]
        ax.hlines(y, XMIN, 1000, color=S.LIGHT, lw=0.4, zorder=0)
        lo = max(d["elo"], XMIN)
        ax.errorbar(d["em"], y - 0.2, xerr=[[d["em"] - lo], [d["ehi"] - d["em"]]], fmt="D", ms=3.0, mfc="white", mec=S.MUTED, color=S.MUTED,
                    elinewidth=0.7, capsize=1.4, capthick=0.7)
        if d["elo"] < XMIN:
            ax.plot(XMIN, y - 0.2, "<", ms=3.0, color=S.MUTED, clip_on=False)
        ax.plot(d["obs"], y + 0.2, S.MK[g], ms=3.8, color=S.COL[g])
        ax.text(1.04, y, f"{d['pairs']:,}", transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=6)
    ax.text(1.04, len(order) - 0.4, "Pairs", transform=ax.get_yaxis_transform(), va="bottom", ha="left", fontsize=6, color=S.MUTED)
    ax.set_xscale("log"); ax.set_xlim(XMIN, 500)
    ax.set_xticks([1, 10, 100]); ax.set_xticklabels(["1", "10", "100"])
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel("Approved indications (count)")
    S.panel(ax, letter, x=-0.02, title=g)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_yticks(ys)
    ax.set_yticklabels([LABEL[s] for s in order] if letter == "A" else [])

ax = axes[2]
XM = 50
for g, off in (("Hetionet", 0.14), ("PrimeKG", -0.14)):
    for y, s in zip(ys, order):
        d = R[(s, g)]
        hi = min(d["rhi"], XM)
        S.errorbar(ax, d["rm"], y + off, d["rlo"], hi, S.COL[g], S.MK[g], ms=3.4)
        if math.isinf(d["rhi"]):
            ax.plot(XM, y + off, ">", ms=3.0, color=S.COL[g], clip_on=False)
for y in ys:
    ax.hlines(y, 1, XM, color=S.LIGHT, lw=0.4, zorder=0)
ax.axvline(1, color=S.MUTED, lw=0.6, ls=":")
ax.set_xscale("log"); ax.set_xlim(0.9, XM)
ax.set_xticks([1, 2, 5, 10, 20, 50]); ax.set_xticklabels(["1", "2", "5", "10", "20", "50"])
ax.set_ylim(-0.6, len(order) - 0.4); ax.set_yticks(ys); ax.set_yticklabels([]); ax.tick_params(axis="y", length=0)
ax.spines["left"].set_visible(False)
ax.set_xlabel("Observed / expected (log scale)")
S.panel(ax, "C", x=-0.02, title="Ratio and 95% range")

hand = [Line2D([], [], marker=S.MK[g], ls="", color=S.COL[g], ms=3.8, label=f"Observed, {g}") for g in G]
hand.append(Line2D([], [], marker="D", ls="", mfc="white", mec=S.MUTED, ms=3.0, label="Null: mean and 95% range of 1,000 permutations"))
hand.append(Line2D([], [], marker=">", ls="", color=S.MUTED, ms=3.0, label="Limit beyond the axis (lower limit 0 or upper limit \u221e)"))
fig.legend(handles=hand, loc="upper center", ncol=2, bbox_to_anchor=(0.56, 1.02), handletextpad=0.3, columnspacing=1.6, fontsize=6.2)
S.save(fig, C.OUT, "FigS_permutation_null")
print("saved; one-sided P:", {g: [round(p, 4) for p in P[g]] for g in G})
