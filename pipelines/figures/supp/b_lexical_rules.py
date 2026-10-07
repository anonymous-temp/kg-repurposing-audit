#!/usr/bin/env python
"""FigS_lexical_rules: lexical stop-reason rules of software release 0.2.0 (Table S19 in supp_S16_S18.md; error examples in S12 are qualitative).
A  heatmap of precision, recall and F1 against the expert labels, and of Cohen's kappa against the Open Targets classifier
B  expert-labelled positives vs rule positives (log axis) for each rule"""
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import b_common as C

S = C.S
S.setup()
PATH = f"{C.MS}/supp_S16_S18.md"
hdr, rows = C.read_table(PATH, r"Table S19\.")
print(hdr)
R = [dict(rule=r[0], gold=int(C.num(r[1])), rulepos=int(C.num(r[2])), prec=C.num(r[3]), rec=C.num(r[4]), f1=C.num(r[5]), kappa=C.num(r[6])) for r in rows]
assert [d["rule"] for d in R] == ["Efficacy", "Safety", "Operational"]
for d in R:                                                    # F1 consistency check with the rounded precision and recall
    assert abs(2 * d["prec"] * d["rec"] / (d["prec"] + d["rec"]) - d["f1"]) < 0.02, d
note = C.lines_after(PATH, r"Table S19\.", 1)[0]
NK = int(C.num(re.search(r"computed on ([\d,]+) stopped trials", note).group(1)))

cmap = LinearSegmentedColormap.from_list("seq", ["#ffffff", "#bfe9d9", "#1baf7a", "#0a5c40"])
fig = plt.figure(figsize=(6.7, 2.55))
W, H = 6.7, 2.55


def axin(l, b, w, h):
    return fig.add_axes([l / W, b / H, w / W, h / H])


cw, ch = 0.62, 0.55
top = H - 0.75
axA = axin(1.15, top - 3 * ch, 3 * cw, 3 * ch)
axK = axin(1.15 + 3 * cw + 0.12, top - 3 * ch, cw, 3 * ch)
yl = [f"{d['rule']} rule" for d in R]


def heat(ax, cols, labels):
    M = np.array([[d[c] for c in cols] for d in R])
    ax.imshow(M, cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="nearest")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=7, color="white" if M[i, j] > 0.62 else S.INK)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(labels); ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    for k in range(1, M.shape[1]):
        ax.axvline(k - 0.5, color="white", lw=1.2)
    for k in range(1, M.shape[0]):
        ax.axhline(k - 0.5, color="white", lw=1.2)
    ax.xaxis.tick_top()


heat(axA, ["prec", "rec", "f1"], ["Precision", "Recall", "F1"])
heat(axK, ["kappa"], ["Cohen’s κ"])
axA.set_yticks(range(3)); axA.set_yticklabels(yl); axK.set_yticks([])
fig.text(0.1 / W, (top + 0.62) / H, "A", fontsize=9, fontweight="bold", va="bottom")
fig.text(0.3 / W, (top + 0.62) / H, "Agreement of the rules", fontsize=7, va="bottom")
axA.text(0.5, 1.17, "vs expert labels", transform=axA.transAxes, ha="center", fontsize=6.5, color=S.MUTED, va="bottom")
axK.text(0.5, 1.17, "vs Open Targets\nclassifier", transform=axK.transAxes, ha="center", fontsize=6.5, color=S.MUTED, va="bottom", linespacing=1.0)
axA.text(0.5, -0.12, "Scale 0–1; κ on %s stopped trials with a stated reason" % f"{NK:,}", transform=axA.transAxes, ha="left", fontsize=5.8, color=S.MUTED, va="top") if False else None

axB = axin(4.75, top - 3 * ch, 1.75, 3 * ch)
ys = [2, 1, 0]
for y, d in zip(ys, R):
    axB.plot([d["rulepos"], d["gold"]], [y, y], color=S.LIGHT, lw=1.6, zorder=1)
    axB.plot(d["gold"], y, "o", ms=4.6, color=S.MUTED, zorder=3)
    axB.plot(d["rulepos"], y, "o", ms=4.6, mfc="white", mec=S.MUTED, mew=1.0, zorder=3)
    for v, dy in ((d["gold"], 0.30), (d["rulepos"], 0.30)):
        pass
    axB.text(d["gold"] * 1.12 if d["gold"] > d["rulepos"] else d["gold"] / 1.12, y + 0.28, f"{d['gold']:,}", ha="left" if d["gold"] > d["rulepos"] else "right",
             va="center", fontsize=6)
    axB.text(d["rulepos"] * 1.12 if d["rulepos"] > d["gold"] else d["rulepos"] / 1.12, y + 0.28, f"{d['rulepos']:,}",
             ha="left" if d["rulepos"] > d["gold"] else "right", va="center", fontsize=6)
axB.set_xscale("log"); axB.set_xlim(60, 6000)
axB.set_yticks(ys); axB.set_yticklabels([]); axB.tick_params(axis="y", length=0)
axB.set_ylim(-0.5, 2.5)
axB.set_xticks([100, 1000]); axB.set_xticklabels(["100", "1,000"])
axB.set_xlabel("Texts labelled positive (log scale)")
axB.spines["left"].set_visible(False)
axB.legend(handles=[Line2D([], [], marker="o", ls="", color=S.MUTED, ms=4.2, label="Expert labels"),
                    Line2D([], [], marker="o", ls="", mfc="white", mec=S.MUTED, mew=1.0, ms=4.2, label="Rule")],
           loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, handletextpad=0.2, columnspacing=1.0, fontsize=6.2, borderaxespad=0.1)
fig.text(4.15 / W, (top + 0.62) / H, "B", fontsize=9, fontweight="bold", va="bottom")
fig.text(4.35 / W, (top + 0.62) / H, "Positives by expert and rule", fontsize=7, va="bottom")
S.save(fig, C.OUT, "FigS_lexical_rules")
print("saved", R, NK)
