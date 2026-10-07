#!/usr/bin/env python
"""FigS_candidate_export: exported top-100/500/1,000 candidates of the graph head (section S11, supp_S7_S14.md).
A  composition by evidence status (Table S9: approval evidence for the pair / trialled, not approved / no clinical record), 100% stacked bars
B  pair implications assigned to the candidates with a stopped trial (Table S9: defer / qualify / negate; Table S9 has no 'retain' column)
C  the 20 highest-ranked candidates (Table S10): approval evidence, registered trials, stopped trials and implication per rank"""
import re
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import b_common as C

S = C.S
S.setup()
PATH = f"{C.MS}/supp_S7_S14.md"

# ---- Table S9
hdr, rows = C.read_table(PATH, r"Table S9\.")
print(hdr)
D = []
for r in rows:
    n = int(C.num(r[0]))
    appr, trial, none, stopped = (int(C.num(x)) for x in r[1:5])
    df, ql, ng = (int(float(x)) for x in re.findall(r"\d+", r[5]))
    assert appr + trial + none == n, r                                  # the three evidence classes add up to N
    assert df + ql + ng == stopped, r                                   # every candidate with a stopped trial carries one implication
    D.append(dict(n=n, appr=appr, trial=trial, none=none, stopped=stopped, defer=df, qualify=ql, negate=ng, errors=int(C.num(r[6]))))
assert [d["n"] for d in D] == [100, 500, 1000]

# ---- Table S10
_, rows10 = C.read_table(PATH, r"Table S10\.")
T10 = [dict(rank=int(C.num(r[0])), appr=r[3] == "yes", trials=int(C.num(r[4])), stopped=int(C.num(r[5])), impl=r[6]) for r in rows10]
assert len(T10) == 20 and [t["rank"] for t in T10] == list(range(1, 21))

EV = [("appr", "Approval evidence for the same pair", "#1baf7a", S.INK), ("trial", "Registered trial(s), no approval evidence", "#eda100", S.INK),
      ("none", "No clinical record", S.LIGHT, S.INK)]
IM = [("defer", "Defer", "#4d4d4d", "white"), ("qualify", "Qualify", "#eda100", S.INK), ("negate", "Negate", "#e34948", "white")]

fig = plt.figure(figsize=(6.55, 5.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.95], hspace=0.80, wspace=0.34, left=0.125, right=0.975, top=0.93, bottom=0.07)
axA, axB, axC = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])
ys = [2, 1, 0]
for y, d in zip(ys, D):
    left = 0
    for key, lab, col, tc in EV:
        w = 100 * d[key] / d["n"]
        axA.barh(y, w, left=left, height=0.62, color=col, edgecolor="white", linewidth=0.6)
        axA.text(left + w / 2, y, str(d[key]), ha="center", va="center", fontsize=6.3, color=tc)
        left += w
axA.set_yticks(ys); axA.set_yticklabels([f"Top {d['n']:,}" for d in D]); axA.tick_params(axis="y", length=0)
axA.set_xlim(0, 100); axA.set_ylim(-0.5, 2.5); axA.set_xlabel("Candidates (%)")
axA.legend(handles=[Patch(fc=c, ec="none", label=l) for _, l, c, _ in EV], loc="upper left", bbox_to_anchor=(-0.34, -0.27), ncol=1,
           handlelength=1.0, fontsize=6.2, borderaxespad=0)
S.panel(axA, "A", x=-0.22, title="Evidence status of the exported candidates")

for y, d in zip(ys, D):
    left = 0
    for key, lab, col, tc in IM:
        w = 100 * d[key] / d["stopped"]
        axB.barh(y, w, left=left, height=0.62, color=col, edgecolor="white", linewidth=0.6)
        axB.text(left + w / 2, y, str(d[key]), ha="center", va="center", fontsize=6.3, color=tc)
        left += w
axB.set_yticks(ys); axB.set_yticklabels([f"Top {d['n']:,}\n(n = {d['stopped']})" for d in D]); axB.tick_params(axis="y", length=0)
axB.set_xlim(0, 100); axB.set_xticks([0, 20, 40, 60, 80, 100]); axB.set_ylim(-0.5, 2.5)
axB.set_xlabel("Candidates with a stopped trial (%)")
axB.legend(handles=[Patch(fc=c, ec="none", label=l) for _, l, c, _ in IM], loc="upper left", bbox_to_anchor=(0.0, -0.27), ncol=3,
           handlelength=1.0, fontsize=6.2, borderaxespad=0, columnspacing=1.0)
S.panel(axB, "B", x=-0.30, title="Pair implication of the stopped trials")

# C: top 20
rowsC = [("Approval evidence", lambda t: t["appr"], "#1baf7a"), ("Registered trial", lambda t: t["trials"] > 0, "#eda100"),
         ("Stopped trial", lambda t: t["stopped"] > 0, "#4d4d4d")]
for k, (lab, f, col) in enumerate(rowsC):
    y = 3 - k
    for t in T10:
        if f(t):
            axC.plot(t["rank"], y, "o", ms=4.6, color=col)
        else:
            axC.plot(t["rank"], y, "o", ms=4.6, mfc="none", mec=S.LIGHT, mew=0.7)
imcol = {i[1].lower(): (i[2], i[3]) for i in IM}
for t in T10:
    if t["impl"] in imcol:
        c, tc = imcol[t["impl"]]
        axC.plot(t["rank"], 0, "s", ms=5.4, color=c)
    else:
        axC.plot(t["rank"], 0, "s", ms=5.4, mfc="none", mec=S.LIGHT, mew=0.7)
axC.set_yticks([3, 2, 1, 0]); axC.set_yticklabels([r[0] for r in rowsC] + ["Pair implication"]); axC.tick_params(axis="y", length=0)
axC.set_xticks(range(1, 21)); axC.set_xlim(0.4, 20.6); axC.set_ylim(-0.6, 3.6)
axC.set_xlabel("Rank of the candidate")
axC.spines["left"].set_visible(False)
S.panel(axC, "C", x=-0.12, title="The 20 highest-ranked candidates")
axC.legend(handles=[Patch(fc=c, ec="none", label=f"Implication: {l.lower()}") for _, l, c, _ in IM], loc="lower right", bbox_to_anchor=(1.0, 1.02), ncol=3,
           handlelength=1.0, fontsize=6.0, borderaxespad=0, columnspacing=1.0, frameon=False)
S.save(fig, C.OUT, "FigS_candidate_export")
print("saved", D)
