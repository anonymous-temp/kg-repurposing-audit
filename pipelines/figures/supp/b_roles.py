#!/usr/bin/env python
"""FigS_roles: A role of the mapped compound in stopped trials (Table S13, trial-compound combinations, both graphs);
B recorded or approved indications among pairs with a scientific stop before 2015, by role subset (Table S14, Wilson 95% intervals)."""
import re, textwrap
import matplotlib.pyplot as plt
import b_common as C

S = C.S
S.setup()
PATH = f"{C.MS}/supp_S16_S18.md"
G = ("Hetionet", "PrimeKG")

# ---- Table S13
_, rows = C.read_table(PATH, r"Table S13\.")
tot = [r for r in rows if r[0] == "Total"][0]
roles = [r for r in rows if r[0] != "Total"]
N = {"Hetionet": int(C.num(tot[2])), "PrimeKG": int(C.num(tot[3]))}
for j, g in enumerate(G):                                               # the rows must sum to the stated total
    assert sum(int(C.num(r[2 + j])) for r in roles) == N[g], g
A = [dict(role=r[0], Hetionet=int(C.num(r[2])), PrimeKG=int(C.num(r[3]))) for r in roles]

# ---- Table S14
_, rows = C.read_table(PATH, r"Table S14\.")
B = []
for r in rows:
    d = dict(label=r[0])
    for j, g in enumerate(G):
        n = int(C.num(r[1 + 2 * j])); cell = r[2 + 2 * j]
        k = int(C.num(cell)); pct = nums = re.search(r"\((\d+)%\)", cell)
        assert abs(100 * k / n - float(pct.group(1))) <= 0.51
        d[g] = (k, n) + C.wilson(k, n)
    B.append(d)
assert len(B) == 5


def pretty14(lab):
    m = re.match(r"(Same-concept scientific stop|Scientific stop), (.*)", lab)
    head = "Same-concept scientific stop" if m.group(1).startswith("Same") else "Scientific stop"
    tail = {"at least one as investigational agent": "at least one as investigational agent",
            "only as comparator or background therapy": "only as comparator or\nbackground therapy",
            "only other roles": "only other roles"}[m.group(2)]
    return head, tail


fig = plt.figure(figsize=(6.7, 6.9))
gs = fig.add_gridspec(2, 1, height_ratios=[len(A) * 1.0, len(B) * 1.35 + 0.8], hspace=0.30, left=0.30, right=0.80, top=0.965, bottom=0.07)
axA, axB = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

# panel A: grouped bars
ys = list(range(len(A)))[::-1]
H = 0.36
for j, g in enumerate(G):
    for y, d in zip(ys, A):
        pct = 100 * d[g] / N[g]
        yy = y + (H / 2 if j == 0 else -H / 2)
        axA.barh(yy, pct, height=H - 0.04, color=S.COL[g], edgecolor="none", label=f"{g} (n = {N[g]:,})" if y == ys[0] else None)
        axA.text(pct + 0.7, yy, f"{pct:.1f}%  ({d[g]:,})", va="center", ha="left", fontsize=6)
axA.set_yticks(ys); axA.set_yticklabels([d["role"] for d in A]); axA.tick_params(axis="y", length=0)
axA.set_xlim(0, 50); axA.set_ylim(-0.6, len(A) - 0.4)
axA.set_xlabel("Trial–compound combinations (%)")
axA.legend(loc="lower right", handlelength=1.2, bbox_to_anchor=(1.0, 0.02))
S.panel(axA, "A", x=-0.52, title="Role of the mapped compound in stopped trials")

# panel B: dot plot with Wilson intervals
ysb = [4.0, 3.0, 2.0, 0.6, -0.4]            # gap between the two groups of subsets
for g, off in (("Hetionet", 0.17), ("PrimeKG", -0.17)):
    for y, d in zip(ysb, B):
        k, n, lo, hi = d[g]
        S.errorbar(axB, 100 * k / n, y + off, 100 * lo, 100 * hi, S.COL[g], S.MK[g], ms=3.6)
        axB.text(1.03, y + off, f"{100 * k / n:.0f}%  ({k}/{n})", transform=axB.get_yaxis_transform(), va="center", ha="left",
                 fontsize=6, color=S.COL[g] if False else S.INK)
axB.set_yticks(ysb)
axB.set_yticklabels([pretty14(d["label"])[0] + ",\n" + pretty14(d["label"])[1] for d in B])
axB.tick_params(axis="y", length=0)
axB.set_xlim(0, 100); axB.set_ylim(-0.95, 4.6)
axB.set_xlabel("Recorded treatments or approved indications (%, Wilson 95% interval)")
axB.text(1.03, 4.55, "Share (recorded or\napproved / pairs)", transform=axB.get_yaxis_transform(), va="bottom", fontsize=6, color=S.MUTED)
from matplotlib.lines import Line2D
axB.legend(handles=[Line2D([], [], marker=S.MK[g], ls="", color=S.COL[g], ms=3.6, label=g) for g in G], loc="lower right", handletextpad=0.2,
           ncol=2, columnspacing=1.0, bbox_to_anchor=(1.0, 1.0))
S.panel(axB, "B", x=-0.52, title="Indications among pairs with a scientific stop before 2015, by role")
S.save(fig, C.OUT, "FigS_roles")
print("saved", N, [(d["role"], round(100 * d["Hetionet"] / N["Hetionet"], 1), round(100 * d["PrimeKG"] / N["PrimeKG"], 1)) for d in A])
