#!/usr/bin/env python3
"""FigS_measures: the write-back policies move pooled AP, per-disease AP, MRR and Hits@10 in the same direction.
Change relative to no write-back (difference of the means in the tables).
A Hetionet random edge, graph head (Table S4a)      B Hetionet compound disjoint, graph head (Table S4b)
C PrimeKG random edge, label-only MF (Table S16a)   D PrimeKG compound disjoint, graph head (Table S16a)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
import supp_common as C
S = C.S
S.setup()

MEASURES = [("pooled", "Pooled AP", "#6a3d9a", "o"), ("perdis", "Per-disease AP", "#000000", "s"),
            ("mrr", "MRR", "#d6458f", "^"), ("hits10", "Hits@10", "#8a6d3b", "D")]
POLS = ["flat_negative", "typed", "typed_scoped", "typed_role", "typed_role_scoped", "mask_all"]
pkg = C.load_s16a()
PANELS = [("A", "Hetionet, random edge, graph head", C.load_s4("Table S4a."), "graph"),
          ("B", "Hetionet, compound disjoint, graph head", C.load_s4("Table S4b."), "graph"),
          ("C", "PrimeKG, random edge, label-only MF", pkg["random"], "mf"),
          ("D", "PrimeKG, compound disjoint, graph head", pkg["compound"], "graph")]

# Hetionet role policies are not in Table S4: their pooled and per-disease AP changes (graph head) come from Table S18b;
# MRR and Hits@10 are not available for them and are therefore not drawn.
S18 = C.contrasts_hetionet()
ROLE = ("typed_role", "typed_role_scoped")
rows = {0: POLS, 1: POLS}
fig = plt.figure(figsize=(6.7, 5.2))
gs = fig.add_gridspec(2, 2, height_ratios=[len(rows[0]) + 0.9, len(rows[1]) + 0.9], hspace=0.42, wspace=0.62)
offs = np.linspace(-0.27, 0.27, 4)
table = []
for k, (letter, title, tab, sc) in enumerate(PANELS):
    r, c = divmod(k, 2)
    ax = fig.add_subplot(gs[r, c])
    pols = rows[r]
    ys = np.arange(len(pols))[::-1]
    base = tab[("no_writeback", sc)]
    o = "random" if letter == "A" else "compound"
    for m, (mk, lab, col, marker) in enumerate(MEASURES):
        pp, xs = [], []
        for p in pols:
            if (p, sc) in tab:
                pp.append(p); xs.append(tab[(p, sc)][mk][0] - base[mk][0])
            elif letter in "AB" and p in ROLE and mk in ("pooled", "perdis"):
                pp.append(p); xs.append(S18[(o, p, "graph")][mk][0])
        yy = np.array([ys[pols.index(p)] for p in pp])
        ax.plot(xs, yy + offs[::-1][m], marker, color=col, ms=3.4, ls="none", mec=col)
        table += [(letter, p, mk, x) for p, x in zip(pp, xs)]
    if letter in "AB":                                  # note in the empty (lower) slots of the role rows
        for p in ROLE:
            ax.text(0.02, ys[pols.index(p)] - 0.18, "MRR, Hits@10 not available", transform=ax.get_yaxis_transform(), ha="left", va="center",
                    fontsize=5.4, color=S.MUTED, style="italic")
    for y in ys[:-1]:
        ax.axhline(y - 0.5, color=C.S.LIGHT, lw=0.4, zorder=0)
    ax.axvline(0, color=S.MUTED, lw=0.6, ls=(0, (3, 2)), zorder=0)
    ax.set_yticks(ys); ax.set_yticklabels([C.PL[p] for p in pols]); ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.6, len(pols) - 0.4)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 5, 10]))
    ax.set_xlabel("Change vs no write-back")
    S.panel(ax, letter, x=-0.02, title=title)
handles = [Line2D([], [], marker=m, color=col, mfc=col, ls="none", ms=3.6, label=lab) for _, lab, col, m in MEASURES]
fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.01))
fig.subplots_adjust(left=0.2, right=0.985, top=0.95, bottom=0.12)
S.save(fig, C.OUT, "FigS_measures")
print("saved FigS_measures")
# sign agreement of the four measures within each panel (policy cell with non-zero change in all four)
for letter in "ABCD":
    pols = sorted({p for l, p, _, _ in table if l == letter}, key=POLS.index)
    for p in pols:
        v = [x for l, pp, _, x in table if l == letter and pp == p]
        signs = {int(np.sign(round(x, 3))) for x in v}
        print(letter, p, [round(x, 3) for x in v], "same sign" if len(signs - {0}) <= 1 else "MIXED")
