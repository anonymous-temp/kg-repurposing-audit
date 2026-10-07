#!/usr/bin/env python
"""FigS_case_review: review of the 80 indications with the cleanest failure records (Table S21 / case_taxonomy_verified.tsv).
A  confusion matrix of the primary code, first coding (rows) vs blind second coding (columns); % agreement and Cohen's kappa computed here
B  same for the contradiction judgement (yes / unclear / no)
C  secondary codes by graph (number of cases carrying the code in the `secondary` column)
D  primary (final) code in oncology vs other diseases (keyword classification of the disease name), 100% stacked bars"""
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.colors import PowerNorm
import b_common as C

S = C.S
S.setup()
TSV = C.CASES
rows = list(csv.DictReader(open(TSV, encoding="utf-8"), delimiter="\t"))
assert len(rows) == 80
CODES = ["A1", "A2", "B1", "B2", "C", "D", "E", "F"]
LAB = {"A1": "A1 Special population", "A2": "A2 Stage, line or goal", "B1": "B1 Regimen of the drug itself", "B2": "B2 Combination or add-on",
       "C": "C Comparative question", "D": "D Endpoint or subtype", "E": "E Stop not drug-attributable", "F": "F Other"}
CCOL = {"A1": "#1baf7a", "A2": "#008300", "B1": "#e87ba4", "B2": "#4a3aa7", "C": "#eda100", "D": "#e34948", "E": "#4d4d4d", "F": "#d9d9d9"}
ONC = ["cancer", "carcinoma", "tumour", "tumor", "neoplasm", "leukaemia", "leukemia", "lymphoma", "myeloma", "sarcoma", "melanoma", "glioma",
       "blastoma", "malignan"]


def confusion(a, b, labs):
    M = np.zeros((len(labs), len(labs)), int)
    for x, y in zip(a, b):
        M[labs.index(x), labs.index(y)] += 1
    n = M.sum(); po = np.trace(M) / n; pe = (M.sum(0) * M.sum(1)).sum() / n ** 2
    return M, po, (po - pe) / (1 - pe)


M1, po1, k1 = confusion([r["first_coding_primary"] for r in rows], [r["blind_coding_primary"] for r in rows], CODES)
M2, po2, k2 = confusion([r["first_coding_contradicts"] for r in rows], [r["blind_coding_contradicts"] for r in rows], ["yes", "unclear", "no"])
print(f"primary: agreement {100 * po1:.1f}% ({np.trace(M1)}/80), kappa {k1:.3f}; contradiction: agreement {100 * po2:.1f}% ({np.trace(M2)}/80), kappa {k2:.3f}")

for r in rows:
    r["onc"] = any(k in r["disease"].lower() for k in ONC)
cnt = {(g, o): sum(1 for r in rows if r["graph"] == g and r["onc"] == o) for g in ("hetionet", "primekg") for o in (True, False)}
print("oncology / other:", {g: (cnt[(g, True)], cnt[(g, False)]) for g in ("hetionet", "primekg")},
      "pooled", (sum(r["onc"] for r in rows), sum(not r["onc"] for r in rows)))

sec = {g: {c: 0 for c in CODES} for g in ("hetionet", "primekg")}
for r in rows:
    for c in {s.strip() for s in r["secondary"].split(",") if s.strip()}:
        sec[r["graph"]][c] += 1

W, H = 6.7, 7.4
fig = plt.figure(figsize=(W, H))


def axin(l, b, w, h):
    return fig.add_axes([l / W, b / H, w / W, h / H])


def title(x, y, letter, text):
    """Panel letter and title at absolute position (inches from the lower left); the title may have two lines."""
    fig.text(x / W, y / H, letter, fontsize=9, fontweight="bold", va="bottom", ha="left")
    fig.text((x + 0.2) / W, y / H, text, fontsize=7, va="bottom", ha="left")


def heat(ax, M, xlabs, ylabs):
    n = len(xlabs)
    ax.imshow(M, cmap="Greys", norm=PowerNorm(0.5, vmin=0, vmax=M.max() * 1.25), aspect="equal", interpolation="nearest")
    for i in range(n):
        for j in range(n):
            v = M[i, j]
            dark = (v / (M.max() * 1.25)) ** 0.5 > 0.55
            ax.text(j, i, str(v), ha="center", va="center", fontsize=6.5, fontweight="bold" if i == j else "normal",
                    color="white" if dark else (S.INK if v else "#a6a6a6"))
        ax.add_patch(Rectangle((i - 0.5, i - 0.5), 1, 1, fill=False, ec=S.INK, lw=0.9, zorder=3))
    ax.set_xticks(range(n)); ax.set_xticklabels(xlabs)
    ax.set_yticks(range(n)); ax.set_yticklabels(ylabs)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlabel("Blind second coding"); ax.set_ylabel("First coding")


top = H - 0.85
cell = 0.31
axA = axin(1.75, top - 8 * cell, 8 * cell, 8 * cell)
heat(axA, M1, CODES, [LAB[c] for c in CODES])
title(0.1, top + 0.08, "A", f"Primary code: {100 * po1:.1f}% agreement,\nCohen\u2019s \u03ba = {k1:.2f} (n = 80)")
c3 = 0.50
axB = axin(5.05, top - 3 * c3, 3 * c3, 3 * c3)
heat(axB, M2, ["Yes", "Unclear", "No"], ["Yes", "Unclear", "No"])
title(4.0, top + 0.08, "B", f"Contradiction judgement:\n{100 * po2:.1f}% agreement, \u03ba = {k2:.2f}")
axB.text(0.5, -0.62, "Cells outlined in black: agreement", transform=axB.transAxes, ha="center", fontsize=5.8, color=S.MUTED) if False else None

# C: secondary codes by graph
row2_top = top - 8 * cell - 0.72
hgt = 2.35
axC = axin(1.75, row2_top - hgt, 1.3, hgt)
ys = list(range(len(CODES)))[::-1]
hh = 0.36
for g, name, off in (("hetionet", "Hetionet", hh / 2), ("primekg", "PrimeKG", -hh / 2)):
    for y, c in zip(ys, CODES):
        v = sec[g][c]
        axC.barh(y + off, v, height=hh - 0.04, color=S.COL[name], edgecolor="none", label=f"{name}" if c == "A1" else None)
        axC.text(v + 0.3, y + off, str(v), va="center", ha="left", fontsize=6)
axC.set_yticks(ys); axC.set_yticklabels([LAB[c] for c in CODES]); axC.tick_params(axis="y", length=0)
axC.set_xlim(0, 17.5); axC.set_xticks([0, 5, 10, 15]); axC.set_ylim(-0.6, len(CODES) - 0.4)
axC.set_xlabel("Cases (of 40 per graph)")
axC.legend(loc="lower right", handlelength=1.1, bbox_to_anchor=(1.0, 0.0), borderaxespad=0.1)
title(0.1, row2_top + 0.08, "C", "Secondary features by graph")

# D: 100% stacked bars by primary code
axD = axin(4.3, row2_top - hgt, 2.2, hgt)
bars = []
for gname, sel in (("All cases", lambda r: True), ("Hetionet", lambda r: r["graph"] == "hetionet"), ("PrimeKG", lambda r: r["graph"] == "primekg")):
    for onc in (True, False):
        sub = [r for r in rows if sel(r) and r["onc"] == onc]
        bars.append((gname, "Oncology" if onc else "Other diseases", sub))
ypos = [5.55, 4.55, 3.05, 2.05, 0.55, -0.45]
for y, (gname, lab, sub) in zip(ypos, bars):
    left = 0.0
    n = len(sub)
    for c in CODES:
        k = sum(1 for r in sub if r["primary"] == c)
        if k == 0:
            continue
        w = 100 * k / n
        axD.barh(y, w, left=left, height=0.86, color=CCOL[c], edgecolor="white", linewidth=0.6)
        axD.text(left + w / 2, y, str(k), ha="center", va="center", fontsize=6, color=S.INK if c in ("C", "F") else "white")
        left += w
axD.set_yticks(ypos); axD.set_yticklabels([f"{b[1]} (n = {len(b[2])})" for b in bars]); axD.tick_params(axis="y", length=0)
axD.set_xlim(0, 100); axD.set_ylim(-1.0, 6.5)
axD.set_xlabel("Cases (%)")
for y, gname in ((6.15, "All cases"), (3.65, "Hetionet"), (1.15, "PrimeKG")):
    axD.text(-0.02, y, gname, transform=axD.get_yaxis_transform(), ha="right", va="center", fontsize=6.5, fontweight="bold")
title(3.05, row2_top + 0.08, "D", "Primary code, oncology vs other diseases")
fig.legend(handles=[Patch(fc=CCOL[c], ec="none", label=LAB[c] + (" (never final primary)" if c == "B1" else "")) for c in CODES], loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.0),
           handlelength=1.0, columnspacing=1.2, fontsize=6.2, title="Primary code (panel D)", title_fontsize=6.2)
S.save(fig, C.OUT, "FigS_case_review")
