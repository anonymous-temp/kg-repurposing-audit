#!/usr/bin/env python3
"""FigS_full_results: per-disease AP of every scorer under each write-back policy (annotated heatmaps).
A-C Hetionet (Tables S4a, S4b, S4c), D-F PrimeKG (Table S16a): random edge, compound disjoint, external approvals (E2).
Values parsed from the manuscript supplement; policies missing from a table are left grey."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import supp_common as C
S = C.S
S.setup()

POLS = ["no_writeback", "flat_negative", "typed", "typed_scoped", "typed_role", "typed_role_scoped", "mask_all"]
SCORERS_HET = ["degree", "mf", "graph", "hybrid"]
SCORERS_PKG = ["degree", "mf", "graph"]                       # no hybrid scorer was run on PrimeKG
COLLAB = {"degree": "Degree\nreference", "mf": "Label-only\nMF", "graph": "Graph\nhead", "hybrid": "Hybrid"}

het = {"random": C.load_s4("Table S4a."), "compound": C.load_s4("Table S4b."), "e2": C.load_s4("Table S4c.")}
pkg = C.load_s16a()

# Hetionet role-restricted policies are not in Tables S4a-c: cell = no-write-back per-disease AP of the same scorer and outcome set (S4)
# + per-disease AP difference "typed_role(_scoped) - no write-back" of Table S18b (as main Table 4 was built). Grey if S18b lacks the cell.
H18 = C.contrasts_hetionet()
derived = {}
for o in C.OUTCOME_KEYS:
    for p in ("typed_role", "typed_role_scoped"):
        for sc in SCORERS_HET:
            if (o, p, sc) in H18 and ("no_writeback", sc) in het[o] and (p, sc) not in het[o]:
                v = round(het[o][("no_writeback", sc)]["perdis"][0] + H18[(o, p, sc)]["perdis"][0], 3)
                het[o][(p, sc)] = {"perdis": (v, float("nan"))}
                derived[(o, p, sc)] = v
TITLES = {"random": "random edge", "compound": "compound disjoint", "e2": "external approvals"}

# light-to-dark single-hue blue; truncated so that the lightest cell is still distinguishable from the grey of missing cells
base = plt.get_cmap("Blues")
cmap = LinearSegmentedColormap.from_list("blues_t", base(np.linspace(0.06, 0.92, 256)))
cmap.set_bad("#c9c9c9")

fig, axes = plt.subplots(2, 3, figsize=(6.55, 4.7))
missing = []
letters = iter("ABCDEF")
for r, (gname, data, scorers) in enumerate([("Hetionet", het, SCORERS_HET), ("PrimeKG", pkg, SCORERS_PKG)]):
    for c, o in enumerate(C.OUTCOME_KEYS):
        ax = axes[r][c]
        tab = data[o]
        M = np.full((len(POLS), len(scorers)), np.nan)
        for i, p in enumerate(POLS):
            for j, sc in enumerate(scorers):
                d = tab.get((p, sc))
                if d is None:
                    missing.append((gname, o, p, sc))
                else:
                    M[i, j] = d["perdis"][0]
        vmax = np.nanmax(M)
        ax.pcolormesh(np.ma.masked_invalid(M), cmap=cmap, vmin=0, vmax=vmax, edgecolors="white", linewidth=1.0)
        ax.set_xlim(0, len(scorers)); ax.set_ylim(len(POLS), 0)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if not np.isnan(M[i, j]):
                    ax.text(j + 0.5, i + 0.5, f"{M[i, j]:.3f}", ha="center", va="center", fontsize=5.8,
                            color="white" if M[i, j] / vmax > 0.55 else "black")
        ax.set_xticks(np.arange(len(scorers)) + 0.5); ax.set_xticklabels([COLLAB[s] for s in scorers], fontsize=5.8)
        ax.set_yticks(np.arange(len(POLS)) + 0.5)
        ax.set_yticklabels([C.PL[p] for p in POLS] if c == 0 else [])
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        S.panel(ax, next(letters), title=f"{gname}, {TITLES[o]}")

fig.subplots_adjust(left=0.15, right=0.995, top=0.94, bottom=0.075, wspace=0.07, hspace=0.30)
S.save(fig, C.OUT, "FigS_full_results")
print("saved FigS_full_results")
print("missing cells (policy not in source table):")
for g, o, p, s in missing:
    print("  ", g, o, p, s)

print("derived Hetionet role cells (no write-back + S18b difference):", len(derived))
# check: graph-head values against main Table 4 (archive v5)
T4 = os.path.join(C.MS, "bmc_tables_4_5.md")
_, t4 = C.table_dicts(T4, "Table 4*")
col = {"random": "RE per-disease AP (\u0394; 95% CI)", "compound": "CD per-disease AP (\u0394; 95% CI)", "e2": "E2 per-disease AP (\u0394; 95% CI)"}
print("graph head: computed (S4 + S18b) vs Table 4")
nbad = 0
for r in t4:
    try:
        p = C.policy_key(r["Policy"])
    except KeyError:
        continue
    if p not in ("no_writeback", "flat_negative", "typed", "typed_scoped", "typed_role", "typed_role_scoped", "mask_all"):
        continue
    out = []
    for o in C.OUTCOME_KEYS:
        t4v = C.num(r[col[o]].split("(")[0])
        mine = het[o][(p, "graph")]["perdis"][0]
        ok = abs(t4v - mine) < 5e-4
        nbad += not ok
        out.append(f"{o}: {mine:.3f} vs {t4v:.3f} {'ok' if ok else 'DIFF'}")
    print(f"  {C.PL[p]:<22}", " | ".join(out))
print("mismatches:", nbad)
