#!/usr/bin/env python
"""Figure 1: study workflow (vector schematic drawn with matplotlib; no generative image model)."""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = sys.argv[1] if len(sys.argv) > 1 else "work/figures"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "pdf.fonttype": 42, "svg.fonttype": "none"})
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"

fig, ax = plt.subplots(figsize=(7.2, 4.3))
ax.set_xlim(0, 100); ax.set_ylim(-2.5, 57); ax.axis("off")


def box(x, y, w, h, title, body, edge=GRID, title_color=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2", lw=0.8, ec=edge, fc="#fcfcfb"))
    ax.text(x + 0.9, y + h - 1.2, title, ha="left", va="top", fontsize=7.2, fontweight="bold", color=title_color)
    ax.text(x + 0.9, y + h - 4.4, body, ha="left", va="top", fontsize=6.0, color=MUTED, linespacing=1.35)


def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7, lw=0.8, color=MUTED,
                                 shrinkA=0, shrinkB=0))


ax.text(0.5, 56.5, "a  Evidence assembly", fontsize=8, fontweight="bold", va="top")
box(0.5, 38, 22.5, 15.5, "Hetionet v1.0 (2016)", "47,031 nodes, 24 relations\n755 recorded treatments (CtD)\n390 palliative pairs (CpD)")
box(0.5, 18.5, 22.5, 16.5, "Open Targets 26.09", "Clinical reports from\nClinicalTrials.gov, ChEMBL,\nDailyMed, FDA, EMA, PMDA,\nTTD; stop-reason classifier")
box(26.0, 18.5, 24.0, 35, "Mapping and typing",
    "Drugs: DrugBank \u2192 ChEMBL\nDiseases: DOID \u2192 EFO/MONDO\n(or nearest mapped ancestor)\n\nStopped trial (terminated,\nwithdrawn or suspended)\n\u2192 failure type: efficacy,\n    safety, operational, design,\n    uninformative, not reported\n\u2192 scope: same or narrower\n    concept than graph disease", edge=BLUE)
arrow(23.6, 45.5, 25.5, 45.5); arrow(23.6, 26.5, 25.5, 26.5)

ax.text(53, 56.5, "b  Write-back and evaluation", fontsize=8, fontweight="bold", va="top")
box(53, 31.5, 21.5, 22, "Write-back policy",
    "Pairs with trials stopped\nand started before 2015:\n\u2022 no write-back\n\u2022 flat negative\n\u2022 mask all\n\u2022 typed: efficacy/safety\n   negative, others masked\n\u2022 typed + scope: negative\n   only in the same concept", edge=ORANGE)
box(77.5, 31.5, 22, 22, "Scorers",
    "Degree reference\nDisease degree\nLabel-only matrix\n   factorisation (MF)\nGraph head on DistMult\n   embeddings trained with\n   no treatment edges\nHybrid (MF + graph)")
arrow(50.6, 44, 52.5, 44); arrow(75.1, 44, 77.0, 44)
box(53, 13.5, 46.5, 15, "Evaluation on full candidate grids",
    "E1  held-out recorded treatments (random-edge and\n       compound-disjoint tasks, 5 partitions each)\n"
    "E2  approved indications absent from Hetionet\n"
    "E3  approved indications vs later (2017+) scientific failures\n"
    "AP, per-disease AP, MRR, AUROC; disease-cluster bootstrap", edge=AQUA)
arrow(88.5, 31.0, 88.5, 29.0)

ax.text(0.5, 12.2, "c  Evidence handoff", fontsize=8, fontweight="bold", va="top")
box(0.5, -1.5, 99, 11.5, "Candidate record (LinkML schema with Biolink 4.4.5 mappings)",
    "prediction provenance (model, task, policy, rank scope) + dated clinical evidence (approvals, trials)\n+ failure annotations (failure type, scope match) "
    "\u2192 implication: retain (recorded treatment) | negate (same-concept efficacy or\nsafety failure) | qualify (failure in another concept) | defer (operational, design, uninformative or unreported reason)")
arrow(38.0, 18.0, 38.0, 10.6)
fig.savefig(f"{OUT}/Figure1_workflow.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/Figure1_workflow.png", dpi=600, bbox_inches="tight")
print("saved", OUT)
