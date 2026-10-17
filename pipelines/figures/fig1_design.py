#!/usr/bin/env python
"""Figure 1: study design (four questions on one data model; vector schematic drawn with matplotlib).
Usage: fig1_design.py OUT_DIR [NAME]"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

OUT = sys.argv[1] if len(sys.argv) > 1 else "work/figures"
NAME = sys.argv[2] if len(sys.argv) > 2 else "Figure1_design"
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as FS
FS.setup()
INK, MUTED, GRID = "#000000", "#333333", "#8c8c8c"
COLS = ["#2a78d6", "#1baf7a", "#eb6834", "#eda100"]

fig = plt.figure(figsize=(7.2, 4.0))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100); ax.set_ylim(0, 60); ax.axis("off")


def box(x, y, w, h, num, title, question, body, color):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2", lw=0.7, ec=GRID, fc="#ffffff"))
    ax.add_patch(Rectangle((x - 0.3, y + h - 6.2), 1.1, 6.2, color=color, lw=0))
    ax.text(x + 1.8, y + h - 1.0, f"{num}  {title}", ha="left", va="top", fontsize=7.6, fontweight="bold", color=INK)
    ax.text(x + 1.8, y + h - 4.3, question, ha="left", va="top", fontsize=6.4, style="italic", color=INK)
    ax.text(x + 1.8, y + h - 8.4, body, ha="left", va="top", fontsize=6.0, color=MUTED, linespacing=1.38)


def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7, lw=0.8, color=MUTED, shrinkA=0, shrinkB=0))


box(0.8, 39.0, 47.5, 20.0, "1", "Use", "Are stopped-trial pairs used as negative labels?",
    "Europe PMC full texts, searched 17 October 2026\n"
    "• articles that mention repoDB (118 records)\n"
    "• repurposing articles that pair terminated, withdrawn\n   or failed trials with negative samples (134 records)\n"
    "• 197 open-access articles screened with a codebook:\n   negatives for training, evaluation, both, resource or none", COLS[0])
box(51.7, 39.0, 47.5, 20.0, "2", "Validity", "How often are stopped-trial pairs indications?",
    "repoDB: 4,123 failed records from terminated, withdrawn\n   or suspended trials; approvals from DrugCentral and\n   Open Targets; UMLS concepts mapped to disease concepts\n"
    "Hetionet v1.0: 755 recorded treatments; PrimeKG v2:\n   9,388 indications; Open Targets 26.09 clinical reports\n"
    "Stopped trials typed by stop reason (classifier), scope\n   (same or narrower concept) and drug role (trial arms)", COLS[1])
box(0.8, 9.0, 47.5, 22.5, "3", "Mechanism", "Why do stopped pairs overlap with indications?",
    "Operational stops → compare efficacy, safety and other stops\n"
    "Granularity, attribution → restrict scope and drug role\n"
    "Selection by testing →\n"
    "• degree-preserving permutation null (popularity)\n"
    "• overlap by number of registered trials and by phase\n"
    "• logistic models adjusted for testing intensity\n"
    "• review of the 80 cleanest failure records (language-model\n   coding, verified against registry, results and publications)", COLS[2])
box(51.7, 9.0, 47.5, 22.5, "4", "Consequences", "What do stopped-trial negatives do to models?",
    "Label policies: none, flat negative, mask all, typed (+ scope,\n   + role), negate only if ≤2 trials, placebo negatives\n"
    "Scorers: degree reference, label-only MF, graph head on\n   DistMult embeddings, hybrid; registry-history scorers\n"
    "Outcomes: held-out treatments (E1; random-edge and\n   compound-disjoint), external approved indications (E2),\n   later failures (E3); per-disease AP, MRR, AUROC;\n   paired disease-cluster bootstrap", COLS[3])
arrow(48.6, 49.0, 51.4, 49.0); arrow(56.0, 38.7, 44.0, 32.3); arrow(48.6, 20.5, 51.4, 20.5)
ax.add_patch(FancyBboxPatch((0.8, 0.8), 98.4, 5.6, boxstyle="round,pad=0.3,rounding_size=1.2", lw=0.7, ec=GRID, fc="#f4f4f4"))
ax.text(2.6, 3.6, "Outputs:  a label policy (mask, do not negate) • a registry-history baseline for evaluation • model-free checks for any "
        "negative-label source • an evidence-record schema", ha="left", va="center", fontsize=6.3, color=INK)
arrow(24.5, 8.7, 24.5, 6.9); arrow(75.4, 8.7, 75.4, 6.9)
FS.save(fig, OUT, NAME)
print("saved", NAME)
