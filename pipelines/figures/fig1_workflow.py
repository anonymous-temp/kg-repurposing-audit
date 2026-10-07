#!/usr/bin/env python
"""Figure 1: study workflow (vector schematic drawn with matplotlib; no generative image model)."""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = sys.argv[1] if len(sys.argv) > 1 else "work/figures"
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as FS
FS.setup()
INK, MUTED, GRID = "#000000", "#333333", "#8c8c8c"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"

fig, ax = plt.subplots(figsize=(7.2, 5.4))
ax.set_xlim(0, 100); ax.set_ylim(-4.2, 70); ax.axis("off")


def box(x, y, w, h, title, body, edge=GRID, title_color=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2", lw=0.7, ec=GRID, fc="#ffffff"))
    ax.text(x + 0.9, y + h - 1.2, title, ha="left", va="top", fontsize=7.2, fontweight="bold", color=title_color)
    ax.text(x + 0.9, y + h - 4.4, body, ha="left", va="top", fontsize=6.0, color=MUTED, linespacing=1.35)


def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7, lw=0.8, color=MUTED,
                                 shrinkA=0, shrinkB=0))


ax.text(0.5, 69.5, "A   Evidence assembly and model-free tests", fontsize=8, fontweight="bold", va="top")
box(0.5, 47, 22.5, 19.5, "Knowledge graphs", "Hetionet v1.0 (2016): 755\nrecorded treatments (CtD),\n390 palliative pairs\nPrimeKG v2 (2022): 9,388\nindications, 2,568 off-label")
box(0.5, 24.5, 22.5, 19.5, "Trial and approval data", "Open Targets 26.09 clinical\nreports (ClinicalTrials.gov,\nChEMBL, DailyMed, FDA, EMA,\nPMDA, TTD); stop-reason\nclassifier. ClinicalTrials.gov\nAPI: arms of 29,485 trials")
box(26.0, 11.5, 24.0, 55, "Mapping, typing and tests",
    "Drugs: DrugBank \u2192 ChEMBL\nDiseases: DOID or MONDO \u2192\nEFO/MONDO (or nearest\nmapped ancestor)\n\nStopped trial (terminated,\nwithdrawn or suspended)\n\u2192 failure type: efficacy,\n    safety, operational, design,\n    uninformative, not reported\n\u2192 scope: same or narrower\n    concept than graph disease\n\u2192 drug role: investigational,\n    comparator, background\n    therapy or other\n\nTests of selection:\n• indications by number of\n   trials, phase, disease area\n• logistic models adjusted\n   for testing intensity\n• degree-preserving\n   permutation null\n• review of 80 indications\n   with failed trials", edge=BLUE)
arrow(23.6, 56, 25.5, 56); arrow(23.6, 34, 25.5, 34)

ax.text(53, 69.5, "B   Write-back and evaluation", fontsize=8, fontweight="bold", va="top")
box(53, 31.5, 21.5, 35, "Write-back policy",
    "Pairs with trials stopped\nand started before 2015:\n\u2022 no write-back\n\u2022 flat negative\n\u2022 mask all\n\u2022 typed: efficacy/safety\n   negative, others masked\n\u2022 typed + scope\n\u2022 typed + role: negative\n   only for the drug under test\n\u2022 typed + role + scope\n\u2022 typed + testing intensity:\n   negative only if \u22642 trials\n\u2022 placebo negatives: same\n   number, same drug and\n   disease counts", edge=ORANGE)
box(77.5, 31.5, 22, 35, "Scorers",
    "Degree reference\nDisease degree\nLabel-only matrix\n   factorisation (MF)\nGraph head on DistMult\n   embeddings trained with\n   no drug\u2013disease edges\nHybrid (MF + graph;\n   Hetionet only)\n\nRegistry history only:\n   number of trials,\n   number of stopped trials,\n   any stopped trial")
arrow(50.6, 50, 52.5, 50); arrow(75.1, 50, 77.0, 50)
box(53, 11.5, 46.5, 18.0, "Evaluation on full candidate grids",
    "E1  held-out recorded treatments (random-edge and\n       compound-disjoint tasks; 5 or 3 partitions)\n"
    "E2  approved indications absent from the graph\n"
    "E3  approved indications vs later (2017+) scientific failures\n"
    "Tested vs approved: approved, failed, newly tested\n       and untested pairs compared\n"
    "AP, per-disease AP, MRR, AUROC; disease-cluster bootstrap", edge=AQUA)
arrow(88.5, 31.0, 88.5, 30.0)

ax.text(0.5, 9.9, "C   Evidence handoff", fontsize=8, fontweight="bold", va="top")
box(0.5, -3.4, 99, 10.2, "Candidate record (LinkML schema with Biolink 4.4.5 mappings)",
    "prediction provenance (model, task, policy, rank scope) + dated clinical evidence (approvals, trials)\n+ failure annotations (failure type, scope match, drug role) "
    "\u2192 implication: retain (recorded treatment) | negate (same-concept efficacy or\nsafety failure) | qualify (failure in another concept) | defer (operational, design, uninformative or unreported reason)")
arrow(38.0, 11.0, 38.0, 7.3)
FS.save(fig, OUT, "Figure1_workflow")
print("saved", OUT)
