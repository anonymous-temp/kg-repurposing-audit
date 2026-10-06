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

fig, ax = plt.subplots(figsize=(7.2, 4.9))
ax.set_xlim(0, 100); ax.set_ylim(-2.5, 63); ax.axis("off")


def box(x, y, w, h, title, body, edge=GRID, title_color=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2", lw=0.7, ec=GRID, fc="#ffffff"))
    ax.text(x + 0.9, y + h - 1.2, title, ha="left", va="top", fontsize=7.2, fontweight="bold", color=title_color)
    ax.text(x + 0.9, y + h - 4.4, body, ha="left", va="top", fontsize=6.0, color=MUTED, linespacing=1.35)


def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7, lw=0.8, color=MUTED,
                                 shrinkA=0, shrinkB=0))


ax.text(0.5, 62.5, "A   Evidence assembly", fontsize=8, fontweight="bold", va="top")
box(0.5, 41, 22.5, 18.5, "Knowledge graphs", "Hetionet v1.0 (2016): 755\nrecorded treatments (CtD),\n390 palliative pairs\nPrimeKG v2 (2022): 9,388\nindications, 2,568 off-label")
box(0.5, 18.5, 22.5, 19.5, "Trial and approval data", "Open Targets 26.09 clinical\nreports (ClinicalTrials.gov,\nChEMBL, DailyMed, FDA, EMA,\nPMDA, TTD); stop-reason\nclassifier. ClinicalTrials.gov\nAPI: arms of 29,485 trials")
box(26.0, 18.5, 24.0, 41, "Mapping and typing",
    "Drugs: DrugBank \u2192 ChEMBL\nDiseases: DOID or MONDO \u2192\nEFO/MONDO (or nearest\nmapped ancestor)\n\nStopped trial (terminated,\nwithdrawn or suspended)\n\u2192 failure type: efficacy,\n    safety, operational, design,\n    uninformative, not reported\n\u2192 scope: same or narrower\n    concept than graph disease\n\u2192 drug role: investigational,\n    comparator, background\n    therapy or other", edge=BLUE)
arrow(23.6, 50, 25.5, 50); arrow(23.6, 28, 25.5, 28)

ax.text(53, 62.5, "B   Write-back and evaluation", fontsize=8, fontweight="bold", va="top")
box(53, 31.5, 21.5, 28, "Write-back policy",
    "Pairs with trials stopped\nand started before 2015:\n\u2022 no write-back\n\u2022 flat negative\n\u2022 mask all\n\u2022 typed: efficacy/safety\n   negative, others masked\n\u2022 typed + scope\n\u2022 typed + role: negative\n   only for the drug under test\n\u2022 typed + role + scope", edge=ORANGE)
box(77.5, 31.5, 22, 28, "Scorers",
    "Degree reference\nDisease degree\nLabel-only matrix\n   factorisation (MF)\nGraph head on DistMult\n   embeddings trained with\n   no drug\u2013disease edges\nHybrid (MF + graph;\n   Hetionet only)")
arrow(50.6, 46, 52.5, 46); arrow(75.1, 46, 77.0, 46)
box(53, 11.8, 46.5, 17.7, "Evaluation on full candidate grids",
    "E1  held-out recorded treatments (random-edge and\n       compound-disjoint tasks; 5 or 3 partitions)\n"
    "E2  approved indications absent from the graph\n"
    "E3  approved indications vs later (2017+) scientific failures\n"
    "Tested vs approved: approved, failed, newly tested\n       and untested pairs compared\n"
    "AP, per-disease AP, MRR, AUROC; disease-cluster bootstrap", edge=AQUA)
arrow(88.5, 31.0, 88.5, 30.0)

ax.text(0.5, 12.2, "C   Evidence handoff", fontsize=8, fontweight="bold", va="top")
box(0.5, -1.5, 99, 11.5, "Candidate record (LinkML schema with Biolink 4.4.5 mappings)",
    "prediction provenance (model, task, policy, rank scope) + dated clinical evidence (approvals, trials)\n+ failure annotations (failure type, scope match) "
    "\u2192 implication: retain (recorded treatment) | negate (same-concept efficacy or\nsafety failure) | qualify (failure in another concept) | defer (operational, design, uninformative or unreported reason)")
arrow(38.0, 18.0, 38.0, 10.6)
FS.save(fig, OUT, "Figure1_workflow")
print("saved", OUT)
