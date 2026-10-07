#!/usr/bin/env python
"""Draw every figure (main Figures 1-9 and supplementary Figures S1-S16).
Usage: 08_make_figures.py [RESULTS ROLE_RESULTS PRIMEKG_RESULTS TESTEDNESS_RESULTS SELECTION_RESULTS EXTRA_HETIONET EXTRA_PRIMEKG FIGURES_DIR]"""
import subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent / "figures"
res, role, pkg, tested, sel, xhet, xpkg, out = sys.argv[1:9] if len(sys.argv) > 8 else (
    "work/results", "work/results_role", "work/results_primekg", "work/results_testedness", "work/results_selection",
    "work/results_extra", "work/results_extra_primekg", "work/figures")
Path(out).mkdir(parents=True, exist_ok=True)
run = lambda *a: subprocess.run([sys.executable, *map(str, a)], check=True)
run(HERE / "fig1_workflow.py", out)                               # Figure 1
run(HERE / "fig2_stopped.py", res, pkg, role, out)                # Figure 2
run(HERE / "fig_enrichment.py", role, pkg, out, sel)              # Figure 3
run(HERE / "fig_intensity.py", sel, out)                          # Figure 4
run(HERE / "fig_casereview.py", f"{sel}/case_taxonomy_verified.tsv", out)   # Figure 5
run(HERE / "fig_testedness.py", tested, out, sel)                 # Figure 6 (file Figure5_tested_vs_approved)
run(HERE / "fig_policies.py", role, pkg, out, xhet, xpkg)         # Figure 7 (file Figure6_policies)
run(HERE / "fig_negation.py", role, pkg, out, xhet, xpkg, sel)    # Figure 8 (file Figure7_negatives)
run(HERE / "fig_lossmap.py", sel, f"{role}/e1_strata.tsv", out)   # Figure 9 (needs 22_rank_strata.py output in sel)
run(HERE / "fig_results.py", res, out)                            # Figures S3, S10, S11
for f in sorted((HERE / "supp").glob("[bf]*_*.py")):              # Figures S1, S2, S4-S9, S12-S16 (read the markdown tables of 09*)
    if not f.name.endswith("common.py"):
        run(f)
