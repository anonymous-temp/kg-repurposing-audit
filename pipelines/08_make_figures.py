#!/usr/bin/env python
"""Draw every figure.
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
run(HERE / "fig_testedness.py", tested, out, sel)                 # Figure 5
run(HERE / "fig_policies.py", role, pkg, out, xhet, xpkg)         # Figure 6
run(HERE / "fig_negation.py", role, pkg, out, xhet, xpkg, sel)    # Figure 7
run(HERE / "fig_results.py", res, out)                            # Figures S1-S3
