#!/usr/bin/env python
"""Draw every figure. Usage: 08_make_figures.py [RESULTS ROLE_RESULTS PRIMEKG_RESULTS TESTEDNESS_RESULTS FIGURES_DIR]"""
import subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent / "figures"
res, role, pkg, tested, out = sys.argv[1:6] if len(sys.argv) > 5 else ("work/results", "work/results_role", "work/results_primekg",
                                                                        "work/results_testedness", "work/figures")
Path(out).mkdir(parents=True, exist_ok=True)
run = lambda *a: subprocess.run([sys.executable, *map(str, a)], check=True)
run(HERE / "fig1_workflow.py", out)
run(HERE / "fig2_stopped.py", res, pkg, role, out)
run(HERE / "fig_enrichment.py", role, pkg, out)          # Figure 3
run(HERE / "fig_testedness.py", tested, out)             # Figure 4
run(HERE / "fig_policies.py", role, pkg, out)            # Figure 5
run(HERE / "fig_negation.py", role, pkg, out)            # Figure 6
run(HERE / "fig_results.py", res, out)                   # Figures S1-S3
