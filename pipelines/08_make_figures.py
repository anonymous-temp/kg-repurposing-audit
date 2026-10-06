#!/usr/bin/env python
"""Draw every figure. Usage: 08_make_figures.py RESULTS_DIR FIGURES_DIR"""
import subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent / "figures"
res, out = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ("work/results", "work/figures")
subprocess.run([sys.executable, str(HERE / "fig1_workflow.py"), out], check=True)
subprocess.run([sys.executable, str(HERE / "fig2_stopped.py"), res, out], check=True)
subprocess.run([sys.executable, str(HERE / "fig_results.py"), res, out], check=True)
