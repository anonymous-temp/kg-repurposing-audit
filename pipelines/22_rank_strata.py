#!/usr/bin/env python
"""Distribution data for Figure 9D. Percentile rank (graph head and MF fitted without write-back, external grid) of unlabelled written-back pairs with a
scientific stop, by number of registered trials before 2015 (same definitions as pipelines/21), plus approved indications."""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit import writeback as L
g = sys.argv[1]                      # hetionet | primekg
ev = L.load_evidence(); tr = ev["trials"]; pre = tr[tr.start < L.CUTOFF_WRITEBACK]; npre = pre.groupby("pair").report_id.nunique()
if g == "hetionet":
    comps, dises, names, ctd, cpd = L.load_graph(); ctd = set(ctd)
    grid = [(c, d) for c in comps for d in dises if (c, d) not in ctd and (c, d) not in cpd]; pos = {p: i for i, p in enumerate(grid)}
    Z = np.load(os.environ.get("E2_SCORES", "work/results/e2_scores.npz")); key = {"graph": "sci=ignore,other=ignore__w10__graph__s1", "mf": "sci=ignore,other=ignore__w10__mf__s1"}
else:
    G = json.load(open(f"{L.DATA}/graph.json")); ctd = {tuple(p) for p in G["indication"]}
    comps = sorted({c for c, d in ctd}); dises = G["diseases"]
    Z = np.load(os.environ.get("E2_SCORES", "work/results_primekg/e2_scores.npz")); pos = {(comps[c], dises[d]): i for i, (c, d) in enumerate(zip(Z["c"], Z["d"]))}
    key = {"graph": "no_writeback__graph", "mf": "no_writeback__mf"}
ext = {p for p in ev["approved"] if p in pos}; y = np.zeros(len(pos)); y[[pos[p] for p in ext]] = 1
out = {}
for m, k in key.items():
    pct = pd.Series(Z[k].astype(float)).rank(pct=True).values * 100
    rows = [(int(npre.get(p, 0)), float(pct[pos[p]]), int(p in ext)) for p in ev["wb_scientific"] if p in pos]
    out[m] = {"stopped": rows, "approved": [float(v) for v in pct[y > 0]]}
    d = pd.DataFrame(rows, columns=["n", "pct", "ext"]); print(g, m, "le2", d[d.n <= 2].pct.median(), "ge10", d[d.n >= 10].pct.median(), "appr", np.median(pct[y > 0]))
OUT = os.environ.get("OUT", "work/results_selection"); os.makedirs(OUT, exist_ok=True)
json.dump(out, open(f"{OUT}/rank_strata_{g}.json", "w"))
