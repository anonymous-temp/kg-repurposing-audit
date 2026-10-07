#!/usr/bin/env python
"""Supporting numbers for the testing-intensity analyses (one graph per call: hetionet | primekg).

- percentile rank (no write-back, external grid) of unlabelled written-back pairs with a scientific stop, by testing intensity;
- per-disease AP on external approved indications with disease-cluster bootstrap intervals for the scorers fitted without write-back;
- share of external approved indications and of other unlabelled pairs with any registered trial before the cut-off.
Inputs: E2 score files of the main analyses (E2_SCORES, E2_BOOT for PrimeKG). Output: OUT/rank_by_intensity_<g>.json, e2_scorers_<g>.json."""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit import writeback as L, fastboot as F
g = sys.argv[1]; OUT = os.environ.get("OUT", "work/results_selection"); os.makedirs(OUT, exist_ok=True); B = 1000
ev = L.load_evidence(); tr = ev["trials"]; pre = tr[tr.start < L.CUTOFF_WRITEBACK]; npre = pre.groupby("pair").report_id.nunique()
if g == "hetionet":
    comps, dises, names, ctd, cpd = L.load_graph(); ctd = set(ctd)
    grid = [(c, d) for c in comps for d in dises if (c, d) not in ctd and (c, d) not in cpd]; pos = {p: i for i, p in enumerate(grid)}
    di = {d: j for j, d in enumerate(dises)}; gd = np.array([di[d] for c, d in grid])
    Z = np.load(os.environ.get("E2_SCORES", "work/results/e2_scores.npz"))
    key = {"graph": ["sci=ignore,other=ignore__w10__graph__s1"], "mf": [f"sci=ignore,other=ignore__w10__mf__s{s}" for s in (1, 2, 3)],
           "degree": ["sci=ignore,other=ignore__w10__degree__s0"], "hybrid": [f"sci=ignore,other=ignore__w10__hybrid__s{s}" for s in (1, 2, 3)]}
else:
    G = json.load(open(f"{L.DATA}/graph.json")); ctd = {tuple(p) for p in G["indication"]}; cpd = {tuple(p) for p in G["off_label"]} - ctd
    comps = sorted({c for c, d in ctd}); dises = G["diseases"]
    Z = np.load(os.environ.get("E2_SCORES", "work/results_primekg/e2_scores.npz")); gd = Z["d"].astype(int)
    pos = {(comps[c], dises[d]): i for i, (c, d) in enumerate(zip(Z["c"], Z["d"]))}
    key = {m: [f"no_writeback__{m}"] for m in ("graph", "mf", "degree")}
nD = len(dises); ext = {p for p in ev["approved"] if p in pos}
y = np.zeros(len(pos)); y[[pos[p] for p in ext]] = 1
rank = {}
for m, ks in key.items():
    pct = pd.Series(Z[ks[0]].astype(float)).rank(pct=True).values * 100
    d = pd.DataFrame([(npre.get(p, 0), pct[pos[p]]) for p in ev["wb_scientific"] if p in pos], columns=["n", "pct"])
    rank[m] = {"le2_median_pct": float(d[d.n <= 2].pct.median()), "le2_n": int((d.n <= 2).sum()), "ge10_median_pct": float(d[d.n >= 10].pct.median()),
               "ge10_n": int((d.n >= 10).sum()), "approved_median_pct": float(np.median(pct[y > 0]))}
json.dump(rank, open(f"{OUT}/rank_by_intensity_{g}.json", "w"), indent=1)
rng = np.random.default_rng(20261007); W = np.stack([np.bincount(rng.integers(nD, size=nD), minlength=nD) for _ in range(B)]).astype(np.float32)
sc = {}
if g == "primekg":                      # inline bootstrap of the original fit (float64 scores)
    Bz = np.load(os.environ.get("E2_BOOT", "work/results_primekg/e2_boot.npz"))
    for m in key:
        dr = Bz[f"no_writeback__{m}"][1]; sc[m] = {"macroAP": float(Bz[f"no_writeback__{m}__obs"][1]), "ci": [float(np.nanpercentile(dr, 2.5)), float(np.nanpercentile(dr, 97.5))]}
else:
    for m, ks in key.items():
        vals = np.mean([F.per_disease_ap(y, Z[k].astype(float), gd, nD, L.weighted_ap) for k in ks], 0); dr = F.macro_draws(vals, W)
        sc[m] = {"macroAP": float(np.nanmean(vals)), "ci": [float(np.nanpercentile(dr, 2.5)), float(np.nanpercentile(dr, 97.5))]}
tested = {p for p in set(pre.pair) if p in pos}
sc["share_tested"] = {"external": len(ext & tested) / len(ext), "other_unlabelled": len(tested - ext) / (len(pos) - len(ext))}
json.dump(sc, open(f"{OUT}/e2_scorers_{g}.json", "w"), indent=1); print(g, rank, sc)
