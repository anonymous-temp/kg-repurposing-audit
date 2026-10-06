#!/usr/bin/env python
"""Do no-write-back scorers rank pairs that were later tested as highly as pairs that were later approved?

Categories of unlabelled pairs on the external grid (E2 grid):
  A  external approved indication (E2 positive)
  F  later scientific failure (E3: efficacy/safety stop of a trial started 2017+, no earlier stop, not approved)
  T  first clinical trial started 2017 or later, not approved, not F
  N  never in any Open Targets clinical report
AUROCs between categories with a paired disease-cluster bootstrap (B = 1000).
Usage: testedness.py hetionet|primekg
"""
import json, os, sys
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
G = sys.argv[1]
OUTD = "work/results_testedness"
if os.path.exists(f"{OUTD}/{G}.json"):
    sys.exit("skip " + G)
os.environ["DATA"] = "work/mapped" if G == "hetionet" else "work/mapped_primekg"
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit import writeback as L
ev = L.load_evidence()
tp = pd.read_csv(f"{L.DATA}/trial_pairs.tsv", sep="\t", usecols=["compound", "disease"], low_memory=False)
anyrep = set(zip(tp.compound, tp.disease))
first = ev["trials"].groupby("pair").start.min()
later_first = {p for p, t in first.items() if pd.notna(t) and t >= pd.Timestamp("2017-01-01")}
if G == "hetionet":
    comps, dises, names, ctd, cpd = L.load_graph()
    ctd_set = set(ctd)
    grid = [(c, d) for c in comps for d in dises if (c, d) not in ctd_set and (c, d) not in cpd]
    Z = np.load("work/results/e2_scores.npz"); y = Z["labels"].astype(int)
    fail = set(map(tuple, pd.read_csv("work/results/e3_later_failures.tsv", sep="\t").values))
    S = {m: np.mean([Z[k].astype(float) for k in Z.files if k.startswith(f"sci=ignore,other=ignore__w10__{m}__")], 0) for m in ["degree", "mf", "graph"]}
    dis = np.array([d for c, d in grid]); pairs = grid
else:   # PrimeKG: integer pair codes, so that 2.4 million pairs fit in memory
    Gj = json.load(open(f"{L.DATA}/graph.json")); comps = sorted({c for c, d in Gj["indication"]}); dises = Gj["diseases"]
    ci = {c: i for i, c in enumerate(comps)}; di = {d: j for j, d in enumerate(dises)}; nD = len(dises)
    Z = np.load("work/results_primekg/e2_scores.npz"); y = Z["labels"].astype(int)
    code = Z["c"].astype(np.int64) * nD + Z["d"].astype(np.int64)
    enc = lambda P: np.array(sorted(ci[c] * nD + di[d] for c, d in P if c in ci and d in di), np.int64)
    is_fail = Z["later_failure"].astype(bool)
    is_rep = np.isin(code, enc(anyrep)); is_new = np.isin(code, enc(later_first))
    lab = np.full(len(code), "", dtype=object)
    lab[~is_rep] = "N"; lab[is_new] = "T"; lab[is_fail] = "F"; lab[y == 1] = "A"
    S = {m: Z[f"no_writeback__{m}"].astype(float) for m in ["degree", "mf", "graph"]}
    dis = Z["d"].astype(np.int64)
if G == "hetionet":
    lab = np.full(len(pairs), "", dtype=object)
    for i, q in enumerate(pairs):
        if y[i]: lab[i] = "A"
        elif q in fail: lab[i] = "F"
        elif q in later_first: lab[i] = "T"
        elif q not in anyrep: lab[i] = "N"
ud, dcode = np.unique(dis, return_inverse=True)
rng = np.random.default_rng(20261006)
Wd = np.stack([np.bincount(rng.integers(len(ud), size=len(ud)), minlength=len(ud)) for _ in range(1000)])
idx = {k: np.flatnonzero(lab == k) for k in "AFTN"}


def wauc(s, a, b, wa, wb):
    """Weighted AUROC (weights = disease multiplicities)."""
    x = np.r_[s[a], s[b]]; yy = np.r_[np.ones(len(a)), np.zeros(len(b))]; w = np.r_[wa, wb]
    if w[yy == 1].sum() == 0 or w[yy == 0].sum() == 0:
        return np.nan
    return roc_auc_score(yy, x, sample_weight=w)


out = {"n": {k: int(len(v)) for k, v in idx.items()}, "auroc": {}, "median_pct": {}}
for m, s in S.items():
    pct = pd.Series(s).rank(pct=True).values
    out["median_pct"][m] = {k: float(np.median(pct[v])) for k, v in idx.items()}
    for a, b in [("A", "N"), ("F", "N"), ("T", "N"), ("A", "T"), ("A", "F")]:
        ia, ib = idx[a], idx[b]
        if b == "N" and len(ib) > 200000:     # subsample the never-tested reference for the bootstrap (fixed seed)
            ib = np.sort(np.random.default_rng(1).choice(ib, 200000, replace=False))
        est = wauc(s, ia, ib, np.ones(len(ia)), np.ones(len(ib)))
        bs = [wauc(s, ia, ib, Wd[r][dcode[ia]], Wd[r][dcode[ib]]) for r in range(1000)]
        lo, hi = np.nanpercentile(bs, [2.5, 97.5])
        out["auroc"][f"{m}|{a}_vs_{b}"] = {"est": float(est), "lo": float(lo), "hi": float(hi)}
        print(G, m, a, "vs", b, round(est, 3), round(lo, 3), round(hi, 3), flush=True)
json.dump(out, open(f"{OUTD}/{G}.json", "w"), indent=1)
