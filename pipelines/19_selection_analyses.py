#!/usr/bin/env python
"""Model-free tests of the selection explanation, for one graph (argument: hetionet | primekg).

  intensity : logistic models of "pair is a recorded or approved indication" among tested pairs, with and
              without adjustment for testing intensity; paired disease-cluster bootstrap of the odds ratios
  purity    : share of indications among stopped pairs by number of registered trials and by trial phase
  null      : degree-preserving permutation null for the enrichment of external approved indications
              among unlabelled stopped pairs
  history   : scorers built only from a pair's registry history (trial count, any stop) on the external grid
  disease   : (Hetionet) per-disease loss under flat negatives against the share of held-out treatments negated

Writes <OUT>/selection_<graph>.json and TSV tables; every block is skipped if its key exists already.
"""
import json, os, re, sys, time
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit import writeback as L, fastboot as F

GRAPH = sys.argv[1]
OUT = os.environ.get("OUT", "work/results_selection"); os.makedirs(OUT, exist_ok=True)
B = int(os.environ.get("BOOT", 1000)); NPERM = int(os.environ.get("NPERM", 1000))
fo = f"{OUT}/selection_{GRAPH}.json"
res = json.load(open(fo)) if os.path.exists(fo) else {}
ONC = re.compile(r"cancer|carcinoma|neoplasm|lymphoma|leuk|melanoma|sarcoma|glioma|myeloma|tumou?r|blastoma|mesothelioma", re.I)

if GRAPH == "hetionet":
    comps, dises, names, ctd, cpd = L.load_graph(); ctd = set(ctd)
else:
    G = json.load(open(f"{L.DATA}/graph.json"))
    ctd = {tuple(p) for p in G["indication"]}; cpd = {tuple(p) for p in G["off_label"]} - ctd
    comps = sorted({c for c, d in ctd}); dises = G["diseases"]; names = G["names"]
    cpd = {(c, d) for c, d in cpd if c in set(comps) and d in set(dises)}
ci = {c: i for i, c in enumerate(comps)}; di = {d: j for j, d in enumerate(dises)}; nC, nD = len(comps), len(dises)
ev = L.load_evidence()
ind = ctd | ev["approved"]
onc = {d: bool(ONC.search(str(names.get(d, "")))) for d in dises}
inside = lambda S: {q for q in S if q[0] in ci and q[1] in di}


def save():
    json.dump(res, open(fo, "w"), indent=1)


def phase_of(x):
    x = str(x)
    return "4" if "PHASE4" in x else "3" if "PHASE3" in x else "2" if "PHASE2" in x else "1" if "PHASE1" in x else "NA"


def wilson(k, n, z=1.96):
    if n == 0:
        return [float("nan")] * 2
    p = k / n; den = 1 + z * z / n; c = p + z * z / (2 * n); h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return [float((c - h) / den), float((c + h) / den)]


# ----------------------------------------------------------------------------- tested pairs
tr = ev["trials"]
tr = tr[(tr.start < L.CUTOFF_WRITEBACK) & tr.phase.notna() & ~tr.phase.fillna("").str.contains("PHASE4")]
tr = tr[[(p in inside({p})) and p not in cpd for p in tr.pair]].copy()
cats = tr.stop_categories.fillna("")
tr["sci"] = tr.status.isin(L.STOP) & cats.apply(lambda s: bool(set(s.split("|")) & L.SCIENTIFIC))
tr["stop"] = tr.status.isin(L.STOP)
tr["early"] = ~tr.phase.str.contains("PHASE3")
g = tr.groupby("pair")
D = pd.DataFrame({"n": g.report_id.nunique(), "sci": g.sci.any().astype(int), "stop": g.stop.any().astype(int)})
D["n12"] = tr[tr.early].groupby("pair").report_id.nunique().reindex(D.index).fillna(0)
D["ind"] = [int(p in ind) for p in D.index]
D["dis"] = [di[p[1]] for p in D.index]; D["onc"] = [int(onc[p[1]]) for p in D.index]
ncd = D.groupby(D.index.map(lambda p: p[0])).size(); ndc = D.groupby(D.index.map(lambda p: p[1])).size()
D["ln"] = np.log(D.n); D["ln12"] = np.log1p(D.n12)
D["lnc"] = np.log([ncd[p[0]] for p in D.index]); D["lnd"] = np.log([ndc[p[1]] for p in D.index])

SPECS = {"crude": ("all", ["sci"]), "adjusted": ("all", ["ln", "sci"]), "adjusted_any_stop": ("all", ["ln", "sci", "stop"]),
         "adjusted_popularity": ("all", ["ln", "lnc", "lnd", "sci"]), "adjusted_phase12_intensity": ("all", ["ln12", "sci"]),
         "adjusted_oncology": ("onc", ["ln", "sci"]), "adjusted_non_oncology": ("non", ["ln", "sci"]),
         "crude_any_stop": ("all", ["stop"]), "adjusted_any_stop_only": ("all", ["ln", "stop"])}


def logit_or(d, cols, w=None, term=None):
    X = sm.add_constant(d[cols].astype(float)); m = sm.GLM(d.ind.astype(float), X, family=sm.families.Binomial(), freq_weights=w).fit()
    return {c: float(np.exp(m.params[c])) for c in cols}


if "intensity" not in res:
    t0 = time.time(); rng = np.random.default_rng(20261007)
    out = {"n_tested_pairs": int(len(D)), "indication_rate": float(D.ind.mean()), "models": {}}
    subsets = {"all": D, "onc": D[D.onc == 1], "non": D[D.onc == 0]}
    present = np.unique(D.dis)
    drawsW = [np.bincount(present[rng.integers(len(present), size=len(present))], minlength=nD) for _ in range(B)]
    for name, (sub, cols) in SPECS.items():
        d = subsets[sub]; est = logit_or(d, cols); bs = {c: [] for c in cols}
        for wd in drawsW:
            w = wd[d.dis.values].astype(float); keep = w > 0
            try:
                e = logit_or(d[keep], cols, w[keep])
            except Exception:
                continue
            for c in cols:
                bs[c].append(e[c])
        out["models"][name] = {"n": int(len(d)), "terms": {c: {"OR": est[c], "lo": float(np.percentile(bs[c], 2.5)), "hi": float(np.percentile(bs[c], 97.5))} for c in cols}}
        print(GRAPH, name, {c: round(est[c], 2) for c in cols}, f"{time.time() - t0:.0f}s", flush=True)
    bins = [0, 1, 2, 4, 9, 19, 10 ** 6]; labels = ["1", "2", "3-4", "5-9", "10-19", "20+"]
    D["bin"] = pd.cut(D.n, bins, labels=labels)
    rows = []
    for b in labels:
        for s in (0, 1):
            x = D[(D.bin == b) & (D.sci == s)]
            rows.append({"bin": b, "scientific_stop": s, "pairs": int(len(x)), "indications": int(x.ind.sum()),
                         "share": float(x.ind.mean()) if len(x) else float("nan"), "ci": wilson(int(x.ind.sum()), len(x))})
    out["by_trial_count"] = rows
    res["intensity"] = out; save()

# ----------------------------------------------------------------------------- purity of the write-back set
if "purity" not in res:
    allt = ev["trials"]; allt = allt[allt.start < L.CUTOFF_WRITEBACK]
    npre = allt.groupby("pair").report_id.nunique()
    out = {}
    for setname in ["wb_all", "wb_scientific", "wb_scientific_inv_scoped"]:
        S = [p for p in inside(ev[setname]) if p not in cpd]
        w = pd.DataFrame({"n": [npre.get(p, 0) for p in S], "ind": [int(p in ind) for p in S], "onc": [onc[p[1]] for p in S]})
        w["bin"] = pd.cut(w.n, [0, 1, 2, 4, 9, 19, 10 ** 6], labels=["1", "2", "3-4", "5-9", "10-19", "20+"])
        out[setname] = {"by_trial_count": [{"bin": b, "pairs": int((w.bin == b).sum()), "indications": int(w[w.bin == b].ind.sum()),
                                            "share": float(w[w.bin == b].ind.mean()) if (w.bin == b).any() else float("nan"),
                                            "ci": wilson(int(w[w.bin == b].ind.sum()), int((w.bin == b).sum()))} for b in w.bin.cat.categories],
                        "le": {str(k): {"pairs": int((w.n <= k).sum()), "indications": int(w[w.n <= k].ind.sum())} for k in (1, 2, 4)},
                        "oncology": {str(o): {"pairs": int((w.onc == o).sum()), "indications": int(w[w.onc == o].ind.sum())} for o in (True, False)}}
    st = ev["stopped"]; st = st[(st.start < L.CUTOFF_WRITEBACK)]
    sc = st[st.stop_categories.fillna("").apply(lambda s: bool(set(s.split("|")) & L.SCIENTIFIC))].copy()
    sc = sc[[p in inside({p}) and p not in cpd for p in sc.pair]]
    sc["ph"] = sc.phase.map(phase_of)
    ph = []
    for k in ["1", "2", "3", "4"]:
        P = set(sc[sc.ph == k].pair); n_i = sum(p in ind for p in P)
        ph.append({"phase": k, "pairs": len(P), "indications": int(n_i), "share": n_i / len(P) if P else float("nan"), "ci": wilson(n_i, len(P))})
    out["scientific_by_phase"] = ph
    res["purity"] = out; save(); print(GRAPH, "purity done", flush=True)

# ----------------------------------------------------------------------------- degree-preserving null
NULL_SETS = ["wb_all", "wb_scientific", "wb_scientific_scoped", "wb_scientific_inv", "wb_scientific_inv_scoped", "wb_other"]
if any(k not in res.get("null", {}) for k in NULL_SETS):
    t0 = time.time()
    excl = np.zeros((nC, nD), bool)
    for c, d in inside(ctd | cpd):
        excl[ci[c], di[d]] = True
    ext = inside(ev["approved"]) - ctd - cpd
    X = np.zeros((nC, nD), bool)
    for c, d in ext:
        X[ci[c], di[d]] = True
    n_grid = int((~excl).sum()); out = res.get("null", {})
    for setname in [k for k in NULL_SETS if k not in out]:
        S = [(ci[c], di[d]) for c, d in inside(ev[setname]) if not excl[ci[c], di[d]]]
        obs = int(sum(X[a, b] for a, b in S)); n = len(S)
        base = (X.sum() - obs) / (n_grid - n)
        rng = np.random.default_rng(11); null = []
        for r in range(NPERM):
            E = list(S); cur = np.zeros((nC, nD), bool)
            for a, b in E:
                cur[a, b] = True
            for _ in range(10 * n):
                i, j = rng.integers(n, size=2); (a, b), (c, d) = E[i], E[j]
                if a == c or b == d or cur[a, d] or cur[c, b] or excl[a, d] or excl[c, b]:
                    continue
                cur[a, b] = cur[c, d] = False; cur[a, d] = cur[c, b] = True; E[i], E[j] = (a, d), (c, b)
            null.append(int(sum(X[a, b] for a, b in E)))
        null = np.array(null)
        out[setname] = {"pairs": n, "observed": obs, "null_mean": float(null.mean()), "null_lo": float(np.percentile(null, 2.5)), "null_hi": float(np.percentile(null, 97.5)),
                        "base_rate": float(base), "fold_observed": float(obs / n / base), "fold_null": float(null.mean() / n / base),
                        "ratio_observed_to_null": float(obs / null.mean()), "ratio_lo": float(obs / np.percentile(null, 97.5)), "ratio_hi": float(obs / max(np.percentile(null, 2.5), 1e-9)),
                        "p_one_sided": float((np.sum(null >= obs) + 1) / (len(null) + 1))}
        print(GRAPH, setname, out[setname], f"{time.time() - t0:.0f}s", flush=True)
    res["null"] = out; save()

# ----------------------------------------------------------------------------- registry-history scorers
if "history" not in res:
    from sklearn.metrics import roc_auc_score
    gc, gd = np.nonzero(np.ones((nC, nD), bool))
    excl = np.zeros((nC, nD), bool)
    for c, d in inside(ctd | cpd):
        excl[ci[c], di[d]] = True
    keep = ~excl[gc, gd]; gc, gd = gc[keep], gd[keep]
    ext = inside(ev["approved"]) - ctd - cpd
    later_fail = inside(ev["later_scientific"]) - ev["wb_all"] - ev["approved"] - ctd - cpd
    def mat(S, vals=None):
        M = np.zeros((nC, nD))
        for i, (c, d) in enumerate(S):
            if c in ci and d in di:
                M[ci[c], di[d]] = 1.0 if vals is None else vals[i]
        return M
    y = mat(ext)[gc, gd]; fl = mat(later_fail)[gc, gd]
    allt = ev["trials"]; allt = allt[allt.start < L.CUTOFF_WRITEBACK]
    cnt = allt.groupby("pair").report_id.nunique(); scnt = allt[allt.status.isin(L.STOP)].groupby("pair").report_id.nunique()
    scorers = {"trial_count": mat(list(cnt.index), cnt.values)[gc, gd], "stopped_trial_count": mat(list(scnt.index), scnt.values)[gc, gd],
               "any_stop": mat(inside(ev["wb_all"]))[gc, gd]}
    rng = np.random.default_rng(20261007)
    W = np.stack([np.bincount(rng.integers(nD, size=nD), minlength=nD) for _ in range(B)]).astype(np.float32)
    out = {"n_grid": int(len(y)), "n_external": int(y.sum()), "n_later_failures": int(fl.sum())}
    for k, s in scorers.items():
        vals = F.per_disease_ap(y, s, gd, nD, L.weighted_ap); pa = F.pooled_ap_draws(y, s, gd, W); ma = F.macro_draws(vals, W)
        es, fs = s[y > 0], s[fl > 0]
        out[k] = {"pooledAP": L.weighted_ap(y, s), "pooledAP_ci": [float(np.nanpercentile(pa, 2.5)), float(np.nanpercentile(pa, 97.5))],
                  "macroAP": float(np.nanmean(vals)), "macroAP_ci": [float(np.nanpercentile(ma, 2.5)), float(np.nanpercentile(ma, 97.5))],
                  "E3_AUROC": float(roc_auc_score(np.r_[np.ones(len(es)), np.zeros(len(fs))], np.r_[es, fs])),
                  "share_external_with_any_pre2015_trial": float((s[y > 0] > 0).mean())}
        print(GRAPH, k, {a: out[k][a] for a in ("pooledAP", "macroAP", "E3_AUROC")}, flush=True)
    res["history"] = out; save()

# ----------------------------------------------------------------------------- per-disease heterogeneity (Hetionet)
if GRAPH == "hetionet" and "disease" not in res:
    chg = pd.read_csv(os.environ.get("DIS_CHANGES", "work/results_role/e1_disease_changes.tsv"), sep="\t")
    PART = os.environ.get("PART", "work/partitions"); out = {}
    for task in ("random", "compound"):
        tot, neg = {}, {}
        for ps in (42, 1, 2, 3, 4):
            p = json.load(open(f"{PART}/{task}_p{ps}.json"))
            for c, d in map(tuple, p["test_pos"]):
                tot[d] = tot.get(d, 0) + 1; neg[d] = neg.get(d, 0) + int((c, d) in ev["wb_all"])
        x = chg[chg.task == task].copy(); x["share_negated"] = [neg.get(d, 0) / tot[d] if tot.get(d) else np.nan for d in x.disease]
        x = x.dropna(subset=["share_negated"])
        r = spearmanr(x.share_negated, x.flat_minus_none)
        hi = x[x.share_negated >= 0.5]; lo = x[x.share_negated < 0.5]
        out[task] = {"diseases": int(len(x)), "spearman_rho": float(r.statistic), "p": float(r.pvalue),
                     "median_change_share_ge_half": float(hi.flat_minus_none.median()), "n_ge_half": int(len(hi)),
                     "median_change_share_lt_half": float(lo.flat_minus_none.median()), "n_lt_half": int(len(lo)),
                     "oncology_median_change": float(x[x.disease.map(onc)].flat_minus_none.median()), "non_oncology_median_change": float(x[~x.disease.map(onc)].flat_minus_none.median())}
        x.to_csv(f"{OUT}/disease_changes_{task}.tsv", sep="\t", index=False)
    res["disease"] = out; save(); print(out)
print(GRAPH, "all blocks done")
