#!/usr/bin/env python
"""Placebo negatives and testing-intensity-restricted negation on Hetionet (stages: e1 <task>, e2, boot <which>).

Every restricting policy of the main analysis can be written as "mask all stopped pairs, then
negate a set S of them". This script asks two questions about S.

1. Placebo negatives. For each negated set S of the main policies, the same number of
   unlabelled pairs without a stopped trial is negated instead, on top of masking:
     - degree-matched placebo: S rewired by double-edge swaps, so that every compound and every
       disease receives as many negatives as under the real policy;
     - uniform placebo: a uniform random sample of the same size.
   Comparing each with masking separates the harm of negating pairs as such (shifted compound and
   disease biases) from the harm of negating the particular pairs that stopped trials select.
2. Testing intensity. A scientific stop is negated only when the pair had at most k registered
   trials (any status) that started before the write-back cut-off; other stopped pairs are masked.
   k = 2 was fixed in advance from the model-free analysis; k = 1 and k = 4 are sensitivity values.

Baseline scores (no write-back, masking and the main policies) are read from the stored score
files of the main analysis (BASE, BASE_ROLE) so that every contrast is paired on the same partitions
and bootstrap draws. Each stage writes its own outputs and is skipped when they exist.
"""
import collections, json, os, sys, time
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit import writeback as L, placebo as P

PART = os.environ.get("PART", "work/partitions")
EMB = os.environ.get("EMB", "work/kge")
BASE = os.environ.get("BASE", "work/results_role")          # main analysis (grid policies, typed_scoped)
BASE_ROLE = os.environ.get("BASE_ROLE", "work/results_role2")  # role-restricted policies
OUT = os.environ.get("OUT", "work/results_extra")
os.makedirs(OUT, exist_ok=True)
SEEDS, PSEEDS, TASKS = [1, 2, 3], [42, 1, 2, 3, 4], ["random", "compound"]
KS = [1, 2, 4]
SIZE_SETS = ["wb_all", "wb_scientific", "wb_scientific_scoped", "wb_scientific_inv", "wb_scientific_inv_scoped", "wb_lowint2"]
SET_POLICY = {"wb_all": "flat_negative", "wb_scientific": "typed", "wb_scientific_scoped": "typed_scoped",
              "wb_scientific_inv": "typed_role", "wb_scientific_inv_scoped": "typed_role_scoped", "wb_lowint2": "typed_lowint2"}

comps, dises, names, ctd, cpd = L.load_graph()
ci = {c: i for i, c in enumerate(comps)}; di = {d: j for j, d in enumerate(dises)}
nC, nD = len(comps), len(dises)
ev = L.load_evidence()
pre = ev["trials"][ev["trials"].start < L.CUTOFF_WRITEBACK]
n_pre = pre.groupby("pair").report_id.nunique()
for k in KS:
    ev[f"wb_lowint{k}"] = {p for p in ev["wb_scientific"] if n_pre.get(p, 0) <= k}
_emb = {}


def emb(ep):
    if ep not in _emb:
        _emb[ep] = L.load_embeddings(f"{EMB}/pretrain_distmult_typed_s1_ep{ep}.pt", comps, dises)
    return _emb[ep]


def part(task, seed):
    p = json.load(open(f"{PART}/{task}_p{seed}.json"))
    for k in ["fit_pos", "val_pos", "test_pos", "eval_grid", "val_grid"]:
        p[k] = [tuple(x) for x in p[k]]
    return p


def masked_rows(p):
    return () if p["task"] != "compound" else {c for c, d in p["test_pos"]} | {c for c, d in p["val_pos"]}


def take(S, pairs):
    return np.array([S[ci[c], di[d]] for c, d in pairs])


def score(model, cfg, Y, Wt, seed):
    if model == "mf":
        return L.fit_scorer(Y, Wt, mf_dim=cfg["dim"], wd=cfg["wd"], lr=0.02, epochs=(cfg["epochs"],), seed=seed)[cfg["epochs"]]
    if model == "graph":
        return L.fit_scorer(Y, Wt, mf_dim=0, emb=emb(cfg["ckpt"]), wd=cfg["wd"], lr=0.01, epochs=(cfg["epochs"],), seed=seed)[cfg["epochs"]]
    return L.fit_scorer(Y, Wt, mf_dim=cfg["dim"], emb=emb(cfg["ckpt"]), wd=cfg["wd"], lr=0.01, epochs=(cfg["epochs"],), seed=seed)[cfg["epochs"]]


def placebo(S_eff, allowed, kind, rng, swaps_per_edge=30):
    """Placebo set of the same size as S_eff (see kg_audit.placebo); returns (pairs, original pairs never replaced)."""
    if kind == "uniform":
        return {(comps[a], dises[b]) for a, b in P.uniform_sample(allowed, len(S_eff), rng)}, 0
    out, left = P.rewire([(ci[c], di[d]) for c, d in S_eff], allowed, rng, swaps_per_edge)
    return {(comps[a], dises[b]) for a, b in out}, left


def allowed_matrix(fit_pos):
    A = np.ones((nC, nD), bool)
    for c, d in set(fit_pos) | cpd | ev["wb_all"]:
        if c in ci and d in di:
            A[ci[c], di[d]] = False
    return A


def configs(fit_pos, seed_base):
    """[(name, actions, info)]: testing-intensity policies and placebo policies."""
    out = []
    for k in KS:
        acts = {p: ("negate" if p in ev[f"wb_lowint{k}"] else "mask") for p in ev["wb_all"]}
        out.append((f"typed_lowint{k}", acts, {}))
    A = allowed_matrix(fit_pos); fs = set(fit_pos)
    for si, sname in enumerate(SIZE_SETS):
        S_eff = {p for p in ev[sname] if p not in fs and p[0] in ci and p[1] in di}
        for ki, kind in enumerate(["degree", "uniform"]):
            rng = np.random.default_rng(seed_base * 100 + si * 10 + ki)
            P, left = placebo(S_eff, A, kind, rng)
            acts = {p: "mask" for p in ev["wb_all"]}
            acts.update({p: "negate" for p in P})
            out.append((f"placebo_{kind}__{SET_POLICY[sname]}", acts, {"placebo_n": len(P), "placebo_unreplaced": left, "real_n": len(S_eff)}))
    return out


# ----------------------------------------------------------------------------- e1
STAGE = sys.argv[1]
if STAGE == "e1":
    task = sys.argv[2]; sel = json.load(open(f"{BASE}/selection.json"))[task]
    fm = f"{OUT}/e1_{task}_metrics.tsv"
    if os.path.exists(fm):
        sys.exit("skip e1 " + task)
    rows, scores = [], {}
    t0 = time.time()
    for ps in PSEEDS:
        p = part(task, ps); test = set(p["test_pos"]); fs = set(p["fit_pos"])
        y = np.array([1.0 if q in test else 0.0 for q in p["eval_grid"]])
        for name, acts, info in configs(p["fit_pos"], ps):
            Y, Wt = L.training_matrix(comps, dises, p["fit_pos"], cpd, acts, w_neg=10.0, masked_rows=masked_rows(p))
            n_aff = {"negated": sum(1 for q, a in acts.items() if a == "negate" and q not in fs),
                     "test_pos_negated": sum(1 for q in test if acts.get(q) == "negate"),
                     "test_pos_masked": sum(1 for q in test if acts.get(q) == "mask"), **info}
            runs = [("degree", 0, L.degree_reference(comps, dises, Y, Wt)), ("graph", 1, score("graph", sel["graph"], Y, Wt, 1))]
            runs += [("mf", sd, score("mf", sel["mf"], Y, Wt, sd)) for sd in SEEDS]
            if name.startswith("typed_lowint"):
                runs += [("hybrid", sd, score("hybrid", sel["hybrid"], Y, Wt, sd)) for sd in SEEDS]
            for model, sd, S in runs:
                s = take(S, p["eval_grid"])
                rows.append({"task": task, "pseed": ps, "policy": name, "model": model, "seed": sd, **L.grid_metrics(p["eval_grid"], y, s), **n_aff})
                scores[f"p{ps}__{name}__{model}__s{sd}"] = s.astype(np.float32)
        print(task, ps, f"{time.time() - t0:.0f}s", flush=True)
        pd.DataFrame(rows).to_csv(fm + ".partial", sep="\t", index=False)
    np.savez_compressed(f"{OUT}/e1_{task}_scores.npz", **scores)
    pd.DataFrame(rows).to_csv(fm, sep="\t", index=False); os.remove(fm + ".partial")

# ----------------------------------------------------------------------------- e2
if STAGE == "e2":
    sel = json.load(open(f"{BASE}/selection.json"))["random"]
    fm = f"{OUT}/e2_metrics.tsv"
    if os.path.exists(fm):
        sys.exit("skip e2")
    ctd_set = set(ctd)
    grid = [(c, d) for c in comps for d in dises if (c, d) not in ctd_set and (c, d) not in cpd]
    ext = sorted(ev["approved"] - ctd_set - cpd); ext_set = set(ext)
    later_fail = sorted(ev["later_scientific"] - ev["wb_all"] - ev["approved"] - ctd_set - cpd)
    y = np.array([1.0 if q in ext_set else 0.0 for q in grid])
    rows, keep = [], {"labels": y.astype(np.int8)}
    for name, acts, info in configs(ctd, 7):
        Y, Wt = L.training_matrix(comps, dises, ctd, cpd, acts, w_neg=10.0)
        runs = [("degree", 0, L.degree_reference(comps, dises, Y, Wt)), ("graph", 1, score("graph", sel["graph"], Y, Wt, 1))]
        runs += [("mf", sd, score("mf", sel["mf"], Y, Wt, sd)) for sd in SEEDS]
        if name.startswith("typed_lowint"):
            runs += [("hybrid", sd, score("hybrid", sel["hybrid"], Y, Wt, sd)) for sd in SEEDS]
        n_neg = sum(1 for q, a in acts.items() if a == "negate" and q not in ctd_set)
        for model, sd, S in runs:
            s = take(S, grid); es, fl = take(S, ext), take(S, later_fail)
            m = L.grid_metrics(grid, y, s)
            m["E3_AUROC_approved_vs_later_failure"] = float(roc_auc_score(np.r_[np.ones(len(es)), np.zeros(len(fl))], np.r_[es, fl]))
            rows.append({"policy": name, "model": model, "seed": sd, "negated": n_neg, "ext_negated": sum(1 for q in ext if acts.get(q) == "negate"), **m, **info})
            keep[f"{name}__{model}__s{sd}"] = s.astype(np.float16)
        print(name, flush=True)
    np.savez_compressed(f"{OUT}/e2_scores.npz", **keep)
    pd.DataFrame(rows).to_csv(fm, sep="\t", index=False)

# ----------------------------------------------------------------------------- boot
BASE_KEYS = {"no_writeback": "sci=ignore,other=ignore", "flat_negative": "sci=negate,other=negate", "mask_all": "sci=mask,other=mask",
             "typed": "sci=negate,other=mask", "typed_scoped": "typed_scoped", "typed_role": "typed_role", "typed_role_scoped": "typed_role_scoped"}


def sorted_ap_factory(y, s):
    o = np.argsort(-s, kind="stable"); yy = y[o]; ss = s[o]
    ends = np.r_[np.flatnonzero(np.diff(ss)), len(ss) - 1]

    def f(w):
        ww = w[o]
        tp = np.cumsum(ww * yy)[ends]; tot = np.cumsum(ww)[ends]
        return float(np.sum(np.diff(np.r_[0.0, tp]) * (tp / np.maximum(tot, 1e-12))) / tp[-1])
    return f


def macro_vals(y, s, dc):
    vals = np.full(nD, np.nan)
    for j in np.unique(dc):
        m = dc == j
        if 0 < y[m].sum() < m.sum():
            vals[j] = L.weighted_ap(y[m], s[m])
    return vals


if STAGE == "boot":
    which = sys.argv[2]; B = int(os.environ.get("BOOT", 1000))
    fo = f"{OUT}/boot_{which}.json"
    if os.path.exists(fo):
        sys.exit("skip boot " + which)
    MODELS = {"degree": [0], "mf": SEEDS, "graph": [1], "hybrid": SEEDS}
    new_names = [f"typed_lowint{k}" for k in KS] + [f"placebo_{kd}__{SET_POLICY[s]}" for s in SIZE_SETS for kd in ("degree", "uniform")]
    dcode = {d: j for j, d in enumerate(dises)}
    Zb = {}
    if which in TASKS:
        Zn = np.load(f"{OUT}/e1_{which}_scores.npz")
        for f in (f"{BASE}/e1_{which}_scores.npz", f"{BASE_ROLE}/e1_{which}_scores.npz"):
            Zb[f] = np.load(f)
        prefixes = [(ps, f"p{ps}__") for ps in PSEEDS]
    else:
        Zn = np.load(f"{OUT}/e2_scores.npz")
        for f in (f"{BASE}/e2_scores.npz", f"{BASE_ROLE}/e2_scores.npz"):
            Zb[f] = np.load(f)
        prefixes = [(0, "")]

    def base_score(prefix, pol, model, sd):
        key = f"{prefix}{BASE_KEYS[pol]}__w10__{model}__s{sd}"
        for Z in Zb.values():
            if key in Z:
                return Z[key].astype(float)
        raise KeyError(key)

    units = []        # (policy, model, ap_fn, macro_vals, disease codes)
    for ps, pre_ in prefixes:
        if which in TASKS:
            p = part(which, ps); dc = np.array([dcode[d] for c, d in p["eval_grid"]]); y = Zb[f"{BASE}/e1_{which}_scores.npz"][f"p{ps}__labels"].astype(float)
        else:
            ctd_set = set(ctd); grid = [(c, d) for c in comps for d in dises if (c, d) not in ctd_set and (c, d) not in cpd]
            dc = np.array([dcode[d] for c, d in grid]); y = Zn["labels"].astype(float)
        for model, seeds in MODELS.items():
            for sd in seeds:
                for pol in BASE_KEYS:
                    s = base_score(pre_, pol, model, sd)
                    units.append((pol, model, sorted_ap_factory(y, s), macro_vals(y, s, dc), dc))
                for name in new_names:
                    key = f"{pre_}{name}__{model}__s{sd}"
                    if key in Zn:
                        s = Zn[key].astype(float)
                        units.append((name, model, sorted_ap_factory(y, s), macro_vals(y, s, dc), dc))

    def stat(wd):
        acc = collections.defaultdict(list)
        for pol, model, f, mv, dcs in units:
            ok = ~np.isnan(mv)
            acc[(pol, model)].append((f(wd[dcs]), float(np.sum(mv[ok] * wd[ok]) / max(np.sum(wd[ok]), 1e-12))))
        return {k: tuple(np.mean(v, axis=0)) for k, v in acc.items()}

    obs = stat(np.ones(nD))
    present = np.unique(np.concatenate([u[4] for u in units]))
    rng = np.random.default_rng(20261006)
    draws = []
    for b in range(B):
        pick = present[rng.integers(len(present), size=len(present))]
        draws.append(stat(np.bincount(pick, minlength=nD).astype(float)))
    out = {"B": B, "cluster": "disease", "observed": {f"{k[0]}|{k[1]}": {"pooledAP": v[0], "macroAP": v[1]} for k, v in obs.items()}, "contrasts": {}}
    pols = list(BASE_KEYS) + new_names
    for model in MODELS:
        for a in pols:
            for b in ("no_writeback", "mask_all"):
                if a == b or (a, model) not in obs or (b, model) not in obs:
                    continue
                rec = {}
                for i, metric in enumerate(["pooledAP", "macroAP"]):
                    dd = np.array([dr[(a, model)][i] - dr[(b, model)][i] for dr in draws])
                    rec[metric] = {"diff": obs[(a, model)][i] - obs[(b, model)][i], "lo": float(np.percentile(dd, 2.5)), "hi": float(np.percentile(dd, 97.5))}
                out["contrasts"][f"{a} - {b} [{model}]"] = rec
    json.dump(out, open(fo, "w"), indent=1)
    print(which, "contrasts", len(out["contrasts"]))
