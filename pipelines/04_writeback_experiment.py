#!/usr/bin/env python
"""Failure write-back experiment on Hetionet (stages: select, e1, e2, boot).

  select : choose hyper-parameters of the learned scorers on validation AP
           (no write-back), separately for the random-edge and compound-disjoint tasks
  e1     : held-out recorded indications (10 partitions) under every write-back policy
  e2     : models fitted on all recorded indications; external approved indications
           absent from Hetionet (E2) and later scientific failures (E3)
  boot   : paired disease-cluster bootstrap of policy and model contrasts

Every stage writes its own output file and is skipped when that file exists.
"""
import collections, json, os, sys, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit import writeback as L

PART = os.environ.get("PART", "work/partitions")
EMB = os.environ.get("EMB", "work/kge")
OUT = os.environ.get("OUT", "work/results")
EMB_TAG = os.environ.get("EMB_TAG", "pretrain_distmult_typed_s1")
os.makedirs(OUT, exist_ok=True)
STAGE = sys.argv[1]
SEEDS = [1, 2, 3]
PSEEDS = [42, 1, 2, 3, 4]
TASKS = ["random", "compound"]

comps, dises, names, ctd, cpd = L.load_graph()
ci = {c: i for i, c in enumerate(comps)}; di = {d: j for j, d in enumerate(dises)}
ev = L.load_evidence()
_emb_cache = {}


def emb(ep):
    if ep not in _emb_cache:
        _emb_cache[ep] = L.load_embeddings(f"{EMB}/{EMB_TAG}_ep{ep}.pt", comps, dises)
    return _emb_cache[ep]


def part(task, seed):
    p = json.load(open(f"{PART}/{task}_p{seed}.json"))
    for k in ["fit_pos", "val_pos", "test_pos", "eval_grid", "val_grid"]:
        p[k] = [tuple(x) for x in p[k]]
    return p


def masked_rows(p):
    if p["task"] != "compound":
        return ()
    return {c for c, d in p["test_pos"]} | {c for c, d in p["val_pos"]}


def take(S, pairs):
    return np.array([S[ci[c], di[d]] for c, d in pairs])


def score_model(model, cfg, Y, Wt, seed):
    """Return {epoch: score matrix} for one learned scorer configuration."""
    if model == "mf":
        return L.fit_scorer(Y, Wt, mf_dim=cfg["dim"], wd=cfg["wd"], lr=0.02, epochs=tuple(cfg["epochs"]), seed=seed)
    if model == "graph":
        return L.fit_scorer(Y, Wt, mf_dim=0, emb=emb(cfg["ckpt"]), wd=cfg["wd"], lr=0.01, epochs=tuple(cfg["epochs"]), seed=seed)
    if model == "hybrid":
        return L.fit_scorer(Y, Wt, mf_dim=cfg["dim"], emb=emb(cfg["ckpt"]), wd=cfg["wd"], lr=0.01, epochs=tuple(cfg["epochs"]), seed=seed)
    raise ValueError(model)


def no_wb_matrix(p):
    acts = L.policy_sets(ev, L.NAMED["no_writeback"])
    return L.training_matrix(comps, dises, p["fit_pos"], cpd, acts, masked_rows=masked_rows(p))


# ----------------------------------------------------------------------------- select
if STAGE == "select":
    f = f"{OUT}/selection.json"
    if os.path.exists(f):
        print("skip select"); sys.exit(0)
    EPOCHS = [200, 400, 800, 1600]
    rows = []
    t0 = time.time()
    for task in TASKS:
        for ps in PSEEDS:
            p = part(task, ps)
            Y, Wt = no_wb_matrix(p)
            vy = np.array([1.0 if q in set(p["val_pos"]) else 0.0 for q in p["val_grid"]])
            cfgs = [("mf", {"dim": d, "wd": w}) for d in (16, 32) for w in (1e-5, 1e-4, 1e-3)]
            cfgs += [("graph", {"ckpt": e, "wd": w}) for e in (5, 10, 15, 20) for w in (1e-5, 1e-4, 1e-3)]
            for model, cfg in cfgs:
                out = score_model(model, {**cfg, "epochs": EPOCHS}, Y, Wt, seed=1)
                for ep, S in out.items():
                    rows.append({"task": task, "pseed": ps, "model": model, **cfg, "epochs": ep,
                                 "val_AP": L.weighted_ap(vy, take(S, p["val_grid"]))})
            print(task, ps, f"{time.time() - t0:.0f}s", flush=True)
    df = pd.DataFrame(rows)
    sel = {}
    for task in TASKS:
        sel[task] = {}
        for model in ["mf", "graph"]:
            g = df[(df.task == task) & (df.model == model)]
            keys = ["dim", "wd", "epochs"] if model == "mf" else ["ckpt", "wd", "epochs"]
            m = g.groupby(keys).val_AP.mean().reset_index().sort_values("val_AP", ascending=False)
            best = m.iloc[0].to_dict()
            sel[task][model] = {k: (int(best[k]) if k in ("dim", "ckpt", "epochs") else float(best[k])) for k in keys}
            sel[task][model]["val_AP"] = float(best["val_AP"])
    # hybrid: MF dimension and embedding checkpoint from the two selections, weight decay re-selected
    rows2 = []
    for task in TASKS:
        for ps in PSEEDS:
            p = part(task, ps)
            Y, Wt = no_wb_matrix(p)
            vy = np.array([1.0 if q in set(p["val_pos"]) else 0.0 for q in p["val_grid"]])
            for w in (1e-5, 1e-4, 1e-3):
                cfg = {"dim": sel[task]["mf"]["dim"], "ckpt": sel[task]["graph"]["ckpt"], "wd": w}
                out = score_model("hybrid", {**cfg, "epochs": EPOCHS}, Y, Wt, seed=1)
                for ep, S in out.items():
                    rows2.append({"task": task, "pseed": ps, "model": "hybrid", **cfg, "epochs": ep,
                                  "val_AP": L.weighted_ap(vy, take(S, p["val_grid"]))})
    df2 = pd.DataFrame(rows2)
    for task in TASKS:
        m = df2[df2.task == task].groupby(["dim", "ckpt", "wd", "epochs"]).val_AP.mean().reset_index().sort_values("val_AP", ascending=False)
        best = m.iloc[0].to_dict()
        sel[task]["hybrid"] = {"dim": int(best["dim"]), "ckpt": int(best["ckpt"]), "wd": float(best["wd"]), "epochs": int(best["epochs"]),
                               "val_AP": float(best["val_AP"])}
    pd.concat([df, df2]).to_csv(f"{OUT}/selection_grid.tsv", sep="\t", index=False)
    json.dump(sel, open(f, "w"), indent=1)
    print(json.dumps(sel, indent=1))

# ----------------------------------------------------------------------------- e1
POLICIES = L.GRID + ["typed_scoped"]
WSENS = [("flat_negative", 3.0), ("flat_negative", 30.0), ("typed", 3.0), ("typed", 30.0)]


def policy_list():
    if os.environ.get("MAIN_ONLY") == "1":                     # sensitivity analyses: five main policies, weight 10
        return [(L.NAMED[n], 10.0) for n in ["no_writeback", "flat_negative", "mask_all", "typed", "typed_scoped"]]
    out = [(pol, 10.0) for pol in POLICIES]
    out += [(L.NAMED[n], w) for n, w in WSENS]
    return out


if STAGE == "e1":
    sel = json.load(open(f"{OUT}/selection.json"))
    task = sys.argv[2]
    fm = f"{OUT}/e1_{task}_metrics.tsv"
    if os.path.exists(fm):
        print("skip e1", task); sys.exit(0)
    rows, scores = [], {}
    t0 = time.time()
    for ps in PSEEDS:
        p = part(task, ps)
        test = set(p["test_pos"])
        y = np.array([1.0 if q in test else 0.0 for q in p["eval_grid"]])
        scores[f"p{ps}__labels"] = y.astype(np.int8)
        for pol, w in policy_list():
            acts = L.policy_sets(ev, pol)
            Y, Wt = L.training_matrix(comps, dises, p["fit_pos"], cpd, acts, w_neg=w, masked_rows=masked_rows(p))
            n_aff = {"negated": sum(1 for q, a in acts.items() if a == "negate" and q not in set(p["fit_pos"])),
                     "masked": sum(1 for q, a in acts.items() if a == "mask" and q not in set(p["fit_pos"])),
                     "test_pos_negated": sum(1 for q in test if acts.get(q) == "negate"),
                     "test_pos_masked": sum(1 for q in test if acts.get(q) == "mask")}
            runs = [("degree", 0, L.degree_reference(comps, dises, Y, Wt))]
            if pol == L.NAMED["no_writeback"] and w == 10.0:
                runs.append(("disease_degree", 0, L.disease_degree(Y)))
            for model in ["mf", "graph", "hybrid"]:
                cfg = sel[task][model]
                for sd in SEEDS:
                    S = score_model(model, {**cfg, "epochs": [cfg["epochs"]]}, Y, Wt, seed=sd)[cfg["epochs"]]
                    runs.append((model, sd, S))
            for model, sd, S in runs:
                s = take(S, p["eval_grid"])
                m = L.grid_metrics(p["eval_grid"], y, s)
                rows.append({"task": task, "pseed": ps, "policy": pol, "w_neg": w, "model": model, "seed": sd, **m, **n_aff})
                scores[f"p{ps}__{pol}__w{w:g}__{model}__s{sd}"] = s.astype(np.float32)
        print(task, ps, f"{time.time() - t0:.0f}s", flush=True)
        pd.DataFrame(rows).to_csv(fm + ".partial", sep="\t", index=False)
    np.savez_compressed(f"{OUT}/e1_{task}_scores.npz", **scores)
    pd.DataFrame(rows).to_csv(fm, sep="\t", index=False)
    os.remove(fm + ".partial")

# ----------------------------------------------------------------------------- e2
if STAGE == "e2":
    sel = json.load(open(f"{OUT}/selection.json"))["random"]
    fm = f"{OUT}/e2_metrics.tsv"
    if os.path.exists(fm):
        print("skip e2"); sys.exit(0)
    ctd_set = set(ctd)
    ext = sorted(ev["approved"] - ctd_set - cpd)
    grid = [(c, d) for c in comps for d in dises if (c, d) not in ctd_set and (c, d) not in cpd]
    ext_set = set(ext)
    later_fail = sorted(ev["later_scientific"] - ev["wb_all"] - ev["approved"] - ctd_set - cpd)
    y = np.array([1.0 if q in ext_set else 0.0 for q in grid])
    strata = {"in_wb_scientific": ext_set & ev["wb_scientific"], "in_wb_other": ext_set & ev["wb_other"],
              "not_in_wb": ext_set - ev["wb_all"]}
    gi = {q: i for i, q in enumerate(grid)}
    rows, keep_scores = [], {"labels": y.astype(np.int8)}
    for pol, w in policy_list():
        acts = L.policy_sets(ev, pol)
        Y, Wt = L.training_matrix(comps, dises, ctd, cpd, acts, w_neg=w)
        runs = [("degree", 0, L.degree_reference(comps, dises, Y, Wt))]
        if pol == L.NAMED["no_writeback"] and w == 10.0:
            runs.append(("disease_degree", 0, L.disease_degree(Y)))
        for model in ["mf", "graph", "hybrid"]:
            cfg = sel[model]
            for sd in SEEDS:
                runs.append((model, sd, score_model(model, {**cfg, "epochs": [cfg["epochs"]]}, Y, Wt, seed=sd)[cfg["epochs"]]))
        for model, sd, S in runs:
            s = take(S, grid)
            m = L.grid_metrics(grid, y, s)
            pct = pd.Series(s).rank(pct=True).values           # percentile within the external grid
            rec = {"policy": pol, "w_neg": w, "model": model, "seed": sd, **m}
            for k, v in strata.items():
                rec[f"median_pct_{k}"] = float(np.median([pct[gi[q]] for q in v])) if v else float("nan")
                rec[f"n_{k}"] = len(v)
            fail_s = take(S, later_fail); ext_s = take(S, ext)
            from sklearn.metrics import roc_auc_score
            rec["E3_AUROC_approved_vs_later_failure"] = float(roc_auc_score(np.r_[np.ones(len(ext_s)), np.zeros(len(fail_s))], np.r_[ext_s, fail_s]))
            rec["E3_median_pct_later_failure"] = float(np.median([pct[gi[q]] for q in later_fail]))
            rec["n_later_failure"] = len(later_fail)
            rows.append(rec)
            keep_scores[f"{pol}__w{w:g}__{model}__s{sd}"] = s.astype(np.float16)
        print(pol, w, flush=True)
    np.savez_compressed(f"{OUT}/e2_scores.npz", **keep_scores)
    json.dump({"grid_n": len(grid), "external_positives": len(ext), "later_failures": len(later_fail),
               "strata": {k: len(v) for k, v in strata.items()}}, open(f"{OUT}/e2_sets.json", "w"), indent=1)
    pd.DataFrame(ext, columns=["compound", "disease"]).to_csv(f"{OUT}/e2_external_positives.tsv", sep="\t", index=False)
    pd.DataFrame(later_fail, columns=["compound", "disease"]).to_csv(f"{OUT}/e3_later_failures.tsv", sep="\t", index=False)
    pd.DataFrame(rows).to_csv(fm, sep="\t", index=False)

# ----------------------------------------------------------------------------- boot
MAIN = ["no_writeback", "flat_negative", "mask_all", "typed", "typed_scoped"]


def sorted_ap_factory(y, s):
    """Pre-sort once; return f(weights) giving weighted non-interpolated AP (exact ties)."""
    o = np.argsort(-s, kind="stable"); yy = y[o]; ss = s[o]
    ends = np.r_[np.flatnonzero(np.diff(ss)), len(ss) - 1]

    def f(w):
        ww = w[o]
        tp = np.cumsum(ww * yy)[ends]; tot = np.cumsum(ww)[ends]
        return float(np.sum(np.diff(np.r_[0.0, tp]) * (tp / np.maximum(tot, 1e-12))) / tp[-1])
    return f


def macro_factory(y, s, dis_codes, n_dis):
    """Per-disease AP, so a disease bootstrap is a weighted mean of fixed per-disease values."""
    vals = np.full(n_dis, np.nan)
    for j in np.unique(dis_codes):
        m = dis_codes == j
        if 0 < y[m].sum() < m.sum():
            vals[j] = L.weighted_ap(y[m], s[m])
    return vals


if STAGE == "boot":
    which = sys.argv[2]           # random | compound | e2
    B = int(os.environ.get("BOOT", 1000))
    fo = f"{OUT}/boot_{which}.json"
    if os.path.exists(fo):
        print("skip boot", which); sys.exit(0)
    rng = np.random.default_rng(20261006)
    dcode = {d: j for j, d in enumerate(dises)}
    if which in TASKS:
        Z = np.load(f"{OUT}/e1_{which}_scores.npz")
        units = []   # (pseed, policy, model, seed, ap_fn, macro_vals, row_disease_codes)
        for ps in PSEEDS:
            p = part(which, ps)
            y = Z[f"p{ps}__labels"].astype(float)
            dc = np.array([dcode[d] for c, d in p["eval_grid"]])
            for pol in MAIN:
                pname = L.NAMED[pol]
                for model, seeds in [("degree", [0]), ("mf", SEEDS), ("graph", SEEDS), ("hybrid", SEEDS)] + ([("disease_degree", [0])] if pol == "no_writeback" else []):
                    for sd in seeds:
                        s = Z[f"p{ps}__{pname}__w10__{model}__s{sd}"].astype(float)
                        units.append((ps, pol, model, sd, sorted_ap_factory(y, s), macro_factory(y, s, dc, len(dises)), dc))
    else:
        Z = np.load(f"{OUT}/e2_scores.npz")
        y = Z["labels"].astype(float)
        ctd_set = set(ctd)
        grid = [(c, d) for c in comps for d in dises if (c, d) not in ctd_set and (c, d) not in cpd]
        dc = np.array([dcode[d] for c, d in grid])
        units = []
        for pol in MAIN:
            pname = L.NAMED[pol]
            for model, seeds in [("degree", [0]), ("mf", SEEDS), ("graph", SEEDS), ("hybrid", SEEDS)] + ([("disease_degree", [0])] if pol == "no_writeback" else []):
                for sd in seeds:
                    s = Z[f"{pname}__w10__{model}__s{sd}"].astype(float)
                    units.append((0, pol, model, sd, sorted_ap_factory(y, s), macro_factory(y, s, dc, len(dises)), dc))

    def stat(wd):
        """wd: weight per disease (bootstrap multiplicity). Returns {(policy, model): (pooledAP, macroAP)}."""
        acc = collections.defaultdict(list)
        for ps, pol, model, sd, f, mv, dcodes in units:
            ok = ~np.isnan(mv)
            macro = float(np.sum(mv[ok] * wd[ok]) / max(np.sum(wd[ok]), 1e-12))
            acc[(pol, model)].append((f(wd[dcodes]), macro))
        return {k: tuple(np.mean(v, axis=0)) for k, v in acc.items()}

    obs = stat(np.ones(len(dises)))
    present = np.unique(np.concatenate([u[6] for u in units]))     # diseases that occur in the evaluation grid
    draws = []
    for b in range(B):
        pick = present[rng.integers(len(present), size=len(present))]
        wd = np.bincount(pick, minlength=len(dises)).astype(float)
        draws.append(stat(wd))
    contrasts = []
    keys = sorted(obs)
    for pol in MAIN[1:]:
        for model in ["degree", "mf", "graph", "hybrid"]:
            contrasts.append((f"{pol} - no_writeback [{model}]", (pol, model), ("no_writeback", model)))
    for a, b in [("graph", "mf"), ("hybrid", "mf"), ("mf", "degree"), ("graph", "degree"), ("hybrid", "degree"), ("degree", "disease_degree")]:
        contrasts.append((f"{a} - {b} [no_writeback]", ("no_writeback", a), ("no_writeback", b)))
    out = {"B": B, "cluster": "disease", "observed": {f"{k[0]}|{k[1]}": {"pooledAP": v[0], "macroAP": v[1]} for k, v in obs.items()},
           "contrasts": {}}
    for name, ka, kb in contrasts:
        rec = {}
        for i, metric in enumerate(["pooledAP", "macroAP"]):
            d_obs = obs[ka][i] - obs[kb][i]
            dd = np.array([dr[ka][i] - dr[kb][i] for dr in draws])
            rec[metric] = {"diff": d_obs, "lo": float(np.percentile(dd, 2.5)), "hi": float(np.percentile(dd, 97.5)),
                           "boot_mean": float(dd.mean())}
        out["contrasts"][name] = rec
    json.dump(out, open(fo, "w"), indent=1)
    for k, v in out["contrasts"].items():
        print(k, {m: (round(x["diff"], 4), round(x["lo"], 4), round(x["hi"], 4)) for m, x in v.items()})
