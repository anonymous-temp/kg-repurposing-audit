#!/usr/bin/env python
"""Write-back experiment replicated on PrimeKG (stages: select, e1 <task>, e2, contrasts).

Same policies, scorers and metrics as the Hetionet analysis (run_wb.py), with three
reductions for compute: three partition seeds, one initialisation seed and no hybrid model.
Bootstrap draws (disease clusters, B = 1000, shared across policies, models and partitions)
are evaluated inline so that score vectors need not be stored. Each stage writes its own
outputs and is skipped when they exist.
"""
import json, os, sys, time
os.environ.setdefault("DATA", "work/mapped_primekg")
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit import writeback as L, fastboot as F

PART, EMB, OUT = "work/partitions_primekg", "work/kge_primekg", os.environ.get("OUT", "work/results_primekg")
os.makedirs(OUT, exist_ok=True)
THREADS = int(os.environ.get("THREADS", 2)); B = 1000
PSEEDS, TASKS = [42, 1, 2], ["random", "compound"]
POLICIES = ["no_writeback", "flat_negative", "mask_all", "typed", "typed_scoped", "typed_role", "typed_role_scoped"]
G = json.load(open(f"{L.DATA}/graph.json"))
ctd = [tuple(p) for p in G["indication"]]; ctd_set = set(ctd)
comps = sorted({c for c, d in ctd}); dises = G["diseases"]
cpd = {tuple(p) for p in G["off_label"]} - ctd_set
cpd = {(c, d) for c, d in cpd if c in set(comps) and d in set(dises)}   # off-label pairs inside the grid
ci = {c: i for i, c in enumerate(comps)}; di = {d: j for j, d in enumerate(dises)}; nD = len(dises)
ev = L.load_evidence()
_emb = {}


def emb(ep):
    if ep not in _emb:
        _emb[ep] = L.load_embeddings(f"{EMB}/pretrain_distmult_typed_s1_ep{ep}.pt", comps, dises)
    return _emb[ep]


def fit(model, cfg, Y, Wt, epochs):
    if model == "mf":
        return L.fit_scorer(Y, Wt, mf_dim=cfg["dim"], wd=cfg["wd"], lr=0.02, epochs=tuple(epochs), seed=1, threads=THREADS)
    return L.fit_scorer(Y, Wt, mf_dim=0, emb=emb(cfg["ckpt"]), wd=cfg["wd"], lr=0.01, epochs=tuple(epochs), seed=1, threads=THREADS)


def load_part(task, seed):
    """Partition as index arrays (npz written from the json partition files by make_partitions_pkg.py)."""
    z = np.load(f"{PART}/{task}_p{seed}.npz")
    out = {"task": task}
    for k in ("fit_pos", "val_pos", "test_pos"):
        out[k] = [(comps[c], dises[d]) for c, d in zip(z[k + "_c"], z[k + "_d"])]
    for k in ("eval_grid", "val_grid"):
        out[k + "_c"] = z[k + "_c"]; out[k + "_d"] = z[k + "_d"]
    return out


def masked_rows(p):
    return set() if p["task"] != "compound" else {c for c, d in p["test_pos"]} | {c for c, d in p["val_pos"]}


def labels_on(gc, gd, pos):
    m = np.zeros((len(comps), nD), bool)
    for c, d in pos:
        m[ci[c], di[d]] = True
    return m[gc, gd].astype(float)


def rank_metrics(gc, y, s):
    """Filtered reciprocal rank of each positive among the candidates of its compound."""
    o = np.argsort(gc, kind="stable"); cnt = np.bincount(gc, minlength=len(comps)); st = np.r_[0, np.cumsum(cnt)[:-1]]
    rr = np.full(len(y), np.nan)
    for c in np.flatnonzero(cnt):
        idx = o[st[c]:st[c] + cnt[c]]; yy = y[idx] > 0
        if not yy.any():
            continue
        neg = np.sort(s[idx][~yy]); sp = s[idx][yy]
        gt = len(neg) - np.searchsorted(neg, sp, side="right"); eq = np.searchsorted(neg, sp, side="right") - np.searchsorted(neg, sp, side="left")
        rr[idx[yy]] = 1.0 / (1 + gt + 0.5 * eq)
    return rr


def metrics(gc, gd, y, s):
    rr = rank_metrics(gc, y, s); vals = F.per_disease_ap(y, s, gd, nD, L.weighted_ap)
    return {"AP": L.weighted_ap(y, s), "AUROC": float(roc_auc_score(y, s)), "macroAP_disease": float(np.nanmean(vals)),
            "MRR": float(np.nanmean(rr)), "Hits@10": float(np.nanmean(1.0 / rr[~np.isnan(rr)] <= 10)), "n": int(len(y)), "n_pos": int(y.sum())}, rr, vals


def draws(seed):
    rng = np.random.default_rng(seed)
    return np.stack([np.bincount(rng.integers(nD, size=nD), minlength=nD) for _ in range(B)]).astype(np.float32)


STAGE = sys.argv[1]
# ----------------------------------------------------------------------------- select
if STAGE == "select":
    f = f"{OUT}/selection.json"
    if os.path.exists(f):
        sys.exit("skip select")
    EP = [200, 400, 800]; rows = []; sel = {}
    for task in TASKS:
        p = load_part(task, 42)
        Y, Wt = L.training_matrix(comps, dises, p["fit_pos"], cpd, {}, masked_rows=masked_rows(p))
        vy = labels_on(p["val_grid_c"], p["val_grid_d"], p["val_pos"])
        cfgs = [("mf", {"dim": 32, "wd": w}) for w in (1e-5, 1e-4, 1e-3)] + [("graph", {"ckpt": e, "wd": w}) for e in (5, 10) for w in (1e-5, 1e-4, 1e-3)]
        for model, cfg in cfgs:
            t0 = time.time()
            for ep, S in fit(model, cfg, Y, Wt, EP).items():
                rows.append({"task": task, "model": model, **cfg, "epochs": ep, "val_AP": L.weighted_ap(vy, S[p["val_grid_c"], p["val_grid_d"]])})
            print(task, model, cfg, f"{time.time() - t0:.0f}s", rows[-1]["val_AP"], flush=True)
        df = pd.DataFrame(rows); sel[task] = {}
        for model, keys in [("mf", ["dim", "wd", "epochs"]), ("graph", ["ckpt", "wd", "epochs"])]:
            best = df[(df.task == task) & (df.model == model)].sort_values("val_AP", ascending=False).iloc[0]
            sel[task][model] = {k: (int(best[k]) if k in ("dim", "ckpt", "epochs") else float(best[k])) for k in keys} | {"val_AP": float(best["val_AP"])}
    pd.DataFrame(rows).to_csv(f"{OUT}/selection_grid.tsv", sep="\t", index=False); json.dump(sel, open(f, "w"), indent=1); print(sel)

# ----------------------------------------------------------------------------- e1
if STAGE == "e1":
    task = sys.argv[2]; sel = json.load(open(f"{OUT}/selection.json"))[task]
    fm = f"{OUT}/e1_{task}_metrics.tsv"
    if os.path.exists(fm):
        sys.exit("skip e1 " + task)
    W = draws(20261006); rows, strata, boot = [], [], {}
    done = pd.read_csv(fm + ".partial", sep="\t") if os.path.exists(fm + ".partial") else None
    bpart = dict(np.load(f"{OUT}/e1_{task}_boot.partial.npz")) if os.path.exists(f"{OUT}/e1_{task}_boot.partial.npz") else {}
    if done is not None:
        rows = done.to_dict("records"); boot = bpart
        strata = pd.read_csv(f"{OUT}/e1_{task}_strata.partial.tsv", sep="\t").to_dict("records")
    t0 = time.time()
    for ps in PSEEDS:
        if any(r["pseed"] == ps for r in rows):
            continue
        p = load_part(task, ps); gc, gd = p["eval_grid_c"], p["eval_grid_d"]
        y = labels_on(gc, gd, p["test_pos"]); test = set(p["test_pos"]); fitset = set(p["fit_pos"])
        tpos_idx = np.flatnonzero(y)
        strat = np.array(["scientific" if (comps[gc[i]], dises[gd[i]]) in ev["wb_scientific"] else
                          "other" if (comps[gc[i]], dises[gd[i]]) in ev["wb_other"] else "none" for i in tpos_idx])
        for pol in POLICIES:
            acts = L.policy_sets(ev, L.NAMED[pol])
            Y, Wt = L.training_matrix(comps, dises, p["fit_pos"], cpd, acts, masked_rows=masked_rows(p))
            n_aff = {"test_pos_negated": sum(1 for q in test if acts.get(q) == "negate"), "test_pos_masked": sum(1 for q in test if acts.get(q) == "mask"),
                     "negated": sum(1 for q, a in acts.items() if a == "negate" and q not in fitset and q[0] in ci and q[1] in di)}
            runs = [("degree", L.degree_reference(comps, dises, Y, Wt))]
            for model in ["mf", "graph"]:
                runs.append((model, fit(model, sel[model], Y, Wt, [sel[model]["epochs"]])[sel[model]["epochs"]]))
            for model, S in runs:
                s = S[gc, gd].astype(np.float64)
                m, rr, vals = metrics(gc, gd, y, s)
                rows.append({"task": task, "pseed": ps, "policy": pol, "model": model, **m, **n_aff})
                for k in ("scientific", "other", "none"):
                    sel_k = strat == k
                    strata.append({"task": task, "pseed": ps, "policy": pol, "model": model, "stratum": k, "n": int(sel_k.sum()),
                                   "MRR": float(np.mean(rr[tpos_idx][sel_k])) if sel_k.any() else float("nan")})
                boot[f"p{ps}__{pol}__{model}"] = np.stack([F.pooled_ap_draws(y, s, gd, W), F.macro_draws(vals, W)]).astype(np.float32)
                boot[f"p{ps}__{pol}__{model}__obs"] = np.array([m["AP"], m["macroAP_disease"]], np.float32)
            print(task, ps, pol, f"{time.time() - t0:.0f}s", flush=True)
        pd.DataFrame(rows).to_csv(fm + ".partial", sep="\t", index=False)
        pd.DataFrame(strata).to_csv(f"{OUT}/e1_{task}_strata.partial.tsv", sep="\t", index=False)
        np.savez(f"{OUT}/e1_{task}_boot.partial.npz", **boot)
    np.savez(f"{OUT}/e1_{task}_boot.npz", **boot); pd.DataFrame(strata).to_csv(f"{OUT}/e1_{task}_strata.tsv", sep="\t", index=False)
    pd.DataFrame(rows).to_csv(fm, sep="\t", index=False)
    for f in (fm + ".partial", f"{OUT}/e1_{task}_strata.partial.tsv", f"{OUT}/e1_{task}_boot.partial.npz"):
        os.remove(f)

# ----------------------------------------------------------------------------- e2
if STAGE == "e2":
    sel = json.load(open(f"{OUT}/selection.json"))["random"]
    fm = f"{OUT}/e2_metrics.tsv"
    if os.path.exists(fm):
        sys.exit("skip e2")
    gc, gd = np.nonzero(np.ones((len(comps), nD), bool))
    excl = np.zeros((len(comps), nD), bool)
    for c, d in ctd_set | cpd:
        if c in ci and d in di:
            excl[ci[c], di[d]] = True
    keep = ~excl[gc, gd]; gc, gd = gc[keep].astype(np.int32), gd[keep].astype(np.int32)
    inside = lambda S: {q for q in S if q[0] in ci and q[1] in di}
    ext = inside(ev["approved"]) - ctd_set - cpd
    later_fail = inside(ev["later_scientific"]) - ev["wb_all"] - ev["approved"] - ctd_set - cpd
    y = labels_on(gc, gd, ext); fl = labels_on(gc, gd, later_fail)
    W = draws(20261007); rows, boot, keep_scores = [], {}, {"labels": y.astype(np.int8), "later_failure": fl.astype(np.int8), "c": gc, "d": gd}
    for pol in POLICIES:
        acts = L.policy_sets(ev, L.NAMED[pol])
        Y, Wt = L.training_matrix(comps, dises, ctd, cpd, acts)
        runs = [("degree", L.degree_reference(comps, dises, Y, Wt))]
        for model in ["mf", "graph"]:
            runs.append((model, fit(model, sel[model], Y, Wt, [sel[model]["epochs"]])[sel[model]["epochs"]]))
        for model, S in runs:
            s = S[gc, gd].astype(np.float64)
            m, rr, vals = metrics(gc, gd, y, s)
            fs, es = s[fl > 0], s[y > 0]
            m["E3_AUROC_approved_vs_later_failure"] = float(roc_auc_score(np.r_[np.ones(len(es)), np.zeros(len(fs))], np.r_[es, fs]))
            rows.append({"policy": pol, "model": model, **m, "n_later_failure": int(fl.sum())})
            boot[f"{pol}__{model}"] = np.stack([F.pooled_ap_draws(y, s, gd, W), F.macro_draws(vals, W)]).astype(np.float32)
            boot[f"{pol}__{model}__obs"] = np.array([m["AP"], m["macroAP_disease"]], np.float32)
            if pol in ("no_writeback", "mask_all", "flat_negative"):
                keep_scores[f"{pol}__{model}"] = s.astype(np.float16)
        print(pol, flush=True)
    np.savez(f"{OUT}/e2_boot.npz", **boot); np.savez_compressed(f"{OUT}/e2_scores.npz", **keep_scores)
    json.dump({"grid_n": int(len(gc)), "external_positives": int(y.sum()), "later_failures": int(fl.sum())}, open(f"{OUT}/e2_sets.json", "w"), indent=1)
    pd.DataFrame(rows).to_csv(fm, sep="\t", index=False)

# ----------------------------------------------------------------------------- contrasts
if STAGE == "contrasts":
    for which in ["random", "compound", "e2"]:
        fo = f"{OUT}/boot_{which}.json"
        Z = np.load(f"{OUT}/e1_{which}_boot.npz" if which != "e2" else f"{OUT}/e2_boot.npz")
        pre = [f"p{ps}__" for ps in PSEEDS] if which != "e2" else [""]
        get = lambda pol, model: np.mean([Z[f"{q}{pol}__{model}"] for q in pre], 0)            # (2, B): mean over partitions
        obs = lambda pol, model: np.mean([Z[f"{q}{pol}__{model}__obs"] for q in pre], 0)
        out = {"B": B, "cluster": "disease", "observed": {}, "contrasts": {}}
        for pol in POLICIES:
            for model in ["degree", "mf", "graph"]:
                o = obs(pol, model); out["observed"][f"{pol}|{model}"] = {"pooledAP": float(o[0]), "macroAP": float(o[1])}
        pairs = [(f"{pol} - no_writeback [{m}]", (pol, m), ("no_writeback", m)) for pol in POLICIES[1:] for m in ["degree", "mf", "graph"]]
        pairs += [(f"{a} - {b} [{m}]", (a, m), (b, m)) for a, b in [("typed_role", "typed"), ("typed_role", "mask_all"), ("typed_role_scoped", "typed_scoped"),
                                                                    ("typed_role_scoped", "mask_all")] for m in ["degree", "mf", "graph"]]
        pairs += [(f"{a} - {b} [no_writeback]", ("no_writeback", a), ("no_writeback", b)) for a, b in [("graph", "mf"), ("mf", "degree"), ("graph", "degree")]]
        for name, ka, kb in pairs:
            da, db = get(*ka), get(*kb); oa, ob = obs(*ka), obs(*kb); rec = {}
            for i, metric in enumerate(["pooledAP", "macroAP"]):
                dd = da[i] - db[i]
                rec[metric] = {"diff": float(oa[i] - ob[i]), "lo": float(np.nanpercentile(dd, 2.5)), "hi": float(np.nanpercentile(dd, 97.5))}
            out["contrasts"][name] = rec
        json.dump(out, open(fo, "w"), indent=1); print(which, "contrasts written")
