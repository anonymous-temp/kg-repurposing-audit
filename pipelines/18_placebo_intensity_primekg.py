#!/usr/bin/env python
"""Placebo negatives and testing-intensity-restricted negation on PrimeKG (stages: e1 <task>, e2, contrasts).

Same design as 17_placebo_intensity_hetionet.py, with the compute reductions of the PrimeKG
replication (three partitions, one initialisation seed, no hybrid scorer) and placebo sets for four
negated-set sizes (flat, typed, typed + role + scope, testing intensity k = 2). Bootstrap draws are
evaluated inline with the same disease-multiplicity matrices as 14_writeback_primekg.py, so contrasts
with the stored results of that script are paired.
"""
import json, os, sys, time
os.environ.setdefault("DATA", "work/mapped_primekg")
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit import writeback as L, fastboot as F, placebo as P

PART = os.environ.get("PART", "work/partitions_primekg"); EMB = os.environ.get("EMB", "work/kge_primekg")
BASE = os.environ.get("BASE", "work/results_primekg"); OUT = os.environ.get("OUT", "work/results_extra_primekg")
os.makedirs(OUT, exist_ok=True)
THREADS = int(os.environ.get("THREADS", 2)); B = 1000
PSEEDS, TASKS, KS = [42, 1, 2], ["random", "compound"], [1, 2, 4]
SIZE_SETS = ["wb_all", "wb_scientific", "wb_scientific_inv_scoped", "wb_lowint2"]
SET_POLICY = {"wb_all": "flat_negative", "wb_scientific": "typed", "wb_scientific_inv_scoped": "typed_role_scoped", "wb_lowint2": "typed_lowint2"}
G = json.load(open(f"{L.DATA}/graph.json"))
ctd = [tuple(p) for p in G["indication"]]; ctd_set = set(ctd)
comps = sorted({c for c, d in ctd}); dises = G["diseases"]
cpd = {tuple(p) for p in G["off_label"]} - ctd_set
cpd = {(c, d) for c, d in cpd if c in set(comps) and d in set(dises)}
ci = {c: i for i, c in enumerate(comps)}; di = {d: j for j, d in enumerate(dises)}; nC, nD = len(comps), len(dises)
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


def fit(model, cfg, Y, Wt):
    if model == "mf":
        return L.fit_scorer(Y, Wt, mf_dim=cfg["dim"], wd=cfg["wd"], lr=0.02, epochs=(cfg["epochs"],), seed=1, threads=THREADS)[cfg["epochs"]]
    return L.fit_scorer(Y, Wt, mf_dim=0, emb=emb(cfg["ckpt"]), wd=cfg["wd"], lr=0.01, epochs=(cfg["epochs"],), seed=1, threads=THREADS)[cfg["epochs"]]


def load_part(task, seed):
    z = np.load(f"{PART}/{task}_p{seed}.npz"); out = {"task": task}
    for k in ("fit_pos", "val_pos", "test_pos"):
        out[k] = [(comps[c], dises[d]) for c, d in zip(z[k + "_c"], z[k + "_d"])]
    out["eval_grid_c"] = z["eval_grid_c"]; out["eval_grid_d"] = z["eval_grid_d"]
    return out


def masked_rows(p):
    return set() if p["task"] != "compound" else {c for c, d in p["test_pos"]} | {c for c, d in p["val_pos"]}


def labels_on(gc, gd, pos):
    m = np.zeros((nC, nD), bool)
    for c, d in pos:
        if c in ci and d in di:
            m[ci[c], di[d]] = True
    return m[gc, gd].astype(float)


def draws(seed):
    rng = np.random.default_rng(seed)
    return np.stack([np.bincount(rng.integers(nD, size=nD), minlength=nD) for _ in range(B)]).astype(np.float32)


def allowed_matrix(fit_pos):
    A = np.ones((nC, nD), bool)
    for c, d in set(fit_pos) | cpd | ev["wb_all"]:
        if c in ci and d in di:
            A[ci[c], di[d]] = False
    return A


def placebo(S_eff, allowed, kind, rng, swaps_per_edge=30):
    """Placebo set of the same size as S_eff (see kg_audit.placebo); returns (pairs, original pairs never replaced)."""
    if kind == "uniform":
        return {(comps[a], dises[b]) for a, b in P.uniform_sample(allowed, len(S_eff), rng)}, 0
    out, left = P.rewire([(ci[c], di[d]) for c, d in S_eff], allowed, rng, swaps_per_edge)
    return {(comps[a], dises[b]) for a, b in out}, left


def configs(fit_pos, seed_base):
    out = []
    for k in KS:
        out.append((f"typed_lowint{k}", {p: ("negate" if p in ev[f"wb_lowint{k}"] else "mask") for p in ev["wb_all"]}, {}))
    A = allowed_matrix(fit_pos); fs = set(fit_pos)
    for si, sname in enumerate(SIZE_SETS):
        S_eff = {p for p in ev[sname] if p not in fs and p[0] in ci and p[1] in di}
        for ki, kind in enumerate(["degree", "uniform"]):
            P, left = placebo(S_eff, A, kind, np.random.default_rng(seed_base * 100 + si * 10 + ki))
            acts = {p: "mask" for p in ev["wb_all"]}; acts.update({p: "negate" for p in P})
            out.append((f"placebo_{kind}__{SET_POLICY[sname]}", acts, {"placebo_n": len(P), "placebo_unreplaced": left, "real_n": len(S_eff)}))
    return out


STAGE = sys.argv[1]
if STAGE == "e1":
    task = sys.argv[2]; sel = json.load(open(f"{BASE}/selection.json"))[task]
    fm = f"{OUT}/e1_{task}_metrics.tsv"
    if os.path.exists(fm):
        sys.exit("skip e1 " + task)
    W = draws(20261006); rows, boot = [], {}
    if os.path.exists(fm + ".partial"):
        rows = pd.read_csv(fm + ".partial", sep="\t").to_dict("records"); boot = dict(np.load(f"{OUT}/e1_{task}_boot.partial.npz"))
    t0 = time.time()
    for ps in PSEEDS:
        if any(r["pseed"] == ps for r in rows):
            continue
        p = load_part(task, ps); gc, gd = p["eval_grid_c"], p["eval_grid_d"]
        y = labels_on(gc, gd, p["test_pos"]); test = set(p["test_pos"]); fs = set(p["fit_pos"])
        for name, acts, info in configs(p["fit_pos"], ps):
            Y, Wt = L.training_matrix(comps, dises, p["fit_pos"], cpd, acts, masked_rows=masked_rows(p))
            n_aff = {"negated": sum(1 for q, a in acts.items() if a == "negate" and q not in fs and q[0] in ci and q[1] in di),
                     "test_pos_negated": sum(1 for q in test if acts.get(q) == "negate"), "test_pos_masked": sum(1 for q in test if acts.get(q) == "mask"), **info}
            runs = [("degree", L.degree_reference(comps, dises, Y, Wt))] + [(m, fit(m, sel[m], Y, Wt)) for m in ("mf", "graph")]
            for model, S in runs:
                s = S[gc, gd].astype(np.float64); vals = F.per_disease_ap(y, s, gd, nD, L.weighted_ap)
                ap = L.weighted_ap(y, s); mac = float(np.nanmean(vals))
                rows.append({"task": task, "pseed": ps, "policy": name, "model": model, "AP": ap, "macroAP_disease": mac, **n_aff})
                boot[f"p{ps}__{name}__{model}"] = np.stack([F.pooled_ap_draws(y, s, gd, W), F.macro_draws(vals, W)]).astype(np.float32)
                boot[f"p{ps}__{name}__{model}__obs"] = np.array([ap, mac], np.float32)
            print(task, ps, name, f"{time.time() - t0:.0f}s", flush=True)
        pd.DataFrame(rows).to_csv(fm + ".partial", sep="\t", index=False); np.savez(f"{OUT}/e1_{task}_boot.partial.npz", **boot)
    np.savez(f"{OUT}/e1_{task}_boot.npz", **boot); pd.DataFrame(rows).to_csv(fm, sep="\t", index=False)
    for f in (fm + ".partial", f"{OUT}/e1_{task}_boot.partial.npz"):
        os.remove(f)

if STAGE == "e2":
    sel = json.load(open(f"{BASE}/selection.json"))["random"]
    fm = f"{OUT}/e2_metrics.tsv"
    if os.path.exists(fm):
        sys.exit("skip e2")
    gc, gd = np.nonzero(np.ones((nC, nD), bool))
    excl = np.zeros((nC, nD), bool)
    for c, d in ctd_set | cpd:
        if c in ci and d in di:
            excl[ci[c], di[d]] = True
    keep = ~excl[gc, gd]; gc, gd = gc[keep].astype(np.int32), gd[keep].astype(np.int32)
    inside = lambda S: {q for q in S if q[0] in ci and q[1] in di}
    ext = inside(ev["approved"]) - ctd_set - cpd
    later_fail = inside(ev["later_scientific"]) - ev["wb_all"] - ev["approved"] - ctd_set - cpd
    y = labels_on(gc, gd, ext); fl = labels_on(gc, gd, later_fail)
    W = draws(20261007); rows, boot = [], {}
    for name, acts, info in configs(ctd, 7):
        Y, Wt = L.training_matrix(comps, dises, ctd, cpd, acts)
        n_neg = sum(1 for q, a in acts.items() if a == "negate" and q not in ctd_set and q[0] in ci and q[1] in di)
        runs = [("degree", L.degree_reference(comps, dises, Y, Wt))] + [(m, fit(m, sel[m], Y, Wt)) for m in ("mf", "graph")]
        for model, S in runs:
            s = S[gc, gd].astype(np.float64); vals = F.per_disease_ap(y, s, gd, nD, L.weighted_ap)
            ap = L.weighted_ap(y, s); mac = float(np.nanmean(vals)); fs_, es = s[fl > 0], s[y > 0]
            e3 = float(roc_auc_score(np.r_[np.ones(len(es)), np.zeros(len(fs_))], np.r_[es, fs_]))
            rows.append({"policy": name, "model": model, "AP": ap, "macroAP_disease": mac, "E3_AUROC_approved_vs_later_failure": e3,
                         "negated": n_neg, "ext_negated": sum(1 for q in ext if acts.get(q) == "negate"), **info})
            boot[f"{name}__{model}"] = np.stack([F.pooled_ap_draws(y, s, gd, W), F.macro_draws(vals, W)]).astype(np.float32)
            boot[f"{name}__{model}__obs"] = np.array([ap, mac], np.float32)
        print(name, flush=True)
    np.savez(f"{OUT}/e2_boot.npz", **boot); pd.DataFrame(rows).to_csv(fm, sep="\t", index=False)

if STAGE == "contrasts":
    BASEPOL = ["no_writeback", "flat_negative", "mask_all", "typed", "typed_scoped", "typed_role", "typed_role_scoped"]
    NEW = [f"typed_lowint{k}" for k in KS] + [f"placebo_{kd}__{SET_POLICY[s]}" for s in SIZE_SETS for kd in ("degree", "uniform")]
    for which in ["random", "compound", "e2"]:
        Zb = np.load(f"{BASE}/e1_{which}_boot.npz" if which != "e2" else f"{BASE}/e2_boot.npz")
        Zn = np.load(f"{OUT}/e1_{which}_boot.npz" if which != "e2" else f"{OUT}/e2_boot.npz")
        pre_ = [f"p{ps}__" for ps in PSEEDS] if which != "e2" else [""]
        src = lambda pol: Zb if pol in BASEPOL else Zn
        get = lambda pol, m: np.mean([src(pol)[f"{q}{pol}__{m}"] for q in pre_], 0)
        obs = lambda pol, m: np.mean([src(pol)[f"{q}{pol}__{m}__obs"] for q in pre_], 0)
        out = {"B": B, "cluster": "disease", "observed": {}, "contrasts": {}}
        for m in ("degree", "mf", "graph"):
            for pol in BASEPOL + NEW:
                o = obs(pol, m); out["observed"][f"{pol}|{m}"] = {"pooledAP": float(o[0]), "macroAP": float(o[1])}
                for b in ("no_writeback", "mask_all"):
                    if pol == b:
                        continue
                    da, db = get(pol, m), get(b, m); ob = obs(b, m); rec = {}
                    for i, metric in enumerate(["pooledAP", "macroAP"]):
                        dd = da[i] - db[i]
                        rec[metric] = {"diff": float(o[i] - ob[i]), "lo": float(np.nanpercentile(dd, 2.5)), "hi": float(np.nanpercentile(dd, 97.5))}
                    out["contrasts"][f"{pol} - {b} [{m}]"] = rec
        json.dump(out, open(f"{OUT}/boot_{which}.json", "w"), indent=1); print(which, "contrasts written")
