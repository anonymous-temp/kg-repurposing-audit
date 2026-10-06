#!/usr/bin/env python
"""Where does write-back harm fall? Held-out treatments stratified by their stopped trials, and
per-disease changes under flat negatives (graph head)."""
import json, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit import writeback as L

RES = os.environ.get("OUT", "work/results")
PART = os.environ.get("PART", "work/partitions")
comps, dises, names, ctd, cpd = L.load_graph()
ev = L.load_evidence()
NAMED = L.NAMED
ROLE = os.environ.get("POLICY_SET") == "role"
POLS = ["no_writeback", "flat_negative", "mask_all", "typed", "typed_scoped"] + (["typed_role", "typed_role_scoped"] if ROLE else [])
rows, drows = [], []
for task in ["random", "compound"]:
    Z = np.load(f"{RES}/e1_{task}_scores.npz")
    for ps in [42, 1, 2, 3, 4]:
        p = json.load(open(f"{PART}/{task}_p{ps}.json"))
        grid = [tuple(x) for x in p["eval_grid"]]
        y = Z[f"p{ps}__labels"].astype(bool)
        comp = np.array([c for c, d in grid]); dis = np.array([d for c, d in grid])
        if ROLE:     # split scientific stops by the role of the compound (ClinicalTrials.gov arm groups)
            stratum = np.array(["scientific stop, investigational" if q in ev["wb_scientific_inv"] else
                                "scientific stop, other role" if q in ev["wb_scientific"] else
                                ("other stop" if q in ev["wb_other"] else "no stop") for q in grid])
        else:
            stratum = np.array(["scientific stop" if q in ev["wb_scientific"] else ("other stop" if q in ev["wb_other"] else "no stop")
                                for q in grid])
        for pol in POLS:
            for model in ["graph", "mf"]:
                seeds = [1] if model == "graph" else [1, 2, 3]
                for sd in seeds:
                    s = Z[f"p{ps}__{NAMED[pol]}__w10__{model}__s{sd}"].astype(float)
                    df = pd.DataFrame({"c": comp, "d": dis, "y": y, "s": s, "st": stratum})
                    # rank of each held-out treatment among the candidates of its compound (other positives removed)
                    for c, g in df.groupby("c"):
                        neg = g.s.values[~g.y.values]
                        for _, r in g[g.y].iterrows():
                            rank = 1 + (neg > r.s).sum() + 0.5 * (neg == r.s).sum()
                            rows.append({"task": task, "pseed": ps, "policy": pol, "model": model, "seed": sd, "stratum": r.st,
                                         "rank": rank, "rr": 1 / rank, "hit10": rank <= 10})
                    if model == "graph":
                        for d, g in df.groupby("d"):
                            if 0 < g.y.sum() < len(g):
                                drows.append({"task": task, "pseed": ps, "policy": pol, "disease": d, "name": names[d],
                                              "ap": L.weighted_ap(g.y.values.astype(float), g.s.values), "n_pos": int(g.y.sum())})
r = pd.DataFrame(rows)
summ = (r.groupby(["task", "model", "stratum", "policy"]).agg(n=("rr", "size"), MRR=("rr", "mean"), median_rank=("rank", "median"),
                                                              hits10=("hit10", "mean")).reset_index())
summ["n"] = summ["n"] / summ.apply(lambda x: 5 * (1 if x.model == "graph" else 3), axis=1)   # per partition and seed
summ.to_csv(f"{RES}/e1_strata.tsv", sep="\t", index=False)
d = pd.DataFrame(drows)
dd = d.pivot_table(index=["task", "disease", "name"], columns="policy", values="ap", aggfunc="mean").reset_index()
npos = d[d.policy == "no_writeback"].groupby(["task", "disease"]).n_pos.mean().rename("mean_pos").reset_index()
dd = dd.merge(npos, on=["task", "disease"])
dd["flat_minus_none"] = dd["flat_negative"] - dd["no_writeback"]
dd["mask_minus_none"] = dd["mask_all"] - dd["no_writeback"]
dd.to_csv(f"{RES}/e1_disease_changes.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(summ[summ.model == "graph"].pivot_table(index=["task", "stratum"], columns="policy", values="MRR").round(3))
print(summ[summ.model == "graph"].pivot_table(index=["task", "stratum"], columns="policy", values="n").round(1))
for t in ["random", "compound"]:
    x = dd[dd.task == t].sort_values("flat_minus_none")
    print(t, "diseases", len(x), "with loss under flat", int((x.flat_minus_none < 0).sum()), "gain", int((x.flat_minus_none > 0).sum()))
    print(x.head(8)[["name", "mean_pos", "no_writeback", "flat_negative", "flat_minus_none"]].round(3).to_string(index=False))
