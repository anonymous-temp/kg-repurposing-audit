#!/usr/bin/env python
"""Write v0.4.0 result files in the formats allowed in the public tree (idempotent).

- paper_results/scores/role_policies_{e1_random,e1_compound,e2}_scores.tsv.gz: candidate-level scores of the two
  role-restricted policies (columns <model>__<policy>; MF and hybrid averaged over initialisation seeds), same
  layout as the v0.3.0 score files.
- paper_results/role/trial_roles_hetionet.tsv and paper_results/primekg/trial_roles_primekg.tsv: drug roles.
- paper_results/partitions_primekg/<task>_p<seed>.tsv: PrimeKG positives with their split (fit / val / test).
Usage: 04d_export_release_files.py ROLE_ONLY_RESULTS HETIONET_PARTITIONS MAPPED MAPPED_PRIMEKG PARTITIONS_PRIMEKG
"""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit import writeback as L

RES, PART, MAP, MAPP, PARTP = (sys.argv[1:6] if len(sys.argv) > 5 else
                               ("work/results_role_only", "work/partitions", "work/mapped", "work/mapped_primekg", "work/partitions_primekg"))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper_results")
POLS = ["typed_role", "typed_role_scoped"]
MODELS = ["degree", "mf", "graph", "hybrid"]


def model_cols(Z, prefix):
    cols = {}
    for pol in POLS:
        for m in MODELS:
            keys = [k for k in Z.files if k.startswith(f"{prefix}{pol}__w10__{m}__s")]
            if keys:
                cols[f"{m}__{pol}"] = np.round(np.mean([Z[k].astype(float) for k in keys], 0), 5)
    return cols


for task in ["random", "compound"]:
    Z = np.load(f"{RES}/e1_{task}_scores.npz"); frames = []
    for ps in [42, 1, 2, 3, 4]:
        p = json.load(open(f"{PART}/{task}_p{ps}.json"))
        g = pd.DataFrame(p["eval_grid"], columns=["compound", "disease"])
        g.insert(0, "pseed", ps); g.insert(0, "task", task); g["label"] = Z[f"p{ps}__labels"].astype(int)
        for k, v in model_cols(Z, f"p{ps}__").items():
            g[k] = v
        frames.append(g)
    pd.concat(frames).to_csv(f"{OUT}/scores/role_policies_e1_{task}_scores.tsv.gz", sep="\t", index=False)
comps, dises, names, ctd, cpd = L.load_graph()
ctd_set = set(ctd)
grid = pd.DataFrame([(c, d) for c in comps for d in dises if (c, d) not in ctd_set and (c, d) not in cpd], columns=["compound", "disease"])
Z = np.load(f"{RES}/e2_scores.npz"); grid["external_approved"] = Z["labels"].astype(int)
for k, v in model_cols(Z, "").items():
    grid[k] = v
grid.to_csv(f"{OUT}/scores/role_policies_e2_scores.tsv.gz", sep="\t", index=False)
pd.read_csv(f"{MAP}/trial_roles.tsv", sep="\t").to_csv(f"{OUT}/role/trial_roles_hetionet.tsv", sep="\t", index=False)
pd.read_csv(f"{MAPP}/trial_roles.tsv", sep="\t").to_csv(f"{OUT}/primekg/trial_roles_primekg.tsv", sep="\t", index=False)
G = json.load(open(f"{MAPP}/graph.json")); pc = sorted({c for c, d in G["indication"]}); pdz = G["diseases"]
for f in sorted(os.listdir(PARTP)):
    if f.endswith(".npz"):
        z = np.load(f"{PARTP}/{f}")
        rows = [(pc[c], pdz[d], split) for split in ["fit", "val", "test"] for c, d in zip(z[f"{split}_pos_c"], z[f"{split}_pos_d"])]
        pd.DataFrame(rows, columns=["drug", "disease", "split"]).to_csv(f"{OUT}/partitions_primekg/{f.replace('.npz', '.tsv')}", sep="\t", index=False)
print("exported")
