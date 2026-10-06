#!/usr/bin/env python
"""Export the top-ranked novel candidates of the full-label graph model as handoff records.

Candidates are Hetionet pairs that are neither recorded treatments (CtD) nor palliative
pairs (CpD). Each record carries its prediction provenance and the Open Targets clinical
reports for that pair. Output: JSONL records, a summary, and JSON-Schema validation counts.
"""
import collections, hashlib, json, os, sys
import numpy as np
import pandas as pd
import jsonschema

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get("REPO", os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
sys.path.insert(0, HERE)
from kg_audit import writeback as L
from kg_audit.bridge import export_ranked

RES = os.environ.get("OUT", "work/results")
TOPN = [int(x) for x in os.environ.get("TOPN", "100,500,1000").split(",")]
MODEL = os.environ.get("MODEL", "graph")
AS_OF = "2026-10-06"

comps, dises, names, ctd, cpd = L.load_graph()
ctd_set = set(ctd)
grid = [(c, d) for c in comps for d in dises if (c, d) not in ctd_set and (c, d) not in cpd]
Z = np.load(f"{RES}/e2_scores.npz")
s = Z[f"{L.NAMED['no_writeback']}__w10__{MODEL}__s1"].astype(float)
order = np.argsort(-s, kind="stable")
src_hash = hashlib.sha256(open(f"{RES}/e2_scores.npz", "rb").read()).hexdigest()

tp = pd.read_csv(f"{L.DATA}/trial_pairs.tsv", sep="\t", low_memory=False)
dm = pd.read_csv(f"{L.DATA}/disease_map.tsv", sep="\t")
gap = dict(zip(zip(dm.ot_disease, dm.disease), dm.depth_gap))
tp["gap"] = [gap.get((o, d), -1) for o, d in zip(tp.ot_disease, tp.disease)]
by_pair = collections.defaultdict(list)
for r in tp.to_dict("records"):
    by_pair[(r["compound"], r["disease"])].append(r)

schema = json.load(open(os.path.join(REPO, "schemas", "handoff.schema.json")))
validator = jsonschema.Draft7Validator(schema)
summaries = {}
for n in TOPN:
    cands = [(grid[i][0], grid[i][1], s[i], k + 1) for k, i in enumerate(order[:n])]
    prov = {"model": f"{MODEL} head on DistMult embeddings (no write-back), seed 1", "graph": "Hetionet v1.0",
            "task": "all recorded indications as training labels", "candidate_universe": "all compound-disease pairs except CtD and CpD",
            "write_back_policy": "no_writeback", "seed": 1, "rank_scope": "global rank in the candidate universe",
            "source_sha256": src_hash, "universe_size": len(grid)}
    path = f"{RES}/handoff_top{n}.jsonl"
    summ = export_ranked(cands, by_pair, ctd_set, path, prov, AS_OF)
    errs = 0
    with open(path) as fh:
        for line in fh:
            errs += sum(1 for _ in validator.iter_errors(json.loads(line)))
    summ["schema_errors"] = errs
    summaries[n] = summ
    print(n, json.dumps(summ))
json.dump(summaries, open(f"{RES}/handoff_summary.json", "w"), indent=1)
