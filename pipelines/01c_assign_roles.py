#!/usr/bin/env python
"""Role of each mapped compound in each stopped trial (ClinicalTrials.gov arm groups; kg_audit.roles).

Usage: 01c_assign_roles.py MAPPED_DIR OT_DIR ARMS_DIR (hetionet NODES_TSV | primekg GRAPH_JSON)
Writes MAPPED_DIR/trial_roles.tsv (report_id, compound, role, n_arms). Rerunning rewrites identical output.
"""
import glob, json, os, sys
from collections import defaultdict
import pandas as pd, pyarrow.parquet as pq
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit.roles import name_patterns, role_of

MAPPED, OT, ARMS, KIND, NAMES = sys.argv[1:6]
STOP = {"TERMINATED", "WITHDRAWN", "SUSPENDED"}
tp = pd.read_csv(f"{MAPPED}/trial_pairs.tsv", sep="\t", low_memory=False, usecols=["report_id", "compound", "origin", "status"])
tp = tp[(tp.origin == "CLINICAL_TRIAL") & tp.status.isin(STOP)][["report_id", "compound"]].drop_duplicates()
cmap = pd.read_csv(f"{MAPPED}/compound_map.tsv", sep="\t")
dm = pq.read_table(f"{OT}/drug_molecule.parquet", columns=["id", "name", "synonyms", "tradeNames"]).to_pandas().set_index("id")
names = defaultdict(set)
for c, ch in zip(cmap.compound, cmap.chembl_id):
    if ch in dm.index:
        r = dm.loc[ch]
        for x in [r["name"]] + list(r["synonyms"] if r["synonyms"] is not None else []) + list(r["tradeNames"] if r["tradeNames"] is not None else []):
            names[c].add(x["label"] if isinstance(x, dict) else x)
if KIND == "hetionet":
    nodes = pd.read_csv(NAMES, sep="\t"); nodes = nodes[nodes.kind.str.strip() == "Compound"]
    extra = dict(zip(nodes.id.str.replace("Compound::", "", regex=False), nodes.name))
else:
    g = json.load(open(NAMES)); extra = {c: g["names"][c] for c in g["drugs"]}
for c, n in extra.items():
    names[c].add(n)
pat = {c: name_patterns(ns) for c, ns in names.items()}
studies = {}
for f in sorted(glob.glob(f"{ARMS}/batches/*.json")):
    for s in json.load(open(f)):
        p = s["protocolSection"]
        studies[p["identificationModule"]["nctId"].lower()] = p.get("armsInterventionsModule", {})
rows = []
for rid, c in zip(tp.report_id, tp.compound):
    a = studies.get(rid)
    rows.append((rid, c, *(role_of(pat.get(c, []), a) if a is not None else ("not_retrieved", 0))))
out = pd.DataFrame(rows, columns=["report_id", "compound", "role", "n_arms"])
out.to_csv(f"{MAPPED}/trial_roles.tsv", sep="\t", index=False)
print(out.role.value_counts().to_string())
