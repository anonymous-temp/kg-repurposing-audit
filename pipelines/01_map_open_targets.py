#!/usr/bin/env python
"""Map Open Targets 26.09 clinical reports onto Hetionet v1.0 compound/disease nodes.

Outputs (in OUT):
  compound_map.tsv   Hetionet DrugBank id -> ChEMBL ids (OT drug_molecule cross-references)
  disease_map.tsv    OT disease id -> Hetionet DOID (direct xref, curated fix, or nearest mapped ancestor)
  trial_pairs.tsv    one row per (clinical report, Hetionet compound, Hetionet disease)
Idempotent: rewrites its outputs from the same inputs.
"""
import collections, json, os, sys
import pandas as pd, pyarrow.parquet as pq

OT = sys.argv[1] if len(sys.argv) > 1 else "work/data/ot"
HET = sys.argv[2] if len(sys.argv) > 2 else "work/data/hetionet"
OUT = sys.argv[3] if len(sys.argv) > 3 else "work/mapped"
os.makedirs(OUT, exist_ok=True)

nodes = pd.read_csv(f"{HET}/hetionet-v1.0-nodes.tsv", sep="\t")
nodes["kind"] = nodes.kind.str.strip()
comp = set(nodes[nodes.kind == "Compound"].id.str.replace("Compound::", "", regex=False))
dis = set(nodes[nodes.kind == "Disease"].id.str.replace("Disease::", "", regex=False))

# ---- compounds: DrugBank -> ChEMBL (parent and child salt forms) ----
dm = pq.read_table(f"{OT}/drug_molecule.parquet", columns=["id", "crossReferences", "childChemblIds", "parentId"]).to_pandas()
ch2db = collections.defaultdict(set)
for _, r in dm.iterrows():
    dbs = [i for x in (r.crossReferences if r.crossReferences is not None else []) if x["source"] == "drugbank" for i in x["ids"]]
    for db in dbs:
        if db in comp:
            ch2db[r.id].add(db)
            for ch in (r.childChemblIds if r.childChemblIds is not None else []):
                ch2db[ch].add(db)
# parents inherit children's DrugBank ids and vice versa
for _, r in dm.iterrows():
    if r.parentId is not None and r.parentId in ch2db and r.id not in ch2db:
        ch2db[r.id] |= ch2db[r.parentId]
rows = [(db, ch) for ch, s in ch2db.items() for db in s]
pd.DataFrame(rows, columns=["compound", "chembl_id"]).sort_values(["compound", "chembl_id"]).to_csv(f"{OUT}/compound_map.tsv", sep="\t", index=False)

# ---- diseases: OT id -> Hetionet DOID ----
dz = pq.read_table(f"{OT}/disease.parquet", columns=["id", "name", "dbXRefs", "obsoleteXRefs", "ancestors"]).to_pandas()
direct = collections.defaultdict(set)
for _, r in dz.iterrows():
    for x in (r.dbXRefs if r.dbXRefs is not None else []):
        if x in dis:
            direct[r.id].add(x)
# Hetionet terms without a current DOID cross-reference in OT 26.09 (resolved by name / obsolete xrefs)
CURATED = {"DOID:2377": ["EFO_0803536"], "DOID:14227": ["HP_0000027"], "DOID:1595": ["MONDO_0002009"],
           "DOID:175": ["EFO_0003967"], "DOID:9917": ["MONDO_0006294"]}
names = dict(zip(dz.id, dz.name))
for doid, ids in CURATED.items():
    for i in ids:
        if i not in names:
            raise SystemExit(f"curated id missing in OT disease index: {i}")
        direct[i].add(doid)
depth = {r.id: len(r.ancestors) if r.ancestors is not None else 0 for _, r in dz.iterrows()}
anc = {r.id: list(r.ancestors) if r.ancestors is not None else [] for _, r in dz.iterrows()}
mrows = []
for i in dz.id:
    if i in direct:
        for d in sorted(direct[i]):
            mrows.append((i, names[i], d, "direct", 0))
        continue
    mapped = [a for a in anc[i] if a in direct]
    if not mapped:
        continue
    best = max(depth[a] for a in mapped)            # most specific mapped ancestor(s)
    for a in mapped:
        if depth[a] == best:
            for d in sorted(direct[a]):
                mrows.append((i, names[i], d, "ancestor:" + a, depth[i] - depth[a]))
dmap = pd.DataFrame(mrows, columns=["ot_disease", "ot_name", "disease", "route", "depth_gap"]).drop_duplicates(["ot_disease", "disease"])
dmap.to_csv(f"{OUT}/disease_map.tsv", sep="\t", index=False)
o2h = dmap.groupby("ot_disease").disease.apply(set).to_dict()

# ---- clinical reports -> Hetionet pairs ----
cols = ["id", "clinicalStage", "origin", "type", "source", "drugs", "diseases", "trialPhase", "trialOverallStatus", "year",
        "trialStopReasonCategories", "trialWhyStopped", "trialStartDate", "url"]
cr = pq.read_table(f"{OT}/clinical_report.parquet", columns=cols).to_pandas()
out = []
for r in cr.itertuples(index=False):
    cs = {db for x in (r.drugs if r.drugs is not None else []) if x["drugId"] in ch2db for db in ch2db[x["drugId"]]}
    if not cs:
        continue
    ds = {(h, x["diseaseId"]) for x in (r.diseases if r.diseases is not None else []) if x["diseaseId"] in o2h for h in o2h[x["diseaseId"]]}
    if not ds:
        continue
    cats = "|".join(sorted(r.trialStopReasonCategories)) if r.trialStopReasonCategories is not None else ""
    for c in cs:
        for h, od in ds:
            out.append((r.id, c, h, od, r.clinicalStage, r.origin, r.type, r.source, r.trialPhase, r.trialOverallStatus,
                        r.year, cats, r.trialWhyStopped, r.trialStartDate, r.url))
tp = pd.DataFrame(out, columns=["report_id", "compound", "disease", "ot_disease", "stage", "origin", "type", "source", "phase",
                                "status", "year", "stop_categories", "why_stopped", "start_date", "url"])
tp = tp.drop_duplicates(["report_id", "compound", "disease"])
tp.to_csv(f"{OUT}/trial_pairs.tsv", sep="\t", index=False)
print(json.dumps({"hetionet_compounds_mapped": len({c for s in ch2db.values() for c in s}),
                  "ot_diseases_mapped": int(dmap.ot_disease.nunique()), "hetionet_diseases_reached": int(dmap.disease.nunique()),
                  "report_pair_rows": len(tp), "reports": int(tp.report_id.nunique()),
                  "pairs": int(tp[["compound", "disease"]].drop_duplicates().shape[0])}, indent=1))
