#!/usr/bin/env python
"""Map Open Targets 26.09 clinical reports onto PrimeKG (v2, 2022) drug and disease nodes.

Outputs (in OUT): graph.json (drugs, diseases, names, indications, off-label, contraindications),
compound_map.tsv, disease_map.tsv, trial_pairs.tsv -- same layout as the Hetionet mapping.
Disease routes: 'direct' (OT MONDO id is a PrimeKG node or a member of a grouped node),
'xref' (OT term shares a DOID/OMIM/Orphanet/UMLS cross-reference with a directly mapped MONDO
term; depth gap 0), 'ancestor:<id>' (most specific mapped ancestor; depth gap > 0).
"""
import collections, json, os, sys
import numpy as np, pandas as pd, pyarrow.parquet as pq

OT, PK, OUT = "work/data/ot", "work/data/primekg", sys.argv[1] if len(sys.argv) > 1 else "work/mapped_primekg"
os.makedirs(OUT, exist_ok=True)
n = pd.read_csv(f"{PK}/nodes.csv", dtype={"node_id": str})
E = np.load(f"{PK}/edges.npz")
rels = list(E["rels"]); x, y, r = E["x"], E["y"], E["r"]
nid = n.node_id.values; ntype = n.node_type.values
def pairs(rel):
    k = rels.index(rel); m = (r == k) & (ntype[x] == "drug") & (ntype[y] == "disease")
    return sorted({(nid[a], nid[b]) for a, b in zip(x[m], y[m])})
ind, off, contra = pairs("indication"), pairs("off-label use"), pairs("contraindication")
drugs = sorted(n[n.node_type == "drug"].node_id)
dis_all = n[n.node_type == "disease"]
labelled = sorted({d for c, d in ind})
names = dict(zip(n[n.node_type.isin(["drug", "disease"])].node_id, n[n.node_type.isin(["drug", "disease"])].node_name))

# ---- compounds: DrugBank -> ChEMBL ----
dm = pq.read_table(f"{OT}/drug_molecule.parquet", columns=["id", "crossReferences", "childChemblIds", "parentId"]).to_pandas()
dset = set(drugs); ch2db = collections.defaultdict(set)
for _, q in dm.iterrows():
    for db in [i for xx in (q.crossReferences if q.crossReferences is not None else []) if xx["source"] == "drugbank" for i in xx["ids"]]:
        if db in dset:
            ch2db[q.id].add(db)
            for ch in (q.childChemblIds if q.childChemblIds is not None else []):
                ch2db[ch].add(db)
for _, q in dm.iterrows():
    if q.parentId is not None and q.parentId in ch2db and q.id not in ch2db:
        ch2db[q.id] |= ch2db[q.parentId]
pd.DataFrame([(db, ch) for ch, s in ch2db.items() for db in s], columns=["compound", "chembl_id"]).sort_values(["compound", "chembl_id"]).to_csv(f"{OUT}/compound_map.tsv", sep="\t", index=False)

# ---- diseases: OT id -> PrimeKG disease node (only labelled diseases are targets) ----
mondo2node = collections.defaultdict(set)
for i, src in zip(dis_all.node_id, dis_all.node_source):
    for part in i.split("_"):
        mondo2node["MONDO_" + part.zfill(7)].add(i)
lab = set(labelled)
dz = pq.read_table(f"{OT}/disease.parquet", columns=["id", "name", "dbXRefs", "ancestors"]).to_pandas()
onames = dict(zip(dz.id, dz.name))
direct = {i: {v for v in mondo2node.get(i, set()) if v in lab} for i in dz.id}
direct = {k: v for k, v in direct.items() if v}
xr2m = collections.defaultdict(set)
for i, xs in zip(dz.id, dz.dbXRefs):
    if i in direct and xs is not None:
        for v in xs:
            if v.split(":")[0] in ("DOID", "OMIM", "Orphanet", "UMLS"):
                xr2m[v].add(i)
viaxref = {}
for i, xs in zip(dz.id, dz.dbXRefs):
    if i in direct or xs is None:
        continue
    tg = {node for v in xs for m in xr2m.get(v, ()) for node in direct[m]}
    if tg:
        viaxref[i] = tg
depth = {i: (len(a) if a is not None else 0) for i, a in zip(dz.id, dz.ancestors)}
anc = {i: (list(a) if a is not None else []) for i, a in zip(dz.id, dz.ancestors)}
rows = []
for i in dz.id:
    if i in direct:
        rows += [(i, onames[i], d, "direct", 0) for d in sorted(direct[i])]; continue
    if i in viaxref:
        rows += [(i, onames[i], d, "xref", 0) for d in sorted(viaxref[i])]; continue
    mapped = [a for a in anc[i] if a in direct]
    if not mapped:
        continue
    best = max(depth[a] for a in mapped)
    for a in mapped:
        if depth[a] == best:
            rows += [(i, onames[i], d, "ancestor:" + a, depth[i] - depth[a]) for d in sorted(direct[a])]
dmap = pd.DataFrame(rows, columns=["ot_disease", "ot_name", "disease", "route", "depth_gap"]).drop_duplicates(["ot_disease", "disease"])
dmap.to_csv(f"{OUT}/disease_map.tsv", sep="\t", index=False)
o2h = dmap.groupby("ot_disease").disease.apply(set).to_dict()

# ---- clinical reports -> PrimeKG pairs ----
cols = ["id", "clinicalStage", "origin", "type", "source", "drugs", "diseases", "trialPhase", "trialOverallStatus", "year",
        "trialStopReasonCategories", "trialWhyStopped", "trialStartDate", "url"]
cr = pq.read_table(f"{OT}/clinical_report.parquet", columns=cols).to_pandas()
out = []
for q in cr.itertuples(index=False):
    cs = {db for xx in (q.drugs if q.drugs is not None else []) if xx["drugId"] in ch2db for db in ch2db[xx["drugId"]]}
    if not cs:
        continue
    ds = {(h, xx["diseaseId"]) for xx in (q.diseases if q.diseases is not None else []) if xx["diseaseId"] in o2h for h in o2h[xx["diseaseId"]]}
    if not ds:
        continue
    cats = "|".join(sorted(q.trialStopReasonCategories)) if q.trialStopReasonCategories is not None else ""
    for c in cs:
        for h, od in ds:
            out.append((q.id, c, h, od, q.clinicalStage, q.origin, q.type, q.source, q.trialPhase, q.trialOverallStatus,
                        q.year, cats, q.trialWhyStopped, q.trialStartDate, q.url))
tp = pd.DataFrame(out, columns=["report_id", "compound", "disease", "ot_disease", "stage", "origin", "type", "source", "phase",
                                "status", "year", "stop_categories", "why_stopped", "start_date", "url"])
tp = tp.drop_duplicates(["report_id", "compound", "disease"])
# a report mapped to a disease both directly and through an ancestor keeps the smallest depth gap
tp.to_csv(f"{OUT}/trial_pairs.tsv", sep="\t", index=False)
# drugs in the grid: drugs with an indication, plus drugs with any mapped clinical report on a labelled disease
grid_drugs = sorted({c for c, d in ind} | set(tp.compound))
json.dump({"drugs": grid_drugs, "diseases": labelled, "indication": ind, "off_label": off, "contraindication": contra,
           "names": {k: names[k] for k in grid_drugs + labelled}}, open(f"{OUT}/graph.json", "w"))
print(json.dumps({"primekg_drugs": len(drugs), "drugs_mapped_to_chembl": len({c for s in ch2db.values() for c in s}),
                  "labelled_diseases": len(labelled), "indications": len(ind), "off_label": len(off), "contraindications": len(contra),
                  "grid_drugs": len(grid_drugs), "grid_cells": len(grid_drugs) * len(labelled),
                  "ot_diseases_mapped": int(dmap.ot_disease.nunique()), "routes": dmap.route.str.split(":").str[0].value_counts().to_dict(),
                  "report_pair_rows": len(tp), "reports": int(tp.report_id.nunique()),
                  "pairs": int(tp[["compound", "disease"]].drop_duplicates().shape[0])}, indent=1))
