#!/usr/bin/env python
"""Audit of the failed (terminated, withdrawn, suspended) drug-indication pairs of repoDB.

repoDB (Brown and Patel, Sci Data 2017) pairs approved indications from DrugCentral with
"failed" pairs taken from ClinicalTrials.gov trials that were terminated, withdrawn or
suspended, and is distributed as a benchmark of true positives and true negatives. This
script asks how many of its failed pairs are indications and how this depends on the
level at which a failure is attached to a disease, on the stop reason and on testing.

Approval evidence: DrugCentral indications (current release, retrieved through the
DrugCentral API) and Open Targets Platform 26.09 approval-stage reports from non-trial
sources. Levels at which a failed pair (drug, condition) can coincide with an indication:
  exact     the same UMLS concept (DrugCentral only; no ontology mapping needed)
  same      the same Open Targets disease concept
  broader   the condition's concept is a direct child of an indicated concept that has at most
            500 descendants (a stage, line or subtype of the indication)
  narrower  the indication is a descendant of the condition's concept
Testing intensity: Open Targets clinical-trial reports of the drug started before 2017
(repoDB was built from AACT in 2016) on the condition's concept or a narrower concept.

Inputs:  REPODB  repoDB full.csv (figshare 3811674)
         OT      directory with clinical_report, drug_molecule and disease parquet files
         DC      directory with DrugCentral API exports omop_indication.json, id_drugbank.json
Outputs (in OUT): repodb_pairs.tsv (one row per failed pair), repodb_audit.json
Idempotent: rewrites its outputs from the same inputs.
"""
import collections, json, os, re, sys
import numpy as np, pandas as pd, pyarrow.parquet as pq

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from kg_audit.failures import failure_type, SCIENTIFIC  # noqa: E402

REPODB = sys.argv[1] if len(sys.argv) > 1 else "work/data/repodb/full.csv"
OT = sys.argv[2] if len(sys.argv) > 2 else "work/data/ot"
DC = sys.argv[3] if len(sys.argv) > 3 else "work/data/drugcentral"
OUT = sys.argv[4] if len(sys.argv) > 4 else "work/repodb"
os.makedirs(OUT, exist_ok=True)
CUTOFF = pd.Timestamp("2017-01-01")
FAILED = {"Terminated", "Withdrawn", "Suspended"}
STOP = {"TERMINATED", "WITHDRAWN", "SUSPENDED"}
BINS = [(1, 1, "1"), (2, 2, "2"), (3, 4, "3-4"), (5, 9, "5-9"), (10, 19, "10-19"), (20, 10**9, ">=20")]


def wilson(k, n, z=1.96):
    if n == 0:
        return [None, None]
    p, d = k / n, 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(c - h, 4), round(c + h, 4)]


def share(mask_k, mask_n=None):
    k, n = (int(np.sum(mask_k)), int(len(mask_k))) if mask_n is None else (int(np.sum(mask_k & mask_n)), int(np.sum(mask_n)))
    return {"n": n, "k": k, "share": round(k / n, 4) if n else None, "ci95": wilson(k, n)}


def norm(s):
    return re.sub(r"[^a-z0-9]+", " ", str(s).lower()).strip() if isinstance(s, str) else ""


def norm_text(s):
    return re.sub(r"\s+", " ", s).strip().lower() if isinstance(s, str) and s.strip() else ""


rd = pd.read_csv(REPODB)
rd["failed"] = rd.status.isin(FAILED)

# ---- drugs: DrugBank -> ChEMBL (parent and salt forms), as in 01_map_open_targets.py ----
dm = pq.read_table(f"{OT}/drug_molecule.parquet", columns=["id", "crossReferences", "childChemblIds", "parentId"]).to_pandas()
db2ch = collections.defaultdict(set)
for r in dm.itertuples(index=False):
    for x in (r.crossReferences if r.crossReferences is not None else []):
        if x["source"] == "drugbank":
            for db in x["ids"]:
                db2ch[db].add(r.id)
                db2ch[db].update(r.childChemblIds if r.childChemblIds is not None else [])
                if r.parentId is not None:
                    db2ch[db].add(r.parentId)
ch2db = collections.defaultdict(set)
for db, chs in db2ch.items():
    for c in chs:
        ch2db[c].add(db)

# ---- diseases: UMLS CUI (or exact name / synonym) -> Open Targets concepts ----
dz = pq.read_table(f"{OT}/disease.parquet", columns=["id", "name", "exactSynonyms", "dbXRefs", "parents", "descendants"]).to_pandas()
dname = dict(zip(dz.id, dz.name))
cui2ot, name2ot, parents, desc = collections.defaultdict(set), collections.defaultdict(set), {}, {}
for r in dz.itertuples(index=False):
    for x in (r.dbXRefs if r.dbXRefs is not None else []):
        if x.startswith("UMLS:"):
            cui2ot[x[5:]].add(r.id)
    for s in [r.name] + list(r.exactSynonyms if r.exactSynonyms is not None else []):
        name2ot[norm(s)].add(r.id)
    parents[r.id] = set(r.parents) if r.parents is not None else set()
    desc[r.id] = set(r.descendants) if r.descendants is not None else set()


def concepts(cui, name):
    if cui2ot.get(cui):
        return cui2ot[cui], "umls_xref"
    if name2ot.get(norm(name)):
        return name2ot[norm(name)], "name"
    return set(), "unmapped"


GENERAL = 500     # parents with more descendants than this are too general to call the condition a stage or subtype


def up1(ids):
    """Direct parents of the condition's concepts that are specific enough to be an indication of which the
    condition is a stage, line or subtype."""
    p1 = set().union(*[parents.get(i, set()) for i in ids]) if ids else set()
    return {p for p in p1 if len(desc.get(p, ())) <= GENERAL} - set(ids)


# ---- approval evidence in (DrugBank, Open Targets disease) space ----
cols = ["id", "clinicalStage", "origin", "drugs", "diseases", "trialOverallStatus", "trialStopReasonCategories",
        "trialWhyStopped", "trialStartDate"]
cr = pq.read_table(f"{OT}/clinical_report.parquet", columns=cols).to_pandas()
appr_ot = collections.defaultdict(set)
for r in cr[(cr.clinicalStage == "APPROVAL") & (cr.origin != "CLINICAL_TRIAL")].itertuples(index=False):
    for a in (r.drugs if r.drugs is not None else []):
        for db in ch2db.get(a["drugId"], ()):
            for b in (r.diseases if r.diseases is not None else []):
                appr_ot[db].add(b["diseaseId"])
ind = pd.DataFrame(json.load(open(f"{DC}/omop_indication.json")))
ids = pd.DataFrame(json.load(open(f"{DC}/id_drugbank.json")))
s2db = ids.groupby("struct_id").identifier.apply(set).to_dict()
dc_cui, appr_dc = collections.defaultdict(set), collections.defaultdict(set)
for r in ind.itertuples(index=False):
    for db in s2db.get(r.struct_id, ()):
        if isinstance(r.umls_cui, str):
            dc_cui[db].add(r.umls_cui)
            appr_dc[db] |= concepts(r.umls_cui, r.concept_name)[0]
appr = {db: appr_ot.get(db, set()) | appr_dc.get(db, set()) for db in set(appr_ot) | set(appr_dc)}

# ---- trials before the cut-off, and stop-reason types by stop text ----
trials = cr[cr.origin == "CLINICAL_TRIAL"]
early = pd.to_datetime(trials.trialStartDate, errors="coerce") < CUTOFF
drug_trials = collections.defaultdict(lambda: collections.defaultdict(set))   # drugbank -> OT disease -> trial ids
text2type = {}
for r, e in zip(trials.itertuples(index=False), early):
    cats = list(r.trialStopReasonCategories) if r.trialStopReasonCategories is not None else []
    if r.trialOverallStatus in STOP and isinstance(r.trialWhyStopped, str):
        text2type.setdefault(norm_text(r.trialWhyStopped), failure_type(cats))
    if not e:
        continue
    for a in (r.drugs if r.drugs is not None else []):
        for db in ch2db.get(a["drugId"], ()):
            for b in (r.diseases if r.diseases is not None else []):
                drug_trials[db][b["diseaseId"]].add(r.id)

# ---- one row per failed repoDB pair ----
approved_pairs = set(zip(rd[rd.status == "Approved"].drug_id, rd[rd.status == "Approved"].ind_id))
rows = []
for (db, cui), g in rd[rd.failed].groupby(["drug_id", "ind_id"]):
    S, route = concepts(cui, g.ind_name.iloc[0])
    A = appr.get(db, set())
    Aot = appr_ot.get(db, set())
    narrow = set(S).union(*[desc.get(o, set()) for o in S]) if S else set()
    broad = up1(S)
    dt = drug_trials.get(db, {})
    n_trials = len(set().union(*[dt.get(o, set()) for o in narrow])) if S and dt else 0
    types = [text2type.get(norm_text(t), "unmatched") if norm_text(t) else "not_reported" for t in g.DetailedStatus]
    rows.append({
        "drug_id": db, "ind_id": cui, "drug_name": g.drug_name.iloc[0], "ind_name": g.ind_name.iloc[0],
        "approved_in_repodb": (db, cui) in approved_pairs, "failed_records": len(g),
        "statuses": "|".join(sorted(set(g.status))), "phases": "|".join(sorted(set(g.phase.dropna()))),
        "record_failure_types": "|".join(types), "any_scientific": any(t in SCIENTIFIC for t in types),
        "only_operational": all(t == "operational" for t in types),
        "disease_route": route, "drug_in_ot": db in db2ch, "drug_in_drugcentral": db in dc_cui,
        "exact_cui": cui in dc_cui.get(db, set()),
        "same": bool(A & S), "broader": bool(A & broad), "narrower": bool(A & (narrow - set(S))),
        "same_ot_only": bool(Aot & S), "broader_ot_only": bool(Aot & broad),
        "trials_before_2017": n_trials,
        "condition_concepts": "|".join(sorted(dname.get(o, o) for o in S)),
        "matched_indications": "|".join(sorted(dname.get(o, o) for o in (A & (S | broad)))),
    })
fp = pd.DataFrame(rows)
# concept-level edges for the generic checks of 25_audit_negative_source.py: failed (drug, concept) edges with their
# testing counts, and approval evidence (drug, concept) for the same drugs
edges, seen = [], set()
for (db, cui), g in rd[rd.failed].groupby(["drug_id", "ind_id"]):
    if (db, cui) in approved_pairs:
        continue
    S, _ = concepts(cui, g.ind_name.iloc[0])
    dt = drug_trials.get(db, {})
    for o in S:
        if (db, o) not in seen:
            seen.add((db, o))
            n = len(set().union(*[dt.get(x, set()) for x in ({o} | desc.get(o, set()))])) if dt else 0
            edges.append((db, o, n))
pd.DataFrame(edges, columns=["drug", "disease", "n_tests"]).to_csv(f"{OUT}/failed_concept_edges.tsv", sep="\t", index=False)
fdrugs = {e[0] for e in edges}
pd.DataFrame([(db, o) for db in sorted(fdrugs) for o in sorted(x for x in appr.get(db, set()) if x)], columns=["drug", "disease"]).to_csv(
    f"{OUT}/approved_concept_edges.tsv", sep="\t", index=False)
fp["same_or_broader"] = fp.same | fp.broader
fp["any_level"] = fp.same | fp.broader | fp.narrower | fp.exact_cui
fp.to_csv(f"{OUT}/repodb_pairs.tsv", sep="\t", index=False)

# ---- summary ----
res = {"inputs": {"records": int(len(rd)), "approved_pairs": len(approved_pairs), "failed_records": int(rd.failed.sum()),
                  "failed_pairs": int(len(fp)), "failed_pairs_also_approved_in_repodb": int(fp.approved_in_repodb.sum()),
                  "drugcentral_indication_rows": int(len(ind))}}
fp = fp[~fp.approved_in_repodb]
mapped = (fp.disease_route != "unmapped") & (fp.drug_in_ot | fp.drug_in_drugcentral)
res["mapping"] = {"failed_pairs": int(len(fp)), "disease_mapped": int((fp.disease_route != "unmapped").sum()),
                  "by_route": fp.disease_route.value_counts().to_dict(), "mapped_pairs": int(mapped.sum())}
res["records_by_stop_type"] = dict(collections.Counter(t for s in fp.record_failure_types for t in s.split("|")))
res["exact_cui_all_pairs"] = share(fp.exact_cui)
m = fp[mapped]
res["mapped_pairs"] = {lvl: share(m[lvl]) for lvl in ["exact_cui", "same", "broader", "narrower", "same_or_broader", "any_level",
                                                      "same_ot_only", "broader_ot_only"]}
# lookup sensitivity: repoDB approved pairs that the same lookup finds approved (DrugCentral is repoDB's own source,
# so this is checked with Open Targets alone)
chk = []
for (db, cui), g in rd[rd.status == "Approved"].groupby(["drug_id", "ind_id"]):
    S, _ = concepts(cui, g.ind_name.iloc[0])
    if S and db in db2ch:
        chk.append(bool(appr_ot.get(db, set()) & (S | up1(S))))
res["lookup_sensitivity_ot_same_or_parent"] = share(np.array(chk))
res["by_reason"] = {"any_scientific_stop": share(m.same_or_broader.values, m.any_scientific.values),
                    "only_operational_stops": share(m.same_or_broader.values, m.only_operational.values)}
res["by_phase"] = {p: share(m.same_or_broader.values, m.phases.str.contains(p, regex=False).values)
                   for p in ["Phase 1", "Phase 2", "Phase 3"]}
res["by_status"] = {s: share(m.same_or_broader.values, m.statuses.str.contains(s).values)
                    for s in ["Terminated", "Withdrawn", "Suspended"]}
t = m.trials_before_2017
res["by_trials"] = {lab: share(m.same_or_broader.values, ((t >= lo) & (t <= hi)).values) for lo, hi, lab in BINS}
res["by_trials"]["0 (no mapped trial)"] = share(m.same_or_broader.values, (t == 0).values)
res["by_trials_same"] = {lab: share(m.same.values, ((t >= lo) & (t <= hi)).values) for lo, hi, lab in BINS}
res["by_trials_same"]["0 (no mapped trial)"] = share(m.same.values, (t == 0).values)
res["by_reason_same"] = {"any_scientific_stop": share(m.same.values, m.any_scientific.values),
                         "only_operational_stops": share(m.same.values, m.only_operational.values)}
res["by_phase_same"] = {p: share(m.same.values, m.phases.str.contains(p, regex=False).values) for p in ["Phase 1", "Phase 2", "Phase 3"]}
res["by_status_same"] = {s: share(m.same.values, m.statuses.str.contains(s).values) for s in ["Terminated", "Withdrawn", "Suspended"]}
json.dump(res, open(f"{OUT}/repodb_audit.json", "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k not in ("by_trials_same",)}, indent=1))
