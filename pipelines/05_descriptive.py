#!/usr/bin/env python
"""Model-free analyses of stopped trials and of the lexical stop-reason rules.

D1  recorded Hetionet indications that have stopped trials, by stop-reason group
D2  later or concurrent approval of stopped drug-disease pairs, by stop-reason group
    (Hetionet-mapped pairs and, for robustness, all Open Targets drug-disease pairs)
D3  the original lexical cue rules against 3,747 expert-labelled stop reasons
D4  agreement between the lexical rules and the Open Targets classifier on Hetionet trials
Writes JSON/TSV summaries to OUT; rerunning overwrites them with identical content.
"""
import collections, json, os, re, sys
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.metrics import cohen_kappa_score

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit import writeback as L

OT = os.environ.get("OTDIR", "work/data/ot")
GOLD = os.environ.get("GOLD", "work/data/gold/data.json")
OUT = os.environ.get("OUT", "work/results")
os.makedirs(OUT, exist_ok=True)

GROUP = {"Negative": "efficacy", "Safety_Sideeffects": "safety", "Business_Administrative": "operational",
         "Insufficient_Enrollment": "operational", "Logistics_Resources": "operational", "Study_Staff_Moved": "operational",
         "Covid19": "operational", "Another_Study": "operational", "Regulatory": "operational", "Study_Design": "design",
         "Invalid_Reason": "uninformative", "No_Context": "uninformative", "Uncategorised": "uninformative"}
ORDER = ["efficacy", "safety", "operational", "design", "uninformative", "no reason given"]


def groups_of(cats):
    if not isinstance(cats, str) or not cats:
        return {"no reason given"}
    return {GROUP.get(c, "uninformative") for c in cats.split("|")}


comps, dises, names, ctd, cpd = L.load_graph()
ctd_set = set(ctd)
ev = L.load_evidence()
st = ev["stopped"].copy()
st["groups"] = st.stop_categories.apply(groups_of)
res = {}

# ---------------------------------------------------------------- D1
d1 = {"recorded_indications": len(ctd_set),
      "with_any_trial": len(ctd_set & set(ev["trials"].pair)),
      "with_any_stopped_trial": len(ctd_set & set(st.pair)),
      "with_stopped_trial_before_2015": len(ctd_set & ev["wb_all"]),
      "with_scientific_stop_before_2015": len(ctd_set & ev["wb_scientific"])}
for g in ORDER:
    d1[f"with_{g.replace(' ', '_')}_stop_any_date"] = len(ctd_set & set(st[st.groups.apply(lambda s: g in s)].pair))
res["D1"] = d1

# ---------------------------------------------------------------- D2 (Hetionet level)
good = ev["approved"] | ctd_set
rows = []
for g in ORDER:
    sub = st[st.groups.apply(lambda s: g in s)]
    for label, s2 in [("all phases", sub), ("phase 1-3 only", sub[~sub.phase.fillna("").str.contains("PHASE4")])]:
        pairs = set(s2.pair)
        rows.append({"level": "Hetionet", "group": g, "phases": label, "stopped_trials": int(s2.report_id.nunique()),
                     "pairs": len(pairs), "pairs_recorded_or_approved": len(pairs & good),
                     "share": len(pairs & good) / max(len(pairs), 1)})
d2h = pd.DataFrame(rows)

# ---------------------------------------------------------------- D2 (all Open Targets pairs)
cr = pq.read_table(f"{OT}/clinical_report.parquet", columns=["id", "clinicalStage", "origin", "drugs", "diseases", "trialPhase",
                                                             "trialOverallStatus", "trialStopReasonCategories"]).to_pandas()
appr = set()
for r in cr[(cr.clinicalStage == "APPROVAL") & (cr.origin != "CLINICAL_TRIAL")].itertuples(index=False):
    for a in (r.drugs if r.drugs is not None else []):
        for b in (r.diseases if r.diseases is not None else []):
            appr.add((a["drugId"], b["diseaseId"]))
sto = cr[(cr.origin == "CLINICAL_TRIAL") & cr.trialOverallStatus.isin(L.STOP)]
grows = collections.defaultdict(lambda: {"trials": set(), "pairs": set(), "trials_p13": set(), "pairs_p13": set()})
for r in sto.itertuples(index=False):
    cats = "|".join(r.trialStopReasonCategories) if r.trialStopReasonCategories is not None else ""
    p4 = "PHASE4" in (r.trialPhase or "")
    for g in groups_of(cats):
        for a in (r.drugs if r.drugs is not None else []):
            for b in (r.diseases if r.diseases is not None else []):
                q = (a["drugId"], b["diseaseId"])
                grows[g]["trials"].add(r.id); grows[g]["pairs"].add(q)
                if not p4:
                    grows[g]["trials_p13"].add(r.id); grows[g]["pairs_p13"].add(q)
rows = []
for g in ORDER:
    for label, tk, pk in [("all phases", "trials", "pairs"), ("phase 1-3 only", "trials_p13", "pairs_p13")]:
        pairs = grows[g][pk]
        rows.append({"level": "Open Targets", "group": g, "phases": label, "stopped_trials": len(grows[g][tk]), "pairs": len(pairs),
                     "pairs_recorded_or_approved": len(pairs & appr), "share": len(pairs & appr) / max(len(pairs), 1)})
d2 = pd.concat([d2h, pd.DataFrame(rows)])
d2.to_csv(f"{OUT}/D2_stopped_pairs_approval.tsv", sep="\t", index=False)
res["D2"] = d2.to_dict(orient="records")

# ---------------------------------------------------------------- D3 lexical rules vs expert labels
PATTERNS = {   # copied verbatim from kg_audit/registry.py v0.2.0
    "efficacy": r"\bfutil\w*\b|lack of (?:efficacy|effectiveness)|insufficient efficacy|inefficacy|lack of benefit|no (?:clinical )?(?:benefit|efficacy)",
    "safety": r"\bsafety\b|\btoxic\w*\b|\badverse event\w*\b|\btolerab\w*\b",
    "operational": r"recruit\w*|enrol\w*|accrual|funding|financial|business|commercial|sponsor decision|strategic|budget|administrative|staffing|investigator|pandemic|covid",
}
GOLDMAP = {"efficacy": {"Negative"}, "safety": {"Safety_Sideeffects"},
           "operational": {"Business_Administrative", "Insufficient_Enrollment", "Logistics_Resources", "Study_Staff_Moved", "Covid19"}}
gold = [json.loads(l) for l in open(GOLD)]
d3 = {"n_texts": len(gold)}
examples = {}
for k, pat in PATTERNS.items():
    pred = np.array([bool(re.search(pat, g["text"] or "", re.I)) for g in gold])
    true = np.array([bool(set(g["label_descriptions"]) & GOLDMAP[k]) for g in gold])
    tp, fp, fn = int((pred & true).sum()), int((pred & ~true).sum()), int((~pred & true).sum())
    prec = tp / max(tp + fp, 1); rec = tp / max(tp + fn, 1)
    d3[k] = {"gold_positive": int(true.sum()), "rule_positive": int(pred.sum()), "TP": tp, "FP": fp, "FN": fn,
             "precision": prec, "recall": rec, "F1": 2 * prec * rec / max(prec + rec, 1e-12)}
    examples[k] = {"false_positive": [g["text"] for g, p, t in zip(gold, pred, true) if p and not t][:8],
                   "false_negative": [g["text"] for g, p, t in zip(gold, pred, true) if t and not p][:8]}
neg_pat = re.compile(r"\b(no|not|without|unrelated to|not related to|none)\b[^.;]{0,40}\bsafety\b", re.I)
safety_fp = [g["text"] for g in gold if re.search(PATTERNS["safety"], g["text"] or "", re.I) and "Safety_Sideeffects" not in g["label_descriptions"]]
d3["safety_false_positives_with_negation"] = sum(bool(neg_pat.search(t)) for t in safety_fp)
d3["safety_false_positives"] = len(safety_fp)
res["D3"] = d3
json.dump(examples, open(f"{OUT}/D3_rule_error_examples.json", "w"), indent=1)

# ---------------------------------------------------------------- D4 rules vs OT classifier on Hetionet trials
rep = st.drop_duplicates("report_id")
rep = rep[rep.why_stopped.notna()]
d4 = {"n_reports_with_reason": int(len(rep))}
for k, pat in PATTERNS.items():
    rule = rep.why_stopped.apply(lambda t: bool(re.search(pat, t, re.I))).values
    clf = rep.stop_categories.fillna("").apply(lambda s: bool(set(s.split("|")) & GOLDMAP[k])).values
    d4[k] = {"rule_positive": int(rule.sum()), "classifier_positive": int(clf.sum()), "both": int((rule & clf).sum()),
             "kappa": float(cohen_kappa_score(rule, clf))}
res["D4"] = d4

# ---------------------------------------------------------------- D5 why scientific stops coexist with indications
wb = st[st.start < L.CUTOFF_WRITEBACK]
drugs_per_trial = ev["trials"].groupby("report_id").compound.nunique()
npair = ev["trials"].groupby("pair").report_id.nunique()
d5 = {}
for name, sub in [("scientific", wb[wb.scientific]), ("other", wb[~wb.scientific])]:
    reps = sub.report_id.unique()
    d5[name] = {"trials": int(len(reps)),
                "share_trials_listing_2plus_mapped_drugs": float((drugs_per_trial.reindex(reps) >= 2).mean()),
                "median_trials_per_pair": float(npair.reindex(list(set(sub.pair))).median()),
                "pairs": int(sub.pair.nunique()),
                "pairs_recorded_indication": int(len(set(sub.pair) & ctd_set))}
sci_ctd = wb[wb.scientific & wb.pair.isin(ctd_set)]
d5["scientific_stops_on_recorded_indications"] = {
    "pairs": int(sci_ctd.pair.nunique()),
    "pairs_whose_trials_all_list_2plus_drugs": int(sum((drugs_per_trial.reindex(g.report_id.unique()) >= 2).all() for _, g in sci_ctd.groupby("pair"))),
    "pairs_same_concept": int(sci_ctd[sci_ctd.gap == 0].pair.nunique())}
res["D5"] = d5

# ---------------------------------------------------------------- D6 stopped pairs as a plausibility signal
ext = ev["approved"] - ctd_set - cpd
grid_n = len(comps) * len(dises) - len(ctd_set) - len(cpd)
d6 = {}
for name, S in [("wb_all", ev["wb_all"]), ("wb_scientific", ev["wb_scientific"]), ("wb_scientific_same_concept", ev["wb_scientific_scoped"]),
                ("wb_other", ev["wb_other"]), ("stopped_any_date", ev["stopped_any_date"])]:
    S2 = S - ctd_set - cpd
    a = len(S2 & ext); rest = grid_n - len(S2); b = len(ext - S2)
    d6[name] = {"unlabelled_pairs": len(S2), "approved": a, "rate": a / len(S2), "rest_pairs": rest, "rest_approved": b,
                "rest_rate": b / rest, "ratio": (a / len(S2)) / (b / rest)}
res["D6"] = d6
json.dump(res, open(f"{OUT}/descriptive.json", "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k not in ("D2",)}, indent=1))
print(d2.to_string(index=False))
