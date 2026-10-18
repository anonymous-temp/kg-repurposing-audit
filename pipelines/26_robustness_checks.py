#!/usr/bin/env python
"""Robustness checks added in release 0.7.0 (Hetionet unless stated otherwise).

Usage: 26_robustness_checks.py TRIAL_PAIRS RESULTS OUT_DIR [CHEMBL_CACHE]
  TRIAL_PAIRS   trial_pairs.tsv written by 01_map_open_targets.py (one row per clinical report, compound and disease)
  RESULTS       paper_results (partitions, compound_map.tsv, disease_map.tsv, e2_external_positives.tsv, scores/e2_scores.tsv.gz,
                role/trial_roles_hetionet.tsv, case_review/, literature/)
  OUT_DIR       output directory (paper_results/robustness)
  CHEMBL_CACHE  ChEMBL first-approval years of the mapped ChEMBL ids; fetched from the ChEMBL API if the file does not exist

Outputs:
  temporal_order_hetionet.json   when stopped phase 1-3 trials (started before 2015) of indications started, relative to the
                                 drug's first approval (ChEMBL) and to a dated European or Japanese approval for the same disease
  temporal_rows_hetionet.tsv     the trial-pair rows behind the dated comparison (source data of FigS_temporal_order)
  mh_or_hetionet.json            crude and log-adjusted odds ratios of a scientific stop and a stratified Mantel-Haenszel estimate
  approval_sources_hetionet.json enrichment of external approved indications among unlabelled stopped pairs by approval source
  case_review_sensitivity.json   contradiction counts of the 80 reviewed cases under the union of all codings
  literature_second_screen.json  independent second screening of a random sample of excluded articles
"""
import json, os, sys, time, urllib.request
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.contingency_tables import StratifiedTable

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from kg_audit.failures import SCIENTIFIC, STOP_STATUSES, failure_type  # noqa: E402
from kg_audit.label_audit import wilson  # noqa: E402

TP, R, OUT = sys.argv[1:4]
CACHE = sys.argv[4] if len(sys.argv) > 4 else "work/data/chembl_first_approval.tsv"
os.makedirs(OUT, exist_ok=True)
CUTOFF = "2015-01-01"
P13 = {"EARLY_PHASE1", "PHASE1", "PHASE1/PHASE2", "PHASE2", "PHASE2/PHASE3", "PHASE3"}
DATED = ["EMA Human Drugs", "PMDA"]                       # approval-stage reports that carry a year
LABELS = ["DailyMed", "EMA Human Drugs", "EMA", "PMDA", "FDA"]  # regulatory labels and approvals (no ATC or TTD)


def dump(name, obj):
    json.dump(obj, open(f"{OUT}/{name}", "w"), indent=1)
    print(name, "written")


def share(k, n):
    lo, hi = wilson(k, n) if n else (None, None)
    return {"k": int(k), "n": int(n), "share": round(k / n, 4) if n else None, "ci": [round(max(lo, 0.0), 4), round(hi, 4)] if n else None}


# ---------------------------------------------------------------- inputs
tp = pd.read_csv(TP, sep="\t", low_memory=False)
tp["start"] = pd.to_datetime(tp.start_date, errors="coerce")
tp["sy"] = tp.start.dt.year
dmap = pd.read_csv(f"{R}/disease_map.tsv", sep="\t")
gap = {(a, b): g for a, b, g in zip(dmap.ot_disease, dmap.disease, dmap.depth_gap)}
tp["depth_gap"] = [gap.get((a, b)) for a, b in zip(tp.ot_disease, tp.disease)]
part = json.load(open(f"{R}/partitions/random_p42.json"))
rec = {tuple(x) for k in ("fit_pos", "val_pos", "test_pos") for x in part[k]}
ext = {tuple(x) for x in pd.read_csv(f"{R}/e2_external_positives.tsv", sep="\t").itertuples(index=False)}
grid = set(map(tuple, pd.read_csv(f"{R}/scores/e2_scores.tsv.gz", sep="\t", usecols=["compound", "disease"]).values))
assert len(rec) == 755 and len(ext) == 497 and len(grid) == 211479
roles = pd.read_csv(f"{R}/role/trial_roles_hetionet.tsv", sep="\t")
role = dict(zip(zip(roles.report_id, roles.compound), roles.role))
cm = pd.read_csv(f"{R}/compound_map.tsv", sep="\t")

# ---------------------------------------------------------------- ChEMBL first-approval years
if not os.path.exists(CACHE):
    os.makedirs(os.path.dirname(CACHE) or ".", exist_ok=True)
    ids, rows = sorted(set(cm.chembl_id)), []
    for i in range(0, len(ids), 50):
        url = ("https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_chembl_id__in=" + ",".join(ids[i:i + 50]) +
               "&only=molecule_chembl_id,first_approval,max_phase&limit=100")
        for attempt in range(5):
            try:
                rows += json.load(urllib.request.urlopen(url, timeout=120))["molecules"]
                break
            except Exception:
                time.sleep(5 * (attempt + 1))
        else:
            raise SystemExit(f"ChEMBL request failed for batch {i}")
        time.sleep(0.3)
    pd.DataFrame(rows).to_csv(CACHE + ".tmp", sep="\t", index=False)
    os.replace(CACHE + ".tmp", CACHE)
fa = cm.merge(pd.read_csv(CACHE, sep="\t"), left_on="chembl_id", right_on="molecule_chembl_id", how="left")
first = fa.groupby("compound").first_approval.min()

# ---------------------------------------------------------------- 1. temporal order of stopped trials of indications
st = tp[tp.status.isin(STOP_STATUSES) & (tp.start < CUTOFF) & tp.phase.isin(P13)].copy()
st["ftype"] = st.stop_categories.fillna("").map(lambda s: failure_type([c for c in s.split("|") if c]))
st["sci"] = st.ftype.isin(SCIENTIFIC)
st["inv"] = np.array([role.get((r, c)) == "investigational" for r, c in zip(st.report_id, st.compound)])
st["strict"] = st.sci & (st.depth_gap == 0) & st.inv
st["pair"] = list(zip(st.compound, st.disease))
st["group"] = np.select([st.pair.isin(rec), st.pair.isin(ext), st.pair.isin(grid)], ["recorded", "external", "other"], "palliative")
st["first_approval"] = st.compound.map(first)
ap = tp[(tp.stage == "APPROVAL") & tp.source.isin(DATED) & tp.year.notna()]
dated = ap.groupby(["compound", "disease"]).year.min()
st["dated_approval"] = [dated.get(p) for p in st.pair]


def timing(d):
    out = {"pairs": int(d.pair.nunique()), "trial_pair_rows": len(d), "trials": int(d.report_id.nunique())}
    y = d[d.first_approval.notna()]
    out["after_first_approval"] = share(int((y.sy >= y.first_approval).sum()), len(y))
    out["median_years_after_first_approval"] = float(np.median(y.sy - y.first_approval)) if len(y) else None
    z = d[d.dated_approval.notna()]
    out["dated"] = {"pairs": int(z.pair.nunique()), "trials": int(z.report_id.nunique()),
                    "in_or_after": share(int((z.sy >= z.dated_approval).sum()), len(z)),
                    "strictly_after": share(int((z.sy > z.dated_approval).sum()), len(z))}
    if len(z):
        pp = z.assign(after=z.sy >= z.dated_approval).groupby("pair").after.agg(["all", "any"])
        out["dated"]["pairs_all_stops_in_or_after"] = int(pp["all"].sum())
        out["dated"]["pairs_any_stop_before"] = int((~pp["all"]).sum())
        out["dated"]["median_years_relative_to_approval"] = float(np.median(z.sy - z.dated_approval))
    return out


ind = st[st.group.isin(["recorded", "external"])]
temporal = {"definition": "phase 1-3 terminated, withdrawn or suspended trials started before 2015; one row per trial and pair",
            "first_approval": "ChEMBL molecule.first_approval, minimum over the ChEMBL ids mapped to the compound",
            "dated_approval": "earliest approval-stage report from EMA Human Drugs or PMDA (Open Targets 26.09) for the same pair",
            "compounds_with_first_approval": int(first.notna().sum()), "compounds_mapped": int(first.size)}
for label, d in (("indications", ind), ("recorded", st[st.group == "recorded"]), ("external", st[st.group == "external"]),
                 ("other_unlabelled", st[st.group == "other"])):
    temporal[label] = {"all_stops": timing(d), "scientific": timing(d[d.sci]), "strict": timing(d[d.strict])}
dump("temporal_order_hetionet.json", temporal)
rows = ind[ind.dated_approval.notna()][["report_id", "compound", "disease", "group", "phase", "status", "ftype", "sci", "sy",
                                        "dated_approval", "first_approval"]].rename(columns={"sy": "start_year"})
rows.sort_values(["compound", "disease", "report_id"]).to_csv(f"{OUT}/temporal_rows_hetionet.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- 2. Mantel-Haenszel odds ratio of a scientific stop
t = tp[tp.phase.notna() & ~tp.phase.fillna("").str.contains("PHASE4") & (tp.start < CUTOFF)].copy()
t = t[[(c, d) in grid or (c, d) in rec for c, d in zip(t.compound, t.disease)]]          # palliative pairs excluded
t["sci"] = t.status.isin(STOP_STATUSES) & t.stop_categories.fillna("").map(
    lambda s: failure_type([c for c in s.split("|") if c]) in SCIENTIFIC)
g = t.groupby(["compound", "disease"]).agg(n=("report_id", "nunique"), sci=("sci", "any")).reset_index()
g["ind"] = [(c, d) in rec or (c, d) in ext for c, d in zip(g.compound, g.disease)]
y = g.ind.astype(float)
m0 = sm.Logit(y, sm.add_constant(g[["sci"]].astype(float))).fit(disp=0)
m1 = sm.Logit(y, sm.add_constant(pd.DataFrame({"sci": g.sci.astype(float), "logn": np.log(g.n)}))).fit(disp=0)
g["stratum"] = pd.cut(g.n, [0, 1, 2, 4, 9, 19, 10 ** 6], labels=["1", "2", "3-4", "5-9", "10-19", ">=20"])
tabs, strata = [], []
for s, x in g.groupby("stratum", observed=True):
    a, b = int((x.sci & x.ind).sum()), int((x.sci & ~x.ind).sum())
    c, d = int((~x.sci & x.ind).sum()), int((~x.sci & ~x.ind).sum())
    ha = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))                                 # Haldane-Anscombe, for display
    se = np.sqrt(1 / (a + 0.5) + 1 / (b + 0.5) + 1 / (c + 0.5) + 1 / (d + 0.5))
    strata.append({"stratum": str(s), "stop_indication": a, "stop_other": b, "nostop_indication": c, "nostop_other": d,
                   "or_haldane": ha, "ci": [ha * np.exp(-1.96 * se), ha * np.exp(1.96 * se)]})
    tabs.append(np.array([[a, b], [c, d]]))
mh = StratifiedTable(tabs)
dump("mh_or_hetionet.json", {"pairs": len(g), "indication_share": float(y.mean()), "crude_or": float(np.exp(m0.params.sci)),
                             "adjusted_or_log_trials": float(np.exp(m1.params.sci)), "mh_or": float(mh.oddsratio_pooled),
                             "mh_ci": [float(v) for v in mh.oddsratio_pooled_confint()],
                             "homogeneity_p": float(mh.test_equal_odds().pvalue), "strata": strata})

# ---------------------------------------------------------------- 3. approval sources
apx = tp[(tp.stage == "APPROVAL") & (tp.source != "ClinicalTrials.gov")]
stopped = tp[tp.status.isin(STOP_STATUSES) & (tp.start < CUTOFF)]
wb = {p for p in zip(stopped.compound, stopped.disease) if p in grid}
src = {"write_back_unlabelled_pairs": len(wb), "approval_reports_by_source": apx.source.value_counts().to_dict()}
for label, keep in (("all_non_trial_sources", None), ("regulatory_labels", LABELS), ("us_labels", ["DailyMed", "FDA"])):
    a_ = apx if keep is None else apx[apx.source.isin(keep)]
    E = {p for p in zip(a_.compound, a_.disease) if p in grid}
    if keep is None:
        assert E == ext
    a, n1 = len(wb & E), len(wb)
    c, n0 = len(E) - a, len(grid) - n1
    rr, se = (a / n1) / (c / n0), np.sqrt(1 / a - 1 / n1 + 1 / c - 1 / n0)
    src[label] = {"external_approved": len(E), "among_stopped": a, "share_stopped": a / n1, "share_other": c / n0,
                  "rate_ratio": rr, "ci": [rr * np.exp(-1.96 * se), rr * np.exp(1.96 * se)]}
dump("approval_sources_hetionet.json", src)

# ---------------------------------------------------------------- 4. case review under the union of all codings
C = f"{R}/case_review"
fin = pd.read_csv(f"{C}/case_taxonomy_final_v2.tsv", sep="\t").set_index("case_id")
cod = pd.DataFrame({"dossier_first": pd.read_csv(f"{C}/case_taxonomy_first_coding.tsv", sep="\t").set_index("case_id").contradicts,
                    "dossier_second": pd.read_csv(f"{C}/case_taxonomy_blind_second_coding.tsv", sep="\t").set_index("case_id").contradicts,
                    "source_verified": pd.DataFrame([json.loads(x) for x in open(f"{C}/case_taxonomy_source_verified.jsonl")]
                                                    ).set_index("case_id").contradicts,
                    "final": fin.contradicts})
yes, yes_unc = cod.eq("yes").any(axis=1), cod.isin(["yes", "unclear"]).any(axis=1)
dump("case_review_sensitivity.json", {
    "counts_by_coding": {k: cod[k].value_counts().to_dict() for k in cod},
    "final": {"yes": int((fin.contradicts == "yes").sum()), "trials": int(fin[fin.contradicts == "yes"].nct_used.nunique())},
    "any_coding_yes": {"cases": int(yes.sum()), "trials": int(fin.nct_used[yes].nunique())},
    "any_coding_yes_or_unclear": {"cases": int(yes_unc.sum()), "trials": int(fin.nct_used[yes_unc].nunique()),
                                  "case_ids": list(cod.index[yes_unc])}})

# ---------------------------------------------------------------- 5. literature: second screening of excluded articles
L = f"{R}/literature/second_screen_sample.jsonl"
if os.path.exists(L):
    s2 = [json.loads(x) for x in open(L)]
    first_dec = pd.read_csv(f"{R}/literature/literature_screening.tsv", sep="\t").set_index("pmcid").stopped_as_negatives
    k = sum(r["stopped_as_negatives"] != "no" for r in s2)
    dump("literature_second_screen.json", {"sampled_from_excluded": int((first_dec == "no").sum()), "sample": len(s2),
                                           "first_decision_no": int(sum(first_dec.get(r["pmcid"]) == "no" for r in s2)),
                                           "changed_by_second_screen": share(k, len(s2))})
