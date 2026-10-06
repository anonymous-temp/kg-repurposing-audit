#!/usr/bin/env python
"""Model-free analyses on PrimeKG, matching the Hetionet analyses (D1, D2, enrichment, roles)."""
import json, os, sys
os.environ["DATA"] = "work/mapped_primekg"
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit import writeback as L
OUT = "work/results_primekg"; os.makedirs(OUT, exist_ok=True)
G = json.load(open(f"{L.DATA}/graph.json"))
ctd = {tuple(p) for p in G["indication"]}; off = {tuple(p) for p in G["off_label"]} - ctd
comps = sorted({c for c, d in ctd}); dises = G["diseases"]; cset, dset = set(comps), set(dises)
ev = L.load_evidence(); st = ev["stopped"].copy()
GROUP = {"Negative": "efficacy", "Safety_Sideeffects": "safety", "Business_Administrative": "operational", "Insufficient_Enrollment": "operational",
         "Logistics_Resources": "operational", "Study_Staff_Moved": "operational", "Covid19": "operational", "Another_Study": "operational",
         "Regulatory": "operational", "Study_Design": "design", "Invalid_Reason": "uninformative", "No_Context": "uninformative", "Uncategorised": "uninformative"}
ORDER = ["efficacy", "safety", "operational", "design", "uninformative", "no reason given"]
st["groups"] = st.stop_categories.apply(lambda c: {GROUP.get(x, "uninformative") for x in c.split("|")} if isinstance(c, str) and c else {"no reason given"})
res = {}
res["D1"] = {"recorded_indications": len(ctd), "with_any_trial": len(ctd & set(ev["trials"].pair)), "with_any_stopped_trial": len(ctd & set(st.pair)),
             "with_stopped_trial_before_2015": len(ctd & ev["wb_all"]), "with_scientific_stop_before_2015": len(ctd & ev["wb_scientific"]),
             "stopped_trials": int(st.report_id.nunique()), "pairs_with_stopped_trial": int(st.pair.nunique()),
             "wb_trials": int(st[st.start < L.CUTOFF_WRITEBACK].report_id.nunique()), "wb_pairs": len(ev["wb_all"]),
             "wb_scientific": len(ev["wb_scientific"]), "wb_scientific_scoped": len(ev["wb_scientific_scoped"]), "wb_other": len(ev["wb_other"]),
             "wb_scientific_inv": len(ev["wb_scientific_inv"]), "wb_scientific_inv_scoped": len(ev["wb_scientific_inv_scoped"]),
             "approved_pairs": len(ev["approved"]), "approved_and_indication": len(ev["approved"] & ctd)}
for g in ORDER:
    res["D1"][f"with_{g.replace(' ', '_')}_stop_any_date"] = len(ctd & set(st[st.groups.apply(lambda s: g in s)].pair))
good = ev["approved"] | ctd
rows = []
for g in ORDER:
    s2 = st[st.groups.apply(lambda s: g in s) & ~st.phase.fillna("").str.contains("PHASE4")]
    pairs = set(s2.pair)
    rows.append({"level": "PrimeKG", "group": g, "phases": "phase 1-3 only", "pairs": len(pairs), "pairs_recorded_or_approved": len(pairs & good),
                 "share": len(pairs & good) / max(len(pairs), 1)})
pd.DataFrame(rows).to_csv(f"{OUT}/D2_stopped_pairs_approval.tsv", sep="\t", index=False); res["D2"] = rows
# background and enrichment on the evaluation grid (drugs with an indication x diseases with an indication)
ngrid = len(comps) * len(dises)
inside = lambda S: {q for q in S if q[0] in cset and q[1] in dset}
res["background"] = {"grid_pairs": ngrid, "recorded_or_approved": len(inside(good)), "rate": len(inside(good)) / ngrid}
unl = ngrid - len(ctd) - len(inside(off))
ext = inside(ev["approved"]) - ctd - off
wb = inside(ev["wb_all"]) - ctd - off
enr = {"unlabelled_pairs": unl, "external_approved": len(ext), "stopped_unlabelled": len(wb), "stopped_ext": len(wb & ext)}
enr["rate_stopped"] = enr["stopped_ext"] / max(enr["stopped_unlabelled"], 1)
enr["rate_other"] = (len(ext) - enr["stopped_ext"]) / (unl - enr["stopped_unlabelled"])
enr["fold"] = enr["rate_stopped"] / enr["rate_other"]
for k in ["wb_scientific", "wb_scientific_scoped", "wb_scientific_inv", "wb_scientific_inv_scoped"]:
    s = inside(ev[k]) - ctd - off
    enr[k] = {"pairs": len(s), "ext": len(s & ext), "rate": len(s & ext) / max(len(s), 1), "fold": (len(s & ext) / max(len(s), 1)) / enr["rate_other"]}
res["enrichment"] = enr
# roles among write-back scientific stops
roles = pd.read_csv(f"{L.DATA}/trial_roles.tsv", sep="\t")
wbt = st[st.start < L.CUTOFF_WRITEBACK].merge(roles, on=["report_id", "compound"], how="left")
sci = wbt[wbt.scientific]
def by_role(df):
    g = df.groupby("pair").role.apply(set); out = {}
    cls = {"any investigational": g[g.apply(lambda s: "investigational" in s)].index,
           "comparator/backbone only": g[g.apply(lambda s: s <= {"comparator", "backbone", "mixed"})].index,
           "other": g[g.apply(lambda s: "investigational" not in s and not s <= {"comparator", "backbone", "mixed"})].index}
    for k, idx in cls.items():
        idx = set(idx); out[k] = {"pairs": len(idx), "recorded_or_approved": len(idx & good), "pct": round(100 * len(idx & good) / max(len(idx), 1), 1)}
    return out
res["roles_scientific"] = by_role(sci); res["roles_scientific_scoped"] = by_role(sci[sci.gap == 0])
res["trial_compound_roles"] = roles.role.value_counts().to_dict()
json.dump(res, open(f"{OUT}/descriptive.json", "w"), indent=1, default=int)
print(json.dumps({k: res[k] for k in ["D1", "background", "enrichment", "roles_scientific", "roles_scientific_scoped"]}, indent=1, default=int))
print(pd.DataFrame(rows).to_string(index=False))
