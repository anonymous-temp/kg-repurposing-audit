"""Attribution test, descriptive part: are scientific stops in which the compound was the
investigational agent less often recorded treatments / approved indications?"""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit import writeback as L
comps, dises, names, ctd, cpd = L.load_graph(); ev = L.load_evidence()
roles = pd.read_csv(f"{L.DATA}/trial_roles.tsv", sep="\t")
st = ev["stopped"].merge(roles, on=["report_id", "compound"], how="left")
wb = st[st.start < L.CUTOFF_WRITEBACK]
sci = wb[wb.scientific]
pos = set(ctd) | ev["approved"]
def summarise(df, label):
    g = df.groupby("pair").role.apply(set)
    out = {}
    cls = {"any investigational": g[g.apply(lambda s: "investigational" in s)].index,
           "comparator/backbone only": g[g.apply(lambda s: s <= {"comparator", "backbone", "mixed"})].index,
           "other roles only (single-arm combination, head-to-head, unresolved)": g[g.apply(lambda s: "investigational" not in s and not s <= {"comparator", "backbone", "mixed"})].index}
    for k, idx in cls.items():
        idx = set(idx)
        out[k] = {"pairs": len(idx), "recorded_or_approved": len(idx & pos), "pct": round(100 * len(idx & pos) / max(len(idx), 1), 1),
                  "recorded": len(idx & set(ctd))}
    print(label, json.dumps(out, indent=1))
    return out
res = {"trial_compound_roles_all_stopped": roles.role.value_counts().to_dict(),
       "scientific_stop_trial_compound_roles_wb": sci.role.value_counts().to_dict(),
       "wb_scientific": summarise(sci, "scientific stops before 2015"),
       "wb_other": summarise(wb[~wb.pair.isin(set(sci.pair))], "other stops before 2015")}
# scope x role
sco = sci[sci.gap == 0]
res["wb_scientific_scoped"] = summarise(sco, "same-concept scientific stops")
inv_pairs = set(sci[sci.role == "investigational"].pair)
inv_scoped = set(sco[sco.role == "investigational"].pair)
res["n_pairs_inv"] = len(inv_pairs); res["n_pairs_inv_scoped"] = len(inv_scoped)
res["recorded_with_sci_stop"] = len(set(sci.pair) & set(ctd)); res["recorded_with_inv_sci_stop"] = len(inv_pairs & set(ctd))
res["recorded_with_inv_scoped_sci_stop"] = len(inv_scoped & set(ctd))
print({k: v for k, v in res.items() if k.startswith("n_") or k.startswith("recorded")})
json.dump(res, open("work/results_role/role_descriptive.json", "w"), indent=1)
