"""External approved indications among unlabelled stopped pairs of the Hetionet E2 grid, by stop type and drug role
(base: unlabelled pairs without a stop before 2015)."""
import json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit import writeback as L
comps, dises, names, ctd, cpd = L.load_graph(); ev = L.load_evidence()
ctd_set = set(ctd)
grid = {(c, d) for c in comps for d in dises if (c, d) not in ctd_set and (c, d) not in cpd}
ext = (ev["approved"] - ctd_set - cpd) & grid
base = len(ext - ev["wb_all"]) / len(grid - ev["wb_all"])
out = {"grid": len(grid), "ext": len(ext), "rate_other": base}
for k in ["wb_all", "wb_scientific", "wb_scientific_scoped", "wb_scientific_inv", "wb_scientific_inv_scoped"]:
    s = ev[k] & grid
    out[k] = {"pairs": len(s), "ext": len(s & ext), "rate": len(s & ext) / len(s), "fold": len(s & ext) / len(s) / base}
print(json.dumps(out, indent=1))
json.dump(out, open("work/results_role/enrichment_hetionet.json", "w"), indent=1)
