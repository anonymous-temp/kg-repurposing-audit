#!/usr/bin/env python
"""Tables 3-4 and supplementary tables S15-S18 for the role-policy and PrimeKG analyses.
Usage: make_tables_v4.py HET_ROLE_RESULTS PRIMEKG_RESULTS TESTEDNESS_RESULTS OUT_DIR"""
import json, os, sys
import numpy as np
import pandas as pd

HR, PR, TR, OUT = (sys.argv[1:5] if len(sys.argv) > 4 else
                   ("work/results_role", "work/results_primekg", "work/results_testedness", "work/tables"))
os.makedirs(OUT, exist_ok=True)
POLS = ["no_writeback", "flat_negative", "typed", "typed_scoped", "typed_role", "typed_role_scoped", "mask_all"]
PL = {"no_writeback": "No write-back", "flat_negative": "Flat negative", "mask_all": "Mask all", "typed": "Typed", "typed_scoped": "Typed + scope",
      "typed_role": "Typed + role", "typed_role_scoped": "Typed + role + scope"}
NAMED = {"no_writeback": "sci=ignore,other=ignore", "flat_negative": "sci=negate,other=negate", "mask_all": "sci=mask,other=mask",
         "typed": "sci=negate,other=mask", "typed_scoped": "typed_scoped", "typed_role": "typed_role", "typed_role_scoped": "typed_role_scoped"}
ML = {"degree": "Degree reference", "disease_degree": "Disease degree", "mf": "Label-only MF", "graph": "Graph head", "hybrid": "Hybrid"}
purity = json.load(open(f"{HR}/negation_purity.json"))


def md(df):
    cols = list(df.columns)
    return "\n".join(["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"] +
                     ["| " + " | ".join(str(r[c]) for c in cols) + " |" for _, r in df.iterrows()])


def cell(boot, pol, model="graph"):
    base = boot["observed"][f"no_writeback|{model}"]["macroAP"]
    if pol == "no_writeback":
        return f"{base:.3f}"
    c = boot["contrasts"][f"{pol} - no_writeback [{model}]"]["macroAP"]
    return f"{base + c['diff']:.3f} ({c['diff']:+.3f}; {c['lo']:+.3f}, {c['hi']:+.3f})"


def policy_table(res, graph, polkey, models=None):
    models = models or {"random": "graph", "compound": "graph", "e2": "graph"}
    boot = {t: json.load(open(f"{res}/boot_{t}.json")) for t in ["random", "compound", "e2"]}
    e1 = {t: pd.read_csv(f"{res}/e1_{t}_metrics.tsv", sep="\t") for t in ["random", "compound"]}
    e2 = pd.read_csv(f"{res}/e2_metrics.tsv", sep="\t")
    rows = []
    for p in POLS:
        pur = purity[graph].get(p)
        r = {"Policy": PL[p], "Negated pairs (recorded or approved)": "0" if pur is None and p == "no_writeback" else
             ("0" if pur is None else f"{pur['negated_pairs']:,} ({pur['pct']:.0f}%)")}
        for t, tl in [("random", "RE"), ("compound", "CD")]:
            d = e1[t][(e1[t][polkey] == (NAMED[p] if polkey == "policy" and graph == "hetionet" else p)) & (e1[t].model == models[t])]
            r[f"{tl} held-out negated / masked"] = f"{d.test_pos_negated.mean():.1f} / {d.test_pos_masked.mean():.1f}"
            r[f"{tl} per-disease AP (Δ; 95% CI)"] = cell(boot[t], p, models[t])
        r["E2 per-disease AP (Δ; 95% CI)"] = cell(boot["e2"], p, models["e2"])
        d = e2[(e2[polkey] == (NAMED[p] if graph == "hetionet" else p)) & (e2.model == models["e2"])]
        r["E3 AUROC"] = f"{d.E3_AUROC_approved_vs_later_failure.mean():.3f}"
        rows.append(r)
    return pd.DataFrame(rows)


out = {}
out["table3"] = policy_table(HR, "hetionet", "policy")
HAVE_PKG = os.path.exists(f"{PR}/boot_e2.json")
if HAVE_PKG:
    out["table4"] = policy_table(PR, "primekg", "policy", {"random": "mf", "compound": "graph", "e2": "mf"})
# ---------------- S15: PrimeKG selection
sel = json.load(open(f"{PR}/selection.json")) if HAVE_PKG else {}
out["S15_pkg_selection"] = pd.DataFrame([{"Task": t, "Scorer": ML[m], **{k: v for k, v in c.items() if k != "val_AP"}, "Validation AP": f"{c['val_AP']:.3f}"}
                                         for t in sel for m, c in sel[t].items()]).fillna("")
# ---------------- S16: PrimeKG full metrics
rows = []
for t in (["random", "compound"] if HAVE_PKG else []):
    d = pd.read_csv(f"{PR}/e1_{t}_metrics.tsv", sep="\t")
    for (p, m), g in d.groupby(["policy", "model"]):
        rows.append({"Set": {"random": "E1 random edge", "compound": "E1 compound disjoint"}[t], "Policy": PL[p], "Scorer": ML[m],
                     **{k: f"{g[k].mean():.3f} ({g[k].std(ddof=1):.3f})" for k in ["AP", "macroAP_disease", "MRR", "Hits@10", "AUROC"]},
                     "Held-out negated": f"{g.test_pos_negated.mean():.1f}"})
d = pd.read_csv(f"{PR}/e2_metrics.tsv", sep="\t") if HAVE_PKG else pd.DataFrame()
for _, g in d.iterrows():
    rows.append({"Set": "E2 external approvals", "Policy": PL[g.policy], "Scorer": ML[g.model], **{k: f"{g[k]:.3f}" for k in ["AP", "macroAP_disease", "MRR", "Hits@10", "AUROC"]},
                 "Held-out negated": "", "E3 AUROC": f"{g.E3_AUROC_approved_vs_later_failure:.3f}"})
s16 = pd.DataFrame(rows).fillna("") if rows else pd.DataFrame({"Set": [], "Policy": [], "Scorer": []})
s16 = s16.rename(columns={"AP": "Pooled AP", "macroAP_disease": "Per-disease AP"})
s16["_o"] = s16.Policy.map({PL[p]: i for i, p in enumerate(POLS)}); s16 = s16.sort_values(["Set", "_o", "Scorer"], kind="stable").drop(columns="_o")
out["S16_pkg_metrics"] = s16
# contrasts (all scorers) for both graphs, role policies vs typed and mask
rows = []
for gname, res in [("Hetionet", HR)] + ([("PrimeKG", PR)] if HAVE_PKG else []):
    for t in ["random", "compound", "e2"]:
        B = json.load(open(f"{res}/boot_{t}.json"))
        for k, c in B["contrasts"].items():
            if any(k.startswith(x) for x in ["typed_role", "graph - ", "mf - degree", "graph - degree"]) or (gname == "PrimeKG" and " - no_writeback [" in k):
                rows.append({"Graph": gname, "Set": t, "Contrast": k, "Pooled AP Δ (95% CI)": f"{c['pooledAP']['diff']:+.3f} ({c['pooledAP']['lo']:+.3f}, {c['pooledAP']['hi']:+.3f})",
                             "Per-disease AP Δ (95% CI)": f"{c['macroAP']['diff']:+.3f} ({c['macroAP']['lo']:+.3f}, {c['macroAP']['hi']:+.3f})"})
out["S18_contrasts"] = pd.DataFrame(rows)
# ---------------- S17: tested vs approved
rows = []
lab = {"A_vs_N": "Approved vs untested", "F_vs_N": "Later failure vs untested", "T_vs_N": "Newly tested vs untested",
       "A_vs_T": "Approved vs newly tested", "A_vs_F": "Approved vs later failure"}
for g in ["hetionet", "primekg"]:
    f = f"{TR}/{g}.json"
    if not os.path.exists(f):
        continue
    R = json.load(open(f))
    for key, l in lab.items():
        rows.append({"Graph": {"hetionet": "Hetionet", "primekg": "PrimeKG"}[g], "Comparison": l,
                     **{ML[m]: f"{R['auroc'][f'{m}|{key}']['est']:.2f} ({R['auroc'][f'{m}|{key}']['lo']:.2f}, {R['auroc'][f'{m}|{key}']['hi']:.2f})" for m in ["degree", "mf", "graph"]}})
    rows.append({"Graph": {"hetionet": "Hetionet", "primekg": "PrimeKG"}[g], "Comparison": "n (approved / later failure / newly tested / untested)",
                 "Degree reference": f"{R['n']['A']:,} / {R['n']['F']:,} / {R['n']['T']:,} / {R['n']['N']:,}", "Label-only MF": "", "Graph head": ""})
out["S17_tested"] = pd.DataFrame(rows)
for k, v in out.items():
    v.to_csv(f"{OUT}/{k}.tsv", sep="\t", index=False)
    open(f"{OUT}/{k}.md", "w").write(md(v))
print("\n\n".join(f"### {k}\n" + md(v) for k, v in out.items() if k in ("table3", "table4")))
