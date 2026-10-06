#!/usr/bin/env python
"""Markdown tables for the manuscript and supplement from the write-back results."""
import json, os, sys
import numpy as np
import pandas as pd

RES = sys.argv[1] if len(sys.argv) > 1 else "work/results"
OUT = sys.argv[2] if len(sys.argv) > 2 else "work/tables"
os.makedirs(OUT, exist_ok=True)
NAMED = {"no_writeback": "sci=ignore,other=ignore", "flat_negative": "sci=negate,other=negate", "mask_all": "sci=mask,other=mask",
         "typed": "sci=negate,other=mask", "typed_scoped": "typed_scoped"}
POL_LAB = {"no_writeback": "No write-back", "flat_negative": "Flat negative", "mask_all": "Mask all", "typed": "Typed",
           "typed_scoped": "Typed + scope"}
MOD_LAB = {"degree": "Degree reference", "disease_degree": "Disease degree", "mf": "Label-only MF", "graph": "Graph head",
           "hybrid": "Hybrid (MF + graph)"}
e1 = {t: pd.read_csv(f"{RES}/e1_{t}_metrics.tsv", sep="\t") for t in ["random", "compound"]}
e2 = pd.read_csv(f"{RES}/e2_metrics.tsv", sep="\t")
boot = {t: json.load(open(f"{RES}/boot_{t}.json")) for t in ["random", "compound", "e2"]}


def msd(x, nd=3):
    return f"{np.mean(x):.{nd}f} ({np.std(x, ddof=1):.{nd}f})" if len(x) > 1 else f"{np.mean(x):.{nd}f}"


def by_partition(d, col):
    return d.groupby("pseed")[col].mean().values      # mean over initialisation seeds, one value per partition


def md(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


out = {}
# ------------------------------------------------ Table 2: scorers without write-back
rows = []
for m in ["degree", "disease_degree", "mf", "graph", "hybrid"]:
    r = {"Scorer": MOD_LAB[m]}
    for t, tl in [("random", "Random edge"), ("compound", "Compound disjoint")]:
        d = e1[t][(e1[t].policy == NAMED["no_writeback"]) & (e1[t].w_neg == 10) & (e1[t].model == m)]
        r[f"{tl}: pooled AP"] = msd(by_partition(d, "AP"))
        r[f"{tl}: per-disease AP"] = msd(by_partition(d, "macroAP_disease"))
        r[f"{tl}: MRR"] = msd(by_partition(d, "MRR"))
        r[f"{tl}: AUROC"] = msd(by_partition(d, "AUROC"))
    d = e2[(e2.policy == NAMED["no_writeback"]) & (e2.w_neg == 10) & (e2.model == m)]
    r["E2: pooled AP"] = f"{d.AP.mean():.3f}"
    r["E2: per-disease AP"] = f"{d.macroAP_disease.mean():.3f}"
    r["E2: AUROC"] = f"{d.AUROC.mean():.3f}"
    rows.append(r)
t2 = pd.DataFrame(rows)
out["table2"] = t2
# contrasts between scorers
rows = []
for t in ["random", "compound", "e2"]:
    for name in ["graph - mf [no_writeback]", "hybrid - mf [no_writeback]", "mf - degree [no_writeback]", "graph - degree [no_writeback]"]:
        c = boot[t]["contrasts"][name]
        rows.append({"Set": {"random": "E1 random edge", "compound": "E1 compound disjoint", "e2": "E2 external approvals"}[t],
                     "Contrast": name.replace(" [no_writeback]", "").replace("graph", "Graph head").replace("hybrid", "Hybrid")
                     .replace("mf", "MF").replace("degree", "degree reference"),
                     "Pooled AP difference (95% CI)": f"{c['pooledAP']['diff']:+.3f} ({c['pooledAP']['lo']:+.3f}, {c['pooledAP']['hi']:+.3f})",
                     "Per-disease AP difference (95% CI)": f"{c['macroAP']['diff']:+.3f} ({c['macroAP']['lo']:+.3f}, {c['macroAP']['hi']:+.3f})"})
out["table2b"] = pd.DataFrame(rows)

# ------------------------------------------------ Table 3: write-back policies
rows = []
for p in ["no_writeback", "flat_negative", "mask_all", "typed", "typed_scoped"]:
    r = {"Policy": POL_LAB[p]}
    for t, tl in [("random", "Random"), ("compound", "Compound")]:
        d = e1[t][(e1[t].policy == NAMED[p]) & (e1[t].w_neg == 10) & (e1[t].model == "graph")]
        r[f"{tl}: held-out treatments negated / masked"] = f"{d.test_pos_negated.mean():.1f} / {d.test_pos_masked.mean():.1f}"
        base = boot[t]["observed"]["no_writeback|graph"]["macroAP"]
        if p == "no_writeback":
            r[f"{tl}: per-disease AP (Δ, 95% CI)"] = f"{base:.3f}"
        else:
            c = boot[t]["contrasts"][f"{p} - no_writeback [graph]"]["macroAP"]
            r[f"{tl}: per-disease AP (Δ, 95% CI)"] = f"{base + c['diff']:.3f} ({c['diff']:+.3f}; {c['lo']:+.3f}, {c['hi']:+.3f})"
    base = boot["e2"]["observed"]["no_writeback|graph"]["macroAP"]
    if p == "no_writeback":
        r["E2: per-disease AP (Δ, 95% CI)"] = f"{base:.3f}"
    else:
        c = boot["e2"]["contrasts"][f"{p} - no_writeback [graph]"]["macroAP"]
        r["E2: per-disease AP (Δ, 95% CI)"] = f"{base + c['diff']:.3f} ({c['diff']:+.3f}; {c['lo']:+.3f}, {c['hi']:+.3f})"
    d = e2[(e2.policy == NAMED[p]) & (e2.w_neg == 10) & (e2.model == "graph")]
    r["E3: AUROC"] = f"{d.E3_AUROC_approved_vs_later_failure.mean():.3f}"
    rows.append(r)
out["table3"] = pd.DataFrame(rows)

# ------------------------------------------------ Supplementary: all policies, scorers, metrics
for t in ["random", "compound"]:
    d = e1[t]
    g = d.groupby(["policy", "w_neg", "model"]).apply(lambda x: pd.Series({
        "pooled AP": msd(by_partition(x, "AP")), "per-disease AP": msd(by_partition(x, "macroAP_disease")),
        "MRR": msd(by_partition(x, "MRR")), "Hits@10": msd(by_partition(x, "Hits@10")), "AUROC": msd(by_partition(x, "AUROC")),
        "held-out negated": f"{x.test_pos_negated.mean():.1f}", "held-out masked": f"{x.test_pos_masked.mean():.1f}"}), include_groups=False).reset_index()
    out[f"S4_{t}"] = g
g = e2.groupby(["policy", "w_neg", "model"])[["AP", "macroAP_disease", "AUROC", "median_pct_in_wb_scientific", "median_pct_in_wb_other",
                                              "median_pct_not_in_wb", "E3_AUROC_approved_vs_later_failure"]].mean().round(3).reset_index()
out["S4_e2"] = g
# policy contrasts for all scorers with CIs (both metrics)
rows = []
for t in ["random", "compound", "e2"]:
    for k, c in boot[t]["contrasts"].items():
        if "no_writeback [" in k and k.split(" - ")[0] in NAMED:
            rows.append({"Set": t, "Contrast": k, **{f"{m} diff": round(c[m]["diff"], 4) for m in ["pooledAP", "macroAP"]},
                         **{f"{m} 95% CI": f"({c[m]['lo']:+.3f}, {c[m]['hi']:+.3f})" for m in ["pooledAP", "macroAP"]}})
out["S5_contrasts"] = pd.DataFrame(rows)
for k, v in out.items():
    v.to_csv(f"{OUT}/{k}.tsv", sep="\t", index=False)
    with open(f"{OUT}/{k}.md", "w") as fh:
        fh.write(md(v))
print("\n\n".join(f"### {k}\n" + md(v) for k, v in out.items() if k in ("table2", "table2b", "table3")))
