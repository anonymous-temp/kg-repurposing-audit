#!/usr/bin/env python
"""Tables 4-5 (write-back policies incl. the testing-intensity rule) and Additional file 1, Section S22 (placebo
negatives and testing-intensity policies), from the main and the extra (placebo / intensity) results.
Usage: make_tables_v5.py HET_ROLE HET_EXTRA PKG PKG_EXTRA SELECTION OUT_DIR"""
import json, os, sys
import pandas as pd

HR, HX, PR, PX, SEL, OUT = (sys.argv[1:7] if len(sys.argv) > 6 else ("work/results_role", "work/results_extra", "work/results_primekg",
                                                                      "work/results_extra_primekg", "work/results_selection", "work/tables_v5"))
os.makedirs(OUT, exist_ok=True)
POLS = ["no_writeback", "flat_negative", "typed", "typed_scoped", "typed_role", "typed_role_scoped", "typed_lowint2", "mask_all"]
PL = {"no_writeback": "No write-back", "flat_negative": "Flat negative", "mask_all": "Mask all", "typed": "Typed", "typed_scoped": "Typed + scope",
      "typed_role": "Typed + role", "typed_role_scoped": "Typed + role + scope", "typed_lowint1": "Typed + ≤1 trial", "typed_lowint2": "Typed + ≤2 trials",
      "typed_lowint4": "Typed + ≤4 trials"}
NAMED = {"no_writeback": "sci=ignore,other=ignore", "flat_negative": "sci=negate,other=negate", "mask_all": "sci=mask,other=mask",
         "typed": "sci=negate,other=mask", "typed_scoped": "typed_scoped", "typed_role": "typed_role", "typed_role_scoped": "typed_role_scoped"}
ML = {"degree": "Degree reference", "mf": "Label-only MF", "graph": "Graph head", "hybrid": "Hybrid"}
purity = json.load(open(f"{SEL}/negation_purity.json"))
NEW = lambda p: p.startswith(("typed_lowint", "placebo"))


def md(df):
    cols = list(df.columns)
    return "\n".join(["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"] + ["| " + " | ".join(str(r[c]) for c in cols) + " |" for _, r in df.iterrows()])


def load(res, extra, het):
    boot = {t: json.load(open(f"{res}/boot_{t}.json")) for t in ["random", "compound", "e2"]}
    bx = {t: json.load(open(f"{extra}/boot_{t}.json")) for t in ["random", "compound", "e2"]}
    e1 = {}
    for t in ["random", "compound"]:
        m = pd.read_csv(f"{res}/e1_{t}_metrics.tsv", sep="\t")
        if het:
            if os.path.exists(f"{res}2/e1_{t}_metrics.tsv"):
                m = pd.concat([m, pd.read_csv(f"{res}2/e1_{t}_metrics.tsv", sep="\t")]).drop_duplicates(["pseed", "policy", "model", "seed", "w_neg"])
            m = m[m.w_neg == 10]
        e1[t] = pd.concat([m, pd.read_csv(f"{extra}/e1_{t}_metrics.tsv", sep="\t")], ignore_index=True)
    e2 = pd.read_csv(f"{res}/e2_metrics.tsv", sep="\t")
    if het:
        if os.path.exists(f"{res}2/e2_metrics.tsv"):
            e2 = pd.concat([e2, pd.read_csv(f"{res}2/e2_metrics.tsv", sep="\t")]).drop_duplicates(["policy", "model", "seed", "w_neg"])
        e2 = e2[e2.w_neg == 10]
    e2 = pd.concat([e2, pd.read_csv(f"{extra}/e2_metrics.tsv", sep="\t")], ignore_index=True)
    return boot, bx, e1, e2


def cell(boot, bx, pol, model, ref="no_writeback"):
    B = bx if NEW(pol) else boot
    base = B["observed"][f"{ref}|{model}"]["macroAP"]
    if pol == ref:
        return f"{base:.3f}"
    c = B["contrasts"][f"{pol} - {ref} [{model}]"]["macroAP"]
    return f"{base + c['diff']:.3f} ({c['diff']:+.3f}; {c['lo']:+.3f}, {c['hi']:+.3f})"


def key(pol, het):
    return pol if (NEW(pol) or not het) else NAMED[pol]


def policy_table(res, extra, graph, models):
    het = graph == "hetionet"; boot, bx, e1, e2 = load(res, extra, het); rows = []
    for p in POLS:
        pur = purity[graph].get(p)
        r = {"Policy": PL[p], "Negated pairs (recorded or approved)": "0" if pur is None else f"{pur['negated_pairs']:,} ({pur['pct']:.0f}%)"}
        for t, tl in [("random", "RE"), ("compound", "CD")]:
            d = e1[t][(e1[t].policy == key(p, het)) & (e1[t].model == models[t])]
            r[f"{tl} held-out negated / masked"] = f"{d.test_pos_negated.mean():.1f} / {d.test_pos_masked.mean():.1f}"
            r[f"{tl} per-disease AP (Δ; 95% CI)"] = cell(boot[t], bx[t], p, models[t])
        r["E2 per-disease AP (Δ; 95% CI)"] = cell(boot["e2"], bx["e2"], p, models["e2"])
        d = e2[(e2.policy == key(p, het)) & (e2.model == models["e2"])]
        r["E3 AUROC"] = f"{d.E3_AUROC_approved_vs_later_failure.mean():.3f}"
        rows.append(r)
    return pd.DataFrame(rows)


out = {"table4_hetionet": policy_table(HR, HX, "hetionet", {"random": "graph", "compound": "graph", "e2": "graph"}),
       "table5_primekg": policy_table(PR, PX, "primekg", {"random": "mf", "compound": "graph", "e2": "mf"})}

# ---------------------------------------------------------------- S22a: placebo negatives, contrasts with masking
rows = []
for gname, res, extra, het, sets, models in [
        ("Hetionet", HR, HX, True, ["flat_negative", "typed", "typed_scoped", "typed_role", "typed_role_scoped", "typed_lowint2"], ["graph", "mf"]),
        ("PrimeKG", PR, PX, False, ["flat_negative", "typed", "typed_role_scoped", "typed_lowint2"], ["graph", "mf"])]:
    boot, bx, e1, e2 = load(res, extra, het)
    for t, tl in [("random", "Random edge"), ("compound", "Compound disjoint"), ("e2", "External approvals")]:
        for m in models:
            for s in sets:
                rec = {"Graph": gname, "Outcome": tl, "Scorer": ML[m], "Negated set": PL[s]}
                src = e1[t] if t != "e2" else e2
                rec["Negatives added"] = f"{src[(src.policy == f'placebo_degree__{s}') & (src.model == m)].real_n.mean():.0f}"
                for kind, lab in [("real", "Stopped pairs"), ("degree", "Placebo, same drug and disease counts"), ("uniform", "Placebo, uniform")]:
                    name = s if kind == "real" else f"placebo_{kind}__{s}"
                    c = bx[t]["contrasts"].get(f"{name} - mask_all [{m}]")
                    rec[f"{lab}: Δ vs masking (95% CI)"] = "" if c is None else f"{c['macroAP']['diff']:+.3f} ({c['macroAP']['lo']:+.3f}, {c['macroAP']['hi']:+.3f})"
                rows.append(rec)
out["S22a_placebo"] = pd.DataFrame(rows)

# ---------------------------------------------------------------- S22b: testing-intensity policies, all scorers
rows = []
for gname, res, extra, het, models in [("Hetionet", HR, HX, True, ["degree", "mf", "graph", "hybrid"]), ("PrimeKG", PR, PX, False, ["degree", "mf", "graph"])]:
    boot, bx, e1, e2 = load(res, extra, het)
    for k in ("typed_lowint1", "typed_lowint2", "typed_lowint4"):
        pur = purity[gname.lower()][k]
        for t, tl in [("random", "Random edge"), ("compound", "Compound disjoint"), ("e2", "External approvals")]:
            for m in models:
                cn = bx[t]["contrasts"].get(f"{k} - no_writeback [{m}]"); cm = bx[t]["contrasts"].get(f"{k} - mask_all [{m}]")
                if cn is None:
                    continue
                rows.append({"Graph": gname, "Policy": PL[k], "Negated pairs (recorded or approved)": f"{pur['negated_pairs']:,} ({pur['pct']:.0f}%)", "Outcome": tl, "Scorer": ML[m],
                             "Per-disease AP": f"{bx[t]['observed'][f'{k}|{m}']['macroAP']:.3f}",
                             "Δ vs no write-back (95% CI)": f"{cn['macroAP']['diff']:+.3f} ({cn['macroAP']['lo']:+.3f}, {cn['macroAP']['hi']:+.3f})",
                             "Δ vs masking (95% CI)": f"{cm['macroAP']['diff']:+.3f} ({cm['macroAP']['lo']:+.3f}, {cm['macroAP']['hi']:+.3f})"})
out["S22b_intensity"] = pd.DataFrame(rows)
for k, v in out.items():
    v.to_csv(f"{OUT}/{k}.tsv", sep="\t", index=False)
    open(f"{OUT}/{k}.md", "w").write(md(v))
print("\n\n".join(f"### {k}\n" + md(v) for k, v in out.items() if k.startswith("table")))
