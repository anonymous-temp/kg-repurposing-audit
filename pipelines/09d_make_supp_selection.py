#!/usr/bin/env python
"""Additional file 1, Sections S20 (tests of the selection explanation) and S21 (case review), from the analysis outputs."""
import json, sys
import pandas as pd
SEL = sys.argv[1] if len(sys.argv) > 1 else "work/results_selection"
OUT = sys.argv[2] if len(sys.argv) > 2 else "work/manuscript/supp_S20_S21.md"
R = {g: json.load(open(f"{SEL}/selection_{g}.json")) for g in ("hetionet", "primekg")}
RK = {g: json.load(open(f"{SEL}/rank_by_intensity_{g}.json")) for g in ("hetionet", "primekg")}
f2 = lambda x: f"{x:.2f}"
L = ["## S20. Tests of the selection explanation", "",
     "This section gives the full results of the model-free tests summarised in Table 2 and Figs. 3–5 of the main text. "
     "Analyses used the script pipelines/19_selection_analyses.py of the released code.", "",
     "**Table S20a.** Logistic models of being a recorded treatment or approved indication among pairs with at least one phase 1–3 trial started before 2015. "
     "Odds ratios with 95% disease-cluster bootstrap intervals (1,000 draws).", "",
     "| Model | Term | Hetionet OR (95% CI) | PrimeKG OR (95% CI) |", "|---|---|---|---|"]
LAB = {"crude": "Crude", "adjusted": "Adjusted for log number of trials", "adjusted_any_stop": "Adjusted, scientific and any stop",
       "adjusted_popularity": "Adjusted, plus trial counts of drug and disease", "adjusted_phase12_intensity": "Intensity from phase 1–2 trials only",
       "adjusted_oncology": "Adjusted, cancers only", "adjusted_non_oncology": "Adjusted, other diseases only",
       "crude_any_stop": "Crude, any stop", "adjusted_any_stop_only": "Adjusted, any stop"}
TERM = {"sci": "Scientific stop", "stop": "Any stop", "ln": "log number of trials", "ln12": "log(1 + phase 1–2 trials)", "lnc": "log trials of the drug", "lnd": "log trials of the disease"}
for k, lab in LAB.items():
    for t in R["hetionet"]["intensity"]["models"][k]["terms"]:
        h = R["hetionet"]["intensity"]["models"][k]["terms"][t]; p = R["primekg"]["intensity"]["models"][k]["terms"][t]
        L.append(f"| {lab} | {TERM[t]} | {f2(h['OR'])} ({f2(h['lo'])}–{f2(h['hi'])}) | {f2(p['OR'])} ({f2(p['lo'])}–{f2(p['hi'])}) |")
n = {g: R[g]["intensity"]["models"]["crude"]["n"] for g in R}
L += ["", f"Tested pairs: {n['hetionet']:,} (Hetionet) and {n['primekg']:,} (PrimeKG); shares of indications {R['hetionet']['intensity']['indication_rate']:.1%} and {R['primekg']['intensity']['indication_rate']:.1%}.", "",
      "**Table S20b.** Recorded treatments or approved indications among tested pairs, by number of phase 1–3 trials started before 2015 and by the presence of a scientific stop (Fig. 4A–B).", "",
      "| Trials | Hetionet, no scientific stop | Hetionet, scientific stop | PrimeKG, no scientific stop | PrimeKG, scientific stop |", "|---|---|---|---|---|"]
for b in ["1", "2", "3-4", "5-9", "10-19", "20+"]:
    cells = []
    for g in ("hetionet", "primekg"):
        for s in (0, 1):
            r = next(x for x in R[g]["intensity"]["by_trial_count"] if x["bin"] == b and x["scientific_stop"] == s)
            cells.append(f"{r['indications']}/{r['pairs']} ({r['share']:.0%})")
    L.append(f"| {b.replace('-', '–').replace('20+', '≥20')} | " + " | ".join(cells) + " |")
L += ["", "**Table S20c.** Recorded treatments or approved indications among written-back pairs (stopped trial started before 2015; palliative or off-label pairs excluded), by number of registered trials of any status started before 2015 (Fig. 4C).", "",
      "| Trials | Hetionet, all stopped | Hetionet, scientific stop | Hetionet, strictest set | PrimeKG, all stopped | PrimeKG, scientific stop | PrimeKG, strictest set |", "|---|---|---|---|---|---|---|"]
for i, b in enumerate(["1", "2", "3-4", "5-9", "10-19", "20+"]):
    cells = []
    for g in ("hetionet", "primekg"):
        for s in ("wb_all", "wb_scientific", "wb_scientific_inv_scoped"):
            r = R[g]["purity"][s]["by_trial_count"][i]
            cells.append(f"{r['indications']}/{r['pairs']}")
    L.append(f"| {b.replace('-', '–').replace('20+', '≥20')} | " + " | ".join(cells) + " |")
L += ["", "Strictest set: same-concept scientific stop of the investigational agent.", "",
      "**Table S20d.** Pairs with a scientific stop in a trial started before 2015, by phase of the stopped trial (a pair is counted under every phase in which it had a scientific stop; Fig. 4E).", "",
      "| Phase | Hetionet | PrimeKG |", "|---|---|---|"]
for i in range(4):
    h = R["hetionet"]["purity"]["scientific_by_phase"][i]; p = R["primekg"]["purity"]["scientific_by_phase"][i]
    L.append(f"| {h['phase']} | {h['indications']}/{h['pairs']} ({h['share']:.0%}; {h['ci'][0]:.0%}–{h['ci'][1]:.0%}) | {p['indications']}/{p['pairs']} ({p['share']:.0%}; {p['ci'][0]:.0%}–{p['ci'][1]:.0%}) |")
L += ["", "**Table S20e.** Degree-preserving permutation null for external approved indications among unlabelled stopped pairs (1,000 permutations; Fig. 3).", "",
      "| Set | Graph | Pairs | Observed | Expected (95% range) | Observed / expected (95% range) |", "|---|---|---|---|---|---|"]
NAMES = {"wb_all": "All stopped pairs", "wb_scientific": "Scientific stop", "wb_scientific_scoped": "… same concept", "wb_scientific_inv": "… investigational agent",
         "wb_scientific_inv_scoped": "… both", "wb_other": "Other stops only"}
for s, lab in NAMES.items():
    for g, gn in (("hetionet", "Hetionet"), ("primekg", "PrimeKG")):
        v = R[g]["null"][s]; hi = "∞" if v["ratio_hi"] > 1e6 else f"{v['ratio_hi']:.1f}"
        L.append(f"| {lab} | {gn} | {v['pairs']:,} | {v['observed']} | {v['null_mean']:.1f} ({v['null_lo']:.0f}–{v['null_hi']:.0f}) | {v['ratio_observed_to_null']:.1f} ({v['ratio_lo']:.1f}–{hi}) |")
L += ["", "**Table S20f.** Scorers built only from registry history on the external grid (E2). Per-disease and pooled AP with 95% disease-cluster bootstrap intervals; E3 AUROC for approved indications versus later scientific failures.", "",
      "| Scorer | Graph | Per-disease AP | Pooled AP | E3 AUROC |", "|---|---|---|---|---|"]
for k, lab in (("trial_count", "Number of trials before 2015"), ("stopped_trial_count", "Number of stopped trials before 2015"), ("any_stop", "Any stopped trial before 2015")):
    for g, gn in (("hetionet", "Hetionet"), ("primekg", "PrimeKG")):
        v = R[g]["history"][k]
        L.append(f"| {lab} | {gn} | {v['macroAP']:.3f} ({v['macroAP_ci'][0]:.3f}–{v['macroAP_ci'][1]:.3f}) | {v['pooledAP']:.3f} ({v['pooledAP_ci'][0]:.3f}–{v['pooledAP_ci'][1]:.3f}) | {v['E3_AUROC']:.3f} |")
L += ["", f"External approved indications with any registered trial before 2015: {R['hetionet']['history']['trial_count']['share_external_with_any_pre2015_trial']:.0%} (Hetionet) and {R['primekg']['history']['trial_count']['share_external_with_any_pre2015_trial']:.0%} (PrimeKG); other unlabelled pairs: 2.1% and 0.4%.", "",
      "**Table S20g.** Median percentile rank on the external grid (scorers fitted without write-back) of unlabelled written-back pairs with a scientific stop, by number of registered trials before 2015, and of external approved indications.", "",
      "| Scorer | Graph | ≤2 trials (n) | ≥10 trials (n) | External approved indications |", "|---|---|---|---|---|"]
for m, lab in (("graph", "Graph head"), ("mf", "Label-only MF"), ("degree", "Degree reference")):
    for g, gn in (("hetionet", "Hetionet"), ("primekg", "PrimeKG")):
        v = RK[g][m]
        L.append(f"| {lab} | {gn} | {v['le2_median_pct']:.0f} ({v['le2_n']}) | {v['ge10_median_pct']:.0f} ({v['ge10_n']}) | {v['approved_median_pct']:.0f} |")
d = R["hetionet"]["disease"]
L += ["", f"**Per-disease heterogeneity (Hetionet, graph head).** Across diseases, the change in per-disease AP under flat negatives correlated with the share of the disease's held-out treatments that the policy negated "
      f"(Spearman ρ = {d['random']['spearman_rho']:.2f} in the random-edge task, {d['random']['diseases']} diseases; ρ = {d['compound']['spearman_rho']:.2f} in the compound-disjoint task, {d['compound']['diseases']} diseases). "
      f"Diseases in which at least half of the held-out treatments were negated lost a median of {-d['random']['median_change_share_ge_half']:.3f} and {-d['compound']['median_change_share_ge_half']:.3f}, against "
      f"{-d['random']['median_change_share_lt_half']:.3f} and {-d['compound']['median_change_share_lt_half']:.3f} for the other diseases.", ""]
# ---------------------------------------------------------------------------------------------- S21
T = pd.read_csv(f"{SEL}/case_taxonomy.tsv", sep="\t").fillna("")
L += ["## S21. Review of indications with the cleanest failure records", "",
      "Cases were recorded treatments or approved indications with a same-concept efficacy or safety stop of the investigational agent in a trial started before 2015: all 40 in Hetionet (H01–H40) and a random 40 of 101 in PrimeKG (P01–P40; random seed 20261007). "
      "For each trial we retrieved the official title, conditions, phase, allocation, enrolment, arm groups and interventions, primary outcome, eligibility and stop reason from the ClinicalTrials.gov API (66 distinct trials). "
      "The dossiers were pre-coded with a large language model (Claude, Anthropic) instructed to use only the dossier text and the codebook below, and every code was checked against the registry record. "
      "Some trials map to both graphs, so cases are not independent.", "",
      "Codebook (one primary code, up to two secondary codes):", "",
      "- A1 Special population: subgroup defined by age, pregnancy, comorbidity or region.",
      "- A2 Stage, line or goal: different disease stage, line of therapy or treatment goal (refractory or relapsed disease, adjuvant or neoadjuvant, maintenance, prophylaxis, acute or perioperative setting).",
      "- B1 Regimen of the drug itself: different dose, schedule, duration, route or formulation.",
      "- B2 Combination or add-on: the drug was part of a new multi-drug regimen or was added to another therapy.",
      "- C Comparative question: comparison with another active treatment or strategy.",
      "- D Endpoint or subtype: a specific complication, outcome or narrow subtype instead of treatment of the disease.",
      "- E Stop not attributable to the drug's efficacy or safety (for example planned interim analysis without futility, accrual, sponsor decision, external evidence, stop for benefit, harm from a comparator).",
      "- F Other or unclear.",
      "- Contradicts: yes only if the trial would plausibly show that the drug does not work for the disease in general, in the population and setting where it is used.", "",
      "**Table S21.** Codes for the 80 reviewed cases.", "",
      "| Case | Drug | Disease | Primary | Secondary | Contradicts | Trials | Justification |", "|---|---|---|---|---|---|---|---|"]
for r in T.itertuples():
    L.append(f"| {r.case_id} | {r.drug} | {r.disease} | {r.primary} | {r.secondary} | {r.contradicts} | {str(r.nct_used).replace(',', ', ')} | {r.justification} |")
open(OUT, "w").write("\n".join(L) + "\n")
print("written", OUT, len(L))
