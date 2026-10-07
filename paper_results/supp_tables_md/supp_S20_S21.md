## S20. Tests of the selection explanation

This section gives the full results of the model-free tests summarised in Table 2 and Figs. 3–5 of the main text. Analyses used the script pipelines/19_selection_analyses.py of the released code.

**Table S20a.** Logistic models of being a recorded treatment or approved indication among pairs with at least one phase 1–3 trial started before 2015. Odds ratios with 95% disease-cluster bootstrap intervals (1,000 draws).

| Model | Term | Hetionet OR (95% CI) | PrimeKG OR (95% CI) |
|---|---|---|---|
| Crude | Scientific stop | 2.97 (1.80–4.53) | 1.86 (1.49–2.29) |
| Adjusted for log number of trials | log number of trials | 2.86 (2.52–3.25) | 2.41 (2.21–2.66) |
| Adjusted for log number of trials | Scientific stop | 0.36 (0.24–0.53) | 0.45 (0.34–0.58) |
| Adjusted, scientific and any stop | log number of trials | 3.03 (2.63–3.56) | 2.74 (2.46–3.04) |
| Adjusted, scientific and any stop | Scientific stop | 0.39 (0.25–0.59) | 0.58 (0.44–0.76) |
| Adjusted, scientific and any stop | Any stop | 0.73 (0.47–1.07) | 0.51 (0.43–0.62) |
| Adjusted, plus trial counts of drug and disease | log number of trials | 3.48 (3.00–4.09) | 3.38 (3.07–3.73) |
| Adjusted, plus trial counts of drug and disease | log trials of the drug | 0.51 (0.43–0.62) | 0.49 (0.45–0.54) |
| Adjusted, plus trial counts of drug and disease | log trials of the disease | 0.74 (0.55–0.99) | 0.57 (0.52–0.63) |
| Adjusted, plus trial counts of drug and disease | Scientific stop | 0.44 (0.30–0.62) | 0.55 (0.42–0.71) |
| Intensity from phase 1–2 trials only | log(1 + phase 1–2 trials) | 2.27 (1.92–2.70) | 1.55 (1.38–1.74) |
| Intensity from phase 1–2 trials only | Scientific stop | 0.74 (0.48–1.07) | 1.06 (0.82–1.36) |
| Adjusted, cancers only | log number of trials | 3.15 (2.66–3.77) | 2.43 (2.17–2.73) |
| Adjusted, cancers only | Scientific stop | 0.49 (0.32–0.70) | 0.65 (0.45–0.89) |
| Adjusted, other diseases only | log number of trials | 3.75 (3.08–4.60) | 3.28 (2.87–3.82) |
| Adjusted, other diseases only | Scientific stop | 0.49 (0.26–0.89) | 0.56 (0.36–0.87) |
| Crude, any stop | Any stop | 2.64 (1.97–3.54) | 1.42 (1.24–1.66) |
| Adjusted, any stop | log number of trials | 2.78 (2.41–3.24) | 2.65 (2.38–2.95) |
| Adjusted, any stop | Any stop | 0.64 (0.42–0.93) | 0.47 (0.39–0.58) |

Tested pairs: 4,237 (Hetionet) and 10,328 (PrimeKG); shares of indications 17.9% and 21.1%.

**Table S20b.** Recorded treatments or approved indications among tested pairs, by number of phase 1–3 trials started before 2015 and by the presence of a scientific stop (Fig. 4A–B).

| Trials | Hetionet, no scientific stop | Hetionet, scientific stop | PrimeKG, no scientific stop | PrimeKG, scientific stop |
|---|---|---|---|---|
| 1 | 148/2104 (7%) | 0/51 (0%) | 648/5588 (12%) | 7/131 (5%) |
| 2 | 81/630 (13%) | 1/21 (5%) | 314/1543 (20%) | 11/70 (16%) |
| 3–4 | 84/470 (18%) | 2/43 (5%) | 306/1122 (27%) | 15/92 (16%) |
| 5–9 | 107/333 (32%) | 15/56 (27%) | 324/792 (41%) | 30/113 (27%) |
| 10–19 | 82/164 (50%) | 17/51 (33%) | 189/344 (55%) | 35/96 (36%) |
| ≥20 | 116/147 (79%) | 106/167 (63%) | 172/232 (74%) | 129/205 (63%) |

**Table S20c.** Recorded treatments or approved indications among written-back pairs (stopped trial started before 2015; palliative or off-label pairs excluded), by number of registered trials of any status started before 2015 (Fig. 4C).

| Trials | Hetionet, all stopped | Hetionet, scientific stop | Hetionet, strictest set | PrimeKG, all stopped | PrimeKG, scientific stop | PrimeKG, strictest set |
|---|---|---|---|---|---|---|
| 1 | 22/465 | 1/57 | 1/23 | 65/1161 | 8/145 | 1/69 |
| 2 | 21/251 | 1/19 | 0/9 | 76/574 | 10/67 | 3/30 |
| 3–4 | 41/286 | 3/47 | 0/11 | 128/639 | 10/91 | 4/38 |
| 5–9 | 89/345 | 10/65 | 2/20 | 218/720 | 35/138 | 7/45 |
| 10–19 | 82/202 | 18/51 | 6/11 | 194/436 | 36/106 | 14/40 |
| ≥20 | 276/380 | 126/193 | 26/37 | 409/560 | 165/248 | 61/83 |

Strictest set: same-concept scientific stop of the investigational agent.

**Table S20d.** Pairs with a scientific stop in a trial started before 2015, by phase of the stopped trial (a pair is counted under every phase in which it had a scientific stop; Fig. 4E).

| Phase | Hetionet | PrimeKG |
|---|---|---|
| 1 | 34/71 (48%; 37%–59%) | 43/114 (38%; 29%–47%) |
| 2 | 103/293 (35%; 30%–41%) | 147/486 (30%; 26%–34%) |
| 3 | 64/108 (59%; 50%–68%) | 102/207 (49%; 43%–56%) |
| 4 | 21/32 (66%; 48%–80%) | 40/69 (58%; 46%–69%) |

**Table S20e.** Degree-preserving permutation null for external approved indications among unlabelled stopped pairs (1,000 permutations; Fig. 3).

| Set | Graph | Pairs | Observed | Expected (95% range) | Observed / expected (95% range) |
|---|---|---|---|---|---|
| All stopped pairs | Hetionet | 1,553 | 155 | 24.9 (17–33) | 6.2 (4.7–9.1) |
| All stopped pairs | PrimeKG | 3,319 | 319 | 54.1 (42–66) | 5.9 (4.8–7.6) |
| Scientific stop | Hetionet | 313 | 40 | 7.4 (3–12) | 5.4 (3.3–13.3) |
| Scientific stop | PrimeKG | 600 | 69 | 12.3 (7–19) | 5.6 (3.6–9.9) |
| … same concept | Hetionet | 160 | 12 | 1.9 (0–4) | 6.4 (3.0–∞) |
| … same concept | PrimeKG | 480 | 52 | 8.8 (4–14) | 5.9 (3.7–13.0) |
| … investigational agent | Hetionet | 164 | 20 | 4.0 (1–8) | 5.0 (2.5–20.0) |
| … investigational agent | PrimeKG | 318 | 40 | 5.7 (2–10) | 7.0 (4.0–20.0) |
| … both | Hetionet | 80 | 4 | 1.1 (0–3) | 3.7 (1.3–∞) |
| … both | PrimeKG | 246 | 31 | 4.1 (1–8) | 7.5 (3.9–31.0) |
| Other stops only | Hetionet | 1,240 | 115 | 19.3 (12–27) | 6.0 (4.3–9.6) |
| Other stops only | PrimeKG | 2,719 | 250 | 44.3 (33–56) | 5.6 (4.5–7.6) |

**Table S20f.** Scorers built only from registry history on the external grid (E2). Per-disease and pooled AP with 95% disease-cluster bootstrap intervals; E3 AUROC for approved indications versus later scientific failures.

| Scorer | Graph | Per-disease AP | Pooled AP | E3 AUROC |
|---|---|---|---|---|
| Number of trials before 2015 | Hetionet | 0.268 (0.203–0.338) | 0.178 (0.109–0.250) | 0.603 |
| Number of trials before 2015 | PrimeKG | 0.234 (0.199–0.267) | 0.088 (0.064–0.117) | 0.551 |
| Number of stopped trials before 2015 | Hetionet | 0.146 (0.094–0.205) | 0.085 (0.041–0.138) | 0.656 |
| Number of stopped trials before 2015 | PrimeKG | 0.095 (0.074–0.120) | 0.038 (0.023–0.056) | 0.592 |
| Any stopped trial before 2015 | Hetionet | 0.049 (0.033–0.065) | 0.033 (0.019–0.050) | 0.656 |
| Any stopped trial before 2015 | PrimeKG | 0.060 (0.045–0.080) | 0.018 (0.013–0.024) | 0.592 |

External approved indications with any registered trial before 2015: 58% (Hetionet) and 43% (PrimeKG); other unlabelled pairs: 2.1% and 0.4%.

**Table S20g.** Median percentile rank on the external grid (scorers fitted without write-back) of unlabelled written-back pairs with a scientific stop, by number of registered trials before 2015, and of external approved indications.

| Scorer | Graph | ≤2 trials (n) | ≥10 trials (n) | External approved indications |
|---|---|---|---|---|
| Graph head | Hetionet | 50 (75) | 77 (133) | 76 |
| Graph head | PrimeKG | 63 (197) | 93 (200) | 87 |
| Label-only MF | Hetionet | 64 (75) | 96 (133) | 88 |
| Label-only MF | PrimeKG | 69 (197) | 94 (200) | 91 |
| Degree reference | Hetionet | 78 (75) | 91 (133) | 79 |
| Degree reference | PrimeKG | 67 (197) | 85 (200) | 77 |

**Per-disease heterogeneity (Hetionet, graph head).** Across diseases, the change in per-disease AP under flat negatives correlated with the share of the disease's held-out treatments that the policy negated (Spearman ρ = -0.65 in the random-edge task, 72 diseases; ρ = -0.35 in the compound-disjoint task, 73 diseases). Diseases in which at least half of the held-out treatments were negated lost a median of 0.136 and 0.147, against 0.001 and 0.083 for the other diseases.

## S21. Review of indications with the cleanest failure records

Cases were recorded treatments or approved indications with a same-concept efficacy or safety stop of the investigational agent in a trial started before 2015: all 40 in Hetionet (H01–H40) and a random 40 of 101 in PrimeKG (P01–P40; random seed 20261007). For each trial we retrieved the official title, conditions, phase, allocation, enrolment, arm groups and interventions, primary outcome, eligibility and stop reason from the ClinicalTrials.gov API (66 distinct trials). Some trials map to both graphs, so cases are not independent.

**Coding and adjudication.** Each dossier was coded twice and independently against the codebook below; the second coding was made without access to the first. The two codings agreed on the primary code in 82% of cases (Cohen's κ = 0.78) and on the contradiction judgement in 95% (κ = 0.74). Every case was then read again against the full dossier. Where the dossier did not settle the question, the posted results on ClinicalTrials.gov or the primary publication were consulted: for galantamine (NCT00679627) the posted participant flow reports 41 deaths with placebo and 29 with galantamine, and the EVOLVE trial of tobramycin inhalation powder (NCT00125346) was terminated early on positive interim results (Konstan et al., Pediatr Pulmonol 2011; doi:10.1002/ppul.21356). The final codes differ from the first coding in 7 cases; the column Verification gives the reason for each change. Both codings and the adjudication were carried out with a large language model (Claude, Anthropic) instructed to use the dossier and the codebook; the dossiers (paper_results/case_review/dossiers.md), both codings and the final codes are released so that they can be checked by clinical reviewers.

Codebook (one primary code, up to two secondary codes):

- A1 Special population: subgroup defined by age, pregnancy, comorbidity or region.
- A2 Stage, line or goal: different disease stage, line of therapy or treatment goal (refractory or relapsed disease, adjuvant or neoadjuvant, maintenance, prophylaxis, acute or perioperative setting).
- B1 Regimen of the drug itself: different dose, schedule, duration, route or formulation.
- B2 Combination or add-on: the drug was part of a new multi-drug regimen or was added to another therapy, including cases in which the drug was background therapy or a fixed partner of the agent that failed or caused harm.
- C Comparative question: comparison with another active treatment or strategy, including cases in which the drug belonged to the comparator regimen.
- D Endpoint or subtype: a specific complication, outcome or narrow subtype instead of treatment of the disease.
- E Stop not attributable to the drug's efficacy or safety: the stated reason reports no efficacy or safety finding about the tested regimen (planned interim analysis without futility, accrual, ethics, sponsor decision, external evidence, stop for benefit, zero enrolment) or concerns harm from a comparator.
- F Other, including a failure of the drug in its established setting and ontology mismatches.
- Contradicts: yes if the trial tested the drug in its established population, setting and regimen class and stopped for lack of efficacy or for harm attributable to the drug; unclear if the dossier and posted information do not settle this; otherwise no.
- Drug role check: whether the drug was part of the tested regimen, background or control therapy, or a fixed partner of a new agent.

**Table S21.** Final codes for the 80 reviewed cases, with the first and the independent second coding (primary code / contradicts).

| Case | Drug | Disease | Primary | Secondary | Contradicts | Drug role | First coding | Second coding | Trials | Justification | Verification |
|---|---|---|---|---|---|---|---|---|---|---|---|
| H01 | Ramipril | hypertension | E | A1,B1 | no | investigational (tested regimen) | E / no | E / no | NCT00389519 | Planned interim analysis permitted by protocol; no futility or harm stated; paediatric dose-ranging. | agree |
| H02 | Nicotine | nicotine dependence | A1 | B1 | no | investigational (tested regimen) | A1 / no | A1 / no | NCT00046813, NCT00115687 | Schizophrenia (dose comparison) and pregnant smokers. | agree |
| H03 | Pyrimethamine | malaria | A1 | B2 | unclear | investigational (tested regimen) | A1 / unclear | A1 / unclear | NCT00453856 | Early treatment failure of sulfadoxine-pyrimethamine in young children in Gabon; regional resistance; whether this speaks against the indication in general is unclear. | agree |
| H04 | Dapsone | malaria | E | B2,A1 | no | investigational (tested regimen) | E / no | E / no | NCT00361114 | Sulfadoxine-pyrimethamine arms stopped; the dapsone-containing arm was never given (drug unavailable). | agree |
| H05 | Metoprolol | hypertension | D | E,A2 | no | investigational (tested regimen) | D / no | E / no | NCT00491387 | Mechanistic MIBG study in controlled hypertension; stopped because of external reports on beta-blockers as first-line therapy, not trial efficacy. | agree |
| H06 | Metformin | polycystic ovary syndrome | A2 | B2,D | no | investigational (tested regimen) | A2 / no | A2 / no | NCT01208740 | Adjunct to gonadotropins for IVF in poor responders aged 35-45; outcome cycle cancellation. | agree |
| H07 | Metformin | metabolic syndrome X | A1 | A2,D | no | investigational (tested regimen) | A1 / no | A1 / no | NCT01996696 | Prevention of weight gain in men with prostate cancer on androgen deprivation. | agree |
| H08 | Olanzapine | endogenous depression | F | B2 | yes | investigational (tested regimen) | A2 / unclear | B2 / unclear | NCT01687478 | Olanzapine plus fluoxetine versus fluoxetine in treatment-resistant depression is the approved use of the combination; stopped for lack of efficacy at interim. A failure in the established setting. | changed: contradicts unclear->yes; primary A2->F |
| H09 | Gemcitabine | breast cancer | B2 | A2 | no | investigational (tested regimen) | B2 / no | B2 / no | NCT00462865 | Adjuvant gemcitabine-capecitabine-bevacizumab triplet; regimen toxicity and slow accrual. | agree |
| H10 | Epirubicin | stomach cancer | C | E,B2 | no | background or control | C / no | C / no | NCT02076594 | Epirubicin regimen (EOX) is the standard control although both arms are labelled experimental; the experimental low-dose docetaxel regimen failed. | agree (drug role: background or control) |
| H11 | Fluoxetine | endogenous depression | B2 | E,D | no | background or control | B2 / no | B2 / no | NCT01119430 | Fluoxetine was background therapy in both arms; the add-on DU125530 was tested. | agree (drug role: background or control) |
| H12 | Cyclophosphamide | breast cancer | B2 | B1,A2 | no | investigational (tested regimen) | B2 / no | B2 / no | NCT01329627 | Single-arm feasibility of a new sequential regimen; toxicity. | agree |
| H13 | Methotrexate | urinary bladder cancer | A2 | D,B2 | no | investigational (tested regimen) | A2 / no | A2 / no | NCT00005047 | Adjuvant MVAC after cystectomy in p53-altered organ-confined tumours; futility. | agree |
| H14 | Methotrexate | rheumatoid arthritis | B2 | E,A2 | no | background or control | B2 / no | B2 / no | NCT00485589 | Methotrexate was background therapy in all arms; ocrelizumab benefit-risk was unfavourable. | agree (drug role: background or control) |
| H15 | Fluticasone Propionate | asthma | E | D,B2 | no | investigational (tested regimen) | E / no | E / no | NCT00456313 | Withdrawn before enrolment ('lack of data'); mechanistic study. | agree |
| H16 | Niacin | coronary artery disease | E | A1,D | no | investigational (tested regimen) | E / no | E / no | NCT01414166 | Lipid-endpoint study stopped after the external HPS2-THRIVE result. | agree |
| H17 | Galantamine | Alzheimer's disease | E |  | no | investigational (tested regimen) | F / unclear | B1 / yes | NCT00679627 | Posted results: 41 deaths with placebo against 29 with galantamine; the pre-specified mortality imbalance favoured the drug, so the 'safety' category is wrong. | changed: F/unclear -> E/no |
| H18 | Estradiol | prostate cancer | A2 | B1 | no | investigational (tested regimen) | A2 / no | A2 / no | NCT00176644 | Transdermal estradiol in castration-refractory disease after docetaxel. | agree |
| H19 | Estradiol | breast cancer | D | A2,B1 | no | investigational (tested regimen) | D / no | D / no | NCT01083641 | High-dose estradiol in metastatic triple-negative breast cancer. | agree |
| H20 | Diethylpropion | obesity | E | B2,A1 | no | investigational (tested regimen) | E / no | E / no | NCT00115063 | Multicomponent obesity programme stopped for the ethics of an untreated control group. | agree |
| H21 | Salmeterol | asthma | E | D,B2 | no | investigational (tested regimen) | E / no | E / no | NCT00456313 | Same withdrawn mechanistic study as H15. | agree |
| H22 | Exemestane | breast cancer | C | A2,B2 | no | investigational (tested regimen) | C / no | A2 / no | NCT01303679 | Maintenance switch to exemestane plus bevacizumab versus continued paclitaxel plus bevacizumab; no difference. | agree |
| H23 | Doxorubicin | urinary bladder cancer | A2 | D,B2 | no | investigational (tested regimen) | A2 / no | A2 / no | NCT00005047 | Same adjuvant MVAC trial as H13. | agree |
| H24 | Doxorubicin | breast cancer | B2 | A2 | no | investigational (tested regimen) | B2 / no | B2 / no | NCT01329627, NCT02131506 | Sequential triplet feasibility and phase 1b liposomal doxorubicin plus lapatinib; toxicity. | agree |
| H25 | Letrozole | breast cancer | B2 | E,A2 | no | fixed partner of a new agent | B2 / no | B2 / no | NCT01199367 | Phase 1 triplet with KW-2450; no tolerable dose of the new drug in combination. | agree (drug role: fixed partner of a new agent) |
| H26 | Cinacalcet | chronic kidney failure | A1 | D | no | investigational (tested regimen) | A1 / no | A1 / no | NCT01277510 | Children on dialysis; clinical hold after a fatality; the established use is in adults. | agree |
| H27 | Orlistat | obesity | E | B2,A1 | no | investigational (tested regimen) | E / no | E / no | NCT00115063 | Same obesity programme as H20. | agree |
| H28 | Capecitabine | breast cancer | B2 | A2 | no | investigational (tested regimen) | B2 / no | B2 / no | NCT00462865, NCT01200212 | Capecitabine added to taxane plus bevacizumab first-line, and an adjuvant triplet; no benefit or toxicity of the combinations. | agree |
| H29 | Sibutramine | obesity | E | B2,A1 | no | investigational (tested regimen) | E / no | E / no | NCT00115063 | Same obesity programme as H20. | agree |
| H30 | Pioglitazone | type 2 diabetes mellitus | A1 | C,D | no | investigational (tested regimen) | A1 / no | A1 / no | NCT00521820 | Patients with diabetes and heart failure, versus glyburide; more heart-failure hospitalisations (a known label restriction). | agree |
| H31 | Escitalopram | endogenous depression | E | A1 | no | investigational (tested regimen) | E / no | E / no | NCT00387348 | Patients with advanced cancer; DSMB stopped because the placebo arm had more adverse events. | agree |
| H32 | Paclitaxel | breast cancer | B2 | B1,A2 | no | investigational (tested regimen) | B2 / no | B2 / no | NCT01329627 | Same feasibility triplet as H12. | agree |
| H33 | Aripiprazole | schizophrenia | E | A1,A2 | no | investigational (tested regimen) | E / no | E / no | NCT01122927 | Adolescent open-label safety study; stopped because the paediatric investigational plan objective was met. | agree |
| H34 | Docetaxel | prostate cancer | B2 | E,A2 | no | investigational (tested regimen) | B2 / no | E / no | NCT00494338 | Docetaxel plus celecoxib in castration-resistant disease; stopped for celecoxib safety issues. | agree |
| H35 | Docetaxel | stomach cancer | C | B2,B1 | no | investigational (tested regimen) | C / no | C / no | NCT02076594 | Low-dose docetaxel triplet failed against the epirubicin triplet. | agree |
| H36 | Docetaxel | breast cancer | A1 | B1,D | no | investigational (tested regimen) | A1 / no | A1 / no | NCT00104624 | Biweekly docetaxel in women aged 70 years or older; toxicity. | agree |
| H37 | Lapatinib | breast cancer | B2 | E,A2 | no | investigational (tested regimen) | B2 / no | B2 / no | NCT01199367, NCT02131506 | Two phase 1 combinations; one stopped for the new partner drug, one for combination toxicity. | agree |
| H38 | Everolimus | breast cancer | A2 | B2,E | no | investigational (tested regimen) | A2 / no | B2 / no | NCT00674414, NCT02236572 | Neoadjuvant add-on to trastuzumab (accrual) and to an aromatase inhibitor in low-risk tumours (response gate missed). | agree |
| H39 | Degarelix | prostate cancer | A2 | E | no | investigational (tested regimen) | A2 / no | A2 / no | NCT01545882 | Non-metastatic castration-refractory disease after total androgen blockade; five patients. | agree |
| H40 | Fingolimod | multiple sclerosis | D | A2 | no | investigational (tested regimen) | D / no | D / no | NCT00731692 | Primary progressive multiple sclerosis; disability-progression endpoint missed. | agree |
| P01 | Valsartan | congestive heart failure | E | C,D | no | investigational (tested regimen) | E / no | E / no | NCT01035255 | Sacubitril-valsartan versus enalapril stopped early for compelling efficacy; the 'Negative' category is wrong. | agree |
| P02 | Ramipril | hypertensive disorder | E | A1,B1 | no | investigational (tested regimen) | E / no | E / no | NCT00389519 | Same paediatric ramipril trial as H01. | agree |
| P03 | Nicotine | nicotine dependence | A1 | B1 | no | investigational (tested regimen) | A1 / no | A1 / no | NCT00046813, NCT00115687 | Same two trials as H02. | agree |
| P04 | Pyrimethamine | malaria | A1 | B2 | unclear | investigational (tested regimen) | A1 / unclear | A1 / unclear | NCT00453856 | Same trial as H03. | agree |
| P05 | Metoprolol | hypertensive disorder | D | E,A2 | no | investigational (tested regimen) | D / no | E / no | NCT00491387 | Same mechanistic study as H05. | agree |
| P06 | Lidocaine | osteoarthritis | E | C,B1 | no | investigational (tested regimen) | E / no | E / no | NCT00904605 | Lidocaine patch versus celecoxib; stopped for safety concerns about the COX-2 comparator class. | agree |
| P07 | Mitomycin | pancreatic adenocarcinoma | E | D,A2 | no | investigational (tested regimen) | E / no | E / no | NCT00386399 | Withdrawn before enrolment; no consented patient carried the required BRCA2 mutation. | agree |
| P08 | Gefitinib | non-small cell lung carcinoma | A2 | D | no | investigational (tested regimen) | A2 / no | A2 / no | NCT00104728 | Single-agent neoadjuvant gefitinib in unselected resectable early-stage disease; response-rate stopping rule. | agree |
| P09 | Olanzapine | major depressive disorder | F | B2 | yes | investigational (tested regimen) | A2 / unclear | B2 / unclear | NCT01687478 | Same trial as H08: failure of olanzapine plus fluoxetine in its approved setting. | changed: contradicts unclear->yes; primary A2->F |
| P10 | Olanzapine | anxiety disorder | F | A2,B2 | no | investigational (tested regimen) | F / no | F / no | NCT01687478 | Depression trial mapped to an anxiety-disorder concept (ontology mismatch). | agree |
| P11 | Treprostinil | pulmonary arterial hypertension | A1 | B1,E | no | investigational (tested regimen) | A1 / no | A1 / no | NCT00494533 | Intravenous treprostinil in India; stopped for safety problems of outpatient intravenous infusion in that setting. | agree |
| P12 | Sorafenib | liver cancer | B2 | E | no | fixed partner of a new agent | B2 / no | E / no | NCT01334710 | Sorafenib plus OSI-906; stopped because of a safety issue with OSI-906 in another study. | agree (drug role: fixed partner of a new agent) |
| P13 | Sorafenib | hepatocellular carcinoma | B2 | E,A2 | no | investigational (tested regimen) | B2 / no | B2 / no | NCT00971126, NCT01334710 | Sorafenib with thalidomide (dose-limiting toxicity of the combination) and with OSI-906. | agree |
| P14 | Lansoprazole | duodenal ulcer | E | A2,C | no | investigational (tested regimen) | E / no | E / no | NCT00762359 | Lansoprazole superior to gefarnate for ulcer prevention; stopped for success. | agree |
| P15 | Fluoxetine | major depressive disorder | B2 | E,D | no | background or control | B2 / no | B2 / no | NCT01119430 | Same trial as H11. | agree (drug role: background or control) |
| P16 | Ritonavir | HIV infectious disease | A2 | B2,C | no | investigational (tested regimen) | A2 / no | B2 / no | NCT00090779, NCT00420355, NCT00531986, NCT01602822 | Boosted protease-inhibitor strategy trials: immediate versus deferred therapy, monotherapy simplification, PI pairing, post-exposure prophylaxis. | agree |
| P17 | Erlotinib | non-small cell lung carcinoma | B2 | A2,D | unclear | investigational (tested regimen) | B2 / no | B2 / unclear | NCT00281021, NCT00283634, NCT00554775, NCT01115803 | Erlotinib combinations and whole-brain radiotherapy; one trial (erlotinib alone and with bortezomib, relapsed disease) stopped for insufficient efficacy without saying which arm, and erlotinib alone in relapsed disease is the established use. | changed: contradicts no -> unclear (after blind coding) |
| P18 | Cyclophosphamide | plasma cell myeloma | B2 | E,A2 | no | fixed partner of a new agent | B2 / no | E / no | NCT01881789 | Oral cyclophosphamide arm of an oprozomib dose escalation; halted to optimise the oprozomib formulation. | agree (drug role: fixed partner of a new agent) |
| P19 | Propranolol | hemangioma | E | C,A1 | no | investigational (tested regimen) | E / no | E / no | NCT00967226 | Propranolol versus prednisolone in infants; stopped for growth retardation with prednisolone. | agree |
| P20 | Fluticasone propionate | asthma | E | D,B2 | no | investigational (tested regimen) | E / no | E / no | NCT00456313 | Same withdrawn mechanistic study as H15. | agree |
| P21 | Niacin | myocardial infarction | B2 | D,A2 | yes | investigational (tested regimen) | B2 / yes | B2 / yes | NCT00120289 | AIM-HIGH: extended-release niacin added to simvastatin in established vascular disease; stopped for lack of efficacy. Statin background is current practice, so this is a failure in the setting where niacin would be used. | agree |
| P22 | Prednisone | primary cutaneous T-cell non-Hodgkin lymphoma | B2 |  | no | unclear (no arm data) | B2 / unclear | B2 / no | NCT00161239 | Single-arm multi-agent LAMPP regimen stopped for toxicity of the regimen; as for other combination-toxicity cases, this does not speak against the component. | changed: contradicts unclear -> no |
| P23 | Tobramycin | cystic fibrosis | E | B1,A1 | no | investigational (tested regimen) | E / unclear | E / unclear | NCT00125346 | EVOLVE: published report (Konstan et al., Pediatr Pulmonol 2011, doi:10.1002/ppul.21356) states the trial was terminated early on positive interim results; the 'Negative' category is wrong. | changed: contradicts unclear -> no |
| P24 | Estradiol | breast cancer | D | A2,B1 | no | investigational (tested regimen) | D / no | D / no | NCT01083641 | Same triple-negative trial as H19. | agree |
| P25 | Diethylpropion | obesity disorder | E | B2,A1 | no | investigational (tested regimen) | E / no | E / no | NCT00115063 | Same obesity programme as H20. | agree |
| P26 | Acetylsalicylic acid | transient ischemic attack | B2 | C,A2 | no | background or control | B2 / no | B2 / no | NCT01661322 | TARDIS: intensive triple antiplatelet therapy versus guideline therapy, in which aspirin was also used; stopped on reaching a definitive answer. | agree (drug role: background or control) |
| P27 | Methylprednisolone | non-Hodgkin lymphoma | B2 | A2 | no | investigational (tested regimen) | B2 / no | A2 / no | NCT00367497 | Rituximab plus ESHAP in relapsed or refractory aggressive lymphoma; low response rate. | agree |
| P28 | Pioglitazone | type 2 diabetes mellitus | A1 | C,D | no | investigational (tested regimen) | A1 / no | A1 / no | NCT00521820 | Same pioglitazone heart-failure trial as H30. | agree |
| P29 | Escitalopram | major depressive disorder | E | A1 | no | investigational (tested regimen) | E / no | E / no | NCT00387348 | Same escitalopram trial as H31. | agree |
| P30 | Dexamethasone | mantle cell lymphoma | B2 | A2 | no | fixed partner of a new agent | B2 / no | B2 / no | NCT01578343 | Vorinostat added to fludarabine-mitoxantrone-dexamethasone in relapsed or refractory mantle cell lymphoma; low response. | agree |
| P31 | Docetaxel | non-small cell lung carcinoma | B2 | B1,A1 | no | investigational (tested regimen) | B2 / no | B2 / no | NCT00801801 | Metronomic docetaxel plus sorafenib first-line in performance-status-2 patients; funding withdrawn and early efficacy not encouraging. | agree |
| P32 | Docetaxel | breast cancer | A1 | B1,D | no | investigational (tested regimen) | A1 / no | A1 / no | NCT00104624 | Same biweekly docetaxel trial as H36. | agree |
| P33 | Lisdexamfetamine | attention deficit-hyperactivity disorder | E | A1,D | no | investigational (tested regimen) | E / no | E / no | NCT01017263 | Children with ADHD, obesity and glucose intolerance; stopped for screen failures. | agree |
| P34 | Lapatinib | breast cancer | B2 | E,A2 | no | investigational (tested regimen) | B2 / no | B2 / no | NCT01199367, NCT02131506 | Same two phase 1 combinations as H37. | agree |
| P35 | Everolimus | breast cancer | A2 | B2,E | no | investigational (tested regimen) | A2 / no | B2 / no | NCT00674414, NCT02236572 | Same neoadjuvant trials as H38. | agree |
| P36 | Trabectedin | soft tissue sarcoma | C | B2,A2 | no | investigational (tested regimen) | C / no | A2 / no | NCT01104298, NCT01189253 | First-line trabectedin (alone or with doxorubicin) versus doxorubicin; the established use is after anthracycline failure. | agree |
| P37 | Prasugrel | myocardial infarction | E | C,D | no | investigational (tested regimen) | E / no | E / no | NCT00910299 | TRIGGER-PCI: prasugrel versus clopidogrel after elective PCI; stopped for a low event rate. | agree |
| P38 | Cangrelor | acute coronary syndrome | C | A2 | yes | investigational (tested regimen) | C / no | C / unclear | NCT00305162 | CHAMPION PCI: cangrelor versus clopidogrel during PCI, the comparison on which cangrelor was later approved, stopped for insufficient effectiveness. A failure in the established setting (approval followed a later trial). | changed: contradicts no -> yes |
| P39 | Degarelix | prostate cancer | A2 | E | no | investigational (tested regimen) | A2 / no | A2 / no | NCT01545882 | Same degarelix trial as H39. | agree |
| P40 | Rilpivirine | HIV infectious disease | A2 | B1,E | no | investigational (tested regimen) | A2 / no | A2 / no | NCT01049932 | Long-acting injectable rilpivirine as pre-exposure prophylaxis in HIV-negative volunteers; no participant enrolled. | agree |
