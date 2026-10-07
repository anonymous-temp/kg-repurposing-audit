## S16. Drug roles in stopped trials

Arm groups and interventions were retrieved for all 29,485 terminated, withdrawn or suspended trials in Open Targets 26.09 from the ClinicalTrials.gov API version 2 on 6 October 2026. Requests used the endpoint /api/v2/studies in batches of 500 identifiers, with the fields NCTId, ArmGroup, Intervention, Condition, Phase, OverallStatus, WhyStopped, StartDate, DesignAllocation, HasResults and LeadSponsorName. All 29,485 records were returned; 27,232 had arm groups.

Names of a compound were its Hetionet or PrimeKG name and the name, synonyms and trade names of every ChEMBL record mapped to it, including parent and salt forms. Names were lower-cased, non-alphanumeric characters were replaced by spaces, and names shorter than four characters or consisting of generic words (for example "sodium", "tablet" or "placebo") were dropped. A compound matched an intervention when one of its names occurred as a whole phrase in the intervention name or its other names. Interventions whose name contained "placebo" were used only when the name also combined placebo with other agents (for example "Placebo + Paclitaxel + Carboplatin"). A dummy such as "placebo for sorafenib" therefore did not place sorafenib in the control arm. A compound that occurred only in intervention or arm descriptions was classified as "described only" and was not treated as investigational. Matched interventions were linked to arm groups through their arm-group labels.

Table S13. Role of the mapped compound in stopped trials (trial–compound combinations).

| Role | Rule | Hetionet | PrimeKG |
|---|---|---|---|
| Investigational agent | only in experimental arms; or the only drug intervention of a single-arm trial or of a trial without arm data | 3,750 | 7,469 |
| Single-arm combination | one arm with two or more drug interventions | 2,553 | 3,390 |
| Background therapy | in every arm | 2,177 | 3,538 |
| Comparator | only in control arms | 438 | 719 |
| Head-to-head | no experimental arm; in some but not all active arms | 376 | 676 |
| Mixed | in experimental and some control arms | 143 | 222 |
| Described only | only in intervention or arm descriptions | 332 | 595 |
| No arm data | several drug interventions, no arm groups | 666 | 873 |
| Unmatched | no name match | 351 | 692 |
| Total | | 10,786 | 18,174 |

Control arms were active-comparator, placebo-comparator, sham-comparator and no-intervention arms, together with other non-experimental arms in trials that had an experimental arm. In trials without an experimental arm, active-comparator arms were treated as test arms.

A first version of the rules excluded every intervention whose name contained "placebo" and fell back on intervention descriptions when no name matched. A check of 50 randomly sampled assignments among Hetionet scientific stops before 2015 found two systematic errors in that version. Background drugs in arms such as "Placebo + Paclitaxel + Carboplatin" or "placebo + nicotine replacement therapy" were classified as investigational. So was methotrexate mentioned in the description of a rituximab intervention. The rules above correct both errors; all reported results use the corrected rules. After correction, 46 of the 50 sampled assignments agreed with the registry record. Four disagreements remained, and all four placed a drug in a non-investigational class. One registry entry labelled the atorvastatin arms of a factorial trial as placebo comparators. Three interventions named the drug only by a code, an abbreviation or a class (TAS-102 for trifluridine and tipiracil, VPA for valproic acid and "Taxane" for paclitaxel), so the drug was unmatched or described only. None of the 50 assignments placed a comparator or background drug in the investigational class. The check was carried out during the analysis with the AI assistance described in Section S15. It is an internal check; an independent review of a larger sample would be needed to estimate error rates.

Table S14. Recorded treatments or approved indications among pairs with a scientific stop before 2015, by the role of the compound.

| Pairs | Hetionet: pairs | Hetionet: recorded or approved | PrimeKG: pairs | PrimeKG: recorded or approved |
|---|---|---|---|---|
| Scientific stop, at least one as investigational agent | 234 | 85 (36%) | 503 | 114 (23%) |
| Scientific stop, only as comparator or background therapy | 63 | 30 (48%) | 188 | 66 (35%) |
| Scientific stop, only other roles | 162 | 64 (40%) | 319 | 106 (33%) |
| Same-concept scientific stop, at least one as investigational agent | 120 | 40 (33%) | 408 | 101 (25%) |
| Same-concept scientific stop, only as comparator or background therapy | 44 | 24 (55%) | 157 | 58 (37%) |

## S17. Replication in PrimeKG

PrimeKG version 2 (Harvard Dataverse, doi:10.7910/DVN/IXA7BM) was downloaded on 6 October 2026: nodes.csv (129,375 nodes) and edges.csv (8,100,498 rows, each undirected edge in both directions). Each undirected edge was used once in each direction of DistMult training (4,007,618 triples after removal of 42,631 drug–disease edges: indications, contraindications and off-label uses). Training used the Hetionet recipe with a run configured for 15 epochs. The operating system stopped the run after epoch 13 because of memory pressure from other jobs on the host. The checkpoints after 5 and 10 epochs were complete; they are identical to those of a 10-epoch run with the same seed, and the analysis used them. One PrimeKG epoch has about 1.8 times as many triples as a Hetionet epoch, so 10 epochs correspond to about 80 million training triples, against 67 million for the 15-epoch Hetionet checkpoint.

Open Targets conditions were mapped to the 1,363 PrimeKG diseases with at least one indication. Of the condition–disease assignments, 1,604 were direct (MONDO identifier or member of a grouped disease), 47 used a shared DOID, OMIM, Orphanet or UMLS cross-reference and 10,153 used the most specific mapped ancestor. In total, 11,145 Open Targets condition terms reached 1,104 of the 1,363 diseases. Compounds were mapped through DrugBank cross-references (5,400 of 7,957 PrimeKG drugs). The grid comprised the 1,801 drugs and 1,363 diseases with at least one indication. Partitions used seeds 42, 1 and 2. The random-edge task held out 1,878 indications per partition (889–899 test compounds; evaluation grids of 1.21–1.22 million pairs). The compound-disjoint task held out all indications of 360 compounds (1,643–2,003 test positives; 490,318–490,423 pairs).

Hyper-parameters were selected on the validation grid of the seed-42 partition of each task, without write-back (Table S15). MF used dimension 32 with weight decay 10⁻⁵, 10⁻⁴ or 10⁻³. The graph head used the 5- or 10-epoch checkpoint with the same weight decays. Both were evaluated after 200, 400 and 800 epochs. Bootstrap draws (1,000, resampling the 1,363 diseases) were generated once per outcome set and applied to every partition, policy and scorer. Pooled AP under each draw was computed from the sorted scores and the per-disease counts of ranked pairs. This reproduces the stored-score implementation exactly, as checked on simulated data with continuous and tied scores (repository test test_roles_and_fastboot.py).

Table S15. PrimeKG hyper-parameter selection (validation AP on the seed-42 partition, without write-back).

| Task | Scorer | MF dimension | Embedding checkpoint (epochs) | Weight decay | Training epochs | Validation AP |
|---|---|---|---|---|---|---|
| Random edge | Label-only MF | 32 | – | 10⁻⁵ | 800 | 0.231 |
| Random edge | Graph head | – | 10 | 10⁻⁵ | 400 | 0.021 |
| Compound disjoint | Label-only MF | 32 | – | 10⁻³ | 200 | 0.016 |
| Compound disjoint | Graph head | – | 10 | 10⁻⁵ | 800 | 0.053 |

Table S16a. PrimeKG results for every policy and scorer: means (SD) over three partitions for E1; single fits for E2. MF, label-only matrix factorisation.

| Set | Policy | Scorer | Pooled AP | Per-disease AP | MRR | Hits@10 | AUROC | Held-out negated | E3 AUROC |
|---|---|---|---|---|---|---|---|---|---|
| E1 compound disjoint | No write-back | Degree reference | 0.017 (0.002) | 0.007 (0.001) | 0.047 (0.006) | 0.086 (0.008) | 0.766 (0.007) | 0.0 |  |
| E1 compound disjoint | No write-back | Graph head | 0.045 (0.018) | 0.144 (0.022) | 0.083 (0.024) | 0.161 (0.054) | 0.850 (0.009) | 0.0 |  |
| E1 compound disjoint | No write-back | Label-only MF | 0.017 (0.001) | 0.007 (0.001) | 0.047 (0.007) | 0.089 (0.010) | 0.758 (0.002) | 0.0 |  |
| E1 compound disjoint | Flat negative | Degree reference | 0.017 (0.002) | 0.007 (0.001) | 0.047 (0.006) | 0.086 (0.008) | 0.766 (0.007) | 148.0 |  |
| E1 compound disjoint | Flat negative | Graph head | 0.040 (0.013) | 0.104 (0.010) | 0.084 (0.015) | 0.165 (0.022) | 0.817 (0.023) | 148.0 |  |
| E1 compound disjoint | Flat negative | Label-only MF | 0.017 (0.001) | 0.007 (0.000) | 0.047 (0.007) | 0.089 (0.010) | 0.737 (0.011) | 148.0 |  |
| E1 compound disjoint | Typed | Degree reference | 0.017 (0.002) | 0.007 (0.001) | 0.047 (0.006) | 0.086 (0.008) | 0.766 (0.007) | 35.7 |  |
| E1 compound disjoint | Typed | Graph head | 0.053 (0.008) | 0.131 (0.006) | 0.096 (0.008) | 0.181 (0.021) | 0.848 (0.010) | 35.7 |  |
| E1 compound disjoint | Typed | Label-only MF | 0.017 (0.001) | 0.007 (0.001) | 0.047 (0.007) | 0.089 (0.010) | 0.755 (0.005) | 35.7 |  |
| E1 compound disjoint | Typed + scope | Degree reference | 0.017 (0.002) | 0.007 (0.001) | 0.047 (0.006) | 0.086 (0.008) | 0.766 (0.007) | 33.7 |  |
| E1 compound disjoint | Typed + scope | Graph head | 0.044 (0.015) | 0.126 (0.009) | 0.087 (0.025) | 0.163 (0.044) | 0.837 (0.016) | 33.7 |  |
| E1 compound disjoint | Typed + scope | Label-only MF | 0.017 (0.001) | 0.007 (0.001) | 0.047 (0.007) | 0.089 (0.010) | 0.755 (0.004) | 33.7 |  |
| E1 compound disjoint | Typed + role | Degree reference | 0.017 (0.002) | 0.007 (0.001) | 0.047 (0.006) | 0.086 (0.008) | 0.766 (0.007) | 13.0 |  |
| E1 compound disjoint | Typed + role | Graph head | 0.057 (0.011) | 0.150 (0.013) | 0.098 (0.009) | 0.181 (0.021) | 0.859 (0.003) | 13.0 |  |
| E1 compound disjoint | Typed + role | Label-only MF | 0.017 (0.001) | 0.007 (0.000) | 0.047 (0.007) | 0.089 (0.010) | 0.757 (0.003) | 13.0 |  |
| E1 compound disjoint | Typed + role + scope | Degree reference | 0.017 (0.002) | 0.007 (0.001) | 0.047 (0.006) | 0.086 (0.008) | 0.766 (0.007) | 13.0 |  |
| E1 compound disjoint | Typed + role + scope | Graph head | 0.050 (0.020) | 0.139 (0.003) | 0.090 (0.014) | 0.169 (0.024) | 0.849 (0.015) | 13.0 |  |
| E1 compound disjoint | Typed + role + scope | Label-only MF | 0.017 (0.001) | 0.007 (0.000) | 0.047 (0.007) | 0.089 (0.010) | 0.757 (0.003) | 13.0 |  |
| E1 compound disjoint | Mask all | Degree reference | 0.017 (0.002) | 0.007 (0.001) | 0.047 (0.006) | 0.086 (0.008) | 0.766 (0.007) | 0.0 |  |
| E1 compound disjoint | Mask all | Graph head | 0.059 (0.009) | 0.151 (0.009) | 0.100 (0.008) | 0.185 (0.020) | 0.858 (0.011) | 0.0 |  |
| E1 compound disjoint | Mask all | Label-only MF | 0.017 (0.001) | 0.007 (0.001) | 0.047 (0.007) | 0.089 (0.010) | 0.759 (0.002) | 0.0 |  |
| E1 random edge | No write-back | Degree reference | 0.018 (0.001) | 0.114 (0.001) | 0.049 (0.003) | 0.092 (0.000) | 0.796 (0.004) | 0.0 |  |
| E1 random edge | No write-back | Graph head | 0.031 (0.007) | 0.130 (0.005) | 0.106 (0.021) | 0.198 (0.030) | 0.899 (0.012) | 0.0 |  |
| E1 random edge | No write-back | Label-only MF | 0.562 (0.006) | 0.620 (0.006) | 0.634 (0.005) | 0.785 (0.002) | 0.930 (0.002) | 0.0 |  |
| E1 random edge | Flat negative | Degree reference | 0.018 (0.001) | 0.114 (0.001) | 0.049 (0.003) | 0.092 (0.000) | 0.796 (0.004) | 152.0 |  |
| E1 random edge | Flat negative | Graph head | 0.031 (0.006) | 0.125 (0.014) | 0.112 (0.009) | 0.213 (0.023) | 0.899 (0.005) | 152.0 |  |
| E1 random edge | Flat negative | Label-only MF | 0.533 (0.005) | 0.609 (0.001) | 0.608 (0.011) | 0.762 (0.005) | 0.926 (0.003) | 152.0 |  |
| E1 random edge | Typed | Degree reference | 0.018 (0.001) | 0.114 (0.001) | 0.049 (0.003) | 0.092 (0.000) | 0.796 (0.004) | 40.0 |  |
| E1 random edge | Typed | Graph head | 0.030 (0.007) | 0.129 (0.011) | 0.107 (0.019) | 0.208 (0.025) | 0.900 (0.009) | 40.0 |  |
| E1 random edge | Typed | Label-only MF | 0.549 (0.004) | 0.616 (0.003) | 0.617 (0.007) | 0.776 (0.002) | 0.929 (0.002) | 40.0 |  |
| E1 random edge | Typed + scope | Degree reference | 0.018 (0.001) | 0.114 (0.001) | 0.049 (0.003) | 0.092 (0.000) | 0.796 (0.004) | 37.7 |  |
| E1 random edge | Typed + scope | Graph head | 0.034 (0.003) | 0.131 (0.010) | 0.115 (0.001) | 0.215 (0.014) | 0.904 (0.002) | 37.7 |  |
| E1 random edge | Typed + scope | Label-only MF | 0.549 (0.005) | 0.616 (0.003) | 0.619 (0.006) | 0.776 (0.003) | 0.929 (0.002) | 37.7 |  |
| E1 random edge | Typed + role | Degree reference | 0.018 (0.001) | 0.114 (0.001) | 0.049 (0.003) | 0.092 (0.000) | 0.796 (0.004) | 13.7 |  |
| E1 random edge | Typed + role | Graph head | 0.034 (0.003) | 0.133 (0.006) | 0.110 (0.012) | 0.214 (0.021) | 0.899 (0.009) | 13.7 |  |
| E1 random edge | Typed + role | Label-only MF | 0.556 (0.004) | 0.620 (0.006) | 0.625 (0.004) | 0.784 (0.002) | 0.930 (0.002) | 13.7 |  |
| E1 random edge | Typed + role + scope | Degree reference | 0.018 (0.001) | 0.114 (0.001) | 0.049 (0.003) | 0.092 (0.000) | 0.796 (0.004) | 13.3 |  |
| E1 random edge | Typed + role + scope | Graph head | 0.035 (0.001) | 0.132 (0.007) | 0.118 (0.003) | 0.222 (0.005) | 0.906 (0.002) | 13.3 |  |
| E1 random edge | Typed + role + scope | Label-only MF | 0.556 (0.004) | 0.620 (0.006) | 0.626 (0.005) | 0.783 (0.002) | 0.930 (0.002) | 13.3 |  |
| E1 random edge | Mask all | Degree reference | 0.018 (0.001) | 0.114 (0.001) | 0.049 (0.003) | 0.092 (0.000) | 0.796 (0.004) | 0.0 |  |
| E1 random edge | Mask all | Graph head | 0.035 (0.002) | 0.132 (0.009) | 0.117 (0.003) | 0.219 (0.008) | 0.905 (0.002) | 0.0 |  |
| E1 random edge | Mask all | Label-only MF | 0.559 (0.005) | 0.621 (0.005) | 0.628 (0.005) | 0.785 (0.002) | 0.931 (0.002) | 0.0 |  |
| E2 external approvals | No write-back | Degree reference | 0.003 | 0.063 | 0.019 | 0.029 | 0.689 |  | 0.589 |
| E2 external approvals | No write-back | Graph head | 0.005 | 0.050 | 0.043 | 0.096 | 0.766 |  | 0.621 |
| E2 external approvals | No write-back | Label-only MF | 0.020 | 0.147 | 0.116 | 0.229 | 0.746 |  | 0.651 |
| E2 external approvals | Flat negative | Degree reference | 0.003 | 0.063 | 0.019 | 0.029 | 0.690 |  | 0.59 |
| E2 external approvals | Flat negative | Graph head | 0.004 | 0.049 | 0.036 | 0.080 | 0.748 |  | 0.631 |
| E2 external approvals | Flat negative | Label-only MF | 0.013 | 0.122 | 0.092 | 0.195 | 0.698 |  | 0.658 |
| E2 external approvals | Typed | Degree reference | 0.003 | 0.063 | 0.019 | 0.029 | 0.689 |  | 0.589 |
| E2 external approvals | Typed | Graph head | 0.005 | 0.050 | 0.042 | 0.094 | 0.764 |  | 0.625 |
| E2 external approvals | Typed | Label-only MF | 0.020 | 0.146 | 0.129 | 0.229 | 0.745 |  | 0.662 |
| E2 external approvals | Typed + scope | Degree reference | 0.003 | 0.063 | 0.019 | 0.029 | 0.689 |  | 0.589 |
| E2 external approvals | Typed + scope | Graph head | 0.005 | 0.051 | 0.045 | 0.098 | 0.776 |  | 0.607 |
| E2 external approvals | Typed + scope | Label-only MF | 0.020 | 0.147 | 0.130 | 0.235 | 0.750 |  | 0.661 |
| E2 external approvals | Typed + role | Degree reference | 0.003 | 0.063 | 0.019 | 0.029 | 0.689 |  | 0.589 |
| E2 external approvals | Typed + role | Graph head | 0.005 | 0.051 | 0.043 | 0.096 | 0.765 |  | 0.623 |
| E2 external approvals | Typed + role | Label-only MF | 0.021 | 0.146 | 0.131 | 0.235 | 0.750 |  | 0.659 |
| E2 external approvals | Typed + role + scope | Degree reference | 0.003 | 0.063 | 0.019 | 0.029 | 0.689 |  | 0.589 |
| E2 external approvals | Typed + role + scope | Graph head | 0.006 | 0.053 | 0.047 | 0.105 | 0.781 |  | 0.602 |
| E2 external approvals | Typed + role + scope | Label-only MF | 0.021 | 0.148 | 0.131 | 0.234 | 0.754 |  | 0.657 |
| E2 external approvals | Mask all | Degree reference | 0.003 | 0.063 | 0.019 | 0.029 | 0.689 |  | 0.589 |
| E2 external approvals | Mask all | Graph head | 0.005 | 0.051 | 0.044 | 0.099 | 0.767 |  | 0.62 |
| E2 external approvals | Mask all | Label-only MF | 0.022 | 0.154 | 0.136 | 0.239 | 0.756 |  | 0.647 |

Table S16b. PrimeKG contrasts with 95% paired disease-cluster bootstrap intervals (1,000 draws).

| Set | Contrast | Pooled AP Δ (95% CI) | Per-disease AP Δ (95% CI) |
|---|---|---|---|
| E1 random edge | flat_negative - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | flat_negative - no_writeback [mf] | −0.030 (−0.039, −0.020) | −0.011 (−0.016, −0.007) |
| E1 random edge | flat_negative - no_writeback [graph] | 0.000 (−0.002, +0.002) | −0.004 (−0.008, +0.001) |
| E1 random edge | mask_all - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | mask_all - no_writeback [mf] | −0.003 (−0.008, 0.000) | +0.001 (0.000, +0.003) |
| E1 random edge | mask_all - no_writeback [graph] | +0.004 (+0.003, +0.006) | +0.002 (0.000, +0.005) |
| E1 random edge | typed - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | typed - no_writeback [mf] | −0.014 (−0.020, −0.009) | −0.004 (−0.007, −0.001) |
| E1 random edge | typed - no_writeback [graph] | −0.001 (−0.003, +0.001) | −0.001 (−0.004, +0.003) |
| E1 random edge | typed_scoped - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | typed_scoped - no_writeback [mf] | −0.014 (−0.019, −0.009) | −0.004 (−0.007, −0.002) |
| E1 random edge | typed_scoped - no_writeback [graph] | +0.003 (+0.001, +0.005) | +0.001 (−0.002, +0.004) |
| E1 random edge | typed_role - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | typed_role - no_writeback [mf] | −0.006 (−0.011, −0.002) | 0.000 (−0.002, +0.002) |
| E1 random edge | typed_role - no_writeback [graph] | +0.003 (+0.001, +0.005) | +0.003 (0.000, +0.007) |
| E1 random edge | typed_role_scoped - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | typed_role_scoped - no_writeback [mf] | −0.006 (−0.011, −0.002) | 0.000 (−0.002, +0.002) |
| E1 random edge | typed_role_scoped - no_writeback [graph] | +0.004 (+0.003, +0.006) | +0.002 (0.000, +0.005) |
| E1 random edge | typed_role - typed [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | typed_role - typed [mf] | +0.007 (+0.003, +0.012) | +0.004 (+0.001, +0.007) |
| E1 random edge | typed_role - typed [graph] | +0.004 (+0.002, +0.007) | +0.004 (0.000, +0.008) |
| E1 random edge | typed_role - mask_all [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | typed_role - mask_all [mf] | −0.003 (−0.005, −0.001) | −0.001 (−0.003, 0.000) |
| E1 random edge | typed_role - mask_all [graph] | −0.001 (−0.002, 0.000) | +0.001 (−0.002, +0.004) |
| E1 random edge | typed_role_scoped - typed_scoped [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | typed_role_scoped - typed_scoped [mf] | +0.007 (+0.003, +0.012) | +0.004 (+0.002, +0.007) |
| E1 random edge | typed_role_scoped - typed_scoped [graph] | +0.001 (0.000, +0.002) | +0.001 (−0.001, +0.004) |
| E1 random edge | typed_role_scoped - mask_all [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | typed_role_scoped - mask_all [mf] | −0.003 (−0.005, −0.001) | −0.001 (−0.003, 0.000) |
| E1 random edge | typed_role_scoped - mask_all [graph] | 0.000 (0.000, +0.001) | 0.000 (−0.002, +0.002) |
| E1 random edge | graph - mf [no_writeback] | −0.532 (−0.568, −0.494) | −0.491 (−0.513, −0.466) |
| E1 random edge | mf - degree [no_writeback] | +0.545 (+0.508, +0.580) | +0.506 (+0.482, +0.529) |
| E1 random edge | graph - degree [no_writeback] | +0.013 (+0.006, +0.023) | +0.015 (+0.004, +0.028) |
| E1 compound disjoint | flat_negative - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | flat_negative - no_writeback [mf] | −0.001 (−0.001, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | flat_negative - no_writeback [graph] | −0.005 (−0.010, 0.000) | −0.040 (−0.047, −0.032) |
| E1 compound disjoint | mask_all - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | mask_all - no_writeback [mf] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | mask_all - no_writeback [graph] | +0.014 (+0.009, +0.021) | +0.007 (+0.003, +0.011) |
| E1 compound disjoint | typed - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed - no_writeback [mf] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed - no_writeback [graph] | +0.008 (+0.002, +0.015) | −0.013 (−0.018, −0.008) |
| E1 compound disjoint | typed_scoped - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_scoped - no_writeback [mf] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_scoped - no_writeback [graph] | −0.002 (−0.005, +0.002) | −0.018 (−0.024, −0.012) |
| E1 compound disjoint | typed_role - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role - no_writeback [mf] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role - no_writeback [graph] | +0.011 (+0.006, +0.019) | +0.006 (+0.001, +0.010) |
| E1 compound disjoint | typed_role_scoped - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role_scoped - no_writeback [mf] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role_scoped - no_writeback [graph] | +0.005 (0.000, +0.011) | −0.005 (−0.011, +0.001) |
| E1 compound disjoint | typed_role - typed [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role - typed [mf] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role - typed [graph] | +0.004 (+0.002, +0.005) | +0.019 (+0.015, +0.023) |
| E1 compound disjoint | typed_role - mask_all [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role - mask_all [mf] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role - mask_all [graph] | −0.003 (−0.004, −0.001) | −0.001 (−0.004, +0.001) |
| E1 compound disjoint | typed_role_scoped - typed_scoped [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role_scoped - typed_scoped [mf] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role_scoped - typed_scoped [graph] | +0.006 (+0.002, +0.011) | +0.013 (+0.008, +0.019) |
| E1 compound disjoint | typed_role_scoped - mask_all [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role_scoped - mask_all [mf] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role_scoped - mask_all [graph] | −0.009 (−0.012, −0.007) | −0.012 (−0.017, −0.007) |
| E1 compound disjoint | graph - mf [no_writeback] | +0.028 (+0.021, +0.038) | +0.137 (+0.127, +0.148) |
| E1 compound disjoint | mf - degree [no_writeback] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | graph - degree [no_writeback] | +0.028 (+0.021, +0.038) | +0.137 (+0.127, +0.148) |
| E2 external approvals | flat_negative - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | flat_negative - no_writeback [mf] | −0.007 (−0.010, −0.004) | −0.025 (−0.037, −0.015) |
| E2 external approvals | flat_negative - no_writeback [graph] | −0.001 (−0.001, 0.000) | −0.001 (−0.004, +0.003) |
| E2 external approvals | mask_all - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | mask_all - no_writeback [mf] | +0.002 (+0.001, +0.004) | +0.007 (0.000, +0.017) |
| E2 external approvals | mask_all - no_writeback [graph] | 0.000 (0.000, 0.000) | 0.000 (0.000, +0.001) |
| E2 external approvals | typed - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | typed - no_writeback [mf] | 0.000 (−0.001, +0.002) | 0.000 (−0.008, +0.009) |
| E2 external approvals | typed - no_writeback [graph] | 0.000 (0.000, 0.000) | −0.001 (−0.002, 0.000) |
| E2 external approvals | typed_scoped - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | typed_scoped - no_writeback [mf] | 0.000 (−0.001, +0.002) | 0.000 (−0.009, +0.010) |
| E2 external approvals | typed_scoped - no_writeback [graph] | 0.000 (0.000, +0.001) | +0.001 (0.000, +0.002) |
| E2 external approvals | typed_role - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | typed_role - no_writeback [mf] | +0.001 (0.000, +0.003) | −0.001 (−0.010, +0.010) |
| E2 external approvals | typed_role - no_writeback [graph] | 0.000 (0.000, 0.000) | 0.000 (0.000, +0.001) |
| E2 external approvals | typed_role_scoped - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | typed_role_scoped - no_writeback [mf] | +0.001 (0.000, +0.003) | +0.001 (−0.008, +0.011) |
| E2 external approvals | typed_role_scoped - no_writeback [graph] | +0.001 (0.000, +0.002) | +0.003 (+0.001, +0.005) |
| E2 external approvals | typed_role - typed [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | typed_role - typed [mf] | +0.001 (0.000, +0.002) | −0.001 (−0.006, +0.004) |
| E2 external approvals | typed_role - typed [graph] | 0.000 (0.000, 0.000) | +0.001 (0.000, +0.002) |
| E2 external approvals | typed_role - mask_all [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | typed_role - mask_all [mf] | −0.001 (−0.002, 0.000) | −0.008 (−0.015, −0.003) |
| E2 external approvals | typed_role - mask_all [graph] | 0.000 (0.000, 0.000) | 0.000 (−0.001, +0.001) |
| E2 external approvals | typed_role_scoped - typed_scoped [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | typed_role_scoped - typed_scoped [mf] | +0.001 (0.000, +0.002) | +0.001 (−0.002, +0.004) |
| E2 external approvals | typed_role_scoped - typed_scoped [graph] | 0.000 (0.000, +0.001) | +0.002 (+0.001, +0.003) |
| E2 external approvals | typed_role_scoped - mask_all [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | typed_role_scoped - mask_all [mf] | −0.001 (−0.001, 0.000) | −0.006 (−0.013, −0.002) |
| E2 external approvals | typed_role_scoped - mask_all [graph] | +0.001 (0.000, +0.002) | +0.003 (+0.001, +0.004) |
| E2 external approvals | graph - mf [no_writeback] | −0.015 (−0.022, −0.009) | −0.096 (−0.122, −0.071) |
| E2 external approvals | mf - degree [no_writeback] | +0.017 (+0.011, +0.025) | +0.083 (+0.056, +0.112) |
| E2 external approvals | graph - degree [no_writeback] | +0.002 (0.000, +0.004) | −0.013 (−0.028, +0.002) |

## S18. Tested versus approved pairs and role-policy details

The tested-versus-approved analysis used the scorers fitted to all recorded treatments without write-back. Newly tested pairs had a first trial of any status that started in 2017 or later and no approval evidence; untested pairs had no Open Targets clinical report. The untested group was represented by a fixed random sample of 200,000 pairs. The bootstrap resampled diseases (1,000 draws) and weighted each pair by the multiplicity of its disease.

Table S17. AUROC between groups of unlabelled pairs, with 95% disease-cluster bootstrap intervals.

| Graph | Comparison | Degree reference | Label-only MF | Graph head |
|---|---|---|---|---|
| Hetionet | Approved vs untested | 0.70 (0.60, 0.78) | 0.74 (0.63, 0.83) | 0.67 (0.59, 0.74) |
| Hetionet | Later failure vs untested | 0.72 (0.64, 0.80) | 0.74 (0.66, 0.81) | 0.52 (0.44, 0.61) |
| Hetionet | Newly tested vs untested | 0.69 (0.64, 0.72) | 0.66 (0.62, 0.70) | 0.54 (0.52, 0.56) |
| Hetionet | Approved vs newly tested | 0.51 (0.44, 0.58) | 0.60 (0.50, 0.69) | 0.62 (0.55, 0.68) |
| Hetionet | Approved vs later failure | 0.50 (0.38, 0.62) | 0.56 (0.43, 0.68) | 0.66 (0.54, 0.75) |
| Hetionet | n (approved / later failure / newly tested / untested) | 497 / 47 / 1,904 / 203,927 |  |  |
| PrimeKG | Approved vs untested | 0.69 (0.61, 0.75) | 0.75 (0.70, 0.80) | 0.77 (0.74, 0.80) |
| PrimeKG | Later failure vs untested | 0.61 (0.55, 0.67) | 0.62 (0.56, 0.68) | 0.65 (0.60, 0.71) |
| PrimeKG | Newly tested vs untested | 0.65 (0.63, 0.68) | 0.62 (0.61, 0.64) | 0.61 (0.59, 0.63) |
| PrimeKG | Approved vs newly tested | 0.55 (0.47, 0.61) | 0.65 (0.60, 0.70) | 0.65 (0.62, 0.68) |
| PrimeKG | Approved vs later failure | 0.59 (0.51, 0.67) | 0.65 (0.59, 0.71) | 0.62 (0.56, 0.68) |
| PrimeKG | n (approved / later failure / newly tested / untested) | 1,735 / 155 / 5,308 / 2,424,284 |  |  |

Table S18a. Hetionet: mean reciprocal rank of held-out treatments among the candidates of their compound, by stopped trials that started before 2015 and by the role of the compound in scientific stops. n is the mean number of held-out treatments per partition.

| Task | Scorer | Stratum | n | No write-back | Flat | Typed | Typed + scope | Typed + role | Typed + role + scope | Mask all |
|---|---|---|---|---|---|---|---|---|---|---|
| Compound disjoint | Graph head | no stop | 74.4 | 0.447 | 0.379 | 0.428 | 0.442 | 0.439 | 0.442 | 0.445 |
| Compound disjoint | Graph head | other stop | 51.2 | 0.505 | 0.285 | 0.465 | 0.486 | 0.487 | 0.49 | 0.507 |
| Compound disjoint | Graph head | scientific stop, investigational | 13.8 | 0.504 | 0.292 | 0.261 | 0.361 | 0.212 | 0.364 | 0.511 |
| Compound disjoint | Graph head | scientific stop, other role | 10.4 | 0.539 | 0.269 | 0.179 | 0.378 | 0.467 | 0.512 | 0.542 |
| Compound disjoint | Label-only MF | no stop | 74.4 | 0.152 | 0.137 | 0.147 | 0.153 | 0.148 | 0.152 | 0.15 |
| Compound disjoint | Label-only MF | other stop | 51.2 | 0.177 | 0.02 | 0.157 | 0.162 | 0.161 | 0.161 | 0.178 |
| Compound disjoint | Label-only MF | scientific stop, investigational | 13.8 | 0.29 | 0.025 | 0.018 | 0.089 | 0.016 | 0.119 | 0.295 |
| Compound disjoint | Label-only MF | scientific stop, other role | 10.4 | 0.26 | 0.021 | 0.017 | 0.177 | 0.224 | 0.253 | 0.263 |
| Random edge | Graph head | no stop | 78.0 | 0.385 | 0.39 | 0.377 | 0.38 | 0.371 | 0.375 | 0.373 |
| Random edge | Graph head | other stop | 49.4 | 0.489 | 0.327 | 0.557 | 0.545 | 0.547 | 0.545 | 0.553 |
| Random edge | Graph head | scientific stop, investigational | 10.4 | 0.456 | 0.271 | 0.277 | 0.387 | 0.27 | 0.408 | 0.563 |
| Random edge | Graph head | scientific stop, other role | 13.2 | 0.422 | 0.281 | 0.251 | 0.324 | 0.571 | 0.572 | 0.579 |
| Random edge | Label-only MF | no stop | 78.0 | 0.368 | 0.361 | 0.358 | 0.354 | 0.367 | 0.362 | 0.361 |
| Random edge | Label-only MF | other stop | 49.4 | 0.41 | 0.182 | 0.452 | 0.447 | 0.453 | 0.458 | 0.456 |
| Random edge | Label-only MF | scientific stop, investigational | 10.4 | 0.406 | 0.154 | 0.143 | 0.268 | 0.138 | 0.301 | 0.441 |
| Random edge | Label-only MF | scientific stop, other role | 13.2 | 0.402 | 0.244 | 0.227 | 0.329 | 0.454 | 0.466 | 0.454 |

Table S18b. Hetionet: contrasts involving the role-restricted policies and between scorers, with 95% paired disease-cluster bootstrap intervals (1,000 draws).

| Set | Contrast | Pooled AP Δ (95% CI) | Per-disease AP Δ (95% CI) |
|---|---|---|---|
| E1 random edge | typed_role - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | typed_role - no_writeback [mf] | +0.001 (−0.016, +0.021) | +0.004 (−0.011, +0.018) |
| E1 random edge | typed_role - no_writeback [graph] | −0.010 (−0.036, +0.025) | +0.005 (−0.019, +0.026) |
| E1 random edge | typed_role - no_writeback [hybrid] | −0.009 (−0.033, +0.020) | +0.005 (−0.020, +0.028) |
| E1 random edge | typed_role_scoped - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | typed_role_scoped - no_writeback [mf] | +0.006 (−0.006, +0.018) | +0.012 (+0.002, +0.021) |
| E1 random edge | typed_role_scoped - no_writeback [graph] | −0.004 (−0.031, +0.032) | +0.011 (−0.011, +0.032) |
| E1 random edge | typed_role_scoped - no_writeback [hybrid] | −0.001 (−0.025, +0.029) | +0.010 (−0.014, +0.033) |
| E1 random edge | typed_role - typed [degree] | 0.000 (−0.001, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | typed_role - typed [mf] | +0.017 (+0.003, +0.034) | +0.022 (+0.012, +0.034) |
| E1 random edge | typed_role - typed [graph] | +0.019 (+0.005, +0.035) | +0.023 (+0.007, +0.042) |
| E1 random edge | typed_role - typed [hybrid] | +0.022 (+0.009, +0.036) | +0.024 (+0.008, +0.042) |
| E1 random edge | typed_role - mask_all [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | typed_role - mask_all [mf] | −0.010 (−0.026, +0.007) | −0.016 (−0.027, −0.006) |
| E1 random edge | typed_role - mask_all [graph] | −0.017 (−0.031, −0.003) | −0.019 (−0.035, −0.006) |
| E1 random edge | typed_role - mask_all [hybrid] | −0.021 (−0.039, −0.003) | −0.018 (−0.032, −0.004) |
| E1 random edge | typed_role_scoped - typed_scoped [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | typed_role_scoped - typed_scoped [mf] | +0.014 (+0.002, +0.031) | +0.023 (+0.011, +0.036) |
| E1 random edge | typed_role_scoped - typed_scoped [graph] | +0.018 (+0.005, +0.032) | +0.014 (−0.004, +0.033) |
| E1 random edge | typed_role_scoped - typed_scoped [hybrid] | +0.018 (+0.005, +0.030) | +0.013 (−0.001, +0.030) |
| E1 random edge | typed_role_scoped - mask_all [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 random edge | typed_role_scoped - mask_all [mf] | −0.006 (−0.016, +0.002) | −0.008 (−0.014, −0.003) |
| E1 random edge | typed_role_scoped - mask_all [graph] | −0.011 (−0.022, −0.001) | −0.013 (−0.026, −0.002) |
| E1 random edge | typed_role_scoped - mask_all [hybrid] | −0.013 (−0.024, −0.003) | −0.013 (−0.027, −0.001) |
| E1 random edge | graph - mf [no_writeback] | +0.079 (−0.016, +0.180) | +0.019 (−0.047, +0.089) |
| E1 random edge | mf - degree [no_writeback] | +0.135 (+0.094, +0.186) | +0.185 (+0.137, +0.237) |
| E1 random edge | graph - degree [no_writeback] | +0.214 (+0.114, +0.313) | +0.204 (+0.150, +0.261) |
| E1 compound disjoint | typed_role - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role - no_writeback [mf] | −0.019 (−0.034, −0.007) | −0.016 (−0.023, −0.009) |
| E1 compound disjoint | typed_role - no_writeback [graph] | −0.050 (−0.084, −0.021) | −0.028 (−0.044, −0.013) |
| E1 compound disjoint | typed_role - no_writeback [hybrid] | −0.048 (−0.083, −0.021) | −0.025 (−0.040, −0.011) |
| E1 compound disjoint | typed_role_scoped - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role_scoped - no_writeback [mf] | −0.011 (−0.018, −0.006) | −0.006 (−0.012, −0.001) |
| E1 compound disjoint | typed_role_scoped - no_writeback [graph] | −0.017 (−0.036, +0.002) | −0.010 (−0.023, +0.002) |
| E1 compound disjoint | typed_role_scoped - no_writeback [hybrid] | −0.017 (−0.036, +0.002) | −0.010 (−0.023, +0.002) |
| E1 compound disjoint | typed_role - typed [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role - typed [mf] | +0.008 (+0.003, +0.013) | +0.002 (0.000, +0.006) |
| E1 compound disjoint | typed_role - typed [graph] | +0.033 (+0.017, +0.049) | +0.030 (+0.011, +0.050) |
| E1 compound disjoint | typed_role - typed [hybrid] | +0.034 (+0.020, +0.049) | +0.036 (+0.016, +0.056) |
| E1 compound disjoint | typed_role - mask_all [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role - mask_all [mf] | −0.019 (−0.033, −0.007) | −0.016 (−0.024, −0.009) |
| E1 compound disjoint | typed_role - mask_all [graph] | −0.050 (−0.090, −0.017) | −0.033 (−0.049, −0.018) |
| E1 compound disjoint | typed_role - mask_all [hybrid] | −0.048 (−0.090, −0.016) | −0.029 (−0.045, −0.015) |
| E1 compound disjoint | typed_role_scoped - typed_scoped [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role_scoped - typed_scoped [mf] | +0.006 (+0.004, +0.008) | +0.005 (0.000, +0.009) |
| E1 compound disjoint | typed_role_scoped - typed_scoped [graph] | +0.020 (+0.009, +0.031) | +0.018 (+0.009, +0.027) |
| E1 compound disjoint | typed_role_scoped - typed_scoped [hybrid] | +0.020 (+0.009, +0.031) | +0.018 (+0.009, +0.027) |
| E1 compound disjoint | typed_role_scoped - mask_all [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E1 compound disjoint | typed_role_scoped - mask_all [mf] | −0.011 (−0.018, −0.005) | −0.007 (−0.012, −0.002) |
| E1 compound disjoint | typed_role_scoped - mask_all [graph] | −0.016 (−0.028, −0.004) | −0.015 (−0.026, −0.005) |
| E1 compound disjoint | typed_role_scoped - mask_all [hybrid] | −0.016 (−0.028, −0.004) | −0.015 (−0.025, −0.005) |
| E1 compound disjoint | graph - mf [no_writeback] | +0.201 (+0.147, +0.262) | +0.276 (+0.235, +0.322) |
| E1 compound disjoint | mf - degree [no_writeback] | −0.005 (−0.007, +0.013) | +0.050 (+0.042, +0.058) |
| E1 compound disjoint | graph - degree [no_writeback] | +0.195 (+0.150, +0.258) | +0.326 (+0.285, +0.374) |
| E2 external approvals | typed_role - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | typed_role - no_writeback [mf] | 0.000 (−0.001, +0.002) | +0.002 (−0.001, +0.005) |
| E2 external approvals | typed_role - no_writeback [graph] | +0.005 (−0.002, +0.014) | +0.006 (−0.008, +0.020) |
| E2 external approvals | typed_role - no_writeback [hybrid] | +0.003 (−0.002, +0.011) | +0.010 (+0.001, +0.022) |
| E2 external approvals | typed_role_scoped - no_writeback [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | typed_role_scoped - no_writeback [mf] | +0.001 (−0.001, +0.003) | +0.005 (+0.001, +0.008) |
| E2 external approvals | typed_role_scoped - no_writeback [graph] | +0.007 (−0.002, +0.017) | +0.008 (−0.006, +0.022) |
| E2 external approvals | typed_role_scoped - no_writeback [hybrid] | +0.006 (−0.002, +0.016) | +0.012 (+0.002, +0.025) |
| E2 external approvals | typed_role - typed [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | typed_role - typed [mf] | 0.000 (−0.002, +0.001) | +0.001 (−0.001, +0.002) |
| E2 external approvals | typed_role - typed [graph] | +0.003 (−0.001, +0.010) | +0.007 (−0.001, +0.023) |
| E2 external approvals | typed_role - typed [hybrid] | +0.003 (−0.001, +0.009) | +0.006 (−0.001, +0.022) |
| E2 external approvals | typed_role - mask_all [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | typed_role - mask_all [mf] | −0.001 (−0.002, +0.001) | −0.002 (−0.004, 0.000) |
| E2 external approvals | typed_role - mask_all [graph] | −0.001 (−0.009, +0.004) | −0.001 (−0.004, +0.001) |
| E2 external approvals | typed_role - mask_all [hybrid] | −0.002 (−0.009, +0.004) | −0.002 (−0.006, +0.002) |
| E2 external approvals | typed_role_scoped - typed_scoped [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | typed_role_scoped - typed_scoped [mf] | 0.000 (−0.001, +0.001) | 0.000 (−0.001, +0.002) |
| E2 external approvals | typed_role_scoped - typed_scoped [graph] | +0.002 (−0.001, +0.007) | +0.006 (−0.002, +0.022) |
| E2 external approvals | typed_role_scoped - typed_scoped [hybrid] | +0.002 (0.000, +0.007) | +0.006 (−0.001, +0.022) |
| E2 external approvals | typed_role_scoped - mask_all [degree] | 0.000 (0.000, 0.000) | 0.000 (0.000, 0.000) |
| E2 external approvals | typed_role_scoped - mask_all [mf] | 0.000 (−0.001, 0.000) | +0.001 (−0.001, +0.003) |
| E2 external approvals | typed_role_scoped - mask_all [graph] | +0.001 (0.000, +0.001) | +0.001 (0.000, +0.002) |
| E2 external approvals | typed_role_scoped - mask_all [hybrid] | +0.001 (0.000, +0.002) | +0.001 (0.000, +0.002) |
| E2 external approvals | graph - mf [no_writeback] | +0.011 (0.000, +0.024) | +0.013 (−0.012, +0.037) |
| E2 external approvals | mf - degree [no_writeback] | +0.008 (+0.002, +0.015) | +0.017 (0.000, +0.035) |
| E2 external approvals | graph - degree [no_writeback] | +0.019 (+0.006, +0.034) | +0.029 (+0.005, +0.054) |

## S19. Tables moved from the main text

Table S19. Lexical stop-reason rules of software release 0.2.0 compared with expert labels, and their agreement with the Open Targets classifier on Hetionet-mapped trials.

| Rule | Expert-labelled positives | Rule positives | Precision | Recall | F1 | Kappa vs classifier (Hetionet trials) |
|---|---|---|---|---|---|---|
| Efficacy | 368 | 100 | 0.89 | 0.24 | 0.38 | 0.48 |
| Safety | 211 | 277 | 0.45 | 0.59 | 0.51 | 0.55 |
| Operational | 2,373 | 1,791 | 0.90 | 0.68 | 0.77 | 0.50 |

Expert categories: efficacy, Negative; safety, Safety_Sideeffects; operational, Business_Administrative, Insufficient_Enrollment, Logistics_Resources, Study_Staff_Moved or Covid19. Kappa was computed on 5,810 stopped trials with a stated reason.

Table S20. Worked examples and the actions assigned by the write-back rule and the checker.

| Example | Evidence encoded | Failure annotation | Rule implication | Checker status |
|---|---|---|---|---|
| Baricitinib / COVID-19 | Graph-derived hypothesis; ACTT-2, COV-BARRIER and RECOVERY; FDA approval (2022) | None | — | Incomplete (exposure and safety domains undocumented) |
| Pimozide / ALS | NCT03272503, status unknown, no results | None | — | Incomplete |
| Evacetrapib / high-risk vascular disease | NCT01687998 terminated for insufficient efficacy; ACCELERATE; REVEAL (anacetrapib) | Stopped trial, efficacy, same concept | Negate (scoped) | Scoped counterevidence |
| Plazomicin / complicated UTI | EPIC; FDA approval (2018); sponsor's Chapter 11 filing (2019) | Corporate event, operational | Defer | Incomplete |

UTI, urinary tract infection.
