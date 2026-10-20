# Changelog

## 0.7.4 (2026-10-20)

- Documentation and design figure describe the literature screening, case coding and label check as automated procedures.

## 0.7.3 (2026-10-19)

- Design figure: the literature box describes the check of excluded articles as a second screening.

## 0.7.2 (2026-10-19)

- repoDB audit (pipeline 23): the DrugCentral endpoint `omop_relationship/relationship_name/indication` also returns contraindication rows (27,755 of 39,986). Only rows with relationship `indication` (12,231) are now used as approval evidence. Corrected repoDB results: same-UMLS-concept matches 18 of 3,260 failed pairs (0.6%); 1,497 mapped pairs, of which 148 (9.9%) are indications of the same concept and 237 (15.8%) of the same concept or a stage, line or subtype; popularity null (pipeline 25) 153 observed against 20.5 expected, ratio 7.4 (5.3 to 11.8). Two sampled stage-or-subtype matches of the label check (R39, R56) rested only on contraindication rows and were removed from `precision_check.tsv` (28 remain at that level).
- repoDB figure labels its percentages from counts.
- Test `tests/test_repodb_outputs.py` checks that the released repoDB outputs use indication rows only.

## 0.7.1 (2026-10-18)

- Case-review sensitivity (pipeline 26) counts registered trials by NCT identifier; some cases list several trials, so the union of contradictory or unclear codes covers 12 trials, not 9 distinct strings. Counts of the final and the union-of-yes codes (4 and 5 trials) are unchanged.
- Timing output adds per-trial shares (all or any of a trial's indication pairs started in or after the dated approval).
- `requirements-reproduction.txt` and `configs/numerical-environment.json` now list the versions of the environment that produced the reported numbers (statsmodels and scipy added).

## 0.7.0 (2026-10-18)

- Robustness checks (pipeline 26, Hetionet; `paper_results/robustness/`): start years of stopped phase 1-3 trials of indications relative to the drug's first approval (ChEMBL) and to the earliest dated European (EMA) or Japanese (PMDA) approval for the same disease; Mantel-Haenszel odds ratio of a scientific stop over six strata of trial counts; enrichment of external approved indications by source of approval evidence; contradiction counts of the case review under the union of all codings; summary of the second literature screening.
- Literature survey: blinded second screening of a random 40 of the 176 excluded articles (`literature/second_screen_sample.jsonl`, codebook `literature/second_screen_codebook.md`).
- New supplementary figure `FigS_temporal_order` (`pipelines/figures/supp/b_temporal_order.py`); study-design figure text updated.
- Tests for the consistency of the robustness outputs (`tests/test_robustness_outputs.py`).
- No earlier analysis or number changed.

## 0.6.0 (2026-10-17)

- repoDB audit (pipeline 23): stop reasons of the failed records (matched to the Open Targets classifier labels of the same text), failed pairs that are indications at three levels (same UMLS concept in DrugCentral; same Open Targets concept; stage, line or subtype of an indication), by testing intensity, phase and stop reason; concept-level edges for the popularity null; a label check of 60 sampled matches (`paper_results/repodb/`).
- Model-free checks for any source of negative labels (`kg_audit.label_audit`, pipeline 25): overlap with known positives, overlap by testing intensity and a degree-preserving popularity null; applied to repoDB.
- Literature survey (pipeline 24): Europe PMC search, passage extraction and the screening of 197 open-access articles for stopped-trial negatives (`paper_results/literature/`).
- Case review: a third, automated source-verified coding of the 80 cases with retrieval of registry records, posted results, publications and labels, with cited sources; adjudication of the 14 disagreements; final codes in `case_review/case_taxonomy_final_v2.tsv` (6 contradicting cases from 4 trials, 1 unclear).
- New figures: study design (`fig1_design.py`), repoDB audit (`fig_repodb.py`); case-review figures updated; `fig2_stopped.py` takes panel letters and an output name.
- No earlier analysis or number changed apart from the case-review codes.

## 0.5.2 (2026-10-07)

- Figures for the final manuscript: Figure 5 (case review, `fig_casereview.py`), Figure 9 (where the loss from negative labels falls, `fig_lossmap.py`, with distribution data from pipeline 22) and supplementary Figures S1, S2, S4-S9 and S12-S16 (`pipelines/figures/supp/`), which plot the supplementary tables instead of listing them.
- The generated supplementary tables used by these scripts are kept in `paper_results/supp_tables_md/`; `paper_results/selection/rank_strata_*.json` holds the percentile ranks behind Figure 9D.
- No analysis or reported number changed.

## 0.5.1 (2026-10-07)

- Case review: independent second coding of the 80 dossiers and adjudication against registry records, posted results and publications (`paper_results/case_review/case_taxonomy_verified.tsv`, with both codings); Section S21 generator reports agreement and changes.

## 0.5.0 (2026-10-07)

- Tests of the selection explanation (pipeline 19): logistic models of being an indication among tested pairs with and without adjustment for testing intensity (disease-cluster bootstrap), share of indications by number of trials, phase and disease area, degree-preserving permutation null for the enrichment of approved indications among stopped pairs, scorers built only from registry history, per-disease heterogeneity.
- Placebo negatives and a testing-intensity write-back policy (pipelines 17 and 18; `kg_audit.placebo`): degree-matched and uniform placebo sets of the same size as each policy's negated set; negation restricted to pairs with at most 1, 2 or 4 registered trials.
- Case review (pipeline 20): dossiers of 80 indications with the cleanest failure records from the ClinicalTrials.gov API v2, codebook and codes.
- Supporting numbers (pipeline 21), Tables 4-5 and Additional file sections S20-S22 (pipelines 09c-09e), seven main figures (pipeline 08 and `pipelines/figures/`).
- Tests for the placebo sets.

## 0.4.0 (2026-10-06)

- Drug roles in stopped trials from ClinicalTrials.gov API v2 arm groups (`kg_audit.roles`; pipelines 01b, 01c): investigational agent, comparator, background therapy, mixed, head-to-head, single-arm combination, described only.
- Two role-restricted write-back policies (`typed_role`, `typed_role_scoped`); negative-set size and purity per policy; enrichment of approved indications among stopped pairs by stop type and role.
- Tested-versus-approved analysis: AUROC between later-approved, later-failed, newly tested and untested pairs for each scorer (pipeline 16).
- PrimeKG v2 replication (pipelines 10-15): Open Targets mapping to MONDO diseases, partitions, DistMult embeddings without drug-disease edges, seven policies, exact inline disease-cluster bootstrap (`kg_audit.fastboot`).
- New figures (tested versus approved; policy effects in both graphs) and tables; tests for role rules, role policies and the inline bootstrap.

## 0.3.0 (2026-10-06)

- New question and analysis: write-back of real stopped-trial evidence (Open Targets 26.09) into Hetionet training labels, with five policies, four scorers and three outcome sets (held-out treatments, approved indications absent from the graph, later scientific failures).
- `kg_audit.failures`: stop-reason categories (Razuvayevskaya et al. 2024) to failure types; scope match from ontology depth; write-back rule (retain / negate / qualify / defer).
- `kg_audit.writeback` replaces the semi-synthetic experiment: per-pair labels and weights on the full compound–disease grid; degree, disease-degree, matrix-factorisation, graph-head and hybrid scorers.
- Evaluation: full candidate grids restricted to diseases with recorded treatments; type-constrained corruption; embeddings pretrained without any treatment edge; five partitions per task; disease-cluster bootstrap with per-disease AP.
- Schema 0.3.0: failure type, stop-reason categories, scope match, event type and Biolink 4.4.5 / PROV-O mappings; JSON Schema generated with LinkML 1.11.1; the checker enforces the write-back rule and normalises scope comparisons.
- `bridge.export_ranked` attaches dated approvals, trials and typed failure annotations to ranked candidates.
- Four clinical examples encoded as schema instances with checker output.
- Removed the MIMIC-IV observational illustration and all patient-data code. Moved the v0.2.0 Hetionet/PrimeKG benchmark, semi-synthetic write-back and registry convenience-sample code to `legacy/v0.2`.
- Added `scripts/audit_v02` (reconstruction of the v0.2.0 evaluation; Supplementary S8).

## 0.2.0

- Align model/reference candidate rows and fitting-graph degree features.
- Separate task definitions, disease-frequency sampling and joint training-degree-bin sampling.
- Add validation-selected PrimeKG checkpoints, run histories and recipe fingerprints.
- Add a separate uncertainty-mask arm to the four-policy write-back experiment.
- Require dated, linked evidence for all eight documentation domains; preserve contextual scope and reject future-evidence use.
- Add a public registry metadata audit without deriving efficacy labels from registry status.
- Add an aggregate-only descriptive observational analysis with explicit landmark and timing checks.
- Provide synthetic examples, input contracts, tests, CI and a tracked-tree release scan.

Database contents, raw candidate manifests and patient-level data are excluded from the public repository.
