# KG Repurposing Audit: failure write-back (v0.5.1)

Code, schema and results for **Stopped clinical trials mark tested hypotheses: implications for negative labels in knowledge-graph drug repurposing**.

The repository asks what a terminated, withdrawn or suspended clinical trial tells a drug-repurposing knowledge graph, and how such trials should change its training labels. It maps Open Targets 26.09 clinical reports onto Hetionet v1.0 and PrimeKG v2, types every stopped trial by stop reason, scope and the role of each drug in the trial arms (ClinicalTrials.gov API v2), tests whether the overlap between stopped pairs and indications reflects testing intensity (logistic models, a degree-preserving permutation null and a review of 80 indications with failed trials), compares write-back policies, placebo negatives and a testing-intensity rule with several scorers on full candidate grids, asks whether scorers rank later-tested pairs as high as later-approved pairs, and provides a LinkML evidence-handoff schema (Biolink 4.4.5 mappings) that records failure type, scope match and write-back implication.

## Contents

| Path | What it is |
|---|---|
| `src/kg_audit/failures.py` | stop-reason categories → failure type; scope match; write-back rule (retain / negate / qualify / defer) |
| `src/kg_audit/writeback.py` | policies (including role-restricted negation), per-pair training labels and weights, scorers (degree, disease degree, MF, graph head, hybrid), metrics |
| `src/kg_audit/roles.py` | role of a drug in a trial (investigational, comparator, background therapy, ...) from ClinicalTrials.gov arm groups |
| `src/kg_audit/placebo.py` | degree-matched (double-edge swap) and uniform placebo negative sets |
| `src/kg_audit/fastboot.py` | exact disease-cluster bootstrap of pooled and per-disease AP without stored score vectors (used for PrimeKG) |
| `src/kg_audit/evidence.py` | record validation and dated documentation checker |
| `src/kg_audit/bridge.py` | exports ranked candidates with Open Targets evidence attached |
| `schemas/handoff.linkml.yaml` | LinkML schema v0.3.0; `handoff.schema.json` is generated from it with LinkML 1.11.1 |
| `examples/` | baricitinib, pimozide, evacetrapib and plazomicin encoded as records, with checker output |
| `pipelines/` | numbered scripts that download the inputs and rebuild every result, table and figure |
| `paper_results/` | aggregate results and candidate-level scores used in the manuscript; v0.4.0 adds `role/` (drug roles, role-restricted policies, the 50-assignment check), `primekg/` and `partitions_primekg/` (replication), `testedness/` and `scores/role_policies_*` (written by `pipelines/04d_export_release_files.py`); v0.5.0 adds `selection/` (testing-intensity models, permutation null, registry-history scorers), `case_review/` (80 dossiers, two independent codings and the adjudicated codes), `extra_hetionet/` and `extra_primekg/` (placebo negatives and testing-intensity policies) and `tables_v5/` |
| `scripts/audit_v02/` | scripts that reconstruct and audit the v0.2.0 evaluation (Supplementary S8) |
| `legacy/v0.2/` | the v0.2.0 benchmark and semi-synthetic code, kept for traceability; not used by v0.3.0 |

## Quick start

```bash
python -m venv .venv && source .venv/bin/activate
python -m pip install -e '.[test]'
python -m pytest            # write-back tests are skipped without torch
python -m pip install -e '.[benchmark]'
python -m pytest            # full suite
```

## Reproduce the manuscript

All inputs are public; nothing patient-level is used.

```bash
IN=$PWD/work/data; W=$PWD/work
pipelines/00_download_data.sh $IN                          # Hetionet v1.0, Open Targets 26.09, expert stop reasons
python pipelines/01_map_open_targets.py $IN/ot $IN/hetionet $W/mapped
python pipelines/02_make_partitions.py $IN/hetionet $W/partitions
MODE=pretrain MODEL=distmult NEG=typed SEED=1 HET=$IN/hetionet OUT=$W/kge python pipelines/03_train_kge.py
export HET=$IN/hetionet DATA=$W/mapped PART=$W/partitions EMB=$W/kge OUT=$W/results
python pipelines/04_writeback_experiment.py select
python pipelines/04_writeback_experiment.py e1 random
python pipelines/04_writeback_experiment.py e1 compound
python pipelines/04_writeback_experiment.py e2
for x in random compound e2; do python pipelines/04_writeback_experiment.py boot $x; done
OTDIR=$IN/ot GOLD=$IN/gold/data.json python pipelines/05_descriptive.py
python pipelines/07_handoff_export.py
python pipelines/08_make_figures.py $W/results $W/figures
python pipelines/09_make_tables.py $W/results $W/tables
```

Drug roles, role-restricted policies and the tested-versus-approved analysis (v0.4.0):

```bash
python pipelines/01b_fetch_trial_arms.py                  # ClinicalTrials.gov API v2 arms of all stopped trials -> work/data/ctgov
python pipelines/01c_assign_roles.py $W/mapped $IN/ot work/data/ctgov hetionet $IN/hetionet/hetionet-v1.0-nodes.tsv
OUT=work/results_role_only POLICY_SET=role_only python pipelines/04_writeback_experiment.py e1 random   # likewise e1 compound, e2
python pipelines/04c_merge_role_runs.py                    # merge with the main run; then boot with OUT=work/results_role POLICY_SET=role
python pipelines/05b_role_descriptive.py; python pipelines/05c_enrichment_hetionet.py; python pipelines/05d_negation_purity.py
python pipelines/16_tested_vs_approved.py hetionet
```

Tests of the selection explanation, placebo negatives and the testing-intensity rule (v0.5.0):

```bash
python pipelines/19_selection_analyses.py hetionet          # OUT=work/results_selection; likewise with DATA=work/mapped_primekg ... primekg
python pipelines/17_placebo_intensity_hetionet.py e1 random  # likewise e1 compound, e2, then boot random|compound|e2 (BASE=work/results_role)
python pipelines/18_placebo_intensity_primekg.py e1 random   # likewise e1 compound, e2, contrasts (BASE=work/results_primekg)
python pipelines/20_case_review_dossiers.py hetionet && python pipelines/20_case_review_dossiers.py primekg
python pipelines/20_case_review_dossiers.py fetch && python pipelines/20_case_review_dossiers.py dossier   # codes: paper_results/case_review/case_taxonomy.tsv
python pipelines/21_ranks_and_reference_scorers.py hetionet  # likewise primekg
python pipelines/09c_make_tables_v5.py; python pipelines/09d_make_supp_selection.py; python pipelines/09e_make_supp_S22.py
python pipelines/08_make_figures.py                         # Figures 1-7 and S1-S3
```

PrimeKG replication:

```bash
pipelines/10_download_primekg.sh                           # nodes.csv and edges (stored as compact arrays)
python pipelines/11_map_open_targets_primekg.py work/data/ot work/data/primekg work/mapped_primekg
python pipelines/01c_assign_roles.py work/mapped_primekg work/data/ot work/data/ctgov primekg work/mapped_primekg/graph.json
python pipelines/12_make_partitions_primekg.py && python pipelines/13_train_kge_primekg.py
for s in select "e1 random" "e1 compound" e2 contrasts; do python pipelines/14_writeback_primekg.py $s; done
python pipelines/15_descriptive_primekg.py; python pipelines/16_tested_vs_approved.py primekg
```

The end-to-end embedding ablation (Supplementary S7) uses `MODE=e2e PARTITION=$W/partitions/random_p42.json NEG=typed|uniform` with `03_train_kge.py`, summarised by `06_e2e_summary.py`. Every step skips work whose output already exists. On a shared 4-core CPU the Hetionet embedding pretraining took about 30 minutes, the Hetionet write-back experiment about 3 hours, the PrimeKG embeddings about 30 minutes for 10 epochs and the PrimeKG write-back experiment about 1 hour. The PrimeKG run used for the manuscript was configured for 15 epochs and stopped by the operating system after epoch 13; the analysis uses the 5- and 10-epoch checkpoints (`paper_results/primekg/kge_run_note.json`), which `13_train_kge_primekg.py` reproduces with its default of 10 epochs.

## Data and licences

| Resource | Licence | Use |
|---|---|---|
| [Hetionet v1.0](https://github.com/hetio/hetionet) | CC0 | graph and recorded treatments |
| [Open Targets Platform 26.09](https://platform.opentargets.org/downloads) | CC0 | clinical reports, stop-reason categories, drug and disease cross-references |
| [Expert-labelled stop reasons](https://huggingface.co/datasets/opentargets/clinical_trial_reason_to_stop) | Apache-2.0 | evaluation of lexical rules |
| [PrimeKG v2](https://doi.org/10.7910/DVN/IXA7BM) | CC0 1.0 | replication graph |
| [ClinicalTrials.gov API v2](https://clinicaltrials.gov/data-api/api) | public domain (U.S. NLM terms) | arm groups and interventions of stopped trials |

Input checksums of the files used for the manuscript are listed in `paper_results/input_sha256.txt`.

## What the code does not do

A score never fills a clinical field, and a registry status alone never creates a negative label. The write-back rule negates a drug–disease label only for an efficacy or safety stop whose trial condition is the same concept as the graph disease. The analyses in the manuscript show that even such stops, and even stops of the investigational drug, often concern true indications; the records are meant for expert review, not for treatment decisions or automatic negative labels.

## Licence and citation

MIT licence. Please cite the article and the archived release (see `CITATION.cff`).
