# KG Repurposing Audit — failure write-back (v0.3.0)

Code, schema and results for **Stopped clinical trials as negative labels for knowledge-graph drug repurposing: an evaluation of write-back policies**.

The repository asks how terminated, withdrawn or suspended clinical trials should change the training labels of a drug-repurposing knowledge graph. It maps Open Targets 26.09 clinical reports onto Hetionet v1.0, types every stopped trial by stop reason and scope, compares five write-back policies with four scorers on full candidate grids, and provides a LinkML evidence-handoff schema (Biolink 4.4.5 mappings) that records failure type, scope match and write-back implication.

## Contents

| Path | What it is |
|---|---|
| `src/kg_audit/failures.py` | stop-reason categories → failure type; scope match; write-back rule (retain / negate / qualify / defer) |
| `src/kg_audit/writeback.py` | policies, per-pair training labels and weights, scorers (degree, disease degree, MF, graph head, hybrid), metrics |
| `src/kg_audit/evidence.py` | record validation and dated documentation checker |
| `src/kg_audit/bridge.py` | exports ranked candidates with Open Targets evidence attached |
| `schemas/handoff.linkml.yaml` | LinkML schema v0.3.0; `handoff.schema.json` is generated from it with LinkML 1.11.1 |
| `examples/` | baricitinib, pimozide, evacetrapib and plazomicin encoded as records, with checker output |
| `pipelines/` | numbered scripts that download the inputs and rebuild every result, table and figure |
| `paper_results/` | aggregate results and candidate-level scores used in the manuscript |
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

The end-to-end embedding ablation (Supplementary S7) uses `MODE=e2e PARTITION=$W/partitions/random_p42.json NEG=typed|uniform` with `03_train_kge.py`, summarised by `06_e2e_summary.py`. Every step skips work whose output already exists. On a 4-core CPU the embedding pretraining took about 30 minutes and the write-back experiment about 3 hours.

## Data and licences

| Resource | Licence | Use |
|---|---|---|
| [Hetionet v1.0](https://github.com/hetio/hetionet) | CC0 | graph and recorded treatments |
| [Open Targets Platform 26.09](https://platform.opentargets.org/downloads) | CC0 | clinical reports, stop-reason categories, drug and disease cross-references |
| [Expert-labelled stop reasons](https://huggingface.co/datasets/opentargets/clinical_trial_reason_to_stop) | Apache-2.0 | evaluation of lexical rules |

Input checksums of the files used for the manuscript are listed in `paper_results/input_sha256.txt`.

## What the code does not do

A score never fills a clinical field, and a registry status alone never creates a negative label. The write-back rule negates a drug–disease label only for an efficacy or safety stop whose trial condition is the same concept as the graph disease; such identifier-level scope matching is necessary but not sufficient, and the records are meant for expert review, not for treatment decisions.

## Licence and citation

MIT licence. Please cite the article and the archived release (see `CITATION.cff`).
