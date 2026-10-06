# Input and output contracts

## Pair table (`trial_pairs.tsv`, from `01_map_open_targets.py`)

One row per (clinical report, Hetionet compound, Hetionet disease): `report_id`, `compound` (DrugBank id), `disease` (DOID), `ot_disease` (condition term), `stage`, `origin`, `type`, `source`, `phase`, `status`, `year`, `stop_categories` ('|'-joined), `why_stopped`, `start_date`, `url`. `disease_map.tsv` gives the route (direct or ancestor) and the depth gap used for the scope match.

## Partitions (`02_make_partitions.py`)

`{task}_p{seed}.json` with `fit_pos`, `val_pos`, `test_pos`, `eval_grid` and `val_grid` as lists of [compound, disease]. The evaluation grid is every pair of a held-out compound with a disease that has at least one recorded treatment, minus fitting/validation positives and CpD pairs. A grid pair that is not a held-out treatment is unlabelled, not a verified non-indication.

## Candidate-level scores (`paper_results/scores/`)

Tab-separated, gzipped: `task`, `pseed`, `compound`, `disease`, `label`, then one column per policy and scorer (mean over initialisation seeds). Scores are uncalibrated and comparable only within a task, partition and scorer.

## Evidence records

The LinkML and generated JSON Schema files describe the serialised structure; semantic rules are in `kg_audit.evidence` and `kg_audit.failures`. A failure annotation's implication must follow the write-back rule: retain for a recorded treatment, negate only for a same-concept efficacy or safety failure, qualify for a scientific failure in another concept, defer otherwise. A registry status alone never creates a negative label, and evidence dated after the review date is never used.
