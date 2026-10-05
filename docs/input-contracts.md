# Input and output contracts

## Hetionet

`hetionet-v1.0-nodes.tsv` has `id`, `name`, and `kind` columns. The gzipped edge SIF is tab-separated with a header and three columns in source, metaedge, target order. The target relation is `CtD`. Rename the upstream edge file to the documented local filename without editing its contents.

`drugcentral_compound_annot.tsv` must include string columns `drugbank_id` and `smiles`. The `drugbank_id` is the unprefixed identifier used by the mapped Hetionet Compound node. Missing or invalid structures remain unmapped; they are not imputed from drug names. The scaffold grouping is a treatment-label partition and does not remove auxiliary node identities.

Preparation writes the graph-specific `splits.json`. Its keys are `random`, `coldstart`, and `scaffold`; each contains `train_pos`, `train_neg`, `test_pos`, `neg_random`, and `neg_degmatch`, as arrays of compound/disease identifier pairs. The historical key `neg_degmatch` means **disease-frequency sampling**, not joint training-degree-bin matching.

The `scores_{model}_{split}_seed{seed}.tsv` files contain `group`, `compound`, `disease`, and uncalibrated `score`. The summary pipeline exports common candidate manifests and seed-specific fitting/validation records. It will fail on missing model files rather than fabricate a cell.

## PrimeKG

`kg.csv` must contain `relation`, `x_type`, `x_id`, `y_type`, and `y_id`. The pipeline retains drug-protein, disease-protein and protein-protein features. Indication pairs are the labelled target. Contraindication and off-label pairs are excluded from the unlabelled candidate pools.

Within each seed/task, the pipeline holds out task-consistent validation labels before constructing fitting-degree references. Uniform and joint-bin candidate sets share test positives. Joint bins are zero, one, two-to-three, four-to-seven, and successive powers of two. Sampling is with replacement, and the output records both draw counts and unique pairs. Empty support is an explicit error.

Outputs contain fitting positives, separate baseline-fitting unlabelled pairs, validation pairs, test candidate rows, per-model scores, checkpoint state dictionaries, validation histories and a recipe fingerprint. Model scores are not calibrated clinical probabilities.

## User-supplied comparisons

A CSV or TSV has `compound`, `disease`, `label` (0 or 1), and numerical score columns. Supply one model and one reference column for each seed, in matching order. Candidate identities must refer to the same task and sampling design. Repeated identical pairs can encode sampling multiplicity; conflicting labels for the same pair are rejected.

`label=0` means sampled unlabelled in the benchmark. It must not be relabelled as proven clinical inefficacy.

## Evidence records

The LinkML and JSON Schema files describe the serialised structure. Semantic rules are in `kg_audit.evidence`: required evidence domains, valid references, date availability, conflicting observations and contextual scope. Structural validity and readiness for documentation review are separate checks. Unknowns can be omitted; absence never creates a pass.

No trial registry status alone creates a negative drug-disease edge. Scoped counterevidence includes drug, indication, population, regimen, comparator, endpoint and source provenance. A changed comparator or regimen cannot silently supersede another treatment strategy.
