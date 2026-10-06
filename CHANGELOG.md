# Changelog

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
