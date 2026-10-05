# KG Repurposing Audit

Paired evaluation of drug-repurposing knowledge graphs and scoped, dated evidence records for downstream review.

The repository accompanies **Auditing knowledge-graph drug repurposing: task-matched evaluation and failure-typed evidence handoff**. It provides an executable research workflow, a database-free demonstration, unit tests, and data-access instructions. The code does not prescribe treatments or turn trial termination into an efficacy-negative label.

## What is included

- Task-specific evaluation with identical candidate rows for each model/reference comparison.
- Degree references calculated from the fitting graph, paired cluster intervals, and explicit candidate-rank domains.
- Hetionet and PrimeKG training/reanalysis pipelines. PrimeKG checkpoints are selected using task-consistent validation labels, not the test set.
- A four-policy, **semi-synthetic** write-back stress test, including a separate uncertainty-masking arm.
- A dated evidence-handoff implementation that preserves regimen, population, comparator and endpoint scope and does not fill missing evidence with a pass.
- A public ClinicalTrials.gov metadata audit with frozen local responses and request fingerprints.
- LinkML and generated JSON Schema specifications, synthetic examples, tests, and CI.

## Quick start without a database

Python 3.11 or 3.12 is recommended. The core package supports Python 3.9 and later; the recorded manuscript runtime is documented separately.

```bash
git clone https://github.com/anonymous-temp/kg-repurposing-audit.git
cd kg-repurposing-audit
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
kg-audit demo --output outputs/demo
python -m pytest
```

On Windows, activate with `.venv\Scripts\activate`. The demonstration creates explicitly synthetic scores and an incomplete evidence record. It does not reproduce the paper's real-data numbers. CI tests this database-free workflow.

## Data access and licences

**Databases are not bundled.** Keep your data in a separate work directory. Users are responsible for obtaining the source data and meeting each provider's access and licence conditions. Public availability is not a blanket redistribution licence.

| Resource | Access | Use here |
| --- | --- | --- |
| [Hetionet v1.0](https://github.com/hetio/hetionet) | Obtain upstream nodes and edges; retain source/licence information | Static treatment-link benchmark |
| [PrimeKG](https://www.nature.com/articles/s41597-023-01960-3) | Obtain the specified `kg.csv` release from the resource's official distribution | Indication benchmark |
| [DrugCentral](https://drugcentral.org/download) | Obtain and map structure annotations under the provider's terms | Hetionet scaffold annotation |
| [ClinicalTrials.gov API](https://clinicaltrials.gov/data-api/api) | Public API, queried explicitly by the user | Descriptive registry metadata audit |
| [MIMIC-IV 3.1](https://physionet.org/content/mimiciv/3.1/) | Credentialing, training and the applicable data-use agreement | Historical observational illustration only; not needed for primary benchmark or software tests |

No MIMIC records, patient-level derived data, clinical-model weights, credentials, private manuscript files, or database snapshots are included. The manuscript's MIMIC illustration is not a validated causal efficacy finding. There is no public patient-data reproduction shortcut in this repository.

Credentialed users can run the optional descriptive workflow in [docs/observational-illustration.md](docs/observational-illustration.md). It corrects order assignment at the landmark and preserves uncertainty in date-only death records; it does not make the dataset sufficient for causal identification.

Input preparation and exact file contracts are in [docs/data-access.md](docs/data-access.md) and [docs/input-contracts.md](docs/input-contracts.md).

## Reproduce the database-dependent experiments

Install the optional numerical dependencies:

```bash
python -m pip install -e '.[benchmark]'
```

Create a caller-owned work directory with this layout:

```text
work/
  data/
    hetionet/
      hetionet-v1.0-nodes.tsv
      hetionet-edges.sif.gz
      drugcentral_compound_annot.tsv
    primekg/
      kg.csv
  outputs/                       # generated locally, ignored by Git
```

Then run from the cloned repository:

```bash
python scripts/run_pipeline.py hetionet-prepare --workdir /path/to/work
python scripts/run_pipeline.py hetionet-kge --workdir /path/to/work --seeds 1 2 3 --max-epochs 20
python scripts/run_pipeline.py hetionet-rgcn --workdir /path/to/work --seeds 1 2 3 --max-epochs 30
python scripts/run_pipeline.py hetionet-summary --workdir /path/to/work

python scripts/run_pipeline.py primekg-train --workdir /path/to/work --seeds 1 2 3 --max-epochs 40
python scripts/run_pipeline.py primekg-summary --workdir /path/to/work

kg-audit-writeback --splits /path/to/work/outputs/hetionet/raw/splits.json \
  --output /path/to/work/outputs/writeback --seeds 1 2 3 4 5 6 7 8 9 10

kg-audit-registry --cutoff 2026-10-05 --per-status 100 \
  --output /path/to/work/outputs/registry
```

These training commands can take hours on CPU, especially the R-GCN pipeline. Runtime depends on hardware and stopping behaviour. The launcher sets deterministic seed controls. Bitwise identity across hardware/library versions is not guaranteed. A completed cache is not an independent retraining. Use a new work directory when changing inputs or recipes.

The manuscript distinguishes reused, reanalysed Hetionet predictions from new validation-selected PrimeKG fits. Running the entire Hetionet sequence above trains fresh models; it must not be described as proof that every historical checkpoint has been independently reproduced. See [docs/reproducibility.md](docs/reproducibility.md).

## Regenerate plots from published summaries

With the benchmark dependencies installed, run `python scripts/plot_aggregate_results.py --output outputs/reference_figures`. This reads the aggregate tables included in the repository and renders AP plots with seed-SD bars. It verifies the presentation of the reference numbers; reproducing those numbers requires the database-dependent pipelines. Aggregate cluster intervals and validation histories are also supplied under `paper_results/`.

## Use your own candidate scores

```bash
kg-audit compare --scores candidates.tsv \
  --model rotate_s1 rotate_s2 rotate_s3 \
  --reference degree_s1 degree_s2 degree_s3 \
  --clusters disease --bootstrap 5000 --output outputs/comparison.json

kg-audit evidence --record examples/incomplete_strategy.json \
  --as-of 2026-01-01 --output outputs/evidence_review.json
```

`compare` requires aligned `compound`, `disease`, `label`, and score columns. It reports the observed mean seed-specific AP difference and an exploratory paired cluster interval. `evidence` assesses documentation as of a date; even a complete record is ready only for evidence review, not for prescribing or deployment.

## Scientific interpretation

Task splits, candidate sampling and leakage are different concepts. Degree derived from the fitting graph is not automatically an illegitimate feature. Unlabelled pairs are not confirmed treatment failures. Matching bins does not match exact degree values. A sampled AP improvement does not establish a biological mechanism or prospective therapeutic benefit.

The registry audit is a convenience sample of interventional drug trials, not a representative estimate of trial-failure causes and not an external prediction benchmark. Its reason-text cues are lexical descriptors, not expert adjudications. The write-back experiment supplies hypothetical retain/defer annotations to recorded training positives; its benefit is conditional on that information being correct.

## Development and citation

```bash
python -m pytest --cov=kg_audit --cov-report=term-missing
python -m ruff check src tests scripts pipelines
python scripts/check_public_tree.py
```

Use [CITATION.cff](CITATION.cff) to cite the software, including the exact release and commit used. Method references and third-party acknowledgements are in [docs/references.md](docs/references.md). This software is MIT-licensed; source databases retain their own terms.

Funding: 2025 Hebei Provincial Major Science and Technology Support Plan, Innovative Application Scenario Project (252Q0103D).
