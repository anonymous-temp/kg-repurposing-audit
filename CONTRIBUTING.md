# Contributing

Use a small pull request with a clear scientific or software question. Add a regression test when correcting a comparison, sampling, temporal or evidence-scope error. Preserve existing provenance fields and explain any changed result.

Run `python -m pytest`, the fatal-error Ruff check, the synthetic demo and `python scripts/check_public_tree.py`. Biomedical claims require a primary source and explicit population, regimen, comparator and endpoint where relevant.

Use synthetic reproductions in issues. Do not attach source databases, licensed bulk exports, clinical credentials, patient records, patient-level derived data or fitted clinical models. A result on an internal dataset must not be presented as publicly reproducible without an accessible, permitted input path.

Please distinguish model fitting, cached scoring, metric recomputation and external validation. Do not remove inconvenient runs or describe a test-set-selected configuration as validation-selected.
