# Optional observational illustration

This component requires a user's own authorised MIMIC-IV 3.1 files. No patient data are supplied. It is a descriptive sensitivity illustration for evidence handling, not a clinical validation of knowledge-graph predictions or a recommendation to use statins for sepsis.

Install `python -m pip install -e '.[clinical]'` and place the resource's unmodified `hosp/` and `icu/` folders under your work directory's `data/mimic-iv-3.1/` folder. Then run:

```bash
python scripts/run_pipeline.py mimic-extract --workdir /path/to/work
python scripts/run_pipeline.py mimic-severity --workdir /path/to/work
python scripts/run_pipeline.py mimic-assemble --workdir /path/to/work
python -m kg_audit.observational \
  --cohort /path/to/work/outputs/mimic/analytic_full.csv \
  --database /path/to/work/outputs/mimic/mimic.duckdb \
  --output /path/to/work/outputs/mimic/associations.json --bootstrap 1000
```

The extraction stages scan large local tables and can take substantial time and disk space. They retain historical internal column names for compatibility. `prevalent_user` means an observed earlier inpatient order in the same admission, not a documented pre-admission medication washout. `sepsis3` is a simplified operational infection/SOFA flag, not independent validation of a consensus clinical diagnosis.

The analysis recalculates order assignment using exact 48-hour comparisons, including later initiators in the no-order-by-48-hours comparator. It requires hospitalisation and recorded survival at the landmark. Recorded hospital death times take precedence; date-only deaths crossing an eligibility or outcome boundary are flagged and excluded rather than assigned invented hours. The outcome is recorded death through day 28 **from ICU entry**, equivalent to at most 26 days after the landmark.

The primary descriptive rerun does not use the optional code-status exclusion. It reports demographic and extended-covariate overlap-weighted associations, with propensity fitting and imputation repeated in each patient bootstrap. Extended covariates may follow treatment and hospital diagnosis codes lack reliable onset times. Neither measured balance nor a narrow interval establishes a causal effect. The module exports aggregate results only.
