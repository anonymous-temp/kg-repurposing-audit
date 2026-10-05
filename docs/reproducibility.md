# Reproduction levels and interpretation

1. **Software checks:** run the synthetic demo and unit tests without datasets. This verifies code behaviour, not biomedical accuracy.
2. **Metric recomputation:** use the exact candidate score manifests and compute AP, AUROC and paired intervals. This verifies calculations conditional on fitted scores.
3. **Fresh training:** obtain the input snapshots, run the data preparation and training stages, record model selection and regenerate predictions. This is distinct from loading a cache.
4. **External or prospective validation:** use independently sourced or time-valid evidence with overlap auditing. This has not been established by the static graph experiments or by the registry metadata audit.

The paper separates these levels. Hetionet reused prediction files in part; new compound-disjoint fits were added during the earlier revision. The new PrimeKG experiment selects checkpoints on validation AP with a maximum of 40 epochs, evaluations every two epochs after epoch eight, and three non-improving evaluations before stopping. It does not select on test AP and is not a comprehensive hyperparameter search.

Cluster-bootstrap intervals are exploratory and conditional on the supplied graph partition, candidates and fitted seeds. Seed SD is not a confidence interval. Disease and compound clustering address some shared-entity dependence but do not model the full curation process or incomplete therapeutic labels.

The local permutation reference is based on a bipartite edge-swap principle. It is not an exact reproduction of the xswap package, and mixing has not been established theoretically. Zero-degree held-out drugs remain zero degree under every degree-preserving permutation.

The four-policy experiment is semi-synthetic. Its annotation fractions are design parameters, not measured failure-type prevalence. A correctly supplied retain annotation can preserve information lost by blanket masking, but this result does not show that human failure labels will be accurate or that a schema improves clinical decisions.

The ClinicalTrials.gov sample is explicitly non-random: the first 100 returned records by ascending first-posted date within each discontinued-status stratum, first posted in 2015-2025 and updated by the specified cutoff. Descriptive lexical counts do not estimate population prevalence or diagnostic accuracy of a failure taxonomy.

The exact PrimeKG training source before whitespace-only formatting is preserved at commit `f835c95ecbbb116b507da9847cb1138310eb3172`. Its SHA-256 is recorded in `paper_results/primekg_protocol.json`. Formatting was checked for Python abstract-syntax-tree equivalence; scientific recipes were not changed.
