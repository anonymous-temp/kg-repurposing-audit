# Obtain data independently

This is a source-code repository. Do not upload database snapshots or patient-derived data in issues, pull requests, discussions, or release attachments.

## Hetionet

Obtain the official [Hetionet v1.0 TSV distribution](https://github.com/hetio/hetionet/tree/master/hetnet/tsv). Preserve the downloaded version and SHA-256 checksums. The benchmark target is the `CtD` treatment relation. Preserve source licences and attribution; individual upstream resources can impose additional conditions.

For scaffold grouping, obtain licensed DrugCentral structure annotations, map them to DrugBank identifiers, and export the two columns specified in `input-contracts.md`. Without those annotations, the recorded scaffold experiment cannot be reproduced faithfully. Do not substitute guessed structures or rename an alternative split as a scaffold split.

## PrimeKG

Use the [PrimeKG project](https://zitniklab.hms.harvard.edu/projects/PrimeKG/) and its [publication](https://doi.org/10.1038/s41597-023-01960-3) to obtain the intended release. Place the downloaded `kg.csv` in `data/primekg/` of your local work directory. The training recipe records its content hash. The same filename with changed contents is a different input snapshot.

## ClinicalTrials.gov

The public API does not require MIMIC credentials. The registry command saves raw responses locally with query URLs, retrieval times and hashes. It does not commit them. Registry records can change; reproducing a historical snapshot requires that exact saved response. A later live query is a refreshed audit, not an exact reproduction of the original record set.

## MIMIC-IV

MIMIC-IV is a credentialed resource. Follow the [PhysioNet access procedure and agreement](https://physionet.org/content/mimiciv/3.1/). Another person's access does not authorise you to use their downloaded data. Restrictions also matter for patient-level derivatives and trained clinical models.

The primary benchmark, registry audit and software tests do not require MIMIC. No claim in this repository establishes a causal statin effect. The manuscript's observational material is an illustration of evidence limitations, not a supported treatment recommendation.
