# Data access

All inputs of release 0.3.0 are public and can be downloaded with `pipelines/00_download_data.sh DATA_DIR`. Keep data outside the repository; the release tree check (`scripts/check_public_tree.py`) rejects databases and large binary inputs.

## Hetionet v1.0

Nodes (`hetionet-v1.0-nodes.tsv`) and edges (`hetionet-v1.0-edges.sif.gz`) from <https://github.com/hetio/hetionet> (CC0). Hetionet integrates data available in 2016; its Compound–treats–Disease (CtD) edges are the recorded treatments and its Compound–palliates–Disease (CpD) edges are excluded from training losses and evaluation.

## Open Targets Platform 26.09

The `clinical_report`, `clinical_indication`, `clinical_target`, `drug_molecule` and `disease` parquet files from <https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output/> (CC0). Clinical reports integrate ClinicalTrials.gov (via AACT), drug labels, regulatory agencies and curated resources; stopped trials carry the stop-reason categories assigned by the Open Targets classifier. A different platform release is a different input; record its checksums.

## Expert-labelled stop reasons

`data.json` from <https://huggingface.co/datasets/opentargets/clinical_trial_reason_to_stop> (Apache-2.0, doi:10.57967/hf/2600): 3,747 stop-reason texts with expert categories, used only to evaluate the lexical rules of release 0.2.0.

## Not used

Release 0.3.0 uses no patient-level data. The MIMIC-IV illustration and its code were removed; PrimeKG and DrugCentral are not needed.
