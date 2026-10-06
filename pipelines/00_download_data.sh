#!/usr/bin/env bash
# Download the public inputs of the v0.3.0 analysis into a caller-owned data directory.
#   Hetionet v1.0 (CC0), Open Targets Platform 26.09 (CC0), expert-labelled stop reasons (Apache-2.0)
# Usage: pipelines/00_download_data.sh /path/to/data
set -euo pipefail
DATA="${1:?usage: $0 DATA_DIR}"
mkdir -p "$DATA/hetionet" "$DATA/ot" "$DATA/gold"
HET=https://github.com/hetio/hetionet/raw/main/hetnet/tsv
[ -s "$DATA/hetionet/hetionet-v1.0-nodes.tsv" ] || curl -fsSL -o "$DATA/hetionet/hetionet-v1.0-nodes.tsv" "$HET/hetionet-v1.0-nodes.tsv"
[ -s "$DATA/hetionet/hetionet-v1.0-edges.sif.gz" ] || curl -fsSL -o "$DATA/hetionet/hetionet-v1.0-edges.sif.gz" "$HET/hetionet-v1.0-edges.sif.gz"
OT=https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output
for d in clinical_report clinical_indication clinical_target drug_molecule disease; do
  [ -s "$DATA/ot/$d.parquet" ] || curl -fsSL -o "$DATA/ot/$d.parquet" "$OT/$d/00000000.parquet"
done
HF=https://huggingface.co/datasets/opentargets/clinical_trial_reason_to_stop/resolve/main
[ -s "$DATA/gold/data.json" ] || curl -fsSL -o "$DATA/gold/data.json" "$HF/data.json"
( cd "$DATA" && sha256sum hetionet/* ot/* gold/* > SHA256SUMS.txt )
echo "inputs in $DATA; checksums in $DATA/SHA256SUMS.txt"
