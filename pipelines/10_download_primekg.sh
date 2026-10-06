#!/bin/bash
set -euo pipefail
O=work/data/primekg; B=https://dataverse.harvard.edu/api/access/datafile
[ -s $O/nodes.csv ] || { curl -sfL "$B/6180617?format=original" -o $O/nodes.csv.tmp && mv $O/nodes.csv.tmp $O/nodes.csv; }
echo nodes $(wc -l < $O/nodes.csv)
[ -s $O/edges.npz ] || curl -sfL "$B/6180616" | (python3 pipelines/10b_parse_primekg_edges.py)
echo finished
