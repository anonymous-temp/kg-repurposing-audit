#!/usr/bin/env python
"""Model-free checks for any source of negative drug-disease labels (see kg_audit.label_audit).

Usage: 25_audit_negative_source.py NEGATIVES.tsv POSITIVES.tsv OUT.json [N_PERM]
  NEGATIVES.tsv  columns drug, disease and optionally n_tests (how often the pair was tested)
  POSITIVES.tsv  columns drug, disease (known positives on the same identifiers)
Writes contamination, contamination by testing strata and the degree-preserving popularity null.
"""
import json, os, sys
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from kg_audit.label_audit import by_strata, contamination, popularity_null  # noqa: E402

NEG, POS, OUT = sys.argv[1:4]
N_PERM = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
BINS = [(0, 0, "0"), (1, 1, "1"), (2, 2, "2"), (3, 4, "3-4"), (5, 9, "5-9"), (10, 19, "10-19"), (20, 10**9, ">=20")]
neg = pd.read_csv(NEG, sep="\t", dtype={"drug": str, "disease": str})
pos = pd.read_csv(POS, sep="\t", dtype={"drug": str, "disease": str})
negs = list(zip(neg.drug, neg.disease))
poss = set(zip(pos.drug, pos.disease))
res = {"contamination": contamination(negs, poss)}
if "n_tests" in neg:
    counts = dict(zip(negs, neg.n_tests))
    res["by_tests"] = by_strata(negs, poss, counts, BINS)
res["popularity_null"] = popularity_null(negs, poss, n_perm=N_PERM)
json.dump(res, open(OUT, "w"), indent=1)
print(json.dumps(res, indent=1))
