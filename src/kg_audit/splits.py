"""Validation labels follow the prediction task; test labels are never selected here."""

import numpy as np


def validation_partition(non_test_positive, task, seed, fraction=0.1):
    pairs = list(map(tuple, non_test_positive))
    if len(pairs) < 2:
        raise ValueError("At least two non-test positives are required")
    rng = np.random.default_rng(seed)
    if task == "random":
        chosen = set(rng.choice(len(pairs), max(1, int(fraction * len(pairs))), replace=False))
        validation = [p for i, p in enumerate(pairs) if i in chosen]
    elif task == "compound_disjoint":
        drugs = sorted({c for c, d in pairs})
        if len(drugs) < 2:
            raise ValueError("At least two drugs are needed for label-cold-start validation")
        selected = set(rng.choice(drugs, max(1, int(fraction * len(drugs))), replace=False))
        validation = [p for p in pairs if p[0] in selected]
    else:
        raise ValueError("Unsupported validation task")
    valset = set(validation)
    train = [p for p in pairs if p not in valset]
    if not train:
        raise ValueError("Validation holdout leaves no fitting positives")
    return train, validation
