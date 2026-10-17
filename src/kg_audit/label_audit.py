"""Model-free checks for a source of negative drug-disease labels.

A candidate negative-label source (for example pairs with a stopped trial) is checked against known
positives on the same identifiers:

* contamination: the share of candidate negatives that are known positives;
* testing gradient: contamination by the number of times a pair was tested;
* popularity null: the contamination expected if the candidate pairs were rewired at random while every
  drug and every disease keeps its number of candidate pairs (degree-preserving double-edge swaps).

The observed contamination divided by the null expectation measures how much of the overlap is specific
to the pairs, beyond drug and disease popularity. These checks need no model; the model-based check
(matched placebo negatives) is implemented in the write-back pipelines.
"""
import math

import numpy as np


def wilson(k, n, z=1.96):
    """Wilson score interval for k successes out of n."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


def contamination(negatives, positives):
    """Share of candidate negative pairs that are known positives, with a Wilson interval."""
    neg = set(negatives)
    k = len(neg & set(positives))
    lo, hi = wilson(k, len(neg))
    return {"n": len(neg), "k": k, "share": k / len(neg) if neg else float("nan"), "ci95": [lo, hi]}


def by_strata(negatives, positives, counts, bins):
    """Contamination within strata of a testing count. bins: list of (low, high, label), inclusive."""
    pos = set(positives)
    out = {}
    for lo, hi, label in bins:
        sub = [p for p in set(negatives) if lo <= counts.get(p, 0) <= hi]
        k = sum(p in pos for p in sub)
        a, b = wilson(k, len(sub))
        out[label] = {"n": len(sub), "k": k, "share": k / len(sub) if sub else float("nan"), "ci95": [a, b]}
    return out


def rewire(edges, rng, n_swaps=None, forbidden=frozenset()):
    """Degree-preserving double-edge swaps on a bipartite edge list; returns a new list of edges."""
    E = list(edges)
    cur = set(E)
    n = len(E)
    for _ in range(n_swaps if n_swaps is not None else 10 * n):
        i, j = rng.integers(n, size=2)
        (a, b), (c, d) = E[i], E[j]
        if a == c or b == d:
            continue
        if (a, d) in cur or (c, b) in cur or (a, d) in forbidden or (c, b) in forbidden:
            continue
        cur.discard((a, b)); cur.discard((c, d)); cur.add((a, d)); cur.add((c, b))
        E[i], E[j] = (a, d), (c, b)
    return E


def popularity_null(negatives, positives, n_perm=1000, seed=20261017, forbidden=frozenset()):
    """Observed number of candidate negatives that are positives against a degree-preserving null.

    Returns observed and null counts, the fold enrichment over the base rate of positives on the
    drug x disease grid spanned by the candidates, and the ratio of observed to null counts.
    """
    neg = sorted(set(negatives))
    pos = set(positives)
    drugs = {a for a, _ in neg}
    dis = {b for _, b in neg}
    obs = sum(p in pos for p in neg)
    grid_pos = sum(1 for a, b in pos if a in drugs and b in dis)
    base = grid_pos / (len(drugs) * len(dis))
    rng = np.random.default_rng(seed)
    null = np.array([sum(p in pos for p in rewire(neg, rng, forbidden=forbidden)) for _ in range(n_perm)])
    mean = float(null.mean())
    return {"pairs": len(neg), "observed": int(obs), "null_mean": mean,
            "null_lo": float(np.percentile(null, 2.5)), "null_hi": float(np.percentile(null, 97.5)),
            "base_rate": base, "fold_observed": obs / len(neg) / base if base else float("nan"),
            "fold_null": mean / len(neg) / base if base else float("nan"),
            "ratio_observed_to_null": obs / mean if mean else float("nan"),
            "ratio_lo": obs / np.percentile(null, 97.5) if np.percentile(null, 97.5) else float("nan"),
            "ratio_hi": obs / max(np.percentile(null, 2.5), 1e-9),
            "p_one_sided": float((np.sum(null >= obs) + 1) / (len(null) + 1))}
