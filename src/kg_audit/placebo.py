"""Placebo negative sets for write-back experiments.

A placebo replaces a set of negated (compound, disease) pairs by the same number of other pairs:

- degree-matched: double-edge swaps that keep every compound's and every disease's number of negatives;
  a swap is accepted only if both new pairs are allowed (unlabelled, not palliative or off-label, no stopped trial);
- uniform: a uniform random sample of allowed pairs of the same size.

Pairs are integer (row, column) indices into an allowed-pair matrix.
"""
import numpy as np


def uniform_sample(allowed, size, rng):
    """Return a set of `size` distinct allowed (row, col) index pairs drawn uniformly."""
    n_cols = allowed.shape[1]
    idx = rng.choice(np.flatnonzero(allowed.ravel()), size=size, replace=False)
    return {(int(i // n_cols), int(i % n_cols)) for i in idx}


def rewire(edges, allowed, rng, swaps_per_edge=30):
    """Degree-preserving rewiring of `edges` (list of (row, col)) into allowed cells.

    Returns (rewired pairs that are not among the original edges, number of original edges never replaced).
    Original edges that could not be swapped out are dropped and counted.
    """
    E = [tuple(e) for e in edges]
    cur = np.zeros(allowed.shape, bool)
    for a, b in E:
        cur[a, b] = True
    orig = set(E); n = len(E)
    for _ in range(swaps_per_edge * n):
        i, j = rng.integers(n, size=2)
        (a, b), (c, d) = E[i], E[j]
        if a == c or b == d or cur[a, d] or cur[c, b] or not allowed[a, d] or not allowed[c, b]:
            continue
        cur[a, b] = cur[c, d] = False; cur[a, d] = cur[c, b] = True
        E[i], E[j] = (a, d), (c, b)
    left = sum(1 for e in E if e in orig)
    return {e for e in E if e not in orig}, left
