import numpy as np
from kg_audit import placebo as P


def _case(seed=0, n_rows=50, n_cols=30, n_edges=200):
    rs = np.random.default_rng(seed)
    allowed = rs.random((n_rows, n_cols)) > 0.25
    edges = sorted({(int(a), int(b)) for a, b in zip(rs.integers(n_rows, size=n_edges), rs.integers(n_cols, size=n_edges))})
    for e in edges:                       # the real negated pairs are stopped pairs, which placebos may not use
        allowed[e] = False
    return edges, allowed


def test_rewire_preserves_row_and_column_counts():
    edges, allowed = _case()
    new, left = P.rewire(edges, allowed, np.random.default_rng(1))
    assert left == 0
    rows0 = np.bincount([a for a, b in edges], minlength=allowed.shape[0]); cols0 = np.bincount([b for a, b in edges], minlength=allowed.shape[1])
    rows1 = np.bincount([a for a, b in new], minlength=allowed.shape[0]); cols1 = np.bincount([b for a, b in new], minlength=allowed.shape[1])
    assert (rows0 == rows1).all() and (cols0 == cols1).all()


def test_rewire_uses_only_allowed_new_cells_and_is_reproducible():
    edges, allowed = _case(2)
    a = P.rewire(edges, allowed, np.random.default_rng(3)); b = P.rewire(edges, allowed, np.random.default_rng(3))
    assert a == b
    assert all(allowed[e] for e in a[0]) and not set(a[0]) & set(edges)


def test_rewire_reports_unreplaced_edges():
    allowed = np.zeros((3, 3), bool)      # nothing allowed: no swap can be accepted
    edges = [(0, 0), (1, 1)]
    new, left = P.rewire(edges, allowed, np.random.default_rng(0), swaps_per_edge=5)
    assert new == set() and left == 2


def test_uniform_sample():
    _, allowed = _case(4)
    s = P.uniform_sample(allowed, 40, np.random.default_rng(0))
    assert len(s) == 40 and all(allowed[e] for e in s)
    assert s == P.uniform_sample(allowed, 40, np.random.default_rng(0))
