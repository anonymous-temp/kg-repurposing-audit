import numpy as np

from kg_audit.label_audit import by_strata, contamination, popularity_null, rewire, wilson


def test_wilson_bounds():
    lo, hi = wilson(5, 10)
    assert 0 < lo < 0.5 < hi < 1


def test_contamination_and_strata():
    neg = [("a", "x"), ("a", "y"), ("b", "x"), ("c", "z")]
    pos = {("a", "x"), ("c", "z"), ("d", "w")}
    c = contamination(neg, pos)
    assert c["n"] == 4 and c["k"] == 2 and c["share"] == 0.5
    s = by_strata(neg, pos, {("a", "x"): 1, ("a", "y"): 1, ("b", "x"): 5, ("c", "z"): 5}, [(1, 1, "1"), (2, 9, "2-9")])
    assert s["1"]["k"] == 1 and s["2-9"]["k"] == 1 and s["2-9"]["n"] == 2


def test_rewire_preserves_degrees():
    rng = np.random.default_rng(1)
    E = [(i, j) for i in range(6) for j in range(4) if (i + j) % 2 == 0]
    R = rewire(E, rng)
    assert len(set(R)) == len(E)
    for k in range(6):
        assert sum(a == k for a, _ in R) == sum(a == k for a, _ in E)
    for k in range(4):
        assert sum(b == k for _, b in R) == sum(b == k for _, b in E)


def test_popularity_null_detects_pair_specific_overlap():
    # every candidate negative is a positive, but positives are spread so that rewiring rarely hits them
    neg = [(i, i) for i in range(30)]
    pos = set(neg)
    r = popularity_null(neg, pos, n_perm=50, seed=3)
    assert r["observed"] == 30 and r["null_mean"] < 5 and r["ratio_observed_to_null"] > 5
