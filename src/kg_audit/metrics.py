"""Candidate and metric contracts shared by the corrected analyses."""

from collections import defaultdict
import math
import numpy as np


def fit_degree_reference(train_positive, train_negative, seed=1):
    """Fit a balanced L2 logistic reference using only supplied training labels.

    Return LIBLINEAR decision scores using explicit feature multiplication;
    rank metrics do not require a probability transform.
    """
    from collections import Counter
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    cdeg = Counter(c for c, d in train_positive)
    ddeg = Counter(d for c, d in train_positive)

    def features(pairs):
        return np.array([[np.log1p(cdeg[c]), np.log1p(ddeg[d])] for c, d in pairs])

    x = features(train_positive + train_negative)
    scaler = StandardScaler().fit(x)
    classifier = LogisticRegression(
        solver="liblinear", class_weight="balanced", C=1.0, max_iter=1000, random_state=seed
    ).fit(scaler.transform(x), np.r_[np.ones(len(train_positive)), np.zeros(len(train_negative))])
    if not np.isfinite(classifier.coef_).all():
        raise ValueError("Non-finite degree baseline coefficients")

    def predict(pairs):
        score = np.sum(scaler.transform(features(pairs)) * classifier.coef_[0], axis=1) + classifier.intercept_[0]
        if not np.isfinite(score).all():
            raise ValueError("Non-finite degree baseline predictions")
        return score

    return predict, cdeg, ddeg


def degree_bin(degree):
    """Zero is separate; positive bins are 1, 2-3, 4-7, ... ."""
    return -1 if degree == 0 else int(math.floor(math.log2(degree)))


def assert_pair_contract(pairs, labels, training_pairs):
    if len(pairs) != len(labels):
        raise ValueError("Pair and label counts differ")
    seen = {}
    for pair, label in zip(pairs, labels):
        pair = tuple(pair)
        if pair in seen and seen[pair] != label:
            raise ValueError("The same pair has conflicting labels")
        seen[pair] = label
    if set(seen).intersection(map(tuple, training_pairs)):
        raise ValueError("Supervised fitting and evaluation pairs overlap")


def ap_function(labels, scores):
    """Build a fast weighted-average-precision function with exact tie handling.

    Repeated cluster draws enter through weights, so multiplicity is preserved.
    The unweighted call agrees with sklearn's non-interpolated average precision.
    """
    labels = np.asarray(labels, dtype=float)
    scores = np.asarray(scores, dtype=float)
    if labels.shape != scores.shape or not np.isfinite(scores).all():
        raise ValueError("Scores and labels must be finite aligned vectors")
    order = np.argsort(-scores, kind="stable")
    yy = labels[order]
    ends = np.r_[np.flatnonzero(np.diff(scores[order])), len(scores) - 1]

    def compute(weights=None):
        ww = np.ones(len(labels)) if weights is None else np.asarray(weights)[order]
        tp = np.cumsum(ww * yy)[ends]
        total = np.cumsum(ww)[ends]
        if tp[-1] <= 0 or total[-1] <= tp[-1]:
            return float("nan")
        precision = np.divide(tp, total, out=np.zeros_like(tp), where=total > 0)
        return float(np.sum(np.diff(np.r_[0.0, tp]) * precision) / tp[-1])

    return compute


def sample_candidates(positives, drugs, diseases, excluded, cdegree, ddegree, scheme, ratio=20, seed=1):
    """Draw evaluation pairs within the declared task universe.

    Uniform draws are uniform over eligible pairs. Matched draws preserve the
    positive pair's two training-degree bins. Both schemes sample with replacement
    and therefore have the same draw count; callers report unique support as well.
    Empty matching strata fail rather than silently switching candidate universes.
    """
    rng = np.random.default_rng(seed)
    excluded = set(map(tuple, excluded))
    drugs, diseases = sorted(set(drugs)), sorted(set(diseases))
    cb, db = defaultdict(list), defaultdict(list)
    for c in drugs:
        cb[degree_bin(cdegree.get(c, 0))].append(c)
    for d in diseases:
        db[degree_bin(ddegree.get(d, 0))].append(d)
    cache = {}
    out = []
    for c, d in positives:
        key = (degree_bin(cdegree.get(c, 0)), degree_bin(ddegree.get(d, 0)))
        if scheme == "uniform":
            # Rejection is uniform over eligible pairs, without a large Cartesian array.
            made, attempts = 0, 0
            while made < ratio:
                attempts += 1
                if attempts > 100000:
                    raise ValueError("Uniform candidate pool is exhausted")
                pair = (drugs[rng.integers(len(drugs))], diseases[rng.integers(len(diseases))])
                if pair not in excluded:
                    out.append(pair)
                    made += 1
        elif scheme == "matched":
            if key not in cache:
                cache[key] = [(a, b) for a in cb[key[0]] for b in db[key[1]] if (a, b) not in excluded]
            pool = cache[key]
            if not pool:
                raise ValueError(f"No unlabelled pair in degree stratum {key}")
            out.extend(pool[i] for i in rng.integers(len(pool), size=ratio))
        else:
            raise ValueError(f"Unknown candidate scheme: {scheme}")
    return out
