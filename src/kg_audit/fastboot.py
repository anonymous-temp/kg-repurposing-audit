"""Disease-cluster bootstrap of pooled and per-disease AP without storing score vectors.

pooled_ap_draws(y, s, dcode, W): W is a (B, n_dis) matrix of disease multiplicities; returns the
pooled AP of every draw, with row weight = multiplicity of the row's disease and exact tie handling
(identical to wb_lib.weighted_ap with those row weights)."""
import numpy as np


def pooled_ap_draws(y, s, dcode, W):
    o = np.argsort(-s, kind="stable"); ss = s[o]; yy = y[o] > 0; dd = dcode[o]
    n = len(ss)
    ends = np.r_[np.flatnonzero(np.diff(ss)), n - 1]
    P = np.flatnonzero(yy)
    e = ends[np.searchsorted(ends, P)]                      # end of each positive's tie group
    nD = W.shape[1]
    C = np.zeros((len(P), nD), np.float32)
    ordd = np.argsort(dd, kind="stable"); cnt = np.bincount(dd, minlength=nD); st = np.r_[0, np.cumsum(cnt)[:-1]]
    for d in np.flatnonzero(cnt):
        pos_d = ordd[st[d]:st[d] + cnt[d]]                   # sorted positions of disease d
        C[:, d] = np.searchsorted(pos_d, e, side="right")
    dpos = dd[P]
    wp = W[:, dpos]                                          # (B, nP)
    cum = np.cumsum(wp, 1)
    last = np.searchsorted(e, e, side="right") - 1           # last positive in the same tie group
    num = cum[:, last]
    den = W.astype(np.float32) @ C.T
    tot = wp.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (wp * num / np.maximum(den, 1e-12)).sum(1) / np.maximum(tot, 1e-12)


def per_disease_ap(y, s, dcode, nD, ap_fn):
    vals = np.full(nD, np.nan)
    o = np.argsort(dcode, kind="stable"); cnt = np.bincount(dcode, minlength=nD); st = np.r_[0, np.cumsum(cnt)[:-1]]
    for d in np.flatnonzero(cnt):
        idx = o[st[d]:st[d] + cnt[d]]
        yy = y[idx]
        if 0 < yy.sum() < len(yy):
            vals[d] = ap_fn(yy, s[idx])
    return vals


def macro_draws(vals, W):
    ok = ~np.isnan(vals)
    return (W[:, ok] * vals[ok]).sum(1) / np.maximum(W[:, ok].sum(1), 1e-12)
