#!/usr/bin/env python
"""Figure 6: size and purity of the negative set of each policy (A), and loss in per-disease AP against the
number of held-out treatments negated per partition (B), in both graphs."""
import json, os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as S
import matplotlib.pyplot as plt

ROLE = sys.argv[1] if len(sys.argv) > 1 else "work/results_role"
PKG = sys.argv[2] if len(sys.argv) > 2 else "work/results_primekg"
OUT = sys.argv[3] if len(sys.argv) > 3 else "work/figures"
S.setup()
POLS = ["flat_negative", "typed", "typed_scoped", "typed_role", "typed_role_scoped", "mask_all"]
SHORT = {"flat_negative": "Flat", "typed": "Typed", "typed_scoped": "+scope", "typed_role": "+role", "typed_role_scoped": "+role+scope", "mask_all": "Mask"}
NAMED = {"flat_negative": "sci=negate,other=negate", "mask_all": "sci=mask,other=mask", "typed": "sci=negate,other=mask",
         "typed_scoped": "typed_scoped", "typed_role": "typed_role", "typed_role_scoped": "typed_role_scoped"}
pur = json.load(open(f"{ROLE}/negation_purity.json"))
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw={"width_ratios": [1, 1.25]})
for g, gname in [("hetionet", "Hetionet"), ("primekg", "PrimeKG")]:
    P = pur[g]; xs = [P[k]["negated_pairs"] for k in POLS[:-1]]; ys = [P[k]["pct"] for k in POLS[:-1]]
    a1.plot(xs, ys, S.MK[gname], color=S.COL[gname], ms=4, label=gname, ls="none")
    for k, x, y in zip(POLS[:-1], xs, ys):
        off = {("hetionet", "typed_scoped"): (-22, 5), ("hetionet", "typed_role"): (5, -2), ("hetionet", "typed_role_scoped"): (4, -10), ("primekg", "typed"): (4, -6), ("primekg", "typed_scoped"): (-24, 5)}.get((g, k), (4, 2))
        a1.annotate(SHORT[k], (x, y), xytext=off, textcoords="offset points", fontsize=5.8, color=S.MUTED)
a1.set_xscale("log"); a1.set_ylim(0, 50); a1.set_xlim(80, 9000)
a1.set_xticks([100, 300, 1000, 3000]); a1.set_xticklabels(["100", "300", "1,000", "3,000"]); a1.minorticks_off()
a1.set_xlabel("Pairs written back as negatives (log scale)"); a1.set_ylabel("Recorded treatments or approved\nindications among them (%)")
a1.legend(loc="lower left"); S.panel(a1, "A")
SETS = [(ROLE, "random", "graph", "Hetionet, random edge", "#2a78d6", "o"), (ROLE, "compound", "graph", "Hetionet, compound disjoint", "#2a78d6", "^"),
        (PKG, "random", "mf", "PrimeKG, random edge (MF)", "#eb6834", "s"), (PKG, "compound", "graph", "PrimeKG, compound disjoint", "#eb6834", "D")]
for res, task, model, lab, col, mk in SETS:
    B = json.load(open(f"{res}/boot_{task}.json")); M = pd.read_csv(f"{res}/e1_{task}_metrics.tsv", sep="\t")
    pk = res == PKG
    xs, ys = [], []
    for k in POLS:
        d = M[(M.policy == (k if pk else NAMED[k])) & (M.model == model)]
        n = d.test_pos_negated.mean(); c = B["contrasts"][f"{k} - no_writeback [{model}]"]["macroAP"]
        S.errorbar(a2, n, c["diff"], c["lo"], c["hi"], col, mk, label=lab if k == POLS[0] else None, horizontal=False, ms=3.4)
        xs.append(n); ys.append(c["diff"])
    order = sorted(range(len(xs)), key=lambda i: xs[i]); a2.plot([xs[i] for i in order], [ys[i] for i in order], color=col, lw=0.5, alpha=0.6)
a2.axhline(0, color=S.LIGHT, lw=0.6)
a2.set_xlabel("Held-out treatments written back as negatives per partition"); a2.set_ylabel("Δ per-disease AP vs no write-back")
a2.legend(loc="lower left", fontsize=6); S.panel(a2, "B")
fig.tight_layout(w_pad=2.5)
S.save(fig, OUT, "Figure6_negation_size")
print("saved")
