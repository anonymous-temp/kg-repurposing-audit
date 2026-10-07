#!/usr/bin/env python
"""Figure 7: what makes negative labels harmful (both graphs).

A  size and purity of the negative set of each policy
B  change in per-disease AP (vs no write-back) against held-out treatments negated per partition
C-F  change in per-disease AP relative to masking when negatives are added on top of masking: negated sets of the
     real policies (stopped pairs) versus degree-matched and uniform placebo sets of the same size
"""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as S
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROLE = sys.argv[1] if len(sys.argv) > 1 else "work/results_role"
ROLE2 = os.environ.get("ROLE2", "work/results_role2")      # role-policy runs, if kept separately from ROLE
PKG = sys.argv[2] if len(sys.argv) > 2 else "work/results_primekg"
OUT = sys.argv[3] if len(sys.argv) > 3 else "work/figures"
XHET = sys.argv[4] if len(sys.argv) > 4 else "work/results_extra"
XPKG = sys.argv[5] if len(sys.argv) > 5 else "work/results_extra_primekg"
SEL = sys.argv[6] if len(sys.argv) > 6 else "work/results_selection"
S.setup()
POLS = ["flat_negative", "typed", "typed_scoped", "typed_role", "typed_role_scoped", "typed_lowint2", "mask_all"]
SHORT = {"flat_negative": "Flat", "typed": "Typed", "typed_scoped": "+scope", "typed_role": "+role", "typed_role_scoped": "+role+scope",
         "typed_lowint2": "≤2 trials", "mask_all": "Mask"}
NAMED = {"flat_negative": "sci=negate,other=negate", "mask_all": "sci=mask,other=mask", "typed": "sci=negate,other=mask",
         "typed_scoped": "typed_scoped", "typed_role": "typed_role", "typed_role_scoped": "typed_role_scoped", "typed_lowint2": "typed_lowint2",
         "no_writeback": "sci=ignore,other=ignore"}
pur = json.load(open(f"{SEL}/negation_purity.json"))

fig = plt.figure(figsize=(7.0, 5.6))
gs = fig.add_gridspec(2, 4, height_ratios=[1.05, 1], hspace=0.62, wspace=0.55)
a1 = fig.add_subplot(gs[0, 0:2]); a2 = fig.add_subplot(gs[0, 2:4])

# ---------------------------------------------------------------- A
for g, gname in [("hetionet", "Hetionet"), ("primekg", "PrimeKG")]:
    P = pur[g]; ks = POLS[:-1]; xs = [P[k]["negated_pairs"] for k in ks]; ys = [P[k]["pct"] for k in ks]
    a1.plot(xs, ys, S.MK[gname], color=S.COL[gname], ms=4, label=gname, ls="none")
    for k, x, y in zip(ks, xs, ys):
        off = {("hetionet", "typed_scoped"): (-22, 5), ("hetionet", "typed_role"): (5, -2), ("hetionet", "typed_role_scoped"): (4, -10),
               ("primekg", "typed"): (4, -6), ("primekg", "typed_scoped"): (-24, 5), ("primekg", "typed_lowint2"): (4, -3)}.get((g, k), (4, 2))
        a1.annotate(SHORT[k], (x, y), xytext=off, textcoords="offset points", fontsize=5.8, color=S.MUTED)
a1.set_xscale("log"); a1.set_ylim(0, 50); a1.set_xlim(40, 9000)
a1.set_xticks([100, 300, 1000, 3000]); a1.set_xticklabels(["100", "300", "1,000", "3,000"]); a1.minorticks_off()
a1.set_xlabel("Pairs written back as negatives (log scale)"); a1.set_ylabel("Recorded or approved\nindications among them (%)")
a1.legend(loc="upper left"); S.panel(a1, "A", title="Size and purity of the negative set")


# ---------------------------------------------------------------- B
def metrics(res, task, extra):
    M = pd.read_csv(f"{res}/e1_{task}_metrics.tsv", sep="\t")
    if res == ROLE:
        if os.path.exists(f"{ROLE2}/e1_{task}_metrics.tsv"):
            M = pd.concat([M, pd.read_csv(f"{ROLE2}/e1_{task}_metrics.tsv", sep="\t")]).drop_duplicates(["pseed", "policy", "model", "seed", "w_neg"])
        M = M[M.w_neg == 10]
    X = pd.read_csv(f"{extra}/e1_{task}_metrics.tsv", sep="\t")
    return pd.concat([M, X], ignore_index=True)


SETS = [(ROLE, XHET, "random", "graph", "Hetionet, random edge", "#2a78d6", "o"), (ROLE, XHET, "compound", "graph", "Hetionet, compound disjoint", "#2a78d6", "^"),
        (PKG, XPKG, "random", "mf", "PrimeKG, random edge (MF)", "#eb6834", "s"), (PKG, XPKG, "compound", "graph", "PrimeKG, compound disjoint", "#eb6834", "D")]
for res, xres, task, model, lab, col, mk in SETS:
    B = json.load(open(f"{res}/boot_{task}.json")); BX = json.load(open(f"{xres}/boot_{task}.json")); M = metrics(res, task, xres)
    pk = res == PKG; xs, ys = [], []
    for k in POLS:
        pol = k if (pk or k == "typed_lowint2") else NAMED[k]
        d = M[(M.policy == pol) & (M.model == model)]
        n = d.test_pos_negated.mean()
        c = (BX if k == "typed_lowint2" else B)["contrasts"][f"{k} - no_writeback [{model}]"]["macroAP"]
        S.errorbar(a2, n, c["diff"], c["lo"], c["hi"], col, mk, label=lab if k == POLS[0] else None, horizontal=False, ms=3.4)
        xs.append(n); ys.append(c["diff"])
    order = sorted(range(len(xs)), key=lambda i: xs[i]); a2.plot([xs[i] for i in order], [ys[i] for i in order], color=col, lw=0.5, alpha=0.6)
a2.axhline(0, color=S.LIGHT, lw=0.6)
a2.set_xlabel("Held-out treatments written back as negatives\nper partition"); a2.set_ylabel("Δ per-disease AP vs no write-back")
a2.legend(loc="lower left", fontsize=6); S.panel(a2, "B", title="Loss against held-out treatments negated")

# ---------------------------------------------------------------- C-F
REAL = {"hetionet": ["flat_negative", "typed", "typed_scoped", "typed_role", "typed_role_scoped", "typed_lowint2"],
        "primekg": ["flat_negative", "typed", "typed_role_scoped", "typed_lowint2"]}
for i, (res, xres, task, model, lab, col, mk) in enumerate(SETS):
    ax = fig.add_subplot(gs[1, i]); g = "primekg" if res == PKG else "hetionet"
    B = json.load(open(f"{res}/boot_{task}.json")); BX = json.load(open(f"{xres}/boot_{task}.json")); M = metrics(res, task, xres)
    for kind, mfc, ls in [("real", col, "-"), ("degree", "white", "--"), ("uniform", None, ":")]:
        xs, ys = [], []
        for k in REAL[g]:
            name = k if kind == "real" else f"placebo_{kind}__{k}"
            pol = name if (res == PKG or name.startswith(("placebo", "typed_lowint"))) else NAMED[k]
            d = M[(M.policy == pol) & (M.model == model)]
            n = d.negated.mean()
            c = BX["contrasts"][f"{name} - mask_all [{model}]"]["macroAP"]
            if kind == "uniform":
                ax.errorbar(n, c["diff"], yerr=[[c["diff"] - c["lo"]], [c["hi"] - c["diff"]]], fmt="x", ms=3.2, color=S.MUTED, elinewidth=0.7, capsize=1.3, capthick=0.7)
            else:
                ax.errorbar(n, c["diff"], yerr=[[c["diff"] - c["lo"]], [c["hi"] - c["diff"]]], fmt=mk, ms=3.2, color=col, mfc=mfc, mec=col, elinewidth=0.7, capsize=1.3, capthick=0.7)
            xs.append(n); ys.append(c["diff"])
        o = np.argsort(xs); ax.plot(np.array(xs)[o], np.array(ys)[o], color=S.MUTED if kind == "uniform" else col, lw=0.6, ls=ls)
    ax.axhline(0, color=S.LIGHT, lw=0.6); ax.set_xscale("log"); ax.minorticks_off()
    ax.set_xticks([100, 1000]); ax.set_xticklabels(["100", "1,000"])
    ax.set_xlabel("Negatives added"); ax.set_title(lab.replace(", ", "\n"), fontsize=6.3, pad=3)
    if i == 0:
        ax.set_ylabel("Δ per-disease AP vs masking")
    S.panel(ax, "CDEF"[i], y=1.16)
hand = [Line2D([], [], marker="o", color="k", mfc="k", ls="-", lw=0.6, ms=3.2, label="Stopped pairs (real policies)"),
        Line2D([], [], marker="o", color="k", mfc="white", ls="--", lw=0.6, ms=3.2, label="Placebo, same drug and disease counts"),
        Line2D([], [], marker="x", color=S.MUTED, ls=":", lw=0.6, ms=3.2, label="Placebo, uniform")]
fig.legend(handles=hand, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.02))
S.save(fig, OUT, "Figure7_negatives")
print("saved")
