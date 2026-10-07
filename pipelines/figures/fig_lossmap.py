#!/usr/bin/env python
"""Figure 9: where the loss from negative labels falls, and why rarely tested negatives do not help.

A, B  Hetionet graph head: per-disease change in AP under flat negatives (vs no write-back) against the share of the
      disease's held-out treatments that the policy negated (random-edge and compound-disjoint tasks)
C     Hetionet graph head, random-edge task: mean reciprocal rank of held-out treatments by their stopped trials, per policy
D     percentile rank on the external grid (graph head fitted without write-back) of unlabelled written-back pairs with
      an efficacy or safety stop, by number of registered trials before 2015; dashed lines: median of external approved
      indications
"""
import json, os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as S
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SEL = sys.argv[1] if len(sys.argv) > 1 else "paper_results/selection"
STRATA = sys.argv[2] if len(sys.argv) > 2 else "paper_results/role/e1_strata.tsv"
OUT = sys.argv[3] if len(sys.argv) > 3 else "work/figures"
S.setup()
fig = plt.figure(figsize=(6.9, 5.2))
gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.42)

# ---------------------------------------------------------------- A, B
for k, (task, title) in enumerate([("random", "Hetionet, random edge"), ("compound", "Hetionet, compound disjoint")]):
    ax = fig.add_subplot(gs[0, k])
    d = pd.read_csv(f"{SEL}/disease_changes_{task}.tsv", sep="\t")
    rho = spearmanr(d.share_negated, d.flat_minus_none).correlation
    ax.axhline(0, color=S.LIGHT, lw=0.8, zorder=0)
    ax.scatter(d.share_negated * 100, d.flat_minus_none, s=6 + 6 * d.mean_pos, color=S.COL["graph"], alpha=0.55, lw=0.4, edgecolor="white")
    # running median over the share negated
    b = pd.cut(d.share_negated, [-0.01, 0.2, 0.4, 0.6, 0.8, 1.0])
    m = d.groupby(b, observed=True).agg(x=("share_negated", "median"), y=("flat_minus_none", "median"))
    ax.plot(m.x * 100, m.y, "-", color=S.INK, lw=1.0)
    ax.set_xlim(-4, 104); ax.set_ylim(-0.85, 0.15)
    ax.set_xlabel("Held-out treatments of the disease negated (%)"); ax.set_ylabel("Δ per-disease AP, flat vs none" if k == 0 else "")
    S.panel(ax, "AB"[k], title=f"{title}: {len(d)} diseases, Spearman ρ = {rho:.2f}".replace("-", "−"))
    if k == 0:
        ax.legend([Line2D([], [], color=S.INK, lw=1.0), Line2D([], [], marker="o", ls="none", color=S.COL["graph"], alpha=0.55, ms=4)],
                  ["Median by bin", "Disease (size: held-out treatments)"], loc="lower left")
    print(task, "rho", round(rho, 3))

# ---------------------------------------------------------------- C
ax = fig.add_subplot(gs[1, 0])
st = pd.read_csv(STRATA, sep="\t")
st = st[(st.task == "random") & (st.model == "graph")]
POL = ["no_writeback", "flat_negative", "typed", "typed_role", "typed_role_scoped", "mask_all"]
PL = ["None", "Flat", "Typed", "+role", "+role\n+scope", "Mask"]
SPEC = [("no stop", "No stop", S.MUTED, "o", None), ("other stop", "Other stops only", "#eda100", "^", 0.530),
        ("scientific stop, investigational", "Sci. stop, investigational", "#4a3aa7", "D", 0.560),
        ("scientific stop, other role", "Sci. stop, other role", "#e87ba4", "s", 0.590)]
for s, lab, col, mk, ylab in SPEC:
    r = st[st.stratum == s].set_index("policy").reindex(POL)
    ax.plot(range(len(POL)), r.MRR, "-", marker=mk, color=col, ms=3.4, lw=0.9)
    yl = r.MRR.iloc[-1] if ylab is None else ylab
    ax.plot([len(POL) - 0.55], [yl], mk, color=col, ms=3.0)
    ax.text(len(POL) - 0.38, yl, lab + f" ({r.n.iloc[0]:.0f})", fontsize=5.8, va="center", color=S.INK)
ax.set_xticks(range(len(POL))); ax.set_xticklabels(PL); ax.set_xlim(-0.3, len(POL) + 2.4)
ax.set_ylim(0.2, 0.62); ax.set_ylabel("Mean reciprocal rank of held-out treatments")
ax.set_xlabel("Write-back policy"); ax.spines["bottom"].set_bounds(0, len(POL) - 1)
S.panel(ax, "C", title="Hetionet graph head, random edge")

# ---------------------------------------------------------------- D
ax = fig.add_subplot(gs[1, 1])
BINS = [(1, 1, "1"), (2, 2, "2"), (3, 4, "3–4"), (5, 9, "5–9"), (10, 19, "10–19"), (20, 10 ** 6, "≥20")]
for gi, g in enumerate(["hetionet", "primekg"]):
    R = json.load(open(f"{SEL}/rank_strata_{g}.json"))["graph"]
    df = pd.DataFrame(R["stopped"], columns=["n", "pct", "ext"])
    gname = "Hetionet" if g == "hetionet" else "PrimeKG"
    pos = np.arange(len(BINS)) + (gi - 0.5) * 0.34
    data = [df[(df.n >= lo) & (df.n <= hi)].pct.values for lo, hi, _ in BINS]
    bp = ax.boxplot(data, positions=pos, widths=0.28, patch_artist=True, showfliers=False, whis=(10, 90),
                    medianprops=dict(color="white", lw=1.0), boxprops=dict(lw=0), whiskerprops=dict(color=S.COL[gname], lw=0.8),
                    capprops=dict(color=S.COL[gname], lw=0.8))
    for p in bp["boxes"]:
        p.set_facecolor(S.COL[gname])
    ax.axhline(np.median(R["approved"]), color=S.COL[gname], ls=(0, (3, 2)), lw=0.8)
    print(g, [len(x) for x in data], [round(float(np.median(x)), 1) for x in data], "approved", round(float(np.median(R["approved"])), 1))
ax.set_xticks(range(len(BINS))); ax.set_xticklabels([b[2] for b in BINS]); ax.set_ylim(0, 100)
ax.set_xlabel("Registered trials before 2015"); ax.set_ylabel("Percentile rank on the external grid")
ax.legend([plt.Rectangle((0, 0), 1, 1, color=S.COL["Hetionet"]), plt.Rectangle((0, 0), 1, 1, color=S.COL["PrimeKG"]),
           Line2D([], [], color=S.MUTED, ls=(0, (3, 2)), lw=0.8)], ["Hetionet", "PrimeKG", "Approved indications (median)"],
          loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, handlelength=1.4)
S.panel(ax, "D", title="Stopped pairs before write-back, graph head")
S.save(fig, OUT, "Figure9_where_loss_falls")
