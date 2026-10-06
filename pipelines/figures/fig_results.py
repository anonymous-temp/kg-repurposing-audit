#!/usr/bin/env python
"""Figures 3-4 and supplementary figures S1-S3 from the write-back results."""
import json, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = sys.argv[1] if len(sys.argv) > 1 else "work/results"
OUT = sys.argv[2] if len(sys.argv) > 2 else "work/figures"
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as FS
FS.setup()
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
COL = {"graph": "#2a78d6", "hybrid": "#eb6834", "mf": "#1baf7a"}            # palette slots 1-3 (validated all-pairs)
MK = {"graph": "o", "hybrid": "s", "mf": "D"}
LAB = {"graph": "Graph head", "hybrid": "Hybrid (MF + graph)", "mf": "Label-only MF"}
POL = ["flat_negative", "mask_all", "typed", "typed_scoped"]
POL_LAB = {"flat_negative": "Flat negative", "mask_all": "Mask all", "typed": "Typed", "typed_scoped": "Typed + scope"}
NAMED = {"no_writeback": "sci=ignore,other=ignore", "flat_negative": "sci=negate,other=negate", "mask_all": "sci=mask,other=mask",
         "typed": "sci=negate,other=mask", "typed_scoped": "typed_scoped"}


def clean(ax):
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    ax.tick_params(length=2, colors=INK)


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def contrast_panel(ax, boot, metric, title, show_legend=False):
    y0 = np.arange(len(POL))[::-1]
    for k, model in enumerate(["graph", "hybrid", "mf"]):
        off = (1 - k) * 0.22
        for yi, pol in zip(y0, POL):
            c = boot["contrasts"][f"{pol} - no_writeback [{model}]"][metric]
            ax.plot([c["lo"], c["hi"]], [yi + off] * 2, color=COL[model], lw=1.2, solid_capstyle="round")
            ax.plot(c["diff"], yi + off, MK[model], ms=4.2, color=COL[model], mec="white", mew=0.8,
                    label=LAB[model] if yi == y0[0] else None)
    ax.axvline(0, color=MUTED, lw=0.7)
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.set_yticks(y0); ax.set_yticklabels([POL_LAB[p] for p in POL])
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.set_title(title, loc="left", fontsize=7.5, fontweight="bold")
    clean(ax)
    if show_legend:
        ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.22), fontsize=6.0, handletextpad=0.3, ncol=1)


# ---------------- Figure 3: E1
boots = {t: json.load(open(f"{RES}/boot_{t}.json")) for t in ["random", "compound"] if os.path.exists(f"{RES}/boot_{t}.json")}
if len(boots) == 2:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.6), sharey=True)
    for j, (t, tl) in enumerate([("random", "Random-edge task"), ("compound", "Compound-disjoint task")]):
        for i, (m, ml) in enumerate([("pooledAP", "pooled AP"), ("macroAP", "per-disease AP")]):
            contrast_panel(axes[i, j], boots[t], m, f"{'abcd'[2 * i + j]}  {tl}: change in {ml}")
            axes[i, j].set_xlabel(f"Difference from no write-back ({ml})", color=MUTED)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.tight_layout(h_pad=1.6, w_pad=2.0, rect=(0, 0.05, 1, 1))
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=6.5, handletextpad=0.3)
    save(fig, "Figure3_E1_policies")

# ---------------- Figure 4: E2 / E3
if os.path.exists(f"{RES}/boot_e2.json") and os.path.exists(f"{RES}/e2_metrics.tsv"):
    b2 = json.load(open(f"{RES}/boot_e2.json"))
    e2 = pd.read_csv(f"{RES}/e2_metrics.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 2.7))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 0.9], wspace=0.55)
    ax = fig.add_subplot(gs[0]); contrast_panel(ax, b2, "macroAP", "a  Approved indications absent\n    from Hetionet: per-disease AP", show_legend=True)
    ax.set_xlabel("Difference from no write-back", color=MUTED)
    ax = fig.add_subplot(gs[1])
    strata = [("in_wb_scientific", "scientific stop before 2015"), ("in_wb_other", "other stop before 2015"), ("not_in_wb", "no stop before 2015")]
    pols = ["no_writeback"] + POL
    y0 = np.arange(len(pols))[::-1]
    sc = {"in_wb_scientific": "#e34948", "in_wb_other": "#eda100", "not_in_wb": "#4a3aa7"}
    smk = {"in_wb_scientific": "o", "in_wb_other": "s", "not_in_wb": "^"}
    for k, (s, sl) in enumerate(strata):
        vals = [e2[(e2.policy == NAMED[p]) & (e2.w_neg == 10) & (e2.model == "graph")][f"median_pct_{s}"].mean() * 100 for p in pols]
        n = int(e2[f"n_{s}"].iloc[0])
        ax.plot(vals, y0 + (1 - k) * 0.22, smk[s], ms=4.2, color=sc[s], mec="white", mew=0.8, label=f"{sl} (n = {n})")
    ax.set_yticks(y0); ax.set_yticklabels(["No write-back"] + [POL_LAB[p] for p in POL])
    ax.set_xlabel("Median percentile rank (graph head)", color=MUTED)
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True); clean(ax)
    ax.set_title("b  Rank of approved indications\n    by write-back stratum", loc="left", fontsize=7.5, fontweight="bold")
    ax.legend(frameon=False, fontsize=5.8, loc="upper center", bbox_to_anchor=(0.45, -0.22), ncol=1, handletextpad=0.3)
    ax = fig.add_subplot(gs[2])
    for k, model in enumerate(["graph", "hybrid", "mf"]):
        vals = [e2[(e2.policy == NAMED[p]) & (e2.w_neg == 10) & (e2.model == model)]["E3_AUROC_approved_vs_later_failure"].mean() for p in pols]
        ax.plot(vals, y0 + (1 - k) * 0.22, MK[model], ms=4.2, color=COL[model], mec="white", mew=0.8, label=LAB[model])
    ax.axvline(0.5, color=MUTED, lw=0.7)
    ax.set_yticks(y0); ax.set_yticklabels([])
    ax.set_xlim(0.45, 0.75)
    ax.set_xlabel("AUROC", color=MUTED)
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True); clean(ax)
    ax.set_title("c  Approved vs later\n    scientific failures", loc="left", fontsize=7.5, fontweight="bold")
    ax.legend(frameon=False, fontsize=6.0, loc="upper center", bbox_to_anchor=(0.45, -0.22), handletextpad=0.3)
    save(fig, "Figure4_E2_E3")

# ---------------- Supplementary S1: 3 x 3 grid heatmaps (graph head, per-disease AP change)
ACT = ["ignore", "negate", "mask"]
mats = []
for t in ["random", "compound"]:
    f = f"{RES}/e1_{t}_metrics.tsv"
    if os.path.exists(f):
        d = pd.read_csv(f, sep="\t"); d = d[(d.w_neg == 10) & (d.model == "graph")]
        g = d.groupby("policy").macroAP_disease.mean()
        mats.append((t, np.array([[g[f"sci={a},other={b}"] - g["sci=ignore,other=ignore"] for b in ACT] for a in ACT])))
if os.path.exists(f"{RES}/e2_metrics.tsv"):
    d = pd.read_csv(f"{RES}/e2_metrics.tsv", sep="\t"); d = d[(d.w_neg == 10) & (d.model == "graph")]
    g = d.groupby("policy").macroAP_disease.mean()
    mats.append(("e2", np.array([[g[f"sci={a},other={b}"] - g["sci=ignore,other=ignore"] for b in ACT] for a in ACT])))
if mats:
    from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("div", ["#e34948", "#f0efec", "#2a78d6"])
    fig, axes = plt.subplots(1, len(mats), figsize=(2.4 * len(mats), 2.3))
    axes = np.atleast_1d(axes)
    titles = {"random": "E1 random-edge", "compound": "E1 compound-disjoint", "e2": "E2 external approvals"}
    for ax, (t, M) in zip(axes, mats):
        lim = max(abs(M).max(), 1e-6)
        ax.imshow(M, cmap=cmap, norm=TwoSlopeNorm(0, -lim, lim))
        for i in range(3):
            for j in range(3):
                ax.text(j, i, f"{M[i, j]:+.3f}", ha="center", va="center", fontsize=6.5, color=INK)
        ax.set_xticks(range(3)); ax.set_xticklabels(ACT); ax.set_yticks(range(3)); ax.set_yticklabels(ACT)
        ax.set_xlabel("Other stopped pairs", color=MUTED); ax.set_ylabel("Scientific stops", color=MUTED)
        ax.set_title(titles[t], fontsize=7.5, fontweight="bold", loc="left")
    fig.tight_layout()
    save(fig, "FigureS1_policy_grid")

# ---------------- Supplementary S2: negative-weight sensitivity
rows = []
for t in ["random", "compound"]:
    f = f"{RES}/e1_{t}_metrics.tsv"
    if os.path.exists(f):
        d = pd.read_csv(f, sep="\t"); d = d[d.model == "graph"]
        base = d[(d.policy == NAMED["no_writeback"])].macroAP_disease.mean()
        for p in ["flat_negative", "typed"]:
            for w in [3.0, 10.0, 30.0]:
                rows.append({"set": t, "policy": p, "w": w, "delta": d[(d.policy == NAMED[p]) & (d.w_neg == w)].macroAP_disease.mean() - base})
if os.path.exists(f"{RES}/e2_metrics.tsv"):
    d = pd.read_csv(f"{RES}/e2_metrics.tsv", sep="\t"); d = d[d.model == "graph"]
    base = d[d.policy == NAMED["no_writeback"]].macroAP_disease.mean()
    for p in ["flat_negative", "typed"]:
        for w in [3.0, 10.0, 30.0]:
            rows.append({"set": "e2", "policy": p, "w": w, "delta": d[(d.policy == NAMED[p]) & (d.w_neg == w)].macroAP_disease.mean() - base})
if rows:
    r = pd.DataFrame(rows)
    sets = list(dict.fromkeys(r.set))
    fig, axes = plt.subplots(1, len(sets), figsize=(2.4 * len(sets), 2.1), sharey=False)
    axes = np.atleast_1d(axes)
    titles = {"random": "E1 random-edge", "compound": "E1 compound-disjoint", "e2": "E2 external approvals"}
    for ax, s in zip(axes, sets):
        for p, c, mk in [("flat_negative", "#2a78d6", "o"), ("typed", "#eb6834", "s")]:
            q = r[(r.set == s) & (r.policy == p)]
            ax.plot(q.w, q.delta, "-", color=c, lw=1.5); ax.plot(q.w, q.delta, mk, color=c, ms=4, mec="white", mew=0.8, label=POL_LAB[p])
        ax.axhline(0, color=MUTED, lw=0.7); ax.set_xscale("log"); ax.set_xticks([3, 10, 30]); ax.set_xticklabels(["3", "10", "30"])
        ax.set_xlabel("Weight of explicit negatives", color=MUTED); ax.set_title(titles[s], fontsize=7.5, fontweight="bold", loc="left")
        ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True); clean(ax)
    axes[0].set_ylabel("Change in per-disease AP (graph head)", color=MUTED); axes[0].legend(frameon=False, fontsize=6.3)
    fig.tight_layout()
    save(fig, "FigureS2_weight_sensitivity")

# ---------------- Supplementary S3: end-to-end ablation curves
if os.path.exists(f"{RES}/e2e_ablation.json"):
    runs = json.load(open(f"{RES}/e2e_ablation.json"))
    fig, axes = plt.subplots(1, 2, figsize=(5.4, 2.1), sharex=True)
    for rr in runs:
        h = pd.DataFrame(rr["history"])
        c, mk, lab = ("#2a78d6", "o", "Type-constrained") if rr["neg"] == "typed" else ("#eb6834", "s", "All-entity (v0.2.0)")
        for ax, col in zip(axes, ["val_AP", "test_AP"]):
            ax.plot(h.epoch, h[col], "-", color=c, lw=1.5); ax.plot(h.epoch, h[col], mk, color=c, ms=3.6, mec="white", mew=0.7, label=lab)
    for ax, t in zip(axes, ["a  Validation AP", "b  Test AP (full grid)"]):
        ax.set_title(t, loc="left", fontsize=7.5, fontweight="bold"); ax.set_xlabel("Epoch", color=MUTED)
        ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True); clean(ax)
    axes[0].legend(frameon=False, fontsize=6.3)
    fig.tight_layout()
    save(fig, "FigureS3_e2e_ablation")
print("figures written to", OUT)
