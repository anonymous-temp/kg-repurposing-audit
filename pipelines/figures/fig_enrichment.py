#!/usr/bin/env python
"""Figure 3: enrichment of approved indications among unlabelled pairs with a stopped trial (rate ratio,
95% CI on the log scale; reference: unlabelled pairs without a stop before 2015), both graphs, with the
enrichment expected from the number of stopped pairs per drug and per disease alone (degree-preserving
permutation of the stopped pairs, 1,000 permutations; mean and 95% range)."""
import json, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as S
import matplotlib.pyplot as plt

ROLE = sys.argv[1] if len(sys.argv) > 1 else "work/results_role"
PKG = sys.argv[2] if len(sys.argv) > 2 else "work/results_primekg"
OUT = sys.argv[3] if len(sys.argv) > 3 else "work/figures"
SEL = sys.argv[4] if len(sys.argv) > 4 else "work/results_selection"
NULL = {g: json.load(open(f"{SEL}/selection_{g.lower()}.json"))["null"] for g in ("Hetionet", "PrimeKG")}
S.setup()
h = json.load(open(f"{ROLE}/enrichment_hetionet.json")); p = json.load(open(f"{PKG}/descriptive.json"))["enrichment"]
# reference groups (unlabelled pairs without a stop before 2015)
href = (h["ext"] - h["wb_all"]["ext"], h["grid"] - h["wb_all"]["pairs"])
pref = (p["external_approved"] - p["stopped_ext"], p["unlabelled_pairs"] - p["stopped_unlabelled"])
ROWS = [("All stopped pairs", "wb_all"), ("Efficacy or safety stop", "wb_scientific"), ("… in the same disease concept", "wb_scientific_scoped"),
        ("… of the investigational drug", "wb_scientific_inv"), ("… both", "wb_scientific_inv_scoped"), ("Other stops only", "other")]


def counts(src, key, graph):
    if graph == "Hetionet":
        if key == "other":
            return h["wb_all"]["ext"] - h["wb_scientific"]["ext"], h["wb_all"]["pairs"] - h["wb_scientific"]["pairs"]
        return h[key]["ext"], h[key]["pairs"]
    if key == "wb_all":
        return p["stopped_ext"], p["stopped_unlabelled"]
    if key == "other":
        return p["stopped_ext"] - p["wb_scientific"]["ext"], p["stopped_unlabelled"] - p["wb_scientific"]["pairs"]
    return p[key]["ext"], p[key]["pairs"]


fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6), sharey=True)
for ax, graph, ref, letter in [(axes[0], "Hetionet", href, "A"), (axes[1], "PrimeKG", pref, "B")]:
    ys = list(range(len(ROWS)))[::-1]
    for y, (lab, key) in zip(ys, ROWS):
        a, n = counts(p, key, graph); b, m = ref
        rr = (a / n) / (b / m); se = math.sqrt(1 / a - 1 / n + 1 / b - 1 / m)
        lo, hi = rr * math.exp(-1.96 * se), rr * math.exp(1.96 * se)
        nk = NULL[graph]["wb_other" if key == "other" else key]
        assert nk["pairs"] == n and nk["observed"] == a, (graph, key, nk["pairs"], n)
        r0 = b / m; nm, nl, nh = (nk[k] / n / r0 for k in ("null_mean", "null_lo", "null_hi"))
        ax.errorbar(nm, y - 0.22, xerr=[[nm - max(nl, 1.0)], [nh - nm]], fmt="D", ms=3.0, mfc="white", mec=S.MUTED, color=S.MUTED,
                    elinewidth=0.7, capsize=1.4, capthick=0.7, label=None)
        S.errorbar(ax, rr, y + 0.1, lo, hi, S.COL[graph], S.MK[graph])
        ax.text(1.02, y, f"{a}/{n:,}", transform=ax.get_yaxis_transform(), fontsize=6, va="center", ha="left", color=S.MUTED)
    ax.text(1.02, len(ROWS) - 0.45, "Approved/pairs", transform=ax.get_yaxis_transform(), fontsize=6, va="bottom", ha="left", color=S.MUTED)
    ax.set_xscale("log"); ax.set_xlim(1, 600)
    ax.set_xticks([1, 3, 10, 30, 100, 300]); ax.set_xticklabels(["1", "3", "10", "30", "100", "300"])
    ax.axvline(1, color=S.LIGHT, lw=0.6)
    ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in ROWS])
    ax.set_xlabel("Rate ratio for approved indication (log scale)")
    S.panel(ax, letter, title=f"{graph} (reference rate {100 * ref[0] / ref[1]:.2g}%)")
    ax.tick_params(axis="y", length=0)
from matplotlib.lines import Line2D
hand = [Line2D([], [], marker=S.MK[g], ls="", color=S.COL[g], ms=3.6, label=f"Observed, {g}") for g in ("Hetionet", "PrimeKG")]
hand.append(Line2D([], [], marker="D", ls="", mfc="white", mec=S.MUTED, ms=3.0, label="Expected from the number of stopped pairs per drug and per disease"))
fig.legend(handles=hand, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.06), handletextpad=0.3, columnspacing=1.4)
fig.tight_layout(w_pad=4.5)
S.save(fig, OUT, "Figure3_enrichment")
print("saved")
