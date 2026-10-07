#!/usr/bin/env python
"""FigS_logistic_models: all logistic models of Table S20a (supp_S20_S21.md), odds ratios with 95% disease-cluster bootstrap intervals
(log scale), Hetionet and PrimeKG.
A  odds ratios of the stop terms (scientific stop; any stop), grouped by term
B  odds ratios of the testing-intensity covariates (per unit of the log-transformed counts)
Every row of the table is plotted exactly once (asserted); values are cross-checked against selection_{graph}.json."""
import json, re
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import b_common as C

S = C.S
S.setup()
G = ("Hetionet", "PrimeKG")
_, rows = C.read_table(f"{C.MS}/supp_S20_S21.md", r"\*\*Table S20a\.\*\*|Table S20a\.")
T = {}
for r in rows:
    d = {}
    for g, cell in zip(G, r[2:4]):
        m = re.match(r"([\d.]+) \(([\d.]+)–([\d.]+)\)", cell)
        d[g] = tuple(float(x) for x in m.groups())
    T[(r[0], r[1])] = d
assert len(T) == 19

# table row (model, term) -> JSON key (model, term)
JK = {("Crude", "Scientific stop"): ("crude", "sci"), ("Adjusted for log number of trials", "log number of trials"): ("adjusted", "ln"),
      ("Adjusted for log number of trials", "Scientific stop"): ("adjusted", "sci"),
      ("Adjusted, scientific and any stop", "log number of trials"): ("adjusted_any_stop", "ln"),
      ("Adjusted, scientific and any stop", "Scientific stop"): ("adjusted_any_stop", "sci"),
      ("Adjusted, scientific and any stop", "Any stop"): ("adjusted_any_stop", "stop"),
      ("Adjusted, plus trial counts of drug and disease", "log number of trials"): ("adjusted_popularity", "ln"),
      ("Adjusted, plus trial counts of drug and disease", "log trials of the drug"): ("adjusted_popularity", "lnc"),
      ("Adjusted, plus trial counts of drug and disease", "log trials of the disease"): ("adjusted_popularity", "lnd"),
      ("Adjusted, plus trial counts of drug and disease", "Scientific stop"): ("adjusted_popularity", "sci"),
      ("Intensity from phase 1–2 trials only", "log(1 + phase 1–2 trials)"): ("adjusted_phase12_intensity", "ln12"),
      ("Intensity from phase 1–2 trials only", "Scientific stop"): ("adjusted_phase12_intensity", "sci"),
      ("Adjusted, cancers only", "log number of trials"): ("adjusted_oncology", "ln"), ("Adjusted, cancers only", "Scientific stop"): ("adjusted_oncology", "sci"),
      ("Adjusted, other diseases only", "log number of trials"): ("adjusted_non_oncology", "ln"),
      ("Adjusted, other diseases only", "Scientific stop"): ("adjusted_non_oncology", "sci"),
      ("Crude, any stop", "Any stop"): ("crude_any_stop", "stop"), ("Adjusted, any stop", "log number of trials"): ("adjusted_any_stop_only", "ln"),
      ("Adjusted, any stop", "Any stop"): ("adjusted_any_stop_only", "stop")}
assert set(JK) == set(T)
J = {g: json.load(open(f"{C.SEL}/selection_{g.lower()}.json"))["intensity"]["models"] for g in G}
for k, v in T.items():                                           # table values equal the JSON (rounded to 2 decimals)
    mk, tk = JK[k]
    for g in G:
        j = J[g][mk]["terms"].get(tk)
        assert j is not None, (g, k)
        assert all(abs(a - b) < 0.0051 for a, b in zip(v[g], (j["OR"], j["lo"], j["hi"]))), (g, k, v[g], j)

MODEL = {"Crude": "Crude", "Adjusted for log number of trials": "Adjusted for log number of trials",
         "Adjusted, scientific and any stop": "Scientific and any stop in one model",
         "Adjusted, plus trial counts of drug and disease": "+ trial counts of drug and disease",
         "Intensity from phase 1–2 trials only": "Intensity from phase 1–2 trials only",
         "Adjusted, cancers only": "Cancers only", "Adjusted, other diseases only": "Other diseases only",
         "Crude, any stop": "Crude", "Adjusted, any stop": "Adjusted for log number of trials"}
A_ROWS = [("Scientific stop", [("Crude", "Scientific stop"), ("Adjusted for log number of trials", "Scientific stop"),
                               ("Adjusted, scientific and any stop", "Scientific stop"),
                               ("Adjusted, plus trial counts of drug and disease", "Scientific stop"),
                               ("Intensity from phase 1–2 trials only", "Scientific stop"),
                               ("Adjusted, cancers only", "Scientific stop"), ("Adjusted, other diseases only", "Scientific stop")]),
          ("Any stop", [("Crude, any stop", "Any stop"), ("Adjusted, any stop", "Any stop"),
                        ("Adjusted, scientific and any stop", "Any stop")])]
B_ROWS = [("Log number of trials", [("Adjusted for log number of trials", "log number of trials"),
                                   ("Adjusted, scientific and any stop", "log number of trials"),
                                   ("Adjusted, plus trial counts of drug and disease", "log number of trials"),
                                   ("Adjusted, cancers only", "log number of trials"), ("Adjusted, other diseases only", "log number of trials"),
                                   ("Adjusted, any stop", "log number of trials")]),
          ("Other intensity terms", [("Adjusted, plus trial counts of drug and disease", "log trials of the drug"),
                                     ("Adjusted, plus trial counts of drug and disease", "log trials of the disease"),
                                     ("Intensity from phase 1–2 trials only", "log(1 + phase 1–2 trials)")])]
used = [k for _, rs in A_ROWS + B_ROWS for k in rs]
assert sorted(used) == sorted(T) and len(set(used)) == 19
MAIN = {("Crude", "Scientific stop"), ("Adjusted for log number of trials", "Scientific stop"),
        ("Adjusted, plus trial counts of drug and disease", "Scientific stop"), ("Intensity from phase 1–2 trials only", "Scientific stop"),
        ("Adjusted, cancers only", "Scientific stop"), ("Adjusted, other diseases only", "Scientific stop")}   # rows of main-text Fig. 4D
OTHER_LAB = {"log trials of the drug": "Log trials of the drug\n(drug and disease model)",
             "log trials of the disease": "Log trials of the disease\n(drug and disease model)",
             "log(1 + phase 1–2 trials)": "Log(1 + phase 1–2 trials)\n(phase 1–2 model)"}
MODEL_B = dict(MODEL, **{"Adjusted for log number of trials": "Scientific-stop model", "Adjusted, any stop": "Any-stop model"})


def layout(groups):
    """y positions (top to bottom) with a header row per group."""
    out, y = [], 0.0
    for head, rs in groups:
        out.append(("head", head, y)); y += 1
        for k in rs:
            out.append(("row", k, y)); y += 1
        y += 0.5
    return out, y - 0.5


def draw(ax, groups, labfun, letter, title):
    items, total = layout(groups)
    for kind, v, y in items:
        yy = total - y
        if kind == "head":
            ax.text(-0.02, yy, v, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.8, fontweight="bold")
            continue
        for g, off in (("Hetionet", 0.21), ("PrimeKG", -0.21)):
            o, lo, hi = T[v][g]
            S.errorbar(ax, o, yy + off, lo, hi, S.COL[g], S.MK[g], ms=3.2)
            ax.text(1.02, yy + off, f"{o:.2f} ({lo:.2f}–{hi:.2f})", transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=5.8,
                    color=S.COL[g] if False else S.INK)
        ax.text(-0.02, yy, labfun(v), transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.5)
    ax.axvline(1, color=S.MUTED, lw=0.6, ls=":")
    ax.set_xscale("log"); ax.set_xlim(0.2, 6)
    ax.set_xticks([0.2, 0.5, 1, 2, 5]); ax.set_xticklabels(["0.2", "0.5", "1", "2", "5"])
    ax.set_ylim(0.1, total + 0.7); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    S.panel(ax, letter, x=-0.62, y=1.0 + 0.03, title=title)
    ymin, ymax = ax.get_ylim()
    return {v: (total - y - ymin) / (ymax - ymin) for kind, v, y in items if kind == "head"}


fig = plt.figure(figsize=(6.7, 7.1))
nA = sum(1 + len(r) for _, r in A_ROWS) + 0.5 * (len(A_ROWS) - 1)
nB = sum(1 + len(r) for _, r in B_ROWS) + 0.5 * (len(B_ROWS) - 1)
gs = fig.add_gridspec(2, 1, height_ratios=[nA + 1.3, nB + 1.3], hspace=0.28, left=0.40, right=0.80, top=0.965, bottom=0.075)
axA, axB = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
headA = draw(axA, A_ROWS, lambda k: MODEL[k[0]], "A", "Odds ratio of a stop for being an indication")
draw(axB, B_ROWS, lambda k: MODEL_B[k[0]] if k[1] == "log number of trials" else OTHER_LAB[k[1]], "B", "Odds ratio of the testing-intensity terms")
axA.set_xlabel("Odds ratio (log scale)"); axB.set_xlabel("Odds ratio per unit of the log count (log scale)")
for ax in (axA, axB):
    ax.text(1.02, 1.0, "OR (95% CI)", transform=ax.transAxes, fontsize=6, color=S.MUTED, va="bottom", ha="left")
hand = [Line2D([], [], marker=S.MK[g], ls="", color=S.COL[g], ms=3.6, label=g) for g in G]
axA.legend(handles=hand, loc="center right", bbox_to_anchor=(0.995, headA["Any stop"]), handletextpad=0.2, fontsize=6.3, ncol=2, columnspacing=1.0)
S.save(fig, C.OUT, "FigS_logistic_models")
print("saved; rows:", len(used), "of which in main Fig. 4D:", len(MAIN))
