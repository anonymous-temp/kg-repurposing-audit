"""Shared journal style for all figures: sans-serif 7 pt, outward ticks, no top/right spines,
bold upper-case panel letters, colour-blind-safe palette (validated), vector PDF + 600 dpi PNG."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, MUTED, LIGHT = "#000000", "#4d4d4d", "#d9d9d9"
COL = {"graph": "#2a78d6", "hybrid": "#eb6834", "mf": "#1baf7a", "degree": "#eda100",
       "Hetionet": "#2a78d6", "PrimeKG": "#eb6834"}
MK = {"graph": "o", "hybrid": "s", "mf": "D", "degree": "^", "Hetionet": "o", "PrimeKG": "s"}
LAB = {"graph": "Graph head", "hybrid": "Hybrid", "mf": "Label-only MF", "degree": "Degree reference"}


def setup():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK, "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.dpi": 600, "figure.dpi": 150, "lines.linewidth": 0.9})


def panel(ax, letter, x=-0.02, y=1.04, title=None):
    """Bold upper-case panel letter at the top left of the axes, optional short plain title after it."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right")
    if title:
        ax.text(x + 0.015, y, title, transform=ax.transAxes, fontsize=7, va="bottom", ha="left")


def errorbar(ax, x, y, lo, hi, color, marker, label=None, horizontal=True, ms=3.6):
    kw = dict(fmt=marker, ms=ms, color=color, mfc=color, mec=color, elinewidth=0.8, capsize=1.6, capthick=0.8, label=label)
    if horizontal:
        ax.errorbar(x, y, xerr=[[x - lo], [hi - x]], **kw)
    else:
        ax.errorbar(x, y, yerr=[[y - lo], [hi - y]], **kw)


def save(fig, out, name):
    fig.savefig(f"{out}/{name}.pdf", bbox_inches="tight")
    fig.savefig(f"{out}/{name}.png", dpi=600, bbox_inches="tight")
