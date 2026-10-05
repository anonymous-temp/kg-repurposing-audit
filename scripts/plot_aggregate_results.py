"""Render reference plots from the published aggregate tables, without databases."""

import argparse
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path(__file__).resolve().parents[1] / "paper_results")
    parser.add_argument("--output", type=Path, default=Path("outputs/reference_figures"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    models = {
        "degree": ("Degree reference", "#444444"),
        "complex": ("ComplEx", "#0072B2"),
        "distmult": ("DistMult", "#D89000"),
        "rotate": ("RotatE", "#009E73"),
        "rgcn": ("R-GCN", "#CC79A7"),
    }
    plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
    for resource, tasks, samplers, labels in [
        (
            "hetionet",
            ["random", "coldstart", "scaffold"],
            ["neg_random", "neg_degmatch"],
            ["Uniform", "Disease frequency"],
        ),
        ("primekg", ["random", "compound_disjoint"], ["uniform", "matched"], ["Uniform", "Joint degree bins"]),
    ]:
        frame = pd.read_csv(args.results / f"{resource}_table.csv")
        fig, axs = plt.subplots(1, len(tasks), figsize=(7.2, 3.2), sharey=True)
        for ax, task in zip(axs, tasks):
            for model, (label, color) in models.items():
                cell = frame[(frame.task == task) & (frame.model == model)]
                if cell.empty:
                    continue
                cell = cell.set_index("sampler").loc[samplers]
                ax.errorbar(
                    [0, 1], cell.AP_mean, yerr=cell.AP_sd, label=label, color=color, marker="o", capsize=2, markersize=3
                )
            ax.axhline(1 / 21, color=".6", linestyle=":", linewidth=0.8)
            ax.set_title(task.replace("_", " "))
            ax.set_xticks([0, 1], labels, fontsize=8)
            ax.set_xlim(-0.15, 1.15)
        axs[0].set_ylabel("Average precision")
        fig.legend(*axs[0].get_legend_handles_labels(), loc="lower center", ncol=3, frameon=False, fontsize=8)
        fig.tight_layout(rect=[0, 0.16, 1, 1])
        for ext in ["png", "pdf"]:
            fig.savefig(args.output / f"{resource}_reference.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)
    print(f"Reference plots written to {args.output}. Error bars are seed SD, not confidence intervals.")


if __name__ == "__main__":
    main()
