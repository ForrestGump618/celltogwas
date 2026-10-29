"""Plot measured toy calibration/power and real-data disagreement, not discoveries."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot(root: Path, output: Path):
    summary = pd.read_csv(root / "adaptive-benchmark/summary.tsv", sep="\t")
    stable = pd.read_csv(root / "adaptive-stability/combined.tsv", sep="\t")
    magma = pd.read_csv(root / "adaptive-magma/combined.tsv", sep="\t")
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 9.4), layout="constrained")
    colors = ["#155e63", "#bf6b31", "#727879"]
    methods = ["adaptive", "continuous_dense", "binary_dense"]
    labels = ["Adaptive (12 kernels)", "Fixed continuous dense", "Fixed binary dense"]
    nulls = summary[summary.scenario.str.startswith("null_") & (summary.method == "adaptive")]
    ax = axes[0, 0]
    ax.errorbar(
        np.arange(3),
        nulls.rate,
        yerr=[nulls.rate - nulls.ci95_low, nulls.ci95_high - nulls.rate],
        fmt="o",
        color=colors[0],
        capsize=5,
    )
    ax.axhline(0.05, color="black", linestyle="--", linewidth=1)
    ax.set(
        xticks=np.arange(3),
        xticklabels=["Independent", "Correlated", "Nuisance-adjusted"],
        ylim=(0, 0.09),
        ylabel="Type-I rejection probability",
        title="A  Calibration: 1,000 replicates per null",
    )
    ax.text(
        0.02,
        0.94,
        "Nominal alpha = 0.05; Wilson 95% intervals",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
    )
    ax = axes[0, 1]
    scenarios = ["dense_weak", "sparse_strong", "binary_favorable"]
    for j, (method, label, color) in enumerate(zip(methods, labels, colors, strict=False)):
        group = summary.set_index(["scenario", "method"]).loc[[(s, method) for s in scenarios]]
        x = np.arange(3) + (j - 1) * 0.25
        ax.bar(x, group.rate, width=0.23, color=color, label=label)
        ax.errorbar(
            x,
            group.rate,
            yerr=[group.rate - group.ci95_low, group.ci95_high - group.rate],
            fmt="none",
            color="black",
            capsize=2,
            linewidth=0.7,
        )
    ax.set(
        xticks=np.arange(3),
        xticklabels=["Dense weak", "Sparse strong", "Binary-favorable"],
        ylim=(0, 1.2),
        ylabel="Per-cell power at alpha = 0.05",
        title="B  Adaptation improves robustness, not every optimum",
    )
    ax.legend(loc="upper left", fontsize=8, frameon=False, ncol=1)
    ax = axes[1, 0]
    cells = ["Micro.PVM", "Immune"]
    for j, method in enumerate(["standard", "hmagma"]):
        group = stable[stable.method == method].set_index("cell_type").loc[cells]
        ax.bar(
            np.arange(2) + (j - 0.5) * 0.32,
            group.q_by_global48,
            width=0.3,
            color=colors[j],
            label=method,
        )
        for i, value in enumerate(group.q_by_global48):
            ax.text(i + (j - 0.5) * 0.32, value + 0.002, f"{value:.4f}", ha="center", fontsize=9)
    ax.axhline(0.05, color="black", linestyle="--", linewidth=1)
    ax.set(
        xticks=np.arange(2),
        xticklabels=cells,
        ylim=(0, 0.09),
        ylabel="BY q across 48 primary tests",
        title="C  Fixed 99,999-resample stability check",
    )
    ax.legend(frameon=False)
    ax = axes[1, 1]
    for _j, (method, kind) in enumerate(
        [(m, k) for m in ["standard", "hmagma"] for k in ["continuous", "binary"]]
    ):
        group = (
            magma[(magma.method == method) & (magma.kind == kind)].set_index("cell_type").loc[cells]
        )
        ax.plot(np.arange(2), group.q_by_global96, "o-", label=f"{method}, {kind}", linewidth=1)
    ax.axhline(0.05, color="black", linestyle="--", linewidth=1)
    ax.set(
        xticks=np.arange(2),
        xticklabels=cells,
        ylim=(0, 0.3),
        ylabel="BY q across 96 property tests",
        title="D  Actual LD-aware MAGMA property regression",
    )
    ax.legend(fontsize=8, frameon=False)
    fig.suptitle(
        "celltogwas 1.1 alpha | measured power and honest evidence boundaries", fontsize=15
    )
    fig.supxlabel(
        "Toy stationary simulations are not general LD validation. Same AD GWAS reanalysis is not independent replication.",
        fontsize=9,
    )
    fig.savefig(output / "adaptive_validation.png", dpi=180)
    fig.savefig(output / "adaptive_validation.svg")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.root, args.output)
