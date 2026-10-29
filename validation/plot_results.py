"""Render standalone scientific figures from completed validation outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_results(results: Path, output: Path) -> None:
    """Draw rankings, gene-statistic concordance, and sensitivity results."""

    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "svg.fonttype": "none"})
    ranks = pd.read_csv(results / "comparison" / "cell_type_rank_comparison.tsv", sep="\t")
    genes = pd.read_csv(results / "comparison" / "standard_vs_hmagma_genes.tsv.gz", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(15, 8), gridspec_kw={"width_ratios": [1, 1.4, 1.3]})
    columns = [
        "celltogwas_primary_rank",
        "celltogwas_membership_rank",
        "standard_rank",
        "hmagma_rank",
    ]
    image = axes[0].imshow(ranks[columns], aspect="auto", cmap="YlGnBu_r", vmin=1, vmax=len(ranks))
    axes[0].set_yticks(np.arange(len(ranks)), ranks.cell_type)
    axes[0].set_xticks(
        range(4), ["Weighted", "Membership", "MAGMA", "H-MAGMA"], rotation=60, ha="right"
    )
    axes[0].set_title("A  Cell-class rankings")
    fig.colorbar(image, ax=axes[0], shrink=0.5, label="Rank (1 = strongest)")
    hb = axes[1].hexbin(
        genes.standard_z, genes.hmagma_z, gridsize=60, bins="log", mincnt=1, cmap="Blues"
    )
    limit = max(genes.standard_z.max(), genes.hmagma_z.max())
    axes[1].plot([-4, limit], [-4, limit], ls="--", color="#b85c38", lw=1)
    axes[1].set(
        xlabel="Standard MAGMA gene Z",
        ylabel="H-MAGMA gene Z",
        title="B  Gene-statistic concordance",
    )
    fig.colorbar(hb, ax=axes[1], shrink=0.5, label="Gene count (log scale)")
    variants = [
        "primary_weighted",
        "sensitivity_no_mhc",
        "sensitivity_p1e5",
        "sensitivity_all_snps",
    ]
    matrix = pd.DataFrame(
        {
            name: pd.read_csv(results / "celltogwas" / f"{name}.tsv", sep="\t").set_index(
                "cell_type"
            )["z_score"]
            for name in variants
        }
    ).loc[ranks.cell_type]
    bound = max(1, np.abs(matrix.to_numpy()).max())
    im = axes[2].imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-bound, vmax=bound)
    axes[2].set_yticks(np.arange(len(ranks)), ranks.cell_type)
    axes[2].set_xticks(
        range(4), ["Primary", "No MHC", "P < 1e-5", "All SNPs"], rotation=60, ha="right"
    )
    axes[2].set_title("C  celltogwas sensitivity")
    fig.colorbar(im, ax=axes[2], shrink=0.5, label="Circular-shift Z")
    fig.suptitle("AD GWAS GCST90027158 + brainSCOPE DLPFC: external method comparison", fontsize=13)
    fig.text(
        0.5,
        0.015,
        "Rank concordance is descriptive. Gene models and LD assumptions differ; one trait and one tissue were evaluated.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(output / "external_validation.png", dpi=200)
    fig.savefig(output / "external_validation.svg")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot_results(args.results, args.output)
