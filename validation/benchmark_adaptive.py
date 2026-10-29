"""Prespecified end-to-end calibration/power simulations, not trait-P tuning.

Run: python -m validation.benchmark_adaptive --output PATH --replicates 1000
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from celltogwas import file_sha256, run_adaptive_gene_enrichment

SCENARIOS = {
    "null_independent": {"rho": 0.0, "effect": "none"},
    "null_correlated": {"rho": 0.8, "effect": "none"},
    "null_confounder_adjusted": {"rho": 0.6, "effect": "confounder"},
    "dense_weak": {"rho": 0.6, "effect": "dense", "beta": 0.18},
    "sparse_strong": {"rho": 0.6, "effect": "sparse", "beta": 3.5},
    "binary_favorable": {"rho": 0.6, "effect": "binary", "beta": 0.8},
}


def ring_noise(rng: np.random.Generator, n: int, rho: float) -> np.ndarray:
    frequency = np.arange(n // 2 + 1) * 2 * np.pi / n
    transform = np.sqrt(1 - rho**2) / (1 - rho * np.exp(-1j * frequency))
    return np.fft.irfft(np.fft.rfft(rng.normal(size=n)) * transform, n=n)


def wilson(successes: int, total: int) -> tuple[float, float]:
    p, z = successes / total, 1.959963984540054
    center = (p + z * z / (2 * total)) / (1 + z * z / total)
    radius = z * np.sqrt(p * (1 - p) / total + z * z / (4 * total**2)) / (1 + z * z / total)
    return float(center - radius), float(center + radius)


def run_benchmark(output: Path, replicates: int = 1000, permutations: int = 999) -> dict:
    if replicates < 1:
        raise ValueError("replicates must be positive")
    output.mkdir(parents=True, exist_ok=True)
    n, design_seed, outcome_seed = 512, 1301, 20261030
    rng = np.random.default_rng(design_seed)
    raw = np.column_stack([ring_noise(rng, n, 0.65) for _ in range(3)])
    weights = np.argsort(np.argsort(raw, axis=0), axis=0) / (n - 1)
    weights **= 2
    nuisance = 0.5 * (weights[:, 0] - weights[:, 0].mean()) + ring_noise(rng, n, 0.4)
    genes = pd.DataFrame(
        {
            "gene": [f"g{i}" for i in range(n)],
            "chrom": "1",
            "pos": np.arange(n),
            "A": weights[:, 0],
            "B": weights[:, 1],
            "C": weights[:, 2],
            "nuisance": nuisance,
        }
    )
    protocol = {
        "scenarios": SCENARIOS,
        "n_genes": n,
        "replicates_per_scenario": replicates,
        "permutations": permutations,
        "design_seed": design_seed,
        "outcome_seed": outcome_seed,
        "powers": [1, 2, 4],
        "thresholds": [None, 2, 4],
        "binary_fraction": 0.1,
        "alpha": 0.05,
        "n_cell_types": 3,
        "power_target": "A, per-cell alpha=0.05",
        "baseline": "same data, covariates, shifts and empirical calibration; one fixed kernel",
        "module_sha256_before_run": file_sha256(
            Path(__file__).parents[1] / "celltogwas/adaptive.py"
        ),
        "benchmark_sha256_before_run": file_sha256(__file__),
        "limitations": "Toy stationary Gaussian gene-order dependence, not realistic GWAS LD or donor uncertainty.",
    }
    (output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n", encoding="utf-8")
    rows, started = [], time.perf_counter()
    scenario_seeds = np.random.SeedSequence(outcome_seed).spawn(len(SCENARIOS))
    for (name, config), seed in zip(SCENARIOS.items(), scenario_seeds, strict=False):
        rng = np.random.default_rng(seed)
        for iteration in range(replicates):
            z = ring_noise(rng, n, config["rho"])
            effect = config["effect"]
            covariates = []
            if effect == "confounder":
                z += 2 * nuisance
                covariates = ["nuisance"]
            elif effect == "dense":
                z += config["beta"] * (weights[:, 0] - weights[:, 0].mean()) / weights[:, 0].std()
            elif effect == "sparse":
                z += config["beta"] * (weights[:, 0] >= np.quantile(weights[:, 0], 0.98))
            elif effect == "binary":
                z += config["beta"] * (weights[:, 0] >= np.quantile(weights[:, 0], 0.9))
            genes["z"] = z
            result = run_adaptive_gene_enrichment(
                genes,
                cell_types=["A", "B", "C"],
                covariates=covariates,
                n_permutations=permutations,
                random_state=iteration + 9000,
            )
            target = result.table.set_index("cell_type").loc["A"]
            parts = result.components.query("cell_type == 'A'").set_index("kernel")
            rows.append(
                {
                    "scenario": name,
                    "iteration": iteration,
                    "adaptive": target.p_adaptive,
                    "continuous_dense": parts.loc["power1:dense", "p_component_uncorrected"],
                    "binary_dense": parts.loc["binary0.1:dense", "p_component_uncorrected"],
                    "family_holm_any": bool((result.table.p_holm <= 0.05).any()),
                }
            )
        print(
            f"{name}: {replicates} replicates, elapsed={time.perf_counter() - started:.1f}s",
            flush=True,
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(output / "replicates.tsv", sep="\t", index=False, lineterminator="\n")
    summaries = []
    for name, group in frame.groupby("scenario", sort=False):
        for method in ["adaptive", "continuous_dense", "binary_dense", "family_holm_any"]:
            rejected = (
                group[method].to_numpy()
                if method == "family_holm_any"
                else group[method].to_numpy() <= 0.05
            )
            count = int(rejected.sum())
            lo, hi = wilson(count, replicates)
            summaries.append(
                {
                    "scenario": name,
                    "method": method,
                    "rejections": count,
                    "n": replicates,
                    "rate": count / replicates,
                    "ci95_low": lo,
                    "ci95_high": hi,
                }
            )
    summary = pd.DataFrame(summaries)
    summary.to_csv(output / "summary.tsv", sep="\t", index=False, lineterminator="\n")
    null = summary[summary.scenario.str.startswith("null_") & (summary.method == "adaptive")]
    verdict = {
        "elapsed_seconds": time.perf_counter() - started,
        "null_gate": "Each adaptive null Wilson 95% upper bound <= 0.075 (diagnostic, not a proof).",
        "null_gate_passed": bool((null.ci95_high <= 0.075).all()),
        "summary": summaries,
    }
    (output / "summary.json").write_text(json.dumps(verdict, indent=2) + "\n", encoding="utf-8")
    return verdict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--permutations", type=int, default=999)
    args = parser.parse_args()
    run_benchmark(args.output, args.replicates, args.permutations)
