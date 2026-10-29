"""Selection-adjusted, multiscale competitive gene-program enrichment.

This is a gene-order circular randomization test, NOT an LD covariance model.
Validity requires cyclic exchangeability of the residual genomic processes.
See docs/ADAPTIVE.md before interpreting a result as biological evidence.
"""

from __future__ import annotations

import numbers
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .enrichment import benjamini_hochberg
from .interval import normalize_chromosome


@dataclass
class AdaptiveResult:
    """Calibrated cell tests and explicitly uncorrected component diagnostics."""

    table: pd.DataFrame
    components: pd.DataFrame
    diagnostics: dict[str, object]
    null_statistics: np.ndarray | None = None


def _integer(value: int, name: str, minimum: int) -> None:
    if isinstance(value, bool) or not isinstance(value, numbers.Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}.")


def _upper_p(values: np.ndarray) -> np.ndarray:
    """Conservative upper-tail ranks, treating round-off ties symmetrically."""
    output = np.empty_like(values, dtype=float)
    for column in range(values.shape[1]):
        vector = values[:, column]
        tolerance = 1e-10 * max(1.0, float(np.max(np.abs(vector))))
        output[:, column] = (
            len(vector) - np.searchsorted(np.sort(vector), vector - tolerance, side="left")
        ) / len(vector)
    return output


def calibrate_adaptive(statistics: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return adaptive and component P values for EVERY exchangeable row.

    Input is (observed + randomizations, kernels, cell types); row zero is
    observed. Both ranking stages include it. Re-ranking every row symmetrically
    accounts for choosing the smallest component P, unlike reporting min(P).
    Exchangeability must come from the caller's null, not from this function.
    """
    values = np.asarray(statistics, dtype=float)
    if values.ndim != 3 or min(values.shape) < 1 or values.shape[0] < 2:
        raise ValueError("statistics must have shape (at least 2 rows, kernels, cell types).")
    if not np.isfinite(values).all():
        raise ValueError("statistics must be finite.")
    component = _upper_p(values.reshape(len(values), -1)).reshape(values.shape)
    minimum = component.min(axis=1)
    adaptive = _upper_p(-minimum)
    return adaptive, component


def _holm(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values, kind="stable")
    adjusted = np.maximum.accumulate(p_values[order] * np.arange(len(order), 0, -1))
    output = np.empty_like(p_values)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def _residualize(values: np.ndarray, covariates: np.ndarray) -> tuple[np.ndarray, int]:
    centered = covariates - covariates.mean(axis=0, keepdims=True)
    scales = np.sqrt(np.sum(centered**2, axis=0))
    keep = scales > np.finfo(float).eps
    design = np.column_stack((np.ones(len(values)), centered[:, keep] / scales[keep]))
    u, singular, _ = np.linalg.svd(design, full_matrices=False)
    rank = int(np.sum(singular > singular[0] * max(design.shape) * np.finfo(float).eps))
    basis = u[:, :rank]
    return values - basis @ (basis.T @ values), rank


def run_adaptive_gene_enrichment(
    genes: pd.DataFrame,
    *,
    cell_types: Sequence[str],
    covariates: Sequence[str] = (),
    powers: Sequence[float] = (1.0, 2.0, 4.0),
    thresholds: Sequence[float | None] = (None, 2.0, 4.0),
    binary_fraction: float | None = 0.1,
    n_permutations: int = 9_999,
    random_state: int = 13,
    z_cap: float = 10.0,
    return_null_statistics: bool = False,
) -> AdaptiveResult:
    """Test continuous cell weights against dense and sparse gene evidence.

    Required columns: gene, chrom, pos (non-negative integer), z (gene
    association strength, e.g. MAGMA ZSTAT), and each requested cell weight.
    Weights must be non-negative; all background genes must be retained,
    including zero weights. Covariates must be prespecified and finite.

    Each chromosome is sorted by position. Evidence and weight transforms are
    residualized on an intercept and the covariates by SVD. Independent uniform
    circular shifts (including zero) are shared by ALL kernels and cell types.
    FFT cross-correlations retain negative residuals, unlike non-negative SNP
    enrichment. Inference is conditional on this residual shift-null model.
    """
    _integer(n_permutations, "n_permutations", 19)
    _integer(random_state, "random_state", 0)
    cells, covars = list(cell_types), list(covariates)
    if not cells or len(set(cells)) != len(cells):
        raise ValueError("cell_types must be non-empty and unique.")
    if len(set(covars)) != len(covars) or set(covars) & set(cells):
        raise ValueError("covariates must be unique and distinct from cell_types.")
    reserved = {"gene", "chrom", "pos", "z"}
    if reserved & (set(cells) | set(covars)):
        raise ValueError("Weights/covariates cannot use reserved input column names.")
    if not isinstance(genes, pd.DataFrame) or genes.empty or genes.columns.has_duplicates:
        raise ValueError("genes must be a non-empty DataFrame with unique columns.")
    required = ["gene", "chrom", "pos", "z", *cells, *covars]
    if not set(required).issubset(genes.columns):
        raise ValueError(f"Missing input columns: {sorted(set(required) - set(genes.columns))}")
    frame = genes[required].copy()
    if frame["gene"].isna().any() or frame["gene"].astype(str).str.strip().eq("").any():
        raise ValueError("Gene identifiers must be present.")
    frame["gene"] = frame["gene"].astype(str)
    if frame["gene"].duplicated().any():
        raise ValueError(
            "Gene identifiers must be unique; do not count aliases as independent genes."
        )
    if frame["chrom"].isna().any():
        raise ValueError("Chromosomes must be present.")
    frame["chrom"] = frame["chrom"].map(normalize_chromosome)
    if frame["chrom"].eq("").any():
        raise ValueError("Chromosomes must be present.")
    numeric = ["pos", "z", *cells, *covars]
    frame[numeric] = frame[numeric].astype(float)
    if not np.isfinite(frame[numeric].to_numpy()).all():
        raise ValueError("Positions, Z statistics, weights and covariates must be finite.")
    if (frame["pos"] < 0).any() or (frame["pos"] % 1 != 0).any():
        raise ValueError("Positions must be non-negative integers.")
    if (frame[cells] < 0).any().any():
        raise ValueError("Cell weights must be non-negative.")
    powers = tuple(float(p) for p in powers)
    if (
        not powers
        or len(set(powers)) != len(powers)
        or any(not np.isfinite(p) or p <= 0 for p in powers)
    ):
        raise ValueError("powers must be distinct, positive and finite.")
    thresholds = tuple(None if t is None else float(t) for t in thresholds)
    if (
        not thresholds
        or len(set(thresholds)) != len(thresholds)
        or any(t is not None and (not np.isfinite(t) or t < 0) for t in thresholds)
    ):
        raise ValueError("thresholds must be distinct non-negative finite values or None (dense).")
    if not np.isfinite(z_cap) or z_cap <= 0:
        raise ValueError("z_cap must be positive and finite.")
    if any(t is not None and t >= z_cap for t in thresholds):
        raise ValueError("Sparse thresholds must be below z_cap.")
    if binary_fraction is not None and (
        not np.isfinite(binary_fraction) or not 0 < binary_fraction < 1
    ):
        raise ValueError("binary_fraction must be in (0, 1) or None.")

    frame = frame.sort_values(["chrom", "pos", "gene"], kind="stable").reset_index(drop=True)
    weights = frame[cells].to_numpy()
    maxima = weights.max(axis=0)
    weights = np.divide(weights, maxima, out=np.zeros_like(weights), where=maxima > 0)
    weight_arrays = [weights**p for p in powers]
    weight_names = [f"power{p:g}" for p in powers]
    if binary_fraction is not None:
        # Include all cutoff ties; a constant program must remain constant.
        cutoff = np.quantile(weights, 1.0 - binary_fraction, axis=0)
        weight_arrays.append(((weights >= cutoff) & (weights > 0)).astype(float))
        weight_names.append(f"binary{binary_fraction:g}")
    names = [
        f"{w}:{'dense' if t is None else f'z>{t:g}'}" for w in weight_names for t in thresholds
    ]
    statistics = np.zeros((n_permutations + 1, len(names), len(cells)))
    chromosome_diagnostics = []
    groups = list(frame.groupby("chrom", sort=True))
    seeds = np.random.SeedSequence(random_state).spawn(len(groups))
    for (chromosome, group), seed in zip(groups, seeds, strict=False):
        rows = group.index.to_numpy()
        if len(rows) < 8:
            raise ValueError(f"Chromosome {chromosome} has fewer than 8 genes; filter explicitly.")
        nuisance = group[covars].to_numpy()
        z, _ = _residualize(group[["z"]].to_numpy(), nuisance)
        z = np.clip(z[:, 0], -z_cap, z_cap)
        evidence = np.column_stack([z if t is None else np.maximum(z - t, 0) for t in thresholds])
        evidence, rank = _residualize(evidence, nuisance)
        if len(rows) - rank < 4:
            raise ValueError(
                f"Chromosome {chromosome} has insufficient residual degrees of freedom."
            )
        shifts = np.r_[0, np.random.default_rng(seed).integers(0, len(rows), n_permutations)]
        evidence_fft = np.conj(np.fft.rfft(evidence, axis=0))
        for wi, array in enumerate(weight_arrays):
            residual_weights, _ = _residualize(array[rows], nuisance)
            weight_fft = np.fft.rfft(residual_weights, axis=0)
            for ti in range(len(thresholds)):
                cross = np.fft.irfft(weight_fft * evidence_fft[:, ti, None], n=len(rows), axis=0)
                statistics[:, wi * len(thresholds) + ti, :] += cross[shifts]
        chromosome_diagnostics.append(
            {"chrom": str(chromosome), "n_genes": len(rows), "nuisance_rank": rank}
        )

    adaptive, component = calibrate_adaptive(statistics)
    p_values = adaptive[0]
    best = component[0].argmin(axis=0)
    q_bh = benjamini_hochberg(p_values)
    table = pd.DataFrame(
        {
            "cell_type": cells,
            "p_adaptive": p_values,
            "p_holm": _holm(p_values),
            "q_by": np.minimum(q_bh * np.sum(1 / np.arange(1, len(cells) + 1)), 1.0),
            "q_bh_exploratory": q_bh,
            "selected_kernel": [names[i] for i in best],
            "selected_p_uncorrected": component[0, best, np.arange(len(cells))],
        }
    ).sort_values(["p_adaptive", "cell_type"], kind="stable", ignore_index=True)
    components = pd.DataFrame(
        [
            {
                "cell_type": cell,
                "kernel": name,
                "statistic": statistics[0, k, c],
                "p_component_uncorrected": component[0, k, c],
            }
            for c, cell in enumerate(cells)
            for k, name in enumerate(names)
        ]
    )
    return AdaptiveResult(
        table=table,
        components=components,
        diagnostics={
            "method": "adaptive-gene-circular-minP-v1",
            "n_genes": len(frame),
            "n_cell_types": len(cells),
            "n_kernels": len(names),
            "n_permutations": n_permutations,
            "random_state": random_state,
            "powers": list(powers),
            "thresholds": list(thresholds),
            "binary_fraction": binary_fraction,
            "z_cap": z_cap,
            "covariates": covars,
            "chromosomes": chromosome_diagnostics,
            "minimum_p": 1.0 / (n_permutations + 1),
            "null_assumption": "Within-chromosome cyclic exchangeability of residual processes; not arbitrary LD control.",
            "interpretation": "Competitive enrichment, not gene causality or independent replication.",
        },
        null_statistics=statistics if return_null_statistics else None,
    )
