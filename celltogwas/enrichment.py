"""Spatially aware empirical enrichment tests.

The test asks whether trait-associated SNPs receive higher cell-type link
scores than expected when the cell-type annotation is circularly shifted
within each chromosome. Circular shifts preserve the local ordering and broad
clustering of genomic annotations while avoiding a naive independent-SNP null.
"""

from __future__ import annotations

import math
import numbers
import warnings
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .exceptions import EnrichmentError
from .interval import normalize_chromosome
from .io import resolve_column


@dataclass
class EnrichmentResult:
    """Enrichment table plus reproducibility diagnostics."""

    table: pd.DataFrame
    diagnostics: dict[str, object]
    null_scores: np.ndarray | None = None


def benjamini_hochberg(p_values: Sequence[float] | np.ndarray) -> np.ndarray:
    """Return Benjamini-Hochberg adjusted p-values."""

    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1:
        raise ValueError("p_values must be one-dimensional.")
    if len(values) == 0:
        return values.copy()
    if np.any(~np.isfinite(values)) or np.any((values < 0) | (values > 1)):
        raise ValueError("p-values must be finite and lie in [0, 1].")
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def circular_correlation_matrix(
    weights: np.ndarray,
    scores: np.ndarray,
) -> np.ndarray:
    """Compute all circular weighted dot products efficiently with an FFT.

    Entry ``k, c`` is the null numerator for a shift between annotation scores
    and evidence weights on chromosome ``c``:
    ``dot(weights, roll(scores[:, c], -k))``. Entry ``0`` equals
    ``weights @ scores``.
    """

    evidence = np.asarray(weights, dtype=float)
    link_scores = np.asarray(scores, dtype=float)
    if evidence.ndim != 1:
        raise ValueError("weights must be one-dimensional.")
    if link_scores.ndim == 1:
        link_scores = link_scores[:, None]
    if link_scores.ndim != 2 or link_scores.shape[0] != len(evidence):
        raise ValueError("scores must have shape (n_variants, n_cell_types).")
    if len(evidence) == 0:
        return np.empty((0, link_scores.shape[1]), dtype=float)
    transform = np.fft.rfft(link_scores, axis=0) * np.conj(np.fft.rfft(evidence)[:, None])
    correlation = np.fft.irfft(transform, n=len(evidence), axis=0)
    # Cross-correlations of non-negative sequences are non-negative; clip FFT
    # round-off so it cannot leak into null means or fold changes.
    return np.maximum(np.asarray(correlation, dtype=float), 0.0)


def _robust_standardized_magnitudes(
    values: np.ndarray,
    *,
    floor: float = 0.0,
) -> tuple[np.ndarray, float, float]:
    """Standardize absolute effects using the median and scaled MAD."""

    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) == 0:
        raise EnrichmentError("Effect values must be a non-empty one-dimensional array.")
    if not np.isfinite(values).all():
        raise EnrichmentError("Effect values must be finite.")
    magnitudes = np.abs(values)
    center = float(np.median(magnitudes))
    deviations = np.abs(magnitudes - center)
    scale = float(1.4826 * np.median(deviations))
    if scale <= np.finfo(float).eps:
        scale = float(np.mean(deviations))
    if scale <= np.finfo(float).eps:
        evidence = np.zeros_like(magnitudes)
    else:
        evidence = np.maximum((magnitudes - center) / scale, 0.0)
    if floor:
        evidence = np.maximum(evidence - floor, 0.0)
    return evidence, center, scale


def run_enrichment(
    annotated_gwas: pd.DataFrame,
    *,
    cell_types: Sequence[str],
    p_column: str | None = None,
    score_column: str | None = None,
    score_mode: str = "pvalue",
    n_permutations: int = 10_000,
    random_state: int = 13,
    min_shift_fraction: float = 0.01,
    p_upper: float | None = None,
    weight_floor: float = 0.0,
    min_p: float = 1e-300,
    max_evidence: float = 300.0,
    min_chromosome_variants: int = 4,
    return_null_scores: bool = False,
) -> EnrichmentResult:
    """Test cell-type link-score enrichment using circular genomic shifts."""

    if not isinstance(annotated_gwas, pd.DataFrame) or annotated_gwas.empty:
        raise EnrichmentError("annotated_gwas must be a non-empty DataFrame.")
    if not cell_types:
        raise EnrichmentError("At least one cell type is required.")
    if len(set(cell_types)) != len(cell_types):
        raise EnrichmentError("cell_types must be unique.")
    if isinstance(n_permutations, bool) or not isinstance(n_permutations, numbers.Integral):
        raise EnrichmentError("n_permutations must be an integer.")
    if n_permutations < 1:
        raise EnrichmentError("n_permutations must be positive.")
    if not math.isfinite(min_shift_fraction) or not 0 <= min_shift_fraction < 0.5:
        raise EnrichmentError("min_shift_fraction must lie in [0, 0.5).")
    if p_upper is not None and (not math.isfinite(p_upper) or not 0 < p_upper <= 1):
        raise EnrichmentError("p_upper must lie in (0, 1].")
    if not math.isfinite(min_p) or min_p <= 0 or min_p > 1:
        raise EnrichmentError("min_p must lie in (0, 1].")
    if not math.isfinite(max_evidence) or max_evidence <= 0:
        raise EnrichmentError("max_evidence must be positive.")
    if not math.isfinite(weight_floor) or weight_floor < 0:
        raise EnrichmentError("weight_floor must be non-negative.")
    if isinstance(min_chromosome_variants, bool) or not isinstance(
        min_chromosome_variants, numbers.Integral
    ):
        raise EnrichmentError("min_chromosome_variants must be an integer.")
    if min_chromosome_variants < 4:
        raise EnrichmentError("min_chromosome_variants must be at least 4.")
    if isinstance(random_state, bool) or not isinstance(random_state, numbers.Integral):
        raise EnrichmentError("random_state must be a non-negative integer.")
    if random_state < 0:
        raise EnrichmentError("random_state must be a non-negative integer.")
    if score_mode not in {"pvalue", "uniform", "zscore", "abs"}:
        raise EnrichmentError("score_mode must be 'pvalue', 'uniform', 'zscore', or 'abs'.")

    missing_link_columns = set(cell_types).difference(annotated_gwas.columns)
    if missing_link_columns:
        raise EnrichmentError(
            f"Missing cell-type annotation columns: {sorted(missing_link_columns)}"
        )
    reserved_columns = {"__chrom", "__pos0", "__row", "__eligible_row"}
    reserved_collisions = set(cell_types).intersection(reserved_columns)
    if reserved_collisions:
        raise EnrichmentError(
            f"Cell-type names use reserved internal columns: {sorted(reserved_collisions)}"
        )
    chrom_column = resolve_column(
        annotated_gwas.columns,
        ("chr", "chrom", "chromosome", "chr_id", "contig"),
        label="chromosome",
    )
    position_column = resolve_column(
        annotated_gwas.columns,
        (
            "pos",
            "position",
            "bp",
            "base_pair_location",
            "basepairlocation",
            "base_pair",
        ),
        label="position",
    )

    frame = annotated_gwas.copy()
    frame["__chrom"] = frame[chrom_column].map(normalize_chromosome)
    frame["__pos0"] = pd.to_numeric(frame[position_column], errors="coerce") - 1.0
    valid = (
        (frame["__chrom"] != "")
        & frame["__pos0"].notna()
        & np.isfinite(frame["__pos0"])
        & (frame["__pos0"] >= 0)
        & np.equal(frame["__pos0"], np.floor(frame["__pos0"]))
    )

    score_source: str | None = None
    score_center: float | None = None
    score_scale: float | None = None
    if score_mode == "pvalue":
        if p_column is None:
            p_column = resolve_column(
                frame.columns,
                (
                    "p",
                    "p_value",
                    "pvalue",
                    "p_val",
                    "pval",
                    "p-value",
                    "p_bolt_lmm_inf",
                    "p_bolt_lmm",
                ),
                label="P-value",
            )
        if p_column not in frame.columns:
            raise EnrichmentError(f"P-value column {p_column!r} is absent.")
        p_values = pd.to_numeric(frame[p_column], errors="coerce")
        valid &= p_values.notna() & (p_values >= 0) & (p_values <= 1)
        evidence = np.clip(-np.log10(p_values.clip(lower=min_p)), 0.0, max_evidence)
        evidence = np.maximum(evidence.to_numpy(dtype=float) - weight_floor, 0.0)
        if p_upper is not None:
            evidence[p_values.to_numpy(dtype=float) > p_upper] = 0.0
        score_source = p_column
    elif score_mode == "uniform":
        evidence = np.ones(len(frame), dtype=float)
    else:
        if score_column is None:
            score_column = resolve_column(
                frame.columns,
                ("z", "zscore", "z_score", "beta", "effect", "stat"),
                label="effect/statistic",
            )
        if score_column not in frame.columns:
            raise EnrichmentError(f"Score column {score_column!r} is absent.")
        raw_scores = pd.to_numeric(frame[score_column], errors="coerce")
        valid &= raw_scores.notna()
        if score_mode == "zscore":
            evidence = np.zeros(len(frame), dtype=float)
            standardized, score_center, score_scale = _robust_standardized_magnitudes(
                raw_scores.loc[valid].to_numpy(dtype=float),
                floor=weight_floor,
            )
            evidence[valid.to_numpy()] = standardized
        else:
            evidence = np.abs(raw_scores.to_numpy(dtype=float))
            if weight_floor:
                evidence = np.maximum(evidence - weight_floor, 0.0)
        score_source = score_column

    frame = frame.loc[valid].copy()
    evidence = np.asarray(evidence[valid.to_numpy()], dtype=float)
    if frame.empty:
        raise EnrichmentError("No valid GWAS rows remain after filtering.")
    if not np.isfinite(evidence).all():
        raise EnrichmentError("Evidence weights contain non-finite values.")
    if np.sum(evidence) <= 0:
        raise EnrichmentError(
            "All evidence weights are zero. Lower p_upper/weight_floor or use score_mode='uniform'."
        )

    frame = frame.reset_index(drop=True)
    frame["__row"] = np.arange(len(frame), dtype=np.int64)
    link_matrix = frame.loc[:, list(cell_types)].apply(pd.to_numeric, errors="coerce")
    if link_matrix.isna().any().any():
        raise EnrichmentError("Cell-type link scores contain non-numeric values.")
    link_values = link_matrix.to_numpy(dtype=float)
    if not np.isfinite(link_values).all() or np.any(link_values < 0):
        raise EnrichmentError("Cell-type link scores must be finite and non-negative.")

    eligible_parts: list[pd.DataFrame] = []
    skipped_chromosomes: dict[str, int] = {}
    for chrom, indices in frame.groupby("__chrom", sort=False).groups.items():
        size = len(indices)
        if size < min_chromosome_variants:
            skipped_chromosomes[str(chrom)] = int(size)
        else:
            eligible_parts.append(frame.loc[indices])
    if not eligible_parts:
        raise EnrichmentError(
            "No chromosome has enough variants for circular permutation. "
            f"Minimum required: {min_chromosome_variants}."
        )
    eligible = pd.concat(eligible_parts, axis=0).sort_values(
        ["__chrom", "__pos0"], kind="mergesort"
    )
    eligible_rows = eligible["__row"].to_numpy(dtype=np.int64)
    evidence = evidence[eligible_rows]
    link_values = link_values[eligible_rows]
    eligible = eligible.reset_index(drop=True)
    eligible["__eligible_row"] = np.arange(len(eligible), dtype=np.int64)
    if np.sum(evidence) <= 0:
        raise EnrichmentError("Eligible chromosomes have zero total evidence.")

    chromosomes: list[str] = []
    correlations: list[np.ndarray] = []
    weight_totals: list[float] = []
    observed_numerators = np.zeros(len(cell_types), dtype=float)
    allowed_offsets: list[np.ndarray] = []
    for chrom, group in eligible.groupby("__chrom", sort=True):
        source_index = group["__eligible_row"].to_numpy(dtype=np.int64)
        chrom_evidence = evidence[source_index]
        chrom_links = link_values[source_index]
        if np.sum(chrom_evidence) <= 0:
            continue
        correlation = circular_correlation_matrix(chrom_evidence, chrom_links)
        size = len(group)
        minimum_shift = max(1, int(math.ceil(size * min_shift_fraction)))
        offsets = np.arange(
            minimum_shift,
            size - minimum_shift + 1,
            dtype=np.int64,
        )
        if len(offsets) == 0:
            offsets = np.arange(1, size, dtype=np.int64)
        if len(offsets) == 0:
            continue
        chromosomes.append(str(chrom))
        correlations.append(correlation)
        allowed_offsets.append(offsets)
        weight_totals.append(float(np.sum(chrom_evidence)))
        observed_numerators += correlation[0]

    if not correlations:
        raise EnrichmentError("No chromosome produced a valid null distribution.")

    total_weight = float(sum(weight_totals))
    observed_scores = observed_numerators / total_weight
    rng = np.random.default_rng(random_state)
    sampled_offsets = np.column_stack(
        [rng.choice(offsets, size=n_permutations, replace=True) for offsets in allowed_offsets]
    )
    null_numerators = np.zeros((n_permutations, len(cell_types)), dtype=float)
    for column, correlation in enumerate(correlations):
        null_numerators += correlation[sampled_offsets[:, column], :]
    null_scores = null_numerators / total_weight
    null_mean = null_scores.mean(axis=0)
    null_sd = (
        null_scores.std(axis=0, ddof=1)
        if n_permutations > 1
        else np.zeros(len(cell_types), dtype=float)
    )
    tolerance = 1e-12 * np.maximum(1.0, np.abs(observed_scores))
    exceedances = np.sum(null_scores >= observed_scores[None, :] - tolerance, axis=0)
    empirical_p = (1.0 + exceedances) / (n_permutations + 1.0)
    z_scores = np.divide(
        observed_scores - null_mean,
        null_sd,
        out=np.zeros_like(observed_scores),
        where=null_sd > np.finfo(float).eps,
    )
    q_values = benjamini_hochberg(empirical_p)
    null_fold = np.divide(
        observed_scores,
        null_mean,
        out=np.full_like(observed_scores, np.nan),
        where=null_mean > np.finfo(float).eps,
    )

    linked_fraction = np.mean(link_values > 0, axis=0)
    evidence_on_linked = np.divide(
        evidence @ (link_values > 0),
        total_weight,
        out=np.zeros(len(cell_types), dtype=float),
        where=total_weight > 0,
    )
    table = pd.DataFrame(
        {
            "cell_type": list(cell_types),
            "n_variants": len(eligible),
            "n_chromosomes": len(chromosomes),
            "linked_variant_fraction": linked_fraction,
            "weighted_linked_fraction": evidence_on_linked,
            "observed_link_score": observed_scores,
            "null_mean": null_mean,
            "null_sd": null_sd,
            "z_score": z_scores,
            "fold_over_null": null_fold,
            "empirical_p": empirical_p,
            "fdr_q": q_values,
            "n_permutations": n_permutations,
            "random_state": random_state,
        }
    ).sort_values(["empirical_p", "z_score"], ascending=[True, False])
    table = table.reset_index(drop=True)

    if n_permutations < 1000:
        warnings.warn(
            "Fewer than 1000 permutations gives coarse empirical p-values; "
            "use at least 10,000 for publication-grade exploratory inference.",
            RuntimeWarning,
            stacklevel=2,
        )
    diagnostics: dict[str, object] = {
        "method": "within-chromosome circular permutation",
        "score_mode": score_mode,
        "score_source": score_source,
        "score_center": score_center,
        "score_scale": score_scale,
        "n_input_rows": int(len(annotated_gwas)),
        "n_eligible_rows": int(len(eligible)),
        "n_permutations": int(n_permutations),
        "random_state": int(random_state),
        "min_shift_fraction": float(min_shift_fraction),
        "chromosomes": chromosomes,
        "skipped_chromosomes": skipped_chromosomes,
        "total_evidence_weight": total_weight,
        "p_upper": p_upper,
        "n_positive_evidence_rows": int(np.count_nonzero(evidence)),
        "background": "All valid matched variants; p_upper zeros evidence without removing background rows",
        "multiple_testing": "Benjamini-Hochberg across cell types",
        "interpretation_limit": (
            "This is an annotation-shift enrichment test. It does not model LD "
            "directly, partition SNP heritability, or establish causality."
        ),
    }
    table.attrs["diagnostics"] = diagnostics
    return EnrichmentResult(
        table=table,
        diagnostics=diagnostics,
        null_scores=null_scores if return_null_scores else None,
    )


def run_enrichment_by_chromosome(
    annotated_chromosomes: Iterable[pd.DataFrame],
    *,
    cell_types: Sequence[str],
    p_column: str = "P",
    p_upper: float | None = None,
    weight_floor: float = 0.0,
    n_permutations: int = 10_000,
    random_state: int = 13,
) -> EnrichmentResult:
    """Combine independent chromosome shifts while keeping memory bounded.

    Each frame must contain exactly one distinct chromosome. The same statistic
    as ``run_enrichment`` is evaluated on the union of the frames. Independent
    child random seeds are used; samples therefore need not be bit-identical to
    a monolithic call with the same root seed.
    """

    seeds = np.random.SeedSequence(random_state)
    cell_types = list(cell_types)
    totals = np.zeros(len(cell_types))
    null_totals = np.zeros((n_permutations, len(cell_types)))
    coverage = np.zeros(len(cell_types))
    weighted_coverage = np.zeros(len(cell_types))
    total_weight = 0.0
    n_rows = 0
    seen: set[str] = set()
    active_chromosomes: list[str] = []
    details: list[dict[str, object]] = []
    for frame in annotated_chromosomes:
        chrom_column = resolve_column(
            frame.columns, ("chr", "chrom", "chromosome"), label="chromosome"
        )
        chromosomes = frame[chrom_column].map(normalize_chromosome).unique()
        if len(chromosomes) != 1 or chromosomes[0] in seen or not chromosomes[0]:
            raise EnrichmentError(
                "Each stream frame must contain one unique, previously unseen chromosome."
            )
        chrom = str(chromosomes[0])
        seen.add(chrom)
        p = pd.to_numeric(frame[p_column], errors="raise").to_numpy(dtype=float)
        if np.any(~np.isfinite(p)) or np.any((p < 0) | (p > 1)):
            raise EnrichmentError("Streamed P-values must be finite and lie in [0, 1].")
        evidence = np.maximum(
            np.minimum(-np.log10(np.maximum(p, 1e-300)), 300.0) - weight_floor, 0.0
        )
        if p_upper is not None:
            evidence[p > p_upper] = 0.0
        links = frame[cell_types].to_numpy(dtype=float)
        if np.any(~np.isfinite(links)) or np.any(links < 0):
            raise EnrichmentError("Streamed link scores must be finite and non-negative.")
        coverage += (links > 0).sum(axis=0)
        n_rows += len(frame)
        if not evidence.any():
            details.append(
                {"chromosome": chrom, "n_rows": len(frame), "total_evidence_weight": 0.0}
            )
            continue
        child_seed = int(seeds.spawn(1)[0].generate_state(1)[0])
        result = run_enrichment(
            frame,
            cell_types=cell_types,
            p_column=p_column,
            p_upper=p_upper,
            weight_floor=weight_floor,
            n_permutations=n_permutations,
            random_state=child_seed,
            return_null_scores=True,
        )
        table = result.table.set_index("cell_type").loc[cell_types]
        weight = float(result.diagnostics["total_evidence_weight"])
        totals += table["observed_link_score"].to_numpy() * weight
        null_totals += result.null_scores * weight
        weighted_coverage += table["weighted_linked_fraction"].to_numpy() * weight
        total_weight += weight
        active_chromosomes.append(chrom)
        details.append(result.diagnostics)
    if total_weight <= 0:
        raise EnrichmentError("All streamed evidence weights are zero.")
    observed = totals / total_weight
    null = null_totals / total_weight
    mean = null.mean(axis=0)
    sd = null.std(axis=0, ddof=1) if n_permutations > 1 else np.zeros(len(cell_types))
    tolerance = 1e-12 * np.maximum(1.0, np.abs(observed))
    p = (1 + np.sum(null >= observed[None, :] - tolerance, axis=0)) / (n_permutations + 1)
    table = (
        pd.DataFrame(
            {
                "cell_type": cell_types,
                "n_variants": n_rows,
                "n_chromosomes": len(active_chromosomes),
                "linked_variant_fraction": coverage / n_rows,
                "weighted_linked_fraction": weighted_coverage / total_weight,
                "observed_link_score": observed,
                "null_mean": mean,
                "null_sd": sd,
                "z_score": np.divide(
                    observed - mean, sd, out=np.zeros_like(mean), where=sd > np.finfo(float).eps
                ),
                "fold_over_null": np.divide(
                    observed, mean, out=np.full_like(mean, np.nan), where=mean > np.finfo(float).eps
                ),
                "empirical_p": p,
                "fdr_q": benjamini_hochberg(p),
                "n_permutations": n_permutations,
                "random_state": random_state,
            }
        )
        .sort_values(["empirical_p", "z_score"], ascending=[True, False])
        .reset_index(drop=True)
    )
    diagnostics = {
        "method": "within-chromosome circular permutation, chromosome streaming",
        "n_eligible_rows": n_rows,
        "chromosomes": sorted(seen),
        "chromosomes_with_evidence": active_chromosomes,
        "total_evidence_weight": total_weight,
        "n_permutations": n_permutations,
        "random_state": random_state,
        "p_upper": p_upper,
        "weight_floor": weight_floor,
        "chromosome_diagnostics": details,
        "background": "All valid matched variants; p_upper zeros evidence without removing rows",
        "multiple_testing": "Benjamini-Hochberg across cell types",
        "interpretation_limit": "Annotation-shift enrichment; no direct LD model or causal inference.",
    }
    table.attrs["diagnostics"] = diagnostics
    return EnrichmentResult(table=table, diagnostics=diagnostics)


def circular_permutation_enrichment(
    annotated_gwas: pd.DataFrame,
    **kwargs: object,
) -> pd.DataFrame:
    """Convenience wrapper returning the enrichment dataframe.

    Diagnostics are attached to ``DataFrame.attrs['diagnostics']``.
    """

    return run_enrichment(annotated_gwas, **kwargs).table
