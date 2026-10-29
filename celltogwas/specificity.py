"""Cell-type expression-specificity scores.

The implementation combines a gene-level specificity index with the gene's
relative expression in each cell type. The relative-expression factor weights the cell types within each gene.
Use ``min_expression`` to exclude genes with negligible absolute signal.
"""

from __future__ import annotations

import math
from typing import Literal

import numpy as np
import pandas as pd

from .exceptions import InputFormatError

SpecificityMethod = Literal["tau", "entropy", "zscore"]


def _prepare_expression(expression: pd.DataFrame, pseudocount: float) -> pd.DataFrame:
    if not isinstance(expression, pd.DataFrame):
        raise InputFormatError("expression must be a pandas DataFrame.")
    if expression.empty:
        raise InputFormatError("Expression matrix is empty.")
    if expression.index.has_duplicates:
        raise InputFormatError("Expression matrix contains duplicate gene identifiers.")
    if expression.columns.has_duplicates:
        raise InputFormatError("Expression matrix contains duplicate cell-type columns.")
    try:
        numeric = expression.apply(pd.to_numeric, errors="raise").astype(float)
    except (TypeError, ValueError) as error:
        raise InputFormatError("Expression matrix must be numeric.") from error
    if not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise InputFormatError("Expression matrix contains NaN or infinite values.")
    if float(numeric.to_numpy().min()) < 0:
        raise InputFormatError(
            "Expression values must be non-negative. Use normalized counts, "
            "CPM/TPM, or a non-negative imputed matrix."
        )
    if not math.isfinite(pseudocount) or pseudocount <= 0:
        raise InputFormatError("pseudocount must be positive.")
    return numeric


def _tau_specificity(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    maxima = values.max(axis=1)
    safe = np.divide(
        values,
        maxima[:, None],
        out=np.zeros_like(values),
        where=maxima[:, None] > 0,
    )
    if values.shape[1] == 1:
        tau = np.ones(values.shape[0], dtype=float)
    else:
        tau = np.sum(1.0 - safe, axis=1) / (values.shape[1] - 1)
    return tau, safe


def _entropy_specificity(values: np.ndarray, pseudocount: float) -> tuple[np.ndarray, np.ndarray]:
    maxima = values.max(axis=1)
    relative = np.divide(
        values,
        maxima[:, None],
        out=np.zeros_like(values),
        where=maxima[:, None] > 0,
    )
    adjusted = values + pseudocount
    probabilities = adjusted / adjusted.sum(axis=1, keepdims=True)
    if values.shape[1] == 1:
        return np.ones(values.shape[0], dtype=float), relative
    entropy = -np.sum(probabilities * np.log(probabilities), axis=1)
    specificity = 1.0 - entropy / math.log(values.shape[1])
    return np.clip(specificity, 0.0, 1.0), relative


def _zscore_weights(
    values: np.ndarray,
    pseudocount: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    maxima = values.max(axis=1)
    relative = np.divide(
        values,
        maxima[:, None],
        out=np.zeros_like(values),
        where=maxima[:, None] > 0,
    )
    if values.shape[1] == 1:
        specificity = np.ones(values.shape[0], dtype=float)
        return specificity, relative, relative.copy()

    medians = np.median(values, axis=1, keepdims=True)
    deviations = np.abs(values - medians)
    mad = np.median(deviations, axis=1, keepdims=True)
    scale = 1.4826 * mad
    # A zero MAD is common when most cell types have zero expression. In that
    # case the mean absolute deviation about the median retains relative scale.
    zero_mad = scale <= pseudocount
    mean_absolute_deviation = deviations.mean(axis=1, keepdims=True)
    scale = np.where(zero_mad, mean_absolute_deviation, scale)
    z_scores = np.divide(
        values - medians,
        scale,
        out=np.zeros_like(values),
        where=scale > pseudocount,
    )
    # A bounded exponential transform keeps weights interpretable on [0, 1).
    cell_weights = 1.0 - np.exp(-np.maximum(z_scores, 0.0))
    specificity = cell_weights.mean(axis=1)
    return np.clip(specificity, 0.0, 1.0), relative, cell_weights


def compute_specificity(
    expression: pd.DataFrame,
    *,
    method: SpecificityMethod = "tau",
    min_expression: float = 0.0,
    pseudocount: float = 1e-9,
    power: float = 1.0,
    top_n: int | None = None,
    min_weight: float = 0.0,
) -> pd.DataFrame:
    """Compute long-form cell-type marker weights.

    Parameters
    ----------
    expression:
        Genes in rows and cell types in columns. Values must be non-negative,
        such as TPM, CPM, or normalized pseudobulk counts.
    method:
        ``tau`` uses the Tau specificity index, ``entropy`` uses normalized
        Shannon entropy, and ``zscore`` uses a median/MAD cross-cell-type
        z-score with a bounded transform.
    min_expression:
        Drop genes whose maximum expression is below this value.
    power:
        Exponent applied to the final marker weight.
    top_n:
        Retain at most this many genes per cell type.
    min_weight:
        Drop markers below this final weight.
    """

    if method not in {"tau", "entropy", "zscore"}:
        raise InputFormatError("method must be 'tau', 'entropy', or 'zscore'.")
    if not math.isfinite(power) or power <= 0:
        raise InputFormatError("power must be positive.")
    if not math.isfinite(min_expression) or min_expression < 0:
        raise InputFormatError("min_expression must be finite and non-negative.")
    if top_n is not None and top_n < 1:
        raise InputFormatError("top_n must be positive when provided.")
    if not math.isfinite(min_weight) or min_weight < 0:
        raise InputFormatError("min_weight must be non-negative.")

    values = _prepare_expression(expression, pseudocount).to_numpy(dtype=float)
    maxima = values.max(axis=1)
    keep = (maxima > 0) & (maxima >= min_expression)
    values = values[keep]
    maxima = maxima[keep]
    genes = expression.index.to_numpy()[keep]
    cell_types = np.asarray(expression.columns, dtype=object)
    if len(values) == 0:
        raise InputFormatError("No genes pass min_expression.")

    if method == "tau":
        specificity, relative = _tau_specificity(values)
        cell_weights = relative * specificity[:, None]
    elif method == "entropy":
        specificity, relative = _entropy_specificity(values, pseudocount)
        cell_weights = relative * specificity[:, None]
    else:
        specificity, relative, cell_weights = _zscore_weights(values, pseudocount)

    cell_weights = np.power(np.clip(cell_weights, 0.0, 1.0), power)
    output = pd.DataFrame(
        {
            "cell_type": np.tile(cell_types.astype(str), len(genes)),
            "gene": np.repeat(genes.astype(str), len(cell_types)),
            "weight": cell_weights.ravel(order="C"),
            "specificity": np.repeat(specificity, len(cell_types)),
            "relative_expression": relative.ravel(order="C"),
            "max_expression": np.repeat(maxima, len(cell_types)),
        }
    )
    output = output.loc[(output["weight"] >= min_weight) & (output["weight"] > 0)].reset_index(
        drop=True
    )
    if top_n is not None:
        output = (
            output.sort_values(["cell_type", "weight"], ascending=[True, False])
            .groupby("cell_type", sort=False, as_index=False)
            .head(top_n)
            .reset_index(drop=True)
        )
    if output.empty:
        raise InputFormatError("No positive marker weights remain after filtering.")
    return output.sort_values(
        ["cell_type", "weight", "gene"], ascending=[True, False, True]
    ).reset_index(drop=True)


def specificity_to_gene_sets(
    specificity: pd.DataFrame,
    *,
    min_weight: float = 0.0,
    top_n: int | None = None,
) -> dict[str, dict[str, float]]:
    """Convert a long specificity table to weighted gene sets."""

    if not math.isfinite(min_weight) or min_weight < 0:
        raise InputFormatError("min_weight must be finite and non-negative.")
    if top_n is not None and top_n < 1:
        raise InputFormatError("top_n must be positive when provided.")
    required = {"cell_type", "gene", "weight"}
    missing = required.difference(specificity.columns)
    if missing:
        raise InputFormatError(f"Specificity table is missing columns: {sorted(missing)}")
    numeric_weight = pd.to_numeric(specificity["weight"], errors="coerce")
    present = numeric_weight.notna()
    if np.any(~np.isfinite(numeric_weight.loc[present].to_numpy(dtype=float))):
        raise InputFormatError("Specificity weights must be finite.")
    frame = specificity.loc[(numeric_weight >= min_weight) & (numeric_weight > 0)].copy()
    frame["weight"] = numeric_weight.loc[frame.index].astype(float)
    if top_n is not None:
        frame = (
            frame.sort_values(["cell_type", "weight"], ascending=[True, False])
            .groupby("cell_type", sort=False, as_index=False)
            .head(top_n)
        )
    output: dict[str, dict[str, float]] = {}
    for cell_type, group in frame.groupby("cell_type", sort=False):
        weights: dict[str, float] = {}
        for gene, weight in zip(
            group["gene"].astype(str),
            pd.to_numeric(group["weight"], errors="coerce"),
            strict=False,
        ):
            if pd.notna(weight) and float(weight) > 0:
                weights[gene] = max(weights.get(gene, 0.0), float(weight))
        if weights:
            output[str(cell_type)] = weights
    if not output:
        raise InputFormatError("No gene sets remain after marker filtering.")
    return output


def specificity_matrix(specificity: pd.DataFrame) -> pd.DataFrame:
    """Pivot a long specificity table into genes by cell types."""

    required = {"cell_type", "gene", "weight"}
    missing = required.difference(specificity.columns)
    if missing:
        raise InputFormatError(f"Specificity table is missing columns: {sorted(missing)}")
    return (
        specificity.pivot(index="gene", columns="cell_type", values="weight")
        .fillna(0.0)
        .sort_index()
        .sort_index(axis=1)
    )


def specificity_summary(specificity: pd.DataFrame) -> pd.DataFrame:
    """Summarize marker count and weight distribution by cell type."""

    required = {"cell_type", "gene", "weight"}
    missing = required.difference(specificity.columns)
    if missing:
        raise InputFormatError(f"Specificity table is missing columns: {sorted(missing)}")
    return (
        specificity.groupby("cell_type", sort=False)
        .agg(
            n_markers=("gene", "nunique"),
            total_weight=("weight", "sum"),
            mean_weight=("weight", "mean"),
            max_weight=("weight", "max"),
        )
        .reset_index()
    )
