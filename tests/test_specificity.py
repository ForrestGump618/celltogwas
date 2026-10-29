from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from celltogwas.specificity import (
    compute_specificity,
    specificity_matrix,
    specificity_to_gene_sets,
)


def expression_fixture() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "T_cell": [10.0, 5.0, 0.0, 1.0],
            "B_cell": [0.0, 5.0, 6.0, 1.0],
            "Myeloid": [0.0, 5.0, 2.0, 8.0],
        },
        index=["CD3D", "HOUSEKEEPING", "MS4A1", "LYZ"],
    )


def test_tau_specificity_preserves_cell_gene_alignment():
    result = compute_specificity(expression_fixture(), method="tau")
    matrix = specificity_matrix(result)

    assert matrix.loc["CD3D", "T_cell"] == 1.0
    assert matrix.loc["CD3D", "B_cell"] == 0.0
    assert matrix.loc["MS4A1", "B_cell"] > matrix.loc["MS4A1", "T_cell"]
    assert "HOUSEKEEPING" not in matrix.index


def test_entropy_and_top_n_produce_weighted_gene_sets():
    result = compute_specificity(
        expression_fixture(),
        method="entropy",
        min_expression=0.1,
        top_n=2,
    )
    gene_sets = specificity_to_gene_sets(result, min_weight=0.01)

    assert "T_cell" in gene_sets
    assert "CD3D" in gene_sets["T_cell"]
    assert len(gene_sets["Myeloid"]) <= 2
    assert np.isfinite(result["weight"]).all()


def test_zscore_weights_are_bounded():
    result = compute_specificity(expression_fixture(), method="zscore")
    assert result["weight"].between(0.0, 1.0).all()


def test_expression_filtering_keeps_max_expression_aligned():
    expression = pd.DataFrame(
        {
            "T_cell": [10.0, 0.0, 3.0],
            "B_cell": [0.0, 0.0, 1.0],
        },
        index=["KEPT", "ZERO", "BELOW_THRESHOLD"],
    )

    result = compute_specificity(
        expression,
        method="tau",
        min_expression=5.0,
    )

    assert set(result["gene"]) == {"KEPT"}
    assert result.loc[result["gene"] == "KEPT", "max_expression"].tolist() == [10.0]


def test_zscore_uses_median_scale_when_most_cell_types_are_zero():
    expression = pd.DataFrame(
        {"T_cell": [10.0], "B_cell": [0.0], "Myeloid": [0.0]},
        index=["CD3D"],
    )

    result = compute_specificity(expression, method="zscore")

    assert result.loc[0, "cell_type"] == "T_cell"
    assert result.loc[0, "weight"] == pytest.approx(1.0 - np.exp(-3.0))
    assert result.loc[0, "weight"] > 0.9


def test_specificity_gene_set_threshold_is_inclusive_and_positive():
    table = pd.DataFrame(
        {
            "cell_type": ["T_cell", "T_cell", "B_cell"],
            "gene": ["CD3D", "CD3E", "MS4A1"],
            "weight": [0.5, 0.49, 0.0],
        }
    )

    gene_sets = specificity_to_gene_sets(table, min_weight=0.5)

    assert gene_sets == {"T_cell": {"CD3D": 0.5}}


def test_specificity_to_gene_sets_rejects_nonfinite_weights():
    table = pd.DataFrame(
        {
            "cell_type": ["T_cell"],
            "gene": ["CD3D"],
            "weight": [float("inf")],
        }
    )

    with pytest.raises(ValueError, match="finite"):
        specificity_to_gene_sets(table)
