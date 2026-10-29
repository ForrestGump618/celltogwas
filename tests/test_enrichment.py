from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from celltogwas.enrichment import (
    _robust_standardized_magnitudes,
    benjamini_hochberg,
    circular_correlation_matrix,
    run_enrichment,
)


def test_fft_circular_correlation_matches_observed_and_roll():
    rng = np.random.default_rng(7)
    weights = rng.random(20)
    scores = rng.random((20, 2))

    correlation = circular_correlation_matrix(weights, scores)

    np.testing.assert_allclose(correlation[0], weights @ scores, rtol=1e-12, atol=1e-12)
    for shift in (1, 4, 11):
        expected = np.array(
            [
                np.dot(weights, np.roll(scores[:, column], -shift))
                for column in range(scores.shape[1])
            ]
        )
        np.testing.assert_allclose(correlation[shift], expected, rtol=1e-12, atol=1e-12)


def test_benjamini_hochberg_adjustment():
    adjusted = benjamini_hochberg([0.01, 0.04, 0.03, 0.2])
    np.testing.assert_allclose(adjusted, [0.04, 0.0533333333, 0.0533333333, 0.2])


def test_robust_effect_standardization_handles_zero_mad():
    evidence, center, scale = _robust_standardized_magnitudes(np.array([0.0, 0.0, 0.0, 10.0]))

    assert center == 0.0
    assert scale == 2.5
    np.testing.assert_allclose(evidence, [0.0, 0.0, 0.0, 4.0])


def test_zscore_score_mode_records_center_and_scale():
    gwas = pd.DataFrame(
        {
            "CHR": "1",
            "POS": np.arange(1, 21),
            "EFFECT": [0.1] * 19 + [5.0],
            "CellA": 1.0,
        }
    )

    result = run_enrichment(
        gwas,
        cell_types=["CellA"],
        score_mode="zscore",
        score_column="EFFECT",
        n_permutations=20,
        random_state=3,
    )

    assert result.diagnostics["score_center"] == pytest.approx(0.1)
    assert result.diagnostics["score_scale"] > 0


def test_spatial_enrichment_detects_colocalized_signal():
    positions = np.arange(1, 201)
    cell_a = ((positions >= 50) & (positions <= 65)).astype(float)
    cell_b = ((positions >= 150) & (positions <= 190)).astype(float)
    p_values = np.full(len(positions), 0.5)
    p_values[cell_a.astype(bool)] = 1e-8
    gwas = pd.DataFrame(
        {
            "CHR": "1",
            "POS": positions,
            "P": p_values,
            "CellA": cell_a,
            "CellB": cell_b,
        }
    )

    result = run_enrichment(
        gwas,
        cell_types=["CellA", "CellB"],
        n_permutations=399,
        random_state=4,
        min_shift_fraction=0.01,
    )

    table = result.table.set_index("cell_type")
    assert table.loc["CellA", "empirical_p"] < table.loc["CellB", "empirical_p"]
    assert table.loc["CellA", "fdr_q"] >= table.loc["CellA", "empirical_p"]
    assert result.diagnostics["chromosomes"] == ["1"]
    assert "heritability" in result.diagnostics["interpretation_limit"].lower()


def test_enrichment_rejects_all_zero_evidence():
    gwas = pd.DataFrame(
        {
            "CHR": "1",
            "POS": np.arange(1, 21),
            "P": 1.0,
            "CellA": 1.0,
        }
    )
    with pytest.raises(ValueError, match="All evidence weights are zero"):
        run_enrichment(
            gwas,
            cell_types=["CellA"],
            weight_floor=1.0,
            n_permutations=10,
        )


def test_underflowed_p_zero_is_capped_instead_of_dropped():
    gwas = pd.DataFrame(
        {
            "CHR": "1",
            "POS": np.arange(1, 21),
            "P": [0.0, 0.0, 0.0] + [0.5] * 17,
            "CellA": 1.0,
        }
    )
    result = run_enrichment(
        gwas,
        cell_types=["CellA"],
        n_permutations=20,
        random_state=1,
    )
    assert result.diagnostics["n_input_rows"] == 20
    assert result.diagnostics["n_eligible_rows"] == 20


def test_p_upper_retains_background_and_zeros_only_evidence():
    frame = pd.DataFrame(
        {
            "CHR": "1",
            "POS": np.arange(1, 101),
            "P": [1e-10] + [0.5] * 99,
            "CellA": [1.0] + [0.0] * 99,
        }
    )
    result = run_enrichment(
        frame, cell_types=["CellA"], p_upper=5e-8, n_permutations=1000, return_null_scores=True
    )
    assert result.diagnostics["n_eligible_rows"] == 100
    assert result.diagnostics["n_positive_evidence_rows"] == 1
    assert result.table.loc[0, "observed_link_score"] == pytest.approx(1)
    assert result.table.loc[0, "empirical_p"] < 0.01
    assert result.null_scores.shape == (1000, 1)


def test_streaming_observed_statistic_matches_pooled_and_rejects_repeated_chromosomes():
    from celltogwas.enrichment import run_enrichment_by_chromosome

    first = pd.DataFrame(
        {
            "CHR": "1",
            "POS": np.arange(1, 101),
            "P": [1e-10] + [0.5] * 99,
            "CellA": [1.0] + [0.0] * 99,
        }
    )
    second = first.assign(CHR="2", P=[1e-9] + [0.5] * 99)
    options = dict(cell_types=["CellA"], p_upper=5e-8, n_permutations=1000)
    streamed = run_enrichment_by_chromosome([first, second], **options)
    pooled = run_enrichment(pd.concat([first, second]), **options)
    assert streamed.table.loc[0, "observed_link_score"] == pytest.approx(
        pooled.table.loc[0, "observed_link_score"]
    )
    assert streamed.table.loc[0, "empirical_p"] < 0.01
    with pytest.raises(ValueError, match="previously unseen"):
        run_enrichment_by_chromosome([first, first], **options)


def test_fft_ties_in_constant_annotation_have_p_one():
    frame = pd.DataFrame(
        {"CHR": "1", "POS": np.arange(1, 102), "P": np.linspace(1e-12, 0.99, 101), "CellA": 0.3}
    )
    result = run_enrichment(frame, cell_types=["CellA"], n_permutations=1000)
    assert result.table.loc[0, "empirical_p"] == 1


def test_streamed_null_mean_matches_independent_brute_force_shifts():
    from celltogwas.enrichment import run_enrichment_by_chromosome

    weights = [np.array([5.0, 0, 2, 0, 0, 0, 0, 1]), np.array([0.0, 3, 0, 0, 1, 0, 0, 0])]
    scores = [np.array([1.0, 1, 0, 0, 0, 0, 0, 0]), np.array([0.0, 1, 1, 1, 0, 0, 0, 0])]
    frames = [
        pd.DataFrame({"CHR": str(i + 1), "POS": np.arange(1, 9), "P": 10**-w, "CellA": v})
        for i, (w, v) in enumerate(zip(weights, scores, strict=True))
    ]
    expected = sum(
        np.mean([np.dot(w, np.roll(v, -k)) for k in range(1, 8)])
        for w, v in zip(weights, scores, strict=True)
    ) / sum(w.sum() for w in weights)
    result = run_enrichment_by_chromosome(frames, cell_types=["CellA"], n_permutations=10000)
    assert result.table.loc[0, "null_mean"] == pytest.approx(expected, abs=0.01)
