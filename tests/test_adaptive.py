from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from celltogwas.adaptive import (
    _holm,
    _residualize,
    calibrate_adaptive,
    run_adaptive_gene_enrichment,
)
from celltogwas.cli import main


def fixture_genes(n=80):
    rng = np.random.default_rng(5)
    return pd.DataFrame(
        {
            "gene": [f"g{i}" for i in range(n)],
            "chrom": np.repeat([1, 2], n // 2),
            "pos": np.tile(np.arange(n // 2), 2),
            "z": rng.normal(size=n),
            "A": rng.random(n),
            "B": rng.random(n),
            "nuisance": rng.normal(size=n),
        }
    )


def test_calibration_matches_independent_brute_force_with_ties():
    rng = np.random.default_rng(9)
    stats = rng.integers(-4, 5, (31, 4, 2)).astype(float)
    adaptive, components = calibrate_adaptive(stats)
    brute = np.empty_like(stats)
    for row in range(31):
        brute[row] = (stats >= stats[row]).mean(axis=0)
    minimum = brute.min(axis=1)
    expected = np.stack([(minimum <= minimum[i]).mean(axis=0) for i in range(31)])
    np.testing.assert_allclose(components, brute)
    np.testing.assert_allclose(adaptive, expected)
    assert np.all(adaptive >= minimum)
    assert np.all(adaptive <= np.minimum(1, 4 * minimum) + 1e-12)


def test_row_exchangeability_and_duplicate_kernel_invariance():
    rng = np.random.default_rng(19)
    stats = rng.normal(size=(71, 3, 4))
    p, _ = calibrate_adaptive(stats)
    order = rng.permutation(71)
    shuffled, _ = calibrate_adaptive(stats[order])
    duplicated, _ = calibrate_adaptive(np.concatenate((stats, stats), axis=1))
    np.testing.assert_array_equal(shuffled, p[order])
    np.testing.assert_array_equal(duplicated, p)


def test_fft_signed_residuals_match_direct_rolls_and_shared_shifts():
    genes = fixture_genes()
    result = run_adaptive_gene_enrichment(
        genes,
        cell_types=["A", "B"],
        powers=[1],
        thresholds=[None],
        binary_fraction=None,
        n_permutations=39,
        return_null_statistics=True,
    )
    expected = np.zeros((40, 1, 2))
    seeds = np.random.SeedSequence(13).spawn(2)
    weights = genes[["A", "B"]].to_numpy(copy=True)
    weights /= weights.max(axis=0)
    for chrom, seed in zip([1, 2], seeds, strict=False):
        rows = np.flatnonzero(genes.chrom == chrom)
        z = genes.z.to_numpy()[rows].copy()
        z -= z.mean()
        w = weights[rows] - weights[rows].mean(axis=0)
        shifts = np.r_[0, np.random.default_rng(seed).integers(0, len(rows), 39)]
        for i, shift in enumerate(shifts):
            expected[i, 0] += z @ np.roll(w, -shift, axis=0)
    np.testing.assert_allclose(result.null_statistics, expected, atol=1e-12)
    assert expected.min() < 0


def test_reproducible_row_order_and_collinear_covariates():
    genes = fixture_genes()
    genes["twice"] = 2 * genes.nuisance
    kwargs = dict(cell_types=["A", "B"], covariates=["nuisance", "twice"], n_permutations=99)
    first = run_adaptive_gene_enrichment(genes, **kwargs)
    second = run_adaptive_gene_enrichment(genes.sample(frac=1, random_state=3), **kwargs)
    pd.testing.assert_frame_equal(first.table, second.table)
    assert all(c["nuisance_rank"] == 2 for c in first.diagnostics["chromosomes"])
    assert np.all(first.table.p_holm >= first.table.p_adaptive)
    assert np.all(first.table.q_by >= first.table.q_bh_exploratory)


def test_constant_zero_weights_and_roundoff_are_not_significant():
    genes = fixture_genes()
    genes["A"] = 1.0
    genes["B"] = 0.0
    result = run_adaptive_gene_enrichment(genes, cell_types=["A", "B"], n_permutations=39)
    assert (result.table.p_adaptive == 1).all()
    assert (result.components.p_component_uncorrected == 1).all()
    p, _ = calibrate_adaptive(np.ones((30, 2, 1)) + np.arange(30)[:, None, None] * 1e-14)
    assert (p == 1).all()


def test_projection_removes_covariates_and_holm_known_answer():
    rng = np.random.default_rng(11)
    x = rng.normal(size=(50, 2))
    residual, rank = _residualize(rng.normal(size=(50, 3)) + x @ np.ones((2, 3)), x)
    assert rank == 3
    np.testing.assert_allclose(x.T @ residual, 0, atol=1e-12)
    np.testing.assert_allclose(residual.mean(axis=0), 0, atol=1e-12)
    np.testing.assert_allclose(_holm(np.array([0.01, 0.04, 0.03, 0.9])), [0.04, 0.09, 0.09, 0.9])


@pytest.mark.parametrize(
    "change",
    [
        {"powers": []},
        {"powers": [0]},
        {"powers": [1, 1]},
        {"powers": [np.inf]},
        {"thresholds": []},
        {"thresholds": [-1]},
        {"thresholds": [10]},
        {"thresholds": [None, None]},
        {"binary_fraction": 1},
        {"z_cap": 0},
        {"n_permutations": True},
        {"n_permutations": 1},
        {"random_state": -1},
        {"cell_types": ["A", "A"]},
        {"cell_types": ["z"]},
        {"covariates": ["A"]},
    ],
)
def test_invalid_options(change):
    kwargs = {"cell_types": ["A"], "n_permutations": 39, **change}
    with pytest.raises(ValueError):
        run_adaptive_gene_enrichment(fixture_genes(), **kwargs)


@pytest.mark.parametrize(
    "column,value",
    [
        ("z", np.nan),
        ("A", -1),
        ("pos", -1),
        ("pos", 1.5),
        ("gene", None),
        ("chrom", None),
    ],
)
def test_invalid_rows(column, value):
    genes = fixture_genes()
    if column == "pos":
        genes[column] = genes[column].astype(float)
    genes.loc[0, column] = value
    with pytest.raises(ValueError):
        run_adaptive_gene_enrichment(genes, cell_types=["A"], n_permutations=39)


def test_duplicate_ids_missing_columns_small_chromosomes_and_saturated_design():
    genes = fixture_genes()
    genes.loc[0, "gene"] = genes.loc[1, "gene"]
    with pytest.raises(ValueError, match="unique"):
        run_adaptive_gene_enrichment(genes, cell_types=["A"])
    with pytest.raises(ValueError, match="Missing"):
        run_adaptive_gene_enrichment(fixture_genes().drop(columns="z"), cell_types=["A"])
    with pytest.raises(ValueError, match="fewer than 8"):
        run_adaptive_gene_enrichment(fixture_genes(12), cell_types=["A"])
    genes = fixture_genes(16)
    for i in range(8):
        genes[f"x{i}"] = np.tile(np.eye(8)[:, i], 2)
    with pytest.raises(ValueError, match="degrees of freedom"):
        run_adaptive_gene_enrichment(
            genes, cell_types=["A"], covariates=[f"x{i}" for i in range(8)]
        )


@pytest.mark.parametrize("bad", [np.ones((5, 2)), np.ones((1, 2, 3)), np.full((3, 2, 1), np.nan)])
def test_invalid_calibration(bad):
    with pytest.raises(ValueError):
        calibrate_adaptive(bad)


def test_cli_runs_and_records_input_hash(tmp_path):
    path = tmp_path / "genes.tsv"
    fixture_genes().to_csv(path, sep="\t", index=False)
    out, components, report = [
        tmp_path / name for name in ["out.tsv", "components.tsv", "report.json"]
    ]
    assert (
        main(
            [
                "adaptive-gene",
                "--genes",
                str(path),
                "--cell-types",
                "A",
                "B",
                "--covariates",
                "nuisance",
                "--permutations",
                "39",
                "--output",
                str(out),
                "--components-output",
                str(components),
                "--report",
                str(report),
            ]
        )
        == 0
    )
    assert len(pd.read_csv(out, sep="\t")) == 2
    assert len(pd.read_csv(components, sep="\t")) == 24
    assert '"input_sha256"' in report.read_text()
