from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from validation.benchmark_adaptive import ring_noise, wilson
from validation.run_adaptive_real import prepare_inputs


def inputs():
    n, cells = 120, ["A", "B", "C"]
    ids = [f"id{i}" for i in range(n)]
    symbols = [f"g{i}" for i in range(n)]
    annotation = pd.DataFrame(
        {
            "gene_id": [f"{i}.1" for i in ids],
            "gene_name": symbols,
            "chrom": "1",
            "start": np.arange(n) * 1000,
            "end": np.arange(n) * 1000 + 100,
        }
    )
    expression = pd.DataFrame(np.random.default_rng(7).random((n, 3)), index=symbols, columns=cells)
    detection, donors = expression * 0 + 0.5, expression * 0 + 10
    result = pd.DataFrame(
        {"GENE": ids, "ZSTAT": np.arange(n) / n, "P": 0.5, "NSNPS": 10, "NPARAM": 4, "N": 1000}
    )
    return annotation, expression, detection, donors, {"standard": result, "hmagma": result.copy()}


def test_common_background_no_missingness_imputation_and_ambiguous_aliases():
    annotation, expression, detection, donors, results = inputs()
    donors.loc["g0", "A"] = 0
    annotation.loc[2, "gene_name"] = "g1"
    results["hmagma"] = results["hmagma"].iloc[1:-1]
    frames, cells, manifest = prepare_inputs(annotation, expression, detection, donors, results)
    assert cells == ["A", "B", "C"]
    assert manifest["common_genes"] == 116
    assert not set(frames["standard"].gene) & {"id0", "id1", "id2", "id119"}
    pd.testing.assert_frame_equal(frames["standard"], frames["hmagma"])
    assert np.isfinite(frames["standard"][["z", *cells]].to_numpy()).all()


@pytest.mark.parametrize("which", [1, 2, 3])
def test_nonfinite_profiles_are_not_silently_set_to_zero(which):
    args = inputs()
    args[which].iloc[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        prepare_inputs(*args)


def test_misaligned_rows_duplicate_ids_and_too_small_background_fail():
    args = list(inputs())
    args[2] = args[2].iloc[::-1]
    with pytest.raises(ValueError, match="rows"):
        prepare_inputs(*args)
    args = inputs()
    args[4]["hmagma"].loc[0, "GENE"] = "id1"
    with pytest.raises(ValueError, match="unique"):
        prepare_inputs(*args)
    args = inputs()
    args[4]["hmagma"] = args[4]["hmagma"].iloc[:10]
    with pytest.raises(ValueError, match="Too few"):
        prepare_inputs(*args)


def test_detection_and_donor_domains():
    args = inputs()
    args[2].iloc[0, 0] = 1.1
    with pytest.raises(ValueError, match="fractions"):
        prepare_inputs(*args)


def test_simulation_helpers():
    assert wilson(0, 1000)[0] == pytest.approx(0)
    assert wilson(1000, 1000)[1] == pytest.approx(1)
    rng = np.random.default_rng(14)
    samples = np.stack([ring_noise(rng, 512, 0.8) for _ in range(500)])
    assert samples.var() == pytest.approx(1, abs=0.04)
    assert np.corrcoef(samples.ravel(), np.roll(samples, 1, axis=1).ravel())[0, 1] == pytest.approx(
        0.8, abs=0.02
    )
