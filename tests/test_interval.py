from __future__ import annotations

import numpy as np
import pandas as pd

from celltogwas.interval import (
    WeightedIntervalIndex,
    check_value_in_intervals,
    merge_intervals,
    normalize_chromosome,
)


def test_normalize_chromosome_and_merge_intervals():
    assert normalize_chromosome("chr1") == "1"
    assert normalize_chromosome("CHR M") == "MT"
    assert normalize_chromosome(7) == "7"
    assert normalize_chromosome(pd.NA) == ""
    assert normalize_chromosome(float("inf")) == ""
    assert merge_intervals([(0, 5), (5, 10), (20, 30), (25, 28)]) == [
        [0, 10],
        [20, 30],
    ]
    assert check_value_in_intervals([[0, 10], [20, 30]], 10) is False
    assert check_value_in_intervals([[0, 10], [20, 30]], 25) is True


def test_weighted_interval_index_uses_maximum_active_weight():
    index = WeightedIntervalIndex.from_weighted_intervals(
        [(0, 10, 0.5), (5, 15, 1.0), (8, 12, 0.25)]
    )

    query = index.query(np.array([0, 4, 5, 14, 15]))

    np.testing.assert_allclose(query, [0.5, 0.5, 1.0, 1.0, 0.0])
    assert index.covered_bp == 15
    assert index.weighted_covered_bp == 12.5
