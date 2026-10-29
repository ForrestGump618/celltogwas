"""Genomic interval primitives.

All intervals in this module use the 0-based, half-open convention ``[start,
end)``. GTF/GFF records are converted at input time. Keeping one coordinate
system avoids the off-by-one and boundary errors common in GWAS annotation
code.
"""

from __future__ import annotations

import heapq
import math
import numbers
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd


def normalize_chromosome(value: object) -> str:
    """Normalize chromosome labels without changing contig semantics.

    Examples: ``chr1`` -> ``1``, ``CHR X`` -> ``X``, ``chrM`` -> ``MT``.
    """

    if value is None:
        return ""
    try:
        if bool(pd.isna(value)):
            return ""
    except (TypeError, ValueError):
        pass
    if isinstance(value, numbers.Real) and not math.isfinite(float(value)):
        return ""
    if isinstance(value, numbers.Integral):
        return str(int(value))
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "<na>", "nat"}:
        return ""
    if text.lower().startswith("chr"):
        text = text[3:].strip()
    if text.endswith(".0") and text[:-2].isdigit():
        text = text[:-2]
    upper = text.upper()
    aliases = {
        "M": "MT",
        "MITO": "MT",
        "MITOCHONDRIA": "MT",
        "MITOCHONDRION": "MT",
    }
    return aliases.get(upper, upper)


def prefix_chromosome(chrom: object) -> str:
    """Return a normalized chromosome with a UCSC ``chr`` prefix."""

    normalized = normalize_chromosome(chrom)
    return f"chr{normalized}" if normalized else ""


def merge_intervals(intervals: Iterable[Sequence[int]]) -> list[list[int]]:
    """Merge overlapping or directly adjacent half-open intervals.

    Adjacent intervals are merged because their union has no gaps.
    """

    cleaned: list[tuple[int, int]] = []
    for interval in intervals:
        if len(interval) < 2:
            raise ValueError(f"Interval must contain start and end: {interval!r}")
        start, end = int(interval[0]), int(interval[1])
        if start < 0:
            raise ValueError(f"Interval start must be non-negative: {interval!r}")
        if end <= start:
            raise ValueError(f"Interval end must exceed start: {interval!r}")
        cleaned.append((start, end))
    if not cleaned:
        return []
    cleaned.sort()
    merged = [list(cleaned[0])]
    for start, end in cleaned[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def check_value_in_intervals(intervals: Sequence[Sequence[int]], value: int) -> bool:
    """Return whether a position lies in sorted, non-overlapping intervals."""

    if not intervals:
        return False
    starts = [int(interval[0]) for interval in intervals]
    index = int(np.searchsorted(starts, int(value), side="right")) - 1
    if index < 0:
        return False
    return int(intervals[index][0]) <= int(value) < int(intervals[index][1])


@dataclass(frozen=True)
class WeightedIntervalIndex:
    """Sorted, non-overlapping intervals with a maximum-coverage weight.

    Overlapping source intervals are resolved by taking the maximum active
    weight. This lets gene bodies, promoters and optional regulatory elements
    contribute to one interpretable cell-type score without double counting.
    """

    starts: np.ndarray
    ends: np.ndarray
    weights: np.ndarray
    chromosomes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        starts = np.asarray(self.starts, dtype=np.int64)
        ends = np.asarray(self.ends, dtype=np.int64)
        weights = np.asarray(self.weights, dtype=np.float64)
        if not (len(starts) == len(ends) == len(weights)):
            raise ValueError("Interval starts, ends and weights must have equal length.")
        if np.any(ends <= starts):
            raise ValueError("Every interval end must exceed its start.")
        if np.any(weights <= 0) or np.any(~np.isfinite(weights)):
            raise ValueError("Interval weights must be finite and positive.")
        if len(starts) > 1 and np.any(starts[1:] < ends[:-1]):
            raise ValueError("Intervals must be sorted, non-overlapping segments.")
        object.__setattr__(self, "starts", starts)
        object.__setattr__(self, "ends", ends)
        object.__setattr__(self, "weights", weights)
        object.__setattr__(self, "chromosomes", tuple(sorted(set(self.chromosomes))))

    @classmethod
    def empty(cls, chromosomes: Iterable[str] = ()) -> WeightedIntervalIndex:
        """Create an empty index."""

        return cls(
            starts=np.empty(0, dtype=np.int64),
            ends=np.empty(0, dtype=np.int64),
            weights=np.empty(0, dtype=np.float64),
            chromosomes=tuple(sorted(set(chromosomes))),
        )

    @classmethod
    def from_weighted_intervals(
        cls,
        intervals: Iterable[tuple[int, int, float]],
        chromosomes: Iterable[str] = (),
    ) -> WeightedIntervalIndex:
        """Build a maximum-coverage index from possibly overlapping intervals."""

        source: list[tuple[int, int, float]] = []
        for start_raw, end_raw, weight_raw in intervals:
            start, end = int(start_raw), int(end_raw)
            weight = float(weight_raw)
            if not math.isfinite(weight):
                raise ValueError("Interval weight must be finite.")
            if weight <= 0:
                continue
            if start < 0:
                raise ValueError("Interval start must be non-negative.")
            if end <= start:
                continue
            source.append((start, end, weight))
        if not source:
            return cls.empty(chromosomes)

        events: list[tuple[int, int, float]] = []
        for start, end, weight in source:
            events.append((start, 1, weight))
            events.append((end, -1, weight))
        events.sort(key=lambda event: event[0])

        active: dict[float, int] = {}
        max_heap: list[float] = []
        segments: list[tuple[int, int, float]] = []
        previous = events[0][0]
        cursor = 0
        while cursor < len(events):
            position = events[cursor][0]
            while max_heap and active.get(-max_heap[0], 0) <= 0:
                heapq.heappop(max_heap)
            if previous < position and max_heap:
                segments.append((previous, position, -max_heap[0]))
            while cursor < len(events) and events[cursor][0] == position:
                _, delta, weight = events[cursor]
                count = active.get(weight, 0) + delta
                if count == 0:
                    active.pop(weight, None)
                else:
                    active[weight] = count
                if delta > 0:
                    heapq.heappush(max_heap, -weight)
                cursor += 1
            previous = position

        merged: list[tuple[int, int, float]] = []
        for start, end, weight in segments:
            if merged and merged[-1][1] == start and merged[-1][2] == weight:
                merged[-1] = (merged[-1][0], end, weight)
            else:
                merged.append((start, end, weight))

        return cls(
            starts=np.asarray([segment[0] for segment in merged], dtype=np.int64),
            ends=np.asarray([segment[1] for segment in merged], dtype=np.int64),
            weights=np.asarray([segment[2] for segment in merged], dtype=np.float64),
            chromosomes=tuple(sorted(set(chromosomes))),
        )

    def __len__(self) -> int:
        return int(len(self.starts))

    @property
    def covered_bp(self) -> int:
        """Number of genomic bases represented by any positive-weight interval."""

        return int(np.sum(self.ends - self.starts, dtype=np.int64))

    @property
    def weighted_covered_bp(self) -> float:
        """Integral of interval weights over genomic coordinates."""

        return float(np.sum((self.ends - self.starts) * self.weights, dtype=np.float64))

    def query(self, positions: Sequence[int] | np.ndarray) -> np.ndarray:
        """Look up weights for 0-based genomic positions."""

        query = np.asarray(positions, dtype=np.int64)
        output = np.zeros(query.shape, dtype=np.float64)
        if len(self) == 0 or query.size == 0:
            return output
        indices = np.searchsorted(self.starts, query, side="right") - 1
        valid = indices >= 0
        safe_indices = np.where(valid, indices, 0)
        valid &= query < self.ends[safe_indices]
        output[valid] = self.weights[safe_indices[valid]]
        return output


def build_chromosome_indexes(
    records: pd.DataFrame,
    chrom_col: str = "chrom",
    start_col: str = "start",
    end_col: str = "end",
    weight_col: str = "weight",
) -> dict[str, WeightedIntervalIndex]:
    """Build one weighted interval index per chromosome."""

    indexes: dict[str, WeightedIntervalIndex] = {}
    if records.empty:
        return indexes
    for chrom, frame in records.groupby(chrom_col, sort=False):
        normalized = normalize_chromosome(chrom)
        if not normalized:
            continue
        records_for_chrom = frame.loc[:, [start_col, end_col, weight_col]].itertuples(
            index=False, name=None
        )
        indexes[normalized] = WeightedIntervalIndex.from_weighted_intervals(
            records_for_chrom,
            chromosomes=(normalized,),
        )
    return indexes
