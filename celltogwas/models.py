"""Shared data models and configuration objects."""

from __future__ import annotations

import math
import numbers
from dataclasses import asdict, dataclass, field
from typing import Literal

import numpy as np
import pandas as pd

from .exceptions import LinkageError
from .interval import normalize_chromosome

KernelName = Literal["uniform", "exponential"]


@dataclass(frozen=True)
class GeneAnnotation:
    """Gene records normalized to the 0-based, half-open UCSC convention.

    Required columns are ``gene_id``, ``gene_name``, ``chrom``, ``start``,
    ``end``, ``strand`` and ``biotype``. Coordinates are internal, so callers
    must not treat ``start`` as the one-based GTF start.
    """

    genes: pd.DataFrame
    source_path: str | None = None
    source_sha256: str | None = None

    def __post_init__(self) -> None:
        required = {
            "gene_id",
            "gene_name",
            "chrom",
            "start",
            "end",
            "strand",
            "biotype",
        }
        missing = required.difference(self.genes.columns)
        if missing:
            raise ValueError(f"Gene annotation is missing columns: {sorted(missing)}")
        if self.genes.empty:
            raise ValueError("Gene annotation contains no records.")
        frame = self.genes.copy().reset_index(drop=True)
        start_numeric = pd.to_numeric(frame["start"], errors="raise")
        end_numeric = pd.to_numeric(frame["end"], errors="raise")
        start_values = start_numeric.to_numpy(dtype=float)
        end_values = end_numeric.to_numpy(dtype=float)
        if not (
            np.isfinite(start_values).all()
            and np.isfinite(end_values).all()
            and np.equal(start_values, np.floor(start_values)).all()
            and np.equal(end_values, np.floor(end_values)).all()
        ):
            raise ValueError("Gene coordinates must be finite integer base positions.")
        frame["start"] = start_numeric.astype("int64")
        frame["end"] = end_numeric.astype("int64")
        if (frame["start"] < 0).any() or (frame["end"] <= frame["start"]).any():
            raise ValueError("Gene coordinates must be non-negative, 0-based, half-open intervals.")
        frame["chrom"] = frame["chrom"].map(normalize_chromosome)
        if (frame["chrom"] == "").any():
            raise ValueError("Gene annotation contains an empty chromosome.")
        for column in ("gene_id", "gene_name", "strand", "biotype"):
            frame[column] = frame[column].astype(str)
        object.__setattr__(self, "genes", frame)

    @property
    def n_genes(self) -> int:
        """Number of normalized gene records."""

        return int(len(self.genes))

    @property
    def chromosomes(self) -> tuple[str, ...]:
        """Chromosomes represented in the annotation."""

        return tuple(sorted(self.genes["chrom"].astype(str).unique()))

    def as_legacy_frame(self) -> pd.DataFrame:
        """Return one-based coordinates using the historical column names."""

        output = self.genes.loc[:, ["gene_name", "chrom", "start", "end"]].copy()
        output["start"] = output["start"].astype(int) + 1
        output["end"] = output["end"].astype(int)
        return output.rename(
            columns={
                "gene_name": "GENE",
                "chrom": "CHR",
                "start": "START",
                "end": "END",
            }
        )

    def summary(self) -> dict[str, object]:
        """Return a JSON-serializable annotation summary."""

        return {
            "source_path": self.source_path,
            "source_sha256": self.source_sha256,
            "n_genes": self.n_genes,
            "n_unique_gene_ids": int(self.genes["gene_id"].nunique()),
            "n_unique_gene_names": int(self.genes["gene_name"].nunique()),
            "n_chromosomes": len(self.chromosomes),
            "chromosomes": list(self.chromosomes),
        }


@dataclass(frozen=True)
class LinkConfig:
    """Configuration for mapping a SNP to cell-type-linked genomic domains.

    The exponential kernel is evaluated in discrete distance bands. This keeps
    interval lookup logarithmic while retaining a distance-decayed SNP-to-gene
    score. Set ``kernel='uniform'`` for the historical hard window.
    """

    window: int = 100_000
    promoter_window: int = 2_000
    promoter_weight: float = 1.0
    gene_body_weight: float = 1.0
    kernel: KernelName = "exponential"
    decay_scale: int | None = None
    kernel_levels: int = 6
    max_weight: float = 1.0

    def __post_init__(self) -> None:
        for name, value in (
            ("window", self.window),
            ("promoter_window", self.promoter_window),
            ("kernel_levels", self.kernel_levels),
        ):
            if isinstance(value, bool) or not isinstance(value, numbers.Integral):
                raise LinkageError(f"{name} must be an integer.")
        if self.window < 0:
            raise LinkageError("window must be non-negative.")
        if self.promoter_window < 0:
            raise LinkageError("promoter_window must be non-negative.")
        if self.kernel not in {"uniform", "exponential"}:
            raise LinkageError("kernel must be 'uniform' or 'exponential'.")
        if self.kernel_levels < 1:
            raise LinkageError("kernel_levels must be at least 1.")
        if self.decay_scale is not None:
            if isinstance(self.decay_scale, bool) or not isinstance(
                self.decay_scale, numbers.Integral
            ):
                raise LinkageError("decay_scale must be an integer when provided.")
            if self.decay_scale <= 0:
                raise LinkageError("decay_scale must be positive when provided.")
        for name, value in (
            ("promoter_weight", self.promoter_weight),
            ("gene_body_weight", self.gene_body_weight),
            ("max_weight", self.max_weight),
        ):
            if not math.isfinite(float(value)):
                raise LinkageError(f"{name} must be finite.")
            if value < 0:
                raise LinkageError(f"{name} must be non-negative.")
        if self.max_weight == 0:
            raise LinkageError("max_weight must be positive.")

    @property
    def effective_decay_scale(self) -> float:
        """Decay constant used by the exponential kernel."""

        if self.decay_scale is not None:
            return float(self.decay_scale)
        return max(1.0, self.window / 3.0)

    def as_dict(self) -> dict[str, object]:
        """Return a JSON-serializable configuration."""

        return asdict(self)


@dataclass
class CellTypeBuildDiagnostic:
    """Diagnostics collected while building one cell-type domain."""

    cell_type: str
    requested_genes: int
    matched_genes: int
    matched_gene_records: int
    source_intervals: int
    merged_intervals: int
    covered_bp: int
    weighted_covered_bp: float
    missing_genes: list[str] = field(default_factory=list)
    ambiguous_genes: list[str] = field(default_factory=list)

    @property
    def match_fraction(self) -> float:
        """Fraction of requested gene symbols represented in the annotation."""

        if self.requested_genes == 0:
            return 0.0
        return self.matched_genes / self.requested_genes

    def as_dict(self) -> dict[str, object]:
        """Return a JSON-serializable diagnostic."""

        payload = asdict(self)
        payload["match_fraction"] = self.match_fraction
        return payload
