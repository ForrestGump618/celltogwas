"""Backward-compatible functions from celltogwas 0.1.

These wrappers intentionally keep the historical dataframe shapes while using
the corrected parser and interval engine internally.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from .exceptions import InputFormatError
from .interval import (
    check_value_in_intervals,
    merge_intervals,
    normalize_chromosome,
)
from .io import load_gene_annotation


def extract_gene_info(gtf_file: str | Path) -> pd.DataFrame:
    """Return legacy ``GENE/CHR/START/END`` one-based gene records.

    Unlike the 0.1 implementation, this function retains gene-level records
    rather than matching the first transcript-like attribute block it sees.
    """

    return load_gene_annotation(gtf_file, feature_types=("gene",)).as_legacy_frame()


def celltype_genome_range(
    cell_gene_dict: Mapping[str, Sequence[str]],
    gene_info: pd.DataFrame,
    windowsize: int,
) -> dict[str, dict[str, list[list[int]]]]:
    """Build legacy 0-based, half-open merged ranges per cell type."""

    if windowsize < 0:
        raise InputFormatError("windowsize must be non-negative.")
    required = {"GENE", "CHR", "START", "END"}
    normalized = _normalize_gene_info_columns(gene_info)
    missing = required.difference(normalized.columns)
    if missing:
        raise InputFormatError(f"gene_info is missing columns: {sorted(missing)}")
    normalized["START"] = pd.to_numeric(normalized["START"], errors="coerce")
    normalized["END"] = pd.to_numeric(normalized["END"], errors="coerce")
    normalized = normalized.dropna(subset=["START", "END"])

    output: dict[str, dict[str, list[list[int]]]] = {}
    for cell_type, genes in cell_gene_dict.items():
        selected = normalized.loc[normalized["GENE"].isin(list(genes))]
        by_chrom: dict[str, list[tuple[int, int]]] = {}
        for record in selected.itertuples(index=False):
            chrom = normalize_chromosome(record.CHR)
            if not chrom:
                continue
            start = max(0, int(record.START) - int(windowsize) - 1)
            end = int(record.END) + int(windowsize)
            if end > start:
                by_chrom.setdefault(f"chr{chrom}", []).append((start, end))
        cell_range: dict[str, list[list[int]]] = {}
        for chrom, intervals in by_chrom.items():
            cell_range[chrom] = merge_intervals(intervals)
        output[str(cell_type)] = cell_range
    return output


def gwas_annotation(
    gwas: pd.DataFrame,
    celltype_range: Mapping[str, Mapping[str, Sequence[Sequence[int]]]],
) -> pd.DataFrame:
    """Annotate the first two GWAS columns with legacy string labels."""

    if not isinstance(gwas, pd.DataFrame) or gwas.shape[1] < 2:
        raise InputFormatError("gwas must contain at least chromosome and position columns.")
    chrom_values = gwas.iloc[:, 0].map(normalize_chromosome).to_numpy()
    positions_one_based = pd.to_numeric(gwas.iloc[:, 1], errors="coerce").to_numpy()
    positions_zero_based = positions_one_based - 1
    output = gwas.copy()
    for cell_type, chromosome_ranges in celltype_range.items():
        labels = np.full(len(gwas), "0", dtype=object)
        for chrom, intervals in chromosome_ranges.items():
            mask = chrom_values == normalize_chromosome(chrom)
            if not np.any(mask):
                continue
            labels[mask] = [
                "1" if _position_in_intervals(intervals, position) else "0"
                for position in positions_zero_based[mask]
            ]
        output[str(cell_type)] = labels
    return output


def _position_in_intervals(
    intervals: Sequence[Sequence[int]],
    position: float,
) -> bool:
    if not np.isfinite(position):
        return False
    return check_value_in_intervals(intervals, int(position))


def _normalize_gene_info_columns(frame: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(frame, pd.DataFrame):
        raise InputFormatError("gene_info must be a pandas DataFrame.")
    rename = {str(column).strip().upper(): column for column in frame.columns}
    return frame.rename(columns={source: target for target, source in rename.items()})


# Preserve the misspelling used in the original README and downstream scripts.
gwas_annotaion = gwas_annotation


def run() -> None:
    """Compatibility no-op retained for legacy imports."""

    return None
