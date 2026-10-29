"""Regulatory-aware SNP-to-cell-type mapping."""

from __future__ import annotations

import math
import numbers
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

from ._version import __version__
from .exceptions import InputFormatError, LinkageError
from .interval import (
    WeightedIntervalIndex,
    normalize_chromosome,
)
from .io import (
    load_gene_annotation,
    load_gene_sets,
    load_regulatory_bed,
    resolve_column,
)
from .models import CellTypeBuildDiagnostic, GeneAnnotation, LinkConfig

AnnotationMode = Literal["weight", "binary", "both"]


@dataclass
class GenomeLinker:
    """Map GWAS positions to weighted cell-type genomic domains."""

    annotation: GeneAnnotation
    gene_sets: dict[str, dict[str, float]]
    config: LinkConfig
    indexes: dict[str, dict[str, WeightedIntervalIndex]]
    diagnostics: pd.DataFrame
    aliases: dict[str, list[int]]
    regulatory: pd.DataFrame | None = None

    @classmethod
    def from_paths(
        cls,
        annotation_path: str | Path,
        gene_sets_path: str | Path,
        *,
        config: LinkConfig | None = None,
        regulatory_path: str | Path | None = None,
        regulatory_cell_type: str | None = None,
        feature_types: Sequence[str] | None = ("gene",),
        strict_annotation: bool = False,
    ) -> GenomeLinker:
        """Build a linker from a GTF/GFF file and weighted gene sets."""

        annotation = load_gene_annotation(
            annotation_path,
            feature_types=feature_types,
            strict=strict_annotation,
        )
        gene_sets = load_gene_sets(gene_sets_path)
        regulatory = (
            load_regulatory_bed(regulatory_path, cell_type=regulatory_cell_type)
            if regulatory_path is not None
            else None
        )
        return cls.from_annotation(
            annotation,
            gene_sets,
            config=config,
            regulatory=regulatory,
        )

    @classmethod
    def from_annotation(
        cls,
        annotation: GeneAnnotation,
        gene_sets: Mapping[str, Mapping[str, float] | Sequence[str]],
        *,
        config: LinkConfig | None = None,
        regulatory: pd.DataFrame | None = None,
    ) -> GenomeLinker:
        """Build a linker from in-memory annotation and gene sets."""

        active_config = config or LinkConfig()
        normalized_sets = _normalize_gene_sets(gene_sets)
        normalized_regulatory = _normalize_regulatory(regulatory)
        if normalized_regulatory is not None:
            for cell_type in normalized_regulatory["cell_type"].unique():
                normalized_sets.setdefault(str(cell_type), {})
        if not normalized_sets:
            raise LinkageError("At least one non-empty gene set is required.")

        aliases = _build_alias_map(annotation.genes)
        indexes: dict[str, dict[str, WeightedIntervalIndex]] = {}
        diagnostics: list[CellTypeBuildDiagnostic] = []
        for cell_type, weights in normalized_sets.items():
            cell_regulatory = None
            if normalized_regulatory is not None:
                cell_regulatory = normalized_regulatory.loc[
                    normalized_regulatory["cell_type"] == cell_type
                ]
            cell_index, diagnostic = _build_cell_type_index(
                cell_type,
                weights,
                annotation,
                aliases,
                active_config,
                cell_regulatory,
            )
            indexes[cell_type] = cell_index
            diagnostics.append(diagnostic)

        return cls(
            annotation=annotation,
            gene_sets=normalized_sets,
            config=active_config,
            indexes=indexes,
            diagnostics=pd.DataFrame([diagnostic.as_dict() for diagnostic in diagnostics]),
            aliases=aliases,
            regulatory=normalized_regulatory,
        )

    @property
    def cell_types(self) -> tuple[str, ...]:
        """Cell types represented by the linker."""

        return tuple(self.indexes)

    def annotate(
        self,
        gwas: pd.DataFrame,
        *,
        mode: AnnotationMode = "weight",
        threshold: float = 0.05,
        chunk_size: int = 500_000,
    ) -> pd.DataFrame:
        """Annotate GWAS rows with cell-type link weights or binary labels.

        ``mode='weight'`` adds one continuous [0, 1] score per cell type.
        ``mode='binary'`` adds 0/1 columns. ``mode='both'`` adds the binary
        column under the cell-type name and its continuous score as
        ``<cell_type>_weight``.
        """

        if mode not in {"weight", "binary", "both"}:
            raise InputFormatError("mode must be 'weight', 'binary', or 'both'.")
        if not math.isfinite(threshold) or not 0 <= threshold <= self.config.max_weight:
            raise InputFormatError("threshold must lie between 0 and LinkConfig.max_weight.")
        if isinstance(chunk_size, bool) or not isinstance(chunk_size, numbers.Integral):
            raise InputFormatError("chunk_size must be an integer.")
        if chunk_size < 1:
            raise InputFormatError("chunk_size must be positive.")
        if not isinstance(gwas, pd.DataFrame) or gwas.empty:
            raise InputFormatError("gwas must be a non-empty DataFrame.")

        chrom_column = resolve_column(
            gwas.columns,
            ("chr", "chrom", "chromosome", "chr_id", "contig"),
            label="chromosome",
        )
        position_column = resolve_column(
            gwas.columns,
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
        chromosome = gwas[chrom_column].map(normalize_chromosome)
        position = pd.to_numeric(gwas[position_column], errors="coerce")
        position_zero_based = position - 1.0

        output_names = _annotation_output_names(self.cell_types, mode)
        if len(set(output_names)) != len(output_names):
            raise InputFormatError(
                "Cell-type names produce duplicate output columns in mode='both'. "
                "Rename names ending in '_weight'."
            )
        collisions = set(output_names).intersection(map(str, gwas.columns))
        if collisions:
            raise InputFormatError(
                "Cell-type annotation names collide with GWAS columns: "
                f"{sorted(collisions)}. Rename the GWAS columns before annotation."
            )

        chunks: list[pd.DataFrame] = []
        for start in range(0, len(gwas), chunk_size):
            stop = min(start + chunk_size, len(gwas))
            chunk = gwas.iloc[start:stop].copy()
            chunk_chrom = chromosome.iloc[start:stop]
            chunk_position = position_zero_based.iloc[start:stop]
            valid_position = (
                chunk_position.notna()
                & np.isfinite(chunk_position)
                & (chunk_position >= 0)
                & np.equal(chunk_position, np.floor(chunk_position))
            )
            for cell_type in self.cell_types:
                values = np.zeros(stop - start, dtype=np.float64)
                cell_indexes = self.indexes[cell_type]
                for chrom in chunk_chrom.loc[valid_position].unique():
                    index = cell_indexes.get(str(chrom))
                    if index is None or len(index) == 0:
                        continue
                    mask = valid_position & (chunk_chrom == chrom)
                    values[mask.to_numpy()] = index.query(
                        chunk_position.loc[mask].to_numpy(dtype=np.int64)
                    )
                if mode == "weight":
                    chunk[cell_type] = values
                elif mode == "binary":
                    chunk[cell_type] = (values >= max(threshold, np.finfo(float).eps)).astype(
                        np.int8
                    )
                else:
                    chunk[cell_type] = (values >= max(threshold, np.finfo(float).eps)).astype(
                        np.int8
                    )
                    chunk[f"{cell_type}_weight"] = values
            chunks.append(chunk)
        return pd.concat(chunks, axis=0)

    def enrich(
        self,
        gwas: pd.DataFrame,
        *,
        p_column: str | None = None,
        score_column: str | None = None,
        score_mode: str = "pvalue",
        n_permutations: int = 10_000,
        random_state: int = 13,
        min_shift_fraction: float = 0.01,
        p_upper: float | None = None,
        weight_floor: float = 0.0,
        chunk_size: int = 500_000,
    ) -> pd.DataFrame:
        """Run spatial circular-permutation enrichment for each cell type."""

        from .enrichment import circular_permutation_enrichment

        weights = self.annotate(
            gwas,
            mode="weight",
            threshold=0.0,
            chunk_size=chunk_size,
        )
        return circular_permutation_enrichment(
            weights,
            cell_types=list(self.cell_types),
            p_column=p_column,
            score_column=score_column,
            score_mode=score_mode,
            n_permutations=n_permutations,
            random_state=random_state,
            min_shift_fraction=min_shift_fraction,
            p_upper=p_upper,
            weight_floor=weight_floor,
        )

    def domain_summary(self) -> pd.DataFrame:
        """Summarize merged genomic coverage for each cell type."""

        rows: list[dict[str, object]] = []
        for cell_type, chromosome_indexes in self.indexes.items():
            n_intervals = sum(len(index) for index in chromosome_indexes.values())
            covered_bp = sum(index.covered_bp for index in chromosome_indexes.values())
            weighted_bp = sum(index.weighted_covered_bp for index in chromosome_indexes.values())
            rows.append(
                {
                    "cell_type": cell_type,
                    "n_intervals": n_intervals,
                    "covered_bp": covered_bp,
                    "weighted_covered_bp": weighted_bp,
                }
            )
        return pd.DataFrame(rows)

    def manifest(self) -> dict[str, object]:
        """Return a reproducible JSON-compatible build manifest."""

        return {
            "package": "celltogwas",
            "version": __version__,
            "annotation": self.annotation.summary(),
            "config": self.config.as_dict(),
            "cell_types": list(self.cell_types),
            "gene_set_sizes": {
                cell_type: len(genes) for cell_type, genes in self.gene_sets.items()
            },
            "diagnostics": self.diagnostics.to_dict(orient="records"),
            "domain_summary": self.domain_summary().to_dict(orient="records"),
            "statistical_scope": (
                "Cell-type domain annotation and spatial permutation enrichment; "
                "not a heritability partition or causal fine-mapping method."
            ),
        }


def _normalize_gene_sets(
    gene_sets: Mapping[str, Mapping[str, float] | Sequence[str]],
) -> dict[str, dict[str, float]]:
    output: dict[str, dict[str, float]] = {}
    for cell_type, genes in gene_sets.items():
        name = str(cell_type).strip()
        if not name:
            raise LinkageError("Cell-type names cannot be empty.")
        if isinstance(genes, str):
            raise LinkageError(
                f"Gene set {name!r} must be a mapping or a sequence of genes, not a single string."
            )
        items = genes.items() if isinstance(genes, Mapping) else ((gene, 1.0) for gene in genes)
        weights: dict[str, float] = {}
        for gene, weight in items:
            gene_name = str(gene).strip()
            value = float(weight)
            if not gene_name:
                continue
            if not math.isfinite(value) or value < 0:
                raise LinkageError(f"Invalid weight {value!r} for {gene_name!r} in {name!r}.")
            if value > 0:
                weights[gene_name] = max(weights.get(gene_name, 0.0), value)
        if weights:
            output[name] = weights
    return output


def _normalize_regulatory(regulatory: pd.DataFrame | None) -> pd.DataFrame | None:
    if regulatory is None:
        return None
    required = {"chrom", "start", "end", "cell_type", "weight"}
    missing = required.difference(regulatory.columns)
    if missing:
        raise LinkageError(f"Regulatory intervals are missing columns: {sorted(missing)}")
    frame = regulatory.loc[:, ["chrom", "start", "end", "cell_type", "weight"]].copy()
    frame["chrom"] = frame["chrom"].map(normalize_chromosome)
    frame["start"] = pd.to_numeric(frame["start"], errors="coerce")
    frame["end"] = pd.to_numeric(frame["end"], errors="coerce")
    frame["weight"] = pd.to_numeric(frame["weight"], errors="coerce")
    frame = frame.dropna(subset=["chrom", "start", "end", "weight"])
    frame = frame.loc[
        (frame["chrom"] != "")
        & (frame["start"] >= 0)
        & np.isfinite(frame["start"])
        & (frame["end"] > frame["start"])
        & np.isfinite(frame["end"])
        & (frame["weight"] > 0)
        & np.isfinite(frame["weight"])
    ].copy()
    if frame.empty:
        return None
    frame["start"] = frame["start"].astype(np.int64)
    frame["end"] = frame["end"].astype(np.int64)
    frame["weight"] = frame["weight"].astype(float)
    frame["cell_type"] = frame["cell_type"].astype(str).str.strip()
    return frame.loc[frame["cell_type"] != ""].reset_index(drop=True)


def _build_alias_map(genes: pd.DataFrame) -> dict[str, list[int]]:
    aliases: dict[str, list[int]] = defaultdict(list)
    for row_index, (gene_id, gene_name) in enumerate(
        zip(genes["gene_id"], genes["gene_name"], strict=False)
    ):
        for alias in {str(gene_id).strip(), str(gene_name).strip()}:
            if not alias:
                continue
            aliases[alias].append(row_index)
            aliases[alias.upper()].append(row_index)
    for alias, indices in aliases.items():
        aliases[alias] = sorted(set(indices))
    return dict(aliases)


def _lookup_gene_rows(aliases: Mapping[str, list[int]], gene: str) -> list[int]:
    return list(aliases.get(gene, aliases.get(gene.upper(), ())))


def _clip_weight(weight: float, config: LinkConfig) -> float:
    return max(0.0, min(float(config.max_weight), float(weight)))


def _window_intervals(
    start: int,
    end: int,
    gene_weight: float,
    config: LinkConfig,
) -> list[tuple[int, int, float]]:
    if config.window == 0:
        return []
    if config.kernel == "uniform":
        output: list[tuple[int, int, float]] = []
        left_start = max(0, start - config.window)
        if start > left_start:
            output.append((left_start, start, gene_weight))
        output.append((end, end + config.window, gene_weight))
        return output

    scale = config.effective_decay_scale
    output: list[tuple[int, int, float]] = []
    for level in range(config.kernel_levels):
        inner = config.window * level / config.kernel_levels
        outer = config.window * (level + 1) / config.kernel_levels
        midpoint = (inner + outer) / 2.0
        band_weight = gene_weight * math.exp(-midpoint / scale)
        left_start = max(0, int(math.floor(start - outer)))
        left_end = max(0, int(math.ceil(start - inner)))
        if left_end > left_start:
            output.append((left_start, left_end, band_weight))
        right_start = int(math.floor(end + inner))
        right_end = int(math.ceil(end + outer))
        if right_end > right_start:
            output.append((right_start, right_end, band_weight))
    return output


def _promoter_interval(start: int, end: int, strand: str, radius: int) -> tuple[int, int]:
    if radius == 0:
        return (0, 0)
    tss = start if strand != "-" else max(start, end - 1)
    return (max(0, tss - radius), tss + radius + 1)


def _build_cell_type_index(
    cell_type: str,
    gene_weights: Mapping[str, float],
    annotation: GeneAnnotation,
    aliases: Mapping[str, list[int]],
    config: LinkConfig,
    regulatory: pd.DataFrame | None,
) -> tuple[dict[str, WeightedIntervalIndex], CellTypeBuildDiagnostic]:
    intervals_by_chrom: dict[str, list[tuple[int, int, float]]] = defaultdict(list)
    matched_genes: list[str] = []
    missing_genes: list[str] = []
    ambiguous_genes: list[str] = []
    resolved_records: dict[int, float] = {}
    source_intervals = 0

    for gene, gene_weight in gene_weights.items():
        row_indices = _lookup_gene_rows(aliases, gene)
        if not row_indices:
            missing_genes.append(gene)
            continue
        matched_genes.append(gene)
        if len(row_indices) > 1:
            ambiguous_genes.append(gene)
        for row_index in row_indices:
            previous = resolved_records.get(row_index)
            if previous is None or gene_weight > previous:
                resolved_records[row_index] = gene_weight

    for row_index, gene_weight in resolved_records.items():
        record = annotation.genes.iloc[row_index]
        chrom = str(record["chrom"])
        start = int(record["start"])
        end = int(record["end"])
        body_weight = _clip_weight(gene_weight * config.gene_body_weight, config)
        if body_weight > 0:
            intervals_by_chrom[chrom].append((start, end, body_weight))
            source_intervals += 1

        promoter_weight = _clip_weight(
            gene_weight * config.promoter_weight,
            config,
        )
        if promoter_weight > 0 and config.promoter_window > 0:
            promoter_start, promoter_end = _promoter_interval(
                start,
                end,
                str(record["strand"]),
                config.promoter_window,
            )
            if promoter_end > promoter_start:
                intervals_by_chrom[chrom].append((promoter_start, promoter_end, promoter_weight))
                source_intervals += 1

        window_intervals = _window_intervals(
            start,
            end,
            _clip_weight(gene_weight, config),
            config,
        )
        intervals_by_chrom[chrom].extend(window_intervals)
        source_intervals += len(window_intervals)

    if regulatory is not None:
        for record in regulatory.itertuples(index=False):
            intervals_by_chrom[str(record.chrom)].append(
                (
                    int(record.start),
                    int(record.end),
                    _clip_weight(float(record.weight), config),
                )
            )
            source_intervals += 1

    indexes: dict[str, WeightedIntervalIndex] = {}
    for chrom, intervals in intervals_by_chrom.items():
        index = WeightedIntervalIndex.from_weighted_intervals(
            intervals,
            chromosomes=(chrom,),
        )
        if len(index) > 0:
            indexes[chrom] = index
    if not indexes:
        raise LinkageError(
            f"Cell type {cell_type!r} has no intervals after gene matching and weight filtering."
        )
    diagnostic = CellTypeBuildDiagnostic(
        cell_type=cell_type,
        requested_genes=len(gene_weights),
        matched_genes=len(matched_genes),
        matched_gene_records=len(resolved_records),
        source_intervals=source_intervals,
        merged_intervals=sum(len(index) for index in indexes.values()),
        covered_bp=sum(index.covered_bp for index in indexes.values()),
        weighted_covered_bp=sum(index.weighted_covered_bp for index in indexes.values()),
        missing_genes=sorted(missing_genes),
        ambiguous_genes=sorted(ambiguous_genes),
    )
    return indexes, diagnostic


def _annotation_output_names(
    cell_types: Sequence[str],
    mode: AnnotationMode,
) -> list[str]:
    if mode == "both":
        return [name for cell_type in cell_types for name in (cell_type, f"{cell_type}_weight")]
    return list(cell_types)
