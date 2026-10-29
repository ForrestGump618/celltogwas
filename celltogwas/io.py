"""Input/output helpers for annotations, gene sets and tabular data."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import math
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pandas as pd

from .exceptions import AnnotationError, InputFormatError
from .interval import normalize_chromosome
from .models import GeneAnnotation


def _path_string(path: str | Path) -> str:
    return str(Path(path).expanduser().resolve())


def _is_gzip(path: str | Path) -> bool:
    with Path(path).open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def open_text(path: str | Path) -> io.TextIOBase:
    """Open plain text or gzip-compressed text by file signature."""

    path = Path(path)
    if _is_gzip(path):
        return gzip.open(path, mode="rt", encoding="utf-8", errors="replace")
    return path.open(mode="rt", encoding="utf-8", errors="replace")


def file_sha256(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Compute a SHA-256 digest for reproducibility manifests."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def parse_attributes(attributes: str, *, gff3: bool | None = None) -> dict[str, str]:
    """Parse GTF2 or GFF3 attribute fields into a dictionary."""

    if gff3 is None:
        gff3 = "=" in attributes and ' "' not in attributes
    if gff3:
        parsed: dict[str, str] = {}
        for token in attributes.strip().split(";"):
            if not token.strip():
                continue
            key, separator, value = token.strip().partition("=")
            if separator:
                parsed[key.strip()] = unquote(value.strip())
        return parsed

    parsed = {}
    for token in attributes.split(";"):
        token = token.strip()
        if not token:
            continue
        key, separator, value = token.partition(" ")
        if not separator:
            key, separator, value = token.partition("\t")
        if not separator:
            continue
        value = value.strip()
        if len(value) >= 2 and value[0] == value[-1] == '"':
            value = value[1:-1]
        parsed[key.strip()] = value
    return parsed


def load_gene_annotation(
    path: str | Path,
    *,
    feature_types: Sequence[str] | None = ("gene",),
    strict: bool = False,
) -> GeneAnnotation:
    """Load gene records from a GTF/GFF3 or gzip-compressed equivalent.

    Parameters
    ----------
    feature_types:
        Feature types retained from column 3. The default ``("gene",)`` avoids
        the historical bug of counting one row per transcript. Use ``None`` to
        retain all records and aggregate them by ``gene_id``.
    strict:
        Raise on malformed non-comment rows instead of warning and skipping.
    """

    path = Path(path).expanduser()
    rows: list[dict[str, object]] = []
    skipped = 0
    retained_features = None if feature_types is None else set(feature_types)
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n\r")
            if not line or line.startswith("#"):
                if line == "##FASTA":
                    break
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                skipped += 1
                if strict:
                    raise AnnotationError(
                        f"{path}:{line_number} has {len(fields)} columns; expected 9."
                    )
                continue
            feature = fields[2]
            if retained_features is not None and feature not in retained_features:
                continue
            try:
                start_one_based = int(fields[3])
                end_one_based = int(fields[4])
            except ValueError:
                skipped += 1
                if strict:
                    raise AnnotationError(
                        f"{path}:{line_number} has non-integer coordinates."
                    ) from None
                continue
            if start_one_based < 1 or end_one_based < start_one_based:
                skipped += 1
                if strict:
                    raise AnnotationError(
                        f"{path}:{line_number} has invalid coordinates "
                        f"{start_one_based}-{end_one_based}."
                    )
                continue
            attributes = parse_attributes(fields[8])
            gene_id = (
                attributes.get("gene_id")
                or attributes.get("ID")
                or attributes.get("gene")
                or attributes.get("Name")
            )
            if not gene_id:
                skipped += 1
                if strict:
                    raise AnnotationError(f"{path}:{line_number} has no gene_id/ID attribute.")
                continue
            gene_name = (
                attributes.get("gene_name")
                or attributes.get("Name")
                or attributes.get("gene")
                or gene_id
            )
            biotype = (
                attributes.get("gene_type")
                or attributes.get("gene_biotype")
                or attributes.get("biotype")
                or ""
            )
            chrom = normalize_chromosome(fields[0])
            if not chrom:
                skipped += 1
                if strict:
                    raise AnnotationError(f"{path}:{line_number} has an empty chromosome.")
                continue
            rows.append(
                {
                    "gene_id": str(gene_id),
                    "gene_name": str(gene_name),
                    "chrom": chrom,
                    "start": start_one_based - 1,
                    "end": end_one_based,
                    "strand": fields[6] if fields[6] in {"+", "-", "."} else ".",
                    "biotype": str(biotype),
                }
            )

    if skipped:
        warnings.warn(
            f"Skipped {skipped} malformed or non-gene GTF/GFF rows from {path}.",
            RuntimeWarning,
            stacklevel=2,
        )
    if not rows:
        raise AnnotationError(
            f"No gene records were parsed from {path}. Check feature_types and attribute naming."
        )

    genes = pd.DataFrame.from_records(rows)
    # Aggregate transcript-level fallback records without losing multi-locus genes.
    genes = (
        genes.groupby(["gene_id", "chrom", "strand"], as_index=False, sort=False)
        .agg(
            gene_name=("gene_name", "first"),
            start=("start", "min"),
            end=("end", "max"),
            biotype=("biotype", "first"),
        )
        .loc[:, ["gene_id", "gene_name", "chrom", "start", "end", "strand", "biotype"]]
    )
    genes["start"] = genes["start"].astype(np.int64)
    genes["end"] = genes["end"].astype(np.int64)
    return GeneAnnotation(
        genes=genes.reset_index(drop=True),
        source_path=_path_string(path),
        source_sha256=file_sha256(path),
    )


def _normalize_gene_weights(
    values: Mapping[str, float] | Sequence[str],
) -> dict[str, float]:
    if isinstance(values, str):
        raise InputFormatError("A gene set must be a mapping or sequence, not a single string.")
    items = values.items() if isinstance(values, Mapping) else ((gene, 1.0) for gene in values)
    output: dict[str, float] = {}
    for gene, weight in items:
        gene_text = str(gene).strip()
        if not gene_text:
            continue
        try:
            numeric = float(weight)
        except (TypeError, ValueError) as error:
            raise InputFormatError(f"Invalid gene weight {weight!r} for {gene_text!r}.") from error
        if not math.isfinite(numeric):
            raise InputFormatError(f"Non-finite gene weight for {gene_text!r}.")
        if numeric < 0:
            raise InputFormatError(f"Negative gene weight for {gene_text!r}.")
        if numeric > 0:
            output[gene_text] = max(output.get(gene_text, 0.0), numeric)
    if not output:
        raise InputFormatError("A gene set cannot be empty.")
    return output


def _read_gene_sets_gmt(path: Path) -> dict[str, dict[str, float]]:
    output: dict[str, dict[str, float]] = {}
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\n\r").split("\t")
            if not fields or not fields[0].strip():
                continue
            name = fields[0].strip()
            genes = fields[2:] if len(fields) >= 2 else []
            if not genes:
                raise InputFormatError(f"{path}:{line_number} contains no genes.")
            weighted: dict[str, float] = {}
            for token in genes:
                gene, separator, weight_text = token.rpartition(":")
                if separator and gene and _is_float(weight_text):
                    numeric = float(weight_text)
                    if numeric < 0:
                        raise InputFormatError(f"{path}:{line_number} has a negative gene weight.")
                    weighted[gene] = max(weighted.get(gene, 0.0), numeric)
                else:
                    weighted[token] = max(weighted.get(token, 0.0), 1.0)
            output[name] = _normalize_gene_weights(weighted)
    if not output:
        raise InputFormatError(f"No gene sets were parsed from {path}.")
    return output


def _is_float(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def read_table(path: str | Path, **kwargs: object) -> pd.DataFrame:
    """Read CSV/TSV (optionally gzip) or parquet data."""

    path = Path(path)
    suffixes = [suffix.lower() for suffix in path.suffixes]
    active_suffix = suffixes[-2] if suffixes[-1:] == [".gz"] and len(suffixes) > 1 else ""
    if not active_suffix and suffixes:
        active_suffix = suffixes[-1]
    if active_suffix in {".csv", ".tsv", ".txt", ".tab"}:
        separator: str | None
        if active_suffix in {".tsv", ".tab"}:
            separator = "\t"
        elif active_suffix == ".txt":
            separator = None
        else:
            separator = ","
        try:
            if separator is None:
                return pd.read_csv(path, sep=None, engine="python", **kwargs)
            return pd.read_csv(path, sep=separator, **kwargs)
        except Exception as error:
            raise InputFormatError(f"Could not read table {path}: {error}") from error
    if active_suffix == ".parquet":
        try:
            return pd.read_parquet(path, **kwargs)
        except ImportError as error:
            raise InputFormatError("Reading parquet requires pyarrow or fastparquet.") from error
    # Delimited files without a conventional extension.
    try:
        return pd.read_csv(path, sep=None, engine="python", **kwargs)
    except Exception as error:
        raise InputFormatError(f"Could not read table {path}: {error}") from error


def load_gene_sets(
    path: str | Path,
    *,
    cell_type_column: str | None = None,
    gene_column: str | None = None,
    weight_column: str | None = None,
) -> dict[str, dict[str, float]]:
    """Load weighted cell-type gene sets from GMT, JSON, CSV/TSV.

    Long tables require cell-type and gene columns; a weight column is
    optional. Wide tables are interpreted as genes in the first column and one
    numeric cell-type column per remaining column.
    """

    path = Path(path)
    name_lower = path.name.lower()
    if name_lower.endswith(".gmt") or name_lower.endswith(".gmt.gz"):
        return _read_gene_sets_gmt(path)
    if name_lower.endswith(".json") or name_lower.endswith(".json.gz"):
        with open_text(path) as handle:
            payload = json.load(handle)
        if not isinstance(payload, Mapping):
            raise InputFormatError("Gene-set JSON must be an object keyed by cell type.")
        return {
            str(cell_type): _normalize_gene_weights(values) for cell_type, values in payload.items()
        }

    frame = read_table(path)
    if cell_type_column is None:
        cell_type_column = _find_column(
            frame.columns, ("cell_type", "celltype", "cell", "celltype_name")
        )
    if gene_column is None:
        gene_column = _find_column(
            frame.columns, ("gene", "gene_name", "gene_symbol", "symbol", "marker")
        )

    output: dict[str, dict[str, float]] = {}
    if cell_type_column is not None and gene_column is not None:
        if weight_column is None:
            weight_column = _find_column(
                frame.columns, ("weight", "score", "specificity", "gene_weight")
            )
        if weight_column is None:
            frame = frame.assign(__celltogwas_weight=1.0)
            weight_column = "__celltogwas_weight"
        for cell_type, group in frame.groupby(cell_type_column, sort=False):
            values: dict[str, float] = {}
            for gene, weight in zip(
                group[gene_column].astype(str),
                pd.to_numeric(group[weight_column], errors="coerce"),
                strict=False,
            ):
                if weight is None or pd.isna(weight):
                    continue
                if not _is_float(weight):
                    raise InputFormatError(f"Non-finite gene weight for {gene!r}.")
                numeric = float(weight)
                if numeric < 0:
                    raise InputFormatError(f"Negative gene weight for {gene!r}.")
                if numeric > 0:
                    values[gene] = max(values.get(gene, 0.0), numeric)
            output[str(cell_type)] = _normalize_gene_weights(values)
        return output

    if frame.shape[1] < 2:
        raise InputFormatError(
            "A wide gene-set table needs a gene column and at least one cell type."
        )
    first_column = frame.columns[0]
    # The wide format always reserves its first column for gene identifiers,
    # even when those identifiers happen to be numeric.
    frame = frame.set_index(first_column)
    for cell_type in frame.columns:
        numeric = pd.to_numeric(frame[cell_type], errors="coerce")
        values: dict[str, float] = {}
        for gene, weight in numeric.items():
            if pd.isna(weight):
                continue
            if not _is_float(weight):
                raise InputFormatError(f"Non-finite gene weight for {gene!r}.")
            numeric_weight = float(weight)
            if numeric_weight < 0:
                raise InputFormatError(f"Negative gene weight for {gene!r}.")
            if numeric_weight > 0:
                gene_text = str(gene)
                values[gene_text] = max(values.get(gene_text, 0.0), numeric_weight)
        if values:
            output[str(cell_type)] = _normalize_gene_weights(values)
    if not output:
        raise InputFormatError(f"No weighted cell-type gene sets were found in {path}.")
    return output


def load_regulatory_bed(
    path: str | Path,
    *,
    cell_type: str | None = None,
    weight: float = 1.0,
) -> pd.DataFrame:
    """Load BED-like regulatory intervals.

    Accepted layouts are:

    * BED3 with the ``cell_type`` and ``weight`` arguments supplied.
    * BED5-like ``chrom, start, end, cell_type, weight`` without a header.
    * A table with a header containing ``chrom``, ``start``, ``end`` and
      ``cell_type``; ``weight`` is optional.
    """

    path = Path(path)
    header_probe = pd.read_csv(path, sep=None, engine="python", nrows=0)
    normalized_headers = {
        str(column).strip().lower().lstrip("#"): str(column) for column in header_probe.columns
    }
    header_aliases = {
        "chrom": ("chrom", "chr", "chromosome", "contig"),
        "start": ("start", "chromstart", "begin"),
        "end": ("end", "chromend", "stop"),
        "cell_type": ("cell_type", "celltype", "cell", "celltype_name"),
        "weight": ("weight", "score", "gene_weight"),
    }
    resolved_headers = {
        canonical: next(
            (normalized_headers[alias] for alias in aliases if alias in normalized_headers),
            None,
        )
        for canonical, aliases in header_aliases.items()
    }
    has_header = all(resolved_headers[name] is not None for name in ("chrom", "start", "end"))
    if has_header:
        frame = read_table(path)
        rename = {
            resolved_headers["chrom"]: "chrom",
            resolved_headers["start"]: "start",
            resolved_headers["end"]: "end",
        }
        if resolved_headers["cell_type"] is not None:
            rename[resolved_headers["cell_type"]] = "cell_type"
        if resolved_headers["weight"] is not None:
            rename[resolved_headers["weight"]] = "weight"
        frame = frame.rename(columns=rename)
    else:
        frame = pd.read_csv(
            path,
            sep=r"\s+",
            comment="#",
            header=None,
        )
        if frame.shape[1] < 3:
            raise InputFormatError("Regulatory BED must contain at least three columns.")
        if frame.shape[1] == 3:
            frame.columns = ["chrom", "start", "end"]
        elif frame.shape[1] == 4:
            if cell_type is not None:
                frame.columns = ["chrom", "start", "end", "weight"]
            else:
                frame.columns = ["chrom", "start", "end", "cell_type"]
        else:
            frame = frame.iloc[:, :5]
            frame.columns = ["chrom", "start", "end", "cell_type", "weight"]

    required = {"chrom", "start", "end"}
    missing = required.difference(frame.columns)
    if missing:
        raise InputFormatError(f"Regulatory BED is missing columns: {sorted(missing)}")
    if cell_type is not None:
        frame["cell_type"] = cell_type
    if "cell_type" not in frame.columns:
        raise InputFormatError("Regulatory BED needs a cell_type column or the cell_type argument.")
    if "weight" not in frame.columns:
        frame["weight"] = weight
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
    frame["start"] = frame["start"].astype(np.int64)
    frame["end"] = frame["end"].astype(np.int64)
    frame["cell_type"] = frame["cell_type"].astype(str).str.strip()
    frame = frame.loc[frame["cell_type"] != ""]
    return frame.loc[:, ["chrom", "start", "end", "cell_type", "weight"]].reset_index(drop=True)


def load_expression_matrix(
    path: str | Path,
    *,
    gene_column: str | None = None,
) -> pd.DataFrame:
    """Load a genes-by-cell-types expression table."""

    frame = read_table(path)
    if gene_column is None:
        gene_column = _find_column(frame.columns, ("gene", "gene_name", "gene_symbol", "symbol"))
    if gene_column is not None:
        if gene_column not in frame.columns:
            raise InputFormatError(f"Expression gene column {gene_column!r} is absent.")
        frame = frame.set_index(gene_column)
    else:
        # Expression matrices are genes-by-cell-types; the first column stays
        # the gene identifier even when all identifiers are numeric.
        frame = frame.set_index(frame.columns[0])
    if frame.empty:
        raise InputFormatError("Expression matrix is empty.")
    try:
        numeric = frame.apply(pd.to_numeric, errors="raise")
    except (TypeError, ValueError) as error:
        raise InputFormatError("Expression matrix contains non-numeric values.") from error
    return numeric


def _find_column(
    columns: Sequence[object] | pd.Index,
    aliases: Sequence[str],
) -> str | None:
    normalized = {str(column).strip().lower().lstrip("#"): str(column) for column in columns}
    for alias in aliases:
        if alias in normalized:
            return normalized[alias]
    return None


def resolve_column(
    columns: Sequence[object] | pd.Index,
    aliases: Sequence[str],
    *,
    label: str,
) -> str:
    """Resolve a column name using case-insensitive aliases."""

    resolved = _find_column(columns, aliases)
    if resolved is None:
        raise InputFormatError(
            f"Could not identify the {label} column. Tried aliases {list(aliases)}; "
            f"available columns are {list(map(str, columns))}."
        )
    return resolved


def write_table(frame: pd.DataFrame, path: str | Path, *, index: bool = False) -> None:
    """Write a dataframe to CSV/TSV (optionally gzip) or parquet."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    suffixes = [suffix.lower() for suffix in path.suffixes]
    if ".parquet" in suffixes:
        frame.to_parquet(path, index=index)
        return
    active = suffixes[-2] if suffixes[-1] == ".gz" and len(suffixes) > 1 else suffixes[-1]
    separator = "\t" if active in {".tsv", ".tab"} else ","
    frame.to_csv(path, sep=separator, index=index)


def write_json(payload: Mapping[str, object], path: str | Path) -> None:
    """Write indented UTF-8 JSON."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(
            json_safe(payload),
            handle,
            ensure_ascii=False,
            indent=2,
            allow_nan=False,
        )


def json_safe(value: object) -> object:
    """Convert common scientific Python values to strict JSON values."""

    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return [json_safe(item) for item in value.tolist()]
    if pd.isna(value):
        return None
    return str(value)


def write_gene_sets(
    gene_sets: Mapping[str, Mapping[str, float]],
    path: str | Path,
) -> None:
    """Write weighted gene sets as a GMT-compatible file."""

    prepared: dict[str, dict[str, float]] = {}
    for cell_type, genes in gene_sets.items():
        name = str(cell_type).strip()
        if not name:
            raise InputFormatError("Gene-set names cannot be empty.")
        prepared[name] = _normalize_gene_weights(genes)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text_for_write(path) as handle:
        for cell_type, genes in prepared.items():
            tokens = [
                gene if weight == 1.0 else f"{gene}:{weight:g}" for gene, weight in genes.items()
            ]
            handle.write(f"{cell_type}\tcelltogwas\t" + "\t".join(tokens) + "\n")


def open_text_for_write(path: str | Path) -> io.TextIOBase:
    """Open plain or gzip text for writing based on suffix."""

    path = Path(path)
    if path.suffix.lower() == ".gz":
        return gzip.open(path, mode="wt", encoding="utf-8")
    return path.open(mode="wt", encoding="utf-8")
