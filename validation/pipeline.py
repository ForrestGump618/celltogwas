"""Real-data validation and MAGMA/H-MAGMA concordance pipeline.

The pipeline is intentionally separate from the installable ``celltogwas``
package: public data, MAGMA binaries, and multi-gigabyte reference panels are
not distributed with the library.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
import subprocess
import time
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from celltogwas import GenomeLinker, LinkConfig
from celltogwas._version import __version__
from celltogwas.enrichment import benjamini_hochberg, run_enrichment_by_chromosome
from celltogwas.io import (
    load_gene_annotation,
    load_gene_sets,
    write_gene_sets,
    write_json,
    write_table,
)
from celltogwas.specificity import compute_specificity, specificity_to_gene_sets

AUTOSOMES = {str(chromosome) for chromosome in range(1, 23)}
PSYCHENCODE_CELL_TYPES = (
    "Astro",
    "Chandelier__Pvalb",
    "Endo__VLMC",
    "Immune",
    "L2.3.IT",
    "L4.IT",
    "L5.6.NP",
    "L5.ET",
    "L5.IT",
    "L6.CT",
    "L6.IT.Car3",
    "L6.IT",
    "L6b",
    "Lamp5.Lhx6",
    "Lamp5",
    "Micro.PVM",
    "OPC",
    "Oligo",
    "PC",
    "Pax6",
    "SMC",
    "Sncg",
    "Sst__Sst.Chodl",
    "Vip",
)
KNOWN_AD_GENES = (
    "ABCA7",
    "ACE",
    "ADAM10",
    "APOE",
    "BIN1",
    "CD33",
    "CLU",
    "CR1",
    "INPP5D",
    "MS4A4A",
    "MS4A6A",
    "PICALM",
    "PLCG2",
    "SORL1",
    "SPI1",
    "TREM2",
)
RSID_PATTERN = re.compile(r"^rs\d+$")


@dataclass(frozen=True)
class ValidationConfig:
    """Analysis parameters fixed before inspecting the validation results."""

    data_root: Path
    output_dir: Path
    top_n: int = 1_000
    min_expression: float = 0.1
    min_detection_fraction: float = 0.1
    link_window: int = 10_000
    primary_p_upper: float = 5e-8
    sensitivity_p_upper: float = 1e-5
    n_permutations: int = 10_000
    random_state: int = 13
    bootstrap_replicates: int = 2_000
    magma_executable: Path | None = None
    magma_reference_prefix: Path | None = None
    gwas_path: Path | None = None
    psychencode_dir: Path | None = None
    gencode_path: Path | None = None
    hmagma_annotation_path: Path | None = None
    standard_annotation_path: Path | None = None

    def __post_init__(self) -> None:
        if self.top_n < 1:
            raise ValueError("top_n must be positive.")
        if not 0 <= self.min_expression < np.inf:
            raise ValueError("min_expression must be finite and non-negative.")
        if not 0 < self.min_detection_fraction <= 1:
            raise ValueError("min_detection_fraction must lie in (0, 1].")
        if self.link_window < 0:
            raise ValueError("link_window must be non-negative.")
        if not 0 < self.primary_p_upper <= 1:
            raise ValueError("primary_p_upper must lie in (0, 1].")
        if not self.primary_p_upper <= self.sensitivity_p_upper <= 1:
            raise ValueError("sensitivity_p_upper must be at least primary_p_upper.")
        if self.n_permutations < 1:
            raise ValueError("n_permutations must be positive.")
        if self.bootstrap_replicates < 1:
            raise ValueError("bootstrap_replicates must be positive.")

    def resolved_paths(self) -> ValidationConfig:
        """Return a copy with conventional data-root defaults resolved."""

        data_root = self.data_root.expanduser().resolve()
        defaults = {
            "magma_executable": data_root / "magma" / "magma_v1.10_win" / "magma.exe",
            "magma_reference_prefix": data_root / "magma" / "g1000_eur" / "g1000_eur",
            "gwas_path": data_root / "gwas" / "GCST90027158_buildGRCh38.tsv.gz",
            "psychencode_dir": data_root / "psychencode",
            "gencode_path": data_root / "gencode" / "gencode.v19.annotation.gtf.gz",
            "hmagma_annotation_path": data_root / "hmagma" / "Adult_brain.genes.annot",
            "standard_annotation_path": data_root / "hmagma" / "MAGMAdefault.genes.annot.gz",
        }
        values = asdict(self)
        for name, default in defaults.items():
            value = values[name]
            values[name] = default if value is None else Path(value).expanduser().resolve()
        values["data_root"] = data_root
        values["output_dir"] = self.output_dir.expanduser().resolve()
        return ValidationConfig(**values)


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp."""

    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: str | Path, *, chunk_size: int = 8 * 1024 * 1024) -> str:
    """Compute a SHA-256 digest with bounded memory."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _require_paths(paths: Iterable[Path | None]) -> None:
    missing = [str(path) for path in paths if path is None or not Path(path).exists()]
    if missing:
        raise FileNotFoundError("Required input files are missing:\n" + "\n".join(missing))


def _ensure_output_tree(output_dir: Path) -> dict[str, Path]:
    directories = {
        name: output_dir / name
        for name in (
            "expression",
            "gene_sets",
            "harmonized",
            "celltogwas",
            "magma",
            "comparison",
            "logs",
        )
    }
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=True)
    return directories


def read_psychencode_cell_type(
    path: str | Path,
) -> tuple[pd.Series, pd.Series, pd.Series, list[str]]:
    """Aggregate one brainSCOPE pseudobulk BED to donor means and detection.

    Rows are genes. The first six columns are genomic metadata and every
    remaining column is one donor. A gene symbol that appears more than once is
    represented by its strongest mean profile, avoiding double counting of
    duplicated symbols while retaining multi-locus genes for annotation.
    """

    path = Path(path)
    frame = pd.read_csv(path, sep="\t", compression="infer")
    required_prefix = ["#chr", "start", "end", "gene", "length", "strand"]
    if list(frame.columns[:6]) != required_prefix:
        raise ValueError(f"{path} does not use the expected six-column brainSCOPE header.")
    donor_ids = [str(column) for column in frame.columns[6:]]
    if not donor_ids or len(set(donor_ids)) != len(donor_ids):
        raise ValueError(f"{path} has empty or duplicate donor columns.")
    genes = frame["gene"].astype(str).str.strip()
    if (genes == "").any():
        raise ValueError(f"{path} contains an empty gene symbol.")
    values = frame.iloc[:, 6:].to_numpy(dtype=np.float64, copy=True)
    if not np.isfinite(values).all():
        raise ValueError(f"{path} contains non-finite expression values.")
    if np.any(values < 0):
        raise ValueError(f"{path} contains negative expression values.")

    means = pd.Series(values.mean(axis=1), index=genes, name="mean_expression")
    detection = pd.Series((values > 0).mean(axis=1), index=genes, name="detection_fraction")
    donors = pd.Series(np.full(len(genes), len(donor_ids), dtype=int), index=genes, name="n_donors")
    aggregated = (
        pd.concat([means, detection, donors], axis=1)
        .groupby(level=0, sort=True)
        .agg(
            {
                "mean_expression": "max",
                "detection_fraction": "max",
                "n_donors": "max",
            }
        )
    )
    return (
        aggregated["mean_expression"].astype(float),
        aggregated["detection_fraction"].astype(float),
        aggregated["n_donors"].astype(int),
        donor_ids,
    )


def build_expression_and_markers(
    config: ValidationConfig,
    directories: Mapping[str, Path],
) -> dict[str, object]:
    """Build cross-cell-type marker weights and matching binary gene sets."""

    psychencode_dir = Path(config.psychencode_dir)
    gencode_path = Path(config.gencode_path)
    _require_paths([psychencode_dir, gencode_path])

    mean_profiles: dict[str, pd.Series] = {}
    detection_profiles: dict[str, pd.Series] = {}
    donor_profiles: dict[str, pd.Series] = {}
    cell_type_rows: list[dict[str, object]] = []
    all_donors: set[str] = set()
    common_donors: set[str] | None = None

    for cell_type in PSYCHENCODE_CELL_TYPES:
        path = psychencode_dir / f"{cell_type}.expr.bed.gz"
        if not path.exists():
            raise FileNotFoundError(f"Missing PsychENCODE matrix for {cell_type}: {path}")
        mean_values, detection_values, donor_counts, donor_ids = read_psychencode_cell_type(path)
        mean_profiles[cell_type] = mean_values
        detection_profiles[cell_type] = detection_values
        donor_profiles[cell_type] = donor_counts
        all_donors.update(donor_ids)
        common_donors = set(donor_ids) if common_donors is None else common_donors & set(donor_ids)
        cell_type_rows.append(
            {
                "cell_type": cell_type,
                "n_donors": len(donor_ids),
                "n_genes": int(len(mean_values)),
                "n_detected_genes": int((mean_values > 0).sum()),
                "mean_expression_median": float(mean_values.median()),
                "detection_fraction_median": float(detection_values.median()),
                "source_file": path.name,
                "sha256": sha256_file(path),
            }
        )

    expression = pd.concat(mean_profiles, axis=1).fillna(0.0)
    detection = pd.concat(detection_profiles, axis=1).fillna(0.0)
    donor_counts = pd.concat(donor_profiles, axis=1).fillna(0).astype(int)
    write_table(
        expression.rename_axis("gene").reset_index(),
        directories["expression"] / "cell_type_mean_expression.tsv.gz",
    )
    write_table(
        detection.rename_axis("gene").reset_index(),
        directories["expression"] / "cell_type_detection_fraction.tsv.gz",
    )
    write_table(
        donor_counts.rename_axis("gene").reset_index(),
        directories["expression"] / "cell_type_donors.tsv.gz",
    )

    annotation = load_gene_annotation(gencode_path, feature_types=("gene",))
    symbol_to_ids: dict[str, set[str]] = {}
    for record in annotation.genes.itertuples(index=False):
        symbol = str(record.gene_name).strip().upper()
        gene_id = str(record.gene_id).strip().split(".", 1)[0]
        if symbol and gene_id:
            symbol_to_ids.setdefault(symbol, set()).add(gene_id)
    gencode_symbols = set(symbol_to_ids)

    common_genes = set.intersection(*(set(profile.index) for profile in mean_profiles.values()))
    complete_expression = expression.loc[expression.index.isin(common_genes)]
    all_weights = compute_specificity(
        complete_expression,
        method="tau",
        min_expression=config.min_expression,
        power=1.0,
        top_n=None,
        min_weight=0.0,
    )
    all_weights["max_expression"] = all_weights["max_expression"].astype(float)
    detection_long = (
        detection.rename_axis("gene")
        .reset_index()
        .melt(
            id_vars="gene",
            var_name="cell_type",
            value_name="detection_fraction",
        )
    )
    donor_long = (
        donor_counts.rename_axis("gene")
        .reset_index()
        .melt(
            id_vars="gene",
            var_name="cell_type",
            value_name="n_donors",
        )
    )
    all_weights = all_weights.merge(
        detection_long, on=["cell_type", "gene"], validate="many_to_one"
    )
    all_weights = all_weights.merge(donor_long, on=["cell_type", "gene"], validate="many_to_one")
    eligible = all_weights.loc[
        (all_weights["weight"] > 0)
        & (all_weights["max_expression"] >= config.min_expression)
        & (all_weights["detection_fraction"] >= config.min_detection_fraction)
        & all_weights["gene"].astype(str).str.strip().str.upper().isin(gencode_symbols)
    ].copy()
    markers = (
        eligible.sort_values(["cell_type", "weight", "gene"], ascending=[True, False, True])
        .groupby("cell_type", sort=False, as_index=False)
        .head(config.top_n)
        .reset_index(drop=True)
    )
    missing_cell_types = set(PSYCHENCODE_CELL_TYPES) - set(markers["cell_type"])
    if missing_cell_types:
        raise ValueError(f"No eligible markers for cell types: {sorted(missing_cell_types)}")
    write_table(markers, directories["expression"] / "marker_weights.tsv.gz")

    marker_summary = (
        markers.groupby("cell_type", sort=False)
        .agg(
            n_markers=("gene", "nunique"),
            min_weight=("weight", "min"),
            median_weight=("weight", "median"),
            max_weight=("weight", "max"),
            median_max_expression_across_cell_types=("max_expression", "median"),
            median_detection_fraction=("detection_fraction", "median"),
        )
        .reset_index()
    )
    write_table(marker_summary, directories["expression"] / "marker_summary.tsv")

    weighted_sets = specificity_to_gene_sets(markers)
    membership_sets = {
        cell_type: dict.fromkeys(group["gene"].astype(str), 1.0)
        for cell_type, group in markers.groupby("cell_type", sort=False)
    }
    write_gene_sets(weighted_sets, directories["gene_sets"] / "markers_weighted.gmt")
    write_gene_sets(membership_sets, directories["gene_sets"] / "markers_membership.gmt")

    ensembl_sets: dict[str, dict[str, float]] = {}
    unmapped: dict[str, list[str]] = {}
    for cell_type, genes in membership_sets.items():
        gene_ids: set[str] = set()
        missing: list[str] = []
        for gene in genes:
            mapped = symbol_to_ids.get(str(gene).strip().upper(), set())
            if mapped:
                gene_ids.update(mapped)
            else:
                missing.append(str(gene))
        if not gene_ids:
            raise ValueError(f"No GENCODE v19 Ensembl IDs map to markers for {cell_type}.")
        ensembl_sets[cell_type] = dict.fromkeys(sorted(gene_ids), 1.0)
        unmapped[cell_type] = missing
    write_gene_sets(ensembl_sets, directories["gene_sets"] / "markers_membership_ensembl.gmt")

    payload = {
        "dataset": "PsychENCODE brainSCOPE adult DLPFC pseudobulk",
        "source": "https://brainscope.gersteinlab.org/output-pseudobulk-expr.html",
        "cell_type_count": len(PSYCHENCODE_CELL_TYPES),
        "n_genes_union": int(expression.shape[0]),
        "n_genes_present_in_all_cell_types": len(common_genes),
        "n_donor_ids_union": len(all_donors),
        "n_donors_present_in_all_cell_types": len(common_donors or set()),
        "cell_types": cell_type_rows,
        "marker_parameters": {
            "method": "tau",
            "power": 1.0,
            "min_expression": config.min_expression,
            "min_detection_fraction": config.min_detection_fraction,
            "top_n": config.top_n,
            "gene_universe": (
                "Genes reported in all 24 matrices and represented by a GENCODE v19 gene symbol; this common "
                "universe is fixed before GWAS analysis."
            ),
        },
        "n_genes_in_gencode_v19_symbol_universe": int(
            expression.index.astype(str).str.strip().str.upper().isin(gencode_symbols).sum()
        ),
        "n_markers_requested": int(len(PSYCHENCODE_CELL_TYPES) * config.top_n),
        "n_markers_selected": int(len(markers)),
        "n_unique_markers_selected": int(markers["gene"].nunique()),
        "ensembl_set_sizes": {cell_type: len(genes) for cell_type, genes in ensembl_sets.items()},
        "marker_summary": marker_summary.to_dict(orient="records"),
        "unmapped_marker_symbols": {key: sorted(value) for key, value in unmapped.items()},
        "interpretation_limits": [
            "Donor columns differ by cell type; means use all donors in which a cell type was observed.",
            "Only one donor is represented in all 24 cell types, so a complete-donor-only profile is not estimable.",
            "Marker weights are derived from expression specificity, not from GWAS association.",
        ],
    }
    write_json(payload, directories["expression"] / "expression_manifest.json")
    return payload


def harmonize_gwas(
    gwas_path: str | Path,
    bim_path: str | Path,
    output_dir: str | Path,
) -> dict[str, object]:
    """Restrict GWAS to unique rsIDs in the GRCh37 1000G EUR reference panel."""

    gwas_path = Path(gwas_path)
    bim_path = Path(bim_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    _require_paths([gwas_path, bim_path])

    # Numeric rsIDs keep a 22-million-variant reference compact in memory.
    reference_parts = []
    n_bim_rows = 0
    for chunk in pd.read_csv(
        bim_path,
        sep=r"\s+",
        header=None,
        usecols=[0, 1, 3],
        names=["chrom", "rsid", "pos"],
        dtype={"rsid": "str"},
        chunksize=1_000_000,
    ):
        n_bim_rows += len(chunk)
        rs = chunk["rsid"].str.strip().str.lower()
        chrom = pd.to_numeric(chunk["chrom"], errors="coerce")
        valid = rs.str.fullmatch(RSID_PATTERN) & chrom.between(1, 22)
        part = np.empty(int(valid.sum()), dtype=[("rs", "u8"), ("chrom", "u1"), ("pos", "u4")])
        part["rs"] = rs.loc[valid].str[2:].to_numpy(dtype=np.uint64)
        part["chrom"] = chrom.loc[valid].to_numpy(dtype=np.uint8)
        part["pos"] = chunk.loc[valid, "pos"].to_numpy(dtype=np.uint32)
        reference_parts.append(part)
    reference = np.concatenate(reference_parts)
    del reference_parts
    reference.sort(order="rs")
    duplicate = np.zeros(len(reference), dtype=bool)
    equal = reference["rs"][1:] == reference["rs"][:-1]
    duplicate[1:] |= equal
    duplicate[:-1] |= equal
    duplicate_bim_count = int(duplicate.sum())
    reference = reference[~duplicate]
    n_reference = len(reference)
    reference_ids = np.ascontiguousarray(reference["rs"])
    print(f"Reference loaded: {n_reference:,} unique autosomal rsIDs", flush=True)

    parts = []
    total_gwas_rows = 0
    n_valid_gwas_rows = 0
    all_valid_ids = []
    columns = ["variant_id", "p_value", "n_cases", "n_controls"]
    for chunk in pd.read_csv(gwas_path, sep="\t", usecols=columns, chunksize=500_000):
        total_gwas_rows += len(chunk)
        rs = chunk["variant_id"].fillna("").astype(str).str.strip().str.lower()
        p_values = pd.to_numeric(chunk["p_value"], errors="coerce")
        cases = pd.to_numeric(chunk["n_cases"], errors="coerce")
        controls = pd.to_numeric(chunk["n_controls"], errors="coerce")
        n = cases + controls
        valid = (
            rs.str.fullmatch(RSID_PATTERN)
            & p_values.between(0, 1)
            & np.isfinite(n)
            & (cases >= 0)
            & (controls >= 0)
            & (n > 0)
        )
        ids = rs.loc[valid].str[2:].to_numpy(dtype=np.uint64)
        n_valid_gwas_rows += len(ids)
        all_valid_ids.append(ids)
        locations = np.searchsorted(reference_ids, ids)
        locations = np.minimum(locations, len(reference_ids) - 1)
        keep = reference_ids[locations] == ids
        part = np.empty(
            int(keep.sum()),
            dtype=[("rs", "u8"), ("chrom", "u1"), ("pos", "u4"), ("p", "f8"), ("n", "f8")],
        )
        part["rs"] = ids[keep]
        part["chrom"] = reference["chrom"][locations[keep]]
        part["pos"] = reference["pos"][locations[keep]]
        part["p"] = p_values.loc[valid].to_numpy(dtype=float)[keep]
        part["n"] = n.loc[valid].to_numpy(dtype=float)[keep]
        parts.append(part)
    all_ids = np.concatenate(all_valid_ids)
    unique_ids, counts = np.unique(all_ids, return_counts=True)
    repeated_ids = unique_ids[counts > 1]
    n_duplicate_gwas_rows = int(counts[counts > 1].sum())
    n_unique_gwas_rsids = int((counts == 1).sum())
    del all_valid_ids, all_ids, unique_ids, counts, reference, reference_ids
    matched = np.concatenate(parts)
    del parts
    matched = matched[~np.isin(matched["rs"], repeated_ids)]
    n_below_magma_min_n = int((matched["n"] < 50).sum())
    matched = matched[matched["n"] >= 50]
    matched.sort(order=["chrom", "pos", "rs"])
    if not len(matched):
        raise ValueError("No unique autosomal GWAS rsIDs match the reference with N >= 50.")
    n_zero_p_rows = int((matched["p"] == 0).sum())
    chromosome_counts = {}
    chromosome_dir = output_dir / "chromosomes"
    chromosome_dir.mkdir(exist_ok=True)
    # MAGMA v1.10 expects plain text for the P-value file, too.
    with (
        gzip.open(
            output_dir / "gwas_harmonized.tsv.gz", "wt", encoding="utf-8", newline=""
        ) as combined,
        (output_dir / "gwas_for_magma.tsv").open("w", encoding="utf-8", newline="") as magma,
    ):
        first = True
        for chromosome in range(1, 23):
            subset = matched[matched["chrom"] == chromosome]
            if not len(subset):
                continue
            chromosome_counts[str(chromosome)] = len(subset)
            frame = pd.DataFrame(
                {
                    "SNP": [f"rs{rs}" for rs in subset["rs"]],
                    "CHR": subset["chrom"],
                    "POS": subset["pos"],
                    "P": np.maximum(subset["p"], 1e-300),
                    "N": np.rint(subset["n"]).astype(np.int64),
                }
            )
            frame.to_csv(combined, sep="\t", index=False, header=first)
            frame.rename(columns={"SNP": "rsid"})[["rsid", "P", "N"]].to_csv(
                magma, sep="\t", index=False, header=first
            )
            frame.to_csv(chromosome_dir / f"chr{chromosome}.tsv.gz", sep="\t", index=False)
            first = False
    mhc = (matched["chrom"] == 6) & (matched["pos"] >= 25_000_000) & (matched["pos"] <= 35_000_000)
    payload = {
        "gwas": {
            "study": "Bellenguez Alzheimer's disease GWAS, accession GCST90027158",
            "accession": "GCST90027158",
            "source_build": "GRCh38",
            "analysis_build": "GRCh37 / hg19",
            "coordinate_strategy": "Exact rsID match to 1000G EUR BIM; reference positions; no liftover. P-only analysis, no effect-sign harmonization.",
            "total_rows": total_gwas_rows,
            "valid_rows_before_duplicate_filter": n_valid_gwas_rows,
            "duplicate_rsid_rows_removed": n_duplicate_gwas_rows,
            "unique_valid_rsids": n_unique_gwas_rsids,
            "matched_rsids": len(matched),
            "matched_rsids_below_magma_min_n_removed": n_below_magma_min_n,
            "unmatched_or_low_n_unique_valid_rsids": n_unique_gwas_rsids - len(matched),
            "match_fraction_of_unique_valid": len(matched) / n_unique_gwas_rsids,
            "zero_p_rows": n_zero_p_rows,
            "zero_p_values_clamped_to": 1e-300,
            "sample_size": "N = n_cases + n_controls; N >= 50, matching MAGMA's default",
            "per_snp_n_min": float(matched["n"].min()),
            "per_snp_n_median": float(np.median(matched["n"])),
            "per_snp_n_max": float(matched["n"].max()),
            "chromosome_counts": chromosome_counts,
            "mhc_snps_included_in_primary": int(mhc.sum()),
            "mhc_definition": "GRCh37 chr6:25000000-35000000",
        },
        "reference": {
            "name": "MAGMA 1000 Genomes EUR reference",
            "bim_rows": n_bim_rows,
            "autosomal_unique_rsids": n_reference,
            "duplicate_autosomal_rsid_rows_removed": duplicate_bim_count,
            "reference_individuals": 503,
        },
        "inputs_sha256": {"gwas": sha256_file(gwas_path), "bim": sha256_file(bim_path)},
        "outputs": {
            "harmonized": str((output_dir / "gwas_harmonized.tsv.gz").resolve()),
            "magma_pvalue_file": str((output_dir / "gwas_for_magma.tsv").resolve()),
        },
    }
    write_json(payload, output_dir / "gwas_harmonization.json")
    return payload


def _run_logged_command(
    command: Sequence[str],
    *,
    cwd: Path,
    log_path: Path,
) -> dict[str, object]:
    """Run one external command and persist a complete reproducibility log."""

    started = time.perf_counter()
    result = subprocess.run(
        list(command),
        cwd=cwd,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    elapsed = time.perf_counter() - started
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as handle:
        handle.write(f"command={json.dumps(list(command), ensure_ascii=True)}\n")
        handle.write(f"cwd={cwd}\n")
        handle.write(f"returncode={result.returncode}\n")
        handle.write(f"elapsed_seconds={elapsed:.6f}\n")
        handle.write("\n[stdout]\n")
        handle.write(result.stdout)
        handle.write("\n[stderr]\n")
        handle.write(result.stderr)
    if result.returncode != 0 or re.search(r"(?m)^ERROR\b", result.stdout):
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}; see {log_path}. "
            f"stderr: {result.stderr[-2000:]}"
        )
    return {
        "command": list(command),
        "returncode": result.returncode,
        "elapsed_seconds": elapsed,
        "log": str(log_path.resolve()),
    }


def _ensure_uncompressed_annotation(source: Path, directory: Path) -> Path:
    """MAGMA v1.10 does not read gzip-compressed gene annotation files."""

    if source.suffix.lower() != ".gz":
        return source
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / source.name[:-3]
    temporary = destination.with_name(destination.name + ".part")
    with gzip.open(source, "rb") as input_handle, temporary.open("wb") as output_handle:
        while chunk := input_handle.read(8 * 1024 * 1024):
            output_handle.write(chunk)
    temporary.replace(destination)
    return destination


def run_magma_analyses(
    config: ValidationConfig,
    directories: Mapping[str, Path],
) -> dict[str, object]:
    """Run standard MAGMA and H-MAGMA with the same executable and GWAS input."""

    executable = Path(config.magma_executable)
    reference_prefix = Path(config.magma_reference_prefix)
    standard_annotation = _ensure_uncompressed_annotation(
        Path(config.standard_annotation_path),
        directories["magma"] / "references",
    )
    hmagma_annotation = Path(config.hmagma_annotation_path)
    pvalue_file = directories["harmonized"] / "gwas_for_magma.tsv"
    marker_gmt = directories["gene_sets"] / "markers_membership_ensembl.gmt"
    _require_paths([executable, standard_annotation, hmagma_annotation, pvalue_file, marker_gmt])
    for suffix in (".bed", ".bim", ".fam"):
        reference_file = Path(f"{reference_prefix}{suffix}")
        if not reference_file.exists():
            raise FileNotFoundError(f"MAGMA reference component is missing: {reference_file}")

    version_log = directories["logs"] / "magma_version.log"
    version_result = _run_logged_command(
        [str(executable), "--version"],
        cwd=directories["magma"],
        log_path=version_log,
    )
    version_text = version_log.read_text(encoding="utf-8")
    if "v1.10" not in version_text:
        raise RuntimeError(f"Expected MAGMA v1.10; see {version_log}.")

    analyses: dict[str, dict[str, object]] = {}
    for label, annotation in (
        ("standard", standard_annotation),
        ("hmagma", hmagma_annotation),
    ):
        gene_prefix = directories["magma"] / f"{label}_gene"
        gene_log = directories["logs"] / f"{label}_gene.log"
        gene_command = [
            str(executable),
            "--bfile",
            str(reference_prefix),
            "--pval",
            str(pvalue_file),
            "use=rsid,P",
            "ncol=N",
            "--gene-annot",
            str(annotation),
            "--out",
            str(gene_prefix),
        ]
        gene_run = _run_logged_command(
            gene_command,
            cwd=directories["magma"],
            log_path=gene_log,
        )

        gene_set_prefix = directories["magma"] / f"{label}_geneset"
        gene_set_log = directories["logs"] / f"{label}_geneset.log"
        gene_set_command = [
            str(executable),
            "--gene-results",
            str(gene_prefix) + ".genes.raw",
            "--set-annot",
            str(marker_gmt),
            "--out",
            str(gene_set_prefix),
        ]
        gene_set_run = _run_logged_command(
            gene_set_command,
            cwd=directories["magma"],
            log_path=gene_set_log,
        )
        analyses[label] = {
            "annotation": str(annotation.resolve()),
            "annotation_sha256": sha256_file(annotation),
            "gene_run": gene_run,
            "gene_set_run": gene_set_run,
            "gene_results": str(
                (gene_prefix.parent / f"{gene_prefix.name}.genes.out.txt").resolve()
            ),
            "gene_set_results": str(
                (gene_set_prefix.parent / f"{gene_set_prefix.name}.gsa.out.txt").resolve()
            ),
        }

    payload = {
        "engine": "MAGMA",
        "version": "v1.10 (win/s)",
        "executable_sha256": sha256_file(executable),
        "reference_prefix": str(reference_prefix.resolve()),
        "pvalue_file": str(pvalue_file.resolve()),
        "pvalue_file_sha256": sha256_file(pvalue_file),
        "gene_set_file": str(marker_gmt.resolve()),
        "gene_set_file_sha256": sha256_file(marker_gmt),
        "common_engine_statement": (
            "H-MAGMA is not a separate statistical engine here: both analyses "
            "use the same MAGMA v1.10 executable, GWAS SNP/P/N file, and 1000G "
            "EUR reference. Only the SNP-to-gene annotation differs."
        ),
        "version_run": version_result,
        "analyses": analyses,
    }
    write_json(payload, directories["magma"] / "magma_run_manifest.json")
    return payload


def _run_enrichment_variant(
    linker: GenomeLinker,
    gwas: Path,
    *,
    name: str,
    p_upper: float,
    weight_floor: float,
    n_permutations: int,
    random_state: int,
    output_dir: Path,
    exclude_mhc: bool = False,
) -> dict[str, object]:
    chromosome_paths = sorted(gwas.glob("chr*.tsv.gz"))
    counts = {"input": 0, "selected": 0, "removed_mhc": 0}

    def chromosomes():
        for path in chromosome_paths:
            frame = pd.read_csv(path, sep="\t", usecols=["CHR", "POS", "P"])
            counts["input"] += len(frame)
            counts["selected"] += int((frame["P"] <= p_upper).sum())
            if exclude_mhc:
                mhc = (frame["CHR"] == 6) & frame["POS"].between(25_000_000, 35_000_000)
                counts["removed_mhc"] += int(mhc.sum())
                frame = frame.loc[~mhc].copy()
            print(
                f"{name}: chromosome {int(frame['CHR'].iloc[0])}, {len(frame):,} background SNPs",
                flush=True,
            )
            yield linker.annotate(frame, mode="weight", chunk_size=200_000)

    result = run_enrichment_by_chromosome(
        chromosomes(),
        cell_types=list(linker.cell_types),
        p_column="P",
        p_upper=p_upper,
        weight_floor=weight_floor,
        n_permutations=n_permutations,
        random_state=random_state,
    )
    table_path = output_dir / f"{name}.tsv"
    report_path = output_dir / f"{name}.json"
    write_table(result.table, table_path)
    payload = {
        "name": name,
        "parameters": {
            "p_upper": p_upper,
            "weight_floor": weight_floor,
            "n_permutations": n_permutations,
            "random_state": random_state,
            "min_shift_fraction": 0.01,
            "exclude_mhc": exclude_mhc,
        },
        "n_input_harmonized_snps": counts["input"],
        "n_snps_selected_by_p": counts["selected"],
        "n_mhc_snps_removed": counts["removed_mhc"],
        "n_snps_annotated": counts["input"] - counts["removed_mhc"],
        "diagnostics": result.diagnostics,
        "table": str(table_path.resolve()),
    }
    write_json(payload, report_path)
    return payload


def run_celltogwas_analyses(
    config: ValidationConfig,
    directories: Mapping[str, Path],
) -> dict[str, object]:
    """Run weighted and membership-aligned celltogwas enrichment analyses."""

    gwas_path = directories["harmonized"] / "gwas_harmonized.tsv.gz"
    weighted_gmt = directories["gene_sets"] / "markers_weighted.gmt"
    membership_gmt = directories["gene_sets"] / "markers_membership.gmt"
    _require_paths([gwas_path, weighted_gmt, membership_gmt, config.gencode_path])
    gwas = directories["harmonized"] / "chromosomes"
    annotation = load_gene_annotation(Path(config.gencode_path), feature_types=("gene",))
    link_config = LinkConfig(
        window=config.link_window,
        promoter_window=0,
        promoter_weight=0.0,
        gene_body_weight=1.0,
        kernel="uniform",
        max_weight=1.0,
    )
    weighted_linker = GenomeLinker.from_annotation(
        annotation,
        load_gene_sets(weighted_gmt),
        config=link_config,
    )
    membership_linker = GenomeLinker.from_annotation(
        annotation,
        load_gene_sets(membership_gmt),
        config=link_config,
    )
    primary_floor = -math.log10(config.primary_p_upper)
    sensitivity_floor = -math.log10(config.sensitivity_p_upper)
    output_dir = directories["celltogwas"]
    variants = {
        "sensitivity_all_snps": _run_enrichment_variant(
            weighted_linker,
            gwas,
            name="sensitivity_all_snps",
            p_upper=1.0,
            weight_floor=0.0,
            n_permutations=config.n_permutations,
            random_state=config.random_state,
            output_dir=output_dir,
        ),
        "primary_weighted": _run_enrichment_variant(
            weighted_linker,
            gwas,
            name="primary_weighted",
            p_upper=config.primary_p_upper,
            weight_floor=primary_floor,
            n_permutations=config.n_permutations,
            random_state=config.random_state,
            output_dir=output_dir,
        ),
        "membership_aligned": _run_enrichment_variant(
            membership_linker,
            gwas,
            name="membership_aligned",
            p_upper=config.primary_p_upper,
            weight_floor=primary_floor,
            n_permutations=config.n_permutations,
            random_state=config.random_state,
            output_dir=output_dir,
        ),
        "sensitivity_p1e5": _run_enrichment_variant(
            weighted_linker,
            gwas,
            name="sensitivity_p1e5",
            p_upper=config.sensitivity_p_upper,
            weight_floor=sensitivity_floor,
            n_permutations=config.n_permutations,
            random_state=config.random_state,
            output_dir=output_dir,
        ),
        "sensitivity_no_mhc": _run_enrichment_variant(
            weighted_linker,
            gwas,
            name="sensitivity_no_mhc",
            p_upper=config.primary_p_upper,
            weight_floor=primary_floor,
            n_permutations=config.n_permutations,
            random_state=config.random_state,
            output_dir=output_dir,
            exclude_mhc=True,
        ),
    }
    payload = {
        "package_version": __version__,
        "genomic_link_definition": {
            "annotation": "GENCODE v19 gene bodies plus bilateral 10 kb",
            "window": config.link_window,
            "promoter_window": 0,
            "promoter_weight": 0.0,
            "gene_body_weight": 1.0,
            "kernel": "uniform",
            "overlap_rule": "maximum active weight per cell type",
        },
        "weighted_linker_diagnostics": weighted_linker.diagnostics.to_dict(orient="records"),
        "membership_linker_diagnostics": membership_linker.diagnostics.to_dict(orient="records"),
        "variants": variants,
        "interpretation_limit": (
            "celltogwas uses circular annotation shifts and does not model LD. "
            "Agreement with MAGMA is rank concordance, not numerical equivalence "
            "or causal validation."
        ),
    }
    write_json(payload, output_dir / "celltogwas_run_manifest.json")
    return payload


def _column_name(frame: pd.DataFrame, aliases: Sequence[str]) -> str:
    normalized = {str(column).strip().upper(): str(column) for column in frame.columns}
    for alias in aliases:
        if alias.upper() in normalized:
            return normalized[alias.upper()]
    raise KeyError(f"None of {list(aliases)} are present in columns {list(frame.columns)}")


def read_magma_gene_results(path: str | Path) -> pd.DataFrame:
    """Read a MAGMA ``.genes.out`` table into canonical columns."""

    frame = pd.read_csv(path, sep=r"\s+", comment="#")
    gene = _column_name(frame, ("GENE",))
    zstat = _column_name(frame, ("ZSTAT", "Z"))
    pvalue = _column_name(frame, ("P",))
    output = pd.DataFrame(
        {
            "GENE": frame[gene].astype(str),
            "ZSTAT": pd.to_numeric(frame[zstat], errors="coerce"),
            "P": pd.to_numeric(frame[pvalue], errors="coerce"),
        }
    )
    try:
        output["NSNPS"] = pd.to_numeric(
            frame[_column_name(frame, ("NSNPS",))],
            errors="coerce",
        )
    except KeyError:
        output["NSNPS"] = np.nan
    return output.loc[np.isfinite(output["ZSTAT"]) & np.isfinite(output["P"])].reset_index(
        drop=True
    )


def read_magma_gene_set_results(path: str | Path) -> pd.DataFrame:
    """Read a MAGMA ``.gsa.out`` table into canonical columns."""

    frame = pd.read_csv(path, sep=r"\s+", comment="#")
    set_column = _column_name(
        frame,
        ("FULL_NAME", "VARIABLE", "SET", "GENESET", "GENE_SET"),
    )
    z_column = None
    for candidate in ("ZSTAT", "Z"):
        try:
            z_column = _column_name(frame, (candidate,))
            break
        except KeyError:
            continue
    output = pd.DataFrame({"cell_type": frame[set_column].astype(str)})
    if z_column is not None:
        output["magma_z"] = pd.to_numeric(frame[z_column], errors="coerce")
    else:
        beta = pd.to_numeric(frame[_column_name(frame, ("BETA",))], errors="coerce")
        standard_error = pd.to_numeric(frame[_column_name(frame, ("SE",))], errors="coerce")
        output["magma_z"] = beta / standard_error.replace(0.0, np.nan)
    try:
        output["magma_p"] = pd.to_numeric(frame[_column_name(frame, ("P",))], errors="coerce")
    except KeyError:
        output["magma_p"] = np.nan
    try:
        output["n_genes"] = pd.to_numeric(
            frame[_column_name(frame, ("NGENES", "N_GENES"))],
            errors="coerce",
        ).astype("Int64")
    except KeyError:
        output["n_genes"] = pd.Series(pd.NA, index=output.index, dtype="Int64")
    return output.loc[np.isfinite(output["magma_z"])].reset_index(drop=True)


def spearman_rank(
    first: Sequence[float] | np.ndarray, second: Sequence[float] | np.ndarray
) -> float:
    """Compute Spearman correlation from average ranks without SciPy."""

    left = pd.Series(np.asarray(first, dtype=float)).rank(method="average")
    right = pd.Series(np.asarray(second, dtype=float)).rank(method="average")
    if left.isna().any() or right.isna().any():
        raise ValueError("Spearman inputs must be finite.")
    left_values = left.to_numpy(dtype=float)
    right_values = right.to_numpy(dtype=float)
    if len(left_values) < 2:
        return float("nan")
    left_centered = left_values - left_values.mean()
    right_centered = right_values - right_values.mean()
    denominator = math.sqrt(
        float(left_centered @ left_centered) * float(right_centered @ right_centered)
    )
    if denominator <= np.finfo(float).eps:
        return float("nan")
    return float(np.clip(left_centered @ right_centered / denominator, -1.0, 1.0))


def bootstrap_spearman_ci(
    first: Sequence[float] | np.ndarray,
    second: Sequence[float] | np.ndarray,
    *,
    replicates: int,
    random_state: int,
) -> tuple[float, float]:
    """Return a percentile bootstrap CI for a rank correlation."""

    left = pd.Series(np.asarray(first, dtype=float)).rank(method="average").to_numpy()
    right = pd.Series(np.asarray(second, dtype=float)).rank(method="average").to_numpy()
    if len(left) != len(right) or len(left) < 3:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(random_state)
    values = np.empty(replicates, dtype=float)
    for index in range(replicates):
        sample = rng.integers(0, len(left), size=len(left))
        values[index] = spearman_rank(left[sample], right[sample])
    finite = values[np.isfinite(values)]
    if not len(finite):
        return (float("nan"), float("nan"))
    lower, upper = np.quantile(finite, [0.025, 0.975])
    return (float(lower), float(upper))


def top_overlap(
    first: Sequence[tuple[str, float]],
    second: Sequence[tuple[str, float]],
    *,
    k: int,
) -> tuple[int, float, list[str], list[str]]:
    """Compare the top-k ranked identifiers from two score tables."""

    first_top = [
        identifier for identifier, _ in sorted(first, key=lambda item: (-item[1], item[0]))[:k]
    ]
    second_top = [
        identifier for identifier, _ in sorted(second, key=lambda item: (-item[1], item[0]))[:k]
    ]
    first_set = set(first_top)
    second_set = set(second_top)
    intersection = first_set & second_set
    union = first_set | second_set
    jaccard = len(intersection) / len(union) if union else 1.0
    return len(intersection), float(jaccard), first_top, second_top


def _symbol_to_versionless_ids(annotation_path: Path) -> dict[str, list[str]]:
    annotation = load_gene_annotation(annotation_path, feature_types=("gene",))
    mapping: dict[str, set[str]] = {}
    for record in annotation.genes.itertuples(index=False):
        symbol = str(record.gene_name).strip().upper()
        gene_id = str(record.gene_id).strip().split(".", 1)[0]
        if symbol and gene_id:
            mapping.setdefault(symbol, set()).add(gene_id)
    return {symbol: sorted(ids) for symbol, ids in mapping.items()}


def compare_magma_gene_results(
    standard_path: str | Path,
    hmagma_path: str | Path,
    *,
    output_dir: Path,
    gencode_path: Path,
    bootstrap_replicates: int,
    random_state: int,
) -> dict[str, object]:
    """Compare standard MAGMA and H-MAGMA gene-level statistics."""

    standard = read_magma_gene_results(standard_path).rename(
        columns={"ZSTAT": "standard_z", "P": "standard_p", "NSNPS": "standard_nsnps"}
    )
    hmagma = read_magma_gene_results(hmagma_path).rename(
        columns={"ZSTAT": "hmagma_z", "P": "hmagma_p", "NSNPS": "hmagma_nsnps"}
    )
    merged = standard.merge(hmagma, on="GENE", how="inner", validate="one_to_one")
    merged = merged.loc[
        np.isfinite(merged["standard_z"])
        & np.isfinite(merged["hmagma_z"])
        & (merged["standard_nsnps"].fillna(0) > 0)
        & (merged["hmagma_nsnps"].fillna(0) > 0)
    ].copy()
    rho_z = spearman_rank(merged["standard_z"], merged["hmagma_z"])
    rho_p = spearman_rank(merged["standard_p"], merged["hmagma_p"])
    ci_z = bootstrap_spearman_ci(
        merged["standard_z"],
        merged["hmagma_z"],
        replicates=bootstrap_replicates,
        random_state=random_state,
    )
    merged["standard_rank"] = merged["standard_z"].rank(ascending=False, method="min").astype(int)
    merged["hmagma_rank"] = merged["hmagma_z"].rank(ascending=False, method="min").astype(int)
    top_metrics: dict[str, object] = {}
    native_top_metrics: dict[str, object] = {}
    for k in (100, 500, 1_000):
        overlap, jaccard, _, _ = top_overlap(
            list(zip(merged["GENE"], merged["standard_z"], strict=False)),
            list(zip(merged["GENE"], merged["hmagma_z"], strict=False)),
            k=k,
        )
        top_metrics[str(k)] = {
            "overlap": overlap,
            "jaccard": jaccard,
            "fraction_of_available_top_k": overlap / min(k, len(merged)),
        }
        overlap, jaccard, _, _ = top_overlap(
            list(zip(standard["GENE"], standard["standard_z"], strict=False)),
            list(zip(hmagma["GENE"], hmagma["hmagma_z"], strict=False)),
            k=k,
        )
        native_top_metrics[str(k)] = {
            "overlap": overlap,
            "jaccard": jaccard,
            "fraction_of_available_top_k": overlap / min(k, len(standard), len(hmagma)),
        }

    mapping = _symbol_to_versionless_ids(gencode_path)
    reverse_mapping: dict[str, list[str]] = {}
    for symbol, gene_ids in mapping.items():
        for gene_id in gene_ids:
            reverse_mapping.setdefault(gene_id, []).append(symbol)
    known_ids = {gene_id for symbol in KNOWN_AD_GENES for gene_id in mapping.get(symbol, [])}
    known = merged.loc[merged["GENE"].isin(known_ids)].copy()
    known["gene_symbol"] = known["GENE"].map(
        lambda gene_id: ",".join(sorted(reverse_mapping.get(gene_id, [])))
    )
    known = known.sort_values("standard_z", ascending=False)
    write_table(known, output_dir / "known_ad_genes_magma.tsv")
    write_table(merged, output_dir / "standard_vs_hmagma_genes.tsv.gz")
    payload = {
        "n_genes_standard": int(len(standard)),
        "n_genes_hmagma": int(len(hmagma)),
        "n_matched_genes_with_snps": int(len(merged)),
        "zstat_spearman": rho_z,
        "zstat_bootstrap_95ci": list(ci_z),
        "p_spearman": rho_p,
        "top_gene_overlap": native_top_metrics,
        "top_gene_overlap_within_common_genes": top_metrics,
        "top_gene_overlap_background": "Primary top-k lists use each model's full gene universe; common-universe results are also retained.",
        "known_ad_gene_rows": int(len(known)),
        "known_ad_gene_output": str((output_dir / "known_ad_genes_magma.tsv").resolve()),
        "bootstrap_note": (
            "Percentile bootstrap resamples genes and is descriptive; shared LD "
            "and annotation dependence mean it must not be interpreted as an "
            "independent-replicate confidence interval."
        ),
    }
    write_json(payload, output_dir / "magma_gene_concordance.json")
    return payload


def _load_celltogwas_table(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    required = {"cell_type", "z_score", "empirical_p", "fdr_q"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    return frame


def compare_cell_type_methods(
    standard_gsa_path: str | Path,
    hmagma_gsa_path: str | Path,
    celltogwas_dir: Path,
    *,
    output_dir: Path,
    bootstrap_replicates: int,
    random_state: int,
) -> dict[str, object]:
    """Compare cell-type rankings across celltogwas, MAGMA, and H-MAGMA."""

    standard = read_magma_gene_set_results(standard_gsa_path).rename(
        columns={"magma_z": "standard_z", "magma_p": "standard_p", "n_genes": "standard_n_genes"}
    )
    hmagma = read_magma_gene_set_results(hmagma_gsa_path).rename(
        columns={"magma_z": "hmagma_z", "magma_p": "hmagma_p", "n_genes": "hmagma_n_genes"}
    )
    standard["standard_q"] = benjamini_hochberg(standard["standard_p"].to_numpy())
    hmagma["hmagma_q"] = benjamini_hochberg(hmagma["hmagma_p"].to_numpy())
    primary = _load_celltogwas_table(celltogwas_dir / "primary_weighted.tsv").rename(
        columns={
            "z_score": "celltogwas_primary_z",
            "empirical_p": "celltogwas_primary_p",
            "fdr_q": "celltogwas_primary_q",
        }
    )
    membership = _load_celltogwas_table(celltogwas_dir / "membership_aligned.tsv").rename(
        columns={
            "z_score": "celltogwas_membership_z",
            "empirical_p": "celltogwas_membership_p",
            "fdr_q": "celltogwas_membership_q",
        }
    )
    merged = primary.merge(
        membership.loc[
            :,
            [
                "cell_type",
                "celltogwas_membership_z",
                "celltogwas_membership_p",
                "celltogwas_membership_q",
            ],
        ],
        on="cell_type",
        how="inner",
        validate="one_to_one",
    )
    all_snps = _load_celltogwas_table(celltogwas_dir / "sensitivity_all_snps.tsv").rename(
        columns={
            "z_score": "celltogwas_all_snps_z",
            "empirical_p": "celltogwas_all_snps_p",
            "fdr_q": "celltogwas_all_snps_q",
        }
    )
    merged = merged.merge(
        all_snps[
            ["cell_type", "celltogwas_all_snps_z", "celltogwas_all_snps_p", "celltogwas_all_snps_q"]
        ],
        on="cell_type",
        how="inner",
        validate="one_to_one",
    )
    merged = merged.merge(standard, on="cell_type", how="inner", validate="one_to_one")
    merged = merged.merge(hmagma, on="cell_type", how="inner", validate="one_to_one")
    score_columns = (
        "celltogwas_primary_z",
        "celltogwas_membership_z",
        "celltogwas_all_snps_z",
        "standard_z",
        "hmagma_z",
    )
    for column in score_columns:
        rank_column = column.replace("_z", "_rank")
        merged[rank_column] = merged[column].rank(ascending=False, method="min").astype(int)
    merged = merged.sort_values("celltogwas_primary_rank")
    write_table(merged, output_dir / "cell_type_rank_comparison.tsv")

    comparisons = (
        ("primary_weighted_vs_standard", "celltogwas_primary_z", "standard_z"),
        ("primary_weighted_vs_hmagma", "celltogwas_primary_z", "hmagma_z"),
        ("membership_aligned_vs_standard", "celltogwas_membership_z", "standard_z"),
        ("membership_aligned_vs_hmagma", "celltogwas_membership_z", "hmagma_z"),
        ("all_snps_sensitivity_vs_standard", "celltogwas_all_snps_z", "standard_z"),
        ("all_snps_sensitivity_vs_hmagma", "celltogwas_all_snps_z", "hmagma_z"),
        ("standard_vs_hmagma", "standard_z", "hmagma_z"),
        ("primary_vs_membership_aligned", "celltogwas_primary_z", "celltogwas_membership_z"),
    )
    rows: list[dict[str, object]] = []
    for name, left_column, right_column in comparisons:
        finite = merged.loc[np.isfinite(merged[left_column]) & np.isfinite(merged[right_column])]
        rho = spearman_rank(finite[left_column], finite[right_column])
        ci = bootstrap_spearman_ci(
            finite[left_column],
            finite[right_column],
            replicates=bootstrap_replicates,
            random_state=random_state,
        )
        row: dict[str, object] = {
            "comparison": name,
            "n_cell_types": int(len(finite)),
            "spearman_rho": rho,
            "bootstrap_ci_lower": ci[0],
            "bootstrap_ci_upper": ci[1],
        }
        for k in (3, 5):
            overlap, jaccard, left_top, right_top = top_overlap(
                list(zip(finite["cell_type"], finite[left_column], strict=False)),
                list(zip(finite["cell_type"], finite[right_column], strict=False)),
                k=k,
            )
            row[f"top{k}_overlap"] = overlap
            row[f"top{k}_jaccard"] = jaccard
            row[f"top{k}_{name}_left"] = ";".join(left_top)
            row[f"top{k}_{name}_right"] = ";".join(right_top)
        rows.append(row)
    comparison_table = pd.DataFrame(rows)
    write_table(comparison_table, output_dir / "method_concordance.tsv")
    top_cell_types = {
        method: merged.sort_values(column, ascending=False)
        .head(5)
        .loc[:, ["cell_type", column]]
        .to_dict(orient="records")
        for method, column in (
            ("celltogwas_primary", "celltogwas_primary_z"),
            ("celltogwas_membership", "celltogwas_membership_z"),
            ("celltogwas_all_snps_sensitivity", "celltogwas_all_snps_z"),
            ("standard_magma", "standard_z"),
            ("hmagma", "hmagma_z"),
        )
    }
    payload = {
        "n_common_cell_types": int(len(merged)),
        "fdr_significant_cell_types": {
            label: sorted(merged.loc[merged[column] < 0.05, "cell_type"].tolist())
            for label, column in (
                ("celltogwas_primary", "celltogwas_primary_q"),
                ("celltogwas_membership", "celltogwas_membership_q"),
                ("celltogwas_all_snps_sensitivity", "celltogwas_all_snps_q"),
                ("standard_magma", "standard_q"),
                ("hmagma", "hmagma_q"),
            )
        },
        "comparisons": comparison_table.to_dict(orient="records"),
        "top_cell_types": top_cell_types,
        "rank_interpretation": (
            "Cell-type ranking agreement assesses numerical consistency and "
            "biological plausibility only; it is not causal evidence."
        ),
    }
    write_json(payload, output_dir / "cell_type_method_concordance.json")
    return payload


def compare_all(
    config: ValidationConfig,
    directories: Mapping[str, Path],
) -> dict[str, object]:
    """Run gene-level and cell-type-level concordance analyses."""

    standard_genes = directories["magma"] / "standard_gene.genes.out.txt"
    hmagma_genes = directories["magma"] / "hmagma_gene.genes.out.txt"
    standard_gsa = directories["magma"] / "standard_geneset.gsa.out.txt"
    hmagma_gsa = directories["magma"] / "hmagma_geneset.gsa.out.txt"
    _require_paths([standard_genes, hmagma_genes, standard_gsa, hmagma_gsa])
    gene_metrics = compare_magma_gene_results(
        standard_genes,
        hmagma_genes,
        output_dir=directories["comparison"],
        gencode_path=Path(config.gencode_path),
        bootstrap_replicates=config.bootstrap_replicates,
        random_state=config.random_state,
    )
    cell_type_metrics = compare_cell_type_methods(
        standard_gsa,
        hmagma_gsa,
        directories["celltogwas"],
        output_dir=directories["comparison"],
        bootstrap_replicates=config.bootstrap_replicates,
        random_state=config.random_state,
    )
    payload = {
        "standard_vs_hmagma_genes": gene_metrics,
        "cell_type_comparisons": cell_type_metrics,
        "limits": [
            "Rank agreement is not numerical equivalence.",
            "MAGMA models LD; celltogwas uses circular annotation shifts and does not estimate LD.",
            "MAGMA and H-MAGMA use different SNP-to-gene models and therefore different gene-level backgrounds.",
            "GWAS is GRCh38; all methods use rsID-matched GRCh37 1000G EUR positions without liftover.",
            "scRNA-seq donor coverage differs by cell type; only one donor occurs in all 24 brainSCOPE classes.",
        ],
    }
    write_json(payload, directories["comparison"] / "validation_summary.json")
    return payload
