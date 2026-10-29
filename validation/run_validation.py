"""Command-line entry point for the celltogwas external validation."""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

from celltogwas._version import __version__
from celltogwas.io import write_json
from validation.pipeline import (
    ValidationConfig,
    build_expression_and_markers,
    compare_all,
    harmonize_gwas,
    run_celltogwas_analyses,
    run_magma_analyses,
    sha256_file,
    utc_now,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run real scRNA-seq/GWAS validation and MAGMA/H-MAGMA "
            "method-consistency analyses for celltogwas."
        )
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        required=True,
        help="Directory containing gwas/, psychencode/, gencode/, magma/, and hmagma/.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for machine-readable validation results and logs.",
    )
    parser.add_argument("--top-n", type=int, default=1_000)
    parser.add_argument("--min-expression", type=float, default=0.1)
    parser.add_argument("--min-detection-fraction", type=float, default=0.1)
    parser.add_argument("--link-window", type=int, default=10_000)
    parser.add_argument("--primary-p-upper", type=float, default=5e-8)
    parser.add_argument("--sensitivity-p-upper", type=float, default=1e-5)
    parser.add_argument("--permutations", type=int, default=10_000)
    parser.add_argument("--bootstrap", type=int, default=2_000)
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--magma-executable", type=Path)
    parser.add_argument("--magma-reference-prefix", type=Path)
    parser.add_argument("--gwas", type=Path)
    parser.add_argument("--psychencode-dir", type=Path)
    parser.add_argument("--gencode", type=Path)
    parser.add_argument("--hmagma-annotation", type=Path)
    parser.add_argument("--standard-annotation", type=Path)
    parser.add_argument(
        "--resume-prepared",
        action="store_true",
        help="Reuse existing expression and GWAS preparation outputs.",
    )
    return parser


def _validate_prepared_inputs(
    config: ValidationConfig,
    directories: dict[str, Path],
) -> tuple[dict[str, object], dict[str, object]]:
    expression_manifest = directories["expression"] / "expression_manifest.json"
    gwas_manifest = directories["harmonized"] / "gwas_harmonization.json"
    required = [
        expression_manifest,
        gwas_manifest,
        directories["gene_sets"] / "markers_weighted.gmt",
        directories["gene_sets"] / "markers_membership.gmt",
        directories["gene_sets"] / "markers_membership_ensembl.gmt",
        directories["harmonized"] / "gwas_harmonized.tsv.gz",
        directories["harmonized"] / "gwas_for_magma.tsv",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Prepared outputs are incomplete:\n" + "\n".join(missing))
    expression = json.loads(expression_manifest.read_text(encoding="utf-8"))
    gwas = json.loads(gwas_manifest.read_text(encoding="utf-8"))
    parameters = expression["marker_parameters"]
    for name in ("top_n", "min_expression", "min_detection_fraction"):
        if parameters[name] != getattr(config, name):
            raise ValueError(
                f"Prepared marker parameter {name} differs from the requested configuration."
            )
    expected = {"gwas": Path(config.gwas_path), "bim": Path(f"{config.magma_reference_prefix}.bim")}
    for name, path in expected.items():
        if gwas["inputs_sha256"][name] != sha256_file(path):
            raise ValueError(f"Prepared GWAS input fingerprint changed: {name}")
    return expression, gwas


def _input_hashes(config: ValidationConfig) -> dict[str, str]:
    files = {
        "gwas": Path(config.gwas_path),
        "gencode": Path(config.gencode_path),
        "magma_executable": Path(config.magma_executable),
        "standard_annotation": Path(config.standard_annotation_path),
        "hmagma_annotation": Path(config.hmagma_annotation_path),
        "magma_reference_bed": Path(f"{config.magma_reference_prefix}.bed"),
        "magma_reference_bim": Path(f"{config.magma_reference_prefix}.bim"),
        "magma_reference_fam": Path(f"{config.magma_reference_prefix}.fam"),
    }
    return {label: sha256_file(path) for label, path in files.items()}


def run(args: argparse.Namespace) -> dict[str, object]:
    started = time.perf_counter()
    config = ValidationConfig(
        data_root=args.data_root,
        output_dir=args.output_dir,
        top_n=args.top_n,
        min_expression=args.min_expression,
        min_detection_fraction=args.min_detection_fraction,
        link_window=args.link_window,
        primary_p_upper=args.primary_p_upper,
        sensitivity_p_upper=args.sensitivity_p_upper,
        n_permutations=args.permutations,
        random_state=args.seed,
        bootstrap_replicates=args.bootstrap,
        magma_executable=args.magma_executable,
        magma_reference_prefix=args.magma_reference_prefix,
        gwas_path=args.gwas,
        psychencode_dir=args.psychencode_dir,
        gencode_path=args.gencode,
        hmagma_annotation_path=args.hmagma_annotation,
        standard_annotation_path=args.standard_annotation,
    ).resolved_paths()
    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
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

    manifest_path = output_dir / "validation_manifest.json"
    manifest: dict[str, object] = {
        "package": "celltogwas",
        "package_version": __version__,
        "started_at_utc": utc_now(),
        "completed_at_utc": None,
        "status": "running",
        "runtime": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "analysis_code_sha256": {
            str(path.relative_to(Path(__file__).resolve().parents[1])).replace(
                "\\", "/"
            ): sha256_file(path)
            for folder in ("celltogwas", "validation")
            for path in sorted((Path(__file__).resolve().parents[1] / folder).glob("*.py"))
        },
        "config": {
            key: str(value) if isinstance(value, Path) else value
            for key, value in asdict(config).items()
        },
        "prespecified_analysis": {
            "primary_p_upper": config.primary_p_upper,
            "primary_weight_floor": -np.log10(config.primary_p_upper),
            "sensitivity_p_upper": config.sensitivity_p_upper,
            "n_permutations": config.n_permutations,
            "random_state": config.random_state,
            "marker_top_n_per_cell_type": config.top_n,
            "marker_min_expression": config.min_expression,
            "marker_min_detection_fraction": config.min_detection_fraction,
            "celltogwas_link_definition": "GENCODE v19 gene body plus bilateral 10 kb",
        },
        "controls": {
            "primary_evidence_threshold": config.primary_p_upper,
            "report_significance_at_fdr_q": 0.05,
            "multiple_testing": "Benjamini-Hochberg across 24 cell types",
            "mhc_sensitivity": "Exclude GRCh37 chr6:25000000-35000000",
            "gwas_harmonization": "Exact rsID match to 1000G EUR BIM; no liftover",
        },
    }
    write_json(manifest, manifest_path)

    try:
        print("[1/5] Fingerprinting inputs", file=sys.stderr)
        manifest["input_sha256"] = _input_hashes(config)
        write_json(manifest, manifest_path)

        if (
            args.resume_prepared
            and (directories["expression"] / "expression_manifest.json").exists()
        ):
            print("[2/5] Reusing prepared expression and GWAS tables", file=sys.stderr)
            expression_manifest, gwas_manifest = _validate_prepared_inputs(config, directories)
        else:
            print("[2/5] Building expression profiles and marker sets", file=sys.stderr)
            expression_manifest = build_expression_and_markers(config, directories)
            print("[3/5] Harmonizing GWAS by rsID to the GRCh37 reference", file=sys.stderr)
            gwas_manifest = harmonize_gwas(
                Path(config.gwas_path),
                Path(f"{config.magma_reference_prefix}.bim"),
                directories["harmonized"],
            )
        manifest["expression"] = expression_manifest
        manifest["gwas_harmonization"] = gwas_manifest
        write_json(manifest, manifest_path)

        print("[4/5] Running MAGMA, H-MAGMA, and celltogwas analyses", file=sys.stderr)
        with ThreadPoolExecutor(max_workers=2) as executor:
            magma = executor.submit(run_magma_analyses, config, directories)
            celltogwas = executor.submit(run_celltogwas_analyses, config, directories)
            manifest["magma"] = magma.result()
            write_json(manifest, manifest_path)
            manifest["celltogwas"] = celltogwas.result()

        print("[5/5] Computing method concordance", file=sys.stderr)
        manifest["comparisons"] = compare_all(config, directories)
        manifest["status"] = "completed"
        manifest["completed_at_utc"] = utc_now()
        manifest["elapsed_seconds"] = time.perf_counter() - started
        write_json(manifest, manifest_path)
    except Exception as error:
        manifest["status"] = "failed"
        manifest["failed_at_utc"] = utc_now()
        manifest["error"] = f"{type(error).__name__}: {error}"
        manifest["elapsed_seconds"] = time.perf_counter() - started
        write_json(manifest, manifest_path)
        raise
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        manifest = run(args)
    except Exception as error:
        print(f"validation failed: {type(error).__name__}: {error}", file=sys.stderr)
        return 1
    print(
        "Validation completed: "
        f"{manifest['package_version']} with "
        f"{manifest['expression']['cell_type_count']} cell types and "
        f"{manifest['gwas_harmonization']['gwas']['matched_rsids']:,} matched SNPs.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
