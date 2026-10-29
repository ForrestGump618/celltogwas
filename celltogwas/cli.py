"""Command-line interface for celltogwas."""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from ._version import __version__
from .adaptive import run_adaptive_gene_enrichment
from .exceptions import CellToGWASError
from .io import (
    file_sha256,
    load_expression_matrix,
    read_table,
    resolve_column,
    write_gene_sets,
    write_json,
    write_table,
)
from .linker import GenomeLinker
from .models import LinkConfig
from .specificity import (
    compute_specificity,
    specificity_summary,
    specificity_to_gene_sets,
)


def build_parser() -> argparse.ArgumentParser:
    """Build the public CLI parser."""

    parser = argparse.ArgumentParser(
        prog="celltogwas",
        description=(
            "Regulatory-aware mapping of GWAS loci to cell types using "
            "weighted gene programs and spatial permutation tests."
        ),
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    annotate = subparsers.add_parser(
        "annotate",
        help="Map GWAS positions to weighted cell-type domains.",
    )
    _add_linker_arguments(annotate, gwas_required=True)
    annotate.add_argument(
        "--mode",
        choices=("weight", "binary", "both"),
        default="weight",
        help="Continuous weights, binary labels, or both.",
    )
    annotate.add_argument(
        "--threshold",
        type=float,
        default=0.05,
        help="Positive-weight threshold used by binary mode.",
    )

    enrich = subparsers.add_parser(
        "enrich",
        help="Run within-chromosome circular-permutation cell-type enrichment.",
    )
    _add_linker_arguments(enrich, gwas_required=True)
    enrich.add_argument("--p-column", help="P-value column name.")
    enrich.add_argument("--score-column", help="Effect/statistic column name.")
    enrich.add_argument(
        "--score-mode",
        choices=("pvalue", "uniform", "zscore", "abs"),
        default="pvalue",
        help="Evidence-weight construction.",
    )
    enrich.add_argument(
        "--permutations",
        type=int,
        default=10_000,
        help="Number of circular-shift permutations.",
    )
    enrich.add_argument("--seed", type=int, default=13, help="Random seed.")
    enrich.add_argument(
        "--min-shift-fraction",
        type=float,
        default=0.01,
        help="Minimum circular offset as a fraction of chromosome variants.",
    )
    enrich.add_argument(
        "--p-upper",
        type=float,
        help="Give evidence only to SNPs at or below this P; retain all valid background SNPs.",
    )
    enrich.add_argument(
        "--weight-floor",
        type=float,
        default=0.0,
        help="Subtract this value from evidence weights before testing.",
    )
    enrich.add_argument(
        "--snps-output",
        help="Optional output with per-SNP cell-type link weights.",
    )

    specificity = subparsers.add_parser(
        "specificity",
        help="Derive marker weights from a genes-by-cell-types expression matrix.",
    )
    specificity.add_argument("--expression", required=True)
    specificity.add_argument("--output", required=True)
    specificity.add_argument(
        "--method",
        choices=("tau", "entropy", "zscore"),
        default="tau",
    )
    specificity.add_argument("--min-expression", type=float, default=0.0)
    specificity.add_argument("--power", type=float, default=1.0)
    specificity.add_argument("--top-n", type=int)
    specificity.add_argument("--min-weight", type=float, default=0.0)
    specificity.add_argument("--gmt-output")
    specificity.add_argument("--summary-output")
    specificity.add_argument("--report")
    specificity.add_argument("--quiet", action="store_true")

    validate = subparsers.add_parser(
        "validate",
        help="Validate annotations, gene matching and optional GWAS columns.",
    )
    _add_linker_arguments(validate, gwas_required=False)
    validate.add_argument(
        "--gwas-rows",
        type=int,
        default=10_000,
        help="Rows used for optional GWAS validation.",
    )
    adaptive = subparsers.add_parser(
        "adaptive-gene",
        help="Selection-adjusted multiscale gene-program enrichment (experimental).",
    )
    adaptive.add_argument(
        "--genes", required=True, help="Gene table: gene,chrom,pos,z and cell weights."
    )
    adaptive.add_argument("--cell-types", required=True, nargs="+")
    adaptive.add_argument("--covariates", nargs="+", default=[])
    adaptive.add_argument("--permutations", type=int, default=9_999)
    adaptive.add_argument("--seed", type=int, default=13)
    adaptive.add_argument("--output", required=True)
    adaptive.add_argument("--components-output")
    adaptive.add_argument("--report")
    return parser


def _add_linker_arguments(
    parser: argparse.ArgumentParser,
    *,
    gwas_required: bool,
) -> None:
    parser.add_argument("--gtf", required=True, help="GTF/GFF3 gene annotation.")
    parser.add_argument("--gene-sets", required=True, help="GMT/JSON/CSV gene sets.")
    parser.add_argument("--gwas", required=gwas_required, help="GWAS summary table.")
    if gwas_required:
        parser.add_argument("--output", required=True, help="Output table path.")
    parser.add_argument("--regulatory", help="Optional BED-like regulatory intervals.")
    parser.add_argument("--regulatory-cell-type")
    parser.add_argument("--window", type=int, default=100_000)
    parser.add_argument("--promoter-window", type=int, default=2_000)
    parser.add_argument("--promoter-weight", type=float, default=1.0)
    parser.add_argument("--gene-body-weight", type=float, default=1.0)
    parser.add_argument(
        "--kernel",
        choices=("uniform", "exponential"),
        default="exponential",
    )
    parser.add_argument("--decay-scale", type=int)
    parser.add_argument("--kernel-levels", type=int, default=6)
    parser.add_argument("--max-weight", type=float, default=1.0)
    parser.add_argument("--chunk-size", type=int, default=500_000)
    parser.add_argument(
        "--all-features",
        action="store_true",
        help="Retain all GTF/GFF features rather than only feature_type=gene.",
    )
    parser.add_argument("--strict-annotation", action="store_true")
    parser.add_argument("--report", help="Write a JSON run/build report.")
    parser.add_argument("--quiet", action="store_true")


def _config_from_args(args: argparse.Namespace) -> LinkConfig:
    return LinkConfig(
        window=args.window,
        promoter_window=args.promoter_window,
        promoter_weight=args.promoter_weight,
        gene_body_weight=args.gene_body_weight,
        kernel=args.kernel,
        decay_scale=args.decay_scale,
        kernel_levels=args.kernel_levels,
        max_weight=args.max_weight,
    )


def _runtime_payload() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }


def _build_linker(args: argparse.Namespace) -> GenomeLinker:
    return GenomeLinker.from_paths(
        args.gtf,
        args.gene_sets,
        config=_config_from_args(args),
        regulatory_path=args.regulatory,
        regulatory_cell_type=args.regulatory_cell_type,
        feature_types=None if args.all_features else ("gene",),
        strict_annotation=args.strict_annotation,
    )


def _report_payload(
    args: argparse.Namespace,
    *,
    command: str,
    started: float,
    linker: GenomeLinker,
    output: str | None = None,
    extra: dict[str, object] | None = None,
) -> dict[str, object]:
    inputs = [args.gtf, args.gene_sets]
    if getattr(args, "gwas", None):
        inputs.append(args.gwas)
    if getattr(args, "regulatory", None):
        inputs.append(args.regulatory)
    payload: dict[str, object] = {
        "package": "celltogwas",
        "version": __version__,
        "command": command,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.perf_counter() - started,
        "runtime": _runtime_payload(),
        "inputs": {str(Path(path).resolve()): file_sha256(path) for path in inputs if path},
        "output": str(Path(output).resolve()) if output else None,
        "linker": linker.manifest(),
    }
    if extra:
        payload.update(extra)
    return payload


def _run_annotate(args: argparse.Namespace) -> int:
    started = time.perf_counter()
    gwas = read_table(args.gwas)
    linker = _build_linker(args)
    annotated = linker.annotate(
        gwas,
        mode=args.mode,
        threshold=args.threshold,
        chunk_size=args.chunk_size,
    )
    write_table(annotated, args.output)
    if args.report:
        write_json(
            _report_payload(
                args,
                command="annotate",
                started=started,
                linker=linker,
                output=args.output,
                extra={
                    "annotation_mode": args.mode,
                    "threshold": args.threshold,
                    "n_rows": int(len(annotated)),
                },
            ),
            args.report,
        )
    if not args.quiet:
        print(
            f"Wrote {len(annotated):,} annotated rows and "
            f"{len(linker.cell_types)} cell-type columns to {args.output}",
            file=sys.stderr,
        )
    return 0


def _run_enrich(args: argparse.Namespace) -> int:
    started = time.perf_counter()
    gwas = read_table(args.gwas)
    linker = _build_linker(args)
    result = linker.enrich(
        gwas,
        p_column=args.p_column,
        score_column=args.score_column,
        score_mode=args.score_mode,
        n_permutations=args.permutations,
        random_state=args.seed,
        min_shift_fraction=args.min_shift_fraction,
        p_upper=args.p_upper,
        weight_floor=args.weight_floor,
        chunk_size=args.chunk_size,
    )
    write_table(result, args.output)
    if args.snps_output:
        write_table(linker.annotate(gwas, chunk_size=args.chunk_size), args.snps_output)
    if args.report:
        write_json(
            _report_payload(
                args,
                command="enrich",
                started=started,
                linker=linker,
                output=args.output,
                extra={
                    "enrichment": result.attrs.get("diagnostics", {}),
                    "n_rows": int(len(result)),
                },
            ),
            args.report,
        )
    if not args.quiet:
        preview = result.loc[
            :,
            ["cell_type", "observed_link_score", "z_score", "empirical_p", "fdr_q"],
        ].head(10)
        print(preview.to_string(index=False), file=sys.stderr)
    return 0


def _run_specificity(args: argparse.Namespace) -> int:
    started = time.perf_counter()
    expression = load_expression_matrix(args.expression)
    result = compute_specificity(
        expression,
        method=args.method,
        min_expression=args.min_expression,
        power=args.power,
        top_n=args.top_n,
        min_weight=args.min_weight,
    )
    write_table(result, args.output)
    if args.gmt_output:
        write_gene_sets(
            specificity_to_gene_sets(result),
            args.gmt_output,
        )
    summary_table = specificity_summary(result)
    if args.summary_output:
        write_table(summary_table, args.summary_output)
    if args.report:
        write_json(
            {
                "package": "celltogwas",
                "version": __version__,
                "command": "specificity",
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "elapsed_seconds": time.perf_counter() - started,
                "runtime": _runtime_payload(),
                "input": {str(Path(args.expression).resolve()): file_sha256(args.expression)},
                "output": str(Path(args.output).resolve()),
                "method": args.method,
                "parameters": {
                    "min_expression": args.min_expression,
                    "power": args.power,
                    "top_n": args.top_n,
                    "min_weight": args.min_weight,
                },
                "outputs": {
                    "markers": str(Path(args.output).resolve()),
                    "gmt": (
                        str(Path(args.gmt_output).resolve())
                        if args.gmt_output is not None
                        else None
                    ),
                    "summary": (
                        str(Path(args.summary_output).resolve())
                        if args.summary_output is not None
                        else None
                    ),
                },
                "n_cell_types": int(expression.shape[1]),
                "n_genes": int(expression.shape[0]),
                "n_markers": int(len(result)),
                "summary": summary_table.to_dict(orient="records"),
            },
            args.report,
        )
    if not args.quiet:
        print(
            f"Wrote {len(result):,} marker weights for "
            f"{expression.shape[1]} cell types to {args.output}",
            file=sys.stderr,
        )
    return 0


def _run_validate(args: argparse.Namespace) -> int:
    started = time.perf_counter()
    if args.gwas_rows < 1:
        raise ValueError("gwas_rows must be positive.")
    linker = _build_linker(args)
    inputs = [args.gtf, args.gene_sets, args.regulatory, args.gwas]
    payload = {
        "package": "celltogwas",
        "version": __version__,
        "command": "validate",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": None,
        "runtime": _runtime_payload(),
        "valid": True,
        "inputs": {
            str(Path(path).resolve()): file_sha256(path) for path in inputs if path is not None
        },
        "linker": linker.manifest(),
    }
    if args.gwas:
        gwas = read_table(args.gwas, nrows=args.gwas_rows)
        resolve_column(
            gwas.columns,
            ("chr", "chrom", "chromosome", "chr_id", "contig"),
            label="chromosome",
        )
        resolve_column(
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
        annotated = linker.annotate(gwas, chunk_size=args.chunk_size)
        coverage = {
            cell_type: float((annotated[cell_type] > 0).mean()) for cell_type in linker.cell_types
        }
        payload["gwas_validation"] = {
            "path": str(Path(args.gwas).resolve()),
            "n_rows_checked": int(len(gwas)),
            "link_coverage": coverage,
        }
    payload["elapsed_seconds"] = time.perf_counter() - started
    if args.report:
        write_json(payload, args.report)
    elif not args.quiet:
        print(json.dumps(payload, ensure_ascii=False, indent=2, default=str))
    return 0


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""

    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "adaptive-gene":
            result = run_adaptive_gene_enrichment(
                read_table(args.genes),
                cell_types=args.cell_types,
                covariates=args.covariates,
                n_permutations=args.permutations,
                random_state=args.seed,
            )
            write_table(result.table, args.output)
            if args.components_output:
                write_table(result.components, args.components_output)
            if args.report:
                write_json(
                    {
                        "version": __version__,
                        "command": "adaptive-gene",
                        "input_sha256": file_sha256(args.genes),
                        "runtime": _runtime_payload(),
                        **result.diagnostics,
                    },
                    args.report,
                )
            return 0
        if args.command == "annotate":
            return _run_annotate(args)
        if args.command == "enrich":
            return _run_enrich(args)
        if args.command == "specificity":
            return _run_specificity(args)
        if args.command == "validate":
            return _run_validate(args)
        parser.error(f"Unknown command: {args.command}")
    except (CellToGWASError, OSError, ValueError) as error:
        print(f"celltogwas: error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
