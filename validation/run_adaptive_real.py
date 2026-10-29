"""Frozen multiscale reanalysis of the existing public AD/brainSCOPE evidence.

This is NOT an independent replication dataset. Vendor gene results are inputs;
MAGMA gene-level LD calculations are not replaced by the circular test.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from celltogwas import (
    benjamini_hochberg,
    file_sha256,
    load_gene_annotation,
    run_adaptive_gene_enrichment,
)
from celltogwas.adaptive import _holm

COVARIATES = [
    "log_length",
    "log_nsnp",
    "log_nparam",
    "log_n",
    "mean_expression",
    "log_length_sq",
    "log_nsnp_sq",
    "mean_expression_sq",
]


def prepare_inputs(annotation, expression, detection, donors, gene_results):
    """One-to-one ID/symbol mapping and a shared fully measured gene universe."""
    cells = list(expression.columns)
    if expression.index.has_duplicates or expression.columns.has_duplicates:
        raise ValueError("Expression gene symbols/cell types must be unique.")
    if not detection.index.equals(expression.index) or not donors.index.equals(expression.index):
        raise ValueError("Detection/donor rows must match expression rows exactly.")
    if list(detection.columns) != cells or list(donors.columns) != cells:
        raise ValueError("Detection/donor columns must match expression columns exactly.")
    if len(cells) < 2:
        raise ValueError("At least two cell types are required for Tau specificity.")
    for name, matrix in [("expression", expression), ("detection", detection), ("donors", donors)]:
        if not np.isfinite(matrix.to_numpy()).all() or (matrix < 0).any().any():
            raise ValueError(f"{name} must be finite and non-negative; no implicit imputation.")
    if (detection > 1).any().any() or (donors % 1 != 0).any().any():
        raise ValueError("Detection fractions must be <=1 and donor counts must be integers.")
    background = expression.loc[(donors > 0).all(axis=1) & (expression.max(axis=1) >= 0.1)]
    if not np.isfinite(background.to_numpy()).all() or (background < 0).any().any():
        raise ValueError("Measured expression profiles must be finite and non-negative.")
    annot = annotation.copy()
    annot["gene"] = annot.gene_id.str.split(".").str[0]
    annot = annot[annot.chrom.astype(str).isin([str(i) for i in range(1, 23)])]
    ambiguous = annot.gene_name.duplicated(keep=False) | annot.gene.duplicated(keep=False)
    annot = annot.loc[~ambiguous & annot.gene_name.isin(background.index)].copy()
    common = set(annot.gene)
    for frame in gene_results.values():
        if frame.GENE.duplicated().any():
            raise ValueError("MAGMA gene IDs must be unique.")
        common &= set(frame.GENE)
    annot = annot[annot.gene.isin(common)].sort_values("gene").reset_index(drop=True)
    if len(annot) < 100:
        raise ValueError("Too few common measured genes.")
    profiles = background.loc[annot.gene_name].to_numpy()
    relative = np.divide(
        profiles,
        profiles.max(axis=1, keepdims=True),
        out=np.zeros_like(profiles),
        where=profiles.max(axis=1, keepdims=True) > 0,
    )
    tau = (1 - relative).sum(axis=1) / (len(cells) - 1)
    weights = relative * tau[:, None]
    weights *= detection.loc[annot.gene_name, cells].to_numpy() >= 0.1
    inputs = {}
    for method, result in gene_results.items():
        aligned = result.set_index("GENE").loc[annot.gene]
        frame = pd.DataFrame(
            {
                "gene": annot.gene.to_numpy(),
                "symbol": annot.gene_name.to_numpy(),
                "chrom": annot.chrom.to_numpy(),
                "pos": annot.start.to_numpy(),
                "z": aligned.ZSTAT.to_numpy(),
                "gene_p": aligned.P.to_numpy(),
                "log_length": np.log1p((annot.end - annot.start).to_numpy()),
                "log_nsnp": np.log1p(aligned.NSNPS.to_numpy()),
                "log_nparam": np.log1p(aligned.NPARAM.to_numpy()),
                "log_n": np.log1p(aligned.N.to_numpy()),
                "mean_expression": profiles.mean(axis=1),
            }
        )
        for column in ["log_length", "log_nsnp", "mean_expression"]:
            frame[f"{column}_sq"] = frame[column] ** 2
        frame[cells] = weights
        inputs[method] = frame
    return (
        inputs,
        cells,
        {
            "measured_symbols_before_mapping": len(background),
            "common_genes": len(annot),
            "excluded_ambiguous_autosomal_annotation_records": int(ambiguous.sum()),
            "mapping": "one-to-one GENCODEv19 stable ID/symbol; intersection of both gene-result universes",
            "missingness": "require donor count >0 in all 24 classes; do not impute unmeasured profiles",
        },
    )


def run(data: Path, results: Path, output: Path, benchmark: Path, permutations=19999):
    output.mkdir(parents=True, exist_ok=True)
    files = {
        "gtf": data / "gencode/gencode.v19.annotation.gtf.gz",
        "expression": results / "expression/cell_type_mean_expression.tsv.gz",
        "detection": results / "expression/cell_type_detection_fraction.tsv.gz",
        "donors": results / "expression/cell_type_donors.tsv.gz",
        "standard": results / "magma/standard_batched.genes.out.txt",
        "hmagma": results / "magma/hmagma_batched.genes.out.txt",
        "benchmark_summary": benchmark / "summary.json",
    }
    benchmark_summary = json.loads(files["benchmark_summary"].read_text())
    if not benchmark_summary["null_gate_passed"]:
        raise ValueError("Simulation calibration gate failed; stop before real-data reanalysis.")
    source = Path(__file__).parents[1] / "celltogwas/adaptive.py"
    benchmark_protocol = json.loads((benchmark / "protocol.json").read_text())
    if benchmark_protocol["module_sha256_before_run"] != file_sha256(source):
        raise ValueError("Algorithm changed after calibration; rerun the fixed benchmark.")
    protocol = {
        "analysis_date_local": "2026-10-30",
        "status": "exploratory reanalysis, not replication",
        "fixed_before_reanalysis": True,
        "n_permutations": permutations,
        "seed": 13,
        "powers": [1, 2, 4],
        "thresholds": [None, 2, 4],
        "binary_fraction": 0.1,
        "covariates": COVARIATES,
        "primary_family": "24 classes x 2 mappings = 48 tests; global BY/Holm",
        "sensitivity_families": ["exclude MHC chr6:25-35Mb", "exclude chromosome 19"],
        "sensitivity_role": "robustness only; not replacements for primary or independent replications",
        "algorithm_sha256": file_sha256(source),
        "runner_sha256": file_sha256(__file__),
        "inputs": {
            key: {"path": str(path), "sha256": file_sha256(path)} for key, path in files.items()
        },
    }
    (output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n", encoding="utf-8")
    expression, detection, donors = [
        pd.read_csv(files[key], sep="\t", index_col=0)
        for key in ["expression", "detection", "donors"]
    ]
    gene_results = {key: pd.read_csv(files[key], sep=r"\s+") for key in ["standard", "hmagma"]}
    inputs, cells, diagnostics = prepare_inputs(
        load_gene_annotation(files["gtf"]).genes,
        expression,
        detection,
        donors,
        gene_results,
    )
    tables, prioritization = [], []
    for method, frame in inputs.items():
        frame.to_csv(output / f"{method}.input.tsv.gz", sep="\t", index=False)
        for sensitivity in ["primary", "no_mhc", "no_chr19"]:
            eligible = frame.copy()
            if sensitivity == "no_mhc":
                eligible = eligible[
                    ~((eligible.chrom == "6") & eligible.pos.between(25_000_000, 35_000_000))
                ]
            elif sensitivity == "no_chr19":
                eligible = eligible[eligible.chrom != "19"]
            result = run_adaptive_gene_enrichment(
                eligible,
                cell_types=cells,
                covariates=COVARIATES,
                n_permutations=permutations,
                random_state=13,
            )
            stem = f"{method}.{sensitivity}"
            result.table.to_csv(output / f"{stem}.tsv", sep="\t", index=False, lineterminator="\n")
            result.components.to_csv(
                output / f"{stem}.components.tsv", sep="\t", index=False, lineterminator="\n"
            )
            (output / f"{stem}.json").write_text(
                json.dumps(result.diagnostics, indent=2) + "\n", encoding="utf-8"
            )
            table = result.table.assign(method=method, analysis=sensitivity)
            tables.append(table)
            print(f"{stem}: {len(eligible)} genes; min P={table.p_adaptive.min():.6g}", flush=True)
        for cell in cells:
            support = frame[["gene", "symbol", "chrom", "pos", "z", "gene_p", cell]].copy()
            support["descriptive_support"] = support[cell] * np.maximum(support.z, 0)
            support = support.sort_values("descriptive_support", ascending=False).head(20)
            prioritization.append(
                support.rename(columns={cell: "cell_weight"}).assign(
                    method=method,
                    cell_type=cell,
                    status="descriptive prioritization, NOT causal/novel-gene significance",
                )
            )
    combined = pd.concat(tables, ignore_index=True)
    for _, group in combined.groupby("analysis"):
        p = group.p_adaptive.to_numpy()
        combined.loc[group.index, "p_holm_global48"] = _holm(p)
        combined.loc[group.index, "q_by_global48"] = np.minimum(
            benjamini_hochberg(p) * np.sum(1 / np.arange(1, len(p) + 1)),
            1,
        )
    combined.to_csv(output / "combined.tsv", sep="\t", index=False, lineterminator="\n")
    pd.concat(prioritization, ignore_index=True).to_csv(
        output / "gene_prioritization.tsv", sep="\t", index=False, lineterminator="\n"
    )
    diagnostics["input_table_sha256"] = {
        key: file_sha256(output / f"{key}.input.tsv.gz") for key in inputs
    }
    (output / "mapping.json").write_text(json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=19999)
    args = parser.parse_args()
    run(args.data, args.results, args.output, args.benchmark, args.permutations)
