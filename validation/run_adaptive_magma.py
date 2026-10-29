"""Actual LD-aware MAGMA continuous/binary regressions on the shared universe."""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

from celltogwas import benjamini_hochberg, file_sha256
from celltogwas.adaptive import _holm
from validation.run_adaptive_real import COVARIATES


def run(executable: Path, gene_results: Path, adaptive: Path, output: Path):
    output.mkdir(parents=True, exist_ok=True)
    tables, runs = [], []
    for method in ["standard", "hmagma"]:
        frame = pd.read_csv(adaptive / f"{method}.input.tsv.gz", sep="\t")
        metadata = {"gene", "symbol", "chrom", "pos", "z", "gene_p", *COVARIATES}
        cells = [column for column in frame if column not in metadata]
        covariates = frame[["gene", *COVARIATES]].copy()
        for cell in cells:
            weights = frame[cell].to_numpy()
            covariates[f"continuous_{cell}"] = weights
            covariates[f"binary_{cell}"] = (
                (weights >= np.quantile(weights, 0.9)) & (weights > 0)
            ).astype(int)
        covariate_path = output / f"{method}.covariates.tsv"
        covariates.to_csv(covariate_path, sep="\t", index=False, lineterminator="\n")
        raw = gene_results / f"{method}_batched.genes.raw"
        prefix = output / f"{method}.properties"
        command = [
            str(executable),
            "--gene-results",
            str(raw),
            "--gene-covar",
            str(covariate_path),
            "missing-values=drop",
            "--model",
            "direction-covar=greater",
            "condition-residualize=" + ",".join(COVARIATES),
            "--out",
            str(prefix),
        ]
        manifest = {
            "command": command,
            "inputs": {str(p): file_sha256(p) for p in [executable, raw, covariate_path]},
            "same_universe_as_adaptive": len(frame),
            "purpose": "Prespecified continuous and fixed top-10% binary gene-property comparators, not best-kernel selection.",
        }
        manifest_path = output / f"{method}.receipt.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        started = time.perf_counter()
        with (output / f"{method}.console.log").open("w", encoding="utf-8") as stream:
            completed = subprocess.run(
                command, stdout=stream, stderr=subprocess.STDOUT, check=False
            )
        manifest["returncode"] = completed.returncode
        manifest["elapsed_seconds"] = time.perf_counter() - started
        if (
            completed.returncode != 0
            or "End time is" not in prefix.with_suffix(".properties.log").read_text()
        ):
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
            raise RuntimeError(f"MAGMA {method} failed; inspect console log.")
        table = pd.read_csv(Path(str(prefix) + ".gsa.out.txt"), sep=r"\s+", comment="#")
        if len(table) != 2 * len(cells) or not np.isfinite(table.P).all():
            raise RuntimeError("Incomplete gene-property result.")
        table["kind"] = table.VARIABLE.str.split("_").str[0]
        table["cell_type"] = table.VARIABLE.str.replace(r"^(continuous|binary)_", "", regex=True)
        tables.append(table.assign(method=method))
        manifest["completed"] = True
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        runs.append(manifest)
        print(
            f"{method}: {len(table)} actual MAGMA property tests, {manifest['elapsed_seconds']:.1f}s",
            flush=True,
        )
    combined = pd.concat(tables, ignore_index=True)
    p = combined.P.to_numpy()
    combined["p_holm_global96"] = _holm(p)
    combined["q_by_global96"] = np.minimum(
        benjamini_hochberg(p) * np.sum(1 / np.arange(1, len(p) + 1)), 1
    )
    combined.to_csv(output / "combined.tsv", sep="\t", index=False, lineterminator="\n")
    (output / "runs.json").write_text(json.dumps(runs, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--gene-results", type=Path, required=True)
    parser.add_argument("--adaptive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(
        args.executable.resolve(),
        args.gene_results.resolve(),
        args.adaptive.resolve(),
        args.output.resolve(),
    )
