"""One fixed higher-precision resampling check; never iterate until significant."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from celltogwas import benjamini_hochberg, file_sha256, run_adaptive_gene_enrichment
from celltogwas.adaptive import _holm
from validation.run_adaptive_real import COVARIATES


def run(source: Path, output: Path):
    output.mkdir(parents=True, exist_ok=True)
    protocol = {
        "reason": "Resolve Monte Carlo uncertainty; one fixed check, no outcome-dependent stopping.",
        "n_permutations": 99999,
        "seed": 20261030,
        "inputs": {
            method: file_sha256(source / f"{method}.input.tsv.gz")
            for method in ["standard", "hmagma"]
        },
    }
    (output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n", encoding="utf-8")
    tables = []
    for method in ["standard", "hmagma"]:
        frame = pd.read_csv(source / f"{method}.input.tsv.gz", sep="\t")
        metadata = {"gene", "symbol", "chrom", "pos", "z", "gene_p", *COVARIATES}
        cells = [c for c in frame if c not in metadata]
        result = run_adaptive_gene_enrichment(
            frame,
            cell_types=cells,
            covariates=COVARIATES,
            n_permutations=99999,
            random_state=20261030,
        )
        tables.append(result.table.assign(method=method))
        print(f"{method}: fixed 99,999-resample check complete", flush=True)
    combined = pd.concat(tables, ignore_index=True)
    p = combined.p_adaptive.to_numpy()
    combined["p_holm_global48"] = _holm(p)
    combined["q_by_global48"] = np.minimum(
        benjamini_hochberg(p) * np.sum(1 / np.arange(1, len(p) + 1)), 1
    )
    combined.to_csv(output / "combined.tsv", sep="\t", index=False, lineterminator="\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.source, args.output)
