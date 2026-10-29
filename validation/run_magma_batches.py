"""Run the official MAGMA batch/merge workflow with identical per-batch inputs."""

from __future__ import annotations

import argparse
import json
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from celltogwas.io import write_json
from validation.pipeline import _run_logged_command, sha256_file


def normalize_batch_filenames(prefix: Path, batches: int) -> None:
    """Keep originals while providing names expected by the Windows merger."""

    for index in range(1, batches + 1):
        stem = f"{prefix}.batch{index}_{batches}"
        plain, windows = Path(stem + ".genes.out"), Path(stem + ".genes.out.txt")
        if windows.exists():
            shutil.copyfile(windows, plain)
        if not plain.is_file() or not Path(stem + ".genes.raw").is_file():
            raise FileNotFoundError(f"Incomplete MAGMA batch {index}/{batches}: {stem}")


def batch_input_hashes(
    executable: Path, reference: Path, pvalues: Path, annotation: Path, gene_sets: Path
) -> dict[str, str]:
    files = {
        "executable": executable,
        "pvalues": pvalues,
        "annotation": annotation,
        "gene_sets": gene_sets,
    }
    for suffix in ("bed", "bim", "fam", "synonyms"):
        path = Path(f"{reference}.{suffix}")
        if suffix != "synonyms" or path.exists():
            files[f"reference_{suffix}"] = path
    return {name: sha256_file(path) for name, path in files.items()}


def verify_batch_inputs(prefix: Path, hashes: dict[str, str], *, resume: bool) -> None:
    receipt = Path(f"{prefix}.inputs.json")
    if resume:
        if not receipt.exists() or json.loads(receipt.read_text(encoding="utf-8")) != hashes:
            raise ValueError("Missing or changed batch input fingerprints; rerun without --resume.")
    else:
        write_json(hashes, receipt)


def run_batches(
    executable: Path,
    reference: Path,
    pvalues: Path,
    annotation: Path,
    gene_sets: Path,
    prefix: Path,
    batches: int = 4,
    resume: bool = False,
) -> None:
    if not 1 <= batches <= 8:
        raise ValueError("Use between 1 and 8 batches to limit memory usage.")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    hashes = batch_input_hashes(executable, reference, pvalues, annotation, gene_sets)
    verify_batch_inputs(prefix, hashes, resume=resume)
    common = [
        str(executable),
        "--bfile",
        str(reference),
        "--pval",
        str(pvalues),
        "use=rsid,P",
        "ncol=N",
        "--gene-annot",
        str(annotation),
        "--out",
        str(prefix),
    ]

    def run_one(index):
        command = common + ["--batch", str(index), str(batches)]
        log_path = Path(f"{prefix}.batch{index}.console.log")
        if resume and log_path.exists():
            text = log_path.read_text(encoding="utf-8")
            previous = json.loads(text.splitlines()[0].partition("=")[2])
            stem = f"{prefix}.batch{index}_{batches}"
            if (
                previous == command
                and "returncode=0\n" in text
                and Path(stem + ".genes.raw").exists()
                and (Path(stem + ".genes.out.txt").exists() or Path(stem + ".genes.out").exists())
            ):
                return {
                    "command": command,
                    "returncode": 0,
                    "log": str(log_path),
                    "resumed": True,
                }
        return _run_logged_command(
            command,
            cwd=prefix.parent,
            log_path=log_path,
        )

    with ThreadPoolExecutor(max_workers=batches) as executor:
        results = list(executor.map(run_one, range(1, batches + 1)))
    normalize_batch_filenames(prefix, batches)
    merged = _run_logged_command(
        [str(executable), "--merge", prefix.name, str(batches), "--out", str(prefix)],
        cwd=prefix.parent,
        log_path=Path(f"{prefix}.merge.console.log"),
    )
    for suffix in (".genes.out.txt", ".genes.raw"):
        if not Path(f"{prefix}{suffix}").is_file():
            raise RuntimeError(f"MAGMA merge did not produce {prefix}{suffix}")
    sets_prefix = Path(f"{prefix}_geneset")
    sets = _run_logged_command(
        [
            str(executable),
            "--gene-results",
            f"{prefix}.genes.raw",
            "--set-annot",
            str(gene_sets),
            "--out",
            str(sets_prefix),
        ],
        cwd=prefix.parent,
        log_path=Path(f"{prefix}.geneset.console.log"),
    )
    write_json(
        {
            "engine": "MAGMA v1.10",
            "inputs_sha256": hashes,
            "batch_count": batches,
            "common_command": common,
            "batches": results,
            "merge": merged,
            "gene_set_run": sets,
            "annotation_sha256": sha256_file(annotation),
            "pvalue_file_sha256": sha256_file(pvalues),
            "gene_set_file_sha256": sha256_file(gene_sets),
            "gene_results": f"{prefix}.genes.out.txt",
            "gene_raw": f"{prefix}.genes.raw",
            "gene_set_results": f"{sets_prefix}.gsa.out.txt",
        },
        Path(f"{prefix}.manifest.json"),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("executable", "reference", "pvalues", "annotation", "gene-sets", "prefix"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--batches", type=int, default=4)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse completed batches with identical command logs and complete outputs.",
    )
    args = parser.parse_args()
    run_batches(
        args.executable,
        args.reference,
        args.pvalues,
        args.annotation,
        args.gene_sets,
        args.prefix,
        args.batches,
        args.resume,
    )
