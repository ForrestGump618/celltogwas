from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from validation.pipeline import (
    bootstrap_spearman_ci,
    compare_magma_gene_results,
    harmonize_gwas,
    read_magma_gene_results,
    read_magma_gene_set_results,
    read_psychencode_cell_type,
    spearman_rank,
    top_overlap,
)


def test_spearman_and_bootstrap() -> None:
    first = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    second = first**2

    assert spearman_rank(first, second) == pytest.approx(1.0)
    assert spearman_rank(first, second[::-1]) == pytest.approx(-1.0)
    lower, upper = bootstrap_spearman_ci(
        first,
        second,
        replicates=100,
        random_state=3,
    )
    assert lower <= 1.0 <= upper


def test_top_overlap_uses_score_order_and_jaccard() -> None:
    first = [("A", 5.0), ("B", 4.0), ("C", 3.0)]
    second = [("C", 9.0), ("B", 8.0), ("D", 7.0)]

    overlap, jaccard, first_top, second_top = top_overlap(first, second, k=2)

    assert overlap == 1
    assert jaccard == pytest.approx(1 / 3)
    assert first_top == ["A", "B"]
    assert second_top == ["C", "B"]


def test_read_psychencode_aggregates_donors_and_duplicate_symbols(tmp_path: Path) -> None:
    path = tmp_path / "Astro.expr.bed.gz"
    text = (
        "#chr\tstart\tend\tgene\tlength\tstrand\tD1\tD2\n"
        "chr1\t1\t2\tGENE1\t1\t+\t0\t2\n"
        "chr1\t3\t4\tGENE1\t1\t+\t4\t4\n"
        "chr1\t5\t6\tGENE2\t1\t+\t0\t0\n"
    )
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)

    mean, detection, donors, donor_ids = read_psychencode_cell_type(path)

    assert donor_ids == ["D1", "D2"]
    assert mean.loc["GENE1"] == pytest.approx(4.0)
    assert detection.loc["GENE1"] == pytest.approx(1.0)
    assert donors.loc["GENE1"] == 2
    assert mean.loc["GENE2"] == 0


def test_harmonize_gwas_uses_bim_positions_and_removes_ambiguous_rsids(
    tmp_path: Path,
) -> None:
    bim = tmp_path / "reference.bim"
    bim.write_text(
        "1\trs1\t0\t100\tA\tG\n1\trs2\t0\t200\tC\tT\n23\trs4\t0\t300\tA\tC\n",
        encoding="utf-8",
    )
    gwas = tmp_path / "gwas.tsv"
    gwas.write_text(
        "variant_id\tp_value\tchromosome\tbase_pair_location\tn_cases\tn_controls\n"
        "rs1\t0\t1\t900\t100\t200\n"
        "rs2\t0.1\t1\t901\t100\t200\n"
        "rs2\t0.2\t1\t902\t100\t200\n"
        "rs4\t0.3\t23\t903\t100\t200\n",
        encoding="utf-8",
    )
    output = tmp_path / "harmonized"

    manifest = harmonize_gwas(gwas, bim, output)
    result = pd.read_csv(output / "gwas_harmonized.tsv.gz", sep="\t")

    assert manifest["gwas"]["matched_rsids"] == 1
    assert manifest["gwas"]["duplicate_rsid_rows_removed"] == 2
    assert result.loc[0, "SNP"] == "rs1"
    assert result.loc[0, "POS"] == 100
    assert result.loc[0, "P"] == pytest.approx(1e-300)
    assert result.loc[0, "N"] == 300


def test_read_magma_gene_results(tmp_path: Path) -> None:
    path = tmp_path / "genes.out"
    path.write_text(
        "GENE CHR START STOP NSNPS NPARAM N ZSTAT P\n"
        "ENSG1 1 1 2 3 2 100 4.2 1e-5\n"
        "ENSG2 1 3 4 2 2 100 1.0 0.1\n",
        encoding="utf-8",
    )

    frame = read_magma_gene_results(path)

    assert list(frame.columns) == ["GENE", "ZSTAT", "P", "NSNPS"]
    assert frame.loc[0, "ZSTAT"] == pytest.approx(4.2)


def test_read_magma_gene_set_results_calculates_z_from_beta_and_se(
    tmp_path: Path,
) -> None:
    path = tmp_path / "sets.gsa.out.txt"
    path.write_text(
        "# TEST_DIRECTION = one-sided\n"
        "VARIABLE TYPE NGENES BETA BETA_STD SE P\n"
        "Astro SET 10 0.2 0.1 0.05 0.001\n",
        encoding="utf-8",
    )

    frame = read_magma_gene_set_results(path)

    assert frame.loc[0, "cell_type"] == "Astro"
    assert frame.loc[0, "magma_z"] == pytest.approx(4.0)
    assert frame.loc[0, "n_genes"] == 10


def test_compare_magma_genes_writes_known_ad_annotations(tmp_path: Path) -> None:
    gencode = tmp_path / "genes.gtf"
    gencode.write_text(
        '1\ts\tgene\t1\t2\t.\t+\t.\tgene_id "ENSG1.1"; gene_name "APOE";\n'
        '1\ts\tgene\t3\t4\t.\t+\t.\tgene_id "ENSG2.1"; gene_name "TREM2";\n',
        encoding="utf-8",
    )
    standard = tmp_path / "standard.genes.out"
    hmagma = tmp_path / "hmagma.genes.out"
    text = (
        "GENE CHR START STOP NSNPS NPARAM N ZSTAT P\n"
        "ENSG1 1 1 2 3 2 100 5.0 1e-7\n"
        "ENSG2 1 3 4 2 2 100 2.0 0.02\n"
    )
    standard.write_text(text + "ENSG3 1 5 6 3 2 100 6.0 1e-9\n", encoding="utf-8")
    hmagma.write_text(text + "ENSG4 1 7 8 3 2 100 6.0 1e-9\n", encoding="utf-8")

    result = compare_magma_gene_results(
        standard,
        hmagma,
        output_dir=tmp_path,
        gencode_path=gencode,
        bootstrap_replicates=20,
        random_state=7,
    )

    known = pd.read_csv(tmp_path / "known_ad_genes_magma.tsv", sep="\t")
    assert result["n_matched_genes_with_snps"] == 2
    assert result["zstat_spearman"] == pytest.approx(1.0)
    assert set(known["gene_symbol"]) == {"APOE", "TREM2"}
    summary = json.loads((tmp_path / "magma_gene_concordance.json").read_text(encoding="utf-8"))
    assert summary["top_gene_overlap"]["100"]["overlap"] == 2
    assert summary["top_gene_overlap"]["100"]["jaccard"] == pytest.approx(0.5)
    assert summary["top_gene_overlap_within_common_genes"]["100"]["jaccard"] == 1.0


def test_external_tool_error_is_detected_even_with_zero_exit_status(tmp_path: Path) -> None:
    import sys

    from validation.pipeline import _run_logged_command

    with pytest.raises(RuntimeError, match="Command failed"):
        _run_logged_command(
            [sys.executable, "-c", "print('ERROR - invalid annotation')"],
            cwd=tmp_path,
            log_path=tmp_path / "tool.log",
        )
    assert "ERROR - invalid annotation" in (tmp_path / "tool.log").read_text()


def test_zip_extraction_rejects_parent_traversal(tmp_path: Path) -> None:
    import zipfile

    from validation.download_data import safe_extract_zip

    archive = tmp_path / "bad.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("../outside.txt", "bad")
    with pytest.raises(ValueError, match="unsafe ZIP member"):
        safe_extract_zip(archive, tmp_path / "extract")
    assert not (tmp_path / "outside.txt").exists()


def test_windows_magma_batch_aliases_preserve_results_and_require_all_raw_files(
    tmp_path: Path,
) -> None:
    from validation.run_magma_batches import normalize_batch_filenames

    prefix = tmp_path / "analysis"
    for index in (1, 2):
        Path(f"{prefix}.batch{index}_2.genes.out.txt").write_text(
            f"GENE P\nGENE{index} 0.5\n", encoding="utf-8"
        )
        Path(f"{prefix}.batch{index}_2.genes.raw").write_bytes(b"test-raw")
    normalize_batch_filenames(prefix, 2)
    for index in (1, 2):
        assert (
            Path(f"{prefix}.batch{index}_2.genes.out").read_bytes()
            == Path(f"{prefix}.batch{index}_2.genes.out.txt").read_bytes()
        )
    Path(f"{prefix}.batch2_2.genes.raw").unlink()
    with pytest.raises(FileNotFoundError, match="Incomplete MAGMA batch 2/2"):
        normalize_batch_filenames(prefix, 2)


def test_batch_resume_refuses_changed_input_fingerprints(tmp_path: Path) -> None:
    from validation.run_magma_batches import verify_batch_inputs

    prefix = tmp_path / "analysis"
    hashes = {"pvalues": "original-hash", "annotation": "annotation-hash"}
    verify_batch_inputs(prefix, hashes, resume=False)
    verify_batch_inputs(prefix, hashes, resume=True)
    with pytest.raises(ValueError, match="changed batch input fingerprints"):
        verify_batch_inputs(prefix, {**hashes, "pvalues": "modified-hash"}, resume=True)
