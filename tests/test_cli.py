from __future__ import annotations

import json

import pandas as pd

from celltogwas.cli import main


def write_inputs(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        "1\tsrc\tgene\t100\t200\t.\t+\t.\t"
        'gene_id "G1"; gene_name "GeneA"; gene_type "protein_coding";\n'
        "1\tsrc\tgene\t1000\t1100\t.\t-\t.\t"
        'gene_id "G2"; gene_name "GeneB"; gene_type "protein_coding";\n',
        encoding="utf-8",
    )
    gene_sets = tmp_path / "sets.gmt"
    gene_sets.write_text("CellA\tdesc\tGeneA\nCellB\tdesc\tGeneB\n", encoding="utf-8")
    gwas = tmp_path / "gwas.csv"
    pd.DataFrame(
        {
            "CHR": [1, 1, 1, 1, 1, 1, 1, 1],
            "POS": [120, 140, 160, 180, 250, 500, 1020, 1080],
            "P": [1e-8, 1e-7, 1e-8, 1e-6, 0.2, 0.5, 0.3, 0.4],
        }
    ).to_csv(gwas, index=False)
    return gtf, gene_sets, gwas


def test_cli_annotate_and_enrich(tmp_path):
    gtf, gene_sets, gwas = write_inputs(tmp_path)
    annotated_path = tmp_path / "annotated.csv"
    annotation_report = tmp_path / "annotation.json"
    exit_code = main(
        [
            "annotate",
            "--gtf",
            str(gtf),
            "--gene-sets",
            str(gene_sets),
            "--gwas",
            str(gwas),
            "--output",
            str(annotated_path),
            "--window",
            "50",
            "--kernel",
            "uniform",
            "--report",
            str(annotation_report),
            "--quiet",
        ]
    )
    assert exit_code == 0
    annotated = pd.read_csv(annotated_path)
    assert {"CellA", "CellB"}.issubset(annotated.columns)
    assert json.loads(annotation_report.read_text(encoding="utf-8"))["version"]

    enrichment_path = tmp_path / "enrichment.csv"
    enrichment_report = tmp_path / "enrichment.json"
    exit_code = main(
        [
            "enrich",
            "--gtf",
            str(gtf),
            "--gene-sets",
            str(gene_sets),
            "--gwas",
            str(gwas),
            "--output",
            str(enrichment_path),
            "--window",
            "50",
            "--kernel",
            "uniform",
            "--permutations",
            "100",
            "--report",
            str(enrichment_report),
            "--quiet",
        ]
    )
    assert exit_code == 0
    enrichment = pd.read_csv(enrichment_path)
    assert {"cell_type", "empirical_p", "fdr_q"}.issubset(enrichment.columns)
    assert "enrichment" in json.loads(enrichment_report.read_text(encoding="utf-8"))


def test_cli_specificity(tmp_path):
    expression = tmp_path / "expression.csv"
    pd.DataFrame(
        {
            "gene": ["CD3D", "MS4A1", "LYZ"],
            "T_cell": [10.0, 0.0, 0.1],
            "B_cell": [0.0, 9.0, 0.1],
            "Myeloid": [0.1, 0.1, 8.0],
        }
    ).to_csv(expression, index=False)
    output = tmp_path / "markers.csv"
    gmt = tmp_path / "markers.gmt"
    report = tmp_path / "specificity.json"

    exit_code = main(
        [
            "specificity",
            "--expression",
            str(expression),
            "--output",
            str(output),
            "--gmt-output",
            str(gmt),
            "--report",
            str(report),
            "--quiet",
        ]
    )

    assert exit_code == 0
    markers = pd.read_csv(output)
    assert {"cell_type", "gene", "weight", "specificity"}.issubset(markers.columns)
    assert gmt.exists()
    assert json.loads(report.read_text(encoding="utf-8"))["method"] == "tau"


def test_cli_validate_writes_input_hashes(tmp_path):
    gtf, gene_sets, gwas = write_inputs(tmp_path)
    report = tmp_path / "validation.json"

    exit_code = main(
        [
            "validate",
            "--gtf",
            str(gtf),
            "--gene-sets",
            str(gene_sets),
            "--gwas",
            str(gwas),
            "--report",
            str(report),
            "--quiet",
        ]
    )

    assert exit_code == 0
    payload = json.loads(report.read_text(encoding="utf-8"))
    assert payload["valid"] is True
    assert len(payload["inputs"]) == 3
    assert payload["gwas_validation"]["n_rows_checked"] == 8
