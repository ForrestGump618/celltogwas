from __future__ import annotations

import gzip
import json

import pandas as pd

from celltogwas.io import (
    load_expression_matrix,
    load_gene_annotation,
    load_gene_sets,
    load_regulatory_bed,
    read_table,
    write_gene_sets,
)

GTF_TEXT = """\
##description: test
1\tsrc\tgene\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "GeneA"; gene_type "protein_coding";
1\tsrc\ttranscript\t100\t150\t.\t+\t.\tgene_id "G1"; gene_name "GeneA"; transcript_id "T1";
1\tsrc\tgene\t400\t500\t.\t-\t.\tgene_id "G2"; gene_name "GeneB"; gene_type "lncRNA";
"""


def test_load_gene_annotation_keeps_gene_records_and_converts_coordinates(tmp_path):
    path = tmp_path / "genes.gtf.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(GTF_TEXT)

    annotation = load_gene_annotation(path)

    assert annotation.n_genes == 2
    frame = annotation.as_legacy_frame()
    gene_a = frame.loc[frame["GENE"] == "GeneA"].iloc[0]
    assert gene_a["CHR"] == "1"
    assert gene_a["START"] == 100
    assert gene_a["END"] == 200
    assert annotation.source_sha256


def test_load_gene_annotation_supports_gff3(tmp_path):
    path = tmp_path / "genes.gff3"
    path.write_text(
        "##gff-version 3\n"
        "chr2\tsrc\tgene\t50\t90\t.\t-\t.\t"
        "ID=ENSG1;Name=Gene%20One;biotype=protein_coding\n",
        encoding="utf-8",
    )

    annotation = load_gene_annotation(path)

    assert annotation.genes.loc[0, "gene_id"] == "ENSG1"
    assert annotation.genes.loc[0, "gene_name"] == "Gene One"
    assert annotation.genes.loc[0, "chrom"] == "2"
    assert annotation.genes.loc[0, "start"] == 49
    assert annotation.genes.loc[0, "end"] == 90


def test_load_gene_sets_from_gmt_json_and_long_table(tmp_path):
    gmt = tmp_path / "sets.gmt"
    gmt.write_text("T_cell\tdesc\tCD3D\tCD3E:0.8\n", encoding="utf-8")
    loaded_gmt = load_gene_sets(gmt)
    assert loaded_gmt["T_cell"] == {"CD3D": 1.0, "CD3E": 0.8}

    json_path = tmp_path / "sets.json"
    json_path.write_text(
        json.dumps({"B_cell": {"MS4A1": 0.9, "CD79A": 0.7}}),
        encoding="utf-8",
    )
    assert load_gene_sets(json_path)["B_cell"]["MS4A1"] == 0.9

    table = tmp_path / "sets.csv"
    pd.DataFrame(
        {
            "cell_type": ["Myeloid", "Myeloid", "T_cell"],
            "gene": ["LYZ", "FCGR3A", "CD3D"],
            "weight": [0.9, 0.6, 0.8],
        }
    ).to_csv(table, index=False)
    loaded_table = load_gene_sets(table)
    assert loaded_table["Myeloid"] == {"LYZ": 0.9, "FCGR3A": 0.6}


def test_load_gene_sets_wide_table_keeps_numeric_gene_ids(tmp_path):
    table = tmp_path / "sets.csv"
    pd.DataFrame(
        {
            "gene": [101, 102],
            "T_cell": [1.0, 0.0],
            "B_cell": [0.0, 0.8],
        }
    ).to_csv(table, index=False)

    loaded = load_gene_sets(table)

    assert loaded == {"T_cell": {"101": 1.0}, "B_cell": {"102": 0.8}}


def test_write_gene_sets_round_trips_through_gmt_reader(tmp_path):
    gene_sets = {
        "T_cell": {"CD3D": 1.0, "CD3E": 0.8},
        "Myeloid": {"LYZ": 0.75},
    }
    path = tmp_path / "sets.gmt.gz"

    write_gene_sets(gene_sets, path)

    assert load_gene_sets(path) == gene_sets


def test_read_table_sniffs_tab_delimited_txt(tmp_path):
    path = tmp_path / "variants.txt"
    path.write_text("CHR\tPOS\tP\n1\t100\t0.01\n", encoding="utf-8")

    frame = read_table(path)

    assert list(frame.columns) == ["CHR", "POS", "P"]
    assert frame.loc[0, "POS"] == 100


def test_load_expression_matrix_keeps_numeric_gene_ids(tmp_path):
    path = tmp_path / "expression.csv"
    path.write_text("gene,T_cell,B_cell\n101,5,0\n102,0,4\n", encoding="utf-8")

    frame = load_expression_matrix(path)

    assert list(frame.index) == [101, 102]
    assert list(frame.columns) == ["T_cell", "B_cell"]


def test_load_regulatory_bed4_with_explicit_cell_type(tmp_path):
    bed = tmp_path / "regulatory.bed"
    bed.write_text("chr1\t10\t20\t0.7\n", encoding="utf-8")

    frame = load_regulatory_bed(bed, cell_type="CellA")

    assert frame.loc[0, "cell_type"] == "CellA"
    assert frame.loc[0, "weight"] == 0.7


def test_load_regulatory_bed_header_aliases(tmp_path):
    bed = tmp_path / "regulatory.bed"
    bed.write_text(
        "chr\tchromStart\tchromEnd\tcelltype\tscore\nchr1\t10\t20\tMicroglia\t0.8\n",
        encoding="utf-8",
    )

    frame = load_regulatory_bed(bed)

    assert frame.loc[0, "chrom"] == "1"
    assert frame.loc[0, "cell_type"] == "Microglia"
    assert frame.loc[0, "weight"] == 0.8
