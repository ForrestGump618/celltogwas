from __future__ import annotations

import gzip

import pandas as pd

from celltogwas.cell2gwas import (
    celltype_genome_range,
    extract_gene_info,
    gwas_annotaion,
    gwas_annotation,
)


def test_legacy_round_trip(tmp_path):
    path = tmp_path / "genes.gtf.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(
            "1\tsrc\tgene\t100\t200\t.\t+\t.\t"
            'gene_id "G1"; gene_name "GeneA"; gene_type "protein_coding";\n'
        )

    gene_info = extract_gene_info(path)
    ranges = celltype_genome_range({"CellA": ["GeneA"]}, gene_info, 10)
    gwas = pd.DataFrame({"CHR": ["chr1", "1", "1", "2"], "POS": [100, 90, 89, 100]})
    annotated = gwas_annotation(gwas, ranges)

    assert ranges["CellA"]["chr1"] == [[89, 210]]
    assert annotated["CellA"].tolist() == ["1", "1", "0", "0"]
    assert gwas_annotaion(gwas, ranges).equals(annotated)
