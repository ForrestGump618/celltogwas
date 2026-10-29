from __future__ import annotations

import pandas as pd
import pytest

from celltogwas.exceptions import LinkageError
from celltogwas.linker import GenomeLinker
from celltogwas.models import GeneAnnotation, LinkConfig


def annotation_fixture() -> GeneAnnotation:
    genes = pd.DataFrame(
        {
            "gene_id": ["G1", "G2"],
            "gene_name": ["GeneA", "GeneB"],
            "chrom": ["1", "1"],
            "start": [100, 1000],
            "end": [200, 1100],
            "strand": ["+", "-"],
            "biotype": ["protein_coding", "protein_coding"],
        }
    )
    return GeneAnnotation(genes=genes)


def config_fixture() -> LinkConfig:
    return LinkConfig(
        window=100,
        promoter_window=10,
        promoter_weight=1.0,
        gene_body_weight=0.8,
        kernel="uniform",
    )


def test_annotate_weight_and_binary_modes():
    linker = GenomeLinker.from_annotation(
        annotation_fixture(),
        {"CellA": {"GeneA": 1.0}, "CellB": ["GeneB"]},
        config=config_fixture(),
    )
    gwas = pd.DataFrame(
        {
            "CHR": ["chr1", "1", "1", "2"],
            "BP": [150, 50, 350, 150],
        }
    )

    weighted = linker.annotate(gwas, mode="weight")
    binary = linker.annotate(gwas, mode="binary", threshold=0.5)

    assert weighted.loc[0, "CellA"] == 0.8
    assert weighted.loc[1, "CellA"] == 1.0
    assert weighted.loc[2, "CellA"] == 0.0
    assert weighted.loc[3, "CellA"] == 0.0
    assert binary.loc[0, "CellA"] == 1
    assert binary.loc[2, "CellA"] == 0


def test_regulatory_intervals_and_diagnostics():
    regulatory = pd.DataFrame(
        {
            "chrom": ["1"],
            "start": [500],
            "end": [600],
            "cell_type": ["RegulatoryOnly"],
            "weight": [0.7],
        }
    )
    linker = GenomeLinker.from_annotation(
        annotation_fixture(),
        {"CellA": {"GeneA": 1.0, "MISSING": 0.5}},
        config=config_fixture(),
        regulatory=regulatory,
    )

    annotated = linker.annotate(pd.DataFrame({"CHR": [1], "POS": [550]}))

    assert annotated.loc[0, "RegulatoryOnly"] == 0.7
    diagnostic = linker.diagnostics.set_index("cell_type")
    assert diagnostic.loc["CellA", "matched_genes"] == 1
    assert diagnostic.loc["CellA", "missing_genes"] == ["MISSING"]


def test_gene_id_and_symbol_resolve_to_one_domain():
    linker = GenomeLinker.from_annotation(
        annotation_fixture(),
        {"CellA": {"GeneA": 0.4, "G1": 0.9}},
        config=config_fixture(),
    )
    annotated = linker.annotate(pd.DataFrame({"CHR": ["1"], "POS": [150]}))
    assert annotated.loc[0, "CellA"] == pytest.approx(0.72)


def test_multilocus_gene_records_are_all_retained():
    genes = pd.DataFrame(
        {
            "gene_id": ["G1", "G1"],
            "gene_name": ["MultiLocus", "MultiLocus"],
            "chrom": ["1", "2"],
            "start": [100, 100],
            "end": [200, 200],
            "strand": ["+", "+"],
            "biotype": ["protein_coding", "protein_coding"],
        }
    )
    linker = GenomeLinker.from_annotation(
        GeneAnnotation(genes=genes),
        {"CellA": {"MultiLocus": 1.0}},
        config=LinkConfig(window=0),
    )

    annotated = linker.annotate(pd.DataFrame({"CHR": ["1", "2"], "POS": [150, 150]}))

    assert annotated["CellA"].tolist() == [1.0, 1.0]
    assert linker.diagnostics.loc[0, "matched_gene_records"] == 2
    assert linker.diagnostics.loc[0, "ambiguous_genes"] == ["MultiLocus"]


def test_link_config_rejects_non_integer_and_non_finite_values():
    with pytest.raises(LinkageError, match="kernel_levels must be an integer"):
        LinkConfig(kernel_levels=2.5)
    with pytest.raises(LinkageError, match="max_weight must be finite"):
        LinkConfig(max_weight=float("nan"))
