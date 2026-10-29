"""Legacy import path.

Existing code can continue to use::

    from celltogwas.cell2gwas import gwas_annotation

New code should prefer :class:`celltogwas.GenomeLinker`.
"""

from .legacy import (
    celltype_genome_range,
    check_value_in_intervals,
    extract_gene_info,
    gwas_annotaion,
    gwas_annotation,
    merge_intervals,
    run,
)

__all__ = [
    "celltype_genome_range",
    "check_value_in_intervals",
    "extract_gene_info",
    "gwas_annotaion",
    "gwas_annotation",
    "merge_intervals",
    "run",
]
