"""celltogwas: regulatory-aware cell-type to GWAS mapping.

Version 1.0 keeps the 0.1 helper functions available while introducing typed
annotation, weighted genomic domains, expression-specificity weights and
spatial permutation enrichment.
"""

from ._version import __version__
from .adaptive import AdaptiveResult, calibrate_adaptive, run_adaptive_gene_enrichment
from .enrichment import (
    EnrichmentResult,
    benjamini_hochberg,
    circular_correlation_matrix,
    circular_permutation_enrichment,
    run_enrichment,
    run_enrichment_by_chromosome,
)
from .exceptions import (
    AnnotationError,
    CellToGWASError,
    EnrichmentError,
    InputFormatError,
    LinkageError,
)
from .interval import (
    WeightedIntervalIndex,
    check_value_in_intervals,
    merge_intervals,
    normalize_chromosome,
)
from .io import (
    file_sha256,
    load_expression_matrix,
    load_gene_annotation,
    load_gene_sets,
    load_regulatory_bed,
    read_table,
    write_gene_sets,
    write_json,
    write_table,
)
from .legacy import (
    celltype_genome_range,
    extract_gene_info,
    gwas_annotaion,
    gwas_annotation,
)
from .linker import GenomeLinker
from .models import GeneAnnotation, LinkConfig
from .specificity import (
    compute_specificity,
    specificity_matrix,
    specificity_summary,
    specificity_to_gene_sets,
)

__all__ = [
    "AdaptiveResult",
    "AnnotationError",
    "CellToGWASError",
    "EnrichmentError",
    "EnrichmentResult",
    "GeneAnnotation",
    "GenomeLinker",
    "InputFormatError",
    "LinkConfig",
    "LinkageError",
    "WeightedIntervalIndex",
    "__version__",
    "benjamini_hochberg",
    "calibrate_adaptive",
    "celltype_genome_range",
    "check_value_in_intervals",
    "circular_correlation_matrix",
    "circular_permutation_enrichment",
    "compute_specificity",
    "extract_gene_info",
    "file_sha256",
    "gwas_annotaion",
    "gwas_annotation",
    "load_expression_matrix",
    "load_gene_annotation",
    "load_gene_sets",
    "load_regulatory_bed",
    "merge_intervals",
    "normalize_chromosome",
    "read_table",
    "run_enrichment",
    "run_adaptive_gene_enrichment",
    "run_enrichment_by_chromosome",
    "specificity_matrix",
    "specificity_summary",
    "specificity_to_gene_sets",
    "write_gene_sets",
    "write_json",
    "write_table",
]
