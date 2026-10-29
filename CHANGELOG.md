# Changelog

## 1.0.0

### Added

- Reproducible real brainSCOPE pseudobulk / AD GWAS validation with 24 cell
  classes and 14,161,368 harmonized SNPs.
- Actual MAGMA v1.10 and adult-brain H-MAGMA gene/gene-set comparisons,
  native/common gene-universe overlap, rank uncertainty, and per-method FDR.
- Official MAGMA batch/merge execution, immutable-input resume receipts,
  Windows filename/prefix compatibility, scientific plots and validation reports.
- Chromosome-streamed enrichment and a full-background P-threshold policy.

- Typed `GeneAnnotation`, `LinkConfig`, and `GenomeLinker` APIs.
- Gene-level GTF/GFF3 parsing with gzip support.
- Weighted promoter, gene-body, distance-kernel, and regulatory intervals.
- Exponentially decayed link scores with configurable distance bands.
- Expression specificity methods: Tau, normalized entropy, and robust z-score.
- Median/MAD standardization with a zero-MAD fallback for robust expression
  and effect weights.
- FFT-accelerated within-chromosome circular permutation enrichment.
- Benjamini-Hochberg FDR correction across cell types.
- CSV/TSV/JSON/GMT input and output helpers.
- CLI commands: `annotate`, `enrich`, `specificity`, and `validate`.
- Reproducible JSON reports with input SHA-256 hashes and diagnostics.
- Test suite, lint configuration, CI workflow, and migration documentation.

### Fixed

- GTF transcript rows no longer duplicate gene intervals.
- Genes with multiple annotation records retain every locus, with ambiguity
  reported in the build manifest.
- One-based GTF/GWAS boundaries are converted consistently.
- Chromosome aliases are normalized without character-stripping bugs.
- Missing chromosomes yield zero annotation rather than `KeyError`.
- Binary legacy output always uses string labels.
- Dataset-wide text files are delimiter-sniffed, and weighted GMT output now
  round-trips through the GMT reader.
- Regulatory BED headers accept common chromosome, coordinate, cell-type and
  score aliases.
- Specificity threshold semantics are consistent between marker computation
  and gene-set export.
- Expression filtering keeps ``max_expression`` aligned with retained genes.
- Wide gene-set tables keep numeric gene identifiers instead of treating them
  as a cell-type column.
- Expression tables keep numeric gene identifiers in the matrix index.
- Non-finite specificity weights are rejected before GMT export.
- Package metadata now follows the repository's actual GPL-3.0-only license.

### Compatibility

- Historical `celltogwas.cell2gwas` imports remain available.
