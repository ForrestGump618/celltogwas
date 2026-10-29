# External validation

This directory contains the reproducible, real-data validation reported for
celltogwas 1.0. The pipeline uses:

- PsychENCODE brainSCOPE adult DLPFC donor-level pseudobulk profiles across 24
  cell classes.
- Bellenguez Alzheimer's disease GWAS `GCST90027158` (the accession metadata cites the 2020 preprint, DOI `10.1101/2020.10.01.20200659`; this is not asserted to be the final 2022 publication's full dataset).
- The official MAGMA 1000 Genomes European reference panel.
- GENCODE human release 19 gene models.
- Standard MAGMA annotation and H-MAGMA adult-brain Hi-C
  annotation, both analyzed with the same MAGMA v1.10 executable.

The repository does not redistribute GWAS, expression, reference-panel, or
MAGMA files. `data_manifest.json` records official URLs, byte sizes, and
SHA-256 checksums. MAGMA's own terms prohibit redistribution of its binary and
source files.

## 1. Download public inputs

```bash
python validation/download_data.py \
  --data-root /path/to/celltogwas-validation-data/v1
```

Add `--verify-only` to verify an existing tree without downloading. ZIP
archives are extracted below `magma/`. The expected layout is:

```text
v1/
  gwas/GCST90027158_buildGRCh38.tsv.gz
  gencode/gencode.v19.annotation.gtf.gz
  hmagma/Adult_brain.genes.annot
  hmagma/MAGMAdefault.genes.annot.gz
  magma/g1000_eur/g1000_eur.{bed,bim,fam}
  magma/magma_v1.10_win/magma.exe
  psychencode/*.expr.bed.gz
```

## 2. Run the analysis

From the repository root:

```bash
python -m validation.run_validation \
  --data-root /path/to/celltogwas-validation-data/v1 \
  --output-dir /path/to/results
```

The analysis parameters are fixed before inspecting results:

- Tau marker specificity with 1,000 markers per cell class, selected from genes reported in all 24 matrices and present in GENCODE v19.
- Minimum mean expression 0.1 and donor detection fraction 0.1.
- Same-genome SNP matching by exact rsID to the GRCh37 1000G EUR BIM; no
  liftover. This is a P-value analysis; effect signs and alleles are not harmonized.
- Per-SNP sample size `N = n_cases + n_controls`, retaining `N >= 50` as required by the MAGMA default.
- All matched SNPs are retained as the shift background. Primary celltogwas evidence: `P <= 5e-8`, weight floor `-log10(5e-8)`,
  10,000 within-chromosome circular shifts, seed 13.
- Benjamini-Hochberg correction across 24 cell classes.
- GENCODE v19 gene body plus bilateral 10 kb as the celltogwas genomic link.
- Sensitivity analyses at `P <= 1e-5` and excluding the MHC
  (`chr6:25,000,000-35,000,000`, GRCh37).

All analyses use autosomal SNPs. H-MAGMA is treated correctly as a Hi-C
SNP-to-gene annotation used inside the MAGMA statistical engine, not as an
independent statistical method.

## 3. Main outputs

- `validation_manifest.json`: runtime, parameters, input hashes, and summary.
- `harmonized/gwas_harmonization.json`: GWAS/reference matching counts.
- `expression/marker_summary.tsv`: marker counts and detection diagnostics.
- `celltogwas/*.tsv`: primary and sensitivity enrichment tables.
- `magma/*.genes.out.txt`: standard MAGMA and H-MAGMA gene-level results.
- `magma/*.gsa.out.txt`: cell-set association results.
- `comparison/*.json`: rank concordance, top-gene overlap, and method limits.

## Interpretation boundary

MAGMA performs an LD-aware gene and gene-set regression. celltogwas uses a
maximum interval annotation followed by circular spatial shifts and does not
model LD. Numeric values are therefore not expected to be identical. The
validation evaluates rank and membership consistency, biological plausibility,
and sensitivity to prespecified choices. It does not establish causality.


The supplied conventional MAGMA annotation and GENCODE v19 intervals are different gene models (the supplied conventional annotation uses strand-dependent upstream/downstream windows). This is a method comparison, not an identical-interval test. All-SNP evidence is also evaluated as a sensitivity analysis.

## Figures

Install the optional plotting dependency with `python -m pip install ".[validation]"`,
then render PNG and editable SVG figures:

```bash
python -m validation.plot_results --results /path/to/results --output /path/to/figures
```

The downloaded conventional MAGMA and adult-brain H-MAGMA annotation files use
older gene models than some brainSCOPE symbols. Marker selection uses genes
present in all 24 expression matrices and resolves symbols to GENCODE v19 IDs;
unmapped symbols, different donor coverage, and missing genes are reported.
The source values are log-normalized expression, not reconstructed raw counts.

## Official MAGMA batch execution

For larger input tables, `python -m validation.run_magma_batches` runs 1–8
independent batches using MAGMA's documented `--batch i n`, then its official
`--merge`. Every batch receives the same reference, P/N file, annotation, and
output prefix. All required batch outputs must exist before merging; the
merged `.genes.raw` retains the correlations required for competitive gene-set
analysis. Supply `--executable`, `--reference`, `--pvalues`, `--annotation`,
`--gene-sets`, `--prefix`, and optionally `--batches 4`. Run each annotation
model in a separate output prefix. The complete commands and input hashes are
written next to the merged output.

On Windows MAGMA v1.10 may emit `.genes.out.txt` while its merger looks for
`.genes.out`. The batch runner preserves originals and adds byte-identical
aliases, and calls the merger from the output directory with a relative prefix.
`--resume` reuses complete batches only when successful command logs and all input SHA-256 fingerprints match.
