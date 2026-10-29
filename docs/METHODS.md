# Methods

This document defines the algorithms used by celltogwas 1.0 and the
interpretation boundary of each result.

## 1. Coordinates and gene records

GTF/GFF3 records are one-based and end-inclusive. They are converted once to
the 0-based, half-open convention `[start, end)`:

```text
start_0 = start_1 - 1
end_0 = end_1
```

Circulating chromosome labels are normalized by removing a leading `chr` and
mapping mitochondrial aliases to `MT`. The GWAS position is converted as
`position_0 = position_1 - 1`.

By default, only `gene` feature rows are retained. This is intentional:
counting transcript or exon rows as independent genes duplicates intervals and
distorts both domain size and cell-type coverage.

## 2. Weighted genomic domains

For each cell type `c`, a gene set is a weighted mapping
`G_c = {(g, w_g)}`. Each resolved gene record contributes:

- gene body interval with weight `w_g * gene_body_weight`
- directional TSS promoter interval with weight
  `w_g * promoter_weight`
- optional cis distance bands with weight `w_g * K(d)`
- optional regulatory intervals with their own weights

The cell-type link score is the maximum active weight at a SNP:

```text
S_c(x) = max_{(I,w) in intervals_c, x in I} w
```

Using the maximum rather than a sum avoids inflating a score when several
transcripts or overlapping annotations describe the same locus.

### Interval merge algorithm

The implementation sweeps over interval starts and ends while maintaining a
maximum heap with lazy deletion over active weights. It emits non-overlapping
segments whose weight is the maximum active source weight. Lookup is then
`searchsorted(start, x) - 1`, followed by an end-boundary check, giving
`O(log n)` query time after an `O(n log n)` build.

### Exponential kernel

The default soft kernel is:

```text
K(d) = exp(-d / gamma),  gamma = window / 3 unless specified
```

It is discretized into `kernel_levels` distance bands. For a band
`[a, b]`, the representative distance is `(a + b) / 2`. This is an
explicit approximation, documented in diagnostics by the kernel and band
count. Use `kernel="uniform"` when a hard window is required.

## 3. Expression specificity

Expression rows are genes and columns are cell types. Values must be
non-negative. The default Tau-style specificity index is:

```text
tau_g = sum_c (1 - x_gc / max_c(x_g)) / (C - 1)
```

For one cell type, `tau_g = 1`. Housekeeping profiles have `tau` close to zero;
strictly cell-type-specific profiles have `tau` close to one.

The relative expression and final marker weight are:

```text
relative_gc = x_gc / max_c(x_g)
marker_gc = (tau_g * relative_gc) ** power
```

The entropy alternative uses:

```text
p_gc = (x_gc + epsilon) / sum_c (x_gc + epsilon)
specificity_g = 1 - H(p_g) / log(C)
```

The robust z-score alternative centers each gene by its cross-cell-type
median:

```text
z_gc = (x_gc - median_c(x_gc)) / s_g
s_g = 1.4826 * MAD_c(x_gc)
marker_gc = 1 - exp(-max(z_gc, 0))
```

When `MAD = 0`, which is common when most cell types have zero expression for
a gene, `s_g` falls back to the mean absolute deviation about the median. A
constant profile therefore receives zero marker weight instead of an
undefined or arbitrary score.

Gene-level specificity alone is insufficient: a gene with a highly specific
but very low expression profile can be biologically uninformative. Relative expression down-weights cell types below that gene's maximum; it is invariant to scaling the whole gene profile. Use `min_expression` for absolute-expression filtering.

## 4. Spatial enrichment

Let `w_i` be a non-negative GWAS evidence weight and `S_ic` the link score of
SNP `i` to cell type `c`. The observed statistic is:

```text
T_c = (sum_i w_i S_ic) / (sum_i w_i)
```

For `-log10(P)` scoring, a SNP with `P=1` contributes no evidence. The
`weight_floor` option subtracts a constant before clipping at zero, which can
be used to focus on genome-wide-significant-like signals.

Evidence modes are explicit:

- `pvalue`: `w_i = max(-log10(P_i) - weight_floor, 0)`
- `abs`: `w_i = max(|effect_i| - weight_floor, 0)`
- `zscore`: robustly standardizes `|effect_i|` by its median and scaled MAD,
  clamps values below the median to zero, then applies `weight_floor`
- `uniform`: every eligible SNP receives weight one

The fitted center and scale for `zscore` mode are recorded in diagnostics.

### Circular null

Within each chromosome, the annotation vector `S_c` is circularly shifted
relative to `w`:

```text
T_c^(b) = (sum_i w_i S_c[(i + k_b) mod n]) / (sum_i w_i)
```

Each permutation samples one offset per chromosome. Offsets closer than
`min_shift_fraction * n` to zero are excluded where possible to reduce
contamination by the observed alignment.

The shift is applied to the ordered matched-variant vector, not to a fixed
base-pair grid. It preserves the ordered score sequence, but does not preserve physical base-pair distances across shifted positions. It is not a fixed-distance genomic rotation.

The full set of circular weighted dot products is computed with:

```text
cross_correlation = IFFT(FFT(S_c) * conjugate(FFT(w)))
```

This replaces a naive `O(n * n_permutations)` loop with an FFT transform plus
one lookup per sampled offset.

### Empirical p-value and FDR

The one-sided empirical p-value is:

```text
p_c = (1 + #{b: T_c^(b) >= T_c}) / (n_permutations + 1)
```

Benjamini-Hochberg adjusted values are computed across cell types:

```text
q_(i) = min_{j >= i} p_(j) * m / j
```

The test is exploratory. It accounts for local ordering and clustering through
within-chromosome shifts, but it does not explicitly estimate linkage
disequilibrium or heritability. Use S-LDSC, MAGMA/H-MAGMA, or a suitable
fine-mapping framework when those are the inferential targets.

## 5. Quality controls

The package reports:

- parsed gene count and source SHA-256
- requested versus matched gene symbols per cell type
- unmatched genes
- gene symbols resolving to multiple annotation records
- source and merged interval counts
- covered bases and weighted covered bases
- number of eligible GWAS rows and chromosomes
- skipped small chromosomes
- permutation count, random seed, evidence-weight mode, and multiple-testing
  method

## References

- de Leeuw CA, Mooij JM, Heskes T, Posthuma D. MAGMA: generalized gene-set
  analysis of GWAS data. PLoS Computational Biology. 2015.
  https://doi.org/10.1371/journal.pcbi.1004219
- Sey NYA et al. A computational tool for improved prediction of brain-disorder
  risk genes by incorporating brain chromatin interaction profiles. Nature
  Neuroscience. 2020. https://doi.org/10.1038/s41593-020-0603-0
- Jagadeesh KA et al. Identifying disease-critical cell types and cellular
  processes by integrating single-cell RNA sequencing and human genetics.
  Nature Genetics. 2022. https://doi.org/10.1038/s41588-022-01187-9
- Finucane HK et al. Heritability enrichment of specifically expressed genes
  identifies disease-relevant tissues and cell types. Nature Genetics. 2018.
  https://doi.org/10.1038/s41588-018-0081-4
- Timshel PN, Thompson JJ, Pers TH. Genetic mapping of etiologic brain cell
  types for obesity. eLife. 2020. https://doi.org/10.7554/eLife.55851
- Yanai I et al. Genome-wide midrange transcription profiles reveal expression
  level relationships in human tissue specification. Bioinformatics. 2005.
  https://doi.org/10.1093/bioinformatics/bti042
- McLean CY et al. GREAT improves functional interpretation of cis-regulatory
  regions. Nature Biotechnology. 2010. https://doi.org/10.1038/nbt.1630

`p_upper` assigns zero evidence above the threshold while retaining all valid SNPs as the circular-shift background. `run_enrichment_by_chromosome` combines independent per-chromosome shifts with bounded memory and child seeds derived from the root seed.
