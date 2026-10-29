# Adaptive gene-program enrichment (experimental 1.1)

## Scientific target

Increase power across unknown sparse/dense association architectures without
optimizing the target trait's significant-hit count. The unit of inference is
a **cell program**, not a causal gene or a particular variant. This method
tests competitive enrichment on a fully measured background. Do not filter the
input down to GWAS hits or markers; zero-weight genes belong in the background.

Evolutionary/genetic search against the validation GWAS P values is deliberately
not used. Any future learned weights must be trained on disjoint data, then
frozen and evaluated on independent traits/cohorts. An adaptive test corrects
for a fixed kernel family; it does not correct for an analyst repeatedly changing
the family, dataset, covariates or phenotype after inspecting results.

## Fixed transforms

For cell c, normalize non-negative gene weights by their genome-wide maximum.
Use weight powers 1, 2, 4 and a binary top-10% program (including cutoff ties,
excluding zero weights). These span diffuse and concentrated expression
programs without a hard top-1,000-only truncation. No expression ranking uses GWAS.

For each chromosome, construct the nuisance design X from an intercept and
prespecified covariates. Center/scale columns and compute its numerical column
space Q using rank-revealing SVD. Redundant columns are dropped, not inverted.
Let M = I - Q Q'. First residualize the gene association-strength Z on X, and
clip M Z to [-10, 10]. Define evidence transforms

```
e_dense = clip(M Z, -10, 10)
e_2     = max(clip(M Z, -10, 10) - 2, 0)
e_4     = max(clip(M Z, -10, 10) - 4, 0)
```

The sparse cutoffs concern **residual association strength**, not SNP P values
or signed allelic effects. For each weight transform w and evidence e, use

```
T[k,c](s) = sum_over_chromosomes (M e)' roll(M w_c, -s_chromosome)
```

This gives 4 x 3 = 12 kernels per cell. All use the SAME independent uniform
chromosome offsets in each resample, including shift zero. The observed row uses
zero on every chromosome. Shifts operate in sorted gene order, not uniform
base-pair distance. FFT computes all cross-correlations for a chromosome, then
samples their offsets. Signed residual correlations MUST NOT be clipped at zero.

## Selection correction

Let B be the number of randomizations and b=0 the observed row. For every row,
kernel and cell, compute conservative upper-tail ranks:

```
p[b,k,c] = count_j(T[j,k,c] >= T[b,k,c]) / (B + 1)
m[b,c]   = min_k p[b,k,c]
p_adaptive[c] = count_j(m[j,c] <= m[0,c]) / (B + 1)
```

Both stages include the observed row. Near ties are handled symmetrically with
a 1e-10 relative scale tolerance. Every row is processed by the same function;
the null distribution includes the best-kernel selection. Reporting only
`selected_p_uncorrected` is invalid adaptive inference. Duplicating a kernel
does not change the calibrated result. P values cannot be zero.

Given exchangeable rows, the final rank is conservative. This is a statement
about the calibration function conditional on a valid randomization design,
**not a proof that real genomes satisfy the circular null**. Covariate fitting,
non-stationarity and irregular gene spacing can break that design.

## Multiplicity and uncertainty

Return per-cell adaptive P, Holm family-wise adjustment, BY FDR adjustment
(`BH * harmonic_number(m)`, clipped to 1), and explicitly exploratory BH.
BY/Holm do not require independence BETWEEN valid marginal tests. They do NOT
repair a misspecified genomic null or analyst selection outside the fixed family.
Related immune cell labels do not constitute independent replications.

The real analysis corrects 48 primary tests jointly (24 classes x two mappings).
MHC/chr19 exclusions are separate exploratory robustness checks, not additional
opportunities to declare the primary successful. A fixed higher-precision seed
check is reported even when it changes a threshold decision. Do not repeat
resampling until significance is reached. The minimum P is 1/(B+1); results near
the boundary need substantially more resampling and independent follow-up.

## LD and biological limitations

- A MAGMA gene Z input incorporates its original gene-level/reference analysis.
  The adaptive aggregation does NOT use MAGMA's gene-gene covariance matrix.
- Circular shifts preserve gene-order structure but not arbitrary physical LD,
  ancestry mismatch, overlapping-gene dependence or distal Hi-C dependencies.
- Residualization removes the modeled nuisance span, not every nonlinear or
  unmeasured confounder. Per-chromosome fitted residuals need not be stationary.
- Donor composition, cell-state mixtures, exposure and expression uncertainty
  are not represented by the current pseudobulk point weights.
- Real-data gene support tables are descriptive prioritization only. Their
  scores are not gene P values, causal probabilities or evidence of novelty.
- Fine-mapping, allele harmonization, colocalization and independent replication
  remain necessary before claiming a previously unknown biological association.

For an LD-aware comparator, `validation.run_adaptive_magma` actually runs MAGMA
property regression using its cached `.genes.raw` covariance information on the
same measured universe. It includes both continuous and binary weights, controls
the same eight input nuisances plus MAGMA internal covariates, and corrects all
96 tests together. It is a different model, not an exact numerical replica of
the circular test. Agreement is evidence of robustness; disagreement is retained.

## API and CLI

```
gene  chrom  pos  z  Immune  Micro.PVM  mean_expression
g1    1      100  1  0.4     0.3        2.1
...
```

Require one row per gene, finite complete values, non-negative weights, integer
positions, at least 8 genes/chromosome and 4 residual degrees of freedom. Gene
position ordering must be genomic and within one build. Alias duplicates fail.

```python
from celltogwas import run_adaptive_gene_enrichment

result = run_adaptive_gene_enrichment(
    genes,
    cell_types=["Immune", "Micro.PVM"],
    covariates=["mean_expression"],
    n_permutations=19_999,
    random_state=13,
)
print(result.table)
```

```text
celltogwas adaptive-gene --genes genes.tsv --cell-types Immune Micro.PVM \
  --covariates mean_expression --permutations 19999 --seed 13 \
  --output adaptive.tsv --components-output components.tsv --report run.json
```

Backslashes show command continuation only; use one line in PowerShell.
The original SNP-level `enrich` command and its validated 1.0 algorithms remain
unchanged. The core still depends only on NumPy and pandas.

## Sources and scope of inspiration

- [SKAT-O / adaptive sparse-dense testing](https://pmc.ncbi.nlm.nih.gov/articles/PMC3440237/):
  motivation for architecture-adaptive tests, NOT a reimplementation of SKAT-O.
- [MAGMA](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004219):
  gene-level strength, competitive regression and LD-aware comparator.
- [MAGMA official manual](https://cncr.nl/research/magma/), version 1.10: the
  actual local manual governs command flags and default nuisance adjustments.
- [Phipson and Smyth](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf):
  finite randomization tests and inclusion of the observed configuration.
- [Resampling-based minP](https://statistics.berkeley.edu/sites/default/files/tech-reports/633.pdf):
  joint-null calibration; this implementation uses a single symmetric minP
  calibration, NOT the full Westfall-Young step-down procedure.
- [Benjamini and Yekutieli](https://www.math.tau.ac.il/~ybenja/depApr27.pdf):
  dependence-robust multiple testing conditional on valid marginal P values.

ACAT was considered but not used as a universal arbitrary-dependence guarantee;
tail approximations and dependence assumptions would need separate validation.
