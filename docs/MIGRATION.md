# Migration from 0.1

## Compatibility

The historical module path and function names remain valid:

```python
from celltogwas.cell2gwas import (
    extract_gene_info,
    merge_intervals,
    celltype_genome_range,
    check_value_in_intervals,
    gwas_annotation,
    gwas_annotaion,
)
```

The functions also remain available from `celltogwas.legacy`.

## Important corrections

The 0.1 code selected transcript-like rows by matching
`gene_name "..."` followed by `transcript_type`. Version 1.0 selects the
`gene` feature by default, so a gene is not counted once per transcript.

The old interval code added a `chr` prefix by stripping characters from the
front of a chromosome string. Version 1.0 uses explicit chromosome
normalization and supports common aliases such as `chr1`, `CHR1`, and `chrM`.

The old GWAS loop assumed that the first two columns were chromosome and
position and could raise `KeyError` when a chromosome was absent from a
cell-type dictionary. The new API resolves column names, preserves rows, and
assigns zero to missing chromosomes. The legacy wrapper still uses the first
two columns for compatibility, but it no longer raises on a missing chromosome.

Intervals are now consistently 0-based and half-open. GWAS positions remain
one-based at the input boundary and are converted internally. Legacy functions
still expose their historical one-based `START/END` frame and string labels.

## Recommended migration

Replace:

```python
ranges = celltype_genome_range(cell_gene_dict, gene_info, 100000)
annotated = gwas_annotation(gwas, ranges)
```

with:

```python
from celltogwas import GenomeLinker, LinkConfig

linker = GenomeLinker.from_paths(
    "genes.gtf.gz",
    "cell_type_genes.gmt",
    config=LinkConfig(window=100000),
)
annotated = linker.annotate(gwas, mode="weight")
```

The new output uses continuous weights and automatically detected GWAS
columns. Use `mode="binary"` if downstream code needs the historical 0/1
semantics.

## Output differences

| 0.1 | 1.0 |
| --- | --- |
| Hard-window membership | Continuous weighted link score |
| String `"1"` / `"0"` only | `weight`, `binary`, or `both` |
| No expression specificity | Tau, entropy, and robust z-score weights |
| No uncertainty estimate | Spatial circular-permutation p-values and FDR |
| No run metadata | JSON manifests with input hashes and diagnostics |

## Notes for reproducible analyses

Record the exact GTF/GFF release, genome assembly, gene-set source, expression
normalization method, window and kernel settings, permutation count, and random
seed. The CLI report contains these fields plus input hashes; keep it with the
result table.
