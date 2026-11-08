# Identify gene clusters

Identifies genomic clusters of differentially expressed genes.

Differentially expressed genes are mapped to their genomic coordinates
and grouped into clusters based on genomic proximity. Genes separated by
less than `gapwidth` base pairs are assigned to the same cluster.

## Usage

``` r
identify_gene_clusters(
  de_results,
  p_adj_cutoff = 0.01,
  logfc_cutoff = 0,
  gapwidth = 2e+06
)
```

## Arguments

- de_results:

  Differential expression results returned by compare_spatial_clusters()

- p_adj_cutoff:

  Adjusted p-value threshold

- logfc_cutoff:

  Minimum log2 fold-change used to define differentially expressed genes

- gapwidth:

  Maximum genomic distance between genes for inclusion in the same
  cluster

## Value

A list containing cluster assignments for differentially expressed genes
(`cluster_assignments`) and the corresponding genomic ranges (`gr_de`).
