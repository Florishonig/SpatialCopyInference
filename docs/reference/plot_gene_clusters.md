# Plot genomic clusters of differentially expressed genes

Visualises genomic clusters identified using identify_gene_clusters().

Differentially expressed genes are displayed on a genome-wide karyotype
and coloured according to the size of the genomic cluster to which they
belong.

## Usage

``` r
plot_gene_clusters(
  gene_clusters,
  driver_genes = NULL,
  raised_driver_genes = NULL
)
```

## Arguments

- gene_clusters:

  Output from identify_gene_clusters()

- driver_genes:

  Optional vector of driver genes to annotate

- raised_driver_genes:

  Optional subset of driver genes displayed on a separate annotation
  track to reduce label overlap

## Value

An invisible list containing annotated gene ranges, cluster ranges,
driver gene ranges and the karyoploteR object used for visualisation.
