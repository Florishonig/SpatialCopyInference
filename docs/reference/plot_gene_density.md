# Plot genomic density of differentially expressed genes

Visualises the genomic density of differentially expressed genes across
the genome.

Differentially expressed genes are mapped to their genomic coordinates
and summarised using a rolling density score across chromosome arms.
Results are displayed as a genome-wide karyotype density map, with
optional annotation of known driver genes.

## Usage

``` r
plot_gene_density(
  de_results,
  driver_genes = NULL,
  raised_driver_genes = NULL,
  p_adj_cutoff = 0.01,
  logfc_cutoff = 0,
  window_size = 5e+06,
  step_size = 5e+05
)
```

## Arguments

- de_results:

  Differential expression results returned by compare_spatial_clusters()

- driver_genes:

  Optional vector of driver genes to annotate

- raised_driver_genes:

  Optional subset of driver genes displayed above the karyotype to
  reduce label overlap

- p_adj_cutoff:

  Adjusted p-value threshold

- logfc_cutoff:

  Minimum log2 fold-change used to define differentially expressed genes

- window_size:

  Genomic window size used for density estimation

- step_size:

  Step size used for density estimation and smoothing

## Value

An invisible list containing the smoothed genomic signal, differentially
expressed genes and karyoploteR object used for visualisation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Custom driver genes
driver_genes <- c("MYC", "PIK3CA", "CCNE1", "KRAS",
                 "TERT","AKT1", "AKT2", "AKT3", "CCND1",
                 "CCND2", "CCND3", "KRAS", "MECOM")

# Plot gene density with (raised) driver genes
plot_gene_density_obj <- plot_gene_density(
 de_results = cluster1v2$de,
 driver_genes = driver_genes,
 raised_driver_genes = c("MECOM","AKT2"))} # }
```
