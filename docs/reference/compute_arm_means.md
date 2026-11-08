# Compute chromosome arm means

Computes chromosome arm-level expression profiles by averaging
normalised gene expression values across genes assigned to the same
chromosome arm.

Gene annotation is obtained from the packaged dataset hg38_gene_info.

## Usage

``` r
compute_arm_means(xenium.obj, min_genes_per_arm = 5)
```

## Arguments

- xenium.obj:

  Normalised Xenium object.

- min_genes_per_arm:

  Minimum number of genes required for a chromosome arm to be retained.

## Value

A matrix of chromosome arm mean expression values with chromosome arms
as rows and cells as columns.

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute arm means
arm_means_matrix <- compute_arm_means(xenium.obj,
min_genes_per_arm = 5)} # }

```
