# Compute reference

Computes a chromosome arm-level reference expression profile from a set
of reference cells.

Expression values are averaged across selected cells and subsequently
aggregated at the chromosome arm level using the packaged dataset
hg38_gene_info.

## Usage

``` r
compute_reference(
  xenium.obj,
  group_var = NULL,
  group_names = NULL,
  min_genes_per_arm = 5
)
```

## Arguments

- xenium.obj:

  Normalised Xenium object.

- group_var:

  Optional metadata variable used to define reference cells. If NULL,
  all cells are used.

- group_names:

  Groups within group_var to include in the reference profile.

- min_genes_per_arm:

  Minimum number of genes required for a chromosome arm to be retained.

## Value

A data frame containing chromosome arm-level reference expression values
and genomic coordinates.
