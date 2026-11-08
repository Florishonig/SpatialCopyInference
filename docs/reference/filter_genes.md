# Filter genes

Filter lowly expressed genes to reduce noise and improve computational
efficiency.

## Usage

``` r
filter_genes(xenium.obj, min_cells = 50, min_counts = 50)
```

## Arguments

- xenium.obj:

  Imported and cell filtered Xenium object.

- min_cells:

  Minimum number of cells in which a gene must be detected to be
  retained.

- min_counts:

  Minimum total transcript count required for a gene to be retained.

## Value

A filter Seurat object.

## Examples

``` r
if (FALSE) { # \dontrun{
xenium.obj <- filter_genes(xenium.obj,
                          min_cells = 50,
                          min_counts = 50)} # }
```
