# Filter cells

Filter low quality cells based on gene and transcript count thresholds.

## Usage

``` r
filter_cells(xenium.obj, feature_quantile = 0.02, count_quantile = 0.02)
```

## Arguments

- xenium.obj:

  Imported Xenium object.

- feature_quantile:

  Quantile threshold used to filter cells based on number of detected
  features.

- count_quantile:

  Quantile threshold used to filter cells based on number of detected
  transcript counts.

## Value

A filtered Seurat object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Filter low quality cells
xenium.obj <- filter_cells(xenium.obj)} # }
```
