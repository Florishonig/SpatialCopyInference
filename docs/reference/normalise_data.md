# Normalise data

Normalises expression data using Seurat::NormalizeData with the
"LogNormalize" method and subsequently convert the resulting values to
log2 scale.

## Usage

``` r
normalise_data(xenium.obj)
```

## Arguments

- xenium.obj:

  Filtered Xenium object.

## Value

A log2-normalised Xenium object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Perform normalisation
xenium.obj <- normalise_data(xenium.obj)} # }
```
