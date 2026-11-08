# Plot cells

Plot annotated cell types in their spatial locations.

## Usage

``` r
plot_cells(xenium.obj, group.by = "group")
```

## Arguments

- xenium.obj:

  Filtered Xenium object.

- group.by:

  Metadata column used to colour cells.

## Value

A ggplot2 object showing the spatial distribution of annotated cell
types.

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot annotated cells
plot_cells(xenium.obj)} # }
```
