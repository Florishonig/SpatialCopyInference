# Plot clustered log ratio heatmap per cell type

Computes chromosome arm log ratios for each cell type relative to a
reference group and visualises the results as a clustered heatmap.

Cell types are clustered based on their chromosome arm log ratio
profiles, while chromosome arms are displayed in genomic order.

## Usage

``` r
plot_ratio_heatmap(xenium.obj, reference_group, group_var = "group")
```

## Arguments

- xenium.obj:

  Normalised Xenium object

- reference_group:

  Character vector of groups used to construct the reference chromosome
  arm profile.

- group_var:

  Metadata variable defining cell groups.

## Value

A ComplexHeatmap heatmap object with hierarchically clustered cell
types.

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot clustered heatmap
plot_ratio_heatmap(
 xenium.obj,
 reference_group = c("Fallopian Tube Epithelium",
 "Secretory Epithelial Cells"))} # }
```
