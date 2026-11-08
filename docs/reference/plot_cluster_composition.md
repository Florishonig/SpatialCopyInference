# Plot cluster composition

Visualises the composition of cell groups across BANKSY-derived spatial
clusters.

For each cell group, the proportion of cells assigned to each spatial
cluster is calculated and displayed as a stacked bar plot.

Cell groups are ordered by their dominant spatial cluster to highlight
similarities in cluster composition.

## Usage

``` r
plot_cluster_composition(
  xenium.obj,
  cluster_col = "banksy_clusters",
  group_col = "group",
  cluster_palette = NULL
)
```

## Arguments

- xenium.obj:

  Xenium object including BANKSY cluster metadata

- cluster_col:

  Metadata column containing spatial cluster assignments

- group_col:

  Metadata column defining cell groups

- cluster_palette:

  Optional named colour palette for cluster visualisation

## Value

A ggplot object showing cell group composition across spatial clusters.

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot cluster composition with custom palette
plot_cluster_composition(xenium.obj,
                        cluster_palette = c("#8CD17D","#F28E2B",
                        "#4E79A7","#F1CE63"))} # }
```
