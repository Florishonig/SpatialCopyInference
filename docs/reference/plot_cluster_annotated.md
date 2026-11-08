# Plot annotated spatial cluster

Visualises the cellular composition of a BANKSY-derived spatial cluster.

Cells belonging to the selected spatial cluster are displayed in their
spatial context and coloured according to a metadata annotation.

## Usage

``` r
plot_cluster_annotated(
  xenium.obj,
  cluster_id,
  cluster_col = "banksy_clusters",
  group_col = "group",
  image = "fov",
  group_palette = NULL,
  legend_order = NULL
)
```

## Arguments

- xenium.obj:

  Xenium object including BANKSY cluster metadata

- cluster_id:

  Cluster identifier to visualise

- cluster_col:

  Metadata column containing spatial cluster assignments.

- group_col:

  Metadata column used for annotation.

- image:

  Image name used to retrieve spatial coordinates

- group_palette:

  Optional named colour palette for annotation groups

- legend_order:

  Optional legend ordering

## Value

A ggplot object showing the selected spatial cluster coloured by
annotation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot annotated cell types within spatial cluster
plot_cluster_annotated(xenium.obj = xenium.obj,
cluster_id = 1)} # }
```
