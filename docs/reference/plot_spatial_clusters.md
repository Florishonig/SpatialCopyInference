# Plot spatial clusters

Visualises BANKSY-derived spatial clusters within a Xenium sample.

Cluster assignments stored in Xenium metadata are displayed as a
combined spatial map and as separate facet panels for each cluster.

## Usage

``` r
plot_spatial_clusters(
  xenium.obj,
  cluster_col = "banksy_clusters",
  image = "fov",
  cluster_palette = NULL
)
```

## Arguments

- xenium.obj:

  Xenium object including BANKSY cluster metadata

- cluster_col:

  Metadata column containing cluster assignments

- image:

  Image name used to retrieve spatial coordinates

- cluster_palette:

  Optional named colour palette for cluster visualisation

## Value

A list containing a combined spatial cluster plot and a faceted spatial
cluster plot.

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot spatial clusters with custom colour palette
plots <- plot_spatial_clusters(
  xenium.obj,
  cluster_palette = c("1"="#8CD17D","2"="#F28E2B",
  "3"="#4E79A7","4"="#F1CE63"))} # }
```
