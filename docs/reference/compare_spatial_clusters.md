# Compare spatial clusters

Performs differential expression and gene ontology enrichment analysis
between two BANKSY-derived spatial clusters.

Differentially expressed genes are identified using
Seurat::FindMarkers(), followed by Gene Ontology (GO) enrichment
analysis of genes upregulated in the first cluster. GO enrichment
results are visualised using dotplot and cnetplot summaries.

## Usage

``` r
compare_spatial_clusters(
  xenium.obj,
  ident_1,
  ident_2,
  cluster_col = "banksy_clusters",
  assay = "Xenium",
  slot = "data",
  p_adj_cutoff = 0.01,
  show_category = 10
)
```

## Arguments

- xenium.obj:

  Xenium object including BANKSY cluster metadata

- ident_1:

  Cluster identifier of interest

- ident_2:

  Reference cluster identifier

- cluster_col:

  Metadata column containing spatial cluster assignments

- assay:

  Assay used for differential expression

- slot:

  Data slot used for differential expression

- p_adj_cutoff:

  Adjusted p-value threshold

- show_category:

  Number of GO terms displayed in enrichment plots

## Value

A list containing differential expression results, GO enrichment results
and associated visualisations.

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare Cluster 1 vs 2
cluster1v2 <- compare_spatial_clusters(
  xenium.obj,
  ident_1 = "1",
  ident_2 = "2")} # }
```
