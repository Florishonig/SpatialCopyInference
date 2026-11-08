
# SpatialCopyInference

SpatialCopyInference infers chromosome arm-level somatic copy number
aberrations (SCNAs) directly from Xenium spatial transcriptomics data.
The package identifies and spatially maps SCNA-defined cell populations
through spatial clustering and provides downstream tools for
visualisation, differential expression analysis, pathway enrichment and
identification of genes and genomic regions associated with inferred
SCNAs. SpatialCopyInference enables the study of tumour evolution and
early cancer development in situ.

Documentation, including the entire workflow using a human ovarian
cancer dataset, is available
[here](https://florishonig.github.io/SpatialCopyInference/).

You can read more about SpatialCopyInference in our
[preprint](https://www.biorxiv.org/content/10.64898/2026.10.23.761433v1.article-info).

## Installation

You can install SpatialCopyInference with the following command:

``` r
remotes::install_github("florishonig/SpatialCopyInference")
```

## Example workflow

A minimum workflow is shown below:

``` r
# Import data
xenium.obj <- import_xenium(data_dir, annotation_file)

# Quality control
xenium.obj <- filter_cells(xenium.obj)
xenium.obj <- filter_genes(xenium.obj)

# Normalisation
xenium.obj <- normalise_data(xenium.obj)

# SCNA inference
arm_means <- compute_arm_means(xenium.obj)
arm_ratios <- compute_arm_ratios(arm_means, xenium.obj)

# Spatial clustering
BANKSY.obj <- compute_spatial_clusters(arm_ratios, xenium.obj)

xenium.obj <- add_banksy_clusters(xenium.obj, BANKSY.obj)

# Visualisation
plot_cluster_heatmaps(BANKSY.obj)
plot_spatial_clusters(xenium.obj)

# Downstream analyses
cluster_results <- compare_spatial_clusters(xenium.obj, ident_1 = "1", ident_2 = "2")
plot_gene_density(de_results = cluster_results$de)
gene_clusters <- identify_gene_clusters(de_results = cluster_results$de)
plot_gene_clusters(gene_clusters = gene_clusters)
```
