#' @title Compare spatial clusters
#'
#' @description
#' Performs differential expression and gene ontology
#' enrichment analysis between two BANKSY-derived
#' spatial clusters.
#'
#' Differentially expressed genes are identified using
#' Seurat::FindMarkers(), followed by Gene Ontology (GO)
#' enrichment analysis of genes upregulated in the first
#' cluster. GO enrichment results are visualised using
#' dotplot and cnetplot summaries.
#'
#' @param xenium.obj Xenium object including BANKSY
#' cluster metadata
#' @param ident_1 Cluster identifier of interest
#' @param ident_2 Reference cluster identifier
#' @param cluster_col Metadata column containing spatial
#' cluster assignments
#' @param assay Assay used for differential expression
#' @param slot Data slot used for differential
#' expression
#' @param p_adj_cutoff Adjusted p-value threshold
#' @param show_category Number of GO terms displayed in
#' enrichment plots
#'
#' @returns A list containing differential expression
#' results, GO enrichment results and associated
#' visualisations.
#' @examples
#' \dontrun{
#' # Compare Cluster 1 vs 2
#' cluster1v2 <- compare_spatial_clusters(
#'   xenium.obj,
#'   ident_1 = "1",
#'   ident_2 = "2")}
#'
#' @export
compare_spatial_clusters <- function(
    xenium.obj,
    ident_1,
    ident_2,
    cluster_col = "banksy_clusters",
    assay = "Xenium",
    slot = "data",
    p_adj_cutoff = 0.01,
    show_category = 10
) {

  # Set identities
  Seurat::Idents(xenium.obj) <- xenium.obj@meta.data[[cluster_col]]

  # Differential expression
  de <- Seurat::FindMarkers(
    xenium.obj,
    ident.1 = ident_1,
    ident.2 = ident_2,
    assay = assay,
    slot = slot,
    test.use = "wilcox"
  )

  # Format output
  de_out <- de %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      log10_padj = -log10(p_val_adj)
    )

  # Significant genes upregulated in ident_1
  de_up <- de %>%
    dplyr::filter(
      p_val_adj < p_adj_cutoff,
      avg_log2FC >= 0
    )

  # Convert SYMBOL -> ENTREZ
  genes_up <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = rownames(de_up),
    columns = "ENTREZID",
    keytype = "SYMBOL"
  ) %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::pull(ENTREZID)

  if (length(genes_up) == 0) {
    stop(
      "No significantly upregulated genes identified."
    )
  }

  # GO enrichment
  ego <-  clusterProfiler::enrichGO(
    gene = genes_up,
    keyType = "ENTREZID",
    OrgDb =   org.Hs.eg.db::org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {

    warning(
      "No significantly enriched GO terms identified."
    )

    return(
      list(
        de = de,
        de_out = de_out,
        ego = ego,
        go_dotplot = NULL,
        go_cnetplot = NULL
      )
    )

  }

  # GO dotplot
  go_dotplot <- enrichplot::dotplot(
    ego,
    showCategory = show_category
  ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 9),
      axis.text.y = ggplot2::element_text(colour = "black", size = 9),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(colour = "black", size = 9),
      axis.title.x = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(size = 9)
    ) +
    ggplot2::ggtitle(
      paste0(
        "GO: Cluster ",
        ident_1,
        " vs ",
        ident_2
      )
    ) +
    ggplot2::scale_y_discrete(
      labels = function(x)
        stringr::str_wrap(x, width = 100)
    )

  # GO cnetplot
  go_cnetplot <- enrichplot::cnetplot(
    ego,
    showCategory = show_category
  )

  list(
    de = de,
    de_out = de_out,
    ego = ego,
    go_dotplot = go_dotplot,
    go_cnetplot = go_cnetplot
  )
}
