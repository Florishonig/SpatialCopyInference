#' @title Plot spatial clusters
#'
#' @description
#' Visualises BANKSY-derived spatial clusters within a
#' Xenium sample.
#'
#' Cluster assignments stored in Xenium metadata are
#' displayed as a combined spatial map and as separate
#' facet panels for each cluster.
#'
#' @param xenium.obj Xenium object including BANKSY
#' cluster metadata
#' @param cluster_col Metadata column containing cluster
#' assignments
#' @param image Image name used to retrieve spatial
#' coordinates
#' @param cluster_palette Optional named colour palette
#' for cluster visualisation
#'
#' @returns A list containing a combined spatial cluster
#' plot and a faceted spatial cluster plot.
#'
#' @examples
#' \dontrun{
#' # Plot spatial clusters with custom colour palette
#' plots <- plot_spatial_clusters(
#'   xenium.obj,
#'   cluster_palette = c("1"="#8CD17D","2"="#F28E2B",
#'   "3"="#4E79A7","4"="#F1CE63"))}
#'
#' @export
plot_spatial_clusters <- function(
    xenium.obj,
    cluster_col = "banksy_clusters",
    image = "fov",
    cluster_palette = NULL
) {

  df_plot <- cbind(
    xenium.obj@meta.data,
    Seurat::GetTissueCoordinates(xenium.obj, image = image)
  )

  # All clusters together
  p_clusters <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = x,
      y = y,
      colour = .data[[cluster_col]]
    )
  ) +
    ggplot2::geom_point(size = 0.3, stroke = 0) +
    ggplot2::coord_equal() +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_void() +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = 4),
        ncol = 1
      )
    ) +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_blank(),
      legend.key.height = grid::unit(2, "mm"),
      legend.key.width = grid::unit(1, "mm")
    )

  if (!is.null(cluster_palette)) {
    p_clusters <- p_clusters +
      ggplot2::scale_color_manual(values = cluster_palette)
  }

  # One facet per cluster
  p_facets <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = x,
      y = y,
      colour = .data[[cluster_col]]
    )
  ) +
    ggplot2::geom_point(size = 0.3, stroke = 0) +
    ggplot2::coord_equal() +
    ggplot2::scale_y_reverse() +
    ggplot2::facet_wrap(
      as.formula(paste("~", cluster_col)),
      ncol = 4
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(size = 8)
    )

  if (!is.null(cluster_palette)) {
    p_facets <- p_facets +
      ggplot2::scale_color_manual(values = cluster_palette)
  }

  list(
    spatial_clusters_plot = p_clusters,
    spatial_clusters_plot_facets = p_facets
  )
}
