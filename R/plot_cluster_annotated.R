#' @title Plot annotated spatial cluster
#'
#' @description
#' Visualises the cellular composition of a BANKSY-derived
#' spatial cluster.
#'
#' Cells belonging to the selected spatial cluster are
#' displayed in their spatial context and coloured
#' according to a metadata annotation.
#'
#' @param xenium.obj Xenium object including BANKSY
#' cluster metadata
#' @param cluster_id Cluster identifier to visualise
#' @param cluster_col Metadata column containing spatial
#' cluster assignments.
#' @param group_col Metadata column used for annotation.
#' @param image Image name used to retrieve spatial
#' coordinates
#' @param group_palette Optional named colour palette for
#' annotation groups
#' @param legend_order Optional legend ordering
#'
#' @returns A ggplot object showing the selected spatial
#' cluster coloured by annotation.
#'
#' @examples
#' \dontrun{
#' # Plot annotated cell types within spatial cluster
#' plot_cluster_annotated(xenium.obj = xenium.obj,
#' cluster_id = 1)}
#'
#' @export
plot_cluster_annotated <- function(
    xenium.obj,
    cluster_id,
    cluster_col = "banksy_clusters",
    group_col = "group",
    image = "fov",
    group_palette = NULL,
    legend_order = NULL
) {

  df_plot <- cbind(
    xenium.obj@meta.data,
    Seurat::GetTissueCoordinates(
      xenium.obj,
      image = image
    )
  )

  p <- ggplot2::ggplot(
    df_plot[
      df_plot[[cluster_col]] == cluster_id,
      ,
      drop = FALSE
    ],
    ggplot2::aes(
      x = x,
      y = y
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        colour = .data[[group_col]]
      ),
      size = 0.3
    ) +
    ggplot2::coord_equal() +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_void() +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = 4)
      )
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_blank(),
      legend.key.height = grid::unit(2, "mm"),
      legend.key.width = grid::unit(1, "mm")
    )

  if (!is.null(group_palette)) {

    if (is.null(legend_order)) {

      p <- p +
        ggplot2::scale_color_manual(
          values = group_palette
        )

    } else {

      p <- p +
        ggplot2::scale_color_manual(
          values = group_palette,
          breaks = legend_order
        )

    }

  }

  return(p)

}
