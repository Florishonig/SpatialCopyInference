#' @title Plot cells
#'
#' @description Plot annotated cell types in their spatial locations.
#'
#' @param xenium.obj Filtered Xenium object.
#' @param group.by Metadata column used to colour cells.
#'
#' @returns A ggplot2 object showing the spatial distribution
#' of annotated cell types.
#'
#' @examples
#' \dontrun{
#' # Plot annotated cells
#' plot_cells(xenium.obj)}
#'
#' @export
plot_cells <- function(
    xenium.obj,
    group.by = "group"
) {

  # Create plotting data frame
  plot_df <- cbind(
    xenium.obj@meta.data,
    Seurat::GetTissueCoordinates(
      xenium.obj,
      image = "fov"
    )
  )

  # Tableau 20 colour palette
  cols <- ggthemes::tableau_color_pal("Tableau 20")(20)

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = x,
      y = y,
      colour = .data[[group.by]]
    )
  ) +
    ggplot2::geom_point(
      size = 1/3,
      stroke = 0
    ) +
    ggplot2::coord_equal() +
    ggplot2::scale_y_reverse() +
    ggplot2::scale_color_manual(
      values = cols
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = 4),
        ncol = 1
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_blank(),
      legend.key.height = grid::unit(2, "mm"),
      legend.key.width = grid::unit(1, "mm")
    )
}
