#' @title Plot cluster composition
#'
#' @description
#' Visualises the composition of cell groups across
#' BANKSY-derived spatial clusters.
#'
#' For each cell group, the proportion of cells assigned
#' to each spatial cluster is calculated and displayed as
#' a stacked bar plot.
#'
#' Cell groups are ordered by their dominant spatial
#' cluster to highlight similarities in cluster
#' composition.
#'
#' @param xenium.obj Xenium object including BANKSY
#' cluster metadata
#' @param cluster_col Metadata column containing spatial
#' cluster assignments
#' @param group_col Metadata column defining cell groups
#' @param cluster_palette Optional named colour palette
#' for cluster visualisation
#'
#' @returns A ggplot object showing cell group
#' composition across spatial clusters.
#'
#' @examples
#' \dontrun{
#' # Plot cluster composition with custom palette
#'plot_cluster_composition(xenium.obj,
#'                         cluster_palette = c("#8CD17D","#F28E2B",
#'                         "#4E79A7","#F1CE63"))}
#'
#' @export
plot_cluster_composition <- function(
    xenium.obj,
    cluster_col = "banksy_clusters",
    group_col = "group",
    cluster_palette = NULL
) {

  bar_data <- xenium.obj@meta.data %>%
    dplyr::count(
      .data[[group_col]],
      .data[[cluster_col]]
    ) %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::mutate(
      prop = n / sum(n)
    ) %>%
    dplyr::mutate(
      max_prop = max(prop),
      dominant_cluster = .data[[cluster_col]][which.max(prop)]
    ) %>%
    dplyr::ungroup()

  group_order <- bar_data %>%
    dplyr::distinct(
      .data[[group_col]],
      dominant_cluster,
      max_prop
    ) %>%
    dplyr::arrange(
      dominant_cluster,
      dplyr::desc(max_prop)
    ) %>%
    dplyr::pull(.data[[group_col]])

  bar_data <- bar_data %>%
    dplyr::mutate(
      group = factor(
        .data[[group_col]],
        levels = rev(group_order)
      )
    )

  bar_data <- bar_data %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(
      dplyr::desc(prop),
      .by_group = TRUE
    ) %>%
    dplyr::mutate(
      stack_order = dplyr::row_number()
    ) %>%
    dplyr::ungroup()

  bar_data <- bar_data %>%
    dplyr::mutate(
      stack_id = interaction(
        group,
        stack_order,
        drop = TRUE
      )
    )

  p <- ggplot2::ggplot(
    bar_data,
    ggplot2::aes(
      x = group,
      y = prop,
      fill = factor(.data[[cluster_col]])
    )
  ) +
    ggplot2::geom_col(
      ggplot2::aes(group = stack_id),
      colour = "white"
    ) +
    ggplot2::coord_flip() +
    ggplot2::facet_grid(
      dominant_cluster ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent
    ) +
    ggplot2::labs(
      x = "Cell type",
      y = "Cell type distribution across clusters (%)",
      fill = "Cluster"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        colour = "black",
        size = 9
      ),
      axis.text.y = ggplot2::element_text(
        colour = "black",
        size = 9
      ),
      axis.title.x = ggplot2::element_text(
        colour = "black",
        size = 10
      ),
      axis.title.y = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(
        colour = "black",
        size = 9
      ),
      legend.position = "right"
    )

  if (!is.null(cluster_palette)) {
    p <- p +
      ggplot2::scale_fill_manual(
        values = cluster_palette
      )
  }

  return(p)
}
