#' @title Plot clustered log ratio heatmap per cell type
#'
#' @description
#' Computes chromosome arm log ratios for each cell type
#' relative to a reference group and visualises the
#' results as a clustered heatmap.
#'
#' Cell types are clustered based on their
#' chromosome arm log ratio profiles, while chromosome
#' arms are displayed in genomic order.
#'
#' @param xenium.obj Normalised Xenium object
#' @param reference_group Character vector of groups used
#' to construct the reference chromosome arm profile.
#' @param group_var Metadata variable defining cell
#' groups.
#'
#' @returns A ComplexHeatmap heatmap object with
#' hierarchically clustered cell types.
#'
#' @examples
#' \dontrun{
#' # Plot clustered heatmap
#' plot_ratio_heatmap(
#'  xenium.obj,
#'  reference_group = c("Fallopian Tube Epithelium",
#'  "Secretory Epithelial Cells"))}
#'
#' @export
plot_ratio_heatmap <- function(
    xenium.obj,
    reference_group,
    group_var = "group"
) {

  # Compute arm means for every cell type
  all_group_refs <- purrr::map_dfr(
    unique(xenium.obj@meta.data[[group_var]]),
    ~ compute_reference(
      xenium.obj = xenium.obj,
      group_var = group_var,
      group_names = .x
    ) %>%
      dplyr::mutate(group = .x)
  )

  # Compute log ratios
  ref_df <- all_group_refs %>%
    dplyr::filter(group %in% reference_group) %>%
    dplyr::group_by(chrom, arm) %>%
    dplyr::summarise(
      ref_mean = mean(mean_arm_expr),
      .groups = "drop"
    )

  arm_logratio <- all_group_refs %>%
    dplyr::left_join(
      ref_df,
      by = c("chrom", "arm")
    ) %>%
    dplyr::mutate(
      logratio = mean_arm_expr - ref_mean,
      chr_arm = paste0(chrom, arm)
    )

  # Add a single reference row
  reference_row <- ref_df %>%
    dplyr::mutate(
      group = "Reference",
      logratio = 0,
      chr_arm = paste0(chrom, arm)
    )

  arm_logratio <- arm_logratio %>%
    dplyr::select(
      group,
      chrom,
      arm,
      chr_arm,
      logratio
    ) %>%
    dplyr::bind_rows(
      reference_row %>%
        dplyr::select(
          group,
          chrom,
          arm,
          chr_arm,
          logratio
        )
    )

  # Create matrix
  mat <- arm_logratio %>%
    dplyr::select(group, chr_arm, logratio) %>%
    tidyr::pivot_wider(
      names_from = chr_arm,
      values_from = logratio
    ) %>%
    as.data.frame()

  rownames(mat) <- mat$group
  mat <- as.matrix(mat[, -1])

  # Order chromosome arms
  arm_order <- arm_logratio %>%
    dplyr::distinct(chrom, arm) %>%
    dplyr::mutate(
      chr_arm = paste0(chrom, arm),
      arm = factor(arm, levels = c("p", "q")),
      chrom_num = match(
        chrom,
        c(paste0("chr", 1:22), "chrX")
      )
    ) %>%
    dplyr::arrange(chrom_num, arm) %>%
    dplyr::pull(chr_arm)

  ComplexHeatmap::Heatmap(
    mat[, arm_order],
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_names_gp = grid::gpar(fontsize = 5.5),
    row_names_gp = grid::gpar(fontsize = 5.5),
    col = circlize::colorRamp2(
      c(-1, 0, 1),
      c("blue", "white", "red")
    ),
    column_names_rot = 60,
    heatmap_legend_param = list(
      title = expression(log[2]*ratio),
      title_gp = grid::gpar(fontsize = 5.5),
      labels_gp = grid::gpar(fontsize = 5.5),
      legend_height = grid::unit(20, "mm"),
      legend_width  = grid::unit(1, "mm"))
  )
}
