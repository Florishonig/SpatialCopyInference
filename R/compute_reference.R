#' @title Compute reference
#'
#' @description
#' Computes a chromosome arm-level reference expression
#' profile from a set of reference cells.
#'
#' Expression values are averaged across selected cells
#' and subsequently aggregated at the chromosome arm
#' level using the packaged dataset hg38_gene_info.
#'
#' @param xenium.obj Normalised Xenium object.
#' @param group_var Optional metadata variable used to define
#' reference cells. If NULL, all cells are used.
#' @param group_names Groups within group_var to include
#' in the reference profile.
#' @param min_genes_per_arm Minimum number of genes
#' required for a chromosome arm to be retained.
#'
#' @returns A data frame containing chromosome arm-level
#' reference expression values and genomic coordinates.
#'
#' @export
compute_reference <- function(
    xenium.obj,
    group_var = NULL,
    group_names = NULL,
    min_genes_per_arm = 5
) {

  # Load built-in hg38 annotation
  gene_info <- hg38_gene_info %>%
    dplyr::filter(
      hgnc_symbol %in% rownames(xenium.obj)
    )

  if (nrow(gene_info) == 0) {

    stop(
      "No matching genes found between hg38_gene_info and xenium.obj."
    )

  }

  genes_per_arm_df <- gene_info %>%
    dplyr::count(
      chrom,
      arm,
      name = "n_genes"
    )

  valid_arms <- genes_per_arm_df %>%
    dplyr::filter(
      n_genes >= min_genes_per_arm
    ) %>%
    dplyr::select(
      chrom,
      arm
    )

  gene_info <- gene_info %>%
    dplyr::semi_join(
      valid_arms,
      by = c("chrom", "arm")
    )

  if (is.null(group_var)) {

    # Use all cells
    cells_use <- colnames(xenium.obj)
    group_names <- "All Cells"

  } else {

    if (is.null(group_names)) {

      group_names <- unique(
        xenium.obj@meta.data[[group_var]]
      )

    }

    cells_use <- rownames(
      xenium.obj@meta.data[
        xenium.obj@meta.data[[group_var]] %in% group_names,
        ,
        drop = FALSE
      ]
    )

  }

  if (length(cells_use) == 0) {

    stop(
      "No cells found for the specified group_names"
    )

  }
  expr_mean <- Matrix::rowMeans(
    xenium.obj@assays$Xenium$data[
      ,
      cells_use,
      drop = FALSE
    ]
  )

  data.frame(
    hgnc_symbol = names(expr_mean),
    mean_expression = as.numeric(expr_mean)
  ) %>%
    merge(
      gene_info,
      by = "hgnc_symbol"
    ) %>%
    dplyr::group_by(chrom, arm) %>%
    dplyr::summarise(
      mean_arm_expr = mean(
        mean_expression,
        na.rm = TRUE
      ),
      p_start = dplyr::first(p_start),
      p_end = dplyr::first(p_end),
      q_start = dplyr::first(q_start),
      q_end = dplyr::first(q_end),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      chrom = factor(
        chrom,
        levels = c(
          paste0("chr", 1:22),
          "chrX"
        )
      )
    ) %>%
    dplyr::arrange(chrom, arm) %>%
    dplyr::mutate(
      group = paste(
        group_names,
        collapse = " + "
      )
    )
}
