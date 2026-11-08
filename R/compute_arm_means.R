#' @title Compute chromosome arm means
#'
#' @description
#' Computes chromosome arm-level expression profiles by
#' averaging normalised gene expression values across
#' genes assigned to the same chromosome arm.
#'
#' Gene annotation is obtained from the packaged
#' dataset hg38_gene_info.
#'
#' @param xenium.obj Normalised Xenium object.
#' @param min_genes_per_arm Minimum number of genes
#' required for a chromosome arm to be retained.
#'
#' @returns A matrix of chromosome arm mean expression
#' values with chromosome arms as rows and cells as
#' columns.
#'
#' @examples
#' \dontrun{
#' # Compute arm means
#' arm_means_matrix <- compute_arm_means(xenium.obj,
#' min_genes_per_arm = 5)}
#'
#'
#' @export
compute_arm_means <- function(
    xenium.obj,
    min_genes_per_arm = 5
) {

  # Load built-in hg38 annotation
  gene_info <- hg38_gene_info %>%
    dplyr::filter(
      hgnc_symbol %in% rownames(xenium.obj)
    )

  # Count genes per chromosome arm
  genes_per_arm_df <- gene_info %>%
    dplyr::count(chrom, arm, name = "n_genes") %>%
    dplyr::mutate(
      chrom = factor(
        chrom,
        levels = c(
          paste0("chr", 1:22),
          "chrX"
        )
      )
    ) %>%
    dplyr::arrange(chrom, arm)

  # Keep only arms with enough genes
  valid_arms <- genes_per_arm_df %>%
    dplyr::filter(
      n_genes >= min_genes_per_arm
    ) %>%
    dplyr::transmute(
      chrom,
      arm
    )

  gene_info <- gene_info %>%
    dplyr::semi_join(
      valid_arms,
      by = c("chrom", "arm")
    )

  # Match genes to matrix
  log_mat <- xenium.obj@assays$Xenium$data[
    gene_info$hgnc_symbol,
    ,
    drop = FALSE
  ]

  # Create chromosome arm groups
  groups <- split(
    seq_len(nrow(log_mat)),
    paste(
      gene_info$chrom,
      gene_info$arm,
      sep = "_"
    )
  )

  # Compute arm means
  arm_means_matrix <- do.call(
    rbind,
    lapply(
      names(groups),
      function(g) {

        matrix(
          Matrix::colMeans(
            log_mat[
              groups[[g]],
              ,
              drop = FALSE
            ]
          ),
          nrow = 1,
          dimnames = list(
            g,
            colnames(log_mat)
          )
        )
      }
    )
  )

  # Sort chromosome arms
  arm_order <- paste(
    valid_arms$chrom,
    valid_arms$arm,
    sep = "_"
  )

  arm_means_matrix <- arm_means_matrix[
    arm_order,
    ,
    drop = FALSE
  ]

  return(arm_means_matrix)
}
