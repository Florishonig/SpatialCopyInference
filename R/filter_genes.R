#' @title Filter genes
#'
#' @description
#' Filter lowly expressed genes to reduce noise
#' and improve computational efficiency.
#'
#' @param xenium.obj Imported and cell filtered Xenium object.
#' @param min_cells Minimum number of cells in which a
#' gene must be detected to be retained.
#' @param min_counts Minimum total transcript count
#' required for a gene to be retained.
#'
#' @returns A filter Seurat object.
#'
#' @examples
#' \dontrun{
#` # Filter lowly expressed genes
#'xenium.obj <- filter_genes(xenium.obj,
#'                           min_cells = 50,
#'                           min_counts = 50)}
#' @export
filter_genes <- function(
    xenium.obj,
    min_cells = 50,
    min_counts = 50
) {

  gene_counts <- Matrix::rowSums(
    xenium.obj@assays$Xenium$counts
  )

  keep_genes <- names(gene_counts)[
    Matrix::rowSums(
      xenium.obj@assays$Xenium$counts > 0
    ) >= min_cells &
      gene_counts >= min_counts
  ]

  xenium.obj@misc$removed_genes <- setdiff(
    rownames(xenium.obj),
    keep_genes
  )

  xenium.obj <- subset(
    x = xenium.obj,
    features = keep_genes
  )

  return(xenium.obj)
}
