#' @title Import Xenium data
#'
#' @description
#' Imports Xenium spatial transcriptomics data using
#' Seurat::LoadXenium() and optionally adds cell type
#' annotations from a CSV file.
#'
#' @param data_dir Path to Xenium directory containing default Xenium output files.
#' @param annotation_file Optional path to a CSV file containing cell type annotations.
#'
#' @returns A Seurat object containing Xenium data and associated metadata.
#'
#' @examples
#' \dontrun{
#' # Import Xenium data only
#' xenium.obj <- import_xenium(
#'   data_dir = "path/to/Xenium_outputs/outs",
#'   annotation_file = NULL)
#'
#' # Import Xenium data and cell type annotations
#' xenium.obj <- import_xenium(
#'   data_dir = "path/to/Xenium_outputs/outs",
#'   annotation_file = "path/to/cell_type_annotations.csv")}
#'
#' @export
import_xenium <- function(data_dir, annotation_file = NULL) {

  # import Xenium data
  xenium.obj <- Seurat::LoadXenium(
    data.dir = data_dir
  )

  # import annotations if provided
  if (!is.null(annotation_file)) {

    annotation <- utils::read.csv(
      annotation_file,
      row.names = "cell_id"
    )

    xenium.obj <- Seurat::AddMetaData(
      xenium.obj,
      metadata = annotation
    )
  }

  return(xenium.obj)
}
