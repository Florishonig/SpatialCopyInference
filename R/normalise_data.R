#' @title Normalise data
#'
#' @description
#' Normalises expression data using Seurat::NormalizeData
#' with the "LogNormalize" method and subsequently convert
#' the resulting values to log2 scale.
#'
#' @param xenium.obj Filtered Xenium object.
#'
#' @returns A log2-normalised Xenium object.
#'
#' @examples
#' \dontrun{
#' # Perform normalisation
#' xenium.obj <- normalise_data(xenium.obj)}
#'
#' @export
normalise_data <- function(
    xenium.obj
) {

  # Perform log normalisation
  xenium.obj <- Seurat::NormalizeData(
    xenium.obj,
    normalization.method = "LogNormalize"
  )

  # Convert log1p(x) to log2(x)
  xenium.obj@assays$Xenium$data <-
    xenium.obj@assays$Xenium$data / log(2)

  return(xenium.obj)
}
