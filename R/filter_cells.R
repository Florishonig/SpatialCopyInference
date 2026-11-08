#' @title Filter cells
#'
#' @description
#' Filter low quality cells based on gene
#' and transcript count thresholds.
#'
#' @param xenium.obj Imported Xenium object.
#' @param feature_quantile Quantile threshold used to filter cells based
#' on number of detected features.
#' @param count_quantile Quantile threshold used to filter cells based
#' on number of detected transcript counts.
#'
#' @returns A filtered Seurat object.
#'
#' @examples
#' \dontrun{
#' # Filter low quality cells
#' xenium.obj <- filter_cells(xenium.obj)}
#'
#' @export
filter_cells <- function(
    xenium.obj,
    feature_quantile = 0.02,
    count_quantile = 0.02
) {

  xenium.obj <- subset(
    xenium.obj,
    subset =
      nFeature_Xenium >= stats::quantile(
        xenium.obj$nFeature_Xenium,
        feature_quantile
      ) &
      nCount_Xenium >= stats::quantile(
        xenium.obj$nCount_Xenium,
        count_quantile
      )
  )

  return(xenium.obj)
}
