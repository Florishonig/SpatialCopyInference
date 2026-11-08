# Adapted from SeuratWrappers::RunBanksy()
# Modified to use SeuratObject::SetAssayData()
run_banksy_internal <- function (object, lambda, assay = "RNA", slot = "data", use_agf = FALSE,
          dimx = NULL, dimy = NULL, dimz = NULL, ndim = 2, features = "variable",
          group = NULL, split.scale = TRUE, k_geom = 15, n = 2, sigma = 1.5,
          alpha = 0.05, k_spatial = 10, spatial_mode = "kNN_median",
          assay_name = "BANKSY", M = NULL, verbose = TRUE)
{
  check_package_internal(package = "data.table", repository = "CRAN")
  check_package_internal(package = "Matrix", repository = "CRAN")
  check_package_internal(package = "Banksy", repository = "github")
  if (lambda < 0 || lambda > 1)
    stop("Lambda must be between 0 and 1")
  data_own <- get_data_internal(object, assay, slot, features, verbose)
  locs <- get_locs_internal(object, dimx, dimy, dimz, ndim, data_own,
                   group, verbose)
  if (!is.null(group)) {
    object <- Seurat::AddMetaData(object, metadata = locs, col.name = paste0("staggered_",
                                                                     colnames(locs)))
  }
  knn_list <- lapply(k_geom, function(kg) {
    computeNeighbors_internal(locs, spatial_mode = spatial_mode,
                              k_geom = kg, n = n, sigma = sigma, alpha = alpha,
                              k_spatial = k_spatial, verbose = verbose)
  })
  M <- seq(0, max(getM_internal(use_agf, M)))
  center <- rep(TRUE, length(M))
  center[1] <- FALSE
  har <- Map(function(knn_df, M, center) {
    x <- computeHarmonics_internal(data_own, knn_df, M, center,
                                   verbose)
    rownames(x) <- paste0(rownames(x), ".m", M)
    x
  }, knn_list, M, center)
  lambdas <- getLambdas_internal(lambda, n_harmonics = length(har))
  if (verbose)
    message("Creating Banksy matrix")
  data_banksy <- c(list(data_own), har)
  if (verbose)
    message("Scaling BANKSY matrix. Do not call ScaleData on assay ",
            assay_name)
  data_scaled <- lapply(data_banksy, fast_scaler_internal, object, group,
                        split.scale, verbose)
  data_banksy <- Map(function(lam, mat) lam * mat, lambdas,
                     data_banksy)
  data_scaled <- Map(function(lam, mat) lam * mat, lambdas,
                     data_scaled)
  data_banksy <- do.call(rbind, data_banksy)
  data_scaled <- do.call(rbind, data_scaled)
  if (grepl(pattern = "counts", x = slot)) {
    banksy_assay <- Seurat::CreateAssayObject(counts = data_banksy)
  }
  else {
    banksy_assay <- Seurat::CreateAssayObject(data = data_banksy)
  }
  if (verbose)
    message("Setting default assay to ", assay_name)
  object[[assay_name]] <- banksy_assay
  Seurat::DefaultAssay(object) <- assay_name
  object <- SeuratObject::SetAssayData(object, layer = "scale.data", new.data = data_scaled,
                         assay = assay_name)
  object <- Seurat::LogSeuratCommand(object = object)
  return(object)
}
