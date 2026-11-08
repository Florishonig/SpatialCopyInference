# Internal BANKSY helper functions
# Copied to avoid reliance on unexported package internals
check_package_internal <- function (package, repository, ...)
{
  if (!requireNamespace(package = basename(path = package),
                        quietly = TRUE)) {
    if (interactive()) {
      message("Package ", package, " is not yet installed")
      message("Install now?")
      choice <- menu(choices = c("yes", "no"))
      if (choice == 1) {
        repository <- match.arg(arg = tolower(x = repository),
                                choices = c("github", "bioconductor", "cran"))
        switch(EXPR = repository, github = remotes::install_github(repo = package,
                                                                   ...), bioconductor = BiocManager::install(pkgs = package,
                                                                                                             ...), cran = install.packages(pkgs = package,
                                                                                                                                           ...), stop("Unknown repository ", repository,
                                                                                                                                                      call. = FALSE))
        return(invisible(x = NULL))
      }
    }
    stop("Unable to find package ", package, ", please install",
         call. = FALSE)
  }
}

get_data_internal <- function (object, assay, slot, features, verbose)
{
  if (verbose)
    message("Fetching data from slot ", slot, " from assay ",
            assay)
  data_own <- Seurat::GetAssayData(object = object, assay = assay,
                                   layer = slot)
  if (features[1] != "all") {
    if (verbose)
      message("Subsetting by features")
    if (features[1] == "variable") {
      feat <- Seurat::VariableFeatures(object)
      if (length(feat) == 0) {
        warning("No variable features found. Running Seurat::FindVariableFeatures")
        object <- Seurat::FindVariableFeatures(object)
        feat <- Seurat::VariableFeatures(object)
      }
    }
    else {
      feat <- features[which(rownames(object) %in% features)]
      if (length(feat) == 0)
        stop("None of the specified features found. Check if features in Seurat object")
    }
    data_own <- data_own[feat, , drop = FALSE]
  }
  data_own <- as.matrix(x = data_own)
  return(data_own)
}


get_locs_internal <- function (object, dimx, dimy, dimz, ndim, data_own, group, verbose)
{
  if (!is.null(dimx) & !is.null(dimy)) {
    locs <- data.frame(sdimx = unlist(object[[dimx]]), sdimy = unlist(object[[dimy]]))
    rownames(locs) <- colnames(object)
    if (!is.null(dimz))
      locs$sdimz = object[[dimz]]
    obj_samples <- colnames(data_own)
    locs_samples <- rownames(locs)
    if (any(is.na(match(obj_samples, locs_samples)))) {
      na_id <- which(is.na(match(obj_samples, locs_samples)))
      warning("No centroids found for samples: ", paste(obj_samples[na_id],
                                                        collapse = ", "), ". Dropping samples.")
      data_own <- data_own[, -na_id, drop = FALSE]
    }
    locs <- locs[match(obj_samples, locs_samples), , drop = FALSE]
  }
  else {
    locs <- Seurat::GetTissueCoordinates(object)[, seq_len(ndim)]
  }
  dim_names <- paste0("sdim", c("x", "y", "z"))
  colnames(locs) <- dim_names[seq_len(ncol(locs))]
  if (!is.null(group)) {
    if (verbose)
      message("Staggering locations by ", group)
    locs[, 1] = locs[, 1] + abs(min(locs[, 1]))
    max_x = max(locs[, 1]) * 2
    n_groups = length(unique(unlist(object[[group]])))
    shift = seq(from = 0, length.out = n_groups, by = max_x)
    locs[, 1] = locs[, 1] + rep(shift, table(object[[group]]))
  }
  return(locs)
}

fast_scaler_internal <- function (data, object, group, split.scale, verbose)
{
    if (!is.null(group) & split.scale) {
        groups = unlist(object[[group]])
        ugroups = unique(groups)
        for (curr_group in ugroups) {
            if (verbose)
                message("Scaling group: ", curr_group)
            curr_group_id <- which(curr_group == groups)
            data[, curr_group_id] <- Seurat::FastRowScale(data[,
                curr_group_id])
        }
    }
    else {
        data <- Seurat::FastRowScale(data)
    }
    data
}

getM_internal <- function (use_agf, M)
{
  if (is.null(M))
    M <- as.numeric(use_agf)
  sort(M)
}

getLambdas_internal <- function (lambda, n_harmonics)
{
  weights <- lambda * (2^-seq(0, n_harmonics - 1))
  weights <- weights/sum(2^-seq(0, n_harmonics - 1))
  lam <- c(1 - sum(weights), weights)
  sqrt(lam)
}

computeNeighbors_internal <- function (locs, spatial_mode = "kNN_median", k_geom = 15, n = 2,
          sigma = 1.5, alpha = 0.05, k_spatial = 100, sample_size = NULL,
          sample_renorm = TRUE, seed = NULL, dimensions = "all", verbose = FALSE)
{
  from <- to <- phi <- NULL
  locs <- as.matrix(locs)
  kernelRadius <- sqrt(-ncol(locs) * log(alpha))
  if (verbose)
    message("Computing neighbors...")
  knnDF <- switch(spatial_mode, rNN_gauss = withRNNgauss_internal(locs = locs,
                                                         sigma = sigma, k_spatial = k_spatial, kernelRadius = kernelRadius,
                                                         verbose = verbose), kNN_rank = withKNNrank_internal(locs = locs,
                                                                                                    k_geom = k_geom, verbose = verbose), kNN_r = withKNNr_internal(locs = locs,
                                                                                                                                                          k_geom = k_geom, verbose = verbose), kNN_rn = withKNNrn_internal(locs,
                                                                                                                                                                                                                  k_geom = k_geom, n = n, verbose = verbose), kNN_unif = withKNNunif_internal(locs = locs,
                                                                                                                                                                                                                                                                                     k_geom = k_geom, verbose = verbose), kNN_median = withKNNmedian_internal(locs = locs,
                                                                                                                                                                                                                                                                                                                                                     k_geom = k_geom, verbose = verbose))
  knnDF[, `:=`(phi, getPhi_internal(locs, from, to)), by = from][]
  if (!is.null(sample_size)) {
    if (verbose)
      message("Subsampling to ", sample_size, " neighbors")
    knnDF <- subsampler_internal(knnDF, sample_size = sample_size,
                        sample_renorm = sample_renorm, seed = seed)
  }
  if (verbose)
    message("Done")
  return(knnDF)
}

computeHarmonics_internal <- function (gcm, knn_df, M, center, verbose, chunk_size = NULL,
          parallel = FALSE, num_cores = NULL, row_limit_factor = 0.75)
{
  from <- to <- weight <- phi <- .N <- count <- . <- NULL
  j <- sqrt(as.complex(-1))
  mean_k <- round(mean(knn_df[, .(count = .N), by = from]$count),
                  1)
  total_rows <- as.double(nrow(gcm)) * ncol(gcm)
  max_rows <- (2^31 - 1) * row_limit_factor
  if (total_rows > max_rows || !is.null(chunk_size)) {
    if (verbose)
      message("Computing neighborhood matrices in chunks...")
    max_chunk_size <- floor(max_rows/ncol(gcm))
    if (!is.null(chunk_size)) {
      if (chunk_size > max_chunk_size) {
        stop("Specified chunk_size too large. Must be smaller than ",
             floor(max_rows/ncol(gcm)))
      }
      max_chunk_size <- chunk_size
    }
  }
  else {
    max_chunk_size <- nrow(gcm)
  }
  num_chunks <- ceiling(nrow(gcm)/max_chunk_size)
  if (verbose && num_chunks > 1) {
    message("Processing in ", num_chunks, " chunks of max ",
            max_chunk_size, " genes each")
    if (parallel) {
      if (.Platform$OS.type == "windows") {
        message("Parallel processing not supported on Windows - using sequential processing")
      }
      else {
        if (!is.null(num_cores)) {
          message("Using parallel processing with ",
                  num_cores, " cores")
        }
        else {
          message("Using parallel processing with default backend")
        }
      }
    }
  }
  if (parallel) {
    if (.Platform$OS.type == "windows") {
      if (verbose)
        message("Parallel processing not supported on Windows. Using sequential processing.")
      apply_fun <- lapply
    }
    else {
      if (!requireNamespace("BiocParallel", quietly = TRUE)) {
        stop("BiocParallel package is required for parallel processing. Please install it with: BiocManager::install('BiocParallel')")
      }
      if (!is.null(num_cores)) {
        bp_param <- BiocParallel::MulticoreParam(workers = num_cores)
        apply_fun <- function(x, fun) BiocParallel::bplapply(x,
                                                             fun, BPPARAM = bp_param)
      }
      else {
        apply_fun <- BiocParallel::bplapply
      }
    }
  }
  else {
    apply_fun <- lapply
  }
  process_chunk <- function(chunk) {
    start_idx <- (chunk - 1) * max_chunk_size + 1
    end_idx <- min(chunk * max_chunk_size, nrow(gcm))
    if (verbose && num_chunks > 1 && (!parallel || .Platform$OS.type ==
                                      "windows")) {
      message("Processing chunk ", chunk, "/", num_chunks,
              " (rows ", start_idx, " to ", end_idx, ")")
    }
    gcm_chunk <- gcm[start_idx:end_idx, , drop = FALSE]
    if (center) {
      if (verbose && chunk == 1 && (!parallel || .Platform$OS.type ==
                                    "windows"))
        message("Centering")
      chunk_aggr <- knn_df[, abs(fscale_internal(gcm_chunk[, to,
                                                  drop = FALSE]) %*% (weight * exp(j * M * phi))),
                           by = from]
    }
    else {
      chunk_aggr <- knn_df[, abs(gcm_chunk[, to, drop = FALSE] %*%
                                   (weight * exp(j * M * phi))), by = from]
    }
    result_vector <- chunk_aggr$V1
    rm(gcm_chunk, chunk_aggr)
    return(result_vector)
  }
  chunk_results <- apply_fun(1:num_chunks, process_chunk)
  ncm <- matrix(0, nrow = nrow(gcm), ncol = ncol(gcm))
  for (i in seq_along(chunk_results)) {
    start_idx <- (i - 1) * max_chunk_size + 1
    end_idx <- min(i * max_chunk_size, nrow(gcm))
    ncm[start_idx:end_idx, ] <- matrix(chunk_results[[i]],
                                       nrow = end_idx - start_idx + 1, ncol = ncol(gcm))
  }
  rm(chunk_results)
  rownames(ncm) <- rownames(gcm)
  colnames(ncm) <- colnames(gcm)
  if (verbose)
    message("Done")
  return(ncm)
}

withRNNgauss_internal <- function (locs, sigma, k_spatial, kernelRadius, verbose)
{
  if (verbose)
    message("Spatial mode is rNN gaussian")
  if (verbose)
    message("Parameters: sigma=", sigma, ", k_spatial=",
            k_spatial)
  tryCatch({
    knn <- dbscan::kNN(x = locs, k = k_spatial)
  }, error = function(cond) {
    stop("Not enough neighbours at k_spatial=", k_spatial)
  })
  medianDist <- median(as.vector(knn$dist[, 1]))
  knnDF <- data.table(from = rep(seq_len(nrow(knn$id)), k_spatial),
                      to = as.vector(knn$id), weight = dnorm(as.vector(knn$dist),
                                                             mean = 0, sd = medianDist * sigma), distance = as.vector(knn$dist))
  distance <- norm.weight <- weight <- from <- to <- NULL
  knnDF <- knnDF[distance < sigma * kernelRadius * medianDist,
  ]
  setDT(knnDF)[, `:=`(norm.weight, weight/sum(weight)), by = from]
  knnDF <- knnDF[, -3, with = FALSE]
  iso <- setdiff(seq_len(nrow(locs)), unique(knnDF$from))
  isomat <- c(rep(iso, 2), rep(0, 2 * length(iso)))
  isomat <- data.table(matrix(isomat, ncol = ncol(knnDF)))
  knnDF <- rbindlist(list(knnDF, isomat), use.names = FALSE)
  knnDF <- knnDF[order(from, to)]
  setnames(knnDF, "norm.weight", "weight")
  return(knnDF)
}

withKNNrank_internal <- function (locs, k_geom, verbose)
{
  if (verbose)
    message("Spatial mode is kNN_rank")
  if (verbose)
    message("Parameters: k_geom=", k_geom)
  tryCatch({
    knn <- dbscan::kNN(x = locs, k = k_geom)
    unnormWt <- exp(-seq(1, k_geom, 1)^2/(2 * (k_geom/1.5)^2))
    normWt <- unnormWt/sum(unnormWt)
    weightMatrix <- t(matrix(normWt, nrow = k_geom, ncol = nrow(knn$id)))
    knnDF <- data.table(from = rep(seq_len(nrow(knn$id)),
                                   k_geom), to = as.vector(knn$id), weight = as.vector(weightMatrix),
                        distance = as.vector(knn$dist))
  }, error = function(cond) {
    stop("Not enough neighbours at k_geom=", k_geom)
  })
  return(knnDF)
}

withKNNr_internal <- function (locs, k_geom, verbose)
{
  if (verbose)
    message("Spatial mode is kNN_r")
  if (verbose)
    message("Parameters: k_geom=", k_geom)
  tryCatch({
    knn <- dbscan::kNN(x = locs, k = k_geom)
  }, error = function(cond) {
    stop("Not enough neighbours at k_geom=", k_geom)
  })
  norm.weight <- weight <- from <- NULL
  knnDF <- data.table(from = rep(seq_len(nrow(knn$id)), k_geom),
                      to = as.vector(knn$id), weight = 1/as.vector(knn$dist),
                      distance = as.vector(knn$dist))
  knnDF[, `:=`(norm.weight, weight/sum(weight)), by = from]
  knnDF <- knnDF[, -3, with = FALSE]
  setnames(knnDF, "norm.weight", "weight")
  return(knnDF)
}

withKNNrn_internal <- function (locs, k_geom, n, verbose)
{
  if (verbose)
    message("Spatial mode is kNN_rn")
  if (verbose)
    message("Parameters: k_geom=", k_geom, ", n=", n)
  tryCatch({
    knn <- dbscan::kNN(x = locs, k = k_geom)
  }, error = function(cond) {
    stop("Not enough neighbours at k_geom=", k_geom)
  })
  norm.weight <- weight <- from <- NULL
  knnDF <- data.table(from = rep(seq_len(nrow(knn$id)), k_geom),
                      to = as.vector(knn$id), weight = 1/(as.vector(knn$dist)^n),
                      distance = as.vector(knn$dist))
  knnDF[, `:=`(norm.weight, weight/sum(weight)), by = from]
  knnDF <- knnDF[, -3, with = FALSE]
  setnames(knnDF, "norm.weight", "weight")
  return(knnDF)
}

withKNNunif_internal <- function (locs, k_geom, verbose)
{
  if (verbose)
    message("Spatial mode is kNN_unif")
  if (verbose)
    message("Parameters: k_geom = ", k_geom)
  tryCatch({
    knn <- dbscan::kNN(x = locs, k = k_geom)
  }, error = function(cond) {
    stop("Not enough neighbours at k_geom=", k_geom)
  })
  norm.weight <- weight <- from <- NULL
  knnDF <- data.table(from = rep(seq_len(nrow(knn$id)), k_geom),
                      to = as.vector(knn$id), weight = 1, distance = as.vector(knn$dist))
  knnDF[, `:=`(norm.weight, weight/sum(weight)), by = from]
  knnDF <- knnDF[, -3, with = FALSE]
  setnames(knnDF, "norm.weight", "weight")
  return(knnDF)
}

withKNNmedian_internal <- function (locs, k_geom, verbose)
{
    if (verbose)
        message("Spatial mode is kNN_median")
    if (verbose)
        message("Parameters: k_geom=", k_geom)
    tryCatch({
        knn <- dbscan::kNN(x = locs, k = k_geom)
    }, error = function(cond) {
        stop("Not enough neighbours at k_geom=", k_geom)
    })
    norm.weight <- weight <- from <- distance <- NULL
    knnDF <- data.table(from = rep(seq_len(nrow(knn$id)), k_geom),
        to = as.vector(knn$id), distance = as.vector(knn$dist))
    knnDF[, `:=`(weight, exp(-distance^2/median(distance)^2)),
        by = from]
    knnDF[, `:=`(norm.weight, weight/sum(weight)), by = from]
    knnDF <- knnDF[, -c(3, 4), with = FALSE]
    setnames(knnDF, "norm.weight", "weight")
    return(knnDF)
}

getPhi_internal <- function (locs, from, to)
{
  out <- sweep(locs[to, , drop = FALSE], 2, locs[from, , drop = FALSE],
               "-")
  phi <- atan2(out[, 2, drop = FALSE], out[, 1, drop = FALSE])
  phi + as.integer(phi < 0) * 2 * pi
}

subsampler_internal <- function (knnDF, sample_size = NULL, sample_renorm = TRUE, seed = NULL)
{
  verbose.seed_internal(seed)
  from <- weight <- NULL
  x <- knnDF[, .SD[sample(.N, min(sample_size, .N), replace = FALSE)],
             by = from]
  if (sample_renorm)
    x[, `:=`(weight, weight/sum(weight)), by = from]
  data.table(x)
}

fscale_internal <- function (x)
{
  rm <- rowMeans(x)
  x <- (x - rm)
  return(x)
}

verbose.seed_internal <- function (seed)
{
  if (!is.null(seed)) {
    set.seed(seed)
    message("Using seed=", seed)
  }
}
