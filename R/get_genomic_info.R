#' @title Get genomic information
#'
#' @description
#' Retrieves genomic coordinates and chromosome arm
#' annotations for a set of genes from Ensembl and
#' combines them with chromosome size and centromere
#' information.
#'
#' The function returns gene level and chromosome level
#' annotation tables that can be used for chromosome arm
#' spatial copy number inference.
#'
#' This function was used to generate the packaged
#' datasets hg38_gene_info and hg38_chrom_info.
#'
#' @param xenium.obj Xenium object containing genes to annotate.
#' @param genes Character vector of HGNC gene symbols.
#' Either `xenium.obj` or `genes` must be supplied.
#' @param chrom_sizes_file Path to chromosome size file.
#' @param centromere_file Path to centromere annotation file.
#'
#' @returns A list containing gene and chromosome
#' annotation information.
#'
#' @export
get_genomic_info <- function(
    xenium.obj = NULL,
    genes = NULL,
    chrom_sizes_file = system.file(
      "extdata",
      "hg38.chrom.sizes.txt",
      package = "SpatialCopyInference"
    ),
    centromere_file = system.file(
      "extdata",
      "centromeres_hg38.tsv",
      package = "SpatialCopyInference"
    )
) {
  # Determine genes to use
  if (!is.null(xenium.obj)) {

    genes <- rownames(xenium.obj)

  } else if (!is.null(genes)) {

    genes <- unique(genes)

  } else {

    stop(
      "Provide either xenium.obj or genes."
    )

  }

  # Obtain gene coordinates from Ensembl
  mart <- biomaRt::useEnsembl(
    "ensembl",
    dataset = "hsapiens_gene_ensembl",
    mirror = "useast"
  )

  coords <- biomaRt::getBM(
    attributes = c(
      "hgnc_symbol",
      "ensembl_gene_id",
      "chromosome_name",
      "start_position",
      "end_position",
      "band"
    ),
    filters = "hgnc_symbol",
    values = genes,
    mart = mart
  )

  # Keep autosomes and chromosome X
  coords <- coords[
    coords$chromosome_name %in% c(
      as.character(1:22),
      "X"
    ),
  ]

  # Change chromosome naming
  coords$chrom <- paste0(
    "chr",
    coords$chromosome_name
  )

  coords$mid_position <-
    (coords$start_position +
       coords$end_position) / 2

  coords$arm <- substr(
    coords$band,
    1,
    1
  )

  # Chromosome sizes
  hg38_chrom_sizes <- utils::read.delim(
    chrom_sizes_file,
    header = FALSE,
    col.names = c(
      "chrom",
      "length"
    )
  )

  hg38_chrom_sizes <- hg38_chrom_sizes[
    hg38_chrom_sizes$chrom %in%
      c(paste0("chr", 1:22), "chrX"),
  ]

  # Centromeres
  hg38_centromeres <- utils::read.delim(
    centromere_file,
    header = FALSE,
    col.names = c(
      "chrom",
      "centromere_start",
      "centromere_end",
      "other"
    )
  )

  hg38_centromeres <- hg38_centromeres %>%
    dplyr::group_by(chrom) %>%
    dplyr::summarise(
      centromere_start = min(centromere_start),
      centromere_end = max(centromere_end),
      .groups = "drop"
    )

  # Chromosome information
  hg38_chrom_info <- merge(
    hg38_chrom_sizes,
    hg38_centromeres
  )

  # Sort chromosomes
  hg38_chrom_info$chrom <- factor(
    hg38_chrom_info$chrom,
    levels = c(
      paste0("chr", 1:22),
      "chrX"
    )
  )

  hg38_chrom_info <- hg38_chrom_info[
    order(hg38_chrom_info$chrom),
  ]

  rownames(hg38_chrom_info) <- NULL

  hg38_chrom_info$cum_length <-
    c(
      0,
      utils::head(
        cumsum(as.numeric(hg38_chrom_info$length)),
        -1
      )
    )

  hg38_chrom_info$cum_centromere_start <-
    hg38_chrom_info$centromere_start +
    hg38_chrom_info$cum_length

  hg38_chrom_info$cum_centromere_end <-
    hg38_chrom_info$centromere_end +
    hg38_chrom_info$cum_length

  hg38_chrom_info$cum_centromere_mid <-
    (
      (hg38_chrom_info$centromere_start +
         hg38_chrom_info$centromere_end) / 2
    ) +
    hg38_chrom_info$cum_length

  hg38_chrom_info$cum_chromosome_mid <-
    (hg38_chrom_info$length / 2) +
    hg38_chrom_info$cum_length

  hg38_chrom_info$p_start <- hg38_chrom_info$cum_length
  hg38_chrom_info$p_end <- hg38_chrom_info$cum_centromere_start
  hg38_chrom_info$q_start <- hg38_chrom_info$cum_centromere_end
  hg38_chrom_info$q_end <- cumsum(as.numeric(hg38_chrom_info$length))

  # Merge gene and chromosome information
  genomic_info <- merge(
    coords,
    hg38_chrom_info,
    by = "chrom"
  )

  # Calculate cumulative positions
  genomic_info$cum_start <-
    genomic_info$cum_length +
    genomic_info$start_position

  genomic_info$cum_end <-
    genomic_info$cum_length +
    genomic_info$end_position

  genomic_info$cum_mid <-
    genomic_info$cum_length +
    genomic_info$mid_position

  return(
    list(
      genomic_info = genomic_info,
      hg38_chrom_info = hg38_chrom_info
    )
  )
}
