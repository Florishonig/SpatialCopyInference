#' hg38 gene information
#'
#' Gene level information for the human reference genome
#' GRCh38 (hg38).
#'
#' The dataset was generated from Ensembl gene annotations
#' using the SpatialCopyInference functions
#' `get_genomic_info()` and `create_hg38_annotation()`.
#'
#' It contains genomic coordinates, chromosome
#' assignments and chromosome arm annotations for all
#' HGNC genes annotated in Ensembl.
#'
#' @format A data frame with 41,480 rows and 24 columns
#' containing Ensembl gene identifiers, HGNC symbols,
#' genomic coordinates, chromosome assignments and
#' chromosome arm annotations.
#'
#' @source Ensembl
"hg38_gene_info"


#' hg38 chromosome information
#'
#' Chromosome and chromosome arm information for the
#' human reference genome GRCh38 (hg38).
#'
#' The dataset was generated from Ensembl chromosome and
#' centromere coordinates using the SpatialCopyInference
#' functions `get_genomic_info()` and
#' `create_hg38_annotation()`.
#'
#' @format A data frame with 23 rows and 13 columns
#' containing chromosome lengths, centromere positions,
#' chromosome arm boundaries and cumulative genomic
#' coordinates.
#'
#' @source Ensembl
"hg38_chrom_info"
