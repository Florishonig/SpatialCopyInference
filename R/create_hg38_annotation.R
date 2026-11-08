#' @title Create hg38 annotation
#'
#' @description
#' Generates genome-wide gene and chromosome annotation
#' datasets for the human reference genome GRCh38 (hg38)
#' using Ensembl.
#'
#' This function retrieves all HGNC genes from Ensembl,
#' annotates their genomic positions and chromosome arm
#' assignments using `get_genomic_info()`, and saves the
#' resulting datasets as RDS files.
#'
#' This function was used to generate the packaged
#' datasets hg38_gene_info and hg38_chrom_info.
#'
#' @param gene_info_file Output file for gene annotation.
#' @param chrom_info_file Output file for chromosome annotation.
#'
#' @returns A list containing the generated hg38_gene_info and
#' hg38_chrom_info files.
#'
#' @export
create_hg38_annotation <- function(
    gene_info_file = "hg38_gene_info.rds",
    chrom_info_file = "hg38_chrom_info.rds"
) {

  # Connect to Ensembl
  mart <- biomaRt::useEnsembl(
    "ensembl",
    dataset = "hsapiens_gene_ensembl",
    mirror = "useast"
  )

  # Retrieve all hgnc symbols
  all_genes <- biomaRt::getBM(
    attributes = "hgnc_symbol",
    mart = mart
  )

  # Clean gene list
  all_genes <- all_genes %>%
    dplyr::filter(hgnc_symbol != "") %>%
    dplyr::distinct(hgnc_symbol)

  message(
    "Retrieved ",
    nrow(all_genes),
    " genes."
  )

  # Run annotation pipeline
  annotation <- get_genomic_info(
    genes = all_genes$hgnc_symbol
  )

  # Separate objects
  hg38_gene_info <- annotation$genomic_info
  hg38_chrom_info <- annotation$hg38_chrom_info

  # Save separately
  saveRDS(
    hg38_gene_info,
    file = gene_info_file
  )

  saveRDS(
    hg38_chrom_info,
    file = chrom_info_file
  )

  message(
    "Saved ",
    nrow(hg38_gene_info),
    " annotated genes to ",
    gene_info_file
  )

  message(
    "Saved chromosome information to ",
    chrom_info_file
  )

  return(
    list(
      hg38_gene_info = hg38_gene_info,
      hg38_chrom_info = hg38_chrom_info
    )
  )
}
