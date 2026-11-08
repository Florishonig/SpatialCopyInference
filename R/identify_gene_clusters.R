#' @title Identify gene clusters
#'
#' @description
#' Identifies genomic clusters of differentially
#' expressed genes.
#'
#' Differentially expressed genes are mapped to their
#' genomic coordinates and grouped into clusters based
#' on genomic proximity. Genes separated by less than
#' `gapwidth` base pairs are assigned to the same
#' cluster.
#'
#' @param de_results Differential expression results
#' returned by compare_spatial_clusters()
#' @param p_adj_cutoff Adjusted p-value threshold
#' @param logfc_cutoff Minimum log2 fold-change used to
#' define differentially expressed genes
#' @param gapwidth Maximum genomic distance between genes
#' for inclusion in the same cluster
#'
#' @importFrom GenomicRanges seqnames
#' @importFrom S4Vectors queryHits subjectHits
#'
#' @returns A list containing cluster assignments for
#' differentially expressed genes (`cluster_assignments`)
#' and the corresponding genomic ranges (`gr_de`).
#'
#' @export
identify_gene_clusters <- function(
    de_results,
    p_adj_cutoff = 0.01,
    logfc_cutoff = 0,
    gapwidth = 2e6
) {
  # Create genomic ranges for all panel genes
  gr <- GenomicRanges::GRanges(
    seqnames = hg38_gene_info$chrom,
    ranges = IRanges::IRanges(
      start = hg38_gene_info$start_position,
      end = hg38_gene_info$end_position
    ),
    gene = hg38_gene_info$hgnc_symbol
  )

  # Significant upregulated genes
  de_genes <- de_results %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(
      p_val_adj < p_adj_cutoff,
      avg_log2FC >= logfc_cutoff
    ) %>%
    dplyr::pull(gene)

  gr_de <- gr[gr$gene %in% de_genes]

  chroms <- paste0("chr", c(1:22, "X"))

  cluster_assignments <- list()

  for (chr in chroms) {

    gr_chr <- gr_de[seqnames(gr_de) == chr]

    if (length(gr_chr) == 0) {
      next
    }

    gr_chr <- sort(gr_chr)

    clusters <- GenomicRanges::reduce(
      gr_chr,
      min.gapwidth = gapwidth
    )

    hits <- GenomicRanges::findOverlaps(
      gr_chr,
      clusters
    )

    if (length(hits) > 0) {

      df <- data.frame(
        gene = gr_chr$gene[S4Vectors::queryHits(hits)],
        cluster_id = S4Vectors::subjectHits(hits)
      )

      cluster_sizes <- table(df$cluster_id)

      valid <- names(
        cluster_sizes[
          cluster_sizes >= 2
        ]
      )

      if (length(valid) > 0) {

        df <- df[
          df$cluster_id %in% valid,
        ]

        if (nrow(df) > 0) {

          df$gene_count <- cluster_sizes[
            as.character(df$cluster_id)
          ]

          df$cluster_size_cat <- ifelse(
            df$gene_count >= 5,
            "5+",
            as.character(df$gene_count)
          )

          df$cluster_size_cat <- factor(
            df$cluster_size_cat,
            levels = c(
              "2",
              "3",
              "4",
              "5+"
            )
          )

          df$cluster_global <- paste(
            chr,
            df$cluster_id,
            sep = "_"
          )

          cluster_assignments[[chr]] <- df
        }
      }
    }
  }

  cluster_assignments <- dplyr::bind_rows(
    cluster_assignments
  )

  list(
    cluster_assignments = cluster_assignments,
    gr_de = gr_de
  )
}
