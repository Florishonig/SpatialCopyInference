#' @title Plot genomic clusters of differentially expressed genes
#'
#' @description
#' Visualises genomic clusters identified using
#' identify_gene_clusters().
#'
#' Differentially expressed genes are displayed on a
#' genome-wide karyotype and coloured according to the
#' size of the genomic cluster to which they belong.
#'
#' @param gene_clusters Output from
#' identify_gene_clusters()
#' @param driver_genes Optional vector of driver genes
#' to annotate
#' @param raised_driver_genes Optional subset of driver
#' genes displayed on a separate annotation track to
#' reduce label overlap
#'
#' @importFrom GenomicRanges seqnames
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors mcols mcols<-
#'
#' @returns An invisible list containing annotated gene
#' ranges, cluster ranges, driver gene ranges and the
#' karyoploteR object used for visualisation.
#'
#' @export
plot_gene_clusters <- function(
    gene_clusters,
    driver_genes = NULL,
    raised_driver_genes = NULL
) {

  gr_de <- gene_clusters$gr_de
  cluster_df <- gene_clusters$cluster_assignments

  # Create genomic ranges panel
  gr <- GenomicRanges::GRanges(
    seqnames = hg38_gene_info$chrom,
    ranges = IRanges::IRanges(
      start = hg38_gene_info$start_position,
      end = hg38_gene_info$end_position
    ),
    gene = hg38_gene_info$hgnc_symbol
  )

  # Optional driver genes
  if (!is.null(driver_genes)) {

    gr_drivers <- gr[
      gr$gene %in% driver_genes
    ]

  } else {

    gr_drivers <- GenomicRanges::GRanges()

  }

  rownames(cluster_df) <- NULL

  match_idx <- match(
    gr_de$gene,
    cluster_df$gene
  )

  gr_de$cluster <- cluster_df$cluster_global[match_idx]

  gr_de$cluster_size <- cluster_df$cluster_size_cat[match_idx]

  cols <- c(
    "2" = "#8dd3c7",
    "3" = "#80b1d3",
    "4" = "#bebada",
    "5+" = "#fb8072"
  )

  gr_de$col <- "grey80"

  has_cluster <- !is.na(gr_de$cluster_size)

  gr_de$col[has_cluster] <- cols[
    as.character(gr_de$cluster_size[has_cluster])
  ]

  cluster_list <- split(
    gr_de[!is.na(gr_de$cluster)],
    gr_de$cluster[!is.na(gr_de$cluster)]
  )

  cluster_gr <- unlist(
    GenomicRanges::GRangesList(
      lapply(cluster_list, function(gr) {

        cluster_size <- unique(
          gr$cluster_size
        )[1]

        GenomicRanges::GRanges(
          seqnames = seqnames(gr)[1],
          ranges = IRanges::IRanges(
            start = min(start(gr)),
            end = max(end(gr))
          ),
          size = cluster_size
        )
      })
    )
  )

  cluster_cols <- c(
    "2"  = scales::alpha("#8dd3c7", 0.4),
    "3"  = scales::alpha("#80b1d3", 0.4),
    "4"  = scales::alpha("#bebada", 0.4),
    "5+" = scales::alpha("#fb8072", 0.4)
  )

  mcols(cluster_gr)$col_rect <- cluster_cols[
    as.character(mcols(cluster_gr)$size)
  ]

  kp <- karyoploteR::plotKaryotype(
    genome = "hg38",
    chromosomes = paste0("chr", c(1:22, "X")),
    cex = 0.7
  )

  karyoploteR::kpRect(
    kp,
    data = cluster_gr,
    y0 = 0,
    y1 = 0.9,
    col = mcols(cluster_gr)$col_rect,
    border = NA
  )

  karyoploteR::kpPoints(
    kp,
    data = gr_de,
    col = gr_de$col,
    cex = 0.8,
    y = 0.75
  )

  if (length(gr_drivers) > 0) {

    karyoploteR::kpPoints(
      kp,
      data = gr_drivers,
      y = rep(-0.15, length(gr_drivers)),
      col = "#ffed6f",
      cex = 0.8,
      pch = 17
    )

    # Standard driver labels
    karyoploteR::kpPlotMarkers(
      kp,
      data = gr_drivers[
        !gr_drivers$gene %in% raised_driver_genes
      ],
      labels = gr_drivers$gene[
        !gr_drivers$gene %in% raised_driver_genes
      ],
      adjust.label.position = FALSE,
      text.orientation = "horizontal",
      line.color = scales::alpha("white", 0.01),
      label.color = "black",
      r0 = .1,
      r1 = .15,
      cex = 0.5
    )

    # Raised driver labels
    if (!is.null(raised_driver_genes) &&
        length(raised_driver_genes) > 0) {

      karyoploteR::kpPlotMarkers(
        kp,
        data = gr_drivers[
          gr_drivers$gene %in% raised_driver_genes
        ],
        labels = gr_drivers$gene[
          gr_drivers$gene %in% raised_driver_genes
        ],
        adjust.label.position = FALSE,
        text.orientation = "horizontal",
        line.color = scales::alpha("white", 0.01),
        label.color = "black",
        r0 = .1,
        r1 = .8,
        cex = 0.5
      )

    }
  }

  legend(
    "bottomright",
    legend = c(
      "Not clustered",
      "2 genes",
      "3 genes",
      "4 genes",
      "5+ genes",
      "Drivers"
    ),
    col = c(
      "grey80",
      "#8dd3c7",
      "#80b1d3",
      "#bebada",
      "#fb8072",
      "#ffed6f"
    ),
    pch = c(
      16, 16, 16, 16, 16, 17
    ),
    cex = 0.7,
    pt.cex = 0.7,
    bty = "n"
  )

  invisible(
    list(
      gr_de = gr_de,
      gr_drivers = gr_drivers,
      cluster_gr = cluster_gr,
      kp = kp
    )
  )
}
