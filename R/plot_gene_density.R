#' @title Plot genomic density of differentially expressed genes
#'
#' @description
#' Visualises the genomic density of differentially
#' expressed genes across the genome.
#'
#' Differentially expressed genes are mapped to their
#' genomic coordinates and summarised using a rolling
#' density score across chromosome arms. Results are
#' displayed as a genome-wide karyotype density map, with
#' optional annotation of known driver genes.
#'
#' @param de_results Differential expression results
#' returned by compare_spatial_clusters()
#' @param driver_genes Optional vector of driver genes
#' to annotate
#' @param raised_driver_genes Optional subset of driver
#' genes displayed above the karyotype to reduce label
#' overlap
#' @param p_adj_cutoff Adjusted p-value threshold
#' @param logfc_cutoff Minimum log2 fold-change used to
#' define differentially expressed genes
#' @param window_size Genomic window size used for
#' density estimation
#' @param step_size Step size used for density
#' estimation and smoothing
#'
#' @importFrom GenomicRanges seqnames
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors mcols mcols<-
#'
#' @returns An invisible list containing the smoothed
#' genomic signal, differentially expressed genes and
#' karyoploteR object used for visualisation.
#'
#' @examples
#' \dontrun{
#' # Custom driver genes
#'driver_genes <- c("MYC", "PIK3CA", "CCNE1", "KRAS",
#'                  "TERT","AKT1", "AKT2", "AKT3", "CCND1",
#'                  "CCND2", "CCND3", "KRAS", "MECOM")
#'
#' # Plot gene density with (raised) driver genes
#'plot_gene_density_obj <- plot_gene_density(
#'  de_results = cluster1v2$de,
#'  driver_genes = driver_genes,
#'  raised_driver_genes = c("MECOM","AKT2"))}
#'
#' @export
plot_gene_density <- function(
    de_results,
    driver_genes = NULL,
    raised_driver_genes = NULL,
    p_adj_cutoff = 0.01,
    logfc_cutoff = 0,
    window_size = 5e6,
    step_size = 5e5
) {

  # panel genes
  gr <- GenomicRanges::GRanges(
    seqnames = hg38_gene_info$chrom,
    ranges = IRanges::IRanges(
      start = hg38_gene_info$start_position,
      end   = hg38_gene_info$end_position
    ),
    gene = hg38_gene_info$hgnc_symbol
  )

  # DE genes
  de_genes <- de_results %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(
      p_val_adj < p_adj_cutoff,
      avg_log2FC >= logfc_cutoff
    ) %>%
    dplyr::pull(gene)

  gr_de <- gr[gr$gene %in% de_genes]

  # rolling density calculation
  chroms <- paste0("chr", c(1:22, "X"))

  results <- list()

  for(chr in chroms){

    chr_info <- hg38_chrom_info[
      hg38_chrom_info$chrom == chr,
    ]

    if(nrow(chr_info) == 0) next

    chr_len <- chr_info$length

    p_end <- chr_info$centromere_start
    q_start <- chr_info$centromere_end

    windows_list <- list()

    # p arm
    if(p_end > 0){

      if(p_end <= window_size){

        windows_p <- GenomicRanges::GRanges(
          seqnames = chr,
          ranges = IRanges::IRanges(
            start = 1,
            end = p_end
          )
        )

      } else {

        windows_p <- GenomicRanges::GRanges(
          seqnames = chr,
          ranges = IRanges::IRanges(
            start = seq(
              1,
              p_end - window_size,
              by = step_size
            ),
            width = window_size
          )
        )
      }

      windows_list[[length(windows_list)+1]] <- windows_p
    }

    # q arm
    q_length <- chr_len - q_start

    if(q_length > 0){

      if(q_length <= window_size){

        windows_q <- GenomicRanges::GRanges(
          seqnames = chr,
          ranges = IRanges::IRanges(
            start = q_start,
            end = chr_len
          )
        )

      } else {

        windows_q <- GenomicRanges::GRanges(
          seqnames = chr,
          ranges = IRanges::IRanges(
            start = seq(
              q_start,
              chr_len - window_size,
              by = step_size
            ),
            width = window_size
          )
        )
      }

      windows_list[[length(windows_list)+1]] <- windows_q
    }

    windows <- Reduce(c, windows_list)

    panel_chr <- gr[seqnames(gr) == chr]
    de_chr <- gr_de[seqnames(gr_de) == chr]

    panel_counts <- GenomicRanges::countOverlaps(
      windows,
      panel_chr
    )

    de_counts <- GenomicRanges::countOverlaps(
      windows,
      de_chr
    )

    score <- ifelse(
      panel_counts > 0,
      de_counts / sqrt(panel_counts),
      NA_real_
    )

    results[[chr]] <- data.frame(
      chrom = chr,
      start = start(windows),
      end = end(windows),
      score = score
    )
  }

  density_df <- dplyr::bind_rows(results)

  density_df$score_scaled <- scales::rescale(
    density_df$score,
    to = c(0,1),
    na.rm = TRUE
  )

  # smoothing
  window_gr <- GenomicRanges::GRanges(
    seqnames = density_df$chrom,
    ranges = IRanges::IRanges(
      start = density_df$start,
      end = density_df$end
    ),
    score = density_df$score_scaled
  )

  signal_gr <- GenomicRanges::GRanges()

  for(chr in unique(density_df$chrom)) {

    chr_len <- hg38_chrom_info$length[
      hg38_chrom_info$chrom == chr
    ]

    bins <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges = IRanges::IRanges(
        start = seq(
          1,
          chr_len,
          by = step_size
        ),
        width = step_size
      )
    )

    win_chr <- window_gr[
      seqnames(window_gr) == chr
    ]

    hits <- GenomicRanges::findOverlaps(
      bins,
      win_chr
    )

    score <- rep(
      NA_real_,
      length(bins)
    )

    split_hits <- split(
        subjectHits(hits),
        queryHits(hits)
    )

    for(i in names(split_hits)) {

      idx <- as.integer(i)

      score[idx] <- mean(
        win_chr$score[
          split_hits[[i]]
        ],
        na.rm = TRUE
      )
    }

    bins$score <- score

    signal_gr <- c(signal_gr, bins)
  }

  # colours
  cols <- colorRampPalette(
    c(
      "#f7fbff",
      "#deebf7",
      "#c6dbef",
      "#9ecae1",
      "#6baed6",
      "#4292c6",
      "#2171b5",
      "#08519c",
      "#08306b"
    )
  )(100)

  signal_gr$col <- cols[
    cut(
      signal_gr$score,
      breaks = 100,
      include.lowest = TRUE,
      labels = FALSE
    )
  ]

  # driver genes
  if(!is.null(driver_genes)) {

    gr_drivers <- gr[
      gr$gene %in% driver_genes
    ]

  } else {

    gr_drivers <- GenomicRanges::GRanges()
  }

  kp <- karyoploteR::plotKaryotype(
    genome = "hg38",
    chromosomes = paste0("chr", c(1:22,"X")),
    cex = 0.7
  )

  karyoploteR::kpRect(
    kp,
    data = signal_gr,
    y0 = 0,
    y1 = 0.6,
    col = scales::alpha(
      signal_gr$col,
      0.85
    ),
    border = NA
  )

  karyoploteR::kpPoints(
    kp,
    data = gr,
    y = rep(.75,length(gr)),
    col = "grey80",
    cex = 0.8
  )

  karyoploteR::kpPoints(
    kp,
    data = gr_de,
    y = rep(.75,length(gr_de)),
    col = "#E69F00",
    cex = 0.8
  )

  if(length(gr_drivers) > 0) {

    if(is.null(raised_driver_genes)) {
      raised_driver_genes <- character(0)
    }

    # driver genes
    mcols(gr_drivers)$label <- gr_drivers$gene

    # Standard driver labels
    karyoploteR::kpPlotMarkers(
      kp,
      data = gr_drivers[
        !gr_drivers$gene %in% raised_driver_genes
      ],
      adjust.label.position = FALSE,
      labels = gr_drivers$gene[
        !gr_drivers$gene %in% raised_driver_genes
      ],
      text.orientation = "horizontal",
      line.color = scales::alpha("white", 0.01),
      label.color = "black",
      r0 = 0.1,
      r1 = 0.15,
      cex = 0.5
    )

    # Driver gene positions
    karyoploteR::kpPoints(
      kp,
      data = gr_drivers,
      y = rep(-0.15, length(gr_drivers)),
      col = "#ffed6f",
      cex = 0.8,
      pch = 17
    )

    # Raised driver labels
    if(length(raised_driver_genes) > 0) {

      karyoploteR::kpPlotMarkers(
        kp,
        data = gr_drivers[
          gr_drivers$gene %in% raised_driver_genes
        ],
        adjust.label.position = FALSE,
        labels = gr_drivers$gene[
          gr_drivers$gene %in% raised_driver_genes
        ],
        text.orientation = "horizontal",
        line.color = scales::alpha("white", 0.01),
        label.color = "black",
        r0 = 1,
        r1 = 0.5,
        cex = 0.5
      )
    }
  }

  # heatmap scale
  old_xpd <- par("xpd")
  on.exit(par(xpd = old_xpd))
  par(xpd = TRUE)

  legend(
    "bottomright",
    legend = c(
      "Xenium panel",
      paste0(
        "DE genes (n=",
        length(gr_de),
        ")"
      ),
      "Driver genes"
    ),
    col = c(
      "grey80",
      "#E69F00",
      "#ffed6f"
    ),
    pch = c(16,16,17),
    cex = 0.7,
    pt.cex = 0.7,
    bty = "n"
  )

  scale_x0 <- 0.76
  scale_x1 <- 0.79

  scale_y0 <- 0.16
  scale_y1 <- 0.31

  n <- length(cols)

  for(i in seq_len(n)) {

    rect(
      scale_x0,
      scale_y0 + (i-1)/n * (scale_y1-scale_y0),
      scale_x1,
      scale_y0 + i/n * (scale_y1-scale_y0),
      col = cols[i],
      border = NA
    )
  }

  rect(
    scale_x0,
    scale_y0,
    scale_x1,
    scale_y1,
    border = "black",
    lwd = 0.5
  )

  tick_vals <- c(0, 0.25, 0.5, 0.75, 1)

  tick_y <- scale_y0 +
    tick_vals * (scale_y1 - scale_y0)

  segments(
    x0 = scale_x1,
    y0 = tick_y,
    x1 = scale_x1 + 0.004,
    y1 = tick_y,
    xpd = TRUE,
    lwd = 0.5
  )

  text(
    x = scale_x1 + 0.008,
    y = tick_y,
    labels = tick_vals,
    adj = 0,
    cex = 0.7,
    xpd = TRUE
  )

  text(
    scale_x0 - 0.01,
    scale_y1  + 0.03,
    labels = "Density score",
    adj = 0,
    cex = 0.7
  )

  invisible(
    list(
      signal_gr = signal_gr,
      gr_de = gr_de,
      kp = kp
    )
  )
}
