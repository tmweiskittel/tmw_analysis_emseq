#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(optparse)
})

option_list <- list(
  make_option(
  "--sample-contrasts",
  type = "character",
  dest = "sample_contrasts",
  help = "Sample contrast CSV containing sample IDs and cohort assignments"
  ),
  make_option("--mdiff", type = "character", dest = "mdiff"),
  make_option("--tiled-mdiff", type = "character", dest = "tiled_mdiff"),
  make_option("--matrix", type = "character", dest = "matrix"),
  make_option("--annotation", type = "character", dest = "annotation"),
  make_option("--outdir", type = "character", dest = "outdir"),
  make_option("--contrast", type = "character", dest = "contrast"),
  make_option("--qvalue-cutoff", type = "double", dest = "qvalue_cutoff", default = 0.05),
  make_option("--meth-diff-cutoff", type = "double", dest = "meth_diff_cutoff", default = 10),
  make_option("--top-n-heatmap", type = "integer", dest = "top_n_heatmap", default = 500)
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

read_methylkit_bgz <- function(path) {
  x <- fread(
    cmd = paste("zcat", shQuote(path), "| grep -v '^#'"),
    sep = "\t",
    header = TRUE
  )

  if (!all(c("chr", "start", "end", "strand", "pvalue", "qvalue", "meth.diff") %in% names(x))) {
    if (ncol(x) == 7) {
      setnames(x, c("chr", "start", "end", "strand", "pvalue", "qvalue", "meth.diff"))
    } else {
      stop("Unexpected methylDiff format. Columns detected: ", paste(names(x), collapse = ", "))
    }
  }

  x
}

mdiff <- read_methylkit_bgz(opt$mdiff)
tiled_mdiff <- read_methylkit_bgz(opt$tiled_mdiff)
meth_mat <- fread(opt$matrix)
annot <- fread(opt$annotation)

required_cols <- c("chr", "start", "qvalue", "meth.diff")
missing_cols <- setdiff(required_cols, names(mdiff))
if (length(missing_cols) > 0) {
  stop("Missing required methylDiff columns: ", paste(missing_cols, collapse = ", "))
}

sig <- mdiff[qvalue < opt$qvalue_cutoff & abs(meth.diff) >= opt$meth_diff_cutoff]
hyper <- sig[meth.diff > 0]
hypo <- sig[meth.diff < 0]

p <- ggplot(mdiff, aes(x = meth.diff)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  xlab("Methylation difference (%)") +
  ylab("CpG count")

ggsave(
  file.path(opt$outdir, paste0(opt$contrast, "_methylation_difference_histogram.png")),
  p,
  width = 8,
  height = 6,
  dpi = 300
)

p <- ggplot(mdiff, aes(x = meth.diff, y = -log10(qvalue))) +
  geom_point(alpha = 0.35, size = 0.6) +
  theme_bw() +
  xlab("Methylation difference (%)") +
  ylab("-log10(q-value)")

ggsave(
  file.path(opt$outdir, paste0(opt$contrast, "_volcano_plot.png")),
  p,
  width = 8,
  height = 6,
  dpi = 300
)

mdiff[, chr := factor(chr, levels = unique(chr))]

p <- ggplot(mdiff, aes(x = start, y = -log10(qvalue))) +
  geom_point(alpha = 0.35, size = 0.4) +
  facet_wrap(~ chr, scales = "free_x") +
  theme_bw() +
  xlab("Genomic position") +
  ylab("-log10(q-value)")

ggsave(
  file.path(opt$outdir, paste0(opt$contrast, "_manhattan_plot.png")),
  p,
  width = 14,
  height = 8,
  dpi = 300
)

if (all(c("chr", "start", "meth.diff") %in% names(tiled_mdiff))) {
  tiled_mdiff[, chr := factor(chr, levels = unique(chr))]

  p <- ggplot(tiled_mdiff, aes(x = start, y = meth.diff)) +
    geom_point(alpha = 0.35, size = 0.4) +
    facet_wrap(~ chr, scales = "free_x") +
    theme_bw() +
    xlab("Genomic position") +
    ylab("Tiled methylation difference (%)")

  ggsave(
    file.path(opt$outdir, paste0(opt$contrast, "_tiled_methylation_difference.png")),
    p,
    width = 14,
    height = 8,
    dpi = 300
  )
}

# -------------------------------------------------------------------------
# PCA and heatmap
# -------------------------------------------------------------------------

coordinate_cols <- c("chr", "start", "end", "strand")

sample_cols <- setdiff(
  names(meth_mat),
  coordinate_cols
)

# Retain only numeric sample columns.
sample_cols <- sample_cols[
  vapply(
    meth_mat[, ..sample_cols],
    is.numeric,
    logical(1)
  )
]

pca_file <- file.path(
  opt$outdir,
  paste0(opt$contrast, "_PCA.png")
)

if (length(sample_cols) >= 2L) {
  message(
    "Preparing methylation matrix with ",
    length(sample_cols),
    " sample columns"
  )

  mat <- as.matrix(
    meth_mat[, ..sample_cols]
  )

  storage.mode(mat) <- "double"

  rownames(mat) <- paste(
    meth_mat$chr,
    meth_mat$start,
    sep = "_"
  )

  # Convert non-finite values to NA.
  mat[!is.finite(mat)] <- NA_real_

  # Require methylation values for at least half of the samples.
  minimum_observed <- max(
    2L,
    ceiling(ncol(mat) * 0.5)
  )

  observed_per_row <- rowSums(!is.na(mat))

  mat <- mat[
    observed_per_row >= minimum_observed,
    ,
    drop = FALSE
  ]

  message(
    "Rows retained after missingness filter: ",
    format(nrow(mat), big.mark = ",")
  )

  if (nrow(mat) >= 2L) {
    # Impute missing values using each CpG's mean methylation.
    row_means <- rowMeans(
      mat,
      na.rm = TRUE
    )

    missing_index <- which(
      is.na(mat),
      arr.ind = TRUE
    )

    if (nrow(missing_index) > 0L) {
      mat[missing_index] <- row_means[
        missing_index[, "row"]
      ]
    }

    # Remove zero-variance and non-finite rows.
    row_sd <- apply(
      mat,
      1L,
      sd
    )

    keep <- is.finite(row_sd) & row_sd > 0

    mat <- mat[
      keep,
      ,
      drop = FALSE
    ]

    row_sd <- row_sd[keep]

    message(
      "Variable rows retained: ",
      format(nrow(mat), big.mark = ",")
    )

    if (nrow(mat) >= 2L) {
      # PCA does not need millions of CpGs. Use the most variable loci.
      max_pca_rows <- 50000L

      if (nrow(mat) > max_pca_rows) {
        variance_order <- order(
          row_sd,
          decreasing = TRUE
        )

        mat_pca <- mat[
          variance_order[seq_len(max_pca_rows)],
          ,
          drop = FALSE
        ]
      } else {
        mat_pca <- mat
      }

      message(
        "CpGs used for PCA: ",
        format(nrow(mat_pca), big.mark = ",")
      )

      # Center each CpG across samples.
      mat_pca_scaled <- t(
        scale(
          t(mat_pca),
          center = TRUE,
          scale = FALSE
        )
      )

      # Remove any rows that became non-finite.
      finite_rows <- apply(
        mat_pca_scaled,
        1L,
        function(x) all(is.finite(x))
      )

      mat_pca_scaled <- mat_pca_scaled[
        finite_rows,
        ,
        drop = FALSE
      ]

      if (nrow(mat_pca_scaled) >= 2L) {
        pca <- prcomp(
          t(mat_pca_scaled),
          center = FALSE,
          scale. = FALSE
        )

        variance_explained <- (
          pca$sdev^2 /
            sum(pca$sdev^2)
        ) * 100

        pca_df <- data.frame(
          sample = rownames(pca$x),
          PC1 = pca$x[, 1L],
          PC2 = pca$x[, 2L],
          stringsAsFactors = FALSE
        )

        p <- ggplot(
          pca_df,
          aes(
            x = PC1,
            y = PC2,
            label = sample
          )
        ) +
          geom_point(size = 3) +
          geom_text(
            vjust = -0.8,
            size = 3,
            check_overlap = TRUE
          ) +
          theme_bw() +
          xlab(
            sprintf(
              "PC1 (%.1f%%)",
              variance_explained[1L]
            )
          ) +
          ylab(
            sprintf(
              "PC2 (%.1f%%)",
              variance_explained[2L]
            )
          ) +
          ggtitle(
            paste(
              opt$contrast,
              "methylation PCA"
            )
          )

        ggsave(
          filename = pca_file,
          plot = p,
          width = 8,
          height = 7,
          dpi = 300
        )
      }

      # ---------------------------------------------------------------
      # Heatmap
      # ---------------------------------------------------------------

      top <- mdiff[
        is.finite(qvalue)
      ][
        order(qvalue)
      ][
        seq_len(
          min(
            opt$top_n_heatmap,
            .N
          )
        )
      ]

      top_key <- paste(
        as.character(top$chr),
        top$start,
        sep = "_"
      )

      heat <- mat[
        rownames(mat) %in% top_key,
        ,
        drop = FALSE
      ]

      message(
        "Heatmap candidate rows: ",
        nrow(heat)
      )

      if (nrow(heat) >= 2L) {

        # Remove rows with non-finite values before scaling.
        finite_before_scaling <- apply(
          heat,
          1L,
          function(x) all(is.finite(x))
        )

        heat <- heat[
          finite_before_scaling,
          ,
          drop = FALSE
        ]

        # Remove zero-variance rows before scaling.
        if (nrow(heat) >= 2L) {
          heat_row_sd <- apply(
            heat,
            1L,
            sd
          )

          keep_heat_rows <- is.finite(heat_row_sd) &
            heat_row_sd > 0

          heat <- heat[
            keep_heat_rows,
            ,
            drop = FALSE
          ]
        }

        # Row-scale CpGs.
        if (nrow(heat) >= 2L) {
          heat <- t(
            scale(
              t(heat),
              center = TRUE,
              scale = TRUE
            )
          )

          # Safety check after scaling.
          finite_after_scaling <- apply(
            heat,
            1L,
            function(x) all(is.finite(x))
          )

          heat <- heat[
            finite_after_scaling,
            ,
            drop = FALSE
          ]
        }

        # Remove sample columns that cannot be clustered.
        if (nrow(heat) >= 2L) {
          col_sd <- apply(
            heat,
            2L,
            sd,
            na.rm = TRUE
          )

          keep_cols <- is.finite(col_sd) &
            col_sd > 0

          heat <- heat[
            ,
            keep_cols,
            drop = FALSE
          ]
        }

        message(
          "Heatmap rows retained: ",
          nrow(heat),
          "; sample columns retained: ",
          ncol(heat)
        )

        if (nrow(heat) >= 2L && ncol(heat) >= 2L) {

          heatmap_file <- file.path(
            opt$outdir,
            paste0(
              opt$contrast,
              "_top",
              opt$top_n_heatmap,
              "_heatmap.png"
            )
          )

          png(
            filename = heatmap_file,
            width = 2400,
            height = 3000,
            res = 300
          )

          pheatmap(
            heat,
            show_rownames = FALSE,
            clustering_distance_cols = "euclidean",
            clustering_distance_rows = "euclidean",
            main = paste(
              opt$contrast,
              "top differential methylation loci"
            )
          )

          dev.off()
        }
      }
    }
  }
}

# -------------------------------------------------------------------------
# PCA fallback
# -------------------------------------------------------------------------

if (!file.exists(pca_file)) {

  png(
    filename = pca_file,
    width = 2100,
    height = 1800,
    res = 300
  )

  plot.new()

  title("PCA not generated")

  text(
    0.5,
    0.55,
    "The methylation matrix did not contain enough\nfinite, variable CpGs for PCA.",
    cex = 1.1
  )

  text(
    0.5,
    0.40,
    paste(
      "Detected sample columns:",
      length(sample_cols)
    ),
    cex = 0.9
  )

  dev.off()
}
      
annotation_col <- intersect(
  c("annotation", "feature", "gene_annotation", "annot.type"),
  names(annot)
)[1]

if (!is.na(annotation_col)) {
  p <- ggplot(annot, aes(x = .data[[annotation_col]])) +
    geom_bar() +
    coord_flip() +
    theme_bw() +
    xlab("Annotation") +
    ylab("Count")

  ggsave(
    file.path(opt$outdir, paste0(opt$contrast, "_annotation_distribution.png")),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

summary_stats <- data.frame(
  contrast = opt$contrast,
  total_tested = nrow(mdiff),
  significant = nrow(sig),
  hypermethylated = nrow(hyper),
  hypomethylated = nrow(hypo),
  qvalue_cutoff = opt$qvalue_cutoff,
  meth_diff_cutoff = opt$meth_diff_cutoff
)

fwrite(summary_stats, file.path(opt$outdir, "summary_statistics.csv"))
fwrite(sig, file.path(opt$outdir, "significant_DMCs.tsv"), sep = "\t")
fwrite(hyper, file.path(opt$outdir, "hypermethylated_DMCs.tsv"), sep = "\t")
fwrite(hypo, file.path(opt$outdir, "hypomethylated_DMCs.tsv"), sep = "\t")
