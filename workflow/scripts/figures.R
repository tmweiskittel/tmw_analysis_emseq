#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(optparse)
})

option_list <- list(
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

num_cols <- names(meth_mat)[sapply(meth_mat, is.numeric)]
coord_cols <- intersect(c("chr", "start", "end"), names(meth_mat))
sample_cols <- setdiff(num_cols, coord_cols)

if (length(sample_cols) >= 2) {
  mat <- as.matrix(meth_mat[, ..sample_cols])
  rownames(mat) <- paste(meth_mat$chr, meth_mat$start, sep = "_")

  row_sd <- apply(mat, 1, sd, na.rm = TRUE)
  keep <- is.finite(row_sd) & row_sd > 0
  mat <- mat[keep, , drop = FALSE]

  mat_scaled <- t(scale(t(mat)))

  pca <- prcomp(t(mat_scaled), center = TRUE, scale. = FALSE)

  pca_df <- data.frame(
    sample = colnames(mat_scaled),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, label = sample)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.8, size = 3) +
    theme_bw() +
    xlab("PC1") +
    ylab("PC2")

  ggsave(
    file.path(opt$outdir, paste0(opt$contrast, "_PCA.png")),
    p,
    width = 7,
    height = 6,
    dpi = 300
  )

  top <- mdiff[order(qvalue)][1:min(opt$top_n_heatmap, .N)]
  top_key <- paste(top$chr, top$start, sep = "_")

  heat <- mat_scaled[rownames(mat_scaled) %in% top_key, , drop = FALSE]

  if (nrow(heat) >= 2) {
    png(
      filename = file.path(opt$outdir, paste0(opt$contrast, "_top", opt$top_n_heatmap, "_heatmap.png")),
      width = 2400,
      height = 3000,
      res = 300
    )

    pheatmap(
      heat,
      show_rownames = FALSE,
      clustering_distance_cols = "correlation",
      clustering_distance_rows = "euclidean"
    )

    dev.off()
  }
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
