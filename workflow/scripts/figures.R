library(data.table)
library(ggplot2)
library(pheatmap)
library(uwot)

experiment <- "YOUR_EXPERIMENT"
D_OUT <- "PATH_TO_OUTPUT"

FIG_DIR <- file.path(D_OUT, "figures", experiment)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

mdiff <- fread(file.path(D_OUT, "dmr/diff",
                         paste0("methylDiff_", experiment, ".txt.bgz")))

tiled_mdiff <- fread(file.path(D_OUT, "dmr/diff",
                               paste0("methylDiff_", experiment, ".tiled.txt.bgz")))

meth_mat <- fread(file.path(D_OUT, "dmr/diff",
                            paste0(experiment, "_pos_meth.tsv")))

annot <- fread(file.path(D_OUT, "dmr/annotation",
                         paste0(experiment, "_annotated.tsv")))

p <- ggplot(mdiff, aes(x = meth.diff)) +
    geom_histogram(bins = 100) +
    theme_bw() +
    xlab("Methylation difference (%)") +
    ylab("CpG count")

ggsave(
    filename = file.path(FIG_DIR, "methylation_difference_histogram.png"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
)

p <- ggplot(mdiff,
            aes(x = meth.diff,
                y = -log10(qvalue))) +
    geom_point(alpha = 0.35) +
    theme_bw() +
    xlab("Methylation difference (%)") +
    ylab("-log10(q-value)")

ggsave(
    file.path(FIG_DIR, "volcano_plot.png"),
    p,
    width = 8,
    height = 6,
    dpi = 300
)

mdiff[, chr := factor(chr, levels = unique(chr))]

p <- ggplot(mdiff,
            aes(x = start,
                y = -log10(qvalue))) +
    geom_point(alpha = 0.35,
               size = 0.4) +
    facet_wrap(~chr,
               scales = "free_x") +
    theme_bw()

ggsave(
    file.path(FIG_DIR, "manhattan_plot.png"),
    p,
    width = 14,
    height = 8,
    dpi = 300
)

num_cols <- names(meth_mat)[sapply(meth_mat, is.numeric)]

mat <- as.matrix(meth_mat[, ..num_cols])
mat <- t(scale(t(mat)))

pca <- prcomp(t(mat))

pca_df <- data.frame(
    sample = colnames(mat),
    PC1 = pca$x[,1],
    PC2 = pca$x[,2]
)

p <- ggplot(pca_df,
            aes(PC1, PC2,
                label = sample)) +
    geom_point(size = 3) +
    theme_bw()

ggsave(
    file.path(FIG_DIR, "PCA.png"),
    p,
    width = 7,
    height = 6,
    dpi = 300
)

set.seed(1)

embedding <- uwot::umap(t(mat))

umap_df <- data.frame(
    sample = colnames(mat),
    UMAP1 = embedding[,1],
    UMAP2 = embedding[,2]
)

p <- ggplot(umap_df,
            aes(UMAP1,
                UMAP2,
                label = sample)) +
    geom_point(size = 3) +
    theme_bw()

ggsave(
    file.path(FIG_DIR, "UMAP.png"),
    p,
    width = 7,
    height = 6,
    dpi = 300
)

top <- mdiff[order(qvalue)][1:min(500, .N)]

top_key <- paste(top$chr,
                 top$start,
                 sep = "_")

meth_mat[, key := paste(chr,
                        start,
                        sep = "_")]

heat <- meth_mat[key %in% top_key,
                 ..num_cols]

heat <- as.matrix(heat)
heat <- t(scale(t(heat)))

png(
    filename = file.path(FIG_DIR, "Top500_heatmap.png"),
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

if ("annotation" %in% names(annot)) {

    p <- ggplot(annot,
                aes(annotation)) +
        geom_bar() +
        coord_flip() +
        theme_bw()

    ggsave(
        file.path(FIG_DIR, "annotation_distribution.png"),
        p,
        width = 8,
        height = 6,
        dpi = 300
    )
}

summary_stats <- data.frame(
    Total_CpGs = nrow(mdiff),
    Significant = sum(mdiff$qvalue < 0.05),
    Hypermethylated = sum(mdiff$qvalue < 0.05 &
                          mdiff$meth.diff > 10),
    Hypomethylated = sum(mdiff$qvalue < 0.05 &
                         mdiff$meth.diff < -10)
)

write.csv(
    summary_stats,
    file.path(FIG_DIR, "summary_statistics.csv"),
    row.names = FALSE
)

fwrite(
    mdiff[qvalue < 0.05],
    file.path(FIG_DIR, "significant_DMCs.tsv"),
    sep = "\t"
)

fwrite(
    mdiff[qvalue < 0.05 & meth.diff > 10],
    file.path(FIG_DIR, "hypermethylated_DMCs.tsv"),
    sep = "\t"
)

fwrite(
    mdiff[qvalue < 0.05 & meth.diff < -10],
    file.path(FIG_DIR, "hypomethylated_DMCs.tsv"),
    sep = "\t"
)
