#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
})

option_list <- list(
  make_option("--annotation", type = "character", dest = "annotation"),
  make_option("--outdir", type = "character", dest = "outdir"),
  make_option("--contrast", type = "character", dest = "contrast"),
  make_option("--qvalue-cutoff", type = "double", dest = "qvalue_cutoff", default = 0.05),
  make_option("--meth-diff-cutoff", type = "double", dest = "meth_diff_cutoff", default = 10),
  make_option("--show-category", type = "integer", dest = "show_category", default = 20)
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

annot <- fread(opt$annotation)

find_col <- function(x, candidates) {
  found <- intersect(candidates, names(x))
  if (length(found) == 0) return(NA_character_)
  found[1]
}

gene_col <- find_col(
  annot,
  c("gene_name", "gene", "geneSymbol", "symbol", "SYMBOL", "external_gene_name", "gene_id")
)

if (is.na(gene_col)) {
  stop("No gene column found. Available columns: ", paste(names(annot), collapse = ", "))
}

if (!"qvalue" %in% names(annot)) {
  stop("Missing qvalue column in annotation file.")
}

if (!"meth.diff" %in% names(annot)) {
  stop("Missing meth.diff column in annotation file.")
}

clean_genes <- function(x) {
  x <- unique(na.omit(as.character(x)))
  x <- x[x != "" & x != "NA" & x != "."]
  x
}

run_enrichment <- function(genes, label) {
  genes <- clean_genes(genes)

  gene_table <- data.table(
    contrast = opt$contrast,
    label = label,
    gene = genes
  )

  fwrite(
    gene_table,
    file.path(opt$outdir, paste0(opt$contrast, "_", label, "_genes.tsv")),
    sep = "\t"
  )

  if (length(genes) < 5) {
    warning("Skipping ", label, ": fewer than 5 genes.")
    return(invisible(NULL))
  }

  mapped <- suppressMessages(
    bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )

  if (is.null(mapped) || nrow(mapped) < 5) {
    warning("Skipping ", label, ": fewer than 5 mapped Entrez IDs.")
    return(invisible(NULL))
  }

  entrez <- unique(mapped$ENTREZID)

  go_bp <- enrichGO(
    gene = entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  go_mf <- enrichGO(
    gene = entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  go_cc <- enrichGO(
    gene = entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  reactome <- enrichPathway(
    gene = entrez,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  fwrite(as.data.table(mapped), file.path(opt$outdir, paste0(opt$contrast, "_", label, "_gene_mapping.tsv")), sep = "\t")
  fwrite(as.data.table(go_bp), file.path(opt$outdir, paste0(opt$contrast, "_", label, "_GO_BP.tsv")), sep = "\t")
  fwrite(as.data.table(go_mf), file.path(opt$outdir, paste0(opt$contrast, "_", label, "_GO_MF.tsv")), sep = "\t")
  fwrite(as.data.table(go_cc), file.path(opt$outdir, paste0(opt$contrast, "_", label, "_GO_CC.tsv")), sep = "\t")
  fwrite(as.data.table(reactome), file.path(opt$outdir, paste0(opt$contrast, "_", label, "_Reactome.tsv")), sep = "\t")

  plot_enrichment <- function(obj, suffix, title) {
    df <- as.data.frame(obj)
    if (nrow(df) == 0) return(invisible(NULL))

    p <- dotplot(obj, showCategory = opt$show_category) +
      ggtitle(title)

    ggsave(
      file.path(opt$outdir, paste0(opt$contrast, "_", label, "_", suffix, "_dotplot.png")),
      p,
      width = 10,
      height = 7,
      dpi = 300
    )

    p2 <- barplot(obj, showCategory = opt$show_category) +
      ggtitle(title)

    ggsave(
      file.path(opt$outdir, paste0(opt$contrast, "_", label, "_", suffix, "_barplot.png")),
      p2,
      width = 10,
      height = 7,
      dpi = 300
    )
  }

  plot_enrichment(go_bp, "GO_BP", paste(opt$contrast, label, "GO Biological Process"))
  plot_enrichment(go_mf, "GO_MF", paste(opt$contrast, label, "GO Molecular Function"))
  plot_enrichment(go_cc, "GO_CC", paste(opt$contrast, label, "GO Cellular Component"))
  plot_enrichment(reactome, "Reactome", paste(opt$contrast, label, "Reactome"))
}

sig <- annot[qvalue < opt$qvalue_cutoff]
hyper <- sig[meth.diff >= opt$meth_diff_cutoff]
hypo <- sig[meth.diff <= -opt$meth_diff_cutoff]

run_enrichment(sig[[gene_col]], "all_significant")
run_enrichment(hyper[[gene_col]], "hypermethylated")
run_enrichment(hypo[[gene_col]], "hypomethylated")

summary <- data.table(
  contrast = opt$contrast,
  gene_column = gene_col,
  total_annotated_rows = nrow(annot),
  significant_rows = nrow(sig),
  hypermethylated_rows = nrow(hyper),
  hypomethylated_rows = nrow(hypo),
  significant_genes = length(clean_genes(sig[[gene_col]])),
  hypermethylated_genes = length(clean_genes(hyper[[gene_col]])),
  hypomethylated_genes = length(clean_genes(hypo[[gene_col]])),
  qvalue_cutoff = opt$qvalue_cutoff,
  meth_diff_cutoff = opt$meth_diff_cutoff
)

fwrite(
  summary,
  file.path(opt$outdir, paste0(opt$contrast, "_pathway_summary.tsv")),
  sep = "\t"
)
