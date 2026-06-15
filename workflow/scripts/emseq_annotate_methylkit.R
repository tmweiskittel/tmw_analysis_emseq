#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(methylKit)
  library(optparse)
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
})

opts <- parse_args(OptionParser(option_list = list(
  make_option("--db", type = "character"),
  make_option("--out", type = "character"),
  make_option("--gtf", type = "character")
)))

dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)

mdiff <- readMethylDB(opts$db)
df <- as.data.table(getData(mdiff))

message("Rows read from methylDiff: ", nrow(df))

if (nrow(df) == 0) {
  stop("No rows read from methylDiff DB: ", opts$db)
}

gtf <- import(opts$gtf)
gtf <- gtf[gtf$type %in% c("gene", "exon", "CDS", "five_prime_UTR", "three_prime_UTR")]

genes <- gtf[gtf$type == "gene"]
tss <- resize(genes, width = 1, fix = "start")

query <- GRanges(
  seqnames = df$chr,
  ranges = IRanges(start = df$start, end = df$end)
)

hits <- findOverlaps(query, gtf, ignore.strand = TRUE)

df[, feature_type := NA_character_]
df[, gene_id := NA_character_]
df[, gene_name := NA_character_]

if (length(hits) > 0) {
  q <- queryHits(hits)
  s <- subjectHits(hits)

  ann <- data.table(
    row_id = q,
    feature_type = as.character(gtf$type[s]),
    gene_id = as.character(mcols(gtf)$gene_id[s]),
    gene_name = as.character(mcols(gtf)$gene_name[s])
  )

  ann <- ann[
    ,
    .(
      feature_type = paste(unique(na.omit(feature_type)), collapse = ";"),
      gene_id = paste(unique(na.omit(gene_id)), collapse = ";"),
      gene_name = paste(unique(na.omit(gene_name)), collapse = ";")
    ),
    by = row_id
  ]

  df[ann$row_id, feature_type := ann$feature_type]
  df[ann$row_id, gene_id := ann$gene_id]
  df[ann$row_id, gene_name := ann$gene_name]
}

nearest <- distanceToNearest(query, tss, ignore.strand = TRUE)

df[, nearest_tss_gene_id := NA_character_]
df[, nearest_tss_gene_name := NA_character_]
df[, distance_to_tss := NA_integer_]

if (length(nearest) > 0) {
  q <- queryHits(nearest)
  s <- subjectHits(nearest)

  df[q, nearest_tss_gene_id := as.character(mcols(tss)$gene_id[s])]
  df[q, nearest_tss_gene_name := as.character(mcols(tss)$gene_name[s])]
  df[q, distance_to_tss := as.integer(mcols(nearest)$distance)]
}

fwrite(df, opts$out, sep = "\t", quote = FALSE, na = "NA")

message("Wrote annotated methylDiff file: ", opts$out)
