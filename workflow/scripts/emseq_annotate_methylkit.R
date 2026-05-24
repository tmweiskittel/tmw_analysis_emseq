#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(methylKit)
  library(optparse)
  library(data.table)
})

opts <- parse_args(OptionParser(option_list = list(
  make_option("--db", type = "character"),
  make_option("--out", type = "character")
)))

dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)

mdiff <- readMethylDiffDB(opts$db)
df <- as.data.table(getData(mdiff))

# Placeholder annotation fields; reference-based gene/TSS annotation can be added later.
df[, feature_type := NA_character_]
df[, gene_id := NA_character_]
df[, gene_name := NA_character_]
df[, nearest_tss := NA_character_]
df[, distance_to_tss := NA_integer_]

fwrite(df, opts$out, sep = "\t", quote = FALSE, na = "NA")
