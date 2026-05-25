#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(methylKit)
  library(optparse)
})

opts <- parse_args(OptionParser(option_list = list(
  make_option("--lib_db_list", type = "character"),
  make_option("--lib_id_list", type = "character"),
  make_option("--treatment_list", type = "character"),
  make_option("--cores", type = "integer", default = 1),
  make_option("--out_dir", type = "character"),
  make_option("--suffix", type = "character"),
  make_option("--assembly", type = "character", default = "hg38"),
  make_option("--mincov", type = "integer", default = 5),
  make_option("--win_size", type = "integer", default = 1000),
  make_option("--min_per_group", type = "integer", default = 2),
  make_option("--chunk_size", type = "integer", default = 1000000)
)))

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

files <- strsplit(opts$lib_db_list, " ")[[1]]
ids <- strsplit(opts$lib_id_list, " ")[[1]]
tx <- as.integer(strsplit(opts$treatment_list, " ")[[1]])

meth <- methRead(
  location = as.list(files),
  sample.id = as.list(ids),
  assembly = opts$assembly,
  treatment = tx,
  context = "CpG",
  mincov = opts$mincov,
  dbtype = "tabix",
  dbdir = opts$out_dir
)

tiles <- tileMethylCounts(
  meth,
  win.size = opts$win_size,
  step.size = opts$win_size,
  cov.bases = 1,
  mc.cores = opts$cores
)

mbase <- unite(
  tiles,
  destrand = FALSE,
  min.per.group = opts$min_per_group,
  mc.cores = opts$cores,
  save.db = TRUE,
  dbdir = opts$out_dir,
  suffix = paste0(opts$suffix, ".tiled")
)

calculateDiffMeth(
  mbase,
  mc.cores = opts$cores,
  save.db = TRUE,
  dbdir = opts$out_dir,
  suffix = paste0(opts$suffix, ".tiled")
)
