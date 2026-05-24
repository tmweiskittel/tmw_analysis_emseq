#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(methylKit)
  library(optparse)
})

opts <- parse_args(OptionParser(option_list = list(
  make_option("--mbase", type = "character"),
  make_option("--cores", type = "integer", default = 1),
  make_option("--out_dir", type = "character"),
  make_option("--suffix", type = "character"),
  make_option("--chunk_size", type = "integer", default = 1000000)
)))

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

mbase <- readMethylBaseDB(opts$mbase)

calculateDiffMeth(
  mbase,
  mc.cores = opts$cores,
  save.db = TRUE,
  dbdir = opts$out_dir,
  suffix = opts$suffix
)
