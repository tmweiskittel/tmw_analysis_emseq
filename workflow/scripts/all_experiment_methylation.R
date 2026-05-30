#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(methylKit)
  library(optparse)
  library(data.table)
})

opts <- parse_args(OptionParser(option_list = list(
  make_option("--db_file", type = "character"),
  make_option("--out_file", type = "character")
)))

dir.create(dirname(opts$out_file), recursive = TRUE, showWarnings = FALSE)

mbase <- readMethylDB(opts$db_file)
df <- as.data.table(getData(mbase))

coverage_cols <- grep("\\.coverage$", names(df), value = TRUE)
numCs_cols <- sub("\\.coverage$", ".numCs", coverage_cols)

out <- df[, .(chr, start, end, strand)]

for (i in seq_along(coverage_cols)) {
  sample <- sub("\\.coverage$", "", coverage_cols[i])
  out[[sample]] <- 100 * df[[numCs_cols[i]]] / df[[coverage_cols[i]]]
}

fwrite(out, opts$out_file, sep = "\t", quote = FALSE, na = "NA")
