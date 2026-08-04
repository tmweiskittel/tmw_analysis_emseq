#!/usr/bin/env Rscript

mbase <- readMethylDB(opts$db_file)

tmp_file <- tempfile(fileext = ".txt")
on.exit(unlink(tmp_file), add = TRUE)

system2(
    "bgzip",
    c("-dc", shQuote(opts$db_file)),
    stdout = tmp_file
)

df <- fread(
    tmp_file,
    sep = "\t",
    header = FALSE,
    skip = length(Rsamtools::headerTabix(opts$db_file)$header)
)

coverage_cols <- grep("\\.coverage$", names(df), value = TRUE)
numCs_cols <- sub("\\.coverage$", ".numCs", coverage_cols)

out <- df[, .(chr, start, end, strand)]

for (i in seq_along(coverage_cols)) {
  sample <- sub("\\.coverage$", "", coverage_cols[i])
  out[[sample]] <- 100 * df[[numCs_cols[i]]] / df[[coverage_cols[i]]]
}

fwrite(out, opts$out_file, sep = "\t", quote = FALSE, na = "NA")
