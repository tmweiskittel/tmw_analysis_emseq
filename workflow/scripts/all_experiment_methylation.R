#!/usr/bin/env Rscript
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opts <- parse_args(
  OptionParser(
    option_list = list(
      make_option(
        "--db_file",
        type = "character",
        dest = "db_file",
        help = "Input methylBase BGZF database"
      ),
      make_option(
        "--out_file",
        type = "character",
        dest = "out_file",
        help = "Output position-by-sample methylation matrix"
      )
    )
  )
)

if (is.null(opts$db_file) || !nzchar(opts$db_file)) {
  stop("--db_file is required")
}

if (!file.exists(opts$db_file)) {
  stop("Input database does not exist: ", opts$db_file)
}

if (is.null(opts$out_file) || !nzchar(opts$out_file)) {
  stop("--out_file is required")
}

dir.create(
  dirname(opts$out_file),
  recursive = TRUE,
  showWarnings = FALSE
)

message("Reading methylKit metadata")

metadata <- system(
  sprintf(
    "bgzip -dc %s | grep '^#'",
    shQuote(opts$db_file)
  ),
  intern = TRUE
)

if (length(metadata) == 0L) {
  stop("No methylKit metadata lines were found")
}

get_metadata_value <- function(prefix) {
  hit <- metadata[startsWith(metadata, prefix)]

  if (length(hit) == 0L) {
    stop("Missing methylKit metadata field: ", prefix)
  }

  sub(
    paste0("^", prefix),
    "",
    hit[1L]
  )
}

sample_ids <- strsplit(
  get_metadata_value("#SI:"),
  ";",
  fixed = TRUE
)[[1L]]

coverage_indices <- as.integer(
  strsplit(
    get_metadata_value("#CI:"),
    ";",
    fixed = TRUE
  )[[1L]]
)

num_c_indices <- as.integer(
  strsplit(
    get_metadata_value("#NC:"),
    ";",
    fixed = TRUE
  )[[1L]]
)

if (
  length(sample_ids) != length(coverage_indices) ||
  length(sample_ids) != length(num_c_indices)
) {
  stop(
    "Metadata lengths do not match: samples=",
    length(sample_ids),
    ", coverage=",
    length(coverage_indices),
    ", numCs=",
    length(num_c_indices)
  )
}

sample_ids <- sub(
  "\\.hg38\\.bwameth\\.methyldackel$",
  "",
  sample_ids
)

sample_ids <- make.unique(sample_ids)

tmp_file <- tempfile(
  pattern = "methylBase_",
  fileext = ".txt"
)

on.exit(
  unlink(tmp_file, force = TRUE),
  add = TRUE
)

message("Decompressing methylBase database")

status <- system2(
  command = "bgzip",
  args = c("-dc", shQuote(opts$db_file)),
  stdout = tmp_file
)

if (!identical(status, 0L)) {
  stop(
    "Failed to decompress methylBase database; exit status: ",
    status
  )
}

message("Reading decompressed methylBase table")

df <- fread(
  tmp_file,
  sep = "\t",
  header = FALSE,
  skip = length(metadata),
  showProgress = TRUE
)

expected_columns <- max(
  c(
    4L,
    coverage_indices,
    num_c_indices
  )
)

if (ncol(df) < expected_columns) {
  stop(
    "Unexpected methylBase column count. Found ",
    ncol(df),
    " columns, but metadata references column ",
    expected_columns
  )
}

setnames(
  df,
  1:4,
  c("chr", "start", "end", "strand")
)

message(
  "Loaded ",
  format(nrow(df), big.mark = ","),
  " rows and ",
  ncol(df),
  " columns"
)

out <- df[, .(
  chr,
  start,
  end,
  strand
)]

for (i in seq_along(sample_ids)) {
  message(
    "Calculating methylation for ",
    sample_ids[i],
    " (",
    i,
    "/",
    length(sample_ids),
    ")"
  )

  coverage <- df[[coverage_indices[i]]]
  num_c <- df[[num_c_indices[i]]]

  out[, (sample_ids[i]) := fifelse(
    coverage > 0,
    100 * num_c / coverage,
    NA_real_
  )]

  rm(coverage, num_c)

  if (i %% 10L == 0L) {
    invisible(gc())
  }
}

message("Writing output: ", opts$out_file)

fwrite(
  out,
  file = opts$out_file,
  sep = "\t",
  quote = FALSE,
  na = "NA",
  showProgress = TRUE
)

message(
  "Finished: ",
  format(nrow(out), big.mark = ","),
  " rows x ",
  ncol(out),
  " columns"
)
