#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opts <- parse_args(OptionParser(option_list = list(
  make_option("--qc-files", type = "character"),
  make_option("--samples", type = "character"),
  make_option("--groups", type = "character"),
  make_option("--out-aggregate", type = "character"),
  make_option("--out-stats", type = "character")
)))

files <- strsplit(opts$qc_files, " ")[[1]]
samples <- strsplit(opts$samples, " ")[[1]]
groups <- strsplit(opts$groups, " ")[[1]]

if (length(files) != length(samples) || length(samples) != length(groups)) {
  stop("files, samples, and groups must have same length")
}

qc_list <- lapply(seq_along(files), function(i) {
  x <- fread(files[i])
  x[, sample := samples[i]]
  x[, group := groups[i]]
  x
})

qc <- rbindlist(qc_list, fill = TRUE)

setcolorder(qc, c("sample", "group", setdiff(names(qc), c("sample", "group"))))

dir.create(dirname(opts$out_aggregate), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opts$out_stats), recursive = TRUE, showWarnings = FALSE)

fwrite(qc, opts$out_aggregate, sep = "\t", quote = FALSE, na = "NA")

metric_cols <- setdiff(names(qc), c("sample", "group"))
numeric_cols <- metric_cols[sapply(qc[, ..metric_cols], is.numeric)]

results <- rbindlist(lapply(numeric_cols, function(metric) {
  dt <- qc[!is.na(get(metric)), .(value = get(metric), group)]

  if (length(unique(dt$group)) != 2) {
    return(NULL)
  }

  g <- unique(dt$group)
  x <- dt[group == g[1], value]
  y <- dt[group == g[2], value]

  test <- wilcox.test(x, y, exact = FALSE)

  data.table(
    metric = metric,
    group_1 = g[1],
    group_2 = g[2],
    n_group_1 = length(x),
    n_group_2 = length(y),
    median_group_1 = median(x, na.rm = TRUE),
    median_group_2 = median(y, na.rm = TRUE),
    mean_group_1 = mean(x, na.rm = TRUE),
    mean_group_2 = mean(y, na.rm = TRUE),
    difference_mean_group2_minus_group1 = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),
    test = "wilcoxon_rank_sum",
    p_value = test$p.value
  )
}), fill = TRUE)

if (nrow(results) > 0) {
  results[, q_value := p.adjust(p_value, method = "BH")]
}

fwrite(results, opts$out_stats, sep = "\t", quote = FALSE, na = "NA")
