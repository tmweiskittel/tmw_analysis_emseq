#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

option_list <- list(
  make_option(
    "--qc",
    type = "character",
    dest = "qc"
  ),
  make_option(
    "--outdir",
    type = "character",
    dest = "outdir"
  )
)

opts <- parse_args(
  OptionParser(option_list = option_list)
)

if (
  is.null(opts$qc) ||
  is.null(opts$outdir)
) {
  stop("--qc and --outdir are required")
}

dir.create(
  opts$outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

qc <- fread(
  opts$qc,
  na.strings = c("NA", "", "NaN")
)

required <- c(
  "sample",
  "group"
)

missing <- setdiff(
  required,
  names(qc)
)

if (length(missing) > 0) {
  stop(
    "Missing required column(s): ",
    paste(missing, collapse = ", ")
  )
}

# group is the cohort assignment supplied by the analysis pipeline.
qc[, cohort := factor(group)]


# ============================================================================
# Metrics to visualize
# ============================================================================

metric_labels <- c(
  raw_fastq_coverage =
    "Raw sequencing coverage (x)",

  mean_aligned_base_coverage =
    "Mean aligned genome coverage (x)",

  mean_coverage_called_cpgs =
    "Mean called-CpG coverage (x)",

  median_coverage_called_cpgs =
    "Median called-CpG coverage (x)",

  coverage_weighted_methylation_fraction =
    "Coverage-weighted methylation fraction",

  fastp_before_q30_rate =
    "Raw Q30 base fraction",

  fastp_after_q30_rate =
    "Post-filter Q30 base fraction",

  fastp_before_gc_content =
    "Raw GC fraction",

  fastp_after_gc_content =
    "Post-filter GC fraction",

  fastp_duplication_rate =
    "fastp duplication rate",

  fastp_peak_insert_size =
    "Peak insert size (bp)",

  percent_reads_retained =
    "Fraction of reads retained",

  percent_mapped_reads_retained =
    "Fraction of mapped reads retained",

  final_human_mapped_reads =
    "Final human mapped reads",

  cpg_sites_called =
    "Called CpG sites",

  lambda_mean_methylation_fraction =
    "Lambda methylation fraction"
)

metrics <- intersect(
  names(metric_labels),
  names(qc)
)

if (length(metrics) == 0) {
  stop("No recognized QC metrics were found.")
}


# ============================================================================
# Sample-level plots
# ============================================================================

for (metric in metrics) {

  dt <- qc[
    !is.na(get(metric)),
    .(
      sample,
      cohort,
      value = get(metric)
    )
  ]

  if (nrow(dt) == 0) {
    next
  }

  # Sort samples by value so low/high outliers are obvious.
  sample_order <- dt[
    order(value),
    sample
  ]

  dt[
    ,
    sample := factor(
      sample,
      levels = sample_order
    )
  ]

  p <- ggplot(
    dt,
    aes(
      x = sample,
      y = value,
      color = cohort
    )
  ) +
    geom_point(
      size = 3
    ) +
    labs(
      title = metric_labels[[metric]],
      x = "Sample",
      y = metric_labels[[metric]],
      color = "Cohort"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(
      opts$outdir,
      paste0(metric, ".samples.png")
    ),
    plot = p,
    width = 11,
    height = 5,
    dpi = 300
  )
}


# ============================================================================
# Cohort comparison plots
# ============================================================================

for (metric in metrics) {

  dt <- qc[
    !is.na(get(metric)),
    .(
      sample,
      cohort,
      value = get(metric)
    )
  ]

  if (
    nrow(dt) == 0 ||
    uniqueN(dt$cohort) < 2
  ) {
    next
  }

  p <- ggplot(
    dt,
    aes(
      x = cohort,
      y = value,
      color = cohort
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.5
    ) +
    geom_jitter(
      width = 0.12,
      height = 0,
      size = 2.5,
      alpha = 0.8
    ) +
    labs(
      title = metric_labels[[metric]],
      x = "Cohort",
      y = metric_labels[[metric]]
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(
      opts$outdir,
      paste0(metric, ".cohort.png")
    ),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}


# ============================================================================
# Raw sequencing coverage vs aligned coverage
# ============================================================================

if (
  all(
    c(
      "raw_fastq_coverage",
      "mean_aligned_base_coverage"
    ) %in% names(qc)
  )
) {

  dt <- qc[
    !is.na(raw_fastq_coverage) &
    !is.na(mean_aligned_base_coverage)
  ]

  p <- ggplot(
    dt,
    aes(
      x = raw_fastq_coverage,
      y = mean_aligned_base_coverage,
      color = cohort
    )
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    geom_point(
      size = 3
    ) +
    labs(
      title = "Sequencing yield vs aligned genome coverage",
      x = "Raw sequencing coverage (x)",
      y = "Mean aligned genome coverage (x)",
      color = "Cohort"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(
      opts$outdir,
      "raw_vs_aligned_coverage.png"
    ),
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
}


# ============================================================================
# Aligned coverage vs CpG coverage
# ============================================================================

if (
  all(
    c(
      "mean_aligned_base_coverage",
      "mean_coverage_called_cpgs"
    ) %in% names(qc)
  )
) {

  dt <- qc[
    !is.na(mean_aligned_base_coverage) &
    !is.na(mean_coverage_called_cpgs)
  ]

  p <- ggplot(
    dt,
    aes(
      x = mean_aligned_base_coverage,
      y = mean_coverage_called_cpgs,
      color = cohort
    )
  ) +
    geom_point(
      size = 3
    ) +
    labs(
      title = "Aligned genome coverage vs CpG coverage",
      x = "Mean aligned genome coverage (x)",
      y = "Mean called-CpG coverage (x)",
      color = "Cohort"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(
      opts$outdir,
      "aligned_vs_cpg_coverage.png"
    ),
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
}


# ============================================================================
# Coverage metrics together
# ============================================================================

coverage_metrics <- intersect(
  c(
    "raw_fastq_coverage",
    "mean_aligned_base_coverage",
    "mean_coverage_called_cpgs",
    "median_coverage_called_cpgs"
  ),
  names(qc)
)

if (length(coverage_metrics) >= 2) {

  coverage_long <- melt(
    qc,
    id.vars = c(
      "sample",
      "cohort"
    ),
    measure.vars = coverage_metrics,
    variable.name = "coverage_type",
    value.name = "coverage"
  )

  coverage_long <- coverage_long[
    !is.na(coverage)
  ]

  coverage_long[
    ,
    coverage_type := factor(
      coverage_type,
      levels = coverage_metrics,
      labels = unname(
        metric_labels[coverage_metrics]
      )
    )
  ]

  if ("raw_fastq_coverage" %in% names(qc)) {
    sample_order <- qc[
      order(raw_fastq_coverage),
      sample
    ]
  } else {
    sample_order <- unique(qc$sample)
  }

  coverage_long[
    ,
    sample := factor(
      sample,
      levels = sample_order
    )
  ]

  p <- ggplot(
    coverage_long,
    aes(
      x = sample,
      y = coverage,
      color = cohort,
      shape = coverage_type
    )
  ) +
    geom_point(
      size = 2.7,
      position = position_dodge(
        width = 0.5
      )
    ) +
    labs(
      title = "Coverage metrics by sample",
      x = "Sample",
      y = "Coverage (x)",
      color = "Cohort",
      shape = "Metric"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(
      opts$outdir,
      "coverage_metrics_by_sample.png"
    ),
    plot = p,
    width = 12,
    height = 6,
    dpi = 300
  )
}
