#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})


# ============================================================================
# Command-line arguments
# ============================================================================

option_list <- list(
  make_option(
    "--qc-files",
    type = "character",
    dest = "qc_files",
    help = "Whitespace-separated QC summary TSV files."
  ),
  make_option(
    "--samples",
    type = "character",
    dest = "samples",
    help = "Whitespace-separated sample identifiers."
  ),
  make_option(
    "--groups",
    type = "character",
    dest = "groups",
    help = "Whitespace-separated group labels."
  ),
  make_option(
    "--out-aggregate",
    type = "character",
    dest = "out_aggregate",
    help = "Output aggregated QC TSV."
  ),
  make_option(
    "--out-stats",
    type = "character",
    dest = "out_stats",
    help = "Output group-comparison statistics TSV."
  )
)

opts <- parse_args(
  OptionParser(option_list = option_list)
)

required_arguments <- c(
  "qc_files",
  "samples",
  "groups",
  "out_aggregate",
  "out_stats"
)

missing_arguments <- required_arguments[
  vapply(
    required_arguments,
    function(argument_name) {
      value <- opts[[argument_name]]

      is.null(value) ||
        is.na(value) ||
        trimws(value) == ""
    },
    logical(1)
  )
]

if (length(missing_arguments) > 0) {
  stop(
    "Missing required argument(s): ",
    paste(missing_arguments, collapse = ", ")
  )
}


# ============================================================================
# Parse vector arguments
# ============================================================================

split_argument <- function(value) {
  result <- strsplit(
    trimws(value),
    "\\s+"
  )[[1]]

  result[nzchar(result)]
}

files <- split_argument(opts$qc_files)
samples <- split_argument(opts$samples)
groups <- split_argument(opts$groups)

if (
  length(files) != length(samples) ||
  length(samples) != length(groups)
) {
  stop(
    "files, samples, and groups must have the same length. ",
    "files=", length(files),
    ", samples=", length(samples),
    ", groups=", length(groups)
  )
}

if (length(files) == 0) {
  stop("No QC files were provided.")
}

if (anyDuplicated(samples)) {
  duplicate_samples <- unique(
    samples[duplicated(samples)]
  )

  stop(
    "Duplicate sample identifiers were supplied: ",
    paste(duplicate_samples, collapse = ", ")
  )
}

missing_files <- files[!file.exists(files)]

if (length(missing_files) > 0) {
  stop(
    "Missing QC file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}


# ============================================================================
# QC schema definitions
# ============================================================================

# These fields came from the old coverage implementation and should not appear
# in newly reconstructed QC summaries.
obsolete_coverage_fields <- c(
  "total_methylated_counts",
  "total_unmethylated_counts",
  "mean_coverage",
  "median_coverage",
  "mean_methylation_fraction",
  "cpg_sites_ge_1x",
  "cpg_sites_ge_5x",
  "cpg_sites_ge_10x",
  "cpg_sites_ge_20x",
  "cpg_sites_ge_30x",
  "cpg_sites_ge_1x_pct",
  "cpg_sites_ge_5x_pct",
  "cpg_sites_ge_10x_pct",
  "cpg_sites_ge_20x_pct",
  "cpg_sites_ge_30x_pct",
  "reference_bases_ge_1x",
  "reference_bases_ge_5x",
  "reference_bases_ge_10x",
  "reference_bases_ge_20x",
  "reference_bases_ge_30x",
  "reference_bases_ge_1x_pct",
  "reference_bases_ge_5x_pct",
  "reference_bases_ge_10x_pct",
  "reference_bases_ge_20x_pct",
  "reference_bases_ge_30x_pct",
  "estimated_mapped_bases",
  "estimated_genome_depth"
)

# Core fields expected from the corrected coverage implementation.
required_new_coverage_fields <- c(
  "vendor_genome_size_denominator",
  "raw_fastq_coverage",
  "mean_aligned_base_coverage",
  "aligned_reference_size_denominator",
  "minimum_mapping_quality",
  "minimum_base_quality",
  "overlapping_mates_counted_once",
  "excluded_depth_contigs",
  "cpg_sites_called",
  "approx_total_methylated_counts",
  "approx_total_unmethylated_counts",
  "mean_coverage_called_cpgs",
  "median_coverage_called_cpgs",
  "coverage_weighted_methylation_fraction"
)

# Numeric columns that are configuration values, fixed denominators, or
# technical provenance. They belong in the aggregate table but should not be
# treated as sample-level outcomes in group comparisons.
non_test_numeric_fields <- c(
  "vendor_genome_size_denominator",
  "aligned_reference_size_denominator",
  "minimum_mapping_quality",
  "minimum_base_quality"
)

# Character or logical metadata that should never enter numeric tests.
non_metric_fields <- c(
  "sample",
  "group",
  "overlapping_mates_counted_once",
  "excluded_depth_contigs"
)


# ============================================================================
# Read and validate QC files
# ============================================================================

read_qc_file <- function(
  filename,
  expected_sample,
  expected_group
) {
  qc_row <- fread(
    filename,
    na.strings = c("NA", "", "NaN", "nan"),
    check.names = FALSE,
    showProgress = FALSE
  )

  if (nrow(qc_row) != 1) {
    stop(
      "Expected exactly one data row in QC file ",
      filename,
      "; found ",
      nrow(qc_row),
      "."
    )
  }

  # The command-line sample assignment is authoritative for aggregation.
  # Remove a file-provided sample column before replacing it so duplicate
  # columns cannot be created.
  if ("sample" %in% names(qc_row)) {
    file_sample <- as.character(qc_row[["sample"]][1])

    if (
      !is.na(file_sample) &&
      nzchar(file_sample) &&
      file_sample != expected_sample
    ) {
      stop(
        "Sample mismatch for ",
        filename,
        ": command line specifies ",
        expected_sample,
        " but the file contains ",
        file_sample,
        "."
      )
    }

    qc_row[, sample := NULL]
  }

  qc_row[, sample := expected_sample]
  qc_row[, group := expected_group]

  qc_row
}


qc_list <- lapply(
  seq_along(files),
  function(index) {
    read_qc_file(
      filename = files[index],
      expected_sample = samples[index],
      expected_group = groups[index]
    )
  }
)


# ============================================================================
# Prevent mixed old/new schemas
# ============================================================================

file_schemas <- lapply(
  qc_list,
  function(qc_row) {
    setdiff(
      names(qc_row),
      c("sample", "group")
    )
  }
)

reference_schema <- file_schemas[[1]]

schema_mismatches <- which(
  !vapply(
    file_schemas,
    function(current_schema) {
      identical(
        sort(current_schema),
        sort(reference_schema)
      )
    },
    logical(1)
  )
)

if (length(schema_mismatches) > 0) {
  mismatch_descriptions <- vapply(
    schema_mismatches,
    function(index) {
      missing_from_file <- setdiff(
        reference_schema,
        file_schemas[[index]]
      )

      extra_in_file <- setdiff(
        file_schemas[[index]],
        reference_schema
      )

      paste0(
        files[index],
        "\n  missing: ",
        if (length(missing_from_file) > 0) {
          paste(missing_from_file, collapse = ", ")
        } else {
          "none"
        },
        "\n  extra: ",
        if (length(extra_in_file) > 0) {
          paste(extra_in_file, collapse = ", ")
        } else {
          "none"
        }
      )
    },
    character(1)
  )

  stop(
    "QC files do not have a consistent schema. ",
    "Old and new QC summaries should not be aggregated together.\n",
    paste(mismatch_descriptions, collapse = "\n")
  )
}

present_obsolete_fields <- intersect(
  obsolete_coverage_fields,
  reference_schema
)

if (length(present_obsolete_fields) > 0) {
  stop(
    "The QC files still contain obsolete coverage fields:\n",
    paste(present_obsolete_fields, collapse = "\n"),
    "\nRegenerate the QC summaries before running aggregation."
  )
}

missing_new_coverage_fields <- setdiff(
  required_new_coverage_fields,
  reference_schema
)

if (length(missing_new_coverage_fields) > 0) {
  stop(
    "The QC summaries are missing corrected coverage fields:\n",
    paste(missing_new_coverage_fields, collapse = "\n")
  )
}


# ============================================================================
# Aggregate all samples
# ============================================================================

qc <- rbindlist(
  qc_list,
  use.names = TRUE,
  fill = FALSE
)

setcolorder(
  qc,
  c(
    "sample",
    "group",
    setdiff(
      names(qc),
      c("sample", "group")
    )
  )
)

# Preserve the supplied group order rather than relying on whichever group
# happens to appear first after data manipulation.
group_order <- unique(groups)

qc[, group := factor(
  group,
  levels = group_order
)]

dir.create(
  dirname(opts$out_aggregate),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(opts$out_stats),
  recursive = TRUE,
  showWarnings = FALSE
)

fwrite(
  qc,
  opts$out_aggregate,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================================
# Select testable numeric metrics
# ============================================================================

metric_columns <- setdiff(
  names(qc),
  non_metric_fields
)

numeric_columns <- metric_columns[
  vapply(
    qc[, ..metric_columns],
    is.numeric,
    logical(1)
  )
]

testable_numeric_columns <- setdiff(
  numeric_columns,
  non_test_numeric_fields
)

if (length(group_order) != 2) {
  stop(
    "Wilcoxon comparisons require exactly two groups. Found: ",
    paste(group_order, collapse = ", ")
  )
}


# ============================================================================
# Compare groups
# ============================================================================

compare_metric <- function(metric) {
  group_1 <- group_order[1]
  group_2 <- group_order[2]

  values_group_1 <- qc[
    group == group_1 & !is.na(get(metric)),
    get(metric)
  ]

  values_group_2 <- qc[
    group == group_2 & !is.na(get(metric)),
    get(metric)
  ]

  total_group_1 <- qc[group == group_1, .N]
  total_group_2 <- qc[group == group_2, .N]

  missing_group_1 <- total_group_1 - length(values_group_1)
  missing_group_2 <- total_group_2 - length(values_group_2)

  if (
    length(values_group_1) == 0 ||
    length(values_group_2) == 0
  ) {
    return(NULL)
  }

  all_values <- c(
    values_group_1,
    values_group_2
  )

  # A completely constant metric has no meaningful rank separation.
  if (length(unique(all_values)) < 2) {
    p_value <- NA_real_
    test_status <- "constant_metric"
  } else {
    test_result <- tryCatch(
      wilcox.test(
        values_group_1,
        values_group_2,
        exact = FALSE
      ),
      error = function(error_condition) {
        error_condition
      }
    )

    if (inherits(test_result, "error")) {
      p_value <- NA_real_
      test_status <- paste0(
        "test_error: ",
        conditionMessage(test_result)
      )
    } else {
      p_value <- test_result$p.value
      test_status <- "ok"
    }
  }

  data.table(
    metric = metric,
    group_1 = group_1,
    group_2 = group_2,

    n_group_1 = length(values_group_1),
    n_group_2 = length(values_group_2),

    missing_group_1 = missing_group_1,
    missing_group_2 = missing_group_2,

    median_group_1 = median(
      values_group_1,
      na.rm = TRUE
    ),
    median_group_2 = median(
      values_group_2,
      na.rm = TRUE
    ),

    mean_group_1 = mean(
      values_group_1,
      na.rm = TRUE
    ),
    mean_group_2 = mean(
      values_group_2,
      na.rm = TRUE
    ),

    difference_mean_group2_minus_group1 = (
      mean(values_group_2, na.rm = TRUE) -
        mean(values_group_1, na.rm = TRUE)
    ),

    test = "wilcoxon_rank_sum",
    test_status = test_status,
    p_value = p_value
  )
}


result_list <- lapply(
  testable_numeric_columns,
  compare_metric
)

result_list <- Filter(
  Negate(is.null),
  result_list
)

if (length(result_list) > 0) {
  results <- rbindlist(
    result_list,
    use.names = TRUE,
    fill = TRUE
  )

  results[, q_value := p.adjust(
    p_value,
    method = "BH"
  )]

  setorder(
    results,
    q_value,
    p_value,
    metric,
    na.last = TRUE
  )
} else {
  results <- data.table(
    metric = character(),
    group_1 = character(),
    group_2 = character(),
    n_group_1 = integer(),
    n_group_2 = integer(),
    missing_group_1 = integer(),
    missing_group_2 = integer(),
    median_group_1 = numeric(),
    median_group_2 = numeric(),
    mean_group_1 = numeric(),
    mean_group_2 = numeric(),
    difference_mean_group2_minus_group1 = numeric(),
    test = character(),
    test_status = character(),
    p_value = numeric(),
    q_value = numeric()
  )
}

fwrite(
  results,
  opts$out_stats,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)
