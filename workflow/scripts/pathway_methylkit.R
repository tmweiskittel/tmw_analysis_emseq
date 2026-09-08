#!/usr/bin/env Rscript

# ==============================================================================
# Generic methylKit tiled methylation pathway analysis
#
# Usage:
#
#   Rscript pathway_methylkit.R \
#       --contrast recurrcsf_vs_trtcsf
#
# Optional:
#
#   --root /home/jupyter/data/analysis/dmr
#   --primary-delta 20
#   --sensitivity-delta 15
#   --q 0.05
#   --group0 recurrcsf
#   --group1 trtcsf
#
# Expected methylKit input:
#
#   <root>/diff/methylDiff_<contrast>.tiled.txt.bgz
#
# Expected methylKit columns:
#
#   chr
#   start
#   end
#   strand
#   pvalue
#   qvalue
#   meth.diff
#
# methylKit interpretation:
#
#   meth.diff = treatment1 - treatment0
#
# Thus:
#
#   meth.diff > 0 -> group1 hypermethylated
#   meth.diff < 0 -> group0 hypermethylated
#
# The script assumes that a contrast named:
#
#   A_vs_B
#
# was generated with:
#
#   treatment 0 = A
#   treatment 1 = B
#
# This can be overridden explicitly with:
#
#   --group0 A
#   --group1 B
#
# Analyses:
#
#   1. Primary promoter ORA
#      - GO Biological Process
#      - KEGG
#      - Reactome
#
#   2. Sensitivity promoter ORA
#
#   3. Secondary nearest-TSS <= 10 kb ORA
#
#   4. Ranked promoter GSEA
#      - GO BP
#      - KEGG
#      - Reactome
#
# Background:
#   Uses genes represented by tested EM-seq tiles, not all human genes.
#
# Genome:
#   hg38
#
# ==============================================================================


# ==============================================================================
# 1. Argument parser
# ==============================================================================

args <- commandArgs(
    trailingOnly = TRUE
)


get_arg <- function(
    flag,
    default = NULL
) {

    idx <- match(
        flag,
        args
    )

    if (
        is.na(idx) ||
        idx == length(args)
    ) {
        return(default)
    }

    args[idx + 1]
}


has_flag <- function(flag) {

    flag %in% args
}


CONTRAST <- get_arg(
    "--contrast"
)


if (is.null(CONTRAST)) {

    stop(
        paste0(
            "\nRequired argument missing.\n\n",
            "Usage:\n\n",
            "Rscript pathway_methylkit.R ",
            "--contrast recurrcsf_vs_trtcsf\n"
        )
    )
}


ROOT <- get_arg(
    "--root",
    "/home/jupyter/data/analysis/dmr"
)


Q_THRESHOLD <- as.numeric(
    get_arg(
        "--q",
        "0.05"
    )
)


PRIMARY_DELTA <- as.numeric(
    get_arg(
        "--primary-delta",
        "20"
    )
)


SENSITIVITY_DELTA <- as.numeric(
    get_arg(
        "--sensitivity-delta",
        "15"
    )
)


PROMOTER_UPSTREAM <- as.numeric(
    get_arg(
        "--promoter-upstream",
        "2000"
    )
)


PROMOTER_DOWNSTREAM <- as.numeric(
    get_arg(
        "--promoter-downstream",
        "500"
    )
)


NEAREST_TSS_MAX_DISTANCE <- as.numeric(
    get_arg(
        "--nearest-tss-distance",
        "10000"
    )
)


MIN_GS_SIZE <- as.numeric(
    get_arg(
        "--min-gs-size",
        "10"
    )
)


MAX_GS_SIZE <- as.numeric(
    get_arg(
        "--max-gs-size",
        "500"
    )
)


N_SHOW <- as.numeric(
    get_arg(
        "--show-category",
        "20"
    )
)


# ==============================================================================
# 2. Determine group names
# ==============================================================================

contrast_split <- strsplit(
    CONTRAST,
    "_vs_",
    fixed = TRUE
)[[1]]


if (length(contrast_split) != 2) {

    stop(
        paste0(
            "\nContrast must normally have the form:\n\n",
            "  group0_vs_group1\n\n",
            "Example:\n",
            "  recurrcsf_vs_trtcsf\n"
        )
    )
}


GROUP0 <- get_arg(
    "--group0",
    contrast_split[1]
)


GROUP1 <- get_arg(
    "--group1",
    contrast_split[2]
)


# ==============================================================================
# 3. Input/output paths
# ==============================================================================

DIFF_FILE <- file.path(
    ROOT,
    "diff",
    paste0(
        "methylDiff_",
        CONTRAST,
        ".tiled.txt.bgz"
    )
)


ANNOTATION_DIR <- file.path(
    ROOT,
    "annotation"
)


OUTPUT_DIR <- file.path(
    ROOT,
    "pathway",
    CONTRAST
)


dir.create(
    ANNOTATION_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ==============================================================================
# 4. Package checks
# ==============================================================================

required_packages <- c(
    "data.table",
    "GenomicRanges",
    "GenomicFeatures",
    "IRanges",
    "GenomeInfoDb",
    "S4Vectors",
    "AnnotationDbi",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "org.Hs.eg.db",
    "clusterProfiler",
    "ReactomePA",
    "enrichplot",
    "ggplot2"
)


missing_packages <- required_packages[
    !vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]


if (length(missing_packages) > 0) {

    stop(
        paste0(
            "\nMissing packages:\n",
            paste(
                paste0(
                    "  ",
                    missing_packages
                ),
                collapse = "\n"
            ),
            "\n\nInstall with BiocManager before running.\n"
        )
    )
}


suppressPackageStartupMessages({

    library(data.table)

    library(GenomicRanges)

    library(GenomicFeatures)

    library(IRanges)

    library(GenomeInfoDb)

    library(S4Vectors)

    library(AnnotationDbi)

    library(TxDb.Hsapiens.UCSC.hg38.knownGene)

    library(org.Hs.eg.db)

    library(clusterProfiler)

    library(ReactomePA)

    library(enrichplot)

    library(ggplot2)

})


data.table::setDTthreads(0)


set.seed(20260907)


# ==============================================================================
# 5. Utility functions
# ==============================================================================

section <- function(x) {

    message(
        "\n",
        paste(
            rep(
                "=",
                78
            ),
            collapse = ""
        ),
        "\n",
        x,
        "\n",
        paste(
            rep(
                "=",
                78
            ),
            collapse = ""
        )
    )
}


safe_name <- function(x) {

    gsub(
        "[^A-Za-z0-9_.-]",
        "_",
        x
    )
}


add_symbol_column <- function(
    dt,
    id_column,
    symbol_column
) {

    ids <- unique(
        as.character(
            dt[[id_column]]
        )
    )


    ids <- ids[
        !is.na(ids) &
        ids != ""
    ]


    if (length(ids) == 0) {

        dt[
            ,
            (symbol_column) := NA_character_
        ]

        return(dt)
    }


    symbol_map <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = ids,
        keytype = "ENTREZID",
        column = "SYMBOL",
        multiVals = "first"
    )


    dt[
        ,
        (symbol_column) :=
            as.character(
                symbol_map[
                    as.character(
                        get(id_column)
                    )
                ]
            )
    ]


    dt
}


save_enrichment <- function(
    obj,
    filename
) {

    path <- file.path(
        OUTPUT_DIR,
        filename
    )


    if (is.null(obj)) {

        data.table::fwrite(
            data.table(
                message = "Analysis unavailable or failed"
            ),
            path,
            sep = "\t"
        )

        return(
            invisible(NULL)
        )
    }


    result <- as.data.frame(
        obj
    )


    if (nrow(result) == 0) {

        data.table::fwrite(
            data.table(
                message = "No pathways passed enrichment thresholds"
            ),
            path,
            sep = "\t"
        )

        return(
            invisible(NULL)
        )
    }


    data.table::fwrite(
        result,
        path,
        sep = "\t",
        na = "NA"
    )


    invisible(NULL)
}


save_dotplot <- function(
    obj,
    filename,
    title
) {

    if (is.null(obj)) {
        return(
            invisible(NULL)
        )
    }


    result <- tryCatch(
        as.data.frame(obj),
        error = function(e) NULL
    )


    if (
        is.null(result) ||
        nrow(result) == 0
    ) {

        return(
            invisible(NULL)
        )
    }


    p <- enrichplot::dotplot(
        obj,
        showCategory = min(
            N_SHOW,
            nrow(result)
        )
    ) +
        ggplot2::ggtitle(
            title
        ) +
        ggplot2::theme(
            plot.title =
                ggplot2::element_text(
                    face = "bold"
                )
        )


    ggplot2::ggsave(
        filename = file.path(
            OUTPUT_DIR,
            filename
        ),
        plot = p,
        width = 10,
        height = 7,
        units = "in"
    )


    invisible(NULL)
}


# ==============================================================================
# 6. Validate input file
# ==============================================================================

section(
    paste0(
        "Contrast: ",
        CONTRAST
    )
)


message(
    "Group 0: ",
    GROUP0
)


message(
    "Group 1: ",
    GROUP1
)


message(
    "Input: ",
    DIFF_FILE
)


if (!file.exists(DIFF_FILE)) {

    stop(
        "\nInput file not found:\n",
        DIFF_FILE
    )
}


# ==============================================================================
# 7. Read methylKit metadata
# ==============================================================================

section(
    "Reading methylKit metadata"
)


metadata_lines <- system(
    sprintf(
        "zcat %s | grep '^#'",
        shQuote(
            DIFF_FILE
        )
    ),
    intern = TRUE
)


get_metadata <- function(prefix) {

    x <- metadata_lines[
        startsWith(
            metadata_lines,
            prefix
        )
    ]


    if (length(x) == 0) {
        return(NA_character_)
    }


    sub(
        paste0(
            "^",
            prefix
        ),
        "",
        x[1]
    )
}


GENOME <- get_metadata(
    "#AS:"
)


SAMPLE_STRING <- get_metadata(
    "#SI:"
)


TREATMENT_STRING <- get_metadata(
    "#TM:"
)


REGION_TYPE <- get_metadata(
    "#RS:"
)


sample_ids <- strsplit(
    SAMPLE_STRING,
    ";",
    fixed = TRUE
)[[1]]


treatments <- as.integer(
    strsplit(
        TREATMENT_STRING,
        ";",
        fixed = TRUE
    )[[1]]
)


sample_table <- data.table(

    sample = sample_ids,

    treatment = treatments
)


sample_table[
    ,
    group := ifelse(
        treatment == 0,
        GROUP0,
        GROUP1
    )
]


print(
    sample_table
)


data.table::fwrite(
    sample_table,
    file.path(
        OUTPUT_DIR,
        "00_samples.tsv"
    ),
    sep = "\t"
)


if (!identical(
    GENOME,
    "hg38"
)) {

    stop(
        paste0(
            "\nThis pipeline currently uses the hg38 TxDb.\n",
            "Input metadata reports genome: ",
            GENOME,
            "\n"
        )
    )
}


if (
    any(
        !treatments %in% c(
            0,
            1
        )
    )
) {

    stop(
        "Unexpected treatment coding in methylKit metadata."
    )
}


if (
    !any(
        treatments == 0
    ) ||
    !any(
        treatments == 1
    )
) {

    stop(
        "Both treatment 0 and treatment 1 must be represented."
    )
}


# ==============================================================================
# 8. Read methylKit differential tiles
# ==============================================================================

section(
    "Reading tiled differential methylation results"
)


read_command <- sprintf(
    "zcat %s | grep -v '^#'",
    shQuote(
        DIFF_FILE
    )
)


dm <- data.table::fread(
    cmd = read_command,
    header = FALSE,
    sep = "\t",
    col.names = c(
        "chr",
        "start",
        "end",
        "strand",
        "pvalue",
        "qvalue",
        "meth.diff"
    ),
    showProgress = TRUE
)


dm[
    ,
    tile_id := .I
]


data.table::setcolorder(
    dm,
    c(
        "tile_id",
        "chr",
        "start",
        "end",
        "strand",
        "pvalue",
        "qvalue",
        "meth.diff"
    )
)


message(
    "Loaded ",
    format(
        nrow(dm),
        big.mark = ","
    ),
    " tested tiles."
)


# ==============================================================================
# 9. Validate coordinates and tile structure
# ==============================================================================

section(
    "Input QC"
)


dm[
    ,
    width := end - start + 1
]


tile_width_summary <- dm[
    ,
    .N,
    by = width
][
    order(
        -N
    )
]


data.table::fwrite(
    tile_width_summary,
    file.path(
        OUTPUT_DIR,
        "01_tile_widths.tsv"
    ),
    sep = "\t"
)


message(
    "Most common tile width: ",
    tile_width_summary$width[1],
    " bp"
)


message(
    "Genome: ",
    GENOME
)


message(
    "Region type: ",
    REGION_TYPE
)


# ==============================================================================
# 10. Differential methylation QC
# ==============================================================================

GROUP0_HYPER_LABEL <- paste0(
    GROUP0,
    "_hypermethylated"
)


GROUP1_HYPER_LABEL <- paste0(
    GROUP1,
    "_hypermethylated"
)


qc_summary <- data.table(

    metric = c(
        "tested_tiles",
        "mean_meth_diff",
        "median_meth_diff",
        "positive_meth_diff",
        "negative_meth_diff",
        "q_le_threshold",
        "primary_group1_hyper",
        "primary_group0_hyper",
        "sensitivity_group1_hyper",
        "sensitivity_group0_hyper"
    ),

    value = c(

        nrow(dm),

        mean(
            dm$meth.diff,
            na.rm = TRUE
        ),

        median(
            dm$meth.diff,
            na.rm = TRUE
        ),

        sum(
            dm$meth.diff > 0,
            na.rm = TRUE
        ),

        sum(
            dm$meth.diff < 0,
            na.rm = TRUE
        ),

        sum(
            dm$qvalue <= Q_THRESHOLD,
            na.rm = TRUE
        ),

        sum(
            dm$qvalue <= Q_THRESHOLD &
            dm$meth.diff >= PRIMARY_DELTA,
            na.rm = TRUE
        ),

        sum(
            dm$qvalue <= Q_THRESHOLD &
            dm$meth.diff <= -PRIMARY_DELTA,
            na.rm = TRUE
        ),

        sum(
            dm$qvalue <= Q_THRESHOLD &
            dm$meth.diff >= SENSITIVITY_DELTA,
            na.rm = TRUE
        ),

        sum(
            dm$qvalue <= Q_THRESHOLD &
            dm$meth.diff <= -SENSITIVITY_DELTA,
            na.rm = TRUE
        )
    )
)


data.table::fwrite(
    qc_summary,
    file.path(
        OUTPUT_DIR,
        "02_differential_methylation_QC.tsv"
    ),
    sep = "\t"
)


print(
    qc_summary
)


# ==============================================================================
# 11. Significant tiles
# ==============================================================================

section(
    "Selecting significant regions"
)


primary <- dm[
    qvalue <= Q_THRESHOLD &
    abs(meth.diff) >= PRIMARY_DELTA
]


primary[
    ,
    direction :=
        ifelse(
            meth.diff > 0,
            GROUP1_HYPER_LABEL,
            GROUP0_HYPER_LABEL
        )
]


sensitivity <- dm[
    qvalue <= Q_THRESHOLD &
    abs(meth.diff) >= SENSITIVITY_DELTA
]


sensitivity[
    ,
    direction :=
        ifelse(
            meth.diff > 0,
            GROUP1_HYPER_LABEL,
            GROUP0_HYPER_LABEL
        )
]


data.table::fwrite(
    primary,
    file.path(
        OUTPUT_DIR,
        paste0(
            "03_primary_tiles_q",
            Q_THRESHOLD,
            "_abs",
            PRIMARY_DELTA,
            ".tsv.gz"
        )
    ),
    sep = "\t"
)


data.table::fwrite(
    sensitivity,
    file.path(
        OUTPUT_DIR,
        paste0(
            "04_sensitivity_tiles_q",
            Q_THRESHOLD,
            "_abs",
            SENSITIVITY_DELTA,
            ".tsv.gz"
        )
    ),
    sep = "\t"
)


# ==============================================================================
# 12. GRanges
# ==============================================================================

section(
    "Creating genomic intervals"
)


gr_all <- GenomicRanges::GRanges(

    seqnames = dm$chr,

    ranges = IRanges::IRanges(
        start = dm$start,
        end = dm$end
    ),

    strand = "*",

    tile_id = dm$tile_id,

    pvalue = dm$pvalue,

    qvalue = dm$qvalue,

    meth.diff = dm$meth.diff
)


GenomeInfoDb::genome(
    gr_all
) <- GENOME


# ==============================================================================
# 13. hg38 gene/TSS/promoter annotation
# ==============================================================================

section(
    "Building hg38 gene annotation"
)


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


tx_by_gene <- GenomicFeatures::transcriptsBy(
    txdb,
    by = "gene"
)


tx <- unlist(
    tx_by_gene,
    use.names = FALSE
)


tx_gene_id <- rep(
    names(
        tx_by_gene
    ),
    S4Vectors::elementNROWS(
        tx_by_gene
    )
)


# ------------------------------------------------------------------------------
# Promoters
# ------------------------------------------------------------------------------

promoters_tx <- GenomicFeatures::promoters(
    tx,
    upstream = PROMOTER_UPSTREAM,
    downstream = PROMOTER_DOWNSTREAM
)


S4Vectors::mcols(
    promoters_tx
)$gene_id <- tx_gene_id


promoters_by_gene <- split(
    promoters_tx,
    S4Vectors::mcols(
        promoters_tx
    )$gene_id
)


promoter_reduced_list <- GenomicRanges::reduce(
    promoters_by_gene
)


promoters_gene <- unlist(
    promoter_reduced_list,
    use.names = FALSE
)


S4Vectors::mcols(
    promoters_gene
)$gene_id <- rep(
    names(
        promoter_reduced_list
    ),
    S4Vectors::elementNROWS(
        promoter_reduced_list
    )
)


# ------------------------------------------------------------------------------
# TSS
# ------------------------------------------------------------------------------

tss_tx <- GenomicFeatures::promoters(
    tx,
    upstream = 0,
    downstream = 1
)


S4Vectors::mcols(
    tss_tx
)$gene_id <- tx_gene_id


tss_by_gene <- split(
    tss_tx,
    S4Vectors::mcols(
        tss_tx
    )$gene_id
)


tss_reduced_list <- GenomicRanges::reduce(
    tss_by_gene
)


tss_gene <- unlist(
    tss_reduced_list,
    use.names = FALSE
)


S4Vectors::mcols(
    tss_gene
)$gene_id <- rep(
    names(
        tss_reduced_list
    ),
    S4Vectors::elementNROWS(
        tss_reduced_list
    )
)


# Restrict annotations to chromosomes present in dataset.

dataset_chromosomes <- unique(
    dm$chr
)


promoters_gene <- promoters_gene[
    as.character(
        GenomicRanges::seqnames(
            promoters_gene
        )
    ) %in% dataset_chromosomes
]


tss_gene <- tss_gene[
    as.character(
        GenomicRanges::seqnames(
            tss_gene
        )
    ) %in% dataset_chromosomes
]


# ==============================================================================
# 14. Promoter mapping
# ==============================================================================

section(
    "Mapping tiles to promoters"
)


promoter_hits <- GenomicRanges::findOverlaps(
    gr_all,
    promoters_gene,
    ignore.strand = TRUE
)


promoter_map <- data.table(

    tile_id = S4Vectors::queryHits(
        promoter_hits
    ),

    gene_id = as.character(
        S4Vectors::mcols(
            promoters_gene
        )$gene_id[
            S4Vectors::subjectHits(
                promoter_hits
            )
        ]
    )
)


promoter_map <- unique(
    promoter_map,
    by = c(
        "tile_id",
        "gene_id"
    )
)


promoter_universe <- sort(
    unique(
        promoter_map$gene_id
    )
)


message(
    "Promoter background genes: ",
    length(
        promoter_universe
    )
)


promoter_stats <- merge(

    promoter_map,

    dm[
        ,
        .(
            tile_id,
            chr,
            start,
            end,
            pvalue,
            qvalue,
            meth.diff
        )
    ],

    by = "tile_id",

    all.x = TRUE,

    sort = FALSE
)


promoter_stats <- add_symbol_column(
    promoter_stats,
    "gene_id",
    "gene_symbol"
)


# ==============================================================================
# 15. Nearest TSS mapping
# ==============================================================================

section(
    "Mapping tiles to nearest TSS"
)


nearest_hits <- GenomicRanges::distanceToNearest(
    gr_all,
    tss_gene,
    ignore.strand = TRUE
)


nearest_map <- data.table(

    tile_id = S4Vectors::queryHits(
        nearest_hits
    ),

    gene_id = as.character(
        S4Vectors::mcols(
            tss_gene
        )$gene_id[
            S4Vectors::subjectHits(
                nearest_hits
            )
        ]
    ),

    distance_to_tss =
        S4Vectors::mcols(
            nearest_hits
        )$distance
)


nearest_map_10kb <- nearest_map[
    distance_to_tss <=
        NEAREST_TSS_MAX_DISTANCE
]


nearest_universe <- sort(
    unique(
        nearest_map_10kb$gene_id
    )
)


# ==============================================================================
# 16. Write significant tiled annotation
# ==============================================================================

section(
    "Writing tiled annotation"
)


primary_promoter <- promoter_stats[
    qvalue <= Q_THRESHOLD &
    abs(meth.diff) >= PRIMARY_DELTA
]


primary_promoter[
    ,
    direction :=
        ifelse(
            meth.diff > 0,
            GROUP1_HYPER_LABEL,
            GROUP0_HYPER_LABEL
        )
]


data.table::fwrite(
    primary_promoter,
    file.path(
        ANNOTATION_DIR,
        paste0(
            CONTRAST,
            "_tiled_promoter_annotated.tsv.gz"
        )
    ),
    sep = "\t"
)


nearest_primary <- merge(

    primary[
        ,
        .(
            tile_id,
            chr,
            start,
            end,
            strand,
            pvalue,
            qvalue,
            meth.diff,
            direction
        )
    ],

    nearest_map,

    by = "tile_id",

    all.x = TRUE,

    sort = FALSE
)


nearest_primary <- add_symbol_column(
    nearest_primary,
    "gene_id",
    "nearest_tss_gene_symbol"
)


data.table::setnames(
    nearest_primary,
    "gene_id",
    "nearest_tss_gene_id"
)


promoter_summary <- primary_promoter[
    ,
    .(
        promoter_gene_ids =
            paste(
                sort(
                    unique(
                        gene_id
                    )
                ),
                collapse = ";"
            ),

        promoter_gene_symbols =
            paste(
                sort(
                    unique(
                        gene_symbol[
                            !is.na(
                                gene_symbol
                            )
                        ]
                    )
                ),
                collapse = ";"
            )
    ),
    by = tile_id
]


primary_annotation <- merge(

    nearest_primary,

    promoter_summary,

    by = "tile_id",

    all.x = TRUE,

    sort = FALSE
)


primary_annotation[
    ,
    promoter_overlap :=
        !is.na(
            promoter_gene_ids
        ) &
        promoter_gene_ids != ""
]


data.table::setorder(
    primary_annotation,
    chr,
    start
)


data.table::fwrite(
    primary_annotation,
    file.path(
        ANNOTATION_DIR,
        paste0(
            CONTRAST,
            "_tiled_annotated.tsv.gz"
        )
    ),
    sep = "\t",
    na = "NA"
)


# ==============================================================================
# 17. Directional gene set helpers
# ==============================================================================

get_promoter_genes <- function(
    threshold,
    sign
) {

    if (sign == "positive") {

        return(
            unique(
                promoter_stats[
                    qvalue <= Q_THRESHOLD &
                    meth.diff >= threshold,
                    gene_id
                ]
            )
        )
    }


    unique(
        promoter_stats[
            qvalue <= Q_THRESHOLD &
            meth.diff <= -threshold,
            gene_id
        ]
    )
}


nearest_stats <- merge(

    nearest_map_10kb,

    dm[
        ,
        .(
            tile_id,
            qvalue,
            meth.diff
        )
    ],

    by = "tile_id",

    all.x = TRUE,

    sort = FALSE
)


get_nearest_genes <- function(
    threshold,
    sign
) {

    if (sign == "positive") {

        return(
            unique(
                nearest_stats[
                    qvalue <= Q_THRESHOLD &
                    meth.diff >= threshold,
                    gene_id
                ]
            )
        )
    }


    unique(
        nearest_stats[
            qvalue <= Q_THRESHOLD &
            meth.diff <= -threshold,
            gene_id
        ]
    )
}


prom_group1_20 <- get_promoter_genes(
    PRIMARY_DELTA,
    "positive"
)


prom_group0_20 <- get_promoter_genes(
    PRIMARY_DELTA,
    "negative"
)


prom_group1_15 <- get_promoter_genes(
    SENSITIVITY_DELTA,
    "positive"
)


prom_group0_15 <- get_promoter_genes(
    SENSITIVITY_DELTA,
    "negative"
)


near_group1_20 <- get_nearest_genes(
    PRIMARY_DELTA,
    "positive"
)


near_group0_20 <- get_nearest_genes(
    PRIMARY_DELTA,
    "negative"
)


gene_set_summary <- data.table(

    set = c(
        paste0(
            "promoter_",
            GROUP1_HYPER_LABEL,
            "_primary"
        ),
        paste0(
            "promoter_",
            GROUP0_HYPER_LABEL,
            "_primary"
        ),
        paste0(
            "promoter_",
            GROUP1_HYPER_LABEL,
            "_sensitivity"
        ),
        paste0(
            "promoter_",
            GROUP0_HYPER_LABEL,
            "_sensitivity"
        ),
        paste0(
            "nearest10kb_",
            GROUP1_HYPER_LABEL,
            "_primary"
        ),
        paste0(
            "nearest10kb_",
            GROUP0_HYPER_LABEL,
            "_primary"
        )
    ),

    n_genes = c(
        length(
            prom_group1_20
        ),
        length(
            prom_group0_20
        ),
        length(
            prom_group1_15
        ),
        length(
            prom_group0_15
        ),
        length(
            near_group1_20
        ),
        length(
            near_group0_20
        )
    )
)


data.table::fwrite(
    gene_set_summary,
    file.path(
        OUTPUT_DIR,
        "05_gene_set_counts.tsv"
    ),
    sep = "\t"
)


print(
    gene_set_summary
)


# ==============================================================================
# 18. ORA
# ==============================================================================

run_go <- function(
    genes,
    universe
) {

    if (
        length(
            genes
        ) < MIN_GS_SIZE
    ) {

        return(NULL)
    }


    tryCatch(

        clusterProfiler::enrichGO(

            gene = as.character(
                genes
            ),

            universe = as.character(
                universe
            ),

            OrgDb = org.Hs.eg.db,

            keyType = "ENTREZID",

            ont = "BP",

            pAdjustMethod = "BH",

            pvalueCutoff = 0.05,

            qvalueCutoff = 0.05,

            minGSSize = MIN_GS_SIZE,

            maxGSSize = MAX_GS_SIZE,

            readable = TRUE
        ),

        error = function(e) {

            warning(
                "GO failed: ",
                conditionMessage(e)
            )

            NULL
        }
    )
}


run_kegg <- function(
    genes,
    universe
) {

    if (
        length(
            genes
        ) < MIN_GS_SIZE
    ) {

        return(NULL)
    }


    tryCatch(

        clusterProfiler::enrichKEGG(

            gene = as.character(
                genes
            ),

            universe = as.character(
                universe
            ),

            organism = "hsa",

            keyType = "ncbi-geneid",

            pAdjustMethod = "BH",

            pvalueCutoff = 0.05,

            qvalueCutoff = 0.05,

            minGSSize = MIN_GS_SIZE,

            maxGSSize = MAX_GS_SIZE
        ),

        error = function(e) {

            warning(
                "KEGG failed: ",
                conditionMessage(e)
            )

            NULL
        }
    )
}


run_reactome <- function(
    genes,
    universe
) {

    if (
        length(
            genes
        ) < MIN_GS_SIZE
    ) {

        return(NULL)
    }


    tryCatch(

        ReactomePA::enrichPathway(

            gene = as.character(
                genes
            ),

            universe = as.character(
                universe
            ),

            organism = "human",

            pAdjustMethod = "BH",

            pvalueCutoff = 0.05,

            qvalueCutoff = 0.05,

            minGSSize = MIN_GS_SIZE,

            maxGSSize = MAX_GS_SIZE,

            readable = TRUE
        ),

        error = function(e) {

            warning(
                "Reactome failed: ",
                conditionMessage(e)
            )

            NULL
        }
    )
}


run_bundle <- function(
    genes,
    universe,
    label
) {

    label <- safe_name(
        label
    )


    message(
        "\n",
        label,
        ": ",
        length(
            genes
        ),
        " genes"
    )


    go <- run_go(
        genes,
        universe
    )


    kegg <- run_kegg(
        genes,
        universe
    )


    reactome <- run_reactome(
        genes,
        universe
    )


    save_enrichment(
        go,
        paste0(
            label,
            "_GO_BP.tsv"
        )
    )


    save_enrichment(
        kegg,
        paste0(
            label,
            "_KEGG.tsv"
        )
    )


    save_enrichment(
        reactome,
        paste0(
            label,
            "_Reactome.tsv"
        )
    )


    save_dotplot(
        go,
        paste0(
            label,
            "_GO_BP_dotplot.pdf"
        ),
        paste0(
            label,
            " - GO BP"
        )
    )


    save_dotplot(
        kegg,
        paste0(
            label,
            "_KEGG_dotplot.pdf"
        ),
        paste0(
            label,
            " - KEGG"
        )
    )


    save_dotplot(
        reactome,
        paste0(
            label,
            "_Reactome_dotplot.pdf"
        ),
        paste0(
            label,
            " - Reactome"
        )
    )


    invisible(
        list(
            GO = go,
            KEGG = kegg,
            Reactome = reactome
        )
    )
}


# ==============================================================================
# 19. Primary promoter ORA
# ==============================================================================

section(
    "Primary promoter enrichment"
)


primary_group1 <- run_bundle(

    prom_group1_20,

    promoter_universe,

    paste0(
        "primary_promoter_",
        GROUP1_HYPER_LABEL,
        "_abs",
        PRIMARY_DELTA
    )
)


primary_group0 <- run_bundle(

    prom_group0_20,

    promoter_universe,

    paste0(
        "primary_promoter_",
        GROUP0_HYPER_LABEL,
        "_abs",
        PRIMARY_DELTA
    )
)


# ==============================================================================
# 20. Sensitivity promoter ORA
# ==============================================================================

section(
    "Sensitivity promoter enrichment"
)


sensitivity_group1 <- run_bundle(

    prom_group1_15,

    promoter_universe,

    paste0(
        "sensitivity_promoter_",
        GROUP1_HYPER_LABEL,
        "_abs",
        SENSITIVITY_DELTA
    )
)


sensitivity_group0 <- run_bundle(

    prom_group0_15,

    promoter_universe,

    paste0(
        "sensitivity_promoter_",
        GROUP0_HYPER_LABEL,
        "_abs",
        SENSITIVITY_DELTA
    )
)


# ==============================================================================
# 21. Secondary nearest-TSS ORA
# ==============================================================================

section(
    "Nearest-TSS enrichment"
)


nearest_group1 <- run_bundle(

    near_group1_20,

    nearest_universe,

    paste0(
        "secondary_nearest10kb_",
        GROUP1_HYPER_LABEL,
        "_abs",
        PRIMARY_DELTA
    )
)


nearest_group0 <- run_bundle(

    near_group0_20,

    nearest_universe,

    paste0(
        "secondary_nearest10kb_",
        GROUP0_HYPER_LABEL,
        "_abs",
        PRIMARY_DELTA
    )
)


# ==============================================================================
# 22. Ranked promoter-level score
# ==============================================================================

section(
    "Constructing ranked promoter gene statistic"
)

promoter_gene_scores <- promoter_stats[
    ,
    .(
        n_tested_promoter_tiles = uniqueN(tile_id),

        median_meth_diff = median(
            meth.diff,
            na.rm = TRUE
        ),

        mean_meth_diff = mean(
            meth.diff,
            na.rm = TRUE
        ),

        min_qvalue = min(
            qvalue,
            na.rm = TRUE
        ),

        max_abs_meth_diff = max(
            abs(meth.diff),
            na.rm = TRUE
        )
    ),
    by = gene_id
]


data.table::setorder(
    promoter_gene_scores,
    -median_meth_diff,
    -mean_meth_diff,
    gene_id
)


# Add a very small deterministic tie-breaker.
#
# The scale is intentionally tiny relative to methylation differences,
# so it only affects exact or near-exact ties.

tie_rank <- frank(
    -promoter_gene_scores$mean_meth_diff,
    ties.method = "first"
)

epsilon <- 1e-10

promoter_gene_scores[
    ,
    gsea_score :=
        median_meth_diff +
        epsilon * tie_rank
]


gene_rank <- promoter_gene_scores$median_meth_diff


names(
    gene_rank
) <- promoter_gene_scores$gene_id


gene_rank <- gene_rank[
    is.finite(
        gene_rank
    )
]


gene_rank <- sort(
    gene_rank,
    decreasing = TRUE
)


# ==============================================================================
# 23. GSEA
# ==============================================================================

section(
    "Ranked promoter GSEA"
)


gsea_go <- tryCatch(

    clusterProfiler::gseGO(

        geneList = gene_rank,

        OrgDb = org.Hs.eg.db,

        keyType = "ENTREZID",

        ont = "BP",

        minGSSize = MIN_GS_SIZE,

        maxGSSize = MAX_GS_SIZE,

        pvalueCutoff = 0.05,

        pAdjustMethod = "BH",
        
        nPermSimple = 10000,

        verbose = FALSE
    ),

    error = function(e) {

        warning(
            "GO GSEA failed: ",
            conditionMessage(e)
        )

        NULL
    }
)


gsea_kegg <- tryCatch(

    clusterProfiler::gseKEGG(

        geneList = gene_rank,

        organism = "hsa",

        keyType = "ncbi-geneid",

        minGSSize = MIN_GS_SIZE,

        maxGSSize = MAX_GS_SIZE,

        pvalueCutoff = 0.05,

        pAdjustMethod = "BH",
        
        nPermSimple = 10000,
        
        verbose = FALSE
    ),

    error = function(e) {

        warning(
            "KEGG GSEA failed: ",
            conditionMessage(e)
        )

        NULL
    }
)


gsea_reactome <- tryCatch(

    ReactomePA::gsePathway(

        geneList = gene_rank,

        organism = "human",

        minGSSize = MIN_GS_SIZE,

        maxGSSize = MAX_GS_SIZE,

        pvalueCutoff = 0.05,

        pAdjustMethod = "BH",     
        
        nPermSimple = 10000,

        verbose = FALSE
    ),

    error = function(e) {

        warning(
            "Reactome GSEA failed: ",
            conditionMessage(e)
        )

        NULL
    }
)


save_enrichment(
    gsea_go,
    "ranked_promoter_GSEA_GO_BP.tsv"
)


save_enrichment(
    gsea_kegg,
    "ranked_promoter_GSEA_KEGG.tsv"
)


save_enrichment(
    gsea_reactome,
    "ranked_promoter_GSEA_Reactome.tsv"
)


gsea_title_suffix <- paste0(
    "\npositive NES = ",
    GROUP1,
    "; negative NES = ",
    GROUP0
)


save_dotplot(
    gsea_go,
    "ranked_promoter_GSEA_GO_BP_dotplot.pdf",
    paste0(
        "Promoter GSEA - GO BP",
        gsea_title_suffix
    )
)


save_dotplot(
    gsea_kegg,
    "ranked_promoter_GSEA_KEGG_dotplot.pdf",
    paste0(
        "Promoter GSEA - KEGG",
        gsea_title_suffix
    )
)


save_dotplot(
    gsea_reactome,
    "ranked_promoter_GSEA_Reactome_dotplot.pdf",
    paste0(
        "Promoter GSEA - Reactome",
        gsea_title_suffix
    )
)


# ==============================================================================
# 24. Save parameters
# ==============================================================================

parameters <- data.table(

    parameter = c(
        "contrast",
        "group0",
        "group1",
        "genome",
        "input_file",
        "q_threshold",
        "primary_delta",
        "sensitivity_delta",
        "promoter_upstream_bp",
        "promoter_downstream_bp",
        "nearest_tss_max_distance_bp",
        "positive_meth_diff",
        "negative_meth_diff",
        "tested_tiles"
    ),

    value = as.character(
        c(
            CONTRAST,
            GROUP0,
            GROUP1,
            GENOME,
            DIFF_FILE,
            Q_THRESHOLD,
            PRIMARY_DELTA,
            SENSITIVITY_DELTA,
            PROMOTER_UPSTREAM,
            PROMOTER_DOWNSTREAM,
            NEAREST_TSS_MAX_DISTANCE,
            paste0(
                GROUP1,
                "_hypermethylated"
            ),
            paste0(
                GROUP0,
                "_hypermethylated"
            ),
            nrow(
                dm
            )
        )
    )
)


data.table::fwrite(
    parameters,
    file.path(
        OUTPUT_DIR,
        "00_analysis_parameters.tsv"
    ),
    sep = "\t"
)


# ==============================================================================
# 25. README
# ==============================================================================

readme <- c(

    paste0(
        "Pathway analysis: ",
        CONTRAST
    ),

    paste(
        rep(
            "=",
            70
        ),
        collapse = ""
    ),

    "",

    paste0(
        "Group 0: ",
        GROUP0
    ),

    paste0(
        "Group 1: ",
        GROUP1
    ),

    "",

    "methylKit convention:",

    paste0(
        "  meth.diff > 0 = ",
        GROUP1,
        " hypermethylated relative to ",
        GROUP0
    ),

    paste0(
        "  meth.diff < 0 = ",
        GROUP0,
        " hypermethylated relative to ",
        GROUP1
    ),

    "",

    paste0(
        "Primary criterion: q <= ",
        Q_THRESHOLD,
        " and |meth.diff| >= ",
        PRIMARY_DELTA
    ),

    paste0(
        "Sensitivity criterion: q <= ",
        Q_THRESHOLD,
        " and |meth.diff| >= ",
        SENSITIVITY_DELTA
    ),

    "",

    paste0(
        "Promoter: TSS - ",
        PROMOTER_UPSTREAM,
        " bp to TSS + ",
        PROMOTER_DOWNSTREAM,
        " bp"
    ),

    paste0(
        "Secondary nearest-TSS maximum distance: ",
        NEAREST_TSS_MAX_DISTANCE,
        " bp"
    ),

    "",

    "ORA background:",

    "  genes represented by tested methylation tiles.",

    "",

    "Ranked GSEA:",

    "  gene score = median meth.diff across all promoter-overlapping tiles.",

    paste0(
        "  positive NES = ",
        GROUP1,
        "-directed methylation"
    ),

    paste0(
        "  negative NES = ",
        GROUP0,
        "-directed methylation"
    )
)


writeLines(
    readme,
    file.path(
        OUTPUT_DIR,
        "README.txt"
    )
)


# ==============================================================================
# 26. Session information
# ==============================================================================

writeLines(
    capture.output(
        sessionInfo()
    ),
    file.path(
        OUTPUT_DIR,
        "sessionInfo.txt"
    )
)


# ==============================================================================
# 27. Final report
# ==============================================================================

section(
    "ANALYSIS COMPLETE"
)


message(
    "Contrast: ",
    CONTRAST
)


message(
    "Group 0: ",
    GROUP0
)


message(
    "Group 1: ",
    GROUP1
)


message(
    "Input:\n  ",
    DIFF_FILE
)


message(
    "\nAnnotation:\n  ",
    file.path(
        ANNOTATION_DIR,
        paste0(
            CONTRAST,
            "_tiled_annotated.tsv.gz"
        )
    )
)


message(
    "\nPathway output:\n  ",
    OUTPUT_DIR
)


message(
    "\nPrimary promoter gene counts:"
)


message(
    "  ",
    GROUP1_HYPER_LABEL,
    ": ",
    length(
        prom_group1_20
    )
)


message(
    "  ",
    GROUP0_HYPER_LABEL,
    ": ",
    length(
        prom_group0_20
    )
)


message(
    "\nGSEA interpretation:"
)


message(
    "  positive NES = ",
    GROUP1
)


message(
    "  negative NES = ",
    GROUP0
)


message(
    "\nDone."
)
