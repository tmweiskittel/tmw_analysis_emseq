#!/usr/bin/env bash

set -euo pipefail


# ============================================================================
# Configuration
# ============================================================================

REPO_PATH="${REPO_PATH:-/home/jupyter/repos/tmw_analysis_emseq}"

PLOT_SCRIPT="${REPO_PATH}/scripts/visualize_qc.R"

# Search beneath this directory for the aggregate QC files.
ANALYSIS_ROOT="${ANALYSIS_ROOT:-/home/jupyter/data/analysis}"

# Local location for generated QC figures.
LOCAL_OUTPUT_ROOT="${LOCAL_OUTPUT_ROOT:-/home/jupyter/data/qc_visualization}"

# Cloud root containing the three contrast directories.
GCS_ROOT="gs://weiskittel-projects1/radnecrosis/diff_methyl_analysis/differential_methylation"

CONTRASTS=(
    "gbmcsf_metcsf"
    "gbmpla_metpla"
    "csf_pla"
)


# ============================================================================
# Validation
# ============================================================================

for program in Rscript gcloud find; do
    if ! command -v "${program}" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: ${program}" >&2
        exit 1
    fi
done

if [[ ! -s "${PLOT_SCRIPT}" ]]; then
    echo "ERROR: QC plotting script not found:" >&2
    echo "  ${PLOT_SCRIPT}" >&2
    exit 1
fi

mkdir -p "${LOCAL_OUTPUT_ROOT}"


# ============================================================================
# Locate a contrast's aggregate QC file
# ============================================================================

find_qc_aggregate() {
    local contrast="$1"
    local matches
    local count

    matches=$(
        find "${ANALYSIS_ROOT}" \
            -type f \
            -path "*${contrast}*" \
            -name "qc_aggregate.tsv" \
            -print
    )

    count=$(
        printf "%s\n" "${matches}" |
            awk 'NF {n++} END {print n + 0}'
    )

    if [[ "${count}" -eq 0 ]]; then
        echo "ERROR: No qc_aggregate.tsv found for ${contrast}" >&2
        return 1
    fi

    if [[ "${count}" -gt 1 ]]; then
        echo "ERROR: Multiple qc_aggregate.tsv files found for ${contrast}:" >&2
        printf "%s\n" "${matches}" >&2
        return 2
    fi

    printf "%s\n" "${matches}"
}


# ============================================================================
# Process contrasts
# ============================================================================

for contrast in "${CONTRASTS[@]}"; do

    echo
    echo "========================================================================"
    echo "QC visualization: ${contrast}"
    echo "========================================================================"

    if ! QC_AGGREGATE=$(
        find_qc_aggregate "${contrast}"
    ); then
        echo "Skipping ${contrast}."
        continue
    fi

    OUTDIR="${LOCAL_OUTPUT_ROOT}/${contrast}"
    LOG="${OUTDIR}/visualize_qc.log"

    GCS_CONTRAST_DIR="${GCS_ROOT}/${contrast}"
    GCS_QC_DIR="${GCS_CONTRAST_DIR}/qc"

    mkdir -p "${OUTDIR}"

    echo "Input:"
    echo "  ${QC_AGGREGATE}"
    echo
    echo "Local output:"
    echo "  ${OUTDIR}"
    echo
    echo "Cloud output:"
    echo "  ${GCS_QC_DIR}"
    echo


    # ========================================================================
    # Generate plots
    # ========================================================================

    echo "Generating QC plots..."

    Rscript "${PLOT_SCRIPT}" \
        --qc "${QC_AGGREGATE}" \
        --outdir "${OUTDIR}" \
        > "${LOG}" 2>&1


    # ========================================================================
    # Validate output
    # ========================================================================

    plot_count=$(
        find "${OUTDIR}" \
            -maxdepth 1 \
            -type f \
            -name '*.png' \
            | wc -l
    )

    if [[ "${plot_count}" -eq 0 ]]; then
        echo "ERROR: No plots generated for ${contrast}." >&2
        echo "See:" >&2
        echo "  ${LOG}" >&2
        exit 1
    fi

    echo "Generated ${plot_count} plots."


    # ========================================================================
    # Upload plots
    # ========================================================================

    echo "Uploading QC plots..."

    gcloud storage cp \
        "${OUTDIR}"/*.png \
        "${GCS_QC_DIR}/"


    # ========================================================================
    # Upload plotting log
    # ========================================================================

    gcloud storage cp \
        "${LOG}" \
        "${GCS_QC_DIR}/visualize_qc.log"


    # ========================================================================
    # Also upload the aggregate QC and statistics if present
    # ========================================================================

    echo "Uploading aggregate QC table..."

    gcloud storage cp \
        "${QC_AGGREGATE}" \
        "${GCS_QC_DIR}/qc_aggregate.tsv"

    QC_DIR=$(dirname "${QC_AGGREGATE}")
    QC_STATS="${QC_DIR}/qc_group_tests.tsv"

    if [[ -s "${QC_STATS}" ]]; then

        echo "Uploading QC group statistics..."

        gcloud storage cp \
            "${QC_STATS}" \
            "${GCS_QC_DIR}/qc_group_tests.tsv"

    fi


    echo
    echo "Completed ${contrast}"
    echo "Uploaded to:"
    echo "  ${GCS_QC_DIR}"

done


echo
echo "========================================================================"
echo "QC visualization complete"
echo "========================================================================"
echo
echo "Local outputs:"
echo "  ${LOCAL_OUTPUT_ROOT}"
echo
echo "Cloud outputs:"
echo "  ${GCS_ROOT}/gbmcsf_metcsf/qc/"
echo "  ${GCS_ROOT}/gbmpla_metpla/qc/"
echo "  ${GCS_ROOT}/csf_pla/qc/"
