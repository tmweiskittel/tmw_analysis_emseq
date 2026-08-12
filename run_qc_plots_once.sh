#!/usr/bin/env bash

set -euo pipefail


# ============================================================================
# Configuration
# ============================================================================

REPO="/home/jupyter/repos/tmw_analysis_emseq"

PLOT_SCRIPT="${REPO}/scripts/visualize_qc.R"

LOCAL_QC_ROOT="/home/jupyter/data/analysis/qc"
LOCAL_PLOT_ROOT="/home/jupyter/data/analysis/qc_plots"

GCS_ROOT="gs://weiskittel-projects1/radnecrosis/diff_methyl_analysis/differential_methylation"

CONTRASTS=(
    "gbmcsf_metcsf"
    "gbmpla_metpla"
    "csf_pla"
)


# ============================================================================
# Validation
# ============================================================================

for program in Rscript gcloud; do
    if ! command -v "${program}" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: ${program}" >&2
        exit 1
    fi
done

if [[ ! -s "${PLOT_SCRIPT}" ]]; then
    echo "ERROR: Plot script not found:"
    echo "  ${PLOT_SCRIPT}"
    exit 1
fi


# ============================================================================
# Run each contrast
# ============================================================================

for contrast in "${CONTRASTS[@]}"; do

    echo
    echo "======================================================================"
    echo "QC: ${contrast}"
    echo "======================================================================"

    QC_AGGREGATE="${LOCAL_QC_ROOT}/${contrast}/qc_aggregate.tsv"
    QC_STATS="${LOCAL_QC_ROOT}/${contrast}/qc_group_tests.tsv"

    OUTDIR="${LOCAL_PLOT_ROOT}/${contrast}"

    GCS_CONTRAST="${GCS_ROOT}/${contrast}"
    GCS_QC="${GCS_CONTRAST}/qc"

    mkdir -p "${OUTDIR}"

    LOG="${OUTDIR}/visualize_qc.log"


    # ------------------------------------------------------------------------
    # Check aggregate
    # ------------------------------------------------------------------------

    if [[ ! -s "${QC_AGGREGATE}" ]]; then
        echo "ERROR: Missing aggregate QC file:"
        echo "  ${QC_AGGREGATE}"
        exit 1
    fi


    # ------------------------------------------------------------------------
    # Generate plots
    # ------------------------------------------------------------------------

    echo "Generating plots..."

    Rscript "${PLOT_SCRIPT}" \
        --qc "${QC_AGGREGATE}" \
        --outdir "${OUTDIR}" \
        > "${LOG}" 2>&1


    # ------------------------------------------------------------------------
    # Count plots
    # ------------------------------------------------------------------------

    plot_count=$(
        find "${OUTDIR}" \
            -maxdepth 1 \
            -type f \
            -name '*.png' \
            | wc -l
    )

    if [[ "${plot_count}" -eq 0 ]]; then
        echo "ERROR: No plots were generated."
        echo
        echo "Log:"
        cat "${LOG}"
        exit 1
    fi

    echo "Generated ${plot_count} plots."


    # ------------------------------------------------------------------------
    # Upload plots
    # ------------------------------------------------------------------------

    echo "Uploading plots to:"
    echo "  ${GCS_QC}"

    gcloud storage cp \
        "${OUTDIR}"/*.png \
        "${GCS_QC}/"


    # ------------------------------------------------------------------------
    # Upload aggregate table
    # ------------------------------------------------------------------------

    gcloud storage cp \
        "${QC_AGGREGATE}" \
        "${GCS_QC}/qc_aggregate.tsv"


    # ------------------------------------------------------------------------
    # Upload group statistics
    # ------------------------------------------------------------------------

    if [[ -s "${QC_STATS}" ]]; then

        gcloud storage cp \
            "${QC_STATS}" \
            "${GCS_QC}/qc_group_tests.tsv"

    fi


    # ------------------------------------------------------------------------
    # Upload log
    # ------------------------------------------------------------------------

    gcloud storage cp \
        "${LOG}" \
        "${GCS_QC}/visualize_qc.log"


    echo
    echo "Finished ${contrast}"

done


echo
echo "======================================================================"
echo "QC plotting complete"
echo "======================================================================"
echo
echo "Local:"
echo "  ${LOCAL_PLOT_ROOT}"
echo
echo "Cloud:"
echo "  ${GCS_ROOT}/gbmcsf_metcsf/qc/"
echo "  ${GCS_ROOT}/gbmpla_metpla/qc/"
echo "  ${GCS_ROOT}/csf_pla/qc/"
