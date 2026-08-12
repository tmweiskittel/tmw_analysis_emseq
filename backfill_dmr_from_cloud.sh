#!/usr/bin/env bash

set -euo pipefail

D_OUT="/home/jupyter/data/analysis"
BUCKET="weiskittel-projects1"
PREFIX="radnecrosis/diff_methyl_analysis"

GCS_ROOT="gs://${BUCKET}/${PREFIX}/differential_methylation"

mkdir -p \
  "${D_OUT}/dmr/diff" \
  "${D_OUT}/dmr/annotation"

echo "Discovering experiments in ${GCS_ROOT} ..."

EXPERIMENTS=$(
  gcloud storage ls "${GCS_ROOT}/" \
    | sed 's:/$::' \
    | awk -F/ '{print $NF}'
)

for EXPERIMENT in ${EXPERIMENTS}; do

  echo
  echo "============================================================"
  echo "Backfilling ${EXPERIMENT}"
  echo "============================================================"

  SRC="${GCS_ROOT}/${EXPERIMENT}"

  declare -a DOWNLOADS=(

    "methylBase_${EXPERIMENT}.txt.bgz|${D_OUT}/dmr/diff/methylBase_${EXPERIMENT}.txt.bgz"

    "methylBase_${EXPERIMENT}.txt.bgz.tbi|${D_OUT}/dmr/diff/methylBase_${EXPERIMENT}.txt.bgz.tbi"

    "methylDiff_${EXPERIMENT}.txt.bgz|${D_OUT}/dmr/diff/methylDiff_${EXPERIMENT}.txt.bgz"

    "methylDiff_${EXPERIMENT}.txt.bgz.tbi|${D_OUT}/dmr/diff/methylDiff_${EXPERIMENT}.txt.bgz.tbi"

    "methylBase_${EXPERIMENT}.tiled.txt.bgz|${D_OUT}/dmr/diff/methylBase_${EXPERIMENT}.tiled.txt.bgz"

    "methylBase_${EXPERIMENT}.tiled.txt.bgz.tbi|${D_OUT}/dmr/diff/methylBase_${EXPERIMENT}.tiled.txt.bgz.tbi"

    "methylDiff_${EXPERIMENT}.tiled.txt.bgz|${D_OUT}/dmr/diff/methylDiff_${EXPERIMENT}.tiled.txt.bgz"

    "methylDiff_${EXPERIMENT}.tiled.txt.bgz.tbi|${D_OUT}/dmr/diff/methylDiff_${EXPERIMENT}.tiled.txt.bgz.tbi"

    "${EXPERIMENT}_pos_meth.tsv|${D_OUT}/dmr/diff/${EXPERIMENT}_pos_meth.tsv"

    "${EXPERIMENT}_annotated.tsv|${D_OUT}/dmr/annotation/${EXPERIMENT}_annotated.tsv"
  )

  for ITEM in "${DOWNLOADS[@]}"; do

    REMOTE_NAME="${ITEM%%|*}"
    LOCAL_PATH="${ITEM#*|}"

    if [[ -s "${LOCAL_PATH}" ]]; then
      echo "Already exists: ${LOCAL_PATH}"
      continue
    fi

    echo "Downloading ${REMOTE_NAME}"

    mkdir -p "$(dirname "${LOCAL_PATH}")"

    if gcloud storage ls "${SRC}/${REMOTE_NAME}" >/dev/null 2>&1; then

      gcloud storage cp \
        "${SRC}/${REMOTE_NAME}" \
        "${LOCAL_PATH}"

    else
      echo "WARNING: not found in cloud:"
      echo "  ${SRC}/${REMOTE_NAME}"
    fi

  done

done

echo
echo "Backfill complete."
