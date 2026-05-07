#!/usr/bin/env bash
# segment_nodules.sh
#
# Run TotalSegmentator (task=lung_nodules) on every NIfTI in NIFTI_DIR.
# License (free for academic use): totalseg_set_license -l aca_XXXXX  (run once)
#
# Usage: segment_nodules.sh <NIFTI_DIR> <SEG_DIR> [--device gpu|cpu]

set -euo pipefail

NIFTI_DIR=${1:?"usage: $0 <NIFTI_DIR> <SEG_DIR> [--device gpu|cpu]"}
SEG_DIR=${2:?"usage: $0 <NIFTI_DIR> <SEG_DIR> [--device gpu|cpu]"}
DEVICE=${3:-gpu}
DEVICE=${DEVICE#--device}
DEVICE=${DEVICE:-gpu}

LOG_DIR="${SEG_DIR}/../logs"
mkdir -p "$SEG_DIR" "$LOG_DIR"
LOG="${LOG_DIR}/segment_nodules.log"

if ! command -v TotalSegmentator >/dev/null 2>&1; then
    echo "[error] TotalSegmentator not found in PATH" >&2
    exit 127
fi

shopt -s nullglob
for nii in "$NIFTI_DIR"/*.nii.gz; do
    case_id=$(basename "$nii" .nii.gz)
    out="${SEG_DIR}/${case_id}"
    target="${out}/lung_nodules.nii.gz"

    if [[ -f "$target" ]]; then
        echo "[skip] $case_id" | tee -a "$LOG"
        continue
    fi

    mkdir -p "$out"
    echo "[run]  $case_id" | tee -a "$LOG"

    TotalSegmentator -i "$nii" -o "$out" -ta lung_nodules \
        --device "$DEVICE" --quiet >> "$LOG" 2>&1 \
        || echo "[fail] $case_id" | tee -a "$LOG"
done

echo "[done] segmentations in: $SEG_DIR"
