#!/usr/bin/env bash
# convert_dicom_to_nifti.sh
#
# Convert a DICOM cohort to NIfTI using dcm2niix.
#
# Expected layout: DICOM_ROOT/<PatientID>/<StudyDate>/<SeriesUID>/*.dcm
# Output filename: <PatientID>_<StudyDate>_<SeriesUID>.nii.gz
#
# Usage: convert_dicom_to_nifti.sh <DICOM_ROOT> <NIFTI_DIR>

set -euo pipefail

DICOM_ROOT=${1:?"usage: $0 <DICOM_ROOT> <NIFTI_DIR>"}
NIFTI_DIR=${2:?"usage: $0 <DICOM_ROOT> <NIFTI_DIR>"}
LOG_DIR="${NIFTI_DIR}/../logs"

mkdir -p "$NIFTI_DIR" "$LOG_DIR"
LOG="${LOG_DIR}/convert_dicom_to_nifti.log"

if ! command -v dcm2niix >/dev/null 2>&1; then
    echo "[error] dcm2niix not found in PATH" >&2
    exit 127
fi

# Only series directories at depth 3 that contain at least one .dcm
find "$DICOM_ROOT" -mindepth 3 -maxdepth 3 -type d | while read -r series_dir; do
    if ! compgen -G "${series_dir}/*.dcm" > /dev/null; then
        continue
    fi

    rel=$(realpath --relative-to="$DICOM_ROOT" "$series_dir")
    case_id=$(echo "$rel" | tr '/' '_')
    out="${NIFTI_DIR}/${case_id}.nii.gz"

    if [[ -f "$out" ]]; then
        echo "[skip] $case_id" | tee -a "$LOG"
        continue
    fi

    echo "[run]  $case_id" | tee -a "$LOG"
    dcm2niix -z y -m y -i y -f "$case_id" -o "$NIFTI_DIR" "$series_dir" \
        >> "$LOG" 2>&1 || echo "[fail] $case_id" | tee -a "$LOG"
done

echo "[done] NIfTI files in: $NIFTI_DIR"
ls -1 "$NIFTI_DIR"/*.nii.gz 2>/dev/null | wc -l | xargs -I {} echo "[count] {} series converted"
