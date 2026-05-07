# lung-nodule-pipeline — end-to-end orchestration
#
# Override paths on the CLI:
#   make all DICOM_ROOT=/data/dicom DATA_DIR=/data/cohort
#
# Each target only re-runs if its sentinel/output is missing.

SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.ONESHELL:

DICOM_ROOT  ?= ./data/dicom
DATA_DIR    ?= ./data
NIFTI_DIR   ?= $(DATA_DIR)/nifti
SEG_DIR     ?= $(DATA_DIR)/seg
REPORTS_DIR ?= $(DATA_DIR)/reports
LOGS_DIR    ?= $(DATA_DIR)/logs
MIN_VOL_MM3 ?= 14.0

PYTHON ?= python

.PHONY: help install convert segment analyze select-thin print-report diameters final-summary all clean

help:
	@echo "Targets:"
	@echo "  install         pip install -e ."
	@echo "  convert         DICOM -> NIfTI (dcm2niix)"
	@echo "  segment         NIfTI -> lung_nodules masks (TotalSegmentator)"
	@echo "  analyze         masks -> nodules_report.csv + study_summary.csv"
	@echo "  select-thin     keep thinnest reconstruction per study"
	@echo "  print-report    human-readable per-study report"
	@echo "  diameters       full diameter set (Lung-RADS, 3D, eq-sphere)"
	@echo "  final-summary   final consolidated CSV with Lung-RADS"
	@echo "  all             run the full pipeline"
	@echo "  clean           remove derived data (NIfTI, seg, reports, logs)"
	@echo ""
	@echo "Variables: DICOM_ROOT, DATA_DIR, NIFTI_DIR, SEG_DIR, REPORTS_DIR, MIN_VOL_MM3"

install:
	$(PYTHON) -m pip install -e .

convert:
	bash scripts/convert_dicom_to_nifti.sh "$(DICOM_ROOT)" "$(NIFTI_DIR)"

segment:
	bash scripts/segment_nodules.sh "$(NIFTI_DIR)" "$(SEG_DIR)"

analyze:
	$(PYTHON) -m lung_nodule_pipeline.analyze_nodules \
	    --seg-dir "$(SEG_DIR)" \
	    --out-csv "$(REPORTS_DIR)/nodules_report.csv" \
	    --min-vol-mm3 $(MIN_VOL_MM3)

select-thin:
	$(PYTHON) -m lung_nodule_pipeline.select_thin_recon \
	    --reports-dir "$(REPORTS_DIR)"

print-report:
	$(PYTHON) -m lung_nodule_pipeline.print_report \
	    --reports-dir "$(REPORTS_DIR)" \
	    --seg-dir "$(SEG_DIR)" \
	    --min-vol-mm3 $(MIN_VOL_MM3)

diameters:
	$(PYTHON) -m lung_nodule_pipeline.compute_diameters \
	    --reports-dir "$(REPORTS_DIR)" \
	    --seg-dir "$(SEG_DIR)" \
	    --min-vol-mm3 $(MIN_VOL_MM3)

final-summary:
	$(PYTHON) -m lung_nodule_pipeline.final_summary \
	    --reports-dir "$(REPORTS_DIR)"

all: convert segment analyze select-thin diameters final-summary

clean:
	rm -rf "$(NIFTI_DIR)" "$(SEG_DIR)" "$(REPORTS_DIR)" "$(LOGS_DIR)"
