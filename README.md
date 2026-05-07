# lung-nodule-pipeline

Reproducible pipeline for **lung-nodule segmentation and quantification** on chest CT.

The pipeline takes a DICOM cohort, converts it to NIfTI, segments lung nodules
with [TotalSegmentator](https://github.com/wasserth/TotalSegmentator)
(`-ta lung_nodules`), and produces per-nodule metrics (volume, axial / sagittal /
coronal long+short axes, 3D Feret, equivalent-sphere diameter) plus an approximate
**Lung-RADS** category and patient-level summaries.

> Disclaimer — research / educational use only. Not a medical device. The
> Lung-RADS labels are an automated approximation and **not** a clinical diagnosis.

---

## Pipeline at a glance

```
DICOM cohort
    │  (1) scripts/convert_dicom_to_nifti.sh   — dcm2niix
    ▼
NIfTI volumes
    │  (2) scripts/segment_nodules.sh          — TotalSegmentator -ta lung_nodules
    ▼
Per-case lung_nodules.nii.gz
    │  (3) lnp-analyze         — connected components + per-view long axis
    │  (4) lnp-select-thin     — keep thinnest reconstruction per study
    │  (5) lnp-print-report    — human-readable per-study report (+ slice ranges CSV)
    │  (6) lnp-diameters       — full Lung-RADS / 3D / eq-sphere diameters
    │  (7) lnp-final-summary   — final consolidated CSV with Lung-RADS
    ▼
reports/*.csv
```

---

## Requirements

- Python ≥ 3.10
- [`dcm2niix`](https://github.com/rordenlab/dcm2niix) on `PATH`
- [`TotalSegmentator`](https://github.com/wasserth/TotalSegmentator) on `PATH`
  (academic license required for the `lung_nodules` task: run
  `totalseg_set_license -l aca_XXXXX` once)
- A CUDA-capable GPU recommended for segmentation (CPU works but is slow)
- `make` for the orchestrator (optional)

Python dependencies are listed in `requirements.txt` / `pyproject.toml`:
`numpy`, `pandas`, `SimpleITK`, `scikit-image`, `TotalSegmentator`.

---

## Install

```bash
git clone https://github.com/<your-org>/lung-nodule-pipeline.git
cd lung-nodule-pipeline

python -m venv .venv && source .venv/bin/activate
pip install -e .
# or, for a minimal install:  pip install -r requirements.txt
```

This exposes the following CLI entry points:

| Command             | Module                                       |
| ------------------- | -------------------------------------------- |
| `lnp-analyze`       | `lung_nodule_pipeline.analyze_nodules`       |
| `lnp-select-thin`   | `lung_nodule_pipeline.select_thin_recon`     |
| `lnp-print-report`  | `lung_nodule_pipeline.print_report`          |
| `lnp-diameters`     | `lung_nodule_pipeline.compute_diameters`     |
| `lnp-final-summary` | `lung_nodule_pipeline.final_summary`         |

---

## Expected DICOM layout

```
DICOM_ROOT/
├── <PatientID>/
│   └── <StudyDate>/
│       ├── <SeriesUID>/      ← thin recon
│       │   ├── 0001.dcm
│       │   ├── 0002.dcm
│       │   └── ...
│       └── <SeriesUID2>/     ← alternative recon (optional)
│           └── ...
└── ...
```

Each series directory at depth 3 that contains at least one `.dcm` is converted
to one NIfTI named `<PatientID>_<StudyDate>_<SeriesUID>.nii.gz`. A trailing `a`
on a `case_id` is treated as an alternative reconstruction of the same study —
`lnp-select-thin` keeps the variant with the smallest z-spacing.

See [`examples/data_layout.md`](examples/data_layout.md) for a worked example.

---

## Quick start (Makefile)

```bash
# Convert DICOM → NIfTI → segmentation → CSV reports
make all DICOM_ROOT=/path/to/dicom DATA_DIR=/path/to/output

# Or step by step
make convert       DICOM_ROOT=/path/to/dicom DATA_DIR=/path/to/output
make segment       DATA_DIR=/path/to/output
make analyze       DATA_DIR=/path/to/output
make select-thin   DATA_DIR=/path/to/output
make print-report  DATA_DIR=/path/to/output
make diameters     DATA_DIR=/path/to/output
make final-summary DATA_DIR=/path/to/output
```

Variables (override on the command line):

| Variable      | Default               | Description                            |
| ------------- | --------------------- | -------------------------------------- |
| `DICOM_ROOT`  | `./data/dicom`        | Source DICOM cohort                    |
| `DATA_DIR`    | `./data`              | Root for derived outputs               |
| `NIFTI_DIR`   | `$DATA_DIR/nifti`     | Converted NIfTI                        |
| `SEG_DIR`     | `$DATA_DIR/seg`       | TotalSegmentator output                |
| `REPORTS_DIR` | `$DATA_DIR/reports`   | CSV outputs                            |
| `MIN_VOL_MM3` | `14.0`                | Minimum nodule volume to keep          |

Each step skips work it has already done (idempotent), so re-running `make all`
is cheap.

---

## Manual usage (without make)

```bash
bash scripts/convert_dicom_to_nifti.sh $DICOM_ROOT $NIFTI_DIR
bash scripts/segment_nodules.sh        $NIFTI_DIR  $SEG_DIR

lnp-analyze       --seg-dir $SEG_DIR --out-csv $REPORTS_DIR/nodules_report.csv --min-vol-mm3 14
lnp-select-thin   --reports-dir $REPORTS_DIR
lnp-print-report  --reports-dir $REPORTS_DIR --seg-dir $SEG_DIR
lnp-diameters     --reports-dir $REPORTS_DIR --seg-dir $SEG_DIR
lnp-final-summary --reports-dir $REPORTS_DIR
```

---

## Outputs

After a full run, `$REPORTS_DIR/` contains:

| File                                    | Description                                                |
| --------------------------------------- | ---------------------------------------------------------- |
| `nodules_report.csv`                    | One row per nodule, all reconstructions                    |
| `study_summary.csv`                     | Pivoted per-study summary (base / alt_a)                   |
| `nodules_report_thin.csv`               | Same as above, thinnest recon only                         |
| `study_summary_thin.csv`                | Per-study summary, thin only                               |
| `nodules_report_with_slice_ranges.csv`  | Nodules + slice-range first/last per view                  |
| `nodules_full_diameters.csv`            | Long+short axes, 3D max, eq-sphere                         |
| `nodules_FINAL_summary.csv`             | Final per-nodule CSV with Lung-RADS + QC flags             |
| `patient_summary_FINAL.csv`             | Per-study summary with highest Lung-RADS                   |

### Lung-RADS approximation

| Mean axial diameter (mm) | Volume (mm³) | Category        |
| ------------------------ | ------------ | --------------- |
| `≥ 30` or vol > 30 000   | —            | `4X (mass)`     |
| `15–30`                  | —            | `4B`            |
| `8–15`                   | —            | `4A`            |
| `6–8`                    | —            | `3`             |
| `< 6`                    | —            | `2`             |

`qc_flag = ANISOTROPIC_review` when the sagittal or coronal long axis is
≥ 1.5× the axial long axis (likely thick-slice artifact — review manually).

---

## Project layout

```
.
├── Makefile                       end-to-end orchestrator
├── pyproject.toml                 packaging + entry-points + ruff config
├── requirements.txt               minimal pip install
├── scripts/
│   ├── convert_dicom_to_nifti.sh  dcm2niix wrapper
│   └── segment_nodules.sh         TotalSegmentator wrapper
├── src/lung_nodule_pipeline/
│   ├── __init__.py
│   ├── geometry.py                shared 2D/3D Feret + view helpers
│   ├── analyze_nodules.py         step 3 — per-nodule CSV
│   ├── select_thin_recon.py       step 4 — thin-recon filter
│   ├── print_report.py            step 5 — human-readable + slice ranges
│   ├── compute_diameters.py       step 6 — full diameter set
│   └── final_summary.py           step 7 — Lung-RADS final CSV
└── examples/
    └── data_layout.md
```

---

## Privacy / data handling

- **Never commit DICOM, NIfTI, or segmentation files.** `.gitignore` excludes
  `data/`, `dicom/`, `nifti/`, `seg/`, `reports/`, `*.nii*`, `*.dcm`.
- Outputs use the case-ID derived from the directory tree
  (`<PatientID>_<StudyDate>_<SeriesUID>`). Anonymize before committing or
  sharing CSVs.

---

## Citation

If you use TotalSegmentator through this pipeline, please cite:

> Wasserthal et al., *TotalSegmentator: Robust Segmentation of 104 Anatomic
> Structures in CT Images*, Radiology AI 2023.

---

## License

MIT — see [`LICENSE`](LICENSE).
