# Example data layout

A minimal cohort with two patients, one of them with two reconstructions of the
same study:

```
data/dicom/
├── PAT001/
│   └── 20240115/
│       ├── 1.2.840.113619.2.55.3.604688.111/        ← thin recon (1 mm)
│       │   ├── 0001.dcm
│       │   ├── 0002.dcm
│       │   └── ...
│       └── 1.2.840.113619.2.55.3.604688.222/        ← thick recon (3 mm)
│           └── ...
└── PAT002/
    └── 20240220/
        └── 1.2.840.113619.2.55.3.604688.333/
            └── ...
```

After running `make convert`, `data/nifti/` contains:

```
PAT001_20240115_1.2.840.113619.2.55.3.604688.111.nii.gz
PAT001_20240115_1.2.840.113619.2.55.3.604688.222.nii.gz
PAT002_20240220_1.2.840.113619.2.55.3.604688.333.nii.gz
```

After `make segment`, `data/seg/<case_id>/lung_nodules.nii.gz` exists for each
case.

After `make all`, the reports under `data/reports/` are derived from the
**thinnest** reconstruction per `study_id` (`PAT001_20240115_*` keeps the 1 mm
variant; `PAT002_20240220_*` is the only one available).

## Reconstruction variants

If your cohort uses a single `case_id` per study with no thin/thick variants,
the pipeline still works — each study just has one row, with
`recon_variant = base`.

To mark a variant as alternative, append a literal `a` to the directory name
that becomes the `case_id` (matched by the regex
`^(?P<study_id>.+?_\d{8}_[\d.]+?)(?P<variant>a?)$`). Customize that regex in
`src/lung_nodule_pipeline/geometry.py` if your naming convention differs.
