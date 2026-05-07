"""Analyze TotalSegmentator lung_nodules masks → per-nodule and per-study CSVs."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import SimpleITK as sitk
from skimage.measure import label as cc_label

from .geometry import (
    CONNECTIVITY,
    measure_view,
    parse_case_id,
)

log = logging.getLogger(__name__)


def analyze_case(seg_path: Path, case_id: str, min_vol_mm3: float) -> list[dict]:
    img = sitk.ReadImage(str(seg_path))
    arr = sitk.GetArrayFromImage(img).astype(np.uint8)
    sx, sy, sz = img.GetSpacing()
    spacing_zyx = (sz, sy, sx)
    voxel_vol = sx * sy * sz

    study_id, variant = parse_case_id(case_id)
    cc = cc_label(arr, connectivity=CONNECTIVITY)

    rows: list[dict] = []
    for cid in range(1, cc.max() + 1):
        m = cc == cid
        n_vox = int(m.sum())
        vol = n_vox * voxel_vol
        if vol < min_vol_mm3:
            continue

        ax = measure_view(m, "axial", spacing_zyx)
        sg = measure_view(m, "sagittal", spacing_zyx)
        co = measure_view(m, "coronal", spacing_zyx)

        zyx = np.argwhere(m).mean(axis=0)
        cx, cy, cz = img.TransformContinuousIndexToPhysicalPoint(
            (zyx[2], zyx[1], zyx[0])
        )

        rows.append({
            "case_id": case_id,
            "study_id": study_id,
            "recon_variant": variant,
            "nodule_id": int(cid),
            "volume_mm3": round(vol, 2),
            "diam_axial_mm": round(ax["long_mm"], 2),
            "diam_sagittal_mm": round(sg["long_mm"], 2),
            "n_slices_axial": ax["n_slices"],
            "n_slices_sagittal": sg["n_slices"],
            "n_slices_coronal": co["n_slices"],
            "extent_z_mm": round(ax["n_slices"] * sz, 2),
            "centroid_x_mm": round(cx, 2),
            "centroid_y_mm": round(cy, 2),
            "centroid_z_mm": round(cz, 2),
            "spacing_xyz_mm": f"{sx:.2f}/{sy:.2f}/{sz:.2f}",
        })

    rows.sort(key=lambda r: -r["volume_mm3"])
    for new_id, r in enumerate(rows, start=1):
        r["nodule_id"] = new_id
    return rows


def build_study_summary(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    grp = df.groupby(["study_id", "recon_variant"])
    summary = grp.agg(
        n_nodules=("nodule_id", "count"),
        n_nodules_ge_6mm=("diam_axial_mm", lambda s: int((s >= 6).sum())),
        max_diam_axial_mm=("diam_axial_mm", "max"),
        max_diam_sagittal_mm=("diam_sagittal_mm", "max"),
        total_volume_mm3=("volume_mm3", "sum"),
        spacing=("spacing_xyz_mm", "first"),
    ).reset_index()
    pivot = summary.pivot(index="study_id", columns="recon_variant").reset_index()
    pivot.columns = ["study_id"] + [f"{m}__{v}" for m, v in pivot.columns[1:]]
    return pivot


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--seg-dir", required=True, type=Path,
                   help="Directory with <case_id>/lung_nodules.nii.gz")
    p.add_argument("--out-csv", required=True, type=Path,
                   help="Path for nodules_report.csv (study_summary.csv written alongside)")
    p.add_argument("--min-vol-mm3", type=float, default=14.0,
                   help="Minimum component volume to keep (default: 14)")
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    all_rows: list[dict] = []
    for cd in sorted(d for d in args.seg_dir.iterdir() if d.is_dir()):
        seg = cd / "lung_nodules.nii.gz"
        if not seg.exists():
            log.warning("missing %s", seg)
            continue
        try:
            rows = analyze_case(seg, cd.name, args.min_vol_mm3)
            all_rows.extend(rows)
            log.info("%s -> %d nodules", cd.name, len(rows))
        except Exception:
            log.exception("failed: %s", cd.name)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(all_rows)
    df.to_csv(args.out_csv, index=False)
    n_cases = df["case_id"].nunique() if not df.empty else 0
    log.info("wrote %s (%d nodules, %d cases)", args.out_csv, len(df), n_cases)

    summary = build_study_summary(df)
    summary_path = args.out_csv.with_name("study_summary.csv")
    summary.to_csv(summary_path, index=False)
    log.info("wrote %s (%d studies)", summary_path, len(summary))


if __name__ == "__main__":
    main()
