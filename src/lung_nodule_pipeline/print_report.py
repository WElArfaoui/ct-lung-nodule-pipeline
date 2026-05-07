"""Human-readable per-study report; also writes nodules_report_with_slice_ranges.csv."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import SimpleITK as sitk
from skimage.measure import label as cc_label

from .geometry import CONNECTIVITY

log = logging.getLogger(__name__)


def slice_ranges(seg_path: Path, min_vol_mm3: float) -> list[dict]:
    img = sitk.ReadImage(str(seg_path))
    arr = sitk.GetArrayFromImage(img).astype(np.uint8)  # (Z, Y, X)
    cc = cc_label(arr, connectivity=CONNECTIVITY)

    sx, sy, sz = img.GetSpacing()
    voxel_vol = sx * sy * sz
    comps: list[dict] = []
    for cid in range(1, cc.max() + 1):
        m = cc == cid
        vol = float(m.sum() * voxel_vol)
        if vol < min_vol_mm3:
            continue
        zs = np.where(m.any(axis=(1, 2)))[0]
        ys = np.where(m.any(axis=(0, 2)))[0]
        xs = np.where(m.any(axis=(0, 1)))[0]
        comps.append({
            "vol": vol,
            "z_min": int(zs.min()) + 1, "z_max": int(zs.max()) + 1,
            "y_min": int(ys.min()) + 1, "y_max": int(ys.max()) + 1,
            "x_min": int(xs.min()) + 1, "x_max": int(xs.max()) + 1,
        })
    comps.sort(key=lambda c: -c["vol"])
    return comps


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reports-dir", required=True, type=Path)
    p.add_argument("--seg-dir", required=True, type=Path)
    p.add_argument("--min-vol-mm3", type=float, default=14.0)
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    df = pd.read_csv(args.reports_dir / "nodules_report_thin.csv")

    rows_csv: list[dict] = []
    for study_id, group in df.groupby("study_id"):
        case_id = group["case_id"].iloc[0]
        spacing = group["spacing_xyz_mm"].iloc[0]
        print("=" * 100)
        print(f"PATIENT/STUDY: {study_id}")
        print(f"  spacing (x/y/z mm): {spacing}   |   recon: {group['recon_variant'].iloc[0]}")
        print(f"  detected nodules:    {len(group)}")
        print()

        seg_path = args.seg_dir / case_id / "lung_nodules.nii.gz"
        ranges = slice_ranges(seg_path, args.min_vol_mm3)
        if len(ranges) != len(group):
            print(f"  [!] mismatch: ranges={len(ranges)} vs report={len(group)}")

        print(f"  {'#':>2}  {'Vol(mm³)':>10}  {'Diam ax(mm)':>12}  {'Diam sg(mm)':>12}  "
              f"{'Slices ax (Z)':>17}  {'Slices sg (X)':>17}  {'Slices co (Y)':>17}")
        print("  " + "-" * 96)

        for (_, row), rng in zip(group.sort_values("nodule_id").iterrows(), ranges):
            ax = f"{rng['z_min']}-{rng['z_max']} ({row['n_slices_axial']} slices)"
            sg = f"{rng['x_min']}-{rng['x_max']} ({row['n_slices_sagittal']} slices)"
            co = f"{rng['y_min']}-{rng['y_max']} ({row['n_slices_coronal']} slices)"
            print(f"  {int(row['nodule_id']):>2}  "
                  f"{row['volume_mm3']:>10.2f}  "
                  f"{row['diam_axial_mm']:>12.2f}  "
                  f"{row['diam_sagittal_mm']:>12.2f}  "
                  f"{ax:>17}  {sg:>17}  {co:>17}")
            rows_csv.append({
                "study_id": study_id,
                "nodule_id": int(row["nodule_id"]),
                "volume_mm3": row["volume_mm3"],
                "diam_axial_mm": row["diam_axial_mm"],
                "diam_sagittal_mm": row["diam_sagittal_mm"],
                "axial_slice_first": rng["z_min"],
                "axial_slice_last":  rng["z_max"],
                "axial_n_slices":    row["n_slices_axial"],
                "sagittal_slice_first": rng["x_min"],
                "sagittal_slice_last":  rng["x_max"],
                "sagittal_n_slices":    row["n_slices_sagittal"],
                "coronal_slice_first": rng["y_min"],
                "coronal_slice_last":  rng["y_max"],
                "coronal_n_slices":    row["n_slices_coronal"],
                "spacing_xyz_mm": row["spacing_xyz_mm"],
            })
        print()

    out = args.reports_dir / "nodules_report_with_slice_ranges.csv"
    pd.DataFrame(rows_csv).to_csv(out, index=False)
    log.info("wrote %s", out)


if __name__ == "__main__":
    main()
