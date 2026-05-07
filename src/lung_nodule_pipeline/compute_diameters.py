"""Full diameter set: long+short axial (Lung-RADS), sagittal, coronal, 3D max, eq-sphere."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import SimpleITK as sitk
from skimage.measure import label as cc_label

from .geometry import CONNECTIVITY, feret_3d_mm, measure_view

log = logging.getLogger(__name__)


def analyze(case_id: str, study_id: str, variant: str,
            seg_dir: Path, min_vol_mm3: float) -> list[dict]:
    seg_path = seg_dir / case_id / "lung_nodules.nii.gz"
    img = sitk.ReadImage(str(seg_path))
    arr = sitk.GetArrayFromImage(img).astype(np.uint8)
    sx, sy, sz = img.GetSpacing()
    spacing_zyx = (sz, sy, sx)
    voxel_vol = sx * sy * sz

    cc = cc_label(arr, connectivity=CONNECTIVITY)
    comps: list[dict] = []
    for cid in range(1, cc.max() + 1):
        m = cc == cid
        vol = float(m.sum() * voxel_vol)
        if vol < min_vol_mm3:
            continue

        ax = measure_view(m, "axial", spacing_zyx)
        sg = measure_view(m, "sagittal", spacing_zyx)
        co = measure_view(m, "coronal", spacing_zyx)
        d3d = feret_3d_mm(m, spacing_zyx)
        d_eq = (6.0 * vol / np.pi) ** (1 / 3)

        zyx = np.argwhere(m).mean(axis=0)
        cx, cy, cz = img.TransformContinuousIndexToPhysicalPoint(
            (zyx[2], zyx[1], zyx[0])
        )

        comps.append({
            "study_id": study_id,
            "case_id": case_id,
            "recon_variant": variant,
            "volume_mm3": round(vol, 2),
            "diam_long_axial_mm":  round(ax["long_mm"], 2),
            "diam_short_axial_mm": round(ax["short_mm"], 2),
            "diam_mean_axial_mm":  round((ax["long_mm"] + ax["short_mm"]) / 2, 2),
            "diam_long_sagittal_mm":  round(sg["long_mm"], 2),
            "diam_short_sagittal_mm": round(sg["short_mm"], 2),
            "diam_long_coronal_mm":   round(co["long_mm"], 2),
            "diam_short_coronal_mm":  round(co["short_mm"], 2),
            "diam_3d_max_mm":     round(d3d, 2),
            "diam_eq_sphere_mm":  round(d_eq, 2),
            "axial_first":    ax["first"], "axial_last":    ax["last"], "axial_n":    ax["n_slices"],
            "sagittal_first": sg["first"], "sagittal_last": sg["last"], "sagittal_n": sg["n_slices"],
            "coronal_first":  co["first"], "coronal_last":  co["last"], "coronal_n":  co["n_slices"],
            "centroid_x_mm": round(cx, 2),
            "centroid_y_mm": round(cy, 2),
            "centroid_z_mm": round(cz, 2),
            "spacing_xyz_mm": f"{sx:.2f}/{sy:.2f}/{sz:.2f}",
        })

    comps.sort(key=lambda c: -c["volume_mm3"])
    for i, c in enumerate(comps, start=1):
        c["nodule_id"] = i
    return comps


COLUMN_ORDER = (
    ["study_id", "case_id", "recon_variant", "nodule_id", "volume_mm3"]
    + ["diam_long_axial_mm", "diam_short_axial_mm", "diam_mean_axial_mm",
       "diam_long_sagittal_mm", "diam_short_sagittal_mm",
       "diam_long_coronal_mm", "diam_short_coronal_mm",
       "diam_3d_max_mm", "diam_eq_sphere_mm"]
    + ["axial_first", "axial_last", "axial_n",
       "sagittal_first", "sagittal_last", "sagittal_n",
       "coronal_first", "coronal_last", "coronal_n"]
    + ["centroid_x_mm", "centroid_y_mm", "centroid_z_mm", "spacing_xyz_mm"]
)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reports-dir", required=True, type=Path)
    p.add_argument("--seg-dir", required=True, type=Path)
    p.add_argument("--min-vol-mm3", type=float, default=14.0)
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    thin_csv = args.reports_dir / "nodules_report_thin.csv"
    df_thin = pd.read_csv(thin_csv)
    chosen = df_thin[["study_id", "case_id", "recon_variant"]].drop_duplicates()

    all_rows: list[dict] = []
    for _, row in chosen.iterrows():
        rows = analyze(row["case_id"], row["study_id"], row["recon_variant"],
                       args.seg_dir, args.min_vol_mm3)
        all_rows.extend(rows)

    df = pd.DataFrame(all_rows)[list(COLUMN_ORDER)]
    out = args.reports_dir / "nodules_full_diameters.csv"
    df.to_csv(out, index=False)
    log.info("wrote %s: %d nodules / %d studies",
             out, len(df), df["study_id"].nunique())

    for study_id, g in df.groupby("study_id"):
        print()
        print("=" * 110)
        print(f"  {study_id}   spacing={g['spacing_xyz_mm'].iloc[0]}   "
              f"recon={g['recon_variant'].iloc[0]}")
        print("=" * 110)
        print(f"  {'#':>2}  {'Vol mm3':>10}  {'Long ax':>7}  {'Short ax':>8}  {'Mean ax':>7}  "
              f"{'Long sg':>7}  {'Long co':>7}  {'3D max':>7}  {'EqSphere':>8}  "
              f"{'Slices ax':>13}  {'Slices sg':>13}  {'Slices co':>13}")
        for _, r in g.sort_values("nodule_id").iterrows():
            ax = f"{r['axial_first']}-{r['axial_last']} ({r['axial_n']})"
            sg = f"{r['sagittal_first']}-{r['sagittal_last']} ({r['sagittal_n']})"
            co = f"{r['coronal_first']}-{r['coronal_last']} ({r['coronal_n']})"
            print(f"  {int(r['nodule_id']):>2}  {r['volume_mm3']:>10.2f}  "
                  f"{r['diam_long_axial_mm']:>7.2f}  {r['diam_short_axial_mm']:>8.2f}  "
                  f"{r['diam_mean_axial_mm']:>7.2f}  "
                  f"{r['diam_long_sagittal_mm']:>7.2f}  {r['diam_long_coronal_mm']:>7.2f}  "
                  f"{r['diam_3d_max_mm']:>7.2f}  {r['diam_eq_sphere_mm']:>8.2f}  "
                  f"{ax:>13}  {sg:>13}  {co:>13}")


if __name__ == "__main__":
    main()
