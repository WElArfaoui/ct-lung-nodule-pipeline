"""Consolidated final CSV with approximate Lung-RADS classification."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

LUNG_RADS_ORDER = ["2", "3", "4A", "4B", "4X (mass)"]


def lung_rads_category(mean_ax_mm: float, vol_mm3: float) -> str:
    if vol_mm3 > 30000 or mean_ax_mm >= 30:
        return "4X (mass)"
    if mean_ax_mm >= 15:
        return "4B"
    if mean_ax_mm >= 8:
        return "4A"
    if mean_ax_mm >= 6:
        return "3"
    return "2"


def anisotropy_flag(long_ax: float, long_sg: float, long_co: float) -> str:
    if long_ax == 0:
        return ""
    if long_sg / long_ax >= 1.5 or long_co / long_ax >= 1.5:
        return "ANISOTROPIC_review"
    return ""


def highest_lr(series: pd.Series) -> str:
    valid = [c for c in series if c in LUNG_RADS_ORDER]
    if not valid:
        return ""
    return max(valid, key=LUNG_RADS_ORDER.index)


DESIRED_COLS = [
    "patient_id", "study_id", "nodule_id",
    "lung_rads", "qc_flag",
    "volume_mm3",
    "diam_long_axial_mm", "diam_short_axial_mm", "diam_mean_axial_mm",
    "diam_long_sagittal_mm", "diam_short_sagittal_mm",
    "diam_long_coronal_mm", "diam_short_coronal_mm",
    "diam_3d_max_mm", "diam_eq_sphere_mm",
    "slice_range_axial", "axial_n",
    "slice_range_sagittal", "sagittal_n",
    "slice_range_coronal", "coronal_n",
    "centroid_x_mm", "centroid_y_mm", "centroid_z_mm",
    "spacing_xyz_mm",
    "recon_variant", "case_id",
]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reports-dir", required=True, type=Path)
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    df = pd.read_csv(args.reports_dir / "nodules_full_diameters.csv")
    log.info("columns in full_diameters: %s", list(df.columns))

    # Backfill from thin report if needed
    needed_extra = [c for c in ["spacing_xyz_mm", "recon_variant"] if c not in df.columns]
    if needed_extra:
        log.info("missing %s; merging from nodules_report_thin.csv", needed_extra)
        thin = pd.read_csv(args.reports_dir / "nodules_report_thin.csv")
        keep = ["case_id"] + [c for c in needed_extra if c in thin.columns]
        df = df.merge(thin[keep].drop_duplicates("case_id"), on="case_id", how="left")

    df["lung_rads"] = df.apply(
        lambda r: lung_rads_category(r["diam_mean_axial_mm"], r["volume_mm3"]), axis=1)
    df["qc_flag"] = df.apply(
        lambda r: anisotropy_flag(
            r["diam_long_axial_mm"], r["diam_long_sagittal_mm"], r["diam_long_coronal_mm"]),
        axis=1)
    df["patient_id"] = df["study_id"].str.split("_").str[0]
    df["slice_range_axial"]    = df["axial_first"].astype(str)    + "-" + df["axial_last"].astype(str)
    df["slice_range_sagittal"] = df["sagittal_first"].astype(str) + "-" + df["sagittal_last"].astype(str)
    df["slice_range_coronal"]  = df["coronal_first"].astype(str)  + "-" + df["coronal_last"].astype(str)

    cols = [c for c in DESIRED_COLS if c in df.columns]
    missing = [c for c in DESIRED_COLS if c not in df.columns]
    if missing:
        log.warning("columns absent (skipped): %s", missing)

    final = df[cols].sort_values(["patient_id", "study_id", "nodule_id"])
    out = args.reports_dir / "nodules_FINAL_summary.csv"
    final.to_csv(out, index=False)
    log.info("wrote %s: %d rows / %d patients / %d studies",
             out, len(final), final["patient_id"].nunique(), final["study_id"].nunique())

    patient_summary = final.groupby(["patient_id", "study_id"]).agg(
        n_nodules=("nodule_id", "count"),
        n_actionable=("lung_rads",
                      lambda s: int(s.isin(["3", "4A", "4B", "4X (mass)"]).sum())),
        max_mean_diam=("diam_mean_axial_mm", "max"),
        max_3d=("diam_3d_max_mm", "max"),
        total_volume_mm3=("volume_mm3", "sum"),
        highest_lung_rads=("lung_rads", highest_lr),
        qc_flags=("qc_flag", lambda s: ";".join(x for x in s if x)),
    ).reset_index()

    out2 = args.reports_dir / "patient_summary_FINAL.csv"
    patient_summary.to_csv(out2, index=False)
    log.info("wrote %s: %d studies\n%s",
             out2, len(patient_summary), patient_summary.to_string(index=False))


if __name__ == "__main__":
    main()
