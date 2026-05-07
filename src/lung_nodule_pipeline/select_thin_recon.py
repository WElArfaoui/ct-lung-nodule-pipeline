"""Filter nodules_report.csv to the thinnest reconstruction (smallest z-spacing) per study."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)


def select_thin(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = df.copy()
    df["z_spacing"] = df["spacing_xyz_mm"].str.split("/").str[2].astype(float)

    chosen = (
        df.groupby(["study_id", "recon_variant"])["z_spacing"].first()
        .reset_index()
        .sort_values(["study_id", "z_spacing"])
        .drop_duplicates("study_id", keep="first")
        .rename(columns={"recon_variant": "chosen_variant"})
    )

    df_thin = df.merge(chosen[["study_id", "chosen_variant"]], on="study_id")
    df_thin = (
        df_thin[df_thin["recon_variant"] == df_thin["chosen_variant"]]
        .drop(columns=["chosen_variant"])
    )
    return chosen, df_thin


def build_thin_study_summary(df_thin: pd.DataFrame) -> pd.DataFrame:
    return df_thin.groupby("study_id").agg(
        n_nodules=("nodule_id", "count"),
        n_ge_6mm=("diam_axial_mm", lambda s: int((s >= 6).sum())),
        max_diam_axial=("diam_axial_mm", "max"),
        max_diam_sagittal=("diam_sagittal_mm", "max"),
        total_volume_mm3=("volume_mm3", "sum"),
        spacing=("spacing_xyz_mm", "first"),
    ).reset_index()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reports-dir", required=True, type=Path,
                   help="Directory containing nodules_report.csv")
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    src = args.reports_dir / "nodules_report.csv"
    df = pd.read_csv(src)

    chosen, df_thin = select_thin(df)
    log.info("=== Chosen reconstruction per study ===\n%s", chosen.to_string(index=False))

    out_thin = args.reports_dir / "nodules_report_thin.csv"
    df_thin.to_csv(out_thin, index=False)
    log.info("wrote %s: %d nodules / %d studies",
             out_thin, len(df_thin), df_thin["study_id"].nunique())

    study = build_thin_study_summary(df_thin)
    out_summary = args.reports_dir / "study_summary_thin.csv"
    study.to_csv(out_summary, index=False)
    log.info("wrote %s\n%s", out_summary, study.to_string(index=False))


if __name__ == "__main__":
    main()
