"""Geometric measurements on binary nodule masks (2D/3D Feret, slice ranges)."""

from __future__ import annotations

import re

import numpy as np

CONNECTIVITY = 3
MAX_PTS_2D = 5000
MAX_PTS_3D = 3000

# case_id pattern: <PatientID>_<StudyDate>_<SeriesUID>[a]
# The trailing 'a' marks an alternative reconstruction of the same study.
_CASE_RE = re.compile(r"^(?P<study_id>.+?_\d{8}_[\d.]+?)(?P<variant>a?)$")


def parse_case_id(case_id: str) -> tuple[str, str]:
    """Return (study_id, variant) where variant is 'base' or 'alt_a'."""
    m = _CASE_RE.match(case_id)
    if not m:
        return case_id, "base"
    return m.group("study_id"), ("alt_a" if m.group("variant") == "a" else "base")


def _subsample(coords_mm: np.ndarray, max_pts: int) -> np.ndarray:
    if len(coords_mm) > max_pts:
        idx = np.linspace(0, len(coords_mm) - 1, max_pts).astype(int)
        return coords_mm[idx]
    return coords_mm


def feret_max_2d_mm(mask_2d: np.ndarray, s_row_mm: float, s_col_mm: float) -> float:
    """Maximum pairwise distance (Feret) between mask voxels in physical mm."""
    rows, cols = np.where(mask_2d > 0)
    if len(rows) < 2:
        return 0.0
    coords_mm = np.column_stack([rows * s_row_mm, cols * s_col_mm])
    coords_mm = _subsample(coords_mm, MAX_PTS_2D)
    diffs = coords_mm[:, None, :] - coords_mm[None, :, :]
    return float(np.linalg.norm(diffs, axis=2).max())


def feret_2d_perpendicular(coords_mm: np.ndarray) -> tuple[float, float]:
    """Long-axis (max Feret) and short-axis (perpendicular extent) — Lung-RADS style."""
    if len(coords_mm) < 2:
        return 0.0, 0.0
    coords_mm = _subsample(coords_mm, MAX_PTS_2D)

    diffs = coords_mm[:, None, :] - coords_mm[None, :, :]
    dists = np.linalg.norm(diffs, axis=2)
    i, j = np.unravel_index(np.argmax(dists), dists.shape)
    long_axis = float(dists[i, j])
    if long_axis == 0:
        return 0.0, 0.0

    v = coords_mm[j] - coords_mm[i]
    v = v / np.linalg.norm(v)
    v_perp = np.array([-v[1], v[0]])
    proj = coords_mm @ v_perp
    short_axis = float(proj.max() - proj.min())
    return long_axis, short_axis


def feret_3d_mm(mask_3d: np.ndarray, spacing_zyx: tuple[float, float, float]) -> float:
    """Maximum pairwise 3D distance between mask voxels in mm."""
    sz, sy, sx = spacing_zyx
    coords = np.argwhere(mask_3d)
    if len(coords) < 2:
        return 0.0
    coords_mm = coords * np.array([sz, sy, sx])
    coords_mm = _subsample(coords_mm, MAX_PTS_3D)
    diffs = coords_mm[:, None, :] - coords_mm[None, :, :]
    return float(np.linalg.norm(diffs, axis=2).max())


def view_axes(view: str) -> tuple[int, int]:
    """Return (axis_to_drop, axis_layout) for a named anatomical view.

    Volume axes are (Z, Y, X). Returns the axis perpendicular to the view plane.
    """
    if view == "axial":
        return 0, 1  # drop Z
    if view == "coronal":
        return 1, 2  # drop Y
    if view == "sagittal":
        return 2, 3  # drop X
    raise ValueError(f"unknown view: {view}")


def measure_view(
    mask_3d: np.ndarray,
    view: str,
    spacing_zyx: tuple[float, float, float],
) -> dict:
    """Find the largest-area slice in `view` and measure long/short axes there.

    Returns dict with: long_mm, short_mm, n_slices, first, last (1-indexed).
    """
    sz, sy, sx = spacing_zyx
    if view == "axial":
        axis_drop, s_row, s_col = 0, sy, sx
    elif view == "coronal":
        axis_drop, s_row, s_col = 1, sz, sx
    elif view == "sagittal":
        axis_drop, s_row, s_col = 2, sz, sy
    else:
        raise ValueError(f"unknown view: {view}")

    sum_axes = tuple(a for a in (0, 1, 2) if a != axis_drop)
    areas = mask_3d.sum(axis=sum_axes)
    where_present = np.where(areas > 0)[0]
    if where_present.size == 0:
        return {"long_mm": 0.0, "short_mm": 0.0, "n_slices": 0, "first": 0, "last": 0}

    idx_max = int(np.argmax(areas))
    slicer = [slice(None)] * 3
    slicer[axis_drop] = idx_max
    slc = mask_3d[tuple(slicer)].astype(np.uint8)
    rows, cols = np.where(slc > 0)
    coords_mm = np.column_stack([rows * s_row, cols * s_col])
    long_ax, short_ax = feret_2d_perpendicular(coords_mm)

    return {
        "long_mm": long_ax,
        "short_mm": short_ax,
        "n_slices": int(where_present.size),
        "first": int(where_present.min()) + 1,
        "last": int(where_present.max()) + 1,
    }
