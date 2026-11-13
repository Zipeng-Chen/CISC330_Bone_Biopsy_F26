"""Reusable rotation, plane-fitting, and direction-error calculations"""

from dataclasses import dataclass

import numpy as np

from ._validation import points3, unit_vector


@dataclass(frozen=True)
class PlaneFit:
    """Least-squares plane with an unsigned unit normal and RMS distance"""

    center: np.ndarray
    normal: np.ndarray
    rms_distance: float
    singular_values: np.ndarray


def rotation_about_axis(axis, angle_radians):
    """Return a right-handed rotation using Rodrigues' formula"""
    direction = unit_vector(axis)
    angle = np.asarray(angle_radians, dtype=float)
    if angle.ndim != 0 or not np.isfinite(angle):
        raise ValueError("angle_radians must be a finite scalar")
    x, y, z = direction
    skew = np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])
    # Separate the unchanged axial component from the rotating plane
    cosine, sine = np.cos(angle), np.sin(angle)
    return cosine * np.eye(3) + (1 - cosine) * np.outer(direction, direction) + sine * skew


def fit_plane(points):
    """Fit a plane to at least three noncollinear 3D points."""
    cloud = points3(points, "points", minimum=3)
    center = cloud.mean(axis=0)
    centered = cloud - center
    # The least-varying direction is perpendicular to the best-fit plane
    _, singular_values, vt = np.linalg.svd(centered, full_matrices=False)
    if singular_values[0] == 0 or singular_values[1] <= 1e-10 * singular_values[0]:
        raise ValueError("plane points must not be coincident or nearly collinear")
    normal = vt[-1]
    residual = centered @ normal
    return PlaneFit(center, normal, float(np.sqrt(np.mean(residual**2))), singular_values)


def direction_error_degrees(estimated, expected):
    """Return the signed-direction disagreement in [0, 180] degrees"""
    first, second = unit_vector(estimated), unit_vector(expected)
    # atan2 retains accuracy near zero and still detects a reversed axis
    return float(np.degrees(np.arctan2(np.linalg.norm(np.cross(first, second)),
                                      np.dot(first, second))))
