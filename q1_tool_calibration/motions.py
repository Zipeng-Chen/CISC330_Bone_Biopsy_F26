"""Ground-truth constrained motions used only to generate calibration observations"""

from dataclasses import dataclass

import numpy as np

from utils._validation import nonnegative_scalar, points3, unit_vector, vector3
from utils.geometry import rotation_about_axis


@dataclass(frozen=True)
class CalibrationMotion:
    """Simulated observations and true poses for validation, never solver inputs"""

    observations: np.ndarray
    rotations: np.ndarray
    translations: np.ndarray


def _check_setup(model, pose_count, tracker_rotation):
    markers = points3(model.markers, "model.markers", minimum=3)
    tip, axis = vector3(model.tip, "model.tip"), unit_vector(model.axis)
    if (isinstance(pose_count, (bool, np.bool_))
            or not isinstance(pose_count, (int, np.integer)) or pose_count < 3):
        raise ValueError("pose_count must be an integer of at least 3")
    rotation = np.eye(3) if tracker_rotation is None else np.asarray(tracker_rotation, dtype=float)
    if (rotation.shape != (3, 3) or not np.all(np.isfinite(rotation))
            or not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-10, rtol=0)
            or not np.isclose(np.linalg.det(rotation), 1, atol=1e-10, rtol=0)):
        raise ValueError("tracker_rotation must be a proper orthonormal 3x3 rotation")
    return markers, tip, axis, rotation


def _observations(markers, tip, rotations, fixed_tip):
    # R_i * true_tip + t_i = fixed_tip enforces the physical tip constraint
    translations = fixed_tip - np.einsum("nij,j->ni", rotations, tip)
    measured = np.einsum("nij,mj->nmi", rotations, markers) + translations[:, None, :]
    return CalibrationMotion(measured, rotations, translations)


def simulate_pivot_motion(model, *, pose_count=60, cone_half_angle_degrees=40.0,
                          pivot_tracker=(25., -10., 40.), tracker_rotation=None, seed=17):
    """Generate multi-axis tilts inside a cone, with a fixed tip and varied roll."""
    markers, tip, axis, base = _check_setup(model, pose_count, tracker_rotation)
    half_angle = nonnegative_scalar(cone_half_angle_degrees, "cone_half_angle_degrees")
    if not 0 < half_angle <= 90:
        raise ValueError("cone_half_angle_degrees must lie in (0, 90]")
    pivot = vector3(pivot_tracker, "pivot_tracker")
    rng = np.random.default_rng(seed)
    # Build two perpendicular tilt directions for any model-frame shaft axis
    transverse = unit_vector(np.cross(axis, np.eye(3)[np.argmin(np.abs(axis))]))
    other = np.cross(axis, transverse)
    azimuth = rng.uniform(0, 2 * np.pi, pose_count)
    tilt = np.arccos(rng.uniform(np.cos(np.deg2rad(half_angle)), 1, pose_count))
    roll = rng.uniform(0, 2 * np.pi, pose_count)
    rotations = []
    for phi, theta, spin in zip(azimuth, tilt, roll):
        tilt_axis = np.cos(phi) * transverse + np.sin(phi) * other
        rotations.append(base @ rotation_about_axis(tilt_axis, theta)
                         @ rotation_about_axis(axis, spin))
    return _observations(markers, tip, np.array(rotations), pivot)


def simulate_spin_motion(model, *, pose_count=60, arc_degrees=360.0,
                         tip_tracker=(25., -10., 40.), tracker_rotation=None):
    """Spin around the shaft with no tilt, axial sliding, or tip movement."""
    markers, tip, axis, base = _check_setup(model, pose_count, tracker_rotation)
    arc = nonnegative_scalar(arc_degrees, "arc_degrees")
    if not 0 < arc <= 360:
        raise ValueError("arc_degrees must lie in (0, 360]")
    fixed_tip = vector3(tip_tracker, "tip_tracker")
    # Avoid repeating the first pose at the end of a full revolution
    angles = np.linspace(0, np.deg2rad(arc), pose_count, endpoint=arc < 360)
    rotations = np.array([base @ rotation_about_axis(axis, angle) for angle in angles])
    return _observations(markers, tip, rotations, fixed_tip)
