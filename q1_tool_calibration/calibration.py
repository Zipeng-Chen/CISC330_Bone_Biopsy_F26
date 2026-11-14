"""Tool-tip pivot calibration and tool-axis spin calibration from tracked markers.

Inputs are corresponding model markers (N, 3) and tracker observations
(poses, N, 3), in mm. The solvers never receive the ground-truth tip or axis.
"""

from dataclasses import dataclass

import numpy as np

from utils._validation import points3, unit_vector, vector3
from utils.geometry import direction_error_degrees, fit_plane
from q2_registration.registration import register_rigid


@dataclass(frozen=True)
class TipCalibration:
    """Estimated model-frame tip, tracker-frame pivot, and fit diagnostics"""

    tip: np.ndarray
    pivot_tracker: np.ndarray
    rms_residual: float
    condition_number: float


@dataclass(frozen=True)
class AxisCalibration:
    """Insertion direction in both frames and spin-constraint diagnostics"""

    axis: np.ndarray
    axis_tracker: np.ndarray
    plane_rms: np.ndarray
    used_marker_indices: tuple[int, ...]
    direction_scatter_degrees: float


@dataclass(frozen=True)
class CalibratedTool:
    """Model-frame tip and unit axis ready for register_tool"""

    tip: np.ndarray
    axis: np.ndarray
    tip_calibration: TipCalibration
    axis_calibration: AxisCalibration


def _recover_poses(model_markers, observations):
    markers = points3(model_markers, "model_markers", minimum=3)
    measured = np.asarray(observations, dtype=float)
    if (measured.ndim != 3 or measured.shape[1:] != markers.shape
            or len(measured) < 3 or not np.all(np.isfinite(measured))):
        raise ValueError("observations must have shape (poses>=3, markers, 3) and be finite")
    # Reuse registration to estimate model-to-tracker pose in every frame
    poses = [register_rigid(markers, frame) for frame in measured]
    rotations = np.array([pose.rotation for pose in poses])
    translations = np.array([pose.translation for pose in poses])
    return markers, measured, rotations, translations


def calibrate_tip(model_markers, pivot_observations):
    """Solve for a fixed tip pivoted over multiple, well-separated orientations.

    Keep the physical tip fixed in a dimple and tilt about more than one axis.
    Pure translation, repeated poses, or rotation about one axis is insufficient.
    The RMS residual is a constraint diagnostic, not a bound on tip error.
    """
    _, _, rotations, translations = _recover_poses(model_markers, pivot_observations)
    count = len(rotations)
    # Each pose satisfies R_i * tip - pivot = -t_i; both points are unknown
    system = np.concatenate((rotations, np.tile(-np.eye(3), (count, 1, 1))), axis=2)
    system = system.reshape(-1, 6)
    right_hand_side = -translations.reshape(-1)
    solution, _, rank, singular_values = np.linalg.lstsq(system, right_hand_side, rcond=1e-10)
    if rank != 6:
        raise ValueError("pivot motions do not determine all six unknowns; tilt about multiple axes")
    # Report the 3D scatter of the reconstructed tips around their common pivot
    residual = (system @ solution - right_hand_side).reshape(-1, 3)
    rms = float(np.sqrt(np.mean(np.sum(residual**2, axis=1))))
    condition = float(singular_values[0] / singular_values[-1])
    return TipCalibration(solution[:3], solution[3:], rms, condition)


def calibrate_axis(model_markers, spin_observations, *, calibrated_tip):
    """Fit marker-circle planes for a tool spun in a fixed guide sleeve.

    Keep the shaft line fixed, prevent axial sliding, and sample a broad spin
    arc (a full turn is preferred). Markers must lie behind the tip along the
    shaft on average; this physical convention resolves the normal's +/- sign
    using the calibrated tip. It does not use the ground-truth direction.
    """
    markers, measured, rotations, _ = _recover_poses(model_markers, spin_observations)
    tip = vector3(calibrated_tip, "calibrated_tip")
    fits, used = [], []
    for marker_index in range(len(markers)):
        # Each off-axis marker traces a circle in a plane normal to the shaft
        try:
            plane = fit_plane(measured[:, marker_index])
        except ValueError:
            continue  # An on-axis marker has no circle and supplies no normal
        if plane.singular_values[2] > 0.25 * plane.singular_values[1]:
            continue  # A volumetric cloud is not a well-resolved spin plane
        fits.append(plane)
        used.append(marker_index)
    if not fits:
        raise ValueError("spin observations do not contain a resolved marker-circle plane")

    # Align arbitrary SVD signs and weight wider circles more strongly
    weights = np.array([np.sum(plane.singular_values[:2]**2) for plane in fits])
    normals = np.array([plane.normal for plane in fits])
    reference = normals[np.argmax(weights)]
    normals[np.dot(normals, reference) < 0] *= -1
    tracker_axis = unit_vector(np.average(normals, axis=0, weights=weights))
    if np.min(normals @ tracker_axis) < np.cos(np.deg2rad(10)):
        raise ValueError("marker-circle planes are inconsistent with one fixed shaft")

    # Directions rotate back into the model frame; translations do not act on them
    model_directions = np.einsum("nji,j->ni", rotations, tracker_axis)
    axis = unit_vector(model_directions.mean(axis=0))
    scatter = np.array([direction_error_degrees(value, axis) for value in model_directions])
    if np.max(scatter) > 10:
        raise ValueError("spin poses do not preserve a fixed shaft direction")

    # Choose the insertion-facing sign from the calibrated tip, not a known axis
    toward_tip = tip - markers.mean(axis=0)
    alignment = float(np.dot(axis, toward_tip))
    if abs(alignment) <= 1e-10 * max(np.linalg.norm(toward_tip), np.finfo(float).tiny):
        raise ValueError("tip and marker layout do not resolve the insertion-direction sign")
    if alignment < 0:
        axis = -axis
        tracker_axis = -tracker_axis
    return AxisCalibration(axis, tracker_axis, np.array([plane.rms_distance for plane in fits]),
                           tuple(used), float(np.sqrt(np.mean(scatter**2))))


def calibrate_tool(model_markers, pivot_observations, spin_observations):
    """Calibrate the tip first, then use it to orient the spin-calibrated axis."""
    tip_result = calibrate_tip(model_markers, pivot_observations)
    axis_result = calibrate_axis(model_markers, spin_observations, calibrated_tip=tip_result.tip)
    return CalibratedTool(tip_result.tip, axis_result.axis, tip_result, axis_result)
