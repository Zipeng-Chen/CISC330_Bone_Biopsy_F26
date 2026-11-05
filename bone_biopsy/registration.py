"""Least-squares rigid registration with known marker correspondences.

Vectors are mathematically columns: x_destination = R @ x_source + t.
For convenient storage, marker arrays have one point per row, shape (N, 3).
"""

from dataclasses import dataclass

import numpy as np

from ._validation import points3, unit_vector, vector3


@dataclass(frozen=True)
class RigidTransform:
    """Source-to-destination pose and RMS marker-fit residual in mm."""

    rotation: np.ndarray
    translation: np.ndarray
    rms_residual: float

    def transform_point(self, point):
        """Apply rotation and translation to a position."""
        return self.rotation @ vector3(point, "point") + self.translation

    def transform_points(self, points):
        """Transform a row-wise array of positions."""
        return points3(points, "points") @ self.rotation.T + self.translation

    def transform_direction(self, direction):
        """Rotate and normalize a direction; translation never acts on it."""
        return unit_vector(self.rotation @ unit_vector(direction))


@dataclass(frozen=True)
class RegisteredTool:
    """Tool tip and unit insertion axis expressed in the CT frame."""

    tip: np.ndarray
    axis: np.ndarray
    transform: RigidTransform


def _check_noncollinear(centered_points, name):
    singular_values = np.linalg.svd(centered_points, compute_uv=False)
    # Noncollinear triples have rank two; rank three is not required.
    if (singular_values[0] == 0
            or singular_values[1] <= 1e-10 * singular_values[0]):
        raise ValueError(f"{name} must not be coincident or nearly collinear")


def register_rigid(source_markers, destination_markers):
    """Fit a proper rigid transform between corresponding 3D marker sets.

    The i-th source marker must correspond to the i-th destination marker.
    At least three noncollinear markers are required. No scaling or marker
    matching is performed. With noisy measurements, minimize summed squared
    marker residuals. rms_residual measures marker fit, not target accuracy.
    """
    # Validate corresponding marker rows before fitting the pose.
    source = points3(source_markers, "source_markers", minimum=3)
    destination = points3(destination_markers, "destination_markers", minimum=3)
    if source.shape != destination.shape:
        raise ValueError("marker sets must have the same shape and correspondence")

    # Remove translation so the centered points determine rotation.
    source_center = source.mean(axis=0)
    destination_center = destination.mean(axis=0)
    source_zero = source - source_center
    destination_zero = destination - destination_center
    _check_noncollinear(source_zero, "source_markers")
    _check_noncollinear(destination_zero, "destination_markers")

    # Fit the least-squares rotation using H = U S V^T.
    u, _, vt = np.linalg.svd(source_zero.T @ destination_zero)
    # Enforce det(R) = +1 to prevent a reflected coordinate frame.
    correction = np.eye(3)
    correction[2, 2] = 1.0 if np.linalg.det(vt.T @ u.T) >= 0 else -1.0
    rotation = vt.T @ correction @ u.T
    # Recover translation by aligning the rotated source centroid.
    translation = destination_center - rotation @ source_center

    # Measure agreement at the markers; this is not target error.
    residuals = source @ rotation.T + translation - destination
    rms_residual = float(np.sqrt(np.mean(np.sum(residuals ** 2, axis=1))))
    return RigidTransform(rotation, translation, rms_residual)


def register_tool(model_markers, ct_markers, tip_tool, axis_tool):
    """Return calibrated tip and axis in CT coordinates, plus the fitted pose.

    tip_tool and axis_tool must use the same model frame as model_markers.
    Ground-truth values can be supplied for isolated registration tests.
    """
    transform = register_rigid(model_markers, ct_markers)
    # Positions receive rotation and translation; directions receive rotation only.
    return RegisteredTool(
        tip=transform.transform_point(tip_tool),
        axis=transform.transform_direction(axis_tool),
        transform=transform,
    )

