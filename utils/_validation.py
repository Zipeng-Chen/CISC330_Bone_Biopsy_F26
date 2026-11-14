"""Shared checks for finite coordinates, directions, and scalar parameters."""

import numpy as np


def vector3(value, name):
    """Require one finite three-dimensional vector."""
    array = np.asarray(value, dtype=float)
    if array.shape != (3,) or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a finite vector with shape (3,)")
    return array


def points3(value, name, *, minimum=1):
    """Require a finite array with one 3D point per row."""
    array = np.asarray(value, dtype=float)
    if (array.ndim != 2 or array.shape[1] != 3
            or len(array) < minimum or not np.all(np.isfinite(array))):
        raise ValueError(f"{name} must contain at least {minimum} finite 3D points")
    return array


def unit_vector(value, name="direction"):
    """Normalize a finite, nonzero direction without changing its sign."""
    array = vector3(value, name)
    # Scale first to avoid overflow or underflow during normalization.
    scale = np.max(np.abs(array))
    if scale == 0:
        raise ValueError(f"{name} must be nonzero")
    scaled = array / scale
    return scaled / np.linalg.norm(scaled)


def nonnegative_scalar(value, name):
    """Require a finite scalar greater than or equal to zero."""
    array = np.asarray(value, dtype=float)
    if array.ndim != 0 or not np.isfinite(array) or array < 0:
        raise ValueError(f"{name} must be a finite nonnegative scalar")
    return float(array)
