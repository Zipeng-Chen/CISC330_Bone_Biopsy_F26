"""Needle-axis/sphere geometry for the assignment's forward biopsy plan.

All positions, radii, and returned distances use the caller's common length
unit (the project convention is millimetres). Directions are normalized, so
intersection parameters are physical distances rather than arbitrary factors.
The targeting error uses the infinite needle axis, as required by the rubric;
harvestable core and insertion depth use only the forward ray from the tip.
"""

from dataclasses import dataclass
from math import sqrt

import numpy as np

from ._validation import nonnegative_scalar, unit_vector, vector3


@dataclass(frozen=True)
class BiopsyPlan:
    """Geometric plan, without any mechanical model of the coring instrument.

    ``insertion_depth`` is the distance from the current tip to the far end of
    the forward tumor intersection. ``entry_depth`` is clipped at zero when
    the tip is already inside the sphere. A tangent returns its contact depth
    and a zero-length core; it does not represent a successful tissue sample.

    Status ``tangent`` describes a tangent supporting line. A tip starting on
    the sphere and pointing outward has status ``intersects`` but zero core
    and zero insertion depth. ``miss`` and ``behind`` have no forward depths.
    """

    status: str
    targeting_error: float
    tissue_core_length: float
    insertion_depth: float | None
    entry_depth: float | None
    exit_depth: float | None


def _length_tolerance(*lengths: float) -> float:
    """Roundoff allowance in the same units, with no millimetre-sized floor."""
    scale = max(*(abs(value) for value in lengths), np.finfo(float).tiny)
    return float(64.0 * np.finfo(float).eps * scale)


def distance_point_to_line(point, line_point, direction) -> float:
    """Return the shortest distance to the infinite line through ``line_point``.

    For unit direction u and displacement w, the distance is ||w x u||.
    This metric does not indicate whether the point lies ahead of the tip.
    """
    # Normalize the axis so distance is independent of direction-vector length.
    target = vector3(point, "point")
    origin = vector3(line_point, "line_point")
    axis = unit_vector(direction, "direction")
    # The cross-product magnitude is the perpendicular distance to the axis.
    return float(np.linalg.norm(np.cross(target - origin, axis)))


def intersect_line_sphere(line_point, direction, center, radius):
    """Return ordered signed entry/exit distances, or ``None`` for a line miss.

    The line is x(s) = line_point + s*u, where u is the normalized direction.
    Let a = (center-line_point).u and d be the point-to-line distance. Then
    s_near/far = a -/+ sqrt(radius**2 - d**2). Negative distances are retained
    here; ``plan_biopsy`` clips the intersection to the forward ray.

    A tangent (including a zero-radius sphere on the line) returns (a, a).
    Only a small relative floating-point tolerance is used at tangency.
    """
    # Validate the sphere and use a unit axis so signed parameters are in mm.
    origin = vector3(line_point, "line_point")
    axis = unit_vector(direction, "direction")
    target = vector3(center, "center")
    sphere_radius = nonnegative_scalar(radius, "radius")
    # Separate the center displacement into axial and perpendicular distances.
    displacement = target - origin
    along = float(np.dot(displacement, axis))
    distance = float(np.linalg.norm(np.cross(displacement, axis)))
    radial_gap = sphere_radius - distance
    # Include displacement scale to cover roundoff on an oblique axis.
    tolerance = _length_tolerance(sphere_radius, distance, np.linalg.norm(displacement))

    # A negative radial gap misses the sphere; a zero gap is tangent.
    if radial_gap < -tolerance:
        return None
    if abs(radial_gap) <= tolerance:
        return along, along

    # Factoring the squared difference avoids subtraction of two squares.
    half_chord = sqrt(radial_gap) * sqrt(sphere_radius + distance)
    return along - half_chord, along + half_chord


def plan_biopsy(tip, direction, tumor_center, tumor_radius) -> BiopsyPlan:
    """Compute axis error, forward core length, and full-core insertion depth.

    The spherical tumor and straight needle axis are the assignment's idealized
    geometry. The needle's finite diameter, cutter length, and mechanical
    effects are not modeled. A miss or an entirely rearward intersection
    returns zero core and ``None`` for all insertion/contact depths.
    """
    origin = vector3(tip, "tip")
    axis = unit_vector(direction, "direction")
    center = vector3(tumor_center, "tumor_center")
    radius = nonnegative_scalar(tumor_radius, "tumor_radius")
    # Compute the infinite-axis error and its signed sphere intersections.
    error = distance_point_to_line(center, origin, axis)
    intersections = intersect_line_sphere(origin, axis, center, radius)

    if intersections is None:
        return BiopsyPlan("miss", error, 0.0, None, None, None)

    near, far = intersections
    tangent = near == far
    tolerance = _length_tolerance(near, far, radius)
    # Intersections behind the tip cannot be reached by forward insertion.
    if far < -tolerance:
        return BiopsyPlan("behind", error, 0.0, None, None, None)

    # Absorb only floating-point roundoff at the starting surface.
    if abs(far) <= tolerance:
        far = 0.0
    if abs(near) <= tolerance:
        near = 0.0
    # Clip entry at the tip and count only the remaining forward tissue.
    entry = max(0.0, near)
    exit_ = max(0.0, far)
    core = max(0.0, exit_ - entry)
    status = "tangent" if tangent else "intersects"
    return BiopsyPlan(status, error, core, exit_, entry, exit_)
