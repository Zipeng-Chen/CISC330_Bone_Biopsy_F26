"""Ground-truth tool from Figure 2, converted from centimetres to millimetres."""

from dataclasses import dataclass

import numpy as np

from ._validation import nonnegative_scalar


@dataclass(frozen=True)
class ToolModel:
    """Marker positions and tip in mm, with a dimensionless insertion axis."""

    markers: np.ndarray
    tip: np.ndarray
    axis: np.ndarray


def assignment_tool_model(*, mm_per_figure_unit=10.0):
    """Return the tool in mm using the confirmed Figure 2 coordinates in cm.

    The default conversion is 10 mm per cm. Keep this value for the assignment;
    the keyword remains available for custom synthetic-scale experiments.
    Marker order is M1, M2, M3. The insertion axis is dimensionless.
    """
    scale = nonnegative_scalar(mm_per_figure_unit, "mm_per_figure_unit")
    if scale == 0:
        raise ValueError("mm_per_figure_unit must be positive")
    # Figure 2 uses cm: convert positions once, without scaling the direction.
    return ToolModel(
        markers=scale * np.array([[2., 18., 0.], [8., 18., 0.], [2., 24., 0.]]),
        tip=np.zeros(3),
        axis=np.array([0., -1., 0.]),
    )

