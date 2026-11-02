"""Ground-truth tool from assignment Figure 2, with an explicit unit conversion."""

from dataclasses import dataclass

import numpy as np

from ._validation import nonnegative_scalar


@dataclass(frozen=True)
class ToolModel:
    """Marker positions and tip in mm, with a dimensionless insertion axis."""

    markers: np.ndarray
    tip: np.ndarray
    axis: np.ndarray


def assignment_tool_model(*, mm_per_figure_unit):
    """Return the model in mm; the PDF does not label its coordinate units.

    Set mm_per_figure_unit to 1 for mm or 10 for cm. The physical scale must
    be supplied explicitly because it is not stated in Figure 2.
    """
    scale = nonnegative_scalar(mm_per_figure_unit, "mm_per_figure_unit")
    if scale == 0:
        raise ValueError("mm_per_figure_unit must be positive")
    # Convert positions only; the insertion axis is dimensionless.
    return ToolModel(
        markers=scale * np.array([[2., 18., 0.], [8., 18., 0.], [2., 24., 0.]]),
        tip=np.zeros(3),
        axis=np.array([0., -1., 0.]),
    )

