"""Reusable calibration, registration, and planning utilities for CISC/CMPE 330 A2."""

from q1_tool_calibration.calibration import (
    AxisCalibration, CalibratedTool, TipCalibration, calibrate_axis, calibrate_tip, calibrate_tool,
)
from q1_tool_calibration.motions import CalibrationMotion, simulate_pivot_motion, simulate_spin_motion
from utils.noise import add_marker_noise, sigma_for_radial_bound

from utils.models import ToolModel, assignment_tool_model
from q3_biopsy_plan.planning import (
    BiopsyPlan, distance_point_to_line, intersect_line_sphere, plan_biopsy,
)
from q2_registration.registration import RegisteredTool, RigidTransform, register_rigid, register_tool

__all__ = [
    "AxisCalibration", "CalibratedTool", "TipCalibration", "calibrate_axis",
    "calibrate_tip", "calibrate_tool", "CalibrationMotion", "simulate_pivot_motion",
    "simulate_spin_motion", "add_marker_noise", "sigma_for_radial_bound",
    "ToolModel", "assignment_tool_model", "RegisteredTool", "RigidTransform",
    "register_rigid", "register_tool",
    "BiopsyPlan", "distance_point_to_line", "intersect_line_sphere", "plan_biopsy",
]


# Older submodule imports still refer to the same question/shared modules.
from importlib import import_module as _import_module
import sys as _sys

_MODULE_PATHS = {
    "calibration": "q1_tool_calibration.calibration",
    "calibration_simulation": "q1_tool_calibration.motions",
    "registration": "q2_registration.registration",
    "planning": "q3_biopsy_plan.planning",
    "models": "utils.models",
    "geometry": "utils.geometry",
    "noise": "utils.noise",
    "_validation": "utils._validation",
}
for _name, _path in _MODULE_PATHS.items():
    _module = _import_module(_path)
    _sys.modules[f"{__name__}.{_name}"] = _module
    globals()[_name] = _module
del _name, _path, _module, _MODULE_PATHS, _import_module, _sys
