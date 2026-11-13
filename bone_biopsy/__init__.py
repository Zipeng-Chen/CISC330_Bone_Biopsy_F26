"""Reusable calibration, registration, and planning utilities for CISC/CMPE 330 A2."""

from .calibration import (
    AxisCalibration, CalibratedTool, TipCalibration, calibrate_axis, calibrate_tip, calibrate_tool,
)
from .calibration_simulation import CalibrationMotion, simulate_pivot_motion, simulate_spin_motion
from .noise import add_marker_noise, sigma_for_radial_bound

from .models import ToolModel, assignment_tool_model
from .planning import (
    BiopsyPlan, distance_point_to_line, intersect_line_sphere, plan_biopsy,
)
from .registration import RegisteredTool, RigidTransform, register_rigid, register_tool

__all__ = [
    "AxisCalibration", "CalibratedTool", "TipCalibration", "calibrate_axis",
    "calibrate_tip", "calibrate_tool", "CalibrationMotion", "simulate_pivot_motion",
    "simulate_spin_motion", "add_marker_noise", "sigma_for_radial_bound",
    "ToolModel", "assignment_tool_model", "RegisteredTool", "RigidTransform",
    "register_rigid", "register_tool",
    "BiopsyPlan", "distance_point_to_line", "intersect_line_sphere", "plan_biopsy",
]

