"""Check planning cases, the combined calculation, and the Figure 6 metric."""

from pathlib import Path
import sys

# Allow both python -m and direct execution from another working directory.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from numpy.testing import assert_allclose

from bone_biopsy import assignment_tool_model, distance_point_to_line, plan_biopsy, register_tool
from q2_registration.cases import REGISTRATION_CASES
from utils.reporting import run_tests, run_with_log

OUTPUTS = Path(__file__).resolve().parent / "outputs"


def run():
    print("All positions and distances are in mm; comparisons use tolerance 1e-12.")
    print("PLANNING: sphere center (0,0,0), radius 1, direction (0,0,1)")
    # Use a unit sphere to test intersection, tangency, and a complete miss
    plan_cases = [
        ("intersects", [0, 0, -3], (0., 2., 4.)),
        ("tangent", [1, 0, -3], (1., 0., 3.)),
        ("miss", [2, 0, -3], (2., 0., None)),
    ]
    for status, tip, expected in plan_cases:
        plan = plan_biopsy(tip, [0, 0, 1], [0, 0, 0], 1)
        actual = (plan.targeting_error, plan.tissue_core_length, plan.insertion_depth)
        if plan.status != status:
            raise AssertionError(f"Expected {status}, got {plan.status}")
        for computed, predicted in zip(actual, expected):
            if predicted is None:
                if computed is not None:
                    raise AssertionError("A miss must have no insertion depth")
            else:
                assert_allclose(computed, predicted, atol=1e-12, rtol=0)
        print(f"{status}: (error, core, insertion/contact depth)")
        print(f"  expected {expected}; computed {actual}")

    print()
    print("COMBINED REGISTRATION -> PLANNING")
    # A 3 mm offset in a 5 mm sphere gives an 8 mm chord at 100 mm axial depth
    model = assignment_tool_model()
    markers = REGISTRATION_CASES[-1][1]
    last_tool = register_tool(model.markers, markers, model.tip, model.axis)
    plan = plan_biopsy(last_tool.tip, last_tool.axis, [110, -2, 2], 5)
    actual = [plan.targeting_error, plan.entry_depth,
              plan.tissue_core_length, plan.insertion_depth]
    expected = [3., 96., 8., 104.]
    assert_allclose(actual, expected, atol=1e-12, rtol=0)
    print("  (error, entry depth, core length, insertion depth)")
    print(f"  expected {expected}; computed {actual}")
    print()
    print("FIGURE 6 METRIC: displaced target to actual needle axis")
    # Planned target is (0,-100,0); actual tool is shifted 3 mm along x.
    actual_tool = register_tool(
        model.markers, [[23, 180, 0], [83, 180, 0], [23, 240, 0]],
        model.tip, model.axis,
    )
    displaced_target = np.array([0., -104., 0.])
    error = distance_point_to_line(
        displaced_target, actual_tool.tip, actual_tool.axis,
    )
    # The 100 mm endpoint is 5 mm away, but Figure 6 measures the 3 mm normal.
    endpoint = actual_tool.tip + 100 * actual_tool.axis
    endpoint_distance = np.linalg.norm(displaced_target - endpoint)
    assert_allclose([error, endpoint_distance], [3., 5.], atol=1e-12, rtol=0)
    print(f"  perpendicular targeting error: {error:.6f} mm (expected 3)")
    print(f"  endpoint distance, for comparison only: {endpoint_distance:.6f} mm")
    print()
    print("PASS: required examples, combined case, and Figure 6 metric agree.")
    run_tests("q3_biopsy_plan.test_planning", OUTPUTS / "test_report.txt")


def main():
    run_with_log(run, OUTPUTS / "run_log.txt", "python -m q3_biopsy_plan.run")


if __name__ == "__main__":
    main()
