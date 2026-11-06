"""Print and verify the required cases and the clarified Figure 6 error metric.

Run from the repository root: python -m examples.run_required_cases
All positions and distances are in mm; Figure 2 model coordinates start in cm.
"""

import numpy as np
from numpy.testing import assert_allclose

from bone_biopsy import (
    assignment_tool_model, distance_point_to_line, plan_biopsy, register_tool,
)


def vector_text(values):
    # Roundoff is suppressed for readability; assertions use unrounded results
    return "(" + ", ".join(f"{float(x):.6f}" for x in values) + ")"


def main():
    # Convert the confirmed Figure 2 centimetres to millimetres once.
    model = assignment_tool_model()
    print("ASSIGNMENT VALIDATION: positions and distances in mm")
    print("Figure 2 coordinates are cm; the model uses 10 mm per cm.")
    print("Predictions are independent analytic fixtures; tolerance = 1e-12.")
    print()
    print("REGISTRATION")
    # Use explicit marker coordinates from the analytic predictions
    cases = [
        ("translation", [[30, 175, 2], [90, 175, 2], [30, 235, 2]],
         [10, -5, 2], [0, -1, 0]),
        ("rotation", [[-180, 20, 0], [-180, 80, 0], [-240, 20, 0]],
         [0, 0, 0], [1, 0, 0]),
        ("rotation + translation", [[-170, 15, 2], [-170, 75, 2], [-230, 15, 2]],
         [10, -5, 2], [1, 0, 0]),
    ]
    last_tool = None
    for name, markers, expected_tip, expected_axis in cases:
        tool = register_tool(model.markers, markers, model.tip, model.axis)
        # Compare full-precision results before formatting them for display
        assert_allclose(tool.tip, expected_tip, atol=1e-12, rtol=0)
        assert_allclose(tool.axis, expected_axis, atol=1e-12, rtol=0)
        print(f"{name}:")
        print(f"  tip  expected {vector_text(expected_tip)}; computed {vector_text(tool.tip)}")
        print(f"  axis expected {vector_text(expected_axis)}; computed {vector_text(tool.axis)}")
        print(f"  tip difference norm: {np.linalg.norm(tool.tip - expected_tip):.3e} mm")
        print(f"  axis difference norm: {np.linalg.norm(tool.axis - expected_axis):.3e}")
        last_tool = tool

    print()
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


if __name__ == "__main__":
    main()
