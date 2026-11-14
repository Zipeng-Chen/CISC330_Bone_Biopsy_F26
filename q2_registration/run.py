"""Check the three required registration cases against hand predictions."""

from pathlib import Path
import sys

# Allow both python -m and direct execution from another working directory.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from numpy.testing import assert_allclose

from bone_biopsy import assignment_tool_model, register_tool
from q2_registration.cases import REGISTRATION_CASES
from utils.reporting import run_tests, run_with_log

OUTPUTS = Path(__file__).resolve().parent / "outputs"


def vector_text(values):
    # Roundoff is suppressed for readability; assertions use unrounded results
    return "(" + ", ".join(f"{float(x):.6f}" for x in values) + ")"


def run():
    # Convert the confirmed Figure 2 centimetres to millimetres once.
    model = assignment_tool_model()
    print("ASSIGNMENT VALIDATION: positions and distances in mm")
    print("Figure 2 coordinates are cm; the model uses 10 mm per cm.")
    print("Predictions are independent analytic fixtures; tolerance = 1e-12.")
    print()
    print("REGISTRATION")
    # Use explicit marker coordinates from the analytic predictions
    cases = REGISTRATION_CASES
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

    print("\nPASS: translation, rotation, and combined registration agree.")
    run_tests("q2_registration.test_registration", OUTPUTS / "test_report.txt")


def main():
    run_with_log(run, OUTPUTS / "run_log.txt", "python -m q2_registration.run")


if __name__ == "__main__":
    main()
