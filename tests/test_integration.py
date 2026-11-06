"""Check the complete marker-registration-to-biopsy-plan calculation."""

import unittest

from numpy.testing import assert_allclose

from bone_biopsy import assignment_tool_model, plan_biopsy, register_tool


class RegistrationPlanningIntegrationTests(unittest.TestCase):
    def test_registered_tool_produces_expected_off_center_plan(self):
        model = assignment_tool_model()
        ct_markers = [[-170, 15, 2], [-170, 75, 2], [-230, 15, 2]]
        tool = register_tool(model.markers, ct_markers, model.tip, model.axis)
        plan = plan_biopsy(tool.tip, tool.axis, [110, -2, 2], 5)
        # Tip=(10,-5,2), axis=+x. Lateral offset=3, half chord=4.
        self.assertEqual(plan.status, "intersects")
        assert_allclose(
            [plan.targeting_error, plan.entry_depth,
             plan.tissue_core_length, plan.insertion_depth],
            [3, 96, 8, 104], atol=1e-12, rtol=0,
        )

    def test_calibrated_nonzero_tip_is_accepted_by_public_pipeline(self):
        model = assignment_tool_model()
        ct_markers = [[-170, 15, 2], [-170, 75, 2], [-230, 15, 2]]
        # A nonzero calibrated tip checks both its rotation and translation.
        tool = register_tool(model.markers, ct_markers, [1, 2, 3], [0, -2, 0])
        assert_allclose(tool.tip, [8, -4, 5], atol=1e-12, rtol=0)
        plan = plan_biopsy(tool.tip, tool.axis, [108, -4, 5], 5)
        assert_allclose(
            [plan.targeting_error, plan.tissue_core_length, plan.insertion_depth],
            [0, 10, 105], atol=1e-12, rtol=0,
        )

    def test_figure_6_displaced_target_uses_actual_registered_axis(self):
        model = assignment_tool_model()
        # Shift the tool 3 mm sideways; the planned target was (0,-100,0).
        tool = register_tool(
            model.markers, [[23, 180, 0], [83, 180, 0], [23, 240, 0]],
            model.tip, model.axis,
        )
        plan = plan_biopsy(tool.tip, tool.axis, [0, -104, 0], 5)
        self.assertEqual(plan.status, "intersects")
        # At axial distance 104, offset 3 and radius 5 give half-chord 4.
        assert_allclose(
            [plan.targeting_error, plan.entry_depth,
             plan.tissue_core_length, plan.insertion_depth],
            [3, 100, 8, 108], atol=1e-12, rtol=0,
        )


if __name__ == "__main__":
    unittest.main()
