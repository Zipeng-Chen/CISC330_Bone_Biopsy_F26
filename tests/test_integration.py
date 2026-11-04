"""Check the complete marker-registration-to-biopsy-plan calculation."""

import unittest

from numpy.testing import assert_allclose

from bone_biopsy import assignment_tool_model, plan_biopsy, register_tool


class RegistrationPlanningIntegrationTests(unittest.TestCase):
    def test_registered_tool_produces_expected_off_center_plan(self):
        model = assignment_tool_model(mm_per_figure_unit=1)
        ct_markers = [[-8, -3, 2], [-8, 3, 2], [-14, -3, 2]]
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
        model = assignment_tool_model(mm_per_figure_unit=1)
        ct_markers = [[-8, -3, 2], [-8, 3, 2], [-14, -3, 2]]
        # A nonzero calibrated tip checks both its rotation and translation.
        tool = register_tool(model.markers, ct_markers, [1, 2, 3], [0, -2, 0])
        assert_allclose(tool.tip, [8, -4, 5], atol=1e-12, rtol=0)
        plan = plan_biopsy(tool.tip, tool.axis, [108, -4, 5], 5)
        assert_allclose(
            [plan.targeting_error, plan.tissue_core_length, plan.insertion_depth],
            [0, 10, 105], atol=1e-12, rtol=0,
        )


if __name__ == "__main__":
    unittest.main()
