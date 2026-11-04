"""Independent analytic cases for the assignment's geometry requirements."""

import math
import unittest

import numpy as np

from bone_biopsy.planning import (
    distance_point_to_line,
    intersect_line_sphere,
    plan_biopsy,
)


class GeometryUtilitiesTests(unittest.TestCase):
    def test_point_to_infinite_line(self):
        self.assertAlmostEqual(
            distance_point_to_line([3, 4, -20], [0, 0, 2], [0, 0, 8]), 5.0
        )

    def test_intersection_parameters_are_distances(self):
        self.assertEqual(
            intersect_line_sphere([0, 0, -3], [0, 0, 5], [0, 0, 0], 1),
            (2.0, 4.0),
        )

    def test_line_intersection_keeps_negative_parameters(self):
        self.assertEqual(
            intersect_line_sphere([0, 0, 3], [0, 0, 1], [0, 0, 0], 1),
            (-4.0, -2.0),
        )

    def test_tangent_is_a_double_root(self):
        self.assertEqual(
            intersect_line_sphere([1, 0, -3], [0, 0, 1], [0, 0, 0], 1),
            (3.0, 3.0),
        )

    def test_near_tangent_miss_is_not_rounded_to_a_hit(self):
        self.assertIsNone(
            intersect_line_sphere([1 + 1e-10, 0, -3], [0, 0, 1], [0, 0, 0], 1)
        )

    def test_near_tangent_intersection_keeps_nonzero_chord(self):
        roots = intersect_line_sphere(
            [1 - 1e-10, 0, -3], [0, 0, 1], [0, 0, 0], 1
        )
        self.assertIsNotNone(roots)
        self.assertGreater(roots[1] - roots[0], 0)

    def test_geometry_rejects_invalid_inputs(self):
        for direction in ([0, 0, 0], [0, np.nan, 1], [0, 1]):
            with self.subTest(direction=direction), self.assertRaises(ValueError):
                distance_point_to_line([0, 0, 0], [0, 0, 0], direction)
        with self.assertRaises(ValueError):
            intersect_line_sphere([0, 0, 0], [0, 0, 1], [0, 0, 0], -1)


class BiopsyPlanningTests(unittest.TestCase):
    def test_required_intersecting_trajectory(self):
        plan = plan_biopsy([0, 0, -3], [0, 0, 1], [0, 0, 0], 1)
        self.assertEqual(plan.status, "intersects")
        self.assertAlmostEqual(plan.targeting_error, 0)
        self.assertAlmostEqual(plan.tissue_core_length, 2)
        self.assertAlmostEqual(plan.entry_depth, 2)
        self.assertAlmostEqual(plan.exit_depth, 4)
        self.assertAlmostEqual(plan.insertion_depth, 4)

    def test_required_tangent_trajectory(self):
        plan = plan_biopsy([1, 0, -3], [0, 0, 1], [0, 0, 0], 1)
        self.assertEqual(plan.status, "tangent")
        self.assertAlmostEqual(plan.targeting_error, 1)
        self.assertEqual(plan.tissue_core_length, 0)
        self.assertAlmostEqual(plan.insertion_depth, 3)
        self.assertAlmostEqual(plan.entry_depth, 3)
        self.assertAlmostEqual(plan.exit_depth, 3)

    def test_required_missing_trajectory(self):
        plan = plan_biopsy([2, 0, -3], [0, 0, 1], [0, 0, 0], 1)
        self.assertEqual(plan.status, "miss")
        self.assertAlmostEqual(plan.targeting_error, 2)
        self.assertEqual(plan.tissue_core_length, 0)
        self.assertIsNone(plan.insertion_depth)
        self.assertIsNone(plan.entry_depth)
        self.assertIsNone(plan.exit_depth)

    def test_sphere_behind_tip_has_no_forward_core(self):
        plan = plan_biopsy([0, 0, 3], [0, 0, 1], [0, 0, 0], 1)
        self.assertEqual(plan.status, "behind")
        # Axis error alone cannot tell whether the needle points at the target.
        self.assertEqual(plan.targeting_error, 0)
        self.assertEqual(plan.tissue_core_length, 0)
        self.assertIsNone(plan.insertion_depth)

    def test_tip_inside_sphere_only_counts_remaining_forward_tissue(self):
        plan = plan_biopsy([0, 0, 0], [0, 0, 1], [0, 0, 0], 1)
        self.assertEqual(plan.status, "intersects")
        self.assertEqual(plan.entry_depth, 0)
        self.assertAlmostEqual(plan.tissue_core_length, 1)
        self.assertAlmostEqual(plan.insertion_depth, 1)

    def test_surface_tip_pointing_inward(self):
        plan = plan_biopsy([0, 0, -1], [0, 0, 1], [0, 0, 0], 1)
        self.assertEqual(plan.status, "intersects")
        self.assertEqual(plan.entry_depth, 0)
        self.assertAlmostEqual(plan.tissue_core_length, 2)
        self.assertAlmostEqual(plan.insertion_depth, 2)

    def test_surface_tip_pointing_outward_has_no_forward_tissue(self):
        plan = plan_biopsy([0, 0, 1], [0, 0, 1], [0, 0, 0], 1)
        self.assertEqual(plan.status, "intersects")
        self.assertEqual(plan.tissue_core_length, 0)
        self.assertEqual(plan.insertion_depth, 0)
        self.assertEqual(plan.entry_depth, 0)
        self.assertEqual(plan.exit_depth, 0)

    def test_nonunit_direction_does_not_change_depths(self):
        plan = plan_biopsy([0, 0, -3], [0, 0, 10], [0, 0, 0], 1)
        self.assertEqual(plan.insertion_depth, 4)
        self.assertEqual(plan.tissue_core_length, 2)

    def test_reversing_direction_changes_forward_reachability(self):
        plan = plan_biopsy([0, 0, -3], [0, 0, -1], [0, 0, 0], 1)
        self.assertEqual(plan.status, "behind")
        self.assertIsNone(plan.insertion_depth)

    def test_oblique_trajectory_and_translated_sphere(self):
        axis = np.array([1.0, 2.0, 2.0]) / 3
        perpendicular = np.array([2.0, -1.0, 0.0]) / math.sqrt(5)
        center = np.array([5.0, -2.0, 7.0])
        tip = center - 5 * axis + math.sqrt(3) * perpendicular
        # Radius 2, offset sqrt(3): half chord sqrt(4-3)=1.
        plan = plan_biopsy(tip, axis, center, 2)
        self.assertEqual(plan.status, "intersects")
        self.assertAlmostEqual(plan.targeting_error, math.sqrt(3))
        self.assertAlmostEqual(plan.tissue_core_length, 2)
        self.assertAlmostEqual(plan.entry_depth, 4)
        self.assertAlmostEqual(plan.insertion_depth, 6)

    def test_tiny_sphere_is_not_collapsed_by_absolute_tolerance(self):
        plan = plan_biopsy([0, 0, -3e-9], [0, 0, 1], [0, 0, 0], 1e-9)
        self.assertEqual(plan.status, "intersects")
        self.assertAlmostEqual(plan.tissue_core_length, 2e-9, delta=1e-23)
        self.assertAlmostEqual(plan.insertion_depth, 4e-9, delta=1e-23)

    def test_zero_radius_on_axis_is_zero_length_contact(self):
        plan = plan_biopsy([0, 0, -3], [0, 0, 1], [0, 0, 0], 0)
        self.assertEqual(plan.status, "tangent")
        self.assertEqual(plan.tissue_core_length, 0)
        self.assertEqual(plan.insertion_depth, 3)

    def test_zero_radius_off_axis_is_miss(self):
        plan = plan_biopsy([1, 0, -3], [0, 0, 1], [0, 0, 0], 0)
        self.assertEqual(plan.status, "miss")
        self.assertIsNone(plan.insertion_depth)

    def test_zero_radius_on_oblique_axis_tolerates_roundoff(self):
        # The target is on the ray at distance sqrt(2^2 + 3^2 + 7^2).
        plan = plan_biopsy([-2, -3, -7], [2, 3, 7], [0, 0, 0], 0)
        self.assertEqual(plan.status, "tangent")
        self.assertEqual(plan.tissue_core_length, 0)
        self.assertAlmostEqual(plan.insertion_depth, math.sqrt(62))

    def test_tangent_behind_tip(self):
        plan = plan_biopsy([1, 0, 3], [0, 0, 1], [0, 0, 0], 1)
        self.assertEqual(plan.status, "behind")
        self.assertEqual(plan.tissue_core_length, 0)
        self.assertIsNone(plan.insertion_depth)

    def test_inputs_are_not_modified(self):
        tip = np.array([0.0, 0.0, -3.0])
        direction = np.array([0.0, 0.0, 5.0])
        center = np.zeros(3)
        original = [array.copy() for array in (tip, direction, center)]
        plan_biopsy(tip, direction, center, 1)
        for actual, expected in zip((tip, direction, center), original):
            np.testing.assert_array_equal(actual, expected)

    def test_invalid_inputs_raise_useful_errors(self):
        invalid_cases = (
            ([0, 0], [0, 0, 1], [0, 0, 0], 1),
            ([0, 0, np.inf], [0, 0, 1], [0, 0, 0], 1),
            ([0, 0, 0], [0, 0, 0], [0, 0, 0], 1),
            ([0, 0, 0], [0, np.nan, 1], [0, 0, 0], 1),
            ([0, 0, 0], [0, 0, 1], [0, np.inf, 0], 1),
            ([0, 0, 0], [0, 0, 1], [0, 0, 0], -1),
            ([0, 0, 0], [0, 0, 1], [0, 0, 0], np.nan),
            ([0, 0, 0], [0, 0, 1], [0, 0, 0], np.inf),
        )
        for arguments in invalid_cases:
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                plan_biopsy(*arguments)


if __name__ == "__main__":
    unittest.main()
