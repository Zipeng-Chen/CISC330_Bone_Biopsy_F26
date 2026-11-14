"""Calibration checks with independent geometry, exact poses, and noisy observations."""

import unittest

import numpy as np
from numpy.testing import assert_allclose

from bone_biopsy import (
    ToolModel, add_marker_noise, assignment_tool_model, calibrate_axis, calibrate_tip,
    calibrate_tool, plan_biopsy, register_tool, sigma_for_radial_bound,
    simulate_pivot_motion, simulate_spin_motion,
)
from bone_biopsy.geometry import direction_error_degrees, fit_plane, rotation_about_axis


class GeometryTests(unittest.TestCase):
    def test_right_handed_rotation(self):
        rotation = rotation_about_axis([0, 0, 7], np.pi / 2)
        assert_allclose(rotation @ [1, 0, 0], [0, 1, 0], atol=1e-15)
        assert_allclose(rotation.T @ rotation, np.eye(3), atol=1e-15)
        self.assertAlmostEqual(np.linalg.det(rotation), 1)

    def test_plane_recovers_independently_defined_normal(self):
        # Every point satisfies x + 2y + 3z = 6.
        plane = fit_plane([[6, 0, 0], [0, 3, 0], [0, 0, 2], [1, 1, 1]])
        self.assertAlmostEqual(abs(np.dot(plane.normal, np.array([1, 2, 3]) / np.sqrt(14))), 1)
        self.assertLess(plane.rms_distance, 1e-14)

    def test_degenerate_plane_is_rejected(self):
        for cloud in (np.zeros((4, 3)), [[0, 0, 0], [1, 0, 0], [2, 0, 0]]):
            with self.subTest(cloud=cloud), self.assertRaises(ValueError):
                fit_plane(cloud)

    def test_direction_error_preserves_axis_sign(self):
        self.assertEqual(direction_error_degrees([0, -1, 0], [0, 1, 0]), 180)
        self.assertEqual(direction_error_degrees([0, -1, 0], [0, -2, 0]), 0)
        self.assertEqual(direction_error_degrees([0, -1, 0], [1, 0, 0]), 90)

    def test_invalid_rotation_inputs(self):
        for axis, angle in (([0, 0, 0], 1), ([0, 1, 0], np.nan), ([0, 1, 0], [1])):
            with self.subTest(axis=axis, angle=angle), self.assertRaises(ValueError):
                rotation_about_axis(axis, angle)


class MarkerNoiseTests(unittest.TestCase):
    def test_sigma_uses_three_dimensional_radial_quantile(self):
        # Independent chi-square(3) 95th percentile, rather than a 1D 1.96 rule.
        expected = 1 / np.sqrt(7.814727903251179)
        self.assertAlmostEqual(sigma_for_radial_bound(), expected, places=12)
        self.assertAlmostEqual(sigma_for_radial_bound(2), 2 * expected, places=12)

    def test_empirical_radial_probability_and_coordinate_independence(self):
        sigma = sigma_for_radial_bound()
        errors = add_marker_noise(np.zeros((100000, 3)), sigma, rng=np.random.default_rng(41))
        self.assertAlmostEqual(np.mean(np.linalg.norm(errors, axis=1) < 1), .95, delta=.003)
        assert_allclose(errors.mean(axis=0), [0, 0, 0], atol=.004)
        assert_allclose(np.cov(errors.T), np.eye(3) * sigma**2, atol=.002)

    def test_noise_is_reproducible_and_does_not_modify_inputs(self):
        points = np.arange(27, dtype=float).reshape(3, 3, 3)
        original = points.copy()
        first = add_marker_noise(points, .3, rng=np.random.default_rng(8))
        second = add_marker_noise(points, .3, rng=np.random.default_rng(8))
        assert_allclose(first, second, atol=0, rtol=0)
        assert_allclose(points, original, atol=0, rtol=0)
        exact = add_marker_noise(points, 0, rng=np.random.default_rng(8))
        assert_allclose(exact, original, atol=0, rtol=0)

    def test_invalid_noise_parameters(self):
        for bound, probability in ((0, .95), (-1, .95), (1, 0), (1, 1), (1, np.nan)):
            with self.subTest(bound=bound, probability=probability), self.assertRaises(ValueError):
                sigma_for_radial_bound(bound, probability)
        with self.assertRaises(ValueError):
            add_marker_noise([[0, 0, 0]], -1, rng=np.random.default_rng(1))
        with self.assertRaises(ValueError):
            add_marker_noise([[0, np.nan, 0]], 1, rng=np.random.default_rng(1))


class ToolCalibrationTests(unittest.TestCase):
    def setUp(self):
        self.model = assignment_tool_model()
        self.base = rotation_about_axis([1, 2, 3], .4)
        self.pivot = simulate_pivot_motion(self.model, tracker_rotation=self.base)
        self.spin = simulate_spin_motion(self.model, tracker_rotation=self.base)

    def test_simulated_pivot_keeps_tip_fixed_and_excites_multiple_axes(self):
        tips = np.einsum("nij,j->ni", self.pivot.rotations, self.model.tip) + self.pivot.translations
        assert_allclose(tips, np.tile([25, -10, 40], (60, 1)), atol=1e-12)
        system = np.concatenate((self.pivot.rotations, np.tile(-np.eye(3), (60, 1, 1))), axis=2)
        self.assertEqual(np.linalg.matrix_rank(system.reshape(-1, 6)), 6)
        # Rigid motion must preserve the 60 mm separation of M1 and M2.
        assert_allclose(np.linalg.norm(self.pivot.observations[:, 1] - self.pivot.observations[:, 0], axis=1),
                        60, atol=1e-12)

    def test_simulated_spin_keeps_axis_fixed_and_markers_on_circles(self):
        expected_axis = self.base @ [0, -1, 0]
        axes = np.einsum("nij,j->ni", self.spin.rotations, self.model.axis)
        assert_allclose(axes, np.tile(expected_axis, (60, 1)), atol=1e-12)
        # The confirmed model gives radii 20, 80, and 20 mm.
        for index, expected_radius in enumerate([20, 80, 20]):
            cloud = self.spin.observations[:, index]
            centered = cloud - cloud.mean(axis=0)
            assert_allclose(np.linalg.norm(centered, axis=1), expected_radius, atol=1e-12)
            assert_allclose(centered @ expected_axis, 0, atol=1e-12)

    def test_exact_assignment_model_is_recovered(self):
        result = calibrate_tool(self.model.markers, self.pivot.observations, self.spin.observations)
        assert_allclose(result.tip, [0, 0, 0], atol=1e-10, rtol=0)
        assert_allclose(result.axis, [0, -1, 0], atol=1e-12, rtol=0)
        assert_allclose(result.tip_calibration.pivot_tracker, [25, -10, 40], atol=1e-10)
        self.assertLess(result.tip_calibration.rms_residual, 1e-10)
        self.assertLess(max(result.axis_calibration.plane_rms), 1e-10)

    def test_nonzero_tip_and_oblique_model_axis_are_recovered(self):
        # Re-express the model in another frame to catch hard-coded tip/axis values.
        rotation = rotation_about_axis([2, -1, 3], .9)
        offset = np.array([11., -9., 6.])
        model = ToolModel(self.model.markers @ rotation.T + offset, offset,
                          rotation @ self.model.axis)
        pivot = simulate_pivot_motion(model, pivot_tracker=[40, -30, 50], tracker_rotation=self.base)
        spin = simulate_spin_motion(model, tip_tracker=[-10, 20, 30], tracker_rotation=self.base)
        result = calibrate_tool(model.markers, pivot.observations, spin.observations)
        assert_allclose(result.tip, offset, atol=1e-10, rtol=0)
        assert_allclose(result.axis, model.axis, atol=1e-12, rtol=0)

    def test_independent_pivot_fixture_recovers_nonzero_tip(self):
        rotations = np.array([
            np.eye(3), [[1, 0, 0], [0, 0, -1], [0, 1, 0]],
            [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
            [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
        ])
        tip, pivot = np.array([1, 2, 3]), np.array([10, 20, 30])
        measured = np.array([(self.model.markers - tip) @ rotation.T + pivot for rotation in rotations])
        result = calibrate_tip(self.model.markers, measured)
        assert_allclose(result.tip, tip, atol=1e-10, rtol=0)
        assert_allclose(result.pivot_tracker, pivot, atol=1e-10, rtol=0)

    def test_unobservable_pivot_is_rejected(self):
        repeated = np.tile(self.model.markers, (6, 1, 1))
        for observed in (repeated, self.spin.observations):
            with self.subTest(shape=observed.shape), self.assertRaisesRegex(ValueError, "multiple axes"):
                calibrate_tip(self.model.markers, observed)

    def test_stationary_spin_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "plane"):
            calibrate_axis(self.model.markers, np.tile(self.model.markers, (6, 1, 1)),
                           calibrated_tip=[0, 0, 0])

    def test_on_axis_marker_can_be_ignored(self):
        model = ToolModel(np.array([[0., 180., 0.], [80., 180., 0.], [20., 240., 0.]]),
                          np.zeros(3), np.array([0., -1., 0.]))
        spin = simulate_spin_motion(model)
        result = calibrate_axis(model.markers, spin.observations, calibrated_tip=[0, 0, 0])
        self.assertEqual(result.used_marker_indices, (1, 2))
        assert_allclose(result.axis, [0, -1, 0], atol=1e-12)

    def test_axis_sign_requires_tip_separation(self):
        with self.assertRaisesRegex(ValueError, "sign"):
            calibrate_axis(self.model.markers, self.spin.observations,
                           calibrated_tip=self.model.markers.mean(axis=0))

    def test_malformed_observations_are_rejected(self):
        bad = self.pivot.observations.copy()
        bad[0, 0, 0] = np.nan
        for observations in (bad, self.pivot.observations[:2], np.zeros((4, 2, 3)), np.zeros((4, 3))):
            with self.subTest(shape=observations.shape), self.assertRaises(ValueError):
                calibrate_tip(self.model.markers, observations)

    def test_invalid_simulation_setup_is_rejected(self):
        for count in (2, 4.5, True):
            with self.subTest(count=count), self.assertRaises(ValueError):
                simulate_pivot_motion(self.model, pose_count=count)
        with self.assertRaises(ValueError):
            simulate_pivot_motion(self.model, cone_half_angle_degrees=0)
        with self.assertRaises(ValueError):
            simulate_spin_motion(self.model, arc_degrees=0)
        with self.assertRaises(ValueError):
            simulate_spin_motion(self.model, tracker_rotation=np.diag([1, 1, -1]))

    def test_noisy_calibration_is_repeatable_and_has_finite_nonzero_errors(self):
        rng = np.random.default_rng(24)
        sigma = sigma_for_radial_bound()
        pivot = add_marker_noise(self.pivot.observations, sigma, rng=rng)
        spin = add_marker_noise(self.spin.observations, sigma, rng=rng)
        first = calibrate_tool(self.model.markers, pivot, spin)
        second = calibrate_tool(self.model.markers, pivot, spin)
        assert_allclose(first.tip, second.tip, atol=0, rtol=0)
        assert_allclose(first.axis, second.axis, atol=0, rtol=0)
        self.assertTrue(np.isfinite(first.tip).all())
        self.assertAlmostEqual(np.linalg.norm(first.axis), 1)
        self.assertGreater(np.linalg.norm(first.tip), 1e-8)
        self.assertGreater(direction_error_degrees(first.axis, self.model.axis), 1e-8)
        self.assertGreater(first.tip_calibration.rms_residual, 0)

    def test_calibration_flows_into_existing_registration_and_planning(self):
        calibrated = calibrate_tool(self.model.markers, self.pivot.observations, self.spin.observations)
        ct_markers = [[-170, 15, 2], [-170, 75, 2], [-230, 15, 2]]
        tool = register_tool(self.model.markers, ct_markers, calibrated.tip, calibrated.axis)
        plan = plan_biopsy(tool.tip, tool.axis, [110, -2, 2], 5)
        assert_allclose([plan.targeting_error, plan.entry_depth, plan.tissue_core_length,
                         plan.insertion_depth], [3, 96, 8, 104], atol=1e-9, rtol=0)


if __name__ == "__main__":
    unittest.main()
