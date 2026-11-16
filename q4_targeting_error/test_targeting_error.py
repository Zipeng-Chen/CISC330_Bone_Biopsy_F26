"""Independent geometry and probability checks for Q4.

Hand predictions and the purpose of each check are explained in
verification/README.md. Tests do not write simulation output files.
"""

from dataclasses import replace
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose

from q4_targeting_error import simulation as sim
from q4_targeting_error.run import maximum_tolerable_sigma


class TargetingErrorTests(unittest.TestCase):
    def setUp(self):
        self.scene = sim.ideal_scene()

    def test_ideal_scene_has_correct_units_and_zero_error(self):
        assert_allclose(self.scene.ct_markers, [[20, 180, 0], [80, 180, 0], [20, 240, 0]])
        assert_allclose(self.scene.target, [0, -100, 0])
        self.assertLess(sim.exact_error(self.scene), 1e-10)
        # The noise value is a coordinate sigma for a 3D 95% radial bound.
        self.assertAlmostEqual(sim.TRACKER_SIGMA_MM, .3577198742584623, places=12)

    def test_known_marker_shift_uses_perpendicular_not_endpoint_distance(self):
        scene = replace(self.scene, target=np.array([0., -104., 0.]))
        # Shift every CT marker 3 mm sideways; real registration must recover it.
        observations = scene.ct_markers + [3, 0, 0]
        with patch.object(sim, 'add_marker_noise', return_value=observations):
            result = sim.trial(scene, 0, np.random.default_rng(10))
        assert_allclose(result.registered_tip, [3, 0, 0], atol=1e-10)
        assert_allclose(result.registered_axis, [0, -1, 0], atol=1e-12)
        self.assertAlmostEqual(result.targeting_error_mm, 3, places=10)
        endpoint = result.registered_tip + 100 * result.registered_axis
        self.assertAlmostEqual(np.linalg.norm(result.displaced_target - endpoint), 5, places=10)

    def test_known_marker_rotation_gives_hand_predicted_lateral_error(self):
        angle = np.deg2rad(1.)
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                             [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        observations = self.scene.ct_markers @ rotation.T
        with patch.object(sim, 'add_marker_noise', return_value=observations):
            result = sim.trial(self.scene, 0, np.random.default_rng(11))
        # A line tilted one degree from the target direction misses by 100 sin(1 degree).
        self.assertAlmostEqual(result.targeting_error_mm, 100 * np.sin(angle), places=10)

    def test_zero_target_displacement_keeps_instrument_noise_active(self):
        rng = np.random.default_rng(12)
        results = [sim.trial(self.scene, 0, rng) for _ in range(100)]
        assert_allclose([r.displaced_target for r in results],
                        np.tile(self.scene.target, (100, 1)), atol=0, rtol=0)
        self.assertGreater(np.std([r.targeting_error_mm for r in results]), .1)
        self.assertGreater(np.max([r.targeting_error_mm for r in results]), 1.)

    def test_target_only_noise_agrees_with_rayleigh_prediction(self):
        # Only this isolated case has an exact Rayleigh model: fixed true line,
        # with two independent unit-variance lateral target-error components.
        with patch.object(sim, 'TRACKER_SIGMA_MM', 0.):
            errors = sim.simulate(self.scene, 1., trials=4000, seed=13)
        expected_success = 1 - np.exp(-3**2 / 2)
        standard_error = np.sqrt(expected_success * (1 - expected_success) / len(errors))
        self.assertLess(abs(np.mean(errors <= 3) - expected_success), 5 * standard_error)
        self.assertAlmostEqual(np.percentile(errors, 95), np.sqrt(-2*np.log(.05)), delta=.1)

    def test_paired_trials_share_ct_and_target_noise(self):
        exact = SimpleNamespace(tip=self.scene.model.tip, axis=self.scene.model.axis)
        data = sim.calibration_observations(self.scene)
        # Removing calibration error must make each pair identical even with other noise.
        with patch.object(sim, 'calibrate_tool', return_value=exact):
            rng = np.random.default_rng(14)
            for _ in range(20):
                primary, combined = sim.trial_pair(self.scene, 1., rng, data)
                self.assertAlmostEqual(primary, combined, places=12)

    def test_combined_trial_propagates_calibrated_tip_and_axis(self):
        angle = np.deg2rad(1.)
        fixtures = [([3., 0, 0], [0., -1, 0], 3.),
                    ([0., 0, 0], [np.sin(angle), -np.cos(angle), 0], 100*np.sin(angle))]
        data = sim.calibration_observations(self.scene)
        for tip, axis, expected in fixtures:
            with self.subTest(expected=expected):
                calibrated = SimpleNamespace(tip=np.array(tip), axis=np.array(axis))
                with patch.object(sim, 'TRACKER_SIGMA_MM', 0.), patch.object(
                    sim, 'calibrate_tool', return_value=calibrated
                ):
                    primary, combined = sim.trial_pair(self.scene, 0., np.random.default_rng(15), data)
                assert_allclose([primary, combined], [0, expected], atol=1e-10)

    def test_actual_noisy_calibration_is_repeatable(self):
        # Exercise the real calibration solver as well as the controlled fixtures above.
        first = sim.simulate_paired(self.scene, .5, trials=8, seed=16)
        second = sim.simulate_paired(self.scene, .5, trials=8, seed=16)
        assert_allclose(first, second, atol=0, rtol=0)
        self.assertTrue(np.all(np.isfinite(first)))
        self.assertTrue(np.all(np.asarray(first) >= 0))
        self.assertGreater(np.max(np.abs(first[0] - first[1])), 1e-6)

    def test_success_boundary_and_wilson_interval(self):
        stats = sim.summarize([0., 3., 3.01, 4.])
        self.assertEqual(stats['successes'], 2)
        self.assertEqual(stats['success_probability'], .5)
        # For two successes in four trials, the two-sided 95% Wilson interval is:
        assert_allclose([stats['success_ci95_lower'], stats['success_ci95_upper']],
                        [.15003898915214947, .8499610108478506], atol=1e-12)
        self.assertFalse(stats['meets_requirement'])

    def test_no_feasible_sigma_differs_from_zero_passing_sigma(self):
        failed = sim.summarize([4.] * 100)
        borderline = sim.summarize([0.] * 95 + [4.] * 5)
        rows = [{'procedural_sigma_mm': sigma, 'registration_only': failed,
                 'combined_uncertainty': None} for sigma in [0., .25, .5]]
        self.assertIsNone(maximum_tolerable_sigma(rows, 'registration_only'))
        rows[0]['registration_only'] = borderline
        self.assertEqual(maximum_tolerable_sigma(rows, 'registration_only'), 0.)
        # Observing exactly 95% success does not put the confidence interval above 95%.
        self.assertIsNone(maximum_tolerable_sigma(rows, 'registration_only', conservative=True))
        self.assertIsNone(maximum_tolerable_sigma(rows, 'combined_uncertainty'))

    def test_sigma_sweep_changes_target_noise_without_changing_marker_noise(self):
        captured = []
        original_trial = sim.trial

        def record_trial(*args, **kwargs):
            result = original_trial(*args, **kwargs)
            captured.append(result)
            return result

        with patch.object(sim, 'trial', side_effect=record_trial):
            rows = sim.sigma_sweep(self.scene, [0., 1.], trials=6, seed=17)
        self.assertEqual([r['procedural_sigma_mm'] for r in rows], [0., 1.])
        # The same seed reuses the same instrument draws across sigma levels.
        for baseline, moved in zip(captured[:6], captured[6:]):
            assert_allclose(baseline.registered_tip, moved.registered_tip, atol=0, rtol=0)
            assert_allclose(baseline.registered_axis, moved.registered_axis, atol=0, rtol=0)
            self.assertGreater(np.linalg.norm(moved.displaced_target - baseline.displaced_target), 0.)

    def test_invalid_inputs_are_rejected(self):
        for value in [-1., np.nan, np.inf]:
            with self.subTest(sigma=value), self.assertRaises(ValueError):
                sim.trial(self.scene, value, np.random.default_rng(18))
        for values in [[], [-1.], [np.nan], [[1., 2.]]]:
            with self.subTest(errors=values), self.assertRaises(ValueError):
                sim.summarize(values)
        with self.assertRaises(ValueError):
            sim.simulate(self.scene, 0., trials=0)
        with self.assertRaises(ValueError):
            sim.ideal_scene(0)


if __name__ == '__main__':
    unittest.main()
