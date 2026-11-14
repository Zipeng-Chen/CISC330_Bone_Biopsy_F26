"""Independent hand-predicted cases and failure checks for rigid registration."""

import unittest

import numpy as np
from numpy.testing import assert_allclose

from bone_biopsy.models import assignment_tool_model
from bone_biopsy.registration import register_rigid, register_tool


class RegistrationTests(unittest.TestCase):
    def setUp(self):
        # Independent fixture: Figure 2 centimetres converted to millimetres.
        self.markers = np.array([[20., 180., 0.], [80., 180., 0.], [20., 240., 0.]])
        self.tip = np.zeros(3)
        self.axis = np.array([0., -1., 0.])

    def check_pose(self, ct_markers, expected_r, expected_t, expected_tip, expected_axis):
        result = register_tool(self.markers, ct_markers, self.tip, self.axis)
        assert_allclose(result.transform.rotation, expected_r, atol=1e-12, rtol=0)
        assert_allclose(result.transform.translation, expected_t, atol=1e-12, rtol=0)
        assert_allclose(result.tip, expected_tip, atol=1e-12, rtol=0)
        assert_allclose(result.axis, expected_axis, atol=1e-12, rtol=0)
        assert_allclose(result.transform.transform_points(self.markers), ct_markers,
                        atol=1e-12, rtol=0)
        self.assertLess(result.transform.rms_residual, 1e-12)

    def test_required_translation(self):
        # Add (10, -5, 2) to each marker and tip; axis is unchanged.
        self.check_pose([[30, 175, 2], [90, 175, 2], [30, 235, 2]],
                        np.eye(3), [10, -5, 2], [10, -5, 2], [0, -1, 0])

    def test_required_rotation(self):
        # +90 degrees about z maps (x,y,z) to (-y,x,z).
        self.check_pose([[-180, 20, 0], [-180, 80, 0], [-240, 20, 0]],
                        [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
                        [0, 0, 0], [0, 0, 0], [1, 0, 0])

    def test_required_rotation_and_translation(self):
        self.check_pose([[-170, 15, 2], [-170, 75, 2], [-230, 15, 2]],
                        [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
                        [10, -5, 2], [10, -5, 2], [1, 0, 0])

    def test_nonzero_tip_and_out_of_plane_rotation(self):
        # +90 degrees about x maps (x,y,z) to (x,-z,y), then add (3,4,5).
        ct_markers = [[23, 4, 185], [83, 4, 185], [23, 4, 245]]
        result = register_tool(self.markers, ct_markers, [1, 2, 3], [0, -7, 0])
        assert_allclose(result.tip, [4, 1, 7], atol=1e-12, rtol=0)
        assert_allclose(result.axis, [0, 0, -1], atol=1e-12, rtol=0)

    def test_identity_and_reverse_transform(self):
        identity = register_rigid(self.markers, self.markers)
        assert_allclose(identity.rotation, np.eye(3), atol=1e-12, rtol=0)
        ct = np.array([[-170, 15, 2], [-170, 75, 2], [-230, 15, 2]])
        forward = register_rigid(self.markers, ct)
        reverse = register_rigid(ct, self.markers)
        point = np.array([2., 3., 7.])
        assert_allclose(reverse.transform_point(forward.transform_point(point)),
                        point, atol=1e-12, rtol=0)

    def test_noisy_markers_still_produce_proper_rotation(self):
        noisy = self.markers + np.array([[.1, -.1, .05], [-.05, .1, .02], [.03, .04, -.08]])
        transform = register_rigid(self.markers, noisy)
        assert_allclose(transform.rotation.T @ transform.rotation, np.eye(3), atol=1e-12)
        self.assertAlmostEqual(np.linalg.det(transform.rotation), 1.)
        self.assertGreater(transform.rms_residual, 0.)
        # The fitted pose should improve on the identity candidate.
        identity_rms = np.sqrt(np.mean(np.sum((self.markers - noisy) ** 2, axis=1)))
        self.assertLessEqual(transform.rms_residual, identity_rms + 1e-12)

    def test_reflection_is_not_returned_as_rotation(self):
        tetrahedron = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
        reflected = tetrahedron * [-1, 1, 1]
        transform = register_rigid(tetrahedron, reflected)
        self.assertAlmostEqual(np.linalg.det(transform.rotation), 1.)
        self.assertGreater(transform.rms_residual, .1)

    def test_degenerate_and_invalid_marker_sets_rejected(self):
        invalid_sets = [np.zeros((3, 3)), [[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                        [[0, 0, 0], [1, 0, 0], [2, 1e-12, 0]],
                        [[0, 0, 0], [1, 0, 0]],
                        [[0, 0, 0], [1, 0, 0], [0, np.nan, 0]]]
        for invalid in invalid_sets:
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    register_rigid(invalid, invalid)
        with self.assertRaises(ValueError):
            register_rigid(self.markers, np.zeros((4, 3)))
        with self.assertRaises(ValueError):
            register_rigid(self.markers, np.zeros((3, 3)))

    def test_zero_axis_and_bad_tip_rejected(self):
        with self.assertRaises(ValueError):
            register_tool(self.markers, self.markers, self.tip, [0, 0, 0])
        with self.assertRaises(ValueError):
            register_tool(self.markers, self.markers, [np.inf, 0, 0], self.axis)


class ToolModelTests(unittest.TestCase):
    def test_default_converts_confirmed_centimetres_to_millimetres(self):
        model = assignment_tool_model()
        assert_allclose(model.markers, [[20, 180, 0], [80, 180, 0], [20, 240, 0]],
                        atol=0, rtol=0)
        assert_allclose(model.tip, [0, 0, 0], atol=0, rtol=0)
        assert_allclose(model.axis, [0, -1, 0], atol=0, rtol=0)
        self.assertEqual(np.linalg.norm(model.markers[1] - model.markers[0]), 60)

    def test_explicit_conversion_remains_compatible_and_axis_is_not_scaled(self):
        model = assignment_tool_model()
        explicit = assignment_tool_model(mm_per_figure_unit=10)
        synthetic = assignment_tool_model(mm_per_figure_unit=1)
        assert_allclose(explicit.markers, model.markers)
        assert_allclose(model.markers, 10 * synthetic.markers)
        assert_allclose(model.axis, synthetic.axis)

    def test_invalid_model_scale_is_rejected(self):
        for scale in [0, -1, np.inf, np.nan]:
            with self.subTest(scale=scale), self.assertRaises(ValueError):
                assignment_tool_model(mm_per_figure_unit=scale)


if __name__ == "__main__":
    unittest.main()
