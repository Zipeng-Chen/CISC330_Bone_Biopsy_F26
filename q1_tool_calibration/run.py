"""Validate exact calibration and quantify errors from repeated noisy trials.

Run from the repository root: python -m q1_tool_calibration.run
Optional: --plot-dir q1_tool_calibration/outputs
"""

import argparse
import json
from pathlib import Path
import shlex
import sys

# Resolve imports for direct execution as well as python -m.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from numpy.testing import assert_allclose

from bone_biopsy import (
    add_marker_noise, assignment_tool_model, calibrate_axis, calibrate_tip, calibrate_tool, plan_biopsy,
    register_tool, sigma_for_radial_bound, simulate_pivot_motion, simulate_spin_motion,
)
from utils.geometry import direction_error_degrees, rotation_about_axis

from utils.reporting import run_tests, run_with_log

OUTPUTS = Path(__file__).resolve().parent / "outputs"


def summarize(values):
    """Summarize nonnegative error magnitudes, not signed coordinate errors"""
    values = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(values)), "median": float(np.median(values)),
        "rms": float(np.sqrt(np.mean(values**2))),
        "p95": float(np.percentile(values, 95)), "maximum": float(np.max(values)),
    }


def error_record(calibrated, model):
    return {
        "tip_mm": calibrated.tip.tolist(), "axis": calibrated.axis.tolist(),
        "tip_error_mm": float(np.linalg.norm(calibrated.tip - model.tip)),
        "axis_error_degrees": direction_error_degrees(calibrated.axis, model.axis),
        "pivot_rms_mm": calibrated.tip_calibration.rms_residual,
        "pivot_condition_number": calibrated.tip_calibration.condition_number,
        "spin_plane_rms_mm": calibrated.axis_calibration.plane_rms.tolist(),
        "axis_scatter_degrees": calibrated.axis_calibration.direction_scatter_degrees,
    }


def check_hand_calculations(model_markers, calibrated_tip):
    """Check the small pivot and spin examples worked out in q1_tool_calibration/README.md."""
    # These four tracker frames are calculated directly from p=(1,2,3), q=(10,20,30).
    pivot_markers = np.array([
        [[29, 198, 27], [89, 198, 27], [29, 258, 27]],
        [[29, 23, 208], [89, 23, 208], [29, 23, 268]],
        [[7, 198, 11], [7, 198, -49], [7, 258, 11]],
        [[-168, 39, 27], [-168, 99, 27], [-228, 39, 27]],
    ], dtype=float)
    tip_result = calibrate_tip(model_markers, pivot_markers)
    assert_allclose(tip_result.tip, [1, 2, 3], atol=1e-10, rtol=0)
    assert_allclose(tip_result.pivot_tracker, [10, 20, 30], atol=1e-10, rtol=0)

    # Quarter turns around y keep each marker's y coordinate fixed.
    spin_markers = np.array([
        [[20, 180, 0], [80, 180, 0], [20, 240, 0]],
        [[0, 180, -20], [0, 180, -80], [0, 240, -20]],
        [[-20, 180, 0], [-80, 180, 0], [-20, 240, 0]],
        [[0, 180, 20], [0, 180, 80], [0, 240, 20]],
    ], dtype=float)
    axis_result = calibrate_axis(model_markers, spin_markers, calibrated_tip=calibrated_tip)
    assert_allclose(axis_result.axis, [0, -1, 0], atol=1e-12, rtol=0)
    return {
        "tip_expected_mm": [1, 2, 3], "tip_computed_mm": tip_result.tip.tolist(),
        "pivot_expected_mm": [10, 20, 30],
        "pivot_computed_mm": tip_result.pivot_tracker.tolist(),
        "tip_difference_mm": float(np.linalg.norm(tip_result.tip - [1, 2, 3])),
        "axis_expected": [0, -1, 0], "axis_computed": axis_result.axis.tolist(),
        "axis_difference_degrees": direction_error_degrees(axis_result.axis, [0, -1, 0]),
    }


def run_experiment(*, trials=500, poses=60, seed=2026):
    """Repeat calibration with fresh independent noise on fixed, known motions"""
    if (isinstance(trials, (bool, np.bool_))
            or not isinstance(trials, (int, np.integer)) or trials < 1):
        raise ValueError("trials must be a positive integer")
    model = assignment_tool_model()
    base_rotation = rotation_about_axis([1, 2, 3], .4)
    # Generate constrained motions using truth; pass only marker observations to solvers.
    pivot = simulate_pivot_motion(model, pose_count=poses, tracker_rotation=base_rotation, seed=17)
    spin = simulate_spin_motion(model, pose_count=poses, tracker_rotation=base_rotation)
    exact = calibrate_tool(model.markers, pivot.observations, spin.observations)
    assert_allclose(exact.tip, model.tip, atol=1e-10, rtol=0)
    assert_allclose(exact.axis, model.axis, atol=1e-12, rtol=0)
    hand_checks = check_hand_calculations(model.markers, exact.tip)

    sigma = sigma_for_radial_bound(bound_mm=1.0, probability=.95)
    rng = np.random.default_rng(seed)
    tip_errors, axis_errors, tip_estimates = [], [], []
    within_bound = marker_count = 0
    first_noisy = None
    for trial in range(trials):
        # Perturb measurements before registration, independently for both movements
        noisy_pivot = add_marker_noise(pivot.observations, sigma, rng=rng)
        noisy_spin = add_marker_noise(spin.observations, sigma, rng=rng)
        calibrated = calibrate_tool(model.markers, noisy_pivot, noisy_spin)
        record = error_record(calibrated, model)
        if trial == 0:
            first_noisy = record
        tip_errors.append(record["tip_error_mm"])
        axis_errors.append(record["axis_error_degrees"])
        tip_estimates.append(calibrated.tip)
        # Check the radial marker bound using the same draws as the calibration trials
        for observed, truth in ((noisy_pivot, pivot.observations), (noisy_spin, spin.observations)):
            lengths = np.linalg.norm(observed - truth, axis=-1)
            within_bound += int(np.count_nonzero(lengths < 1))
            marker_count += lengths.size

    # Feed the calibrated quantities through the existing registration/planning interface
    tool = register_tool(model.markers, [[-170, 15, 2], [-170, 75, 2], [-230, 15, 2]],
                         exact.tip, exact.axis)
    plan = plan_biopsy(tool.tip, tool.axis, [110, -2, 2], 5)
    pipeline = [plan.targeting_error, plan.entry_depth, plan.tissue_core_length, plan.insertion_depth]
    assert_allclose(pipeline, [3, 96, 8, 104], atol=1e-9, rtol=0)
    report = {
        "configuration": {
            "trials": trials, "pivot_poses": poses, "spin_poses": poses,
            "noise_seed": seed, "motion_seed": 17, "length_unit": "mm",
            "pivot_cone_half_angle_degrees": 40, "spin_arc_degrees": 360,
            "tracker_rotation_axis": [1, 2, 3], "tracker_rotation_angle_radians": .4,
            "pivot_and_spin_tip_tracker_mm": [25, -10, 40],
            "marker_sigma_per_coordinate_mm": sigma,
            "marker_radial_bound_mm": 1.0, "marker_radial_probability": .95,
        },
        "ground_truth": {"markers_mm": model.markers.tolist(), "tip_mm": model.tip.tolist(),
                         "axis": model.axis.tolist()},
        "hand_calculations": hand_checks,
        "exact": error_record(exact, model), "first_noisy_trial": first_noisy,
        "tip_error_mm": summarize(tip_errors), "axis_error_degrees": summarize(axis_errors),
        "tip_bias_vector_mm": (np.mean(tip_estimates, axis=0) - model.tip).tolist(),
        "marker_noise_check": {
            "observations": marker_count, "fraction_below_1_mm": within_bound / marker_count,
            "expected_binomial_standard_error": float(np.sqrt(.95 * .05 / marker_count)),
        },
        "exact_pipeline_error_entry_core_depth_mm": pipeline,
        "trial_tip_errors_mm": tip_errors, "trial_axis_errors_degrees": axis_errors,
    }
    return report, model, pivot, spin


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=500)
    parser.add_argument("--poses", type=int, default=60, help="poses per calibration movement")
    parser.add_argument("--seed", type=int, default=2026, help="marker-noise random seed")
    parser.add_argument("--output", type=Path, default=OUTPUTS / "calibration-results.json",
                        help="JSON results file (default: this question's outputs folder)")
    parser.add_argument("--plot-dir", type=Path, help="optional figures; requires matplotlib")
    args = parser.parse_args()
    command = "python -m q1_tool_calibration.run"
    if sys.argv[1:]:
        command += " " + shlex.join(sys.argv[1:])
    run_with_log(lambda: run(args), OUTPUTS / "run_log.txt", command)


def run(args):
    report, model, pivot, spin = run_experiment(trials=args.trials, poses=args.poses, seed=args.seed)
    print("CALIBRATION: Figure 2 cm converted to mm")
    print(f"Poses: {args.poses} pivot + {args.poses} spin; noisy trials: {args.trials}")
    print(f"Marker sigma per coordinate: {report['configuration']['marker_sigma_per_coordinate_mm']:.9f} mm")
    hand = report["hand_calculations"]
    print("\nHAND-CALCULATED CHECKS (worked derivations in q1_tool_calibration/README.md):")
    print(f"  tip expected {hand['tip_expected_mm']}; computed {np.array(hand['tip_computed_mm'])} mm")
    print(f"  pivot expected {hand['pivot_expected_mm']}; computed {np.array(hand['pivot_computed_mm'])} mm")
    print(f"  tip difference: {hand['tip_difference_mm']:.3e} mm")
    print(f"  axis expected {hand['axis_expected']}; computed {np.array(hand['axis_computed'])}")
    print(f"  axis difference: {hand['axis_difference_degrees']:.3e} degrees")
    for name in ("exact", "first_noisy_trial"):
        result = report[name]
        print(f"\n{name}:")
        print(f"  tip (mm): {np.array(result['tip_mm'])}")
        print(f"  axis: {np.array(result['axis'])}")
        print(f"  tip error: {result['tip_error_mm']:.9g} mm")
        print(f"  axis error: {result['axis_error_degrees']:.9g} degrees")
    print("\nRepeated noisy calibration:")
    for key in ("tip_error_mm", "axis_error_degrees"):
        print(f"  {key}: " + ", ".join(f"{name}={value:.6f}" for name, value in report[key].items()))
    print(f"Marker errors below 1 mm: {report['marker_noise_check']['fraction_below_1_mm']:.4%}")
    print("PASS: exact calibration and calibrated registration/planning match predictions.")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        print(f"Saved results: {args.output.name}")
    if args.plot_dir:
        from q1_tool_calibration.plot_results import save_calibration_figures
        save_calibration_figures(model, pivot, spin, report, args.plot_dir)
        print(f"Saved calibration figures in: {args.plot_dir}")

    run_tests("q1_tool_calibration.test_calibration", OUTPUTS / "test_report.txt")


if __name__ == "__main__":
    main()
