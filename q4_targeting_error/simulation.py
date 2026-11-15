"""Q4: Targeting Error Simulation for CISC/CMPE 330 A2.

Simulates the effect of:
1. CT marker-localization uncertainty (fixed tracker noise).
2. Procedural target displacement (variable Gaussian sigma).
3. Optional uncertainty from noisy Q1 tool calibration.

All distances are in millimetres.

The targeting error is the shortest perpendicular distance from the
displaced target to the registered needle axis, as defined in Figure 6.
"""

from dataclasses import dataclass

import numpy as np

from bone_biopsy import (
    assignment_tool_model,
    add_marker_noise,
    calibrate_tool,
    distance_point_to_line,
    register_tool,
    sigma_for_radial_bound,
    simulate_pivot_motion,
    simulate_spin_motion,
)

from utils.geometry import rotation_about_axis


# ------------------------------------------------------------
# Assignment constants
# ------------------------------------------------------------

INSERTION_DEPTH_MM = 100.0
TARGET_TOLERANCE_MM = 3.0
REQUIRED_SUCCESS_PROBABILITY = 0.95

TRACKER_SIGMA_MM = sigma_for_radial_bound(
    bound_mm=1.0,
    probability=0.95,
)


# ------------------------------------------------------------
# Data structures
# ------------------------------------------------------------

@dataclass(frozen=True)
class IdealScene:
    """Ground-truth biopsy geometry in the CT coordinate frame."""

    model: object
    ct_markers: np.ndarray
    tip: np.ndarray
    axis: np.ndarray
    target: np.ndarray
    insertion_depth_mm: float


@dataclass(frozen=True)
class TrialResult:
    """Results from one simulated biopsy procedure."""

    targeting_error_mm: float
    registered_tip: np.ndarray
    registered_axis: np.ndarray
    displaced_target: np.ndarray
    procedural_sigma_mm: float


# ------------------------------------------------------------
# Ideal scene and exact validation
# ------------------------------------------------------------

def ideal_scene(insertion_depth_mm=INSERTION_DEPTH_MM):
    """Construct the ideal scene described in Figure 6.

    The model and CT frames coincide.

    Therefore:
        R = I
        t = 0

    The planned target is exactly on the ideal needle axis,
    100 mm from the tool tip by default.
    """

    depth = float(insertion_depth_mm)

    if not np.isfinite(depth) or depth <= 0:
        raise ValueError("Insertion depth must be positive and finite.")

    model = assignment_tool_model()

    tip = np.array(model.tip, dtype=float, copy=True)
    axis = np.array(model.axis, dtype=float, copy=True)

    axis /= np.linalg.norm(axis)

    target = tip + depth * axis

    return IdealScene(
        model=model,
        ct_markers=np.array(model.markers, dtype=float, copy=True),
        tip=tip,
        axis=axis,
        target=target,
        insertion_depth_mm=depth,
    )


def exact_error(scene):
    """Verify that the noise-free scene has zero targeting error."""

    registered = register_tool(
        scene.model.markers,
        scene.ct_markers,
        scene.model.tip,
        scene.model.axis,
    )

    return float(
        distance_point_to_line(
            scene.target,
            registered.tip,
            registered.axis,
        )
    )


# ------------------------------------------------------------
# Calibration observations for combined experiment
# ------------------------------------------------------------

def calibration_observations(scene, poses=60):
    """Generate ground-truth pivot and spin observations.

    These observations are generated once and reused.
    Fresh independent Gaussian noise is added in each trial.
    """

    if isinstance(poses, (bool, np.bool_)) or not isinstance(
        poses, (int, np.integer)
    ) or poses < 3:
        raise ValueError("poses must be an integer >= 3.")

    base_rotation = rotation_about_axis([1, 2, 3], 0.4)

    pivot = simulate_pivot_motion(
        scene.model,
        pose_count=poses,
        tracker_rotation=base_rotation,
        seed=17,
    )

    spin = simulate_spin_motion(
        scene.model,
        pose_count=poses,
        tracker_rotation=base_rotation,
    )

    return pivot.observations, spin.observations


# ------------------------------------------------------------
# Individual Monte Carlo trial
# ------------------------------------------------------------

def trial(
    scene,
    procedural_sigma_mm,
    rng,
    *,
    include_calibration_uncertainty=False,
    calibration_data=None,
):
    """Simulate one biopsy targeting-error trial.

    Primary experiment:
        Exact calibration
        + noisy CT markers
        + procedural target displacement.

    Combined experiment:
        Noisy Q1 calibration
        + noisy CT markers
        + procedural target displacement.
    """

    sigma = float(procedural_sigma_mm)

    if not np.isfinite(sigma) or sigma < 0:
        raise ValueError("Procedural sigma must be nonnegative and finite.")

    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng must be a numpy.random.Generator.")

    # --------------------------------------------------------
    # 1. Tool calibration
    # --------------------------------------------------------

    if include_calibration_uncertainty:

        if calibration_data is None:
            calibration_data = calibration_observations(scene)

        pivot_observations, spin_observations = calibration_data

        noisy_pivot = add_marker_noise(
            pivot_observations,
            TRACKER_SIGMA_MM,
            rng=rng,
        )

        noisy_spin = add_marker_noise(
            spin_observations,
            TRACKER_SIGMA_MM,
            rng=rng,
        )

        calibrated = calibrate_tool(
            scene.model.markers,
            noisy_pivot,
            noisy_spin,
        )

        calibrated_tip = calibrated.tip
        calibrated_axis = calibrated.axis

    else:

        # Exact calibration isolates registration uncertainty.
        calibrated_tip = scene.model.tip
        calibrated_axis = scene.model.axis

    # --------------------------------------------------------
    # 2. Independent CT marker-localization error
    # --------------------------------------------------------

    noisy_ct_markers = add_marker_noise(
        scene.ct_markers,
        TRACKER_SIGMA_MM,
        rng=rng,
    )

    # --------------------------------------------------------
    # 3. Propagate marker error through rigid registration
    # --------------------------------------------------------

    registered = register_tool(
        scene.model.markers,
        noisy_ct_markers,
        calibrated_tip,
        calibrated_axis,
    )

    # --------------------------------------------------------
    # 4. Independent procedural target displacement
    # --------------------------------------------------------

    displacement = rng.normal(
        loc=0.0,
        scale=sigma,
        size=3,
    )

    displaced_target = scene.target + displacement

    # --------------------------------------------------------
    # 5. Compound targeting error
    # --------------------------------------------------------

    targeting_error = distance_point_to_line(
        displaced_target,
        registered.tip,
        registered.axis,
    )

    return TrialResult(
        targeting_error_mm=float(targeting_error),
        registered_tip=registered.tip,
        registered_axis=registered.axis,
        displaced_target=displaced_target,
        procedural_sigma_mm=sigma,
    )


# ------------------------------------------------------------
# Paired trial for comparing calibration uncertainty
# ------------------------------------------------------------

def trial_pair(scene, procedural_sigma_mm, rng, calibration_data=None):
    """Run both experiments with shared CT and procedural noise.

    This paired design makes comparisons more reliable.

    The combined experiment introduces additional independent
    calibration noise, while the CT-marker and target-displacement
    samples remain identical across the two experiments.
    """

    sigma = float(procedural_sigma_mm)

    if not np.isfinite(sigma) or sigma < 0:
        raise ValueError("Procedural sigma must be nonnegative and finite.")

    if calibration_data is None:
        calibration_data = calibration_observations(scene)

    pivot_observations, spin_observations = calibration_data

    # Shared CT marker-localization noise.
    noisy_ct_markers = add_marker_noise(
        scene.ct_markers,
        TRACKER_SIGMA_MM,
        rng=rng,
    )

    # Shared procedural displacement.
    displacement = rng.normal(
        loc=0.0,
        scale=sigma,
        size=3,
    )

    displaced_target = scene.target + displacement

    # Experiment 1: exact calibration.
    exact_registered = register_tool(
        scene.model.markers,
        noisy_ct_markers,
        scene.model.tip,
        scene.model.axis,
    )

    registration_error = distance_point_to_line(
        displaced_target,
        exact_registered.tip,
        exact_registered.axis,
    )

    # Experiment 2: noisy calibration.
    noisy_pivot = add_marker_noise(
        pivot_observations,
        TRACKER_SIGMA_MM,
        rng=rng,
    )

    noisy_spin = add_marker_noise(
        spin_observations,
        TRACKER_SIGMA_MM,
        rng=rng,
    )

    calibrated = calibrate_tool(
        scene.model.markers,
        noisy_pivot,
        noisy_spin,
    )

    combined_registered = register_tool(
        scene.model.markers,
        noisy_ct_markers,
        calibrated.tip,
        calibrated.axis,
    )

    combined_error = distance_point_to_line(
        displaced_target,
        combined_registered.tip,
        combined_registered.axis,
    )

    return float(registration_error), float(combined_error)


# ------------------------------------------------------------
# Monte Carlo simulation
# ------------------------------------------------------------

def simulate(
    scene,
    procedural_sigma_mm,
    trials=10000,
    seed=2026,
    *,
    include_calibration_uncertainty=False,
):
    """Run repeated independent targeting-error trials."""

    if isinstance(trials, (bool, np.bool_)) or not isinstance(
        trials, (int, np.integer)
    ) or trials < 1:
        raise ValueError("trials must be a positive integer.")

    rng = np.random.default_rng(seed)

    calibration_data = None

    if include_calibration_uncertainty:
        calibration_data = calibration_observations(scene)

    errors = np.empty(trials, dtype=float)

    for i in range(trials):

        result = trial(
            scene,
            procedural_sigma_mm,
            rng,
            include_calibration_uncertainty=include_calibration_uncertainty,
            calibration_data=calibration_data,
        )

        errors[i] = result.targeting_error_mm

    return errors


def simulate_paired(
    scene,
    procedural_sigma_mm,
    trials=10000,
    seed=2026,
):
    """Run paired Monte Carlo experiments.

    Returns:
        registration_errors
        combined_errors
    """

    if isinstance(trials, (bool, np.bool_)) or not isinstance(
        trials, (int, np.integer)
    ) or trials < 1:
        raise ValueError("trials must be a positive integer.")

    rng = np.random.default_rng(seed)

    calibration_data = calibration_observations(scene)

    registration_errors = np.empty(trials, dtype=float)
    combined_errors = np.empty(trials, dtype=float)

    for i in range(trials):

        registration_error, combined_error = trial_pair(
            scene,
            procedural_sigma_mm,
            rng,
            calibration_data=calibration_data,
        )

        registration_errors[i] = registration_error
        combined_errors[i] = combined_error

    return registration_errors, combined_errors


# ------------------------------------------------------------
# Statistical analysis
# ------------------------------------------------------------

def summarize(errors, tolerance_mm=TARGET_TOLERANCE_MM):
    """Summarize targeting-error magnitudes and clinical success.

    The Wilson confidence interval quantifies uncertainty in
    the estimated probability of satisfying the 3 mm tolerance.
    """

    errors = np.asarray(errors, dtype=float)

    if errors.ndim != 1 or errors.size == 0:
        raise ValueError("errors must be a nonempty 1D array.")

    if not np.all(np.isfinite(errors)) or np.any(errors < 0):
        raise ValueError("errors must be finite and nonnegative.")

    tolerance = float(tolerance_mm)

    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("tolerance must be nonnegative and finite.")

    n = errors.size

    successes = int(np.count_nonzero(errors <= tolerance))

    success_probability = successes / n

    # Two-sided 95% Wilson binomial confidence interval.
    z = 1.959963984540054

    denominator = 1.0 + z**2 / n

    center = (
        success_probability + z**2 / (2.0 * n)
    ) / denominator

    half_width = (
        z
        * np.sqrt(
            success_probability * (1.0 - success_probability) / n
            + z**2 / (4.0 * n**2)
        )
        / denominator
    )

    return {
        "trials": int(n),
        "mean_mm": float(np.mean(errors)),
        "median_mm": float(np.median(errors)),
        "rms_mm": float(np.sqrt(np.mean(errors**2))),
        "p95_mm": float(np.percentile(errors, 95)),
        "maximum_mm": float(np.max(errors)),
        "successes": successes,
        "success_probability": float(success_probability),
        "success_ci95_lower": float(max(0.0, center - half_width)),
        "success_ci95_upper": float(min(1.0, center + half_width)),
        "meets_requirement": bool(
            success_probability >= REQUIRED_SUCCESS_PROBABILITY
        ),
        "conservatively_meets_requirement": bool(
            max(0.0, center - half_width)
            >= REQUIRED_SUCCESS_PROBABILITY
        ),
    }


# ------------------------------------------------------------
# Procedural sigma sweep
# ------------------------------------------------------------

def sigma_sweep(
    scene,
    sigma_values,
    trials=10000,
    seed=2026,
    *,
    include_calibration_uncertainty=False,
):
    """Sweep procedural displacement sigma while tracker noise stays fixed."""

    results = []

    for sigma in sigma_values:

        sigma = float(sigma)

        errors = simulate(
            scene,
            procedural_sigma_mm=sigma,
            trials=trials,
            seed=seed,
            include_calibration_uncertainty=include_calibration_uncertainty,
        )

        stats = summarize(errors)

        results.append({
            "procedural_sigma_mm": sigma,
            "experiment": (
                "combined_uncertainty"
                if include_calibration_uncertainty
                else "registration_only"
            ),
            **stats,
        })

    return results


def paired_sigma_sweep(
    scene,
    sigma_values,
    trials=10000,
    seed=2026,
):
    """Compare both experiments across procedural sigma values.

    Tracker sigma remains fixed throughout the sweep.

    Each sigma uses the same seed to improve reproducibility.
    """

    results = []

    for sigma in sigma_values:

        sigma = float(sigma)

        registration_errors, combined_errors = simulate_paired(
            scene,
            procedural_sigma_mm=sigma,
            trials=trials,
            seed=seed,
        )

        results.append({
            "procedural_sigma_mm": sigma,
            "registration_only": summarize(registration_errors),
            "combined_uncertainty": summarize(combined_errors),
        })

    return results