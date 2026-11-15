
"""Q4: Targeting Error Simulation.

Run from the repository root:

    python -m q4_targeting_error.run

Quick test:

    python -m q4_targeting_error.run --trials 1000 --combined-trials 100 --step 0.5

Main experiment:

    python -m q4_targeting_error.run --trials 10000 --combined-trials 0 --step 0.25

The primary experiment uses exact tool calibration and noisy CT marker
registration. The optional combined experiment additionally includes
uncertainty from Q1 tool calibration.

Instrument/tracker uncertainty remains fixed. Only procedural target
displacement sigma is swept.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from q4_targeting_error.simulation import (
    INSERTION_DEPTH_MM,
    REQUIRED_SUCCESS_PROBABILITY,
    TARGET_TOLERANCE_MM,
    TRACKER_SIGMA_MM,
    exact_error,
    ideal_scene,
    simulate,
    simulate_paired,
    summarize,
)


OUTPUTS = Path(__file__).resolve().parent / "outputs"


# ------------------------------------------------------------
# Maximum tolerable procedural sigma
# ------------------------------------------------------------

def maximum_tolerable_sigma(
    results,
    experiment,
    *,
    conservative=False,
):
    """Return the largest tested sigma meeting the clinical criterion.

    None means no tested sigma meets the criterion.
    """

    criterion = (
        "conservatively_meets_requirement"
        if conservative
        else "meets_requirement"
    )

    passing = [
        row["procedural_sigma_mm"]
        for row in results
        if (
            row[experiment] is not None
            and row[experiment][criterion]
        )
    ]

    return max(passing) if passing else None


# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

def save_json(report, path):
    """Save the complete experiment report as JSON."""

    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as file:
        json.dump(
            report,
            file,
            indent=2,
            allow_nan=False,
        )

    return path


def save_csv(results, path):
    """Save sigma-sweep results in CSV format."""

    path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "procedural_sigma_mm",
        "experiment",
        "trials",
        "mean_mm",
        "median_mm",
        "rms_mm",
        "p95_mm",
        "maximum_mm",
        "successes",
        "success_probability",
        "success_ci95_lower",
        "success_ci95_upper",
        "meets_requirement",
        "conservatively_meets_requirement",
    ]

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames,
        )

        writer.writeheader()

        for row in results:
            sigma = row["procedural_sigma_mm"]

            for experiment in (
                "registration_only",
                "combined_uncertainty",
            ):
                stats = row[experiment]

                if stats is None:
                    continue

                writer.writerow({
                    "procedural_sigma_mm": sigma,
                    "experiment": experiment,
                    **stats,
                })

    return path


# ------------------------------------------------------------
# Print results
# ------------------------------------------------------------

def print_results_table(results):
    """Print targeting-error statistics for each sigma."""

    print("\nPROCEDURAL SIGMA SWEEP")
    print("-" * 92)

    print(
        f"{'Sigma':>8} "
        f"{'Primary p95':>13} "
        f"{'Primary success':>17} "
        f"{'Combined p95':>14} "
        f"{'Combined success':>18}"
    )

    print("-" * 92)

    for row in results:
        sigma = row["procedural_sigma_mm"]

        primary = row["registration_only"]
        combined = row["combined_uncertainty"]

        if combined is None:
            combined_p95 = "N/A"
            combined_success = "N/A"
        else:
            combined_p95 = f"{combined['p95_mm']:.3f}"
            combined_success = (
                f"{combined['success_probability']:.2%}"
            )

        print(
            f"{sigma:8.2f} "
            f"{primary['p95_mm']:>13.3f} "
            f"{primary['success_probability']:>17.2%} "
            f"{combined_p95:>14} "
            f"{combined_success:>18}"
        )


def print_threshold(
    name,
    sigma,
    conservative_sigma,
):
    """Print the estimated maximum tolerable procedural sigma."""

    print(f"\n{name}:")

    if sigma is None:
        print(
            "  No tested procedural sigma satisfies "
            "the 95% clinical requirement."
        )
        print(
            "  Maximum tolerable procedural sigma: "
            "none feasible on the tested grid."
        )
    else:
        print(
            "  Maximum passing grid sigma: "
            f"{sigma:.3f} mm"
        )

    if conservative_sigma is None:
        print(
            "  No tested sigma passes the conservative "
            "95% confidence-bound criterion."
        )
    else:
        print(
            "  Maximum conservatively passing grid sigma: "
            f"{conservative_sigma:.3f} mm"
        )


# ------------------------------------------------------------
# Main experiment
# ------------------------------------------------------------

def run_experiment(
    *,
    trials=10000,
    combined_trials=500,
    sigma_max_mm=3.0,
    sigma_step_mm=0.25,
    seed=2026,
):
    """Run primary and optional combined Monte Carlo experiments."""

    if isinstance(trials, bool) or not isinstance(trials, int) or trials < 1:
        raise ValueError("trials must be a positive integer.")

    if (
        isinstance(combined_trials, bool)
        or not isinstance(combined_trials, int)
        or combined_trials < 0
    ):
        raise ValueError(
            "combined_trials must be a nonnegative integer."
        )

    if not np.isfinite(sigma_max_mm) or sigma_max_mm < 0:
        raise ValueError(
            "sigma_max_mm must be nonnegative and finite."
        )

    if not np.isfinite(sigma_step_mm) or sigma_step_mm <= 0:
        raise ValueError(
            "sigma_step_mm must be positive and finite."
        )

    # Construct the ideal scene.
    scene = ideal_scene()

    # Validate ideal geometry.
    ideal_targeting_error = exact_error(scene)

    if not np.isclose(
        ideal_targeting_error,
        0.0,
        atol=1e-9,
    ):
        raise AssertionError(
            "The ideal scene should have zero targeting error."
        )

    # Build procedural sigma grid.
    sigma_values = np.arange(
        0.0,
        sigma_max_mm + sigma_step_mm * 0.5,
        sigma_step_mm,
    )

    sigma_values = sigma_values[
        sigma_values <= sigma_max_mm + 1e-10
    ]

    if not np.isclose(
        sigma_values[-1],
        sigma_max_mm,
    ):
        sigma_values = np.append(
            sigma_values,
            sigma_max_mm,
        )

    sigma_values = np.unique(
        np.round(sigma_values, 10)
    )

    results = []

    # --------------------------------------------------------
    # Procedural sigma sweep
    # --------------------------------------------------------

    for index, sigma in enumerate(
        sigma_values,
        start=1,
    ):
        sigma = float(sigma)

        print(
            f"\n[{index}/{len(sigma_values)}] "
            f"Procedural sigma = {sigma:.2f} mm",
            flush=True,
        )

        # Primary experiment:
        # Exact calibration + noisy CT registration
        # + procedural target displacement.
        print(
            f"  Running primary experiment "
            f"({trials} trials)...",
            flush=True,
        )

        primary_errors = simulate(
            scene,
            sigma,
            trials=trials,
            seed=seed,
            include_calibration_uncertainty=False,
        )

        primary_stats = summarize(primary_errors)

        # Optional combined experiment:
        # Noisy calibration + noisy CT registration
        # + procedural target displacement.
        if combined_trials > 0:
            print(
                f"  Running combined experiment "
                f"({combined_trials} trials)...",
                flush=True,
            )

            _, combined_errors = simulate_paired(
                scene,
                sigma,
                trials=combined_trials,
                seed=seed,
            )

            combined_stats = summarize(
                combined_errors
            )

        else:
            combined_stats = None

        results.append({
            "procedural_sigma_mm": sigma,
            "registration_only": primary_stats,
            "combined_uncertainty": combined_stats,
        })

        print(
            f"  Primary: "
            f"mean={primary_stats['mean_mm']:.3f} mm, "
            f"p95={primary_stats['p95_mm']:.3f} mm, "
            f"success={primary_stats['success_probability']:.2%}",
            flush=True,
        )

        if combined_stats is not None:
            print(
                f"  Combined: "
                f"mean={combined_stats['mean_mm']:.3f} mm, "
                f"p95={combined_stats['p95_mm']:.3f} mm, "
                f"success={combined_stats['success_probability']:.2%}",
                flush=True,
            )

    # --------------------------------------------------------
    # Maximum tolerable sigma
    # --------------------------------------------------------

    primary_max = maximum_tolerable_sigma(
        results,
        "registration_only",
    )

    primary_conservative = maximum_tolerable_sigma(
        results,
        "registration_only",
        conservative=True,
    )

    if combined_trials > 0:
        combined_max = maximum_tolerable_sigma(
            results,
            "combined_uncertainty",
        )

        combined_conservative = maximum_tolerable_sigma(
            results,
            "combined_uncertainty",
            conservative=True,
        )
    else:
        combined_max = None
        combined_conservative = None

    # --------------------------------------------------------
    # Complete report
    # --------------------------------------------------------

    report = {
        "configuration": {
            "primary_trials_per_sigma": int(trials),
            "combined_trials_per_sigma": int(combined_trials),
            "random_seed": int(seed),
            "length_unit": "mm",
            "insertion_depth_mm": INSERTION_DEPTH_MM,
            "target_tolerance_mm": TARGET_TOLERANCE_MM,
            "required_success_probability": (
                REQUIRED_SUCCESS_PROBABILITY
            ),
            "tracker_sigma_per_coordinate_mm": (
                TRACKER_SIGMA_MM
            ),
            "tracker_radial_bound_mm": 1.0,
            "tracker_radial_probability": 0.95,
            "procedural_sigma_values_mm": (
                sigma_values.tolist()
            ),
            "calibration_poses_per_motion": 60,
        },
        "ideal_scene": {
            "model_markers_mm": (
                scene.model.markers.tolist()
            ),
            "ct_markers_mm": (
                scene.ct_markers.tolist()
            ),
            "tool_tip_mm": (
                scene.tip.tolist()
            ),
            "tool_axis": (
                scene.axis.tolist()
            ),
            "planned_target_mm": (
                scene.target.tolist()
            ),
            "exact_targeting_error_mm": (
                float(ideal_targeting_error)
            ),
        },
        "baseline_sigma_zero": results[0],
        "sigma_sweep": results,
        "maximum_tolerable_sigma_mm": {
            "registration_only": primary_max,
            "combined_uncertainty": combined_max,
        },
        "maximum_conservative_sigma_mm": {
            "registration_only": primary_conservative,
            "combined_uncertainty": combined_conservative,
        },
        "interpretation": (
            "Tracker localization noise remains fixed throughout "
            "the experiment. Only procedural target displacement "
            "sigma is varied. The primary experiment assumes exact "
            "tool calibration. The optional combined experiment "
            "also includes noisy Q1 calibration. "
            "The experiments may have different trial counts. "
            "A null maximum indicates that no tested sigma "
            "satisfies the clinical criterion. "
            "Threshold estimates are limited to the tested grid."
        ),
    }

    return report


# ------------------------------------------------------------
# CLI
# ------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        "--trials",
        type=int,
        default=10000,
        help=(
            "Primary Monte Carlo trials per sigma "
            "(default: 10000)"
        ),
    )

    parser.add_argument(
        "--combined-trials",
        type=int,
        default=500,
        help=(
            "Combined experiment trials per sigma; "
            "0 skips combined experiment (default: 500)"
        ),
    )

    parser.add_argument(
        "--sigma-max",
        type=float,
        default=3.0,
        help=(
            "Maximum procedural sigma in mm "
            "(default: 3.0)"
        ),
    )

    parser.add_argument(
        "--step",
        type=float,
        default=0.25,
        help=(
            "Procedural sigma increment in mm "
            "(default: 0.25)"
        ),
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=2026,
        help=(
            "Random seed (default: 2026)"
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUTS,
        help=(
            "Directory for generated CSV and JSON files"
        ),
    )

    args = parser.parse_args()

    print("=" * 72)
    print("Q4: TARGETING ERROR SIMULATION")
    print("=" * 72)

    print(
        f"\nInsertion depth: "
        f"{INSERTION_DEPTH_MM:.1f} mm"
    )

    print(
        f"Clinical tolerance: "
        f"{TARGET_TOLERANCE_MM:.1f} mm"
    )

    print(
        f"Required success probability: "
        f"{REQUIRED_SUCCESS_PROBABILITY:.1%}"
    )

    print(
        "Fixed tracker sigma per coordinate: "
        f"{TRACKER_SIGMA_MM:.9f} mm"
    )

    print(
        f"Primary trials per sigma: "
        f"{args.trials}"
    )

    print(
        f"Combined trials per sigma: "
        f"{args.combined_trials}"
    )

    print(
        f"Procedural sigma range: "
        f"0 to {args.sigma_max} mm"
    )

    print(
        f"Sigma step: {args.step} mm"
    )

    report = run_experiment(
        trials=args.trials,
        combined_trials=args.combined_trials,
        sigma_max_mm=args.sigma_max,
        sigma_step_mm=args.step,
        seed=args.seed,
    )

    print(
        "\nExact ideal targeting error: "
        f"{report['ideal_scene']['exact_targeting_error_mm']:.3e} mm"
    )

    print_results_table(
        report["sigma_sweep"]
    )

    # --------------------------------------------------------
    # Baseline at procedural sigma = 0
    # --------------------------------------------------------

    baseline = report["baseline_sigma_zero"]

    print(
        "\nTRACKER-ONLY BASELINE "
        "(PROCEDURAL SIGMA = 0)"
    )

    print("-" * 60)

    for name, key in (
        ("Primary", "registration_only"),
        ("Combined", "combined_uncertainty"),
    ):
        stats = baseline[key]

        if stats is None:
            print(f"{name}: not run")
            continue

        print(
            f"{name}: "
            f"mean={stats['mean_mm']:.3f} mm, "
            f"p95={stats['p95_mm']:.3f} mm, "
            f"success={stats['success_probability']:.2%}"
        )

    # --------------------------------------------------------
    # Maximum tolerable sigma
    # --------------------------------------------------------

    print(
        "\nMAXIMUM TOLERABLE PROCEDURAL SIGMA"
    )

    print("-" * 60)

    for name, key in (
        ("Primary experiment", "registration_only"),
        ("Combined experiment", "combined_uncertainty"),
    ):
        if (
            key == "combined_uncertainty"
            and args.combined_trials == 0
        ):
            print(
                "\nCombined experiment: not run"
            )
            continue

        print_threshold(
            name,
            report["maximum_tolerable_sigma_mm"][key],
            report["maximum_conservative_sigma_mm"][key],
        )

    # --------------------------------------------------------
    # Save files
    # --------------------------------------------------------

    output_dir = args.output_dir

    json_path = save_json(
        report,
        output_dir / "targeting-results.json",
    )

    csv_path = save_csv(
        report["sigma_sweep"],
        output_dir / "targeting-sigma-sweep.csv",
    )

    print("\nOUTPUT FILES")
    print("-" * 60)

    print(f"JSON: {json_path}")
    print(f"CSV:  {csv_path}")

    print("\nQ4 experiment completed.")


if __name__ == "__main__":
    main()

