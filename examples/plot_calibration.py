"""Optional calibration figures; computational calibration only requires NumPy"""

from pathlib import Path

import numpy as np


def save_calibration_figures(model, pivot, spin, report, output_dir):
    """Plot the constrained motions and the recorded Monte Carlo error distributions"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError("Install requirements-plots.txt to generate figures.") from error
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    colors = ["#176b87", "#c75b12", "#26794e"]
    fig = plt.figure(figsize=(12, 5.5), layout="constrained")
    for index, (motion, title) in enumerate(
            ((pivot, "Tip calibration: multi-axis pivot"), (spin, "Axis calibration: fixed-sleeve spin")), 1):
        ax = fig.add_subplot(1, 2, index, projection="3d")
        for marker, color in enumerate(colors):
            cloud = motion.observations[:, marker]
            ax.scatter(*cloud.T, s=9, color=color, label=f"M{marker + 1}")
            if index == 2:
                circle = np.vstack([cloud, cloud[0]])
                ax.plot(*circle.T, color=color, lw=1, alpha=.5)
        # Draw a few shaft poses to make the different motion constraints visible
        for pose in np.linspace(0, len(motion.rotations) - 1, 6, dtype=int):
            tip = motion.rotations[pose] @ model.tip + motion.translations[pose]
            axis = motion.rotations[pose] @ model.axis
            ax.plot(*np.array([tip, tip - 180 * axis]).T, color="#65707d", alpha=.4, lw=1)
        ax.scatter(*tip, color="black", s=30, label="Fixed tip")
        cloud = np.vstack([motion.observations.reshape(-1, 3), tip])
        center = (cloud.min(axis=0) + cloud.max(axis=0)) / 2
        half = np.max(np.ptp(cloud, axis=0)) * .56
        ax.set_xlim(center[0] - half, center[0] + half)
        ax.set_ylim(center[1] - half, center[1] + half)
        ax.set_zlim(center[2] - half, center[2] + half)
        ax.set_box_aspect((1, 1, 1))
        ax.set(xlabel="Tracker x (mm)", ylabel="Tracker y (mm)", zlabel="Tracker z (mm)", title=title)
        ax.legend(fontsize=8, loc="upper left")
    fig.savefig(destination / "calibration_motions.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    cases = [
        ("trial_tip_errors_mm", "tip_error_mm", "Tool-tip position error (mm)", colors[0]),
        ("trial_axis_errors_degrees", "axis_error_degrees", "Tool-axis angular error (degrees)", colors[1]),
    ]
    for ax, (samples, summary, label, color) in zip(axes, cases):
        # Plot the actual recorded trials and their empirical 95th percentile.
        ax.hist(report[samples], bins=25, color=color, edgecolor="white", alpha=.85)
        p95 = report[summary]["p95"]
        ax.axvline(p95, color="#222222", linestyle="--", label=f"95th percentile = {p95:.4f}")
        ax.set(xlabel=label, ylabel="Calibration trials")
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(fontsize=9)
    config = report["configuration"]
    fig.suptitle(f"{config['trials']} noisy calibrations: {config['pivot_poses']} pivot + "
                 f"{config['spin_poses']} spin poses; marker sigma = "
                 f"{config['marker_sigma_per_coordinate_mm']:.4f} mm per coordinate")
    fig.savefig(destination / "calibration_errors.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
