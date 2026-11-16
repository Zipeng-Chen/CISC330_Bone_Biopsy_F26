# Verification of Q1-Q4

[Back to the project](../README.md)

From the repository root, run:

```bash
python -m verification.run_all
```

This runs the Q1-Q3 examples and all Q1-Q4 unit tests, stopping if any command
fails. The saved run contains **75 passing tests**: 22 for Q1, 12 for Q2,
26 for Q3, 12 for Q4, and 3 cross-question checks. The runner refreshes Q1-Q3
outputs and the two verification logs below. It does not run the full Q4 sigma
sweep or overwrite Q4's saved simulation results.

To run only Q4's tests, or all tests without regenerating question outputs:

```bash
python -m unittest q4_targeting_error.test_targeting_error -v
python -m unittest discover -s . -p "test_*.py" -v
```

## File guide

| File | Purpose |
| --- | --- |
| [README.md](README.md) | Commands, hand predictions, and validation evidence |
| [run_all.py](run_all.py) | Runs Q1-Q3 examples and Q1-Q4 tests, and saves combined output |
| [test_integration.py](test_integration.py) | Three checks connecting registration, calibrated inputs, planning, and the Figure 6 metric |
| [full_reproduction.txt](full_reproduction.txt) | Actual commands, versions, numerical comparisons, and test output |
| [all_tests.txt](all_tests.txt) | Full list of test names and pass/fail results |
| [Q4 tests](../q4_targeting_error/test_targeting_error.py) | Twelve checks of geometry, noise propagation, statistics, and inputs |
| `__init__.py` | Allows module execution and test discovery |

## Where each requirement is checked

| Requirement | Explanation | Runnable evidence |
| --- | --- | --- |
| Tip and axis calibration | [Q1](../q1_tool_calibration/README.md) | Q1 runner and calibration tests |
| Exact and noisy calibration | [Q1 experiment](../q1_tool_calibration/README.md#exact-and-noisy-calibration-experiments) | Q1 log, JSON trial errors, and figures |
| Translation, rotation, and combined pose | [Q2](../q2_registration/README.md) | Q2 runner and registration tests |
| Intersecting, tangent, and missing paths | [Q3](../q3_biopsy_plan/README.md) | Q3 runner and planning tests |
| Figure 6 distance | [Q3 distance check](../q3_biopsy_plan/README.md#checking-the-meaning-of-figure-6) and Q4 checks below | Q3 runner, planning tests, Q4 tests, and integration tests |
| Components used together | [Q3 combined calculation](../q3_biopsy_plan/README.md#checking-the-complete-calculation) | Q1 and Q3 runners; integration tests |
| Fixed tracker noise and independent target displacement | [Q4](../q4_targeting_error/README.md) | Q1 noise tests; Q4 baseline, target-only, paired-trial, and sweep tests |
| Calibration uncertainty reaching the final error | Q4 checks below | Q4 controlled tip/axis and real noisy-calibration tests |
| Success rate, confidence interval, and feasible sigma | Q4 statistical checks below | Q4 summary and threshold tests; separate Monte Carlo sweep |

## Q4 hand calculations and test coverage

The [Q4 test file](../q4_targeting_error/test_targeting_error.py) uses small cases
whose answers we can predict before running the simulation. Fixed seeds make
the noisy cases repeatable. Controlled observations replace only the input
being isolated; the registration and distance calculations still run.

| Test method (prefix `test_`) | Expected result and reason |
| --- | --- |
| `ideal_scene_has_correct_units_and_zero_error` | Checks markers in mm, target `(0,-100,0)`, tracker coordinate sigma `0.3577198742584623` mm, and exact-scene error below `1e-10` mm. This catches incorrect units or depth. |
| `known_marker_shift_uses_perpendicular_not_endpoint_distance` | Shifts all markers 3 mm in X and places the target at `(0,-104,0)`. The axis is 3 mm from the target, while its 100 mm endpoint is 5 mm away. Targeting error must be 3 mm. |
| `known_marker_rotation_gives_hand_predicted_lateral_error` | Rotates all markers by 1 degree about the tip. The perpendicular miss is `100*sin(1 degree) = 1.745240644` mm. |
| `zero_target_displacement_keeps_instrument_noise_active` | At procedural sigma zero, targets remain fixed but targeting errors still vary because tracker noise remains active. |
| `target_only_noise_agrees_with_rayleigh_prediction` | Temporarily removes tracker noise for this test only. At target sigma 1 mm, predicted success is `98.8891%` and the predicted 95th-percentile error is `2.447746831` mm. |
| `paired_trials_share_ct_and_target_noise` | Supplies exact calibration to the paired trial function. Primary and combined errors must match within each pair, even with CT noise and target displacement. |
| `combined_trial_propagates_calibrated_tip_and_axis` | Supplies a tip shifted 3 mm, then an axis tilted 1 degree. With other noise removed, combined errors must be 3 mm and `1.745240644` mm; primary errors stay zero. |
| `actual_noisy_calibration_is_repeatable` | Runs the real calibration solver twice with the same seed. Results repeat, stay finite and nonnegative, and differ between primary and combined trials. |
| `success_boundary_and_wilson_interval` | Errors `[0,3,3.01,4]` mm give two successes out of four. Exactly 3 mm counts as success. The 95% Wilson interval is approximately `[0.150039,0.849961]`. |
| `no_feasible_sigma_differs_from_zero_passing_sigma` | All failed levels give `None`; only sigma zero passing gives `0.0`. Exactly 95% observed success passes the empirical rule but fails the lower-confidence-bound rule. An unrun combined experiment remains `None`. |
| `sigma_sweep_changes_target_noise_without_changing_marker_noise` | Corresponding trials at sigma zero and 1 mm have identical registered tips and axes but different displaced targets, using the sweep's shared seed. Only procedural uncertainty is swept. |
| `invalid_inputs_are_rejected` | Rejects negative/nonfinite sigma, invalid error arrays, zero trials, and zero insertion depth. |

### Geometry predictions

For the shifted needle, the line is `(3,0,0) + s*(0,-1,0)`. The target's closest
point on that line is `(3,-104,0)`, so the sideways gap is 3 mm. The endpoint
at 100 mm is `(3,-100,0)`, giving `sqrt(3^2 + 4^2) = 5` mm. This checks the
Figure 6 distance definition without using another call to the distance
function to generate the expected answer.

For the rotated needle, the target and tip are 100 mm apart. The perpendicular
gap is the opposite side of the right triangle, so it is `100*sin(1 degree)`.
The tests allow small floating-point roundoff.

### Target-only probability prediction

With a fixed true needle along Y, displacement along Y does not change the
distance to the axis. For independent X and Z displacements of sigma 1 mm,
the error is `sqrt(X^2 + Z^2)`. Its cumulative probability is
`P(error <= r) = 1 - exp(-r^2/2)`. Therefore:

```text
P(error <= 3 mm) = 1 - exp(-9/2) = 0.9888910035
95th percentile = sqrt(-2*ln(0.05)) = 2.447746831 mm
```

The test uses 4,000 trials, allowing five binomial standard errors around the
predicted success probability and 0.1 mm around the predicted percentile.
These tolerances account for sampling variation. This Rayleigh prediction
applies only to the isolated target-only case; it is not an assumed distribution
for the full tracker-plus-target simulation.

### Success statistics and threshold

For two successes out of four trials, `p=0.5`, `n=4`, and `z=1.95996398454`.
The Wilson interval can be calculated independently as:

```text
centre = (p + z^2/(2*n)) / (1 + z^2/n)
half-width = z*sqrt(p*(1-p)/n + z^2/(4*n^2)) / (1 + z^2/n)
interval = centre +/- half-width = [0.150039, 0.849961]
```

The threshold test distinguishes an observed success rate of at least 95%
from the stricter rule that the lower confidence bound also reach 95%. It
also distinguishes no feasible tested sigma from a passing value of zero.
This prevents a failed baseline from being reported as a usable solution.

## Reproducing the full Q4 experiment

Unit tests do not replace the full Monte Carlo experiment. To reproduce the
saved primary and combined sweeps in a separate location without replacing submitted Q4 outputs:

```bash
python -m q4_targeting_error.run --trials 10000 --combined-trials 500 --step 0.25 --seed 2026 --output-dir verification/q4_reproduction
```

This writes JSON and CSV under `verification/q4_reproduction/`, separately
from the verification logs. See the [Q4 report](../q4_targeting_error/README.md#results)
for the reported experiment and its interpretation. The new tests do not
hard-code its observed success rate: a particular simulation percentage is
not independent proof of correctness.

The full sweep was also reproduced separately on 2026-11-16, using 10,000
primary and 500 combined trials at each of the 13 sigma levels, with seed
2026. Every saved success count and threshold agreed; the largest numerical
difference in the summary statistics was `8.88e-16`. All 26 CSV rows agree
with the saved JSON, and both Q4 report tables agree at the displayed precision.
This reproduction is separate from the fast unit-test command above.

The paired tests exercise `trial_pair()` and `simulate_paired()`. They do not
establish that the command-line report compares matched primary and combined
samples: the current runner reports its primary sweep separately and retains
only the combined array from the optional paired run. Use the paired arrays
directly when analyzing trial-by-trial differences.

## Interpreting the saved verification

[all_tests.txt](all_tests.txt) records the 75 passing tests, including the 12 Q4
methods above. [full_reproduction.txt](full_reproduction.txt) also records Q1-Q3
hand predictions, computed values, commands, and Python/NumPy versions. The
two small hand examples printed by Q1 are demonstrations, not extra tests in
this count.

Small last-digit differences can occur with another numerical-library version.
Passing these checks supports the implementation under the stated assumptions.
It does not prove clinical safety or validate unmodeled effects such as needle
bending or correlated tracking errors; see the
[project limitations](../README.md#limitations).
