# Q4: Targeting error simulation

Simulate how CT marker-localization uncertainty and procedural target displacement affect the accuracy of a CT-guided bone biopsy. Determine the maximum tolerable target-displacement standard deviation while maintaining the required clinical targeting accuracy.

[Project overview and units](../README.md) | [Full verification](../verification/README.md)

## Run this question

From the repository root:

```bash
python -m q4_targeting_error.run --trials 10000 --combined-trials 500 --step 0.25 --seed 2026
```

The runner performs 10,000 primary trials and 500 combined trials for each procedural standard deviation from 0 to 3.0 mm, in increments of 0.25 mm. It prints the targeting-error statistics, evaluates the clinical success requirement, and saves the results in `outputs/`.

The primary experiment assumes exact tool calibration. The combined experiment also includes noisy Q1 calibration and is included in the saved results. To run only the primary experiment, use `--combined-trials 0`; this replaces the output files with primary-only results unless a different `--output-dir` is supplied.

## File guide

| File | Purpose |
| --- | --- |
| [README.md](README.md) | Simulation methodology, mathematical model, results, and interpretation |
| [simulation.py](simulation.py) | Ideal scene, Gaussian uncertainty, Monte Carlo trials, and statistical analysis |
| [run.py](run.py) | Runs the sigma sweep, summarizes clinical success, and exports results |
| [test_targeting_error.py](test_targeting_error.py) | Twelve checks of units, hand-calculated geometry, noise propagation, calibration inputs, and statistics |
| [__init__.py](__init__.py) | Makes Q4 importable and enables automatic test discovery |
| [outputs/targeting-results.json](outputs/targeting-results.json) | Simulation configuration, baseline statistics, sigma sweep, and clinical threshold results |
| [outputs/targeting-sigma-sweep.csv](outputs/targeting-sigma-sweep.csv) | Tabular targeting-error statistics for each procedural sigma |

## Geometry and definitions

As specified in Figure 6, assume the tool model frame and CT coordinate frame coincide. The ground-truth marker coordinates are:

```text
M1 = (20, 180, 0) mm
M2 = (80, 180, 0) mm
M3 = (20, 240, 0) mm

Tool tip = (0, 0, 0) mm
Needle axis = (0, -1, 0)
```

The marker coordinates are converted from centimetres to millimetres, following the professor's clarification.

For a 100 mm insertion depth, the planned target lies directly along the ideal needle axis:

```text
Target = Tip + 100 * Axis
       = (0, 0, 0) + 100 * (0, -1, 0)
       = (0, -100, 0) mm
```

Without localization error or target displacement, the needle trajectory passes directly through the target. The expected targeting error is therefore zero.

The simulation confirms this, producing an exact-scene targeting error of approximately `9.887e-14` mm. This is effectively zero, with the small difference caused by floating-point roundoff.

## Tracker localization uncertainty

The assignment specifies that 95% of three-dimensional marker-localization errors must have a magnitude below 1.0 mm.

Each marker observation is perturbed by an independent, zero-mean, three-dimensional Gaussian error:

```text
M_noisy = M_true + e

e ~ N(0, sigma_tracker^2 * I)
```

Here, `I` is the 3-by-3 identity matrix. The Gaussian standard deviation is applied independently to the X, Y, and Z coordinates.

The tracker standard deviation is calculated using the shared `sigma_for_radial_bound()` helper:

```text
Tracker radial bound:         1.0 mm
Required radial probability:  95%
Sigma per coordinate:         0.357719874 mm
```

This distinction is important: the 1.0 mm requirement describes the magnitude of the three-dimensional error, not the standard deviation of each coordinate.

The shared noise model is checked by Q1's `test_empirical_radial_probability_and_coordinate_independence`, using 100,000 samples with seed 41. It checks that the fraction within 1.0 mm is within 0.003 of 0.95, and also checks coordinate means and covariance. Q4's tests separately check that the simulation uses the expected tracker sigma.

## Compound targeting error

For each Monte Carlo trial, the simulation performs the following steps:

1. Begin with the ideal tool geometry and planned target.
2. Add independent Gaussian localization errors to the three CT marker positions.
3. Register the tool using the noisy CT markers to estimate the needle tip and axis.
4. Apply an independent three-dimensional Gaussian displacement to the planned target.
5. Calculate the perpendicular distance between the registered needle trajectory and the displaced target.
6. Record whether the targeting error satisfies the 3.0 mm clinical tolerance.

The procedural target displacement is modeled as:

```text
Displacement ~ N(0, sigma_proc^2 * I)

Displaced target = Planned target + Displacement
```

The tracker standard deviation remains fixed throughout the experiment. Only the procedural displacement standard deviation, `sigma_proc`, is varied.

Let `P` be the registered needle tip, `u` the registered unit needle direction, and `C` the displaced target. Define `w=C-P`.

The compound targeting error is:

```text
error = norm(w - dot(w,u)*u)
```

This is the shortest perpendicular distance from the displaced target to the registered needle axis, consistent with the assignment's definition.

The simulation disregards variations in insertion depth and tissue-core length, as instructed.

## Clinical requirement

The assignment requires that, for a 100 mm needle insertion depth, the targeting error remain within 3.0 mm with at least 95% probability.

A trial is successful when:

```text
targeting_error <= 3.0 mm
```

For each procedural sigma, the estimated success probability is:

```text
success_probability = successful_trials / total_trials
```

The clinical requirement is satisfied when:

```text
success_probability >= 0.95
```

The program also reports the mean, median, root-mean-square, maximum, and 95th-percentile targeting errors. A 95% Wilson confidence interval is calculated for the estimated success probability.

The confidence interval describes uncertainty in the Monte Carlo estimate. It is separate from the assignment's requirement that 95% of procedures satisfy the 3.0 mm targeting tolerance.

## Monte Carlo experiment

The primary experiment uses the following configuration:

| Parameter | Value |
| --- | ---: |
| Insertion depth | 100 mm |
| Targeting tolerance | 3.0 mm |
| Required success probability | 95% |
| Tracker sigma per coordinate | 0.357719874 mm |
| Procedural sigma range | 0.0–3.0 mm |
| Procedural sigma increment | 0.25 mm |
| Trials per sigma | 10,000 |
| Random seed | 2026 |
| Combined calibration trials per sigma | 500 |

The procedural sigma sweep contains 13 values, including the zero-displacement baseline.

At each value, the simulation repeats the registration and targeting-error calculation using randomly generated marker-localization errors and target displacements.

Trials are independent within each sigma level. The random generator restarts with the same seed at the next level, reusing the same underlying draws across levels. This helps compare the effect of changing sigma, but the 13 levels are not 13 independent experiments that can be pooled for a smaller confidence interval.

## Results

The following results were obtained from the saved simulation run.

| Procedural sigma (mm) | Mean error (mm) | 95th-percentile error (mm) | Success probability |
| ---: | ---: | ---: | ---: |
| 0.00 | 2.717 | 5.569 | 62.66% |
| 0.25 | 2.736 | 5.586 | 62.25% |
| 0.50 | 2.793 | 5.686 | 60.66% |
| 0.75 | 2.883 | 5.830 | 58.28% |
| 1.00 | 3.003 | 5.994 | 55.04% |
| 1.25 | 3.149 | 6.257 | 51.60% |
| 1.50 | 3.318 | 6.540 | 48.03% |
| 1.75 | 3.506 | 6.901 | 44.59% |
| 2.00 | 3.710 | 7.273 | 40.89% |
| 2.25 | 3.928 | 7.668 | 37.01% |
| 2.50 | 4.158 | 8.125 | 33.99% |
| 2.75 | 4.399 | 8.567 | 30.91% |
| 3.00 | 4.648 | 9.074 | 28.07% |

### Tracker-only baseline

When procedural sigma is zero, the target remains at its planned position. The remaining targeting error comes from CT marker-localization uncertainty propagated through rigid registration.

The baseline results were:

```text
Trials:                    10000
Mean targeting error:      2.717 mm
Median targeting error:    2.500 mm
RMS targeting error:       3.120 mm
95th-percentile error:     5.569 mm
Successful trials:         6266
Success probability:       62.66%
Required probability:      95.00%
```

The estimated 95% confidence interval for the baseline success probability was approximately 61.71% to 63.60%.

Although the mean targeting error is below 3.0 mm, this does not imply that the clinical requirement is satisfied. The requirement concerns the proportion of trials within the tolerance, not the average targeting error.

Small marker errors can produce a larger miss because the target is far from the marker cluster. Its centroid is `(40,200,0)` mm, while the target is `(0,-100,0)` mm: their separation is about 303 mm. Rotation uncertainty around the cluster therefore has a long distance over which to affect the needle trajectory. For scale, a 1-degree orientation error over 300 mm corresponds to about `300*sin(1 degree) = 5.24` mm of lateral displacement. This is an illustration of sensitivity, not a prediction that every trial has a 1-degree error. The 1 mm tracker bound applies at the markers, not at the target.

### Combined calibration experiment

The saved combined run uses fresh independent noise on 60 pivot and 60 spin observations in each trial. Q1 estimates the tip and axis from those observations, and Q2 registers those estimated values using noisy CT markers. The target displacement is independent of both sets of marker noise.

| Procedural sigma (mm) | Combined 95th-percentile error (mm) | Combined success probability |
| ---: | ---: | ---: |
| 0.00 | 5.838 | 61.40% |
| 0.25 | 5.955 | 58.40% |
| 0.50 | 6.094 | 56.00% |
| 0.75 | 6.106 | 55.00% |
| 1.00 | 6.281 | 50.60% |
| 1.25 | 6.466 | 48.00% |
| 1.50 | 6.721 | 45.40% |
| 1.75 | 7.181 | 41.60% |
| 2.00 | 7.610 | 38.20% |
| 2.25 | 7.904 | 35.20% |
| 2.50 | 8.505 | 32.40% |
| 2.75 | 8.893 | 29.20% |
| 3.00 | 9.252 | 27.40% |

At sigma zero, 307 of 500 combined trials succeed. The mean error is 2.820 mm and the 95% Wilson interval for success is 57.06% to 65.56%. This is much wider than the primary interval because the combined sample is smaller. Neither baseline comes close to the required 95% success rate.

The reported 61.40% and 62.66% rates are separate estimates, not a matched estimate of the effect of calibration. Although `simulate_paired()` computes matched primary and combined errors internally, the command-line runner retains only its combined array and reports a separate 10,000-trial primary run. The difference between the reported rates cannot be attributed entirely to calibration. A precise comparison would retain both paired arrays and use a larger combined sample. These results support the same failure of the clinical requirement, but do not measure a precise calibration penalty.

### Effect of increasing procedural uncertainty

Increasing procedural sigma generally increases the targeting error and decreases the probability of remaining within the clinical tolerance.

At sigma = 1.0 mm, the success probability was 55.04%. At sigma = 2.0 mm, it decreased to 40.89%. At sigma = 3.0 mm, it reached 28.07%.

The 95th-percentile targeting error increased from 5.569 mm at the baseline to 9.074 mm at sigma = 3.0 mm.

These results show that procedural target displacement introduces additional uncertainty beyond the existing registration error.

## Maximum tolerable procedural sigma

The assignment asks for the maximum procedural displacement standard deviation that can be tolerated while maintaining at least 95% targeting success.

The simulation shows that the clinical requirement is not satisfied even when procedural sigma is zero.

```text
Baseline success probability: 62.66%
Required success probability: 95.00%

Maximum tolerable procedural sigma:
None feasible among the tested values.
```

The output JSON records this result as:

```json
"maximum_tolerable_sigma_mm": {
  "registration_only": null,
  "combined_uncertainty": null
}
```

For the primary experiment, `null` indicates that no tested procedural sigma satisfies the clinical requirement.

The combined experiment also ran and failed at every tested sigma, so its `null` value has the same meaning. If `--combined-trials 0` is used instead, a combined `null` means that the combined experiment was not evaluated; check the configuration to distinguish the two cases.

Because the baseline success probability is already below 95%, reducing procedural displacement alone cannot resolve the observed baseline targeting inaccuracy.

## Validation and limits

The ideal-scene calculation produced a targeting error of approximately `9.887e-14` mm, confirming that the noise-free needle trajectory passes through the planned target.

Run the 12 Q4-specific tests from the repository root:

```bash
python -m unittest q4_targeting_error.test_targeting_error -v
```

The tests include two hand-calculated cases: a 3 mm sideways shift must give a 3 mm perpendicular error even when the endpoint distance is 5 mm; a 1-degree tilt at 100 mm depth must give `100*sin(1 degree) = 1.745240644` mm. They also check that tracker noise remains active at zero procedural sigma, calibrated tip and axis errors reach the combined result, and success counts and confidence intervals are correct.

The [verification explanation](../verification/README.md#q4-hand-calculations-and-test-coverage) gives the predictions and purpose of all 12 tests, including an independent target-only probability check. [Saved test results](../verification/all_tests.txt) record their outcomes. These tests are also included in `python -m verification.run_all`; they do not overwrite the Q4 simulation outputs. Q2's tests cover known rigid transformations, and Q1's tests check the shared tracker noise distribution.

The simulation assumes rigid tool markers, exact ground-truth geometry, independent zero-mean Gaussian errors, and a fixed 100 mm insertion depth. The primary experiment also assumes exact tool calibration, so noisy calibration uncertainty from Q1 is not included in the reported primary results.

The saved combined experiment includes calibration uncertainty with 500 trials per sigma. Its smaller sample supports a broad comparison with the requirement, rather than a precise estimate of how much calibration adds to the primary error.

The model does not simulate needle bending, tissue deformation beyond the randomized target displacement, correlated tracking errors, or changes in insertion depth and tissue-core length. Its results therefore describe the specified mathematical uncertainty model rather than a complete clinical procedure.

The [JSON results](outputs/targeting-results.json) contain the full configuration and statistical summaries. The [CSV results](outputs/targeting-sigma-sweep.csv) provide the procedural sigma sweep in tabular form.

The simulation demonstrates that, under the specified marker geometry and uncertainty assumptions, tracker-localization error alone is sufficient to prevent the modeled procedure from meeting the 95% clinical targeting requirement. No tested procedural displacement standard deviation satisfies the requirement.

See the [shared limitations and references](../README.md#limitations) for additional project-wide assumptions.
