# CT-guided Bone Biopsy Planning

CISC/CMPE 330 Assignment 2.

**Start with [the report](docs/report.md).** It explains the methods, works out
the expected answers by hand, compares them with the program output, and
discusses what the results mean. The file guide below shows where each part
is implemented and tested.

## Work included

| Component | Status |
| --- | --- |
| Tool calibration | Tip from pivoting and axis from spinning; exact checks and 500 noisy trials |
| Registration | Model-to-CT pose; translation, rotation, and combined tests |
| Biopsy planning | Targeting error, core length, and insertion depth; hit, tangent, and miss tests |
| Targeting-error simulation | Separate team component; the target-displacement sweep is not included here |

The workflow is **calibration -> registration -> biopsy planning**. Figure 2
positions start in cm and are converted to mm. All later positions and lengths
use mm; directions are unit vectors. Detailed conventions are in report Section 1.

## Run the code

Use Python 3.10 or newer. Run these commands from this repository's root:

```console
python -m pip install -r requirements.txt
python -m examples.run_required_cases
python -m examples.run_calibration
python -m unittest discover -s tests -v
```

The first runner prints the registration/planning predictions and computed
answers. The second prints the small hand-calculated calibration checks,
exact calibration, and the 500-trial noisy experiment. The test suite contains
63 tests. The [saved console output](docs/validation-output.txt) records a
successful run of both examples and the complete test suite.

To save the calibration results again:

```console
python -m examples.run_calibration --trials 500 --poses 60 --seed 2026 --output docs/calibration-results.json
```

The calibration figures are optional outputs. To regenerate them:

```console
python -m pip install -r requirements-plots.txt
python -m examples.run_calibration --trials 500 --poses 60 --seed 2026 --output docs/calibration-results.json --plot-dir docs
```

NumPy is the only dependency needed for the calculations and tests. Matplotlib
is needed only for plotting.

## File guide

### Report and recorded results

| File | What it contains |
| --- | --- |
| [README.md](README.md) | Starting point, run commands, project scope, and this file guide |
| [docs/report.md](docs/report.md) | Main explanation: units, calibration, hand calculations, registration, planning, comparisons, and limitations |
| [docs/validation-output.txt](docs/validation-output.txt) | Actual console output from the two example runners and the 63-test suite |
| [docs/calibration-results.json](docs/calibration-results.json) | Experiment settings, hand-check results, exact/noisy estimates, summary statistics, and all 500 trial errors |
| [docs/calibration_motions.png](docs/calibration_motions.png) | Illustration of the pivot and spin constraints used by calibration |
| [docs/calibration_errors.png](docs/calibration_errors.png) | Histograms of the recorded tip and axis errors |
| [docs/registration_frames.png](docs/registration_frames.png) | Tool geometry before and after the worked rotation/translation |
| [docs/biopsy_trajectories.png](docs/biopsy_trajectories.png) | The intersecting, tangent, and missing unit-sphere examples |

### Implementation

| File | What it does |
| --- | --- |
| [bone_biopsy/models.py](bone_biopsy/models.py) | Defines the shared Figure 2 model and converts cm to mm |
| [bone_biopsy/_validation.py](bone_biopsy/_validation.py) | Checks shapes, finite coordinates, nonnegative inputs, and nonzero directions |
| [bone_biopsy/registration.py](bone_biopsy/registration.py) | Fits a rigid pose from matching markers and transforms a tip or direction |
| [bone_biopsy/planning.py](bone_biopsy/planning.py) | Calculates point-to-axis distance, sphere intersections, core length, and insertion depth |
| [bone_biopsy/calibration.py](bone_biopsy/calibration.py) | Solves pivot-tip and spin-axis calibration using the existing registration function |
| [bone_biopsy/calibration_simulation.py](bone_biopsy/calibration_simulation.py) | Generates the constrained motions and simulated marker observations |
| [bone_biopsy/geometry.py](bone_biopsy/geometry.py) | Provides axis-angle rotation, plane fitting, and angular-error calculations |
| [bone_biopsy/noise.py](bone_biopsy/noise.py) | Converts the radial error bound to coordinate sigma and adds Gaussian marker noise |
| [bone_biopsy/__init__.py](bone_biopsy/__init__.py) | Exposes the functions used by the examples and other team components |

### Examples and tests

| File | What it checks or produces |
| --- | --- |
| [examples/run_required_cases.py](examples/run_required_cases.py) | Runs the six required registration/planning cases, the combined example, and the Figure 6 distance comparison |
| [examples/run_calibration.py](examples/run_calibration.py) | Runs the hand-calculated pivot/spin checks, exact calibration, noisy trials, and calibration-to-planning check |
| [examples/plot_calibration.py](examples/plot_calibration.py) | Draws the motion and error plots from the calibration data; it does not implement another solver |
| [tests/test_registration.py](tests/test_registration.py) | Checks known poses, nonzero tips, reflections, model units, and invalid marker layouts |
| [tests/test_planning.py](tests/test_planning.py) | Checks the three required paths, forward-only insertion, near tangency, and invalid geometry |
| [tests/test_calibration.py](tests/test_calibration.py) | Checks the motion constraints, hand fixture, nonzero/oblique model, noise distribution, degenerate inputs, and integration |
| [tests/test_integration.py](tests/test_integration.py) | Checks registration followed by planning, including the displaced-target example |
| [examples/__init__.py](examples/__init__.py), [tests/__init__.py](tests/__init__.py) | Mark the folders as Python packages; they contain no separate algorithm |

### Setup files

| File | Purpose |
| --- | --- |
| [requirements.txt](requirements.txt) | NumPy dependency for calculations and tests |
| [requirements-plots.txt](requirements-plots.txt) | Adds Matplotlib for optional figures |
| [.gitignore](.gitignore) | Keeps Python caches and local environments out of Git |

The example runners make the worked cases easy to inspect. The unit tests
check additional behavior automatically. Both call the same implementation
modules; neither contains a second copy of a calibration, registration, or
planning algorithm.

## Connecting to the other team component

`calibrate_tool(model_markers, pivot_observations, spin_observations)` returns
`tip` and `axis` in the model frame. Pass them into
`register_tool(model_markers, ct_markers, tip, axis)` to obtain CT-frame values,
then pass those to `plan_biopsy`.

Keep the marker order and mm units consistent. The noise helper
`sigma_for_radial_bound()` returns the per-coordinate sigma for the assignment's
1 mm radial 95% bound. The target-displacement simulation must still be added
and evaluated separately; the calibration results do not establish its limit.

## Submission contents

Include this README, the report and its supporting results/figures, the source
modules, the example runners, the tests, and the dependency files listed above.
The report is the main written explanation; console output and test results
support it.

When preparing a ZIP, leave out `.git`, `.venv`, `__pycache__`, notebook
checkpoints, temporary files, and local backup copies. Keep one current report
and one implementation of each algorithm.
