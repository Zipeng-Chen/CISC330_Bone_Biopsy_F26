# CT-guided Bone Biopsy Planning

CISC/CMPE 330 Assignment 2

Each completed question has its own explanation, code, tests, and saved output.
Read the question's README for the hand calculations and interpretation, then
use its run log to compare the predictions with the computed results.

## Questions

| Question | Folder and explanation | Status |
| --- | --- | --- |
| Q1: Tool calibration | [q1_tool_calibration](q1_tool_calibration/README.md) | Implemented and checked with exact and noisy observations |
| Q2: Registration | [q2_registration](q2_registration/README.md) | Implemented; translation, rotation, and combined cases checked |
| Q3: Biopsy planning | [q3_biopsy_plan](q3_biopsy_plan/README.md) | Implemented; hit, tangent, miss, and combined cases checked |
| Q4: Targeting-error simulation | [q4_targeting_error](q4_targeting_error/) | Implemented and tested with Gaussian tracker noise, target displacement, and 10,000 Monte Carlo trials per sigma |


## Run and reproduce

From the repository root, with Python 3.10 or later:

```bash
python -m pip install -r requirements.txt
python -m verification.run_all
```

This runs Q1-Q3, saves their outputs, and runs all 75 unit tests, including 12 Q4-specific tests. The full Q4 Monte Carlo sweep is run separately using the command below; verification does not overwrite Q4's saved results. The recorded run uses Python 3.12.14 and NumPy 2.3.5;
other numerical-library versions may produce small last-digit differences.
The saved run can be read in
[verification/full_reproduction.txt](verification/full_reproduction.txt), with
all test names in [verification/all_tests.txt](verification/all_tests.txt).

To run one question:

```bash
python -m q1_tool_calibration.run
python -m q2_registration.run
python -m q3_biopsy_plan.run
python -m q4_targeting_error.run --trials 10000 --combined-trials 500 --step 0.25 --seed 2026
```

Q1-Q3 commands print the hand predictions and computed answers, run that
question's tests, and refresh its `outputs/run_log.txt` and
`outputs/test_report.txt`. Q1 also saves all 500 noisy trials to JSON. 
Q4 runs 10,000 primary trials and 500 combined-calibration trials per procedural sigma, from 0 to 3.0 mm in increments of 0.25 mm,
and saves its results to q4_targeting_error/outputs/targeting-results.json 
and q4_targeting_error/outputs/targeting-sigma-sweep.csv.

To run only the tests:

```bash
python -m unittest discover -s . -p "test_*.py" -v
```

NumPy is the only dependency for the calculations and tests. Matplotlib is
optional for regenerating Q1 figures; see the
[Q1 instructions](q1_tool_calibration/README.md#run-this-question).
The four saved figures are included for reading the explanations.

## File guide

The question READMEs explain every file in their folders. Shared files are:

| File or folder | Purpose |
| --- | --- |
| [README.md](README.md) | Question map, run commands, common conventions, limitations, and references |
| [utils/models.py](utils/models.py) | Shared Figure 2 tool model, converted from cm to mm |
| [utils/geometry.py](utils/geometry.py) | Rotation, plane fitting, and angular error calculations |
| [utils/noise.py](utils/noise.py) | Gaussian marker noise and the 3D radial probability calculation |
| [utils/_validation.py](utils/_validation.py) | Shared checks for coordinates, directions, and scalar inputs |
| [utils/reporting.py](utils/reporting.py) | Saves console output and each question's test results |
| [utils/__init__.py](utils/__init__.py) | Makes the shared helpers importable |
| [bone_biopsy/__init__.py](bone_biopsy/__init__.py) | Preserves the team's existing imports; points to the implementations in Q1-Q3 and utils |
| [verification](verification/README.md) | Full reproduction command, cross-question tests, and combined logs |
| [requirements.txt](requirements.txt) | NumPy dependency for calculations and tests |
| [requirements-plots.txt](requirements-plots.txt) | Optional Matplotlib dependency for calibration figures |
| [.gitignore](.gitignore) | Excludes local environments, caches, and compiled Python files |
| [q4_targeting_error/simulation.py](q4_targeting_error/simulation.py) | Monte Carlo simulation of tracker localization error, target displacement, and compound targeting error |
| [q4_targeting_error/run.py](q4_targeting_error/run.py) | Runs the procedural sigma sweep, summarizes targeting success, and saves results |
| [q4_targeting_error/outputs/](q4_targeting_error/outputs/) | Generated JSON and CSV results for the targeting-error simulation |
## How the parts connect

Calibration finds the tip and axis relative to the tool markers. Registration
finds where that calibrated tool is in a new CT image. Planning uses the
registered tip and axis to calculate targeting error, tissue-core length, and
insertion depth. Finally, Q4 simulates how marker-localization uncertainty and target displacement affect targeting accuracy.

The existing public imports still work:

```python
from bone_biopsy import calibrate_tool, register_tool, plan_biopsy
```

`calibrate_tool(model_markers, pivot_observations, spin_observations)` returns
`tip` and `axis` in the model frame. Pass them into
`register_tool(model_markers, ct_markers, tip, axis)` to obtain CT-frame values,
then pass those to `plan_biopsy` with the tumor center and radius in mm.
Existing submodule imports such as `bone_biopsy.registration` also work.
There is one implementation of each algorithm, located in its question folder.

The noise helper `sigma_for_radial_bound()` returns the per-coordinate sigma
for the assignment's 1 mm radial 95% bound.  Q4 reuses this helper and the 
registration functions to simulate compound targeting error from Gaussian marker-localization noise and target displacement. 
It performs 10,000 Monte Carlo trials per procedural sigma, evaluates the 3 mm targeting tolerance at a 100 mm insertion depth, and determines whether the 95% clinical success requirement is satisfied.

## Units and coordinate frames

The coordinates in assignment Figure 2 are in cm, as clarified by the
professor. The code uses mm, so each position is multiplied by 10 once.

| Quantity | Figure 2 | Value used by the code |
| --- | --- | --- |
| M1 | (2,18,0) cm | (20,180,0) mm |
| M2 | (8,18,0) cm | (80,180,0) mm |
| M3 | (2,24,0) cm | (20,240,0) mm |
| Tip | (0,0,0) cm | (0,0,0) mm |
| Insertion axis | (0,-1,0) | (0,-1,0), dimensionless |

A direction has no length unit, so it is not multiplied by 10. Values already
given in mm, such as the 100 mm target depth and the 1 mm marker-error bound,
also stay unchanged. `assignment_tool_model()` applies the cm-to-mm conversion.

There are three frames: the model frame attached to the tool, the optical
tracker frame used during calibration, and the CT frame used during planning.
Every marker array uses the order M1, M2, M3. Matching the wrong rows would
describe a different correspondence problem.

For a point, the model-to-destination transformation is `x = R p + t`. For a
direction, it is `u = R v`. Translation moves a position but cannot change the
direction of a line. The arrays store points as rows, so the equivalent code
for several points is `points @ R.T + t`.

For Q4, the tool model and CT coordinate frames are assumed to coincide, as specified in the assignment. The ideal tool tip is at (0, 0, 0) mm, with the needle axis pointing along (0, -1, 0). At an insertion depth of 100 mm, the planned target is therefore located at (0, -100, 0) mm. Gaussian marker-localization noise is applied in the CT frame, while independent target displacement is added to the planned target position. All targeting errors and displacement values are expressed in mm.

## Limitations

The exact cases check whether the geometry and frame transformations are
implemented correctly. The noisy calibration experiment shows how the stated
marker noise affects one chosen calibration motion design. The motions stay
fixed across the 500 trials; the measurement noise changes independently.

The experiments assume rigid markers, exact model marker positions, independent
Gaussian measurement errors, and correctly constrained pivoting and spinning.
They do not include tip slip, sleeve movement, marker mounting error,
correlated tracking errors, or needle bending. These effects could make the
real errors different from the simulated ones.

The planning model contains a straight needle axis and a sphere. It does not
check nearby anatomy, whether the full needle diameter fits through an access
path, or whether a particular instrument can collect the calculated core.
A geometrically intersecting line is therefore not enough to establish that
a procedure is safe. A small marker-fit or calibration residual also cannot
establish safety by itself.

The Q4 simulation evaluates the compound targeting requirement using independent Gaussian marker-localization errors and target displacement. With 10,000 primary trials per procedural sigma, the baseline targeting success probability was 62.66%, below the required 95%. The saved combined experiment includes noisy Q1 calibration and gives 61.40% baseline success over 500 trials. No tested procedural sigma satisfies the requirement in either experiment. The primary simulation assumes exact calibration; the combined sample is smaller, and its reported difference from the primary run is not a paired estimate of the calibration effect. These findings are limited to the assumed geometry and error distributions; see the Q4 report for confidence intervals and interpretation.

## Submission contents

The four question READMEs form the written report: each explains the method,
hand calculations or statistical predictions, results, and limitations. Keep
the code, tests, question outputs, figures, discussion transcripts, dependency
files, and verification logs together so their relative links work.

Exclude `.git/`, `__pycache__/`, `.venv/`, `.pytest_cache/`, compiled Python
files, and any temporary reproduction folder from an uploaded archive. These
are local working files, not assignment evidence. The `bone_biopsy` folder is
an import interface used by Q4; it is required even though the implementations
live in the question folders.

## References

- *CT-guided Bone Biopsy Planning*, CISC/CMPE 330 Assignment 2: pages 1-2 for
  calibration; pages 3-4 for registration and planning; Figure 6 for targeting
  error; and the grading rubric for the required validation cases.
- Gabor Fichtinger, *Calibration* (2026), slides 5-7 for the least-squares
  pivot equations and slide 11 for spin-axis calibration from marker planes.
- *MathPrimer-Part3*: point/direction transformations, rigid registration,
  and the distinction between marker fit and target accuracy.
- Professor Fichtinger's feedback: Figure 2 uses cm; insertion traverses the
  full available core; and targeting error should be read from Figure 6.
- Wolfram MathWorld, *Kabsch Algorithm*: rigid-body registration using corresponding 3D marker positions, including rotation estimation and least-squares alignment.
  https://mathworld.wolfram.com/KabschAlgorithm.html
- National Institute of Standards and Technology (NIST), *Chi-Square Cumulative Distribution Function*: the chi-square distribution used to convert the tracker's 95% radial localization bound into a per-coordinate Gaussian standard deviation.
  https://itl.nist.gov/div898/software/dataplot/refman2/ch8/chscdf.pdf

