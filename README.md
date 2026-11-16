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

Q4 contains only an empty `.gitkeep` file so Git preserves the folder.
There is no Q4 implementation, explanation, or result yet.

## Run and reproduce

From the repository root, with Python 3.10 or later:

```bash
python -m pip install -r requirements.txt
python -m verification.run_all
```

This runs Q1-Q3, saves their outputs, and runs all 63 unit tests. It explicitly
skips the empty Q4 folder. The recorded run uses Python 3.12.14 and NumPy 2.3.5;
other numerical-library versions may produce small last-digit differences.
The saved run can be read in
[verification/full_reproduction.txt](verification/full_reproduction.txt), with
all test names in [verification/all_tests.txt](verification/all_tests.txt).

To run one question:

```bash
python -m q1_tool_calibration.run
python -m q2_registration.run
python -m q3_biopsy_plan.run
```

Each command prints the hand predictions and computed answers, runs that
question's tests, and refreshes its `outputs/run_log.txt` and
`outputs/test_report.txt`. Q1 also saves all 500 noisy trials to JSON.
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

## How the parts connect

Calibration finds the tip and axis relative to the tool markers. Registration
finds where that calibrated tool is in a new CT image. Planning uses the
registered tip and axis to calculate targeting error, tissue-core length, and
insertion depth.

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
for the assignment's 1 mm radial 95% bound. Q4 can reuse this helper and the
calibration, registration, and planning functions. Its target-displacement
simulation and final threshold analysis still need to be added.

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

The 95th-percentile calibration errors do not prove the separate 100 mm / 3 mm
compound targeting requirement. That requires the remaining simulation with
marker-localization error and displaced targets. No maximum tolerable target
displacement is claimed here.

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

