# Verification of Q1-Q3

[Back to the project](../README.md)

From the repository root, run:

```bash
python -m verification.run_all
```

The runner executes each completed question, then discovers all unit tests.
It stops with a failure status if a question or test fails. Q4 is deliberately
empty and is not included. All 63 tests passed in the saved run: 22 for Q1,
12 for Q2, 26 for Q3, and 3 cross-question checks.

## File guide

| File | Purpose |
| --- | --- |
| [README.md](README.md) | How to reproduce and locate the validation evidence |
| [run_all.py](run_all.py) | Runs Q1-Q3 and saves combined output and all test results |
| [test_integration.py](test_integration.py) | Three checks connecting registration, calibrated inputs, planning, and the Figure 6 metric |
| [full_reproduction.txt](full_reproduction.txt) | Actual commands, version information, numerical comparisons, and test output |
| [all_tests.txt](all_tests.txt) | Full list of test names and pass/fail results |
| `__init__.py` | Allows module execution and test discovery |

## Where each requirement is checked

| Requirement | Explanation and predictions | Runnable evidence |
| --- | --- | --- |
| Tip and axis calibration | [Q1](../q1_tool_calibration/README.md) | Q1 runner and calibration tests |
| Exact and noisy calibration | [Q1 experiment](../q1_tool_calibration/README.md#exact-and-noisy-calibration-experiments) | Q1 run log, JSON trial errors, and figures |
| Translation, rotation, and combined pose | [Q2](../q2_registration/README.md) | Q2 runner and registration tests |
| Intersecting, tangent, and missing paths | [Q3](../q3_biopsy_plan/README.md) | Q3 runner and planning tests |
| Figure 6 distance definition | [Q3 distance check](../q3_biopsy_plan/README.md#checking-the-meaning-of-figure-6) | Q3 runner, planning tests, and integration tests |
| Components used together | [Q3 combined calculation](../q3_biopsy_plan/README.md#checking-the-complete-calculation) | Q1 and Q3 runners; integration tests |

The explanations derive the expected answers from the geometry. The runners
compare full-precision computed values against those predictions. The tests
also cover nonzero tips, oblique axes, invalid marker layouts, and cases with
no forward insertion. The two small hand examples printed by Q1 are additional
demonstrations, not extra tests in the count of 63.

The JSON records seeds, units, pose counts, model coordinates, first-trial
results, summary statistics, and every noisy trial error. Small last-digit
differences can occur with a different numerical-library version. These
checks do not establish the unimplemented Q4 targeting limit; see the
[project limitations](../README.md#limitations).
