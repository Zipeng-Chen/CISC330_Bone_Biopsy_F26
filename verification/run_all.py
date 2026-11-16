"""Reproduce the completed questions and save all verification output."""

from datetime import date
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = Path(__file__).resolve().parent


def main():
    # Use the active Python and repository root, regardless of the caller's folder.
    with (OUTPUTS / "full_reproduction.txt").open("w", encoding="utf-8", newline="\n") as log:
        def record(text):
            print(text, end="", flush=True)
            log.write(text)
            log.flush()

        import numpy as np
        record(f"Run date: {date.today().isoformat()}\n")
        record(f"Python {platform.python_version()}; NumPy {np.__version__}\n")
        record("Q1-Q3 examples and Q1-Q4 unit tests are included.\n")
        record("Run the full Q4 sigma sweep separately; its saved results are not overwritten here.\n")
        for module in ("q1_tool_calibration.run", "q2_registration.run", "q3_biopsy_plan.run"):
            record(f"\n$ python -m {module}\n")
            process = subprocess.run(
                [sys.executable, "-B", "-m", module], cwd=ROOT,
                text=True, encoding="utf-8", errors="replace",
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            )
            record(process.stdout)
            if process.returncode:
                raise SystemExit(process.returncode)

        # Discover the question tests plus the checks connecting the components.
        record('\n$ python -m unittest discover -s . -p "test_*.py" -v\n')
        process = subprocess.run(
            [sys.executable, "-B", "-m", "unittest", "discover", "-s", ".", "-p", "test_*.py", "-v"],
            cwd=ROOT, text=True, encoding="utf-8", errors="replace",
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        (OUTPUTS / "all_tests.txt").write_text(process.stdout, encoding="utf-8")
        record(process.stdout)
        if process.returncode:
            raise SystemExit(process.returncode)
        record("\nPASS: Q1-Q3 examples and all Q1-Q4 unit tests completed.\n")


if __name__ == "__main__":
    main()
