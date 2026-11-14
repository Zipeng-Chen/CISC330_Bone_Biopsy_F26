"""Save real console output and the tests for each question."""

from contextlib import redirect_stdout
from datetime import date
from io import StringIO
from pathlib import Path
import sys
import unittest


class _Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for stream in self.streams:
            stream.write(text)
        return len(text)

    def flush(self):
        for stream in self.streams:
            stream.flush()


def run_with_log(action, output_path, command):
    """Show output while saving the same text beside the question."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="\n") as log:
        with redirect_stdout(_Tee(sys.stdout, log)):
            print(f"Run date: {date.today().isoformat()}")
            print(f"Command: {command}\n")
            action()


def run_tests(module_name, output_path):
    """Keep the test names and results, rather than only a pass count."""
    stream = StringIO()
    suite = unittest.defaultTestLoader.loadTestsFromName(module_name)
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    Path(output_path).write_text(stream.getvalue(), encoding="utf-8")
    if not result.wasSuccessful():
        print(stream.getvalue())
        raise SystemExit(1)
    print(f"PASS: {result.testsRun} unit tests; details in outputs/test_report.txt")
