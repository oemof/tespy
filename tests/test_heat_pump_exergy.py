import os
import runpy
import sys

import pytest

path = os.path.join(
    os.path.dirname(__file__), "..", "tutorial", "heat_pump_exergy"
)
scripts = (
    os.path.join(path, f)
    for f in os.listdir(path)
    if f.endswith(".py") and f != "plots.py"
)


class add_path():
    def __init__(self, path):
        self.path = path

    def __enter__(self):
        sys.path.insert(0, self.path)

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            sys.path.remove(self.path)
        except ValueError:
            pass


@pytest.mark.parametrize('script', scripts)
def test_tutorial_execution(script):
    try:
        with add_path(path):
            runpy.run_path(script)
    except ModuleNotFoundError:
        pytest.skip("Test skipped due to missing dependency.")
