import os
import runpy

import pytest

path = os.path.join(os.path.dirname(__file__), "..", "tutorial", "heat_pump_exergy")
scripts = (os.path.join(path, f) for f in os.listdir(path) if f.endswith(".py") and f != "plots.py")


@pytest.mark.parametrize('script', scripts)
def test_tutorial_execution(script):
    try:
        runpy.run_path(script)
    except ModuleNotFoundError:
        pytest.skip("Test skipped due to missing dependency")