import os
import runpy

import pytest

path = os.path.join(os.path.dirname(__file__), "..", "tutorial", "basics")
scripts = (os.path.join(path, f) for f in os.listdir(path) if f.endswith(".py"))


@pytest.mark.parametrize('script', scripts)
def test_tutorial_execution(script):
    runpy.run_path(script)