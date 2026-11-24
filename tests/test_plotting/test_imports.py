import subprocess
import sys

PLOTTING_LIBRARIES = ("matplotlib", "plotly", "coolplot", "fluprodia")


def test_core_packages_do_not_import_plotting_libraries():
    code = (
        "import sys\n"
        "import tespy, tespy.networks, tespy.components, tespy.connections\n"
        "import tespy.models, tespy.plotting\n"
        f"print([m for m in {PLOTTING_LIBRARIES!r} if m in sys.modules])\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout.strip() == "[]"
