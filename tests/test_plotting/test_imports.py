import subprocess
import sys

import pytest

from tespy.plotting import plot_property_diagram
from tespy.plotting import plot_property_diagram_plotly
from tespy.plotting import plot_QT_diagram_plotly

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


def test_missing_coolplot_names_the_extra(heat_pump, monkeypatch):
    monkeypatch.setitem(sys.modules, "coolplot", None)
    with pytest.raises(ImportError, match=r"tespy\[plot\]"):
        plot_property_diagram(heat_pump, "1")


def test_missing_plotly_names_the_extra(heat_pump, monkeypatch):
    monkeypatch.setitem(sys.modules, "plotly", None)
    monkeypatch.setitem(sys.modules, "plotly.graph_objects", None)
    with pytest.raises(ImportError, match=r"tespy\[interactive\]"):
        plot_QT_diagram_plotly(heat_pump.get_comp("condenser"))
    with pytest.raises(ImportError, match=r"tespy\[interactive\]"):
        plot_property_diagram_plotly(heat_pump, "1")
