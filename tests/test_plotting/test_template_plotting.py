from unittest.mock import MagicMock
from unittest.mock import patch

import pytest
from conftest import build_heat_pump
from matplotlib.axes import Axes

from tespy.models.template import ModelTemplate


class HeatPumpModel(ModelTemplate):

    def _parameter_lookup(self):
        return {}

    def _create_network(self):
        self.nw = build_heat_pump()


@pytest.fixture(scope="module")
def model():
    return HeatPumpModel()


def test_property_diagram(model):
    assert isinstance(model.plot_property_diagram("1", "Ts"), Axes)


def test_property_diagram_plotly(model):
    figure = model.plot_property_diagram_plotly("1")
    assert figure.layout.meta["tespy"]["kind"] == "property_diagram"


def test_QT_diagram(model):
    assert isinstance(model.plot_QT_diagram("condenser"), Axes)


def test_QT_diagram_plotly(model):
    figure = model.plot_QT_diagram_plotly("condenser", logo=False)
    assert len(figure.layout.images) == 0


@pytest.mark.parametrize("method", [
    "plot_Ts_diagram_matplotlib", "plot_logph_diagram_matplotlib"
])
def test_old_cycle_plots_are_deprecated(model, method):
    points = {
        "a": {"p": 1.0, "h": 100.0, "T": 25.0, "s": 1.0, "vol": 0.5},
        "b": {"p": 5.0, "h": 200.0, "T": 50.0, "s": 2.0, "vol": 0.1},
    }
    with (
        patch("tespy.models.template.single_fluid", return_value=None),
        patch("tespy.models.template.get_plotting_data", return_value=({}, points)),
        pytest.warns(FutureWarning, match="0.14"),
    ):
        getattr(model, method)("1")


def test_old_QT_plot_is_deprecated(model):
    with pytest.warns(FutureWarning, match="plot_QT_diagram"):
        model.plot_QT_diagram_matplotlib("condenser")


def test_register_diagram_is_deprecated(model):
    with pytest.warns(FutureWarning, match="0.14"):
        model.register_diagram("R290", MagicMock())
