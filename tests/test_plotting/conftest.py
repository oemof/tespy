import matplotlib
import pytest

matplotlib.use("Agg")

from matplotlib import pyplot as plt  # noqa: E402

from tespy.components import Compressor  # noqa: E402
from tespy.components import CycleCloser  # noqa: E402
from tespy.components import HeatExchanger  # noqa: E402
from tespy.components import MovingBoundaryHeatExchanger  # noqa: E402
from tespy.components import Sink  # noqa: E402
from tespy.components import Source  # noqa: E402
from tespy.components import Valve  # noqa: E402
from tespy.connections import Connection  # noqa: E402
from tespy.networks import Network  # noqa: E402
from tespy.plotting import settings  # noqa: E402


def build_heat_pump(T_evap=0, T_cond=50):
    """R290 heat pump with a moving boundary condenser and a plain
    evaporator, brine and water as heat source and sink."""
    nw = Network(iterinfo=False)
    nw.units.set_defaults(
        temperature="degC", pressure="bar", enthalpy="kJ/kg",
        entropy="kJ/kg/K", heat="kW", power="kW"
    )
    cc = CycleCloser("cycle closer")
    compressor = Compressor("compressor")
    condenser = MovingBoundaryHeatExchanger("condenser")
    valve = Valve("valve")
    evaporator = HeatExchanger("evaporator")
    water_in, water_out = Source("water in"), Sink("water out")
    brine_in, brine_out = Source("brine in"), Sink("brine out")

    c1 = Connection(cc, "out1", compressor, "in1", label="1")
    c2 = Connection(compressor, "out1", condenser, "in1", label="2")
    c3 = Connection(condenser, "out1", valve, "in1", label="3")
    c4 = Connection(valve, "out1", evaporator, "in2", label="4")
    c0 = Connection(evaporator, "out2", cc, "in1", label="0")
    w1 = Connection(water_in, "out1", condenser, "in2", label="w1")
    w2 = Connection(condenser, "out2", water_out, "in1", label="w2")
    b1 = Connection(brine_in, "out1", evaporator, "in1", label="b1")
    b2 = Connection(evaporator, "out1", brine_out, "in1", label="b2")
    nw.add_conns(c1, c2, c3, c4, c0, w1, w2, b1, b2)

    c1.set_attr(fluid={"R290": 1}, m=0.2, td_dew=5, T_dew=T_evap)
    c3.set_attr(td_bubble=3, T_bubble=T_cond)
    w1.set_attr(fluid={"water": 1}, p=2, T=T_cond - 15)
    w2.set_attr(T=T_cond - 5)
    b1.set_attr(fluid={"water": 1}, p=2, T=T_evap + 10)
    b2.set_attr(T=T_evap + 5)
    compressor.set_attr(eta_s=0.75)
    condenser.set_attr(pr1=1, pr2=1)
    evaporator.set_attr(pr1=1, pr2=1)
    nw.solve("design")
    assert nw.status == 0
    return nw


@pytest.fixture(scope="session")
def heat_pump():
    return build_heat_pump()


@pytest.fixture(autouse=True)
def _clean_plot_state():
    settings.reset_defaults()
    yield
    settings.reset_defaults()
    plt.close("all")
