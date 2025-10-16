from pytest import approx
from pytest import fixture

from tespy.components import PolynomialCompressorWithCooling
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network


@fixture
def polynomialcompressornetwork():

    nw = Network()
    nw.units.set_defaults(
        temperature="°C",
        pressure="bar"
    )

    gas_inlet = Source("gas inlet")
    gas_outlet = Sink("gas outlet")

    water_inlet = Source("water cold")
    water_outlet = Sink("water hot")

    compressor = PolynomialCompressorWithCooling("compressor")

    c1 = Connection(gas_inlet, "out1", compressor, "in1", label="c1")
    c2 = Connection(compressor, "out1", gas_outlet, "in1", label="c2")

    b1 = Connection(water_inlet, "out1", compressor, "in2", label="b1")
    b2 = Connection(compressor, "out2", water_outlet, "in1", label="b2")

    nw.add_conns(c1, c2, b1, b2)
    return nw


def test_setup_with_no_heat_dissipation(polynomialcompressornetwork):
    nw = polynomialcompressornetwork
    c1, c2, b1, b2 = nw.get_conn(["c1", "c2", "b1", "b2"])
    compressor = nw.get_comp("compressor")

    c1.set_attr(fluid={"R290": 1}, m=1, T_dew=10, td_dew=10)
    c2.set_attr(T_dew=60, td_dew=25)

    b1.set_attr(fluid={"water": 1}, T=15, p=1)
    b2.set_attr(T=40, p=0.9)
    compressor.set_attr(dissipation_ratio=0.0, eta_recovery=0.9)

    nw.solve("design")
    nw.assert_convergence()

    assert approx(b1.m.val) == 0.0


def test_setup_with_no_recovery(polynomialcompressornetwork):
    nw = polynomialcompressornetwork
    c1, c2, b1, b2 = nw.get_conn(["c1", "c2", "b1", "b2"])
    compressor = nw.get_comp("compressor")

    c1.set_attr(fluid={"R290": 1}, m=1, T_dew=10, td_dew=10)
    c2.set_attr(T_dew=60, td_dew=25)

    b1.set_attr(fluid={"water": 1}, T=15, p=1)
    b2.set_attr(T=40, p=0.9)
    compressor.set_attr(dissipation_ratio=0.1, eta_recovery=0.0)

    nw.solve("design")
    nw.assert_convergence()

    assert approx(b1.m.val) == 0.0


def test_cooling_energy_balance(polynomialcompressornetwork):
    nw = polynomialcompressornetwork
    c1, c2, b1, b2 = nw.get_conn(["c1", "c2", "b1", "b2"])
    compressor = nw.get_comp("compressor")

    c1.set_attr(fluid={"R290": 1}, m=1, T_dew=10, td_dew=10)
    c2.set_attr(T_dew=60, td_dew=25)

    b1.set_attr(fluid={"water": 1}, T=15, p=1)
    b2.set_attr(T=40, p=0.9)
    compressor.set_attr(dissipation_ratio=0.1, eta_recovery=0.8)

    nw.solve("design")
    nw.assert_convergence()

    heat = -b1.m.val_SI * (b2.h.val_SI - b1.h.val_SI)
    heat_other = compressor.Q_diss.val_SI * compressor.eta_recovery.val_SI
    assert approx(heat) == heat_other
