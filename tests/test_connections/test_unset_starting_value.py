import numpy as np
import pytest

from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.networks import Network


@pytest.mark.parametrize("key", ["T0", "x0", "p0", "h0", "m0"])
def test_unset_starting_value_with_none(key):
    nw = Network(iterinfo=False)
    nw.units.set_defaults(temperature="degC", pressure="bar")
    so = Source("source")
    hx = SimpleHeatExchanger("heater")
    si = Sink("sink")
    c1 = Connection(so, "out1", hx, "in1", label="1")
    c2 = Connection(hx, "out1", si, "in1", label="2")
    nw.add_conns(c1, c2)
    c1.set_attr(fluid={"water": 1}, m=1, p=5, T=20)
    hx.set_attr(pr=1, Q=1e6)

    guesses = {"T0": 80, "x0": 0.5, "p0": 4, "h0": 400, "m0": 2}
    c2.set_attr(**{key: guesses[key]})
    c2.set_attr(**{key: None})
    assert np.isnan(c2.get_attr(key[:-1]).val0)

    nw.solve("design")
    nw.assert_convergence()
