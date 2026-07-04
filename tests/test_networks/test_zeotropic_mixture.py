# -*- coding: utf-8

"""Integration tests for zeotropic mixture networks.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_zeotropic_mixture.py

SPDX-License-Identifier: MIT
"""
from CoolProp.CoolProp import get_global_param_string
from pytest import approx
from pytest import mark

from tespy.components import CycleCloser
from tespy.components import Generator
from tespy.components import Motor
from tespy.components import PowerBus
from tespy.components import PowerSink
from tespy.components import Pump
from tespy.components import SectionedHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Turbine
from tespy.connections import Connection
from tespy.connections import PowerConnection
from tespy.networks import Network
from tespy.tools.fluid_properties.functions import T_bubble_p
from tespy.tools.fluid_properties.functions import T_dew_p
from tespy.tools.fluid_properties.functions import p_sat_TQ

_WF = "REFPROP::ISOBUTAN[0.5]&IPENTANE[0.5]|mass"

skipif_no_refprop = mark.skipif(
    get_global_param_string("REFPROP_version") == "n/a",
    reason="This test requires REFPROP.",
)


@skipif_no_refprop
class TestZeotropicORC:
    """Simple ORC with a zeotropic isobutane/isopentane mixture.

    Topology:
        CycleCloser -> Turbine -> SectionedHeatExchanger (condenser)
        -> Pump -> SectionedHeatExchanger (evaporator) -> CycleCloser

    Heat source: air 200 -> ~155 degC through evaporator (td_pinch=10 K).
    Heat sink:   air  10 ->   20 degC through condenser.
    """

    def setup_method(self):
        nw = Network()
        nw.units.set_defaults(
            temperature="degC",
            pressure="bar",
            pressure_difference="bar",
            power="kW",
            heat="kW"
        )

        turbine = Turbine("turbine")
        recuperator = SectionedHeatExchanger("recuperator")
        condenser = SectionedHeatExchanger("condenser")
        pump = Pump("pump")
        preheater = SectionedHeatExchanger("preheater")
        evaporator = SectionedHeatExchanger("evaporator")
        cc = CycleCloser("cc")

        heat_source = Source("heat source")
        heat_outflow = Sink("heat outflow")

        air_source = Source("air source")
        air_sink = Sink("air sink")

        a1 = Connection(heat_source, "out1", evaporator, "in1", label="a1")
        a2 = Connection(evaporator, "out1", preheater, "in1", label="a2")
        a3 = Connection(preheater, "out1", heat_outflow, "in1", label="a3")

        b1 = Connection(cc, "out1", turbine, "in1", label="b1")
        b2 = Connection(turbine, "out1", recuperator, "in1", label="b2")
        b3 = Connection(recuperator, "out1", condenser, "in1", label="b3")
        b4 = Connection(condenser, "out1", pump, "in1", label="b4")
        b5 = Connection(pump, "out1", recuperator, "in2", label="b5")
        b6 = Connection(recuperator, "out2", preheater, "in2", label="b6")
        b7 = Connection(preheater, "out2", evaporator, "in2", label="b7")
        b8 = Connection(evaporator, "out2", cc, "in1", label="b8")

        c1 = Connection(air_source, "out1", condenser, "in2", label="c1")
        c2 = Connection(condenser, "out2", air_sink, "in1", label="c2")

        nw.add_conns(a1, a2, a3, b1, b2, b3, b4, b5, b6, b7, b8, c1, c2)

        generator = Generator("generator")
        motor = Motor("motor")
        power_bus = PowerBus("bus", num_in=1, num_out=2)
        grid = PowerSink("grid")

        e1 = PowerConnection(turbine, "power", generator, "power_in", label="e1")
        e2 = PowerConnection(generator, "power_out", power_bus, "power_in1", label="e2")
        e3 = PowerConnection(power_bus, "power_out1", motor, "power_in", label="e3")
        e4 = PowerConnection(motor, "power_out", pump, "power", label="e4")
        e5 = PowerConnection(power_bus, "power_out2", grid, "power", label="e5")

        nw.add_conns(e1, e2, e3, e4, e5)

        generator.set_attr(eta=0.98)
        motor.set_attr(eta=0.98)

        a1.set_attr(fluid={"air": 1}, T=200, p=1, m=10)

        b1.set_attr(fluid={"REFPROP::ISOBUTAN[0.5]&IPENTANE[0.5]|mass": 1}, x=1, T=150)
        b3.set_attr(td_dew=10)
        b4.set_attr(td_bubble=5)
        b7.set_attr(td_bubble=5)

        c1.set_attr(fluid={"air": 1}, T=10, p=1)
        c2.set_attr(T=20)

        recuperator.set_attr(dp1=0, dp2=0)
        condenser.set_attr(dp1=0, dp2=0)
        preheater.set_attr(dp1=0, dp2=0)
        evaporator.set_attr(dp1=0, dp2=0)

        turbine.set_attr(eta_s=0.8)
        pump.set_attr(eta_s=0.7)

        condenser.set_attr(td_pinch=5)
        evaporator.set_attr(td_pinch=10)

        nw.solve("design")

        self.nw = nw
        self.b1 = b1
        self.b4 = b4

    def test_network_converges(self):
        self.nw.assert_convergence()

    def test_turbine_inlet_pressure_matches_p_sat_TQ(self):
        """Turbine inlet pressure equals the dew-point pressure at T=150 degC (x=1 presolve fix)."""
        T_K = self.b1.T.val_SI
        p_expected = p_sat_TQ(T_K, 1, self.b1.fluid_data)
        assert self.b1.p.val_SI == approx(p_expected, rel=1e-4)

    def test_temperature_glide_at_condenser_pressure(self):
        """Dew temperature exceeds bubble temperature at condenser pressure (zeotropic glide)."""
        assert self.b4.T_dew.val_SI > self.b4.T_bubble.val_SI

    def test_pump_inlet_is_subcooled(self):
        """Pump inlet (b4) is subcooled: temperature below bubble point at that pressure."""
        assert self.b4.T.val_SI < T_bubble_p(self.b4.p.val_SI, self.b4.fluid_data)
