# -*- coding: utf-8 -*-
import numpy as np

from tespy.components import (
    Compressor, CycleCloser, HeatExchanger, Motor,
    MovingBoundaryHeatExchanger, PowerBus, PowerSource,
    Pump, Sink, Source, Valve
)
from tespy.connections import Connection, PowerConnection, Ref
from tespy.models import ModelTemplate
from tespy.networks import Network
from tespy.tools.characteristics import CharLine, load_default_char as ldc


class HeatPumpModel(ModelTemplate):
    """Ground-source heat pump with geothermal collector and heating system.

    Parameters
    ----------
    fluid : str
        Working fluid for the refrigerant cycle, e.g. "NH3" or "R410A".
    Tamb : float
        Ambient (reference) temperature in °C.
    pamb : float
        Ambient (reference) pressure in bar.
    Tgeo : float
        Mean geothermal temperature in °C.
    T_hs_feed : float
        Heating system feed flow temperature in °C.
    T_hs_ret : float
        Heating system return flow temperature in °C.
    """
    def __init__(
        self, fluid,
        Tamb=2.8, pamb=1.013, Tgeo=9.5,
        T_hs_feed=40
    ):
        self.fluid = fluid
        self.Tamb = Tamb
        self.pamb = pamb
        self.Tgeo = Tgeo
        self.T_hs_feed = T_hs_feed
        self._ean = None
        super().__init__()

    # ------------------------------------------------------------------
    # ModelTemplate interface
    # ------------------------------------------------------------------

    def _parameter_lookup(self):
        return {
            "T_geo": {"get": self._get_T_geo, "set": self._set_T_geo},
            "T_hs":  ["Connections", "c23", "T"],
            "Q":     ["Components", "condenser", "Q"],
            "COP":     {"get": self._calc_cop},
            "epsilon": {"get": self._get_epsilon},
        }

    def _create_network(self):
        self.nw = Network()
        self.nw.units.set_defaults(
            temperature="degC", pressure="bar", enthalpy="kJ/kg",
            pressure_difference="bar"
        )

        # components
        cc = CycleCloser("cycle closer")
        cd = MovingBoundaryHeatExchanger("condenser")
        va = Valve("valve")
        ev = HeatExchanger("evaporator")
        cp = Compressor("compressor")
        gh_in = Source("ground heat feed flow")
        gh_out = Sink("ground heat return flow")
        ghp = Pump("ground heat loop pump")
        hs_feed = Sink("heating system feed flow")
        hs_ret = Source("heating system return flow")
        hsp = Pump("heating system pump")

        # refrigerant cycle
        c1 = Connection(cc, "out1", cd, "in1", label="c1")
        c2 = Connection(cd, "out1", va, "in1", label="c2")
        c3 = Connection(va, "out1", ev, "in2", label="c3")
        c4 = Connection(ev, "out2", cp, "in1", label="c4")
        c5 = Connection(cp, "out1", cc, "in1", label="c5")
        self.nw.add_conns(c1, c2, c3, c4, c5)

        # geothermal circuit (boundary labels used in exergy analysis)
        c11 = Connection(gh_in, "out1", ghp, "in1", label="c11")
        c12 = Connection(ghp, "out1", ev, "in1", label="c12")
        c13 = Connection(ev, "out1", gh_out, "in1", label="c13")
        self.nw.add_conns(c11, c12, c13)

        # heating circuit (boundary labels used in exergy analysis)
        c21 = Connection(hs_ret, "out1", hsp, "in1", label="c21")
        c22 = Connection(hsp, "out1", cd, "in2", label="c22")
        c23 = Connection(cd, "out2", hs_feed, "in1", label="c23")
        self.nw.add_conns(c21, c22, c23)

        # component parametrization
        cd.set_attr(
            pr1=0.99, pr2=0.99, Q=-4e3,
            design=["pr2", "td_pinch"], offdesign=["zeta2_d4", "UA_char"],
        )
        UA_char1 = ldc("HeatExchanger", "kA_char1", "DEFAULT", CharLine)
        UA_char2 = ldc("HeatExchanger", "kA_char2", "EVAPORATING FLUID", CharLine)
        ev.set_attr(
            pr1=0.99, pr2=0.99,
            UA_char1=UA_char1, UA_char2=UA_char2,
            design=["pr1", "ttd_l"], offdesign=["zeta1_d4", "UA_char"],
        )
        cp.set_attr(eta_s=0.8,  design=["eta_s"], offdesign=["eta_s_char"])
        hsp.set_attr(eta_s=0.75, design=["eta_s"], offdesign=["eta_s_char"])
        ghp.set_attr(eta_s=0.75, design=["eta_s"], offdesign=["eta_s_char"])

        # connection parametrization
        c1.set_attr(fluid={self.fluid: 1})
        c2.set_attr(x=0, T_bubble=40)
        c4.set_attr(td_dew=3, T_dew=5)
        c11.set_attr(T=self.Tgeo + 1.5, p=1.5, fluid={"water": 1})
        c13.set_attr(T=Ref(c11, 1, -3), p=1.5)
        c21.set_attr(T=Ref(c23, 1, -5), p=2, fluid={"water": 1})
        c23.set_attr(T=self.T_hs_feed, p=2)

        self.nw.solve("design")

        c2.set_attr(T_bubble=None)
        c4.set_attr(T_dew=None)
        cd.set_attr(td_pinch=5)
        ev.set_attr(ttd_l=5)

        # motor + power distribution
        mot_cp = Motor("motor compressor")
        mot_ghp = Motor("motor ground heat pump")
        mot_hsp = Motor("motor heating system pump")

        x = np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
        y = np.array([0.2, 0.86, 0.9, 0.93, 0.95, 0.96, 0.95]) / 0.96
        char = CharLine(x=x, y=y)
        mot_cp.set_attr(eta=0.96,  eta_char=char, design=["eta"], offdesign=["eta_char"])
        mot_ghp.set_attr(eta=0.96, eta_char=char, design=["eta"], offdesign=["eta_char"])
        mot_hsp.set_attr(eta=0.96, eta_char=char, design=["eta"], offdesign=["eta_char"])

        grid = PowerSource("grid")
        power_dist = PowerBus("power distribution", num_in=1, num_out=3)

        e1 = PowerConnection(grid, "power", power_dist, "power_in1", label="e1")
        e2 = PowerConnection(power_dist, "power_out1", mot_cp, "power_in", label="e2")
        e3 = PowerConnection(mot_cp, "power_out", cp, "power", label="e3")
        e4 = PowerConnection(power_dist, "power_out2", mot_ghp, "power_in", label="e4")
        e5 = PowerConnection(mot_ghp, "power_out", ghp, "power", label="e5")
        e6 = PowerConnection(power_dist, "power_out3", mot_hsp, "power_in", label="e6")
        e7 = PowerConnection(mot_hsp, "power_out", hsp, "power", label="e7")
        self.nw.add_conns(e1, e2, e3, e4, e5, e6, e7)

        self.nw.solve("design")

    def _set_T_geo(self, value):
        self.nw.get_conn("c11").set_attr(T=value - 1)

    def _get_T_geo(self):
        return self.nw.get_conn("c11").T.val + 1

    def _calc_cop(self):
        cd = self.nw.get_comp("condenser")
        return abs(cd.Q.val_SI) / self.nw.get_conn("e1").E.val_SI

    def _get_epsilon(self):
        return None if self._ean is None else self._ean.epsilon
