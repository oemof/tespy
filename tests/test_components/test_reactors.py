# -*- coding: utf-8

"""Module for testing components of type reactor.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_reactors.py

SPDX-License-Identifier: MIT
"""
from pytest import approx

from tespy.components import Sink
from tespy.components import Source
from tespy.components import PowerSource
from tespy.components import WaterElectrolyzer
from tespy.connections import Bus
from tespy.connections import Connection
from tespy.connections import PowerConnection
from tespy.networks import Network


class TestWaterElectrolyzer:

    def setup_method(self):
        """Set up network for electrolyzer tests."""
        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "temperature": "degC"
        })
        self.instance = WaterElectrolyzer('electrolyzer')

        fw = Source('feed water')
        cw_in = Source('cooling water')
        o2 = Sink('oxygen sink')
        h2 = Sink('hydrogen sink')
        cw_out = Sink('cooling water sink')

        self.instance.set_attr(pr=0.99, eta=1)

        cw_el = Connection(
            cw_in, 'out1', self.instance, 'in1', fluid={'H2O': 1}, T=20, p=1
       )
        el_cw = Connection(self.instance, 'out1', cw_out, 'in1', T=45)

        self.nw.add_conns(cw_el, el_cw)

        fw_el = Connection(fw, 'out1', self.instance, 'in2', label='h2o')
        el_o2 = Connection(self.instance, 'out2', o2, 'in1')
        el_h2 = Connection(self.instance, 'out3', h2, 'in1', label='h2')

        self.nw.add_conns(fw_el, el_o2, el_h2)

    def test_WaterElectrolyzer(self, tmp_path):
        """Test component properties of water electrolyzer."""
        tmp_path = f'{tmp_path}.json'
        # check bus function:
        # power output on component and bus must be indentical
        self.nw.get_conn('h2o').set_attr(T=25, p=1)
        self.nw.get_conn('h2').set_attr(T=25)
        power = Bus('power')
        power.add_comps({'comp': self.instance, 'param': 'P', 'base': 'bus'})
        power.set_attr(P=2.5e6)
        self.nw.add_busses(power)

        self.nw.solve('design')
        self.nw.assert_convergence()
        assert self.nw.status == 0
        msg = (
            f"Value of power must be {power.P.val}, is {self.instance.P.val}."
        )
        assert round(power.P.val, 1) == round(self.instance.P.val), msg

        # effieciency was set to 100 % with inlet and outlet states of the
        # reaction educts and products beeing identical to reference state
        # therefore Q must be equal to 0
        msg = f"Value of heat must be 0.0, is {self.instance.Q.val}."
        assert round(self.instance.Q.val, 4) == 0.0, msg

        # reset power, change efficiency value and specify heat bus value
        power.set_attr(P=None)
        self.nw.get_conn('h2o').set_attr(T=25, p=1)
        self.nw.get_conn('h2').set_attr(T=50)
        self.instance.set_attr(eta=0.8)
        # check bus function:
        # heat output on component and bus must be indentical
        heat = Bus('heat')
        heat.add_comps({'comp': self.instance, 'param': 'Q'})
        heat.set_attr(P=-8e5)
        self.nw.add_busses(heat)

        self.nw.solve('design')
        self.nw.assert_convergence()
        msg = (
            f"Value of heat must be {heat.P.val}, is {self.instance.Q.val}."
        )
        assert round(heat.P.val, 1) == round(self.instance.Q.val), msg
        self.nw.save(tmp_path)

        # check bus function:
        # heat output on component and bus must identical (offdesign test)
        Q = heat.P.val * 0.9
        heat.set_attr(P=Q)
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw.assert_convergence()
        msg = f"Value of heat must be {Q}, is {self.instance.Q.val}."
        assert round(Q, 1) == round(self.instance.Q.val), msg

        # delete both busses again
        self.nw.del_busses(heat, power)

        # test efficiency vs. specific energy consumption
        self.nw.get_conn('h2').set_attr(m=0.1)
        self.instance.set_attr(eta=0.9, e='var')
        self.nw.solve('design')
        self.nw.assert_convergence()
        msg = (
            f"Value of efficiency must be {self.instance.eta.val}, is "
            f"{self.instance.e0 / self.instance.e.val}."
        )
        eta = round(self.instance.eta.val, 2)
        eta_calc = round(self.instance.e0 / self.instance.e.val, 2)
        assert eta == eta_calc, msg

        # test efficiency value > 1, Q must be larger than 0
        e = 130e6
        self.instance.set_attr(e=None, eta=None)
        self.instance.set_attr(e=e)
        self.nw.solve('design')
        self.nw.assert_convergence()
        # test efficiency
        msg = (
            f"Value of efficiency must be {self.instance.e0 / e}, is "
            f"{self.instance.eta.val}."
        )
        eta = round(self.instance.e0 / e, 2)
        eta_calc = round(self.instance.eta.val, 2)
        assert eta == eta_calc, msg
        # test Q
        msg = (
            f"Value of heat must be larger than zero, is "
            f"{self.instance.Q.val}."
        )
        assert self.instance.Q.val > 0, msg

        # test specific energy consumption
        e = 150e6
        self.instance.set_attr(e=None, eta=None)
        self.instance.set_attr(e=e)
        self.nw.solve('design')
        self.nw.assert_convergence()
        msg = (
            f"Value of specific energy consumption e must be {e}, is "
            f"{self.instance.e.val}."
        )
        assert approx(e) == self.instance.e.val, msg

        # test cooling loop pressure ratio, zeta as variable value
        pr = 0.95
        self.instance.set_attr(
            pr=pr, e=None, eta=None, zeta='var', P=2e7, design=['pr'])
        self.nw.solve('design')
        self.nw.save(tmp_path)
        self.nw.assert_convergence()
        msg = (
            f"Value of pressure ratio must be {pr}, is {self.instance.pr.val}."
        )
        assert approx(pr) == self.instance.pr.val, msg

        # use zeta as offdesign parameter, at design point pressure
        # ratio must not change
        self.instance.set_attr(zeta=None, offdesign=['zeta'])
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw.assert_convergence()
        msg = (
            f"Value of pressure ratio must be {pr}, is {self.instance.pr.val}."
        )
        assert approx(pr) == self.instance.pr.val, msg

        # test heat output specification in offdesign mode
        Q = self.instance.Q.val * 0.9
        self.instance.set_attr(Q=Q, P=None)
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw.assert_convergence()
        msg = f"Value of heat must be {Q}, is {self.instance.Q.val}."
        assert approx(Q) == self.instance.Q.val, msg

    def test_WaterElectrolyzer_PowerConnection(self):
        """test utilization of PowerConnection with WaterElectrolyzer"""
        self.nw.get_conn('h2o').set_attr(T=25, p=1)
        self.nw.get_conn('h2').set_attr(T=25)

        power = PowerSource("grid")
        e1 = PowerConnection(power, "power", self.instance, "power", label="e1")
        self.nw.add_conns(e1)

        e1.set_attr(E=2.5e6)
        self.nw.solve('design')

        assert approx(self.instance.P.val_SI) == e1.E.val_SI
