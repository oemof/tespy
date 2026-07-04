# -*- coding: utf-8

"""Module for testing components of type reactor.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_reactors.py

SPDX-License-Identifier: MIT
"""
from pytest import approx

from tespy.components import FuelCell
from tespy.components import PowerSink
from tespy.components import PowerSource
from tespy.components import Sink
from tespy.components import Source
from tespy.components import WaterElectrolyzer
from tespy.connections import Connection
from tespy.connections import PowerConnection
from tespy.networks import Network


class TestWaterElectrolyzer:

    def setup_method(self):
        """Set up network for electrolyzer tests."""
        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC"
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

    def test_WaterElectrolyzer(self):
        """Test component properties of water electrolyzer."""
        # check PowerConnection:
        # power on component and PowerConnection must be identical
        self.nw.get_conn('h2o').set_attr(T=25, p=1)
        self.nw.get_conn('h2').set_attr(T=25)
        power_source = PowerSource('power source')
        e_power = PowerConnection(
            power_source, 'power', self.instance, 'power', label='e_power'
        )
        self.nw.add_conns(e_power)
        e_power.set_attr(E=2.5e6)

        self.nw.solve('design')
        self.nw.assert_convergence()
        assert self.nw.status == 0
        msg = (
            f"Value of power must be {e_power.E.val}, is {self.instance.P.val_SI}."
        )
        assert approx(e_power.E.val_SI) == self.instance.P.val_SI, msg

        # efficiency was set to 100 % with inlet and outlet states of the
        # reaction educts and products being identical to reference state
        # therefore Q must be equal to 0
        msg = f"Value of heat must be 0.0, is {self.instance.Q.val}."
        assert approx(self.instance.Q.val, abs=1e-4) == 0.0, msg

        # reset power, change efficiency value and specify heat output
        e_power.set_attr(E=None)
        self.nw.get_conn('h2o').set_attr(T=25, p=1)
        self.nw.get_conn('h2').set_attr(T=50)
        self.instance.set_attr(eta=0.8, Q=-8e5)

        self.nw.solve('design')
        self.nw.assert_convergence()
        msg = f"Value of heat must be {-8e5}, is {self.instance.Q.val}."
        assert approx(-8e5) == self.instance.Q.val, msg
        design_state = self.nw.save(as_dict=True)

        # check heat output constraint (offdesign test)
        Q = self.instance.Q.val * 0.9
        self.instance.set_attr(Q=Q)
        self.nw.solve('offdesign', design_path=design_state)
        self.nw.assert_convergence()
        msg = f"Value of heat must be {Q}, is {self.instance.Q.val}."
        assert approx(Q) == self.instance.Q.val, msg
        self.instance.set_attr(Q=None)

        # test efficiency vs. specific energy consumption
        self.nw.get_conn('h2').set_attr(m=0.1)
        self.instance.set_attr(eta=0.9, e='var')
        self.nw.solve('design')
        self.nw.assert_convergence()
        msg = (
            f"Value of efficiency must be {self.instance.eta.val}, is "
            f"{self.instance.e0 / self.instance.e.val}."
        )
        assert approx(self.instance.eta.val) == self.instance.e0 / self.instance.e.val, msg

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
        assert approx(self.instance.e0 / e) == self.instance.eta.val, msg
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

        # test cooling loop pressure ratio
        pr = 0.95
        self.instance.set_attr(pr=pr, e=None, eta=None, P=2e7)
        self.nw.solve('design')
        design_state = self.nw.save(as_dict=True)
        self.nw.assert_convergence()
        msg = (
            f"Value of pressure ratio must be {pr}, is {self.instance.pr.val}."
        )
        assert approx(pr) == self.instance.pr.val, msg

        # use zeta as offdesign parameter, at design point pressure
        # ratio must not change
        self.instance.set_attr(design=['pr'], offdesign=['zeta_d4'])
        self.nw.solve('offdesign', design_path=design_state)
        self.nw.assert_convergence()
        msg = (
            f"Value of pressure ratio must be {pr}, is {self.instance.pr.val}."
        )
        assert approx(pr) == self.instance.pr.val, msg

        # test heat output specification in offdesign mode
        Q = self.instance.Q.val * 0.9
        self.instance.set_attr(Q=Q, P=None)
        self.nw.solve('offdesign', design_path=design_state)
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


class TestFuelCell:

    def setup_method(self):
        """Set up network for fuel cell tests."""
        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC"
        })
        self.instance = FuelCell('fuel cell')

        cw_in = Source('cooling water in')
        cw_out = Sink('cooling water out')
        o2_source = Source('oxygen source')
        h2_source = Source('hydrogen source')
        water_sink = Sink('water sink')

        self.instance.set_attr()

        cw_fc = Connection(
            cw_in, 'out1', self.instance, 'in1', fluid={'H2O': 1}, T=20, p=1
        )
        fc_cw = Connection(self.instance, 'out1', cw_out, 'in1', T=40)
        self.nw.add_conns(cw_fc, fc_cw)

        o2_fc = Connection(
            o2_source, 'out1', self.instance, 'in2', label='o2', T=25, p=1
        )
        h2_fc = Connection(
            h2_source, 'out1', self.instance, 'in3', label='h2', T=25
        )
        fc_h2o = Connection(
            self.instance, 'out2', water_sink, 'in1', label='h2o', T=50
        )
        self.nw.add_conns(o2_fc, h2_fc, fc_h2o)

    def test_FuelCell(self):
        """Test component properties of fuel cell."""
        # test efficiency vs. specific energy consumption
        # eta = e / e0, with both e and e0 negative for the fuel cell
        self.nw.get_conn('h2').set_attr(m=0.01)
        self.instance.set_attr(pr=0.98)
        self.instance.set_attr(eta=0.45)
        self.nw.solve('design')
        self.nw.assert_convergence()
        eta = self.instance.eta.val
        eta_calc = self.instance.e.val / self.instance.e0
        msg = f"Value of efficiency must be {eta}, is {eta_calc}."
        assert approx(eta) == eta_calc, msg

        # at eta < 1 the fuel cell rejects heat to the cooling loop (Q < 0)
        msg = (
            f"Value of heat must be less than zero, is "
            f"{self.instance.Q.val}."
        )
        assert self.instance.Q.val < 0, msg

        # test specific energy consumption specification
        e = self.instance.e0 * 0.5
        self.instance.set_attr(e=None, eta=None)
        self.instance.set_attr(e=e)
        self.nw.solve('design')
        self.nw.assert_convergence()
        msg = (
            f"Value of specific energy consumption e must be {e}, is "
            f"{self.instance.e.val}."
        )
        assert approx(e) == self.instance.e.val, msg

        # test cooling loop pressure ratio
        pr = 0.95
        self.instance.set_attr(pr=pr, e=None, eta=None, P=-2e5)
        self.nw.solve('design')
        design_state = self.nw.save(as_dict=True)
        self.nw.assert_convergence()
        msg = (
            f"Value of pressure ratio must be {pr}, is {self.instance.pr.val}."
        )
        assert approx(pr) == self.instance.pr.val, msg

        # use zeta as offdesign parameter, at design point pressure
        # ratio must not change
        self.instance.set_attr(design=['pr'], offdesign=['zeta_d4'])
        self.nw.solve('offdesign', design_path=design_state)
        self.nw.assert_convergence()
        msg = (
            f"Value of pressure ratio must be {pr}, is {self.instance.pr.val}."
        )
        assert approx(pr) == self.instance.pr.val, msg

        # test heat output specification in offdesign mode
        Q = self.instance.Q.val * 0.9
        self.instance.set_attr(Q=Q, P=None)
        self.nw.solve('offdesign', design_path=design_state)
        self.nw.assert_convergence()
        msg = f"Value of heat must be {Q}, is {self.instance.Q.val}."
        assert approx(Q) == self.instance.Q.val, msg

    def test_FuelCell_PowerConnection(self):
        """Test utilization of PowerConnection with FuelCell."""
        self.instance.set_attr(eta=.5, pr=0.98)
        power_sink = PowerSink("grid")
        e1 = PowerConnection(
            self.instance, "power", power_sink, "power", label="e1"
        )
        self.nw.add_conns(e1)

        e1.set_attr(E=2e5)
        self.nw.solve('design')
        self.nw.assert_convergence()

        assert approx(-self.instance.P.val_SI) == e1.E.val_SI
