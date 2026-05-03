# -*- coding: utf-8

"""Module for testing components of type combustion.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_combustion.py

SPDX-License-Identifier: MIT
"""

import pytest

from tespy.components import CombustionChamber
from tespy.components import CombustionEngine
from tespy.components import DiabaticCombustionChamber
from tespy.components import Motor
from tespy.components import PowerSink
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.connections import PowerConnection
from tespy.networks import Network
from tespy.tools import CharLine
from tespy.tools.helpers import TESPyNetworkError


class TestCombustion:

    def setup_method(self):

        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "temperature": "degC",
            "volumetric_flow": "m3/s"
})
        self.fuel = Source('fuel')
        self.air = Source('ambient air')
        self.fg = Sink('flue gas')

    def setup_CombustionChamber_network(self, instance):

        self.c1 = Connection(self.air, 'out1', instance, 'in1', label="air")
        self.c2 = Connection(self.fuel, 'out1', instance, 'in2', label="fuel")
        self.c3 = Connection(instance, 'out1', self.fg, 'in1', label="fluegas")
        self.nw.add_conns(self.c1, self.c2, self.c3)

    def setup_CombustionEngine_network(self, instance):

        self.cw1_in = Source('cooling water 1 source')
        self.cw2_in = Source('cooling water 2 source')
        self.cw1_out = Sink('cooling water 1 sink')
        self.cw2_out = Sink('cooling water 2 sink')

        self.c1 = Connection(self.air, 'out1', instance, 'in3')
        self.c2 = Connection(self.fuel, 'out1', instance, 'in4')
        self.c3 = Connection(instance, 'out3', self.fg, 'in1')
        self.c4 = Connection(self.cw1_in, 'out1', instance, 'in1')
        self.c5 = Connection(self.cw2_in, 'out1', instance, 'in2')
        self.c6 = Connection(instance, 'out1', self.cw1_out, 'in1')
        self.c7 = Connection(instance, 'out2', self.cw2_out, 'in1')
        self.nw.add_conns(
            self.c1, self.c2, self.c3, self.c4, self.c5, self.c6, self.c7
        )

        motor = Motor("motor")
        grid = PowerSink("grid")
        self.e1 = PowerConnection(instance, "power", motor, "power_in")
        self.e2 = PowerConnection(motor, "power_out", grid, "power")
        self.nw.add_conns(self.e1, self.e2)

    def test_CombustionChamber(self):
        """
        Test component properties of combustion chamber.
        """
        instance = CombustionChamber('combustion chamber')
        self.setup_CombustionChamber_network(instance)

        # connection parameter specification
        air = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129}
        fuel = {'CO2': 0.04, 'CH4': 0.96}
        self.c1.set_attr(fluid=air, p=1, T=30)
        self.c2.set_attr(fluid=fuel, T=30)
        self.c3.set_attr(T=1200)

        # test specified thermal input for CombustionChamber
        instance.set_attr(ti=1e6)
        self.nw.solve('design')
        self.nw.assert_convergence()
        ti = (
            self.c2.m.val_SI * self.c2.fluid.val['CH4']
            * instance.fuels['CH4']['LHV']
        )
        msg = f'Value of thermal input must be {instance.ti.val}, is {ti}.'
        assert round(ti, 1) == round(instance.ti.val, 1), msg

        # test specified lamb for CombustionChamber
        self.c3.set_attr(T=None)
        instance.set_attr(lamb=1)
        self.nw.solve('design')
        self.nw.assert_convergence()
        o2 = round(self.c3.fluid.val['O2'], 4)
        msg = f'Value of oxygen in flue gas must be 0.0, is {o2}.'
        assert 0.0 == o2, msg

    def test_CombustionChamberCarbonMonoxide(self):
        instance = CombustionChamber('combustion chamber')
        self.setup_CombustionChamber_network(instance)

        # connection parameter specification
        air = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129}
        fuel = {'CO': 1}
        self.c1.set_attr(fluid=air, p=1, T=30)
        self.c2.set_attr(fluid=fuel, T=30, m=1)
        instance.set_attr(lamb=3)

        self.nw.solve('design')
        self.nw.assert_convergence()
        assert self.nw.status == 0
        assert instance.fuels["CO"]["LHV"] == pytest.approx(10112000, 1e-3)

        molar_flow = {}
        for c in self.nw.conns["object"]:
            molar_flow[c.label] = {
                key: value["mass_fraction"]
                / value["wrapper"]._molar_mass * c.m.val_SI
                for key, value in c.fluid_data.items()
            }
        o2 = molar_flow["air"]["O2"]
        co2 = molar_flow["fluegas"]["CO2"]
        assert o2 == pytest.approx(co2 * 1.5, 1e-3)

        self.c3.set_attr(T=1500)
        instance.set_attr(lamb=None)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert self.c3.T.val == pytest.approx(1500)

    def test_DiabaticCombustionChamber_H2(self):
        instance = DiabaticCombustionChamber('combustion chamber')
        self.setup_CombustionChamber_network(instance)

        instance.set_attr(pr=0.95, eta=1)

        air = {'O2': 0.21, 'N2': 0.79}
        fuel = {'H2': 1}

        self.c1.set_attr(fluid=air, m=100, p=1, T=20)
        self.c2.set_attr(fluid=fuel, p=1, T=20)
        self.c3.set_attr(T=1000)

        self.nw.solve(mode='design')
        self.nw.assert_convergence()
        assert self.c3.fluid.val["O2"] > 0

    def test_CombustionChamberHighTemperature(self):
        instance = CombustionChamber('combustion chamber')
        self.setup_CombustionChamber_network(instance)

        # connection parameter specification
        air = {'N2': 0.8, 'O2': 0.2}
        fuel = {'CH4': 1}
        self.c1.set_attr(fluid=air, p=1, T=30, m=1)
        self.c2.set_attr(fluid=fuel, T=30)
        self.c3.set_attr(T=1200)
        self.nw.solve('design')
        assert self.nw.status == 0
        # a good guess for outlet fluid composition is necessary
        self.c3.set_attr(T=None)
        instance.set_attr(lamb=1)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert self.c3.T.val_SI == pytest.approx(2110, abs=0.1)

    def test_DiabaticCombustionChamber(self):
        """
        Test component properties of diabatic combustion chamber.
        """
        instance = DiabaticCombustionChamber('combustion chamber')
        self.setup_CombustionChamber_network(instance)

        # connection parameter specification
        air = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129}
        fuel = {'CO2': 0.04, 'CH4': 0.96}
        self.c1.set_attr(fluid=air, p=1.2, T=30)
        self.c2.set_attr(fluid=fuel, T=30, p=1.5)
        self.c3.set_attr(T=1200)

        pr = 0.97
        instance.set_attr(pr=pr, eta=0.95, ti=1e6)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert self.nw.status == 0

        valid = round(self.c1.p.val * pr, 2)
        check = round(self.c3.p.val, 2)

        # test outlet pressure value
        msg = (
            f'Value of outlet pressure must be {valid}, the actual value is '
            f'{check}.'
        )
        assert valid == check, msg

        # test invalid pressure ratio
        instance.set_attr(pr=None)
        self.c1.set_attr(p=1.2)
        self.c2.set_attr(p=1.5)
        self.c3.set_attr(p=1.3)
        self.nw.solve('design')
        self.nw.assert_convergence()

        valid = round(self.c3.p.val / self.c1.p.val, 2)
        check = round(instance.pr.val, 2)
        msg = (
            f'Value of pressure ratio must be {valid}, the actual value is '
            f'{check}.'
        )
        assert valid == check, msg

        # test invalid pressure specification -> leading to linear dependency
        instance.set_attr(pr=pr)
        self.c2.set_attr(p=None)
        self.c3.set_attr(p=1.3)
        with pytest.raises(TESPyNetworkError):
            self.nw.solve('design')

    def test_CombustionEngine(self):
        """Test component properties of combustion engine."""
        instance = CombustionEngine('combustion engine')
        self.setup_CombustionEngine_network(instance)
        eta_char = CharLine(x=[0.5, 0.75, 1.0], y=[0.9, 0.98, 1])
        self.nw.get_comp("motor").set_attr(
            eta=0.98, design=["eta"], offdesign=["eta_char"],
            eta_char=eta_char
        )

        air = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129}
        fuel = {'CO2': 0.04, 'CH4': 0.96}
        water1 = {'H2O': 1}
        water2 = {'H2O': 1}

        # connection parametrisation
        instance.set_attr(
            pr1=0.99, pr2=0.99, lamb=1.0,
            design=['pr1', 'pr2'], offdesign=['zeta1', 'zeta2']
        )
        self.c1.set_attr(p=5, T=30, fluid=air)
        self.c2.set_attr(T=30, fluid=fuel)
        self.c4.set_attr(p=3, T=60, m=50, fluid=water1)
        self.c5.set_attr(p=3, T=80, m=50, fluid=water2)

        # test specified thermal input
        ti = 1e6
        instance.set_attr(ti=ti)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert self.nw.status == 0
        design_state = self.nw.save(as_dict=True)
        # calculate in offdesign mode
        self.nw.solve('offdesign', design_path=design_state)
        self.nw.assert_convergence()
        msg = f'Value of thermal input must be {ti}, is {instance.ti.val}.'
        assert round(ti, 1) == round(instance.ti.val, 1), msg
        instance.set_attr(ti=None)

        # test specified heat output 1
        instance.set_attr(Q1=instance.Q1.val)
        self.nw.solve('offdesign', design_path=design_state)
        self.nw.assert_convergence()
        # heat output is at design point value, thermal input must therefore
        # not have changed
        msg = f'Value of thermal input must be {ti}, is {instance.ti.val}.'
        assert round(ti, 1) == round(instance.ti.val, 1), msg

        # calculate heat output over cooling loop
        heat1 = self.c4.m.val_SI * (self.c6.h.val_SI - self.c4.h.val_SI)
        msg = f'Value of heat output 1 must be {-heat1}, is {instance.Q1.val}.'
        assert round(heat1, 1) == -round(instance.Q1.val, 1), msg
        instance.set_attr(Q1=None)

        # test specified heat output 2
        instance.set_attr(Q2=1.2 * instance.Q2.val)
        self.nw.solve('offdesign', design_path=design_state)
        self.nw.assert_convergence()

        # calculate heat output over cooling loop
        heat2 = self.c5.m.val_SI * (self.c7.h.val_SI - self.c5.h.val_SI)
        msg = f'Value of heat output 2 must be {-heat2}, is {instance.Q2.val}.'
        assert round(heat2, 1) == -round(instance.Q2.val, 1), msg

        # test specified heat loss
        instance.set_attr(Q2=None, Qloss=-1e5)
        self.nw.solve('offdesign', design_path=design_state)
        self.nw.assert_convergence()
        msg = f'Value of heat loss must be {-1e5}, is {instance.Qloss.val}.'
        assert round(-1e5, 1) == round(instance.Qloss.val, 1), msg

        # test connector specification
        instance.set_attr(Qloss=None)
        self.e2.set_attr(E=400e5)
        self.nw.solve('offdesign', design_path=design_state)
        self.nw.assert_convergence()
        msg = (
            f'Value of power must be {-self.e1.E.val}, is '
            f'{instance.P.val}.'
        )
        assert round(-self.e1.E.val, 1) == round(instance.P.val, 1), msg
