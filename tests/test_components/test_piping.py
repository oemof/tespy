# -*- coding: utf-8

"""Module for testing components of the piping module.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_piping.py

SPDX-License-Identifier: MIT
"""
import numpy as np
from pytest import approx

from tespy.components import Pipe
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Valve
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.characteristics import CharLine


class TestPiping:

    def setup_piping_network(self, instance):
        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "temperature": "degC"
        })
        self.source = Source('source')
        self.sink = Sink('sink')
        self.c1 = Connection(self.source, 'out1', instance, 'in1')
        self.c2 = Connection(instance, 'out1', self.sink, 'in1')
        self.nw.add_conns(self.c1, self.c2)

    def test_Valve(self):
        """Test component properties of valves."""
        instance = Valve('valve')
        self.setup_piping_network(instance)

        # parameter specification
        self.c1.set_attr(fluid={'CH4': 1}, m=10, p=10, T=120)
        self.c2.set_attr(p=1)

        self.nw.solve('design')
        self.nw.assert_convergence()
        assert self.nw.status == 0
        pr = round(self.c2.p.val / self.c1.p.val, 2)
        msg = (
            f'Value of pressure ratio must be {pr}, is {instance.pr.val}.'
        )
        assert pr == round(instance.pr.val, 2), msg

        # test variable zeta value
        zeta = round(instance.zeta.val, 0)
        instance.set_attr(zeta='var', pr=None)
        self.nw.solve('design')
        self.nw.assert_convergence()
        msg = (
            f'Value of dimension independent zeta value must be {zeta}, is '
            f'{instance.zeta.val}.'
        )
        assert zeta == round(instance.zeta.val, 0), msg

        # dp char
        x = np.array([8, 9, 10, 11, 12])
        y = np.array([5, 8, 9, 9.5, 9.6]) * 1e5
        dp_char = CharLine(x, y)
        instance.set_attr(zeta=None, dp_char={
            'char_func': dp_char, 'is_set': True})
        m = 11
        self.c1.set_attr(m=m)
        self.c2.set_attr(p=None)
        self.nw.solve('design')
        self.nw.assert_convergence()
        self.nw.print_results()
        dp = round(-dp_char.evaluate(m), 0)
        dp_act = round(self.c2.p.val_SI - self.c1.p.val_SI)
        msg = f"The pressure drop at the valve should be {dp} but is {dp_act}."
        assert dp == dp_act, msg

    def test_Valve_Kv(self):
        instance = Valve('valve')
        self.setup_piping_network(instance)

        # parameter specification
        # mass and volumetric flow 5 at the same time to find temperature that
        # exactly is according to 1000 kg/m3 density
        self.c1.set_attr(fluid={'H2O': 1}, p=5, m=5000 / 3600, v=5 / 3600)
        # one bar pressure loss
        self.c2.set_attr(p=4)
        self.nw.solve("design")

        # volumetric flow must then be exactly equal to Kv
        assert approx(self.c1.v.val * 3600) == instance.Kv.val_SI

        # data from online available resource:
        # https://www.samsongroup.com/document/t80003en.pdf
        kv_data = np.array([
            0.09,0.63,1.1,2.1,3.1,4.2,5.2,6.2,7.2,8.2,9.2,10.3,11.3
        ])
        opening_data = np.array([0,5,10,20,30,40,50,60,70,80,90,100,110]) / 100
        Kv_char = {
            "char_func": CharLine(x=opening_data, y=kv_data) , "is_set": True
        }
        # 5 -> 0.4 + (5 - 4.2) / (5.2 - 4.2) * (0.5 - 0.4)
        instance.set_attr(Kv_char=Kv_char, opening="var")
        self.nw.solve("design")
        assert approx(instance.opening.val) == 0.48

    def test_Valve_Kv_analytical(self):
        instance = Valve('valve')
        self.setup_piping_network(instance)

        def analytical_Kv(opening, param):
            return param * opening

        # parameter specification
        # mass and volumetric flow 5 at the same time to find temperature that
        # exactly is according to 1000 kg/m3 density
        self.c1.set_attr(fluid={'H2O': 1}, p=5, m=5000 / 3600, v=5 / 3600)
        # one bar pressure loss
        self.c2.set_attr(p=4)
        instance.set_attr(
            opening="var",
            Kv_analytical={"method": analytical_Kv, "params": [10]}
        )

        self.nw.solve("design")
        assert approx(instance.opening.val_SI) == 0.5

    def test_Pipe(self):
        """Test component properties of pipe."""
        instance = Pipe('pipe')
        self.setup_piping_network(instance)

        # parameter specification
        self.c1.set_attr(fluid={'H2O': 1}, m=10, p=10, T=220)

        instance.set_attr(
            pr=0.99, Tamb=20, L=1000, D=0.2,
            insulation_thickness=0.1 ,insulation_tc=0.035,
            pipe_thickness=0.002, material='Steel',
            environment_media='air', wind_velocity=2,
        )
        self.nw.solve('design')
        self.nw.assert_convergence()
        Q = -62683.7
        msg = (
            f"The Heat loss of surface pipe should be {Q} but is "
            f"{round(instance.Q.val, 1)}."
        )
        assert Q == round(instance.Q.val, 1), msg

        instance.set_attr(
            wind_velocity=None, environment_media='moist soil', pipe_depth=5
            )
        self.nw.solve('design')
        self.nw.assert_convergence()
        Q = -57081.9
        msg = (
            f"The Heat loss of sub surface pipe should be {Q} but is "
            f"{round(instance.Q.val, 1)}."
        )
        assert Q == round(instance.Q.val, 1), msg
