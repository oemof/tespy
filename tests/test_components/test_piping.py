# -*- coding: utf-8

"""Module for testing components of the piping module.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_piping.py

SPDX-License-Identifier: MIT
"""
import numpy as np

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
        msg = ('The pressure drop at the valve should be ' + str(dp) + ' but '
               'is ' + str(dp_act) + '.')
        assert dp == dp_act, msg

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
