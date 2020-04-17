# -*- coding: utf-8

"""Module for testing components of the piping module.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_piping.py

SPDX-License-Identifier: MIT
"""
import numpy as np

from tespy.components.basics import sink
from tespy.components.basics import source
from tespy.components.piping import pipe
from tespy.components.piping import valve
from tespy.connections import connection
from tespy.networks.networks import network
from tespy.tools.characteristics import char_line
from tespy.tools.data_containers import dc_cc


def convergence_check(lin_dep):
    """Check convergence status of a simulation."""
    msg = 'Calculation did not converge!'
    assert lin_dep is False, msg


class TestPiping:

    def setup_piping_network(self, instance):
        self.nw = network(['CH4'], T_unit='C', p_unit='bar')
        self.source = source('source')
        self.sink = sink('sink')
        self.c1 = connection(self.source, 'out1', instance, 'in1')
        self.c2 = connection(instance, 'out1', self.sink, 'in1')
        self.nw.add_conns(self.c1, self.c2)

    def test_valve(self):
        """Test component properties of valves."""
        instance = valve('valve')
        self.setup_piping_network(instance)

        # parameter specification
        self.c1.set_attr(fluid={'CH4': 1}, m=10, p=10, T=120)
        self.c2.set_attr(p=1)

        # test variable pressure ration
        instance.set_attr(pr='var')
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        pr = round(self.c2.p.val_SI / self.c1.p.val_SI, 2)
        msg = ('Value of pressure ratio must be ' + str(pr) + ', is ' +
               str(round(instance.pr.val, 2)) + '.')
        assert pr == round(instance.pr.val, 2), msg

        # test variable zeta value
        zeta = round(instance.zeta.val, 0)
        instance.set_attr(zeta='var', pr=np.nan)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of dimension independent zeta value must be ' +
               str(zeta) + ', is ' + str(round(instance.zeta.val, 0)) + '.')
        assert zeta == round(instance.zeta.val, 0), msg

        # dp char
        x = np.array([8, 9, 10, 11, 12])
        y = np.array([5, 8, 9, 9.5, 9.6]) * 1e5
        dp_char = char_line(x, y)
        instance.set_attr(zeta=np.nan,
                          dp_char=dc_cc(func=dp_char, is_set=True))
        m = 11
        self.c1.set_attr(m=m)
        self.c2.set_attr(p=np.nan)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        self.nw.print_results()
        dp = round(-dp_char.evaluate(m), 0)
        dp_act = round(self.c2.p.val_SI - self.c1.p.val_SI)
        msg = ('The pressure drop at the valve should be ' + str(dp) + ' but '
               'is ' + str(dp_act) + '.')
        assert dp == dp_act, msg

    def test_pipe(self):
        """Test component properties of pipe."""
        instance = pipe('pipe')
        self.setup_piping_network(instance)

        # NO TEST NEEDED AT THE MOMENT, THE PIPE PROPERTIES ARE IDENTICAL TO
        # THE PROPERTIES OF THE SIMPLE HEAT EXCHANGER. TESTS ARE LOCATED AT
        # heat_exchanger_tests.py
