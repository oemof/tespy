# -*- coding: utf-8

from nose.tools import eq_

from tespy.components.basics import sink, source
from tespy.components.piping import valve, pipe
from tespy.connections import connection
from tespy.networks.networks import network

import numpy as np
import shutil


class piping_tests:

    def setup_piping_network(self, instance):
        self.nw = network(['CH4'], T_unit='C', p_unit='bar')
        self.source = source('source')
        self.sink = sink('sink')
        self.c1 = connection(self.source, 'out1', instance, 'in1')
        self.c2 = connection(instance, 'out1', self.sink, 'in1')
        self.nw.add_conns(self.c1, self.c2)

    def test_valve(self):
        """
        Test component properties of valves.
        """
        instance = valve('valve')
        self.setup_piping_network(instance)

        # parameter specification
        self.c1.set_attr(fluid={'CH4': 1}, m=10, p=10, T=120)
        self.c2.set_attr(p=1)

        # test variable pressure ration
        instance.set_attr(pr='var')
        self.nw.solve('design')
        pr = round(self.c2.p.val_SI / self.c1.p.val_SI, 2)
        msg = ('Value of pressure ratio must be ' + str(pr) + ', is ' +
               str(round(instance.pr.val, 2)) + '.')
        eq_(pr, round(instance.pr.val, 2), msg)

        # test variable zeta value
        zeta = round(instance.zeta.val, 0)
        instance.set_attr(zeta='var', pr=np.nan)
        self.nw.solve('design')
        msg = ('Value of dimension independent zeta value must be ' +
               str(zeta) + ', is ' + str(round(instance.zeta.val, 0)) + '.')
        eq_(zeta, round(instance.zeta.val, 0), msg)

    def test_pipe(self):
        """
        Test component properties of pipe.
        """
        instance = pipe('pipe')
        self.setup_piping_network(instance)

        # NO TEST NEEDED AT THE MOMENT, THE PIPE PROPERTIES ARE IDENTICAL TO
        # THE PROPERTIES OF THE SIMPLE HEAT EXCHANGER. TESTS ARE LOCATED AT
        # heat_exchanger_tests.py
