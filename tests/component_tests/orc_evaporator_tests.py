# -*- coding: utf-8

"""Module for testing components of type orc evaporator.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/component_tests/orc_evaporator_tests.py

SPDX-License-Identifier: MIT
"""

from nose.tools import eq_

from tespy.components.basics import sink, source
from tespy.components.customs import orc_evaporator
from tespy.connections import connection, bus
from tespy.networks.networks import network

import logging

import numpy as np
import shutil


class orc_evaporator_tests:

    def setup(self):

        self.nw = network(['water', 'Isopentane'], T_unit='C', p_unit='bar',
                          h_unit='kJ / kg')
        self.inl1 = source('inlet 1')
        self.outl1 = sink('outlet 1')

        self.inl2 = source('inlet 2')
        self.outl2 = sink('outlet 2')

        self.inl3 = source('inlet 3')
        self.outl3 = sink('outlet 3')

    def setup_orc_evaporator_network(self, instance):

        self.c1 = connection(self.inl1, 'out1', instance, 'in1')
        self.c2 = connection(instance, 'out1', self.outl1, 'in1')
        self.c3 = connection(self.inl2, 'out1', instance, 'in2')
        self.c4 = connection(instance, 'out2', self.outl2, 'in1')
        self.c5 = connection(self.inl3, 'out1', instance, 'in3')
        self.c6 = connection(instance, 'out3', self.outl3, 'in1')

        self.nw.add_conns(self.c1, self.c2, self.c3,
                          self.c4, self.c5, self.c6)

    def test_orc_evap(self):
        """
        Test component properties of orc evaporator.
        """
        instance = orc_evaporator('orc evaporator')
        self.setup_orc_evaporator_network(instance)

        # design specification
        instance.set_attr(pr1=0.93181818, pr2=0.970588, pr3=1,
                          design=['pr1', 'pr2', 'pr3'],
                          offdesign=['zeta1', 'zeta2', 'zeta3', 'kA'])
        self.c1.set_attr(T=146.6, p=4.34, m=20.4, state='g',
                         fluid={'water': 1, 'Isopentane': 0})
        self.c3.set_attr(T=146.6, p=10.2, m=190.8,
                         fluid={'water': 1, 'Isopentane': 0})
        self.c4.set_attr(T=118.6)
        self.c5.set_attr(T=111.6, p=10.8,
                         fluid={'water': 0, 'Isopentane': 1})
        self.nw.solve('design')
        self.nw.save('tmp')

        # Check the state of the steam and working fluid outlet:
        x_outl1 = self.c2.x.val
        x_outl3 = self.c6.x.val

        msg = ('Vapor mass fraction of steam outlet must be ' +
               str(round(x_outl1, 0)) + '.')
        eq_(round(x_outl1, 0), round(instance.x_outl1.val, 0), msg)

        msg = ('Vapor mass fraction of working fluid outlet must be ' +
               str(round(x_outl3, 1)) + '.')
        eq_(round(x_outl3, 1), round(instance.x_outl3.val, 1), msg)

        # test heat transfer
        Q = self.c5.m.val_SI * (self.c6.h.val_SI - self.c5.h.val_SI)
        msg = ('Value of heat flow be ' + str(round(instance.Q.val, 0)) +
               ', is ' + str(round(Q, 0)) + '.')
        eq_(round(Q, 1), round(instance.Q.val, 1), msg)

        shutil.rmtree('./tmp', ignore_errors=True)
