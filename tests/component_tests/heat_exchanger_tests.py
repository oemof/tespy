# -*- coding: utf-8

"""Module for testing components of type heat exchanger.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/component_tests/heat_exchanger_tests.py

SPDX-License-Identifier: MIT
"""

from nose.tools import eq_

from tespy.components.basics import sink, source
from tespy.components.heat_exchangers import (heat_exchanger_simple,
                                              solar_collector,
                                              heat_exchanger, condenser)
from tespy.connections import connection, bus
from tespy.networks.networks import network
from tespy.tools.fluid_properties import T_bp_p

import logging

import numpy as np
import shutil


class heat_exchanger_tests:

    def setup(self):

        self.nw = network(['H2O', 'Ar'], T_unit='C', p_unit='bar',
                          v_unit='m3 / s')
        self.inl1 = source('inlet 1')
        self.outl1 = sink('outlet 1')

    def setup_heat_exchanger_simple_network(self, instance):

        self.c1 = connection(self.inl1, 'out1', instance, 'in1')
        self.c2 = connection(instance, 'out1', self.outl1, 'in1')

        self.nw.add_conns(self.c1, self.c2)

    def setup_heat_exchanger_network(self, instance):

        self.inl2 = source('inlet 2')
        self.outl2 = sink('outlet 2')

        self.c1 = connection(self.inl1, 'out1', instance, 'in1')
        self.c2 = connection(instance, 'out1', self.outl1, 'in1')
        self.c3 = connection(self.inl2, 'out1', instance, 'in2')
        self.c4 = connection(instance, 'out2', self.outl2, 'in1')

        self.nw.add_conns(self.c1, self.c2, self.c3, self.c4)

    def test_heat_ex_simple(self):
        """
        Test component properties of simple heat exchanger.
        """
        instance = heat_exchanger_simple('heat exchanger')
        self.setup_heat_exchanger_simple_network(instance)
        fl = {'Ar': 0, 'H2O': 1}
        self.c1.set_attr(fluid=fl, m=1, p=10, T=100)
        # trigger heat exchanger parameter groups
        instance.set_attr(hydro_group='HW', L=100, ks=100, pr=0.99, Tamb=20)

        # test grouped parameter settings with missing parameters
        instance.hydro_group.is_set = True
        instance.kA_group.is_set = True
        self.nw.solve('design', init_only=True)
        msg = ('Hydro group must no be set, if one parameter is missing!')
        eq_(instance.hydro_group.is_set, False, msg)
        msg = ('kA group must no be set, if one parameter is missing!')
        eq_(instance.kA_group.is_set, False, msg)

        # test diameter calculation from specified dimensions (as pipe)
        # with Hazen-Williams method
        instance.set_attr(hydro_group='HW', D='var', L=100,
                          ks=100, pr=0.99, Tamb=20)
        b = bus('heat', P=-1e5)
        b.add_comps({'c': instance})
        self.nw.add_busses(b)
        self.nw.solve('design')
        pr = round(self.c2.p.val_SI / self.c1.p.val_SI, 3)
        msg = ('Value of pressure ratio must be ' + str(pr) + ', is ' +
               str(instance.pr.val) + '.')
        eq_(pr, round(instance.pr.val, 3), msg)

        # make zeta system variable and use previously calculated diameter
        # to calculate zeta. The value for zeta must not change
        zeta = round(instance.zeta.val, 0)
        instance.set_attr(D=instance.D.val, zeta='var', pr=np.nan)
        instance.D.is_var = False
        self.nw.solve('design')
        msg = ('Value of zeta must be ' + str(zeta) + ', is ' +
               str(round(instance.zeta.val, 0)) + '.')
        eq_(zeta, round(instance.zeta.val, 0), msg)

        # same test with pressure ratio as sytem variable
        pr = round(instance.pr.val, 3)
        instance.set_attr(zeta=np.nan, pr='var')
        self.nw.solve('design')
        msg = ('Value of pressure ratio must be ' + str(pr) +
               ', is ' + str(round(instance.pr.val, 3)) + '.')
        eq_(pr, round(instance.pr.val, 3), msg)

        # test heat transfer coefficient as variable of the system (ambient
        # temperature required)
        instance.set_attr(kA='var', pr=np.nan)
        b.set_attr(P=-5e4)
        self.nw.solve('design')

        # due to heat output being half of reference (for Tamb) kA should be
        # somewhere near to that (actual value is 677)
        msg = ('Value of heat transfer coefficient must be 667, is ' +
               str(instance.kA.val) + '.')
        eq_(677, round(instance.kA.val, 0), msg)

        # test heat transfer as variable of the system
        instance.set_attr(Q='var', kA=np.nan)
        Q = -5e4
        b.set_attr(P=Q)
        self.nw.solve('design')
        msg = ('Value of heat transfer must be ' + str(Q) +
               ', is ' + str(instance.Q.val) + '.')
        eq_(Q, round(instance.Q.val, 0), msg)

    def test_solar_collector(self):
        """
        Test component properties of solar collector.
        """
        instance = solar_collector('solar collector')
        self.setup_heat_exchanger_simple_network(instance)
        fl = {'Ar': 0, 'H2O': 1}
        self.c1.set_attr(fluid=fl, p=10, T=30)
        self.c2.set_attr(T=70)

        # test grouped parameter settings with missing parameters
        instance.hydro_group.is_set = True
        instance.energy_group.is_set = True
        self.nw.solve('design', init_only=True)
        msg = ('Hydro group must no be set, if one parameter is missing!')
        eq_(instance.hydro_group.is_set, False, msg)
        msg = ('Energy group must no be set, if one parameter is missing!')
        eq_(instance.energy_group.is_set, False, msg)

        # test solar collector params as system variables
        instance.set_attr(E=1e3, lkf_lin=1.0, lkf_quad=0.005, A='var',
                          eta_opt=0.9, Q=1e5, Tamb=20, pr=0.99)
        self.nw.solve('design')
        # heat loss must be identical to Q - E * A (internal heat loss
        # calculation)
        T_diff = (self.c2.T.val + self.c1.T.val) / 2 - instance.Tamb.val
        Q_loss = round(instance.A.val * (
            instance.E.val * (1 - instance.eta_opt.val) +
            T_diff * instance.lkf_lin.val +
            T_diff ** 2 * instance.lkf_quad.val), 0)
        msg = ('Value for heat loss of solar collector must be '
               + str(Q_loss) + ', is ' + str(round(instance.Q_loss.val, 0)) +
               '.')
        eq_(Q_loss, round(instance.Q_loss.val, 0), msg)

        # test all parameters of the energy group: E
        instance.set_attr(A=instance.A.val, E='var')
        self.nw.solve('design')
        eq_(Q_loss, round(instance.Q_loss.val, 0), msg)

        # test all parameters of the energy group: eta_opt
        instance.set_attr(E=instance.E.val, eta_opt='var')
        self.nw.solve('design')
        eq_(Q_loss, round(instance.Q_loss.val, 0), msg)

        # test all parameters of the energy group: lkf_lin
        instance.set_attr(eta_opt=instance.eta_opt.val, lkf_lin='var')
        self.nw.solve('design')
        eq_(Q_loss, round(instance.Q_loss.val, 0), msg)

        # test all parameters of the energy group: lkf_quad
        instance.set_attr(lkf_lin=instance.lkf_lin.val, lkf_quad='var')
        self.nw.solve('design')
        eq_(Q_loss, round(instance.Q_loss.val, 0), msg)

        # test all parameters of the energy group: Tamb
        instance.set_attr(lkf_lin=instance.lkf_lin.val, lkf_quad='var')
        self.nw.solve('design')
        eq_(Q_loss, round(instance.Q_loss.val, 0), msg)

    def test_heat_ex(self):
        """
        Test component properties of heat exchanger.
        """
        instance = heat_exchanger('heat exchanger')
        self.setup_heat_exchanger_network(instance)

        # design specification
        instance.set_attr(pr1=0.98, pr2=0.98, ttd_u=5,
                          design=['pr1', 'pr2', 'ttd_u'],
                          offdesign=['zeta1', 'zeta2', 'kA'])
        self.c1.set_attr(T=120, p=3, fluid={'Ar': 0, 'H2O': 1})
        self.c2.set_attr(T=70)
        self.c3.set_attr(T=40, p=5, fluid={'Ar': 1, 'H2O': 0})
        b = bus('heat transfer', P=-80e3)
        b.add_comps({'c': instance})
        self.nw.add_busses(b)
        self.nw.solve('design')
        self.nw.save('tmp')

        # check heat transfer
        Q = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        td_log = ((self.c2.T.val - self.c3.T.val -
                   self.c1.T.val + self.c4.T.val) /
                  np.log((self.c2.T.val - self.c3.T.val) /
                         (self.c1.T.val - self.c4.T.val)))
        kA = round(-Q / td_log, 0)
        msg = ('Value of heat transfer must be ' + str(round(Q, 0)) + ', is ' +
               str(round(instance.Q.val, 0)) + '.')
        eq_(round(Q, 0), round(instance.Q.val, 0), msg)

        # check upper terminal temperature difference
        msg = ('Value of terminal temperature difference must be ' +
               str(round(instance.ttd_u.val, 1)) + ', is ' +
               str(round(self.c1.T.val - self.c4.T.val, 1)) + '.')
        eq_(round(instance.ttd_u.val, 1),
            round(self.c1.T.val - self.c4.T.val, 1), msg)

        # check lower terminal temperature difference
        self.c2.set_attr(T=np.nan)
        instance.set_attr(ttd_l=20)
        self.nw.solve('design')
        msg = ('Value of terminal temperature difference must be ' +
               str(instance.ttd_l.val) + ', is ' +
               str(self.c2.T.val - self.c3.T.val) + '.')
        eq_(round(self.c2.T.val - self.c3.T.val, 1),
            round(instance.ttd_l.val, 1), msg)

        # check specified kA value (by offdesign parameter), reset temperatures
        # to design state
        self.c2.set_attr(T=70)
        instance.set_attr(ttd_l=np.nan)
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of heat flow must be ' + str(instance.Q.val) + ', is ' +
               str(round(Q, 0)) + '.')
        eq_(round(Q, 0), round(instance.Q.val, 0), msg)
        msg = ('Value of heat transfer coefficient must be ' + str(kA) +
               ', is ' + str(round(instance.kA.val, 0)) + '.')
        eq_(kA, round(instance.kA.val, 0), msg)

        # trigger negative lower terminal temperature difference as result
        self.c4.set_attr(T=np.nan)
        self.c2.set_attr(T=30)
        self.nw.solve('design')
        msg = ('Value of upper terminal temperature differences must be '
               'smaller than zero, is ' + str(round(instance.ttd_l.val, 1)) +
               '.')
        eq_(True, instance.ttd_l.val < 0, msg)

        # trigger negative upper terminal temperature difference as result
        self.c4.set_attr(T=100)
        self.c2.set_attr(h=200e3, T=np.nan)
        instance.set_attr(pr1=0.98, pr2=0.98, ttd_u=np.nan,
                          design=['pr1', 'pr2'])
        self.c1.set_attr(h=150e3, T=np.nan)
        self.c3.set_attr(T=40)
        self.nw.solve('design')
        msg = ('Value of upper terminal temperature differences must be '
               'smaller than zero, is ' + str(round(instance.ttd_u.val, 1)) +
               '.')
        eq_(True, instance.ttd_u.val < 0, msg)

        shutil.rmtree('./tmp', ignore_errors=True)

    def test_condenser(self):
        """
        Test component properties of condenser.
        """
        instance = condenser('condenser')
        self.setup_heat_exchanger_network(instance)

        # design specification
        instance.set_attr(pr1=0.98, pr2=0.98, ttd_u=5,
                          offdesign=['zeta2', 'kA'])
        self.c1.set_attr(T=100, p0=0.5, fluid={'Ar': 0, 'H2O': 1})
        self.c3.set_attr(T=30, p=5, fluid={'Ar': 0, 'H2O': 1})
        self.c4.set_attr(T=40)
        instance.set_attr(Q=-80e3)
        self.nw.solve('design')
        self.nw.save('tmp')

        # test heat transfer
        Q = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        msg = ('Value ofheat flow be ' + str(round(instance.Q.val, 0)) +
               ', is ' + str(round(Q, 0)) + '.')
        eq_(round(Q, 1), round(instance.Q.val, 1), msg)

        # test upper terminal temperature difference. For the component
        # condenser the temperature of the condensing fluid is relevant.
        ttd_u = round(T_bp_p(self.c1.to_flow()) - self.c4.T.val_SI, 1)
        p = round(self.c1.p.val_SI, 5)
        kA = instance.kA.val
        msg = ('Value of terminal temperature difference must be ' +
               str(round(instance.ttd_u.val, 1)) + ', is ' +
               str(ttd_u) + '.')
        eq_(ttd_u, round(instance.ttd_u.val, 1), msg)

        # test lower terminal temperature difference
        instance.set_attr(ttd_l=20, ttd_u=np.nan, design=['pr2', 'ttd_l'])
        self.nw.solve('design')
        msg = ('Value of terminal temperature difference must be ' +
               str(instance.ttd_l.val) + ', is ' +
               str(self.c2.T.val - self.c3.T.val) + '.')
        eq_(round(self.c2.T.val - self.c3.T.val, 1),
            round(instance.ttd_l.val, 1), msg)

        # check kA value with condensing pressure in offdesign mode:
        # no changes to design point means: identical pressure
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of condensing pressure be ' + str(p) + ', is ' +
               str(round(self.c1.p.val_SI, 5)) + '.')
        eq_(p, round(self.c1.p.val_SI, 5), msg)
        shutil.rmtree('./tmp', ignore_errors=True)
