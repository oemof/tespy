# -*- coding: utf-8

from nose.tools import eq_

from tespy.components.basics import sink, source
from tespy.components.heat_exchangers import (heat_exchanger_simple,
                                              solar_collector,
                                              heat_exchanger, condenser)
from tespy.connections import connection, bus
from tespy.networks.networks import network

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
        self.c1.set_attr(fluid=fl, m=1, p=10, T=100)

        # trigger solar collector parameter groups
        instance.set_attr(hydro_group='default', L=100, ks=100, pr=0.99,
                          Tamb=20)
        # test grouped parameter settings with missing parameters
        instance.hydro_group.is_set = True
        instance.energy_group.is_set = True
        self.nw.solve('design', init_only=True)
        msg = ('Hydro group must no be set, if one parameter is missing!')
        eq_(instance.hydro_group.is_set, False, msg)
        msg = ('Energy group must no be set, if one parameter is missing!')
        eq_(instance.energy_group.is_set, False, msg)

        # test solar collector params as system variables




    def test_heat_ex(self):
        """
        Test component properties of heat exchanger.
        """
        tesin = sink('TES in')
        tesout = source('TES out')
        hsin = sink('HS in')
        hsout = source('HS out')
        he = heat_exchanger('heat exchanger')
        tes_he = connection(tesout, 'out1', he, 'in2')
        he_tes = connection(he, 'out2', tesin, 'in1')
        hs_he = connection(hsout, 'out1', he, 'in1')
        he_hs = connection(he, 'out1', hsin, 'in1')
        self.nw.add_conns(tes_he, he_tes, hs_he, he_hs)
        # design specification
        he.set_attr(pr1=0.98, pr2=0.98, ttd_u=5,
                    design=['pr1', 'pr2', 'ttd_u'],
                    offdesign=['zeta1', 'zeta2', 'kA'])
        hs_he.set_attr(T=120, p=3, fluid={'Ar': 0, 'H2O': 1})
        he_hs.set_attr(T=70)
        tes_he.set_attr(T=40, p=5, fluid={'Ar': 1, 'H2O': 0})
        b = bus('heat transfer', P=-80e3)
        b.add_comps({'c': he})
        self.nw.add_busses(b)
        self.nw.solve('design')
        # check heat flow
        Q = hs_he.m.val_SI * (he_hs.h.val_SI - hs_he.h.val_SI)
        self.nw.save('tmp')
        msg = ('Value of terminal temperature difference must be ' +
               str(he.ttd_u.val) + ', is ' + str(hs_he.T.val - he_tes.T.val) +
               '.')
        eq_(round(hs_he.T.val - he_tes.T.val, 1), round(he.ttd_u.val, 1), msg)
        # check lower terminal temperature difference
        he_hs.set_attr(T=np.nan)
        he.set_attr(ttd_l=20)
        self.nw.solve('design')
        msg = ('Value of terminal temperature difference must be ' +
               str(he.ttd_l.val) + ', is ' + str(he_hs.T.val - tes_he.T.val) +
               '.')
        eq_(round(he_hs.T.val - tes_he.T.val, 1), round(he.ttd_l.val, 1), msg)
        # check kA value
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of heat flow must be ' + str(he.Q.val) +
               ', is ' + str(Q) + '.')
        eq_(round(Q, 1), round(he.Q.val, 1), msg)
        # trigger errors for negative terminal temperature differences at given
        # kA-value
        he.set_attr(ttd_l=np.nan)
        # ttd_l
        he_hs.set_attr(T=30)
        try:
            self.nw.solve('offdesign', design_path='tmp')
        except ValueError:
            pass
        # ttd_u
        he_hs.set_attr(T=np.nan)
        he_tes.set_attr(T=130)
        try:
            self.nw.solve('offdesign', design_path='tmp')
        except ValueError:
            pass

        # trigger negative lower terminal temperature difference as result
        he_tes.set_attr(T=np.nan)
        he_hs.set_attr(T=30)
        self.nw.solve('design')
        msg = ('Value of upper terminal temperature differences must be '
               'smaller than zero, is ' + str(round(he.ttd_l.val, 1)) + '.')
        eq_(True, he.ttd_l.val < 0, msg)

        # trigger negative upper terminal temperature difference as result
        he_tes.set_attr(T=100)
        hs_he.set_attr(h=200e3, T=np.nan)
        he.set_attr(pr1=0.98, pr2=0.98, design=['pr1', 'pr2'], ttd_u=np.nan)
        he_hs.set_attr(h=150e3, T=np.nan)
        tes_he.set_attr(T=40)
        self.nw.solve('design')
        msg = ('Value of upper terminal temperature differences must be '
               'smaller than zero, is ' + str(round(he.ttd_u.val, 1)) + '.')
        eq_(True, he.ttd_u.val < 0, msg)

        shutil.rmtree('./tmp', ignore_errors=True)

    def test_condenser(self):
        """
        Test component properties of condenser.
        """
        tesin = sink('TES in')
        tesout = source('TES out')
        hsin = sink('Cond in')
        hsout = source('Cond out')
        he = condenser('condenser')
        tes_he = connection(tesout, 'out1', he, 'in2')
        he_tes = connection(he, 'out2', tesin, 'in1')
        hs_he = connection(hsout, 'out1', he, 'in1')
        he_hs = connection(he, 'out1', hsin, 'in1')
        self.nw.add_conns(tes_he, he_tes, hs_he, he_hs)
        # design specification
        he.set_attr(pr1=0.98, pr2=0.98, ttd_u=5, design=['pr2', 'ttd_u'],
                    offdesign=['zeta2', 'kA'])
        hs_he.set_attr(T=100, p0=0.5, fluid={'Ar': 0, 'H2O': 1})
        tes_he.set_attr(T=30, p=5, fluid={'Ar': 0, 'H2O': 1})
        he_tes.set_attr(T=40)
        he.set_attr(Q=-80e3)
        self.nw.solve('design')
        # check heat flow
        Q = hs_he.m.val_SI * (he_hs.h.val_SI - hs_he.h.val_SI)
        msg = ('Value ofheat flow be ' +
               str(he.Q.val) + ', is ' + str(Q) + '.')
        eq_(round(Q, 1), round(he.Q.val, 1), msg)
        self.nw.save('tmp')
        ttd_u = hlp.T_bp_p(hs_he.to_flow()) - he_tes.T.val_SI
        p = hs_he.p.val_SI
        msg = ('Value of terminal temperature difference must be ' +
               str(he.ttd_u.val) + ', is ' + str(ttd_u) + '.')
        eq_(round(ttd_u, 1), round(he.ttd_u.val, 1), msg)
        # check lower terminal temperature difference
        he.set_attr(ttd_l=20, ttd_u=np.nan, design=['pr2', 'ttd_l'])
        self.nw.solve('design')
        msg = ('Value of terminal temperature difference must be ' +
               str(he.ttd_l.val) + ', is ' + str(he_hs.T.val - tes_he.T.val) +
               '.')
        eq_(round(he_hs.T.val - tes_he.T.val, 1), round(he.ttd_l.val, 1), msg)
        # check kA value
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of condensing pressure be ' + str(p) +
               ', is ' + str(hs_he.p.val_SI) + '.')
        eq_(round(p, 1), round(hs_he.p.val_SI, 1), msg)
        shutil.rmtree('./tmp', ignore_errors=True)
