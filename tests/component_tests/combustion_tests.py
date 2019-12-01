# -*- coding: utf-8

from nose.tools import eq_

from tespy.components.basics import sink, source
from tespy.components.combustion import (combustion_chamber, combustion_engine)
from tespy.connections import connection, bus, ref
from tespy.networks.networks import network
from tespy.tools.data_containers import dc_cc, dc_cm
from tespy.tools.fluid_properties import s_mix_ph
import numpy as np
import shutil


class component_tests:

    def setup(self):

        self.nw = network(['H2O', 'N2', 'O2', 'Ar', 'CO2', 'CH4'],
                          T_unit='C', p_unit='bar', v_unit='m3 / s')
        self.fuel = source('fuel')
        self.air = source('ambient air')
        self.fg = sink('flue gas')

    def setup_combustion_chamber_network(self, instance):

        self.c1 = connection(self.air, 'out1', instance, 'in1')
        self.c2 = connection(self.fuel, 'out1', instance, 'in2')
        self.c3 = connection(instance, 'out1', self.fg, 'in1')
        self.nw.add_conns(c1, c2, c3)

    def setup_combustion_engine_network(self, instance):

        self.cw1_in = source('cooling water 1 source')
        self.cw2_in = source('cooling water 2 source')
        self.cw1_out = source('cooling water 1 sink')
        self.cw2_out = source('cooling water 2 sink')

        self.c1 = connection(self.air, 'out1', instance, 'in3')
        self.c2 = connection(self.fuel, 'out1', instance, 'in4')
        self.c3 = connection(instance, 'out3', self.fg, 'in1')
        self.c4 = connection(self.cw1_in, 'out1', instance, 'in1')
        self.c5 = connection(self.cw2_in, 'out1', instance, 'in2')
        self.c6 = connection(instance, 'out1', self.cw1_out, 'in1')
        self.c7 = connection(instance, 'out2', self.cw2_out, 'in1')
        self.nw.add_conns(c1, c2, c3, c4, c5, c6, c7)

    def test_combustion_chamber(self):
        """
        Test component properties of combustion chambers.
        """
        instance = combustion_chamber('combustion chamber')
        self.setup_combustion_chamber_network(instance)

        # connection parameter specification
        air = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129, 'H2O': 0, 'CO2': 0,
               'CH4': 0}
        fuel = {'N2': 0, 'O2': 0, 'Ar': 0, 'H2O': 0, 'CO2': 0.04, 'CH4': 0.96}
        self.c1.set_attr(fluid=air, p=1, T=30)
        self.c2.set_attr(fluid=fuel, T=30)
        self.c3.set_attr(T=1200)

        # test specified bus value on combustion_chamber (must be equal to ti)
        b = bus('thermal input', P=1e6)
        b.add_comps({'c': instance})
        self.nw.add_busses(b)
        self.nw.solve('design')
        msg = ('Value of thermal input must be ' + str(b.P.val) + ', is ' +
               str(instance.ti.val) + '.')
        eq_(round(b.P.val, 1), round(instance.ti.val, 1), msg)

        # test specified thermal input for combustion_chamber
        instance.set_attr(ti=1e6)
        self.nw.solve('design')
        ti = (self.c2.m.val_SI * self.c2.fluid.val['CH4'] *
              instance.fuels['CH4']['LHV'])
        msg = ('Value of thermal input must be ' + str(instance.ti.val) +
               ', is ' + str(ti) + '.')
        eq_(round(ti, 1), round(instance.ti.val, 1), msg)

        # test specified lamd for combustion_chamber

    def test_cogeneration_unit(self):
        """
        Test component properties of cogeneration unit.
        """
        amb = cmp.source('ambient')
        sf = cmp.source('fuel')
        fg = cmp.sink('flue gas outlet')
        cw_in1 = cmp.source('cooling water inlet1')
        cw_in2 = cmp.source('cooling water inlet2')
        cw_out1 = cmp.sink('cooling water outlet1')
        cw_out2 = cmp.sink('cooling water outlet2')
        chp = cmp.cogeneration_unit('cogeneration unit', fuel='CH4')
        amb_comb = connection(amb, 'out1', chp, 'in3')
        sf_comb = connection(sf, 'out1', chp, 'in4')
        comb_fg = connection(chp, 'out3', fg, 'in1')
        self.nw.add_conns(sf_comb, amb_comb, comb_fg)
        cw1_chp1 = connection(cw_in1, 'out1', chp, 'in1')
        cw2_chp2 = connection(cw_in2, 'out1', chp, 'in2')
        self.nw.add_conns(cw1_chp1, cw2_chp2)
        chp1_cw = connection(chp, 'out1', cw_out1, 'in1')
        chp2_cw = connection(chp, 'out2', cw_out2, 'in1')
        self.nw.add_conns(chp1_cw, chp2_cw)

        air = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129, 'INCOMP::DowQ': 0,
               'H2O': 0, 'NH3': 0, 'CO2': 0, 'CH4': 0}
        fuel = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 0, 'H2O': 0,
                'NH3': 0, 'CO2': 0.04, 'CH4': 0.96}
        water1 = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 0, 'H2O': 1,
                  'NH3': 0, 'CO2': 0, 'CH4': 0}
        water2 = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 0, 'H2O': 1,
                  'NH3': 0, 'CO2': 0, 'CH4': 0}

        chp.set_attr(fuel='CH4', pr1=0.99, pr2=0.99, lamb=1.2,
                     design=['pr1', 'pr2'], offdesign=['zeta1', 'zeta2'])
        amb_comb.set_attr(p=5, T=30, fluid=air)
        sf_comb.set_attr(T=30, fluid=fuel)
        cw1_chp1.set_attr(p=3, T=60, m=50, fluid=water1)
        cw2_chp2.set_attr(p=3, T=80, m=50, fluid=water2)

        TI = bus('thermal input')
        Q1 = bus('heat output 1')
        Q2 = bus('heat output 2')
        Q = bus('heat output')
        Qloss = bus('thermal heat loss')

        TI.add_comps({'c': chp, 'p': 'TI'})
        Q1.add_comps({'c': chp, 'p': 'Q1'})
        Q2.add_comps({'c': chp, 'p': 'Q2'})
        Q.add_comps({'c': chp, 'p': 'Q'})
        Qloss.add_comps({'c': chp, 'p': 'Qloss'})

        self.nw.add_busses(TI, Q1, Q2, Q, Qloss)
        ti = 1e6
        TI.set_attr(P=ti)

        # design point
        self.nw.solve('design')
        self.nw.save('tmp')
        self.nw.solve('offdesign', init_path='tmp', design_path='tmp')
        msg = ('Value of thermal input must be ' + str(TI.P.val) + ', is ' +
               str(chp.ti.val) + '.')
        eq_(round(TI.P.val, 1), round(chp.ti.val, 1), msg)
        # ti via component
        TI.set_attr(P=np.nan)
        chp.set_attr(ti=ti)
        self.nw.solve('offdesign', init_path='tmp', design_path='tmp')
        msg = ('Value of thermal input must be ' + str(ti) + ', is ' +
               str(chp.ti.val) + '.')
        eq_(round(ti, 1), round(chp.ti.val, 1), msg)
        chp.set_attr(ti=np.nan)
        # Q1 via bus
        Q1.set_attr(P=chp.Q1.val)
        self.nw.solve('offdesign', init_path='tmp', design_path='tmp')
        msg = ('Value of thermal input must be ' + str(ti) + ', is ' +
               str(chp.ti.val) + '.')
        eq_(round(ti, 1), round(chp.ti.val, 1), msg)
        heat1 = chp1_cw.m.val_SI * (chp1_cw.h.val_SI - cw1_chp1.h.val_SI)
        msg = ('Value of thermal input must be ' + str(heat1) + ', is ' +
               str(chp.Q1.val) + '.')
        eq_(round(heat1, 1), round(chp.Q1.val, 1), msg)
        Q1.set_attr(P=np.nan)
        # Q2 via bus
        Q2.set_attr(P=1.2 * chp.Q2.val)
        self.nw.solve('offdesign', init_path='tmp', design_path='tmp')
        # due to characteristic function Q1 is equal to Q2 for this
        # cogeneration unit
        heat1 = chp1_cw.m.val_SI * (chp1_cw.h.val_SI - cw1_chp1.h.val_SI)
        msg = ('Value of heat output 2 must be ' + str(heat1) + ', is ' +
               str(chp.Q2.val) + '.')
        eq_(round(heat1, 1), round(chp.Q2.val, 1), msg)
        # Q2 via component
        Q2.set_attr(P=np.nan)
        chp.set_attr(Q2=heat1)
        self.nw.solve('offdesign', init_path='tmp', design_path='tmp')
        heat1 = chp1_cw.m.val_SI * (chp1_cw.h.val_SI - cw1_chp1.h.val_SI)
        msg = ('Value of heat output 2 must be ' + str(heat1) + ', is ' +
               str(chp.Q2.val) + '.')
        eq_(round(heat1, 1), round(chp.Q2.val, 1), msg)
        # Q via bus
        chp.set_attr(Q2=np.nan)
        Q.set_attr(P=1.5 * chp.Q1.val)
        self.nw.solve('offdesign', init_path='tmp', design_path='tmp')
        msg = ('Value of total heat output (' + str(Q.P.val) +
               ') must be twice as much as value of heat output 2 (' +
               str(chp.Q2.val) + ').')
        eq_(round(Q.P.val, 1), round(2 * chp.Q2.val, 1), msg)
        # Qloss via bus
        Q.set_attr(P=np.nan)
        Qloss.set_attr(P=1e5)
        self.nw.solve('offdesign', init_path='tmp', design_path='tmp')
        msg = ('Value of heat loss must be ' + str(Qloss.P.val) + ', is ' +
               str(chp.Qloss.val) + '.')
        eq_(round(Qloss.P.val, 1), round(chp.Qloss.val, 1), msg)
        shutil.rmtree('./tmp', ignore_errors=True)
