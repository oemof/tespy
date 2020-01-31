# -*- coding: utf-8

from nose.tools import eq_

from tespy.components.basics import sink, source
from tespy.components.turbomachinery import (turbomachine, compressor, pump,
                                             turbine)
from tespy.connections import connection, bus, ref
from tespy.networks.networks import network
from tespy.tools.data_containers import dc_cc, dc_cm
from tespy.tools.fluid_properties import s_mix_ph
from tespy.tools.characteristics import char_line, compressor_map
from tespy.tools.characteristics import load_default_char as ldc
import numpy as np
import shutil


class turbomachinery_tests:

    def setup_network(self, instance):
        self.nw = network(['INCOMP::DowQ', 'NH3', 'N2', 'O2', 'Ar'],
                          T_unit='C', p_unit='bar', v_unit='m3 / s')
        self.source = source('source')
        self.sink = sink('sink')
        self.c1 = connection(self.source, 'out1', instance, 'in1')
        self.c2 = connection(instance, 'out1', self.sink, 'in1')
        self.nw.add_conns(self.c1, self.c2)

    def test_compressor(self):
        """
        Test component properties of compressors.
        """
        instance = compressor('compressor')
        self.setup_network(instance)

        # compress NH3, other fluids in network are for turbine, pump, ...
        fl = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 0, 'NH3': 1}
        self.c1.set_attr(fluid=fl, v=1, p=5, T=100)
        self.c2.set_attr(p=7)
        instance.set_attr(eta_s=0.8)
        self.nw.solve('design')
        self.nw.save('tmp')

        # test isentropic efficiency value
        eta_s_d = ((instance.h_os('') - self.c1.h.val_SI) /
                   (self.c2.h.val_SI - self.c1.h.val_SI))
        msg = ('Value of isentropic efficiency must be ' + str(eta_s_d) +
               ', is ' + str(instance.eta_s.val) + '.')
        eq_(round(eta_s_d, 3), round(instance.eta_s.val, 3), msg)

        # trigger invalid value for isentropic efficiency
        instance.set_attr(eta_s=1.1)
        self.nw.solve('design')

        # test calculated value
        eta_s = ((instance.h_os('') - self.c1.h.val_SI) /
                 (self.c2.h.val_SI - self.c1.h.val_SI))
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
               ', is ' + str(instance.eta_s.val) + '.')
        eq_(round(eta_s, 3), round(instance.eta_s.val, 3), msg)

        # remove pressure at outlet, use characteristic map for pressure
        # rise calculation
        self.c2.set_attr(p=np.nan)
        instance.set_attr(char_map=dc_cm(func=ldc(
            'compressor', 'char_map', 'DEFAULT', compressor_map),
            is_set=True), eta_s=np.nan)

        # offdesign test, efficiency value should be at design value
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
               ') must be identical to design case (' + str(eta_s) + ').')
        eq_(round(eta_s_d, 2), round(instance.eta_s.val, 2), msg)

        # move to highest available speedline, mass flow below lowest value
        # at that line
        self.c1.set_attr(v=np.nan, m=self.c1.m.val * 0.8, T=30)
        self.nw.solve('offdesign', design_path='tmp')
        # should be value
        eta_s = eta_s_d * instance.char_map.func.z2[6, 0]
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
               ') must be at (' + str(round(eta_s, 4)) + ').')
        eq_(round(eta_s, 4), round(instance.eta_s.val, 4), msg)

        # going below lowest available speedline, above highest mass flow at
        # that line
        self.c1.set_attr(T=300)
        self.nw.solve('offdesign', design_path='tmp', init_path='tmp')
        # should be value
        eta_s = eta_s_d * instance.char_map.func.z2[0, 9]
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
               ') must be at (' + str(round(eta_s, 4)) + ').')
        eq_(round(eta_s, 4), round(instance.eta_s.val, 4), msg)

        # back to design properties, test eta_s_char
        self.c2.set_attr(p=7)
        self.c1.set_attr(v=1, T=100, m=np.nan)

        # test parameter specification for eta_s_char with unset char map
        instance.set_attr(eta_s_char=dc_cc(func=ldc(
            'compressor', 'eta_s_char', 'DEFAULT', char_line),
            is_set=True))
        instance.char_map.is_set = False
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of isentropic efficiency must be ' + str(eta_s_d) +
               ', is ' + str(instance.eta_s.val) + '.')
        eq_(round(eta_s_d, 3), round(instance.eta_s.val, 3), msg)

        # move up in volumetric flow
        self.c1.set_attr(v=1.5)
        self.nw.solve('offdesign', design_path='tmp')
        eta_s = round(eta_s_d * instance.eta_s_char.func.evaluate(
            self.c1.m.val_SI / self.c1.m.design), 3)
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
               ', is ' + str(round(instance.eta_s.val, 3)) + '.')
        eq_(eta_s, round(instance.eta_s.val, 3), msg)

        # test parameter specification for pr
        instance.eta_s_char.set_attr(param='pr')
        self.c1.set_attr(v=1)
        self.c2.set_attr(p=7.5)
        self.nw.solve('offdesign', design_path='tmp')
        expr = (self.c2.p.val_SI * self.c1.p.design /
                (self.c2.p.design * self.c1.p.val_SI))
        eta_s = round(eta_s_d * instance.eta_s_char.func.evaluate(expr), 3)
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
               ', is ' + str(round(instance.eta_s.val, 3)) + '.')
        eq_(eta_s, round(instance.eta_s.val, 3), msg)

        shutil.rmtree('./tmp', ignore_errors=True)

    def test_pump(self):
        """
        Test component properties of pumps.
        """
        instance = pump('pump')
        self.setup_network(instance)
        fl = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 1, 'NH3': 0}
        self.c1.set_attr(fluid=fl, v=1, p=5, T=50)
        self.c2.set_attr(p=7)
        instance.set_attr(eta_s=1)
        self.nw.solve('design')

        # test calculated value for efficiency
        eta_s = ((instance.h_os('') - self.c1.h.val_SI) /
                 (self.c2.h.val_SI - self.c1.h.val_SI))
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
               ', is ' + str(instance.eta_s.val) + '.')
        eq_(eta_s, instance.eta_s.val, msg)

        # isentropic efficiency of 1 means inlet and outlet entropy are
        # identical
        s1 = round(s_mix_ph(self.c1.to_flow()), 4)
        s2 = round(s_mix_ph(self.c2.to_flow()), 4)
        msg = ('Value of entropy must be identical for inlet (' + str(s1) +
               ') and outlet (' + str(s2) +
               ') at 100 % isentropic efficiency.')
        eq_(s1, s2, msg)

        # specify realistic value for efficiency, outlet pressure from flow
        # char
        eta_s_d = 0.8
        instance.set_attr(eta_s=eta_s_d)
        self.nw.solve('design')
        self.nw.save('tmp')
        self.c2.set_attr(p=np.nan)

        # flow char (pressure rise vs. volumetric flow)
        x = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4]
        y = np.array([14, 13.5, 12.5, 11, 9, 6.5, 3.5, 0]) * 1e5
        char = dc_cc(func=char_line(x, y), is_set=True)
        # apply flow char and eta_s char
        instance.set_attr(flow_char=char, eta_s=np.nan,
                          eta_s_char=dc_cc(func=ldc('pump', 'eta_s_char',
                                                    'DEFAULT', char_line),
                                           is_set=True))
        self.nw.solve('offdesign', design_path='tmp')

        # value for difference pressure
        dp = 650000.0
        msg = ('Value of pressure rise must be ' + str(dp) + ', is ' +
               str(self.c2.p.val_SI - self.c1.p.val_SI) + '.')
        eq_(round(self.c2.p.val_SI - self.c1.p.val_SI, 0), dp, msg)

        # test ohter volumetric flow on flow char
        self.c1.set_attr(v=0.9)
        self.nw.solve('offdesign', design_path='tmp')
        dp = 775000.0
        msg = ('Value of pressure rise must be ' + str() + ', is ' +
               str(self.c2.p.val_SI - self.c1.p.val_SI) + '.')
        eq_(self.c2.p.val_SI - self.c1.p.val_SI, dp, msg)

        # test value of isentropic efficiency
        eta_s = round(eta_s_d * instance.eta_s_char.func.evaluate(
            self.c1.v.val_SI / self.c1.v.design), 3)
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
               ', is ' + str(instance.eta_s.val) + '.')
        eq_(eta_s, round(instance.eta_s.val, 3), msg)
        instance.eta_s_char.is_set = False

        # test boundaries of characteristic line:
        # lower boundary
        self.c2.set_attr(T=ref(self.c1, 0, 20))
        self.c1.set_attr(v=-0.1)
        self.nw.solve('design')
        msg = ('Value of power must be ' + str(14e5) + ', is ' +
               str(self.c2.p.val_SI - self.c1.p.val_SI) + '.')
        eq_(self.c2.p.val_SI - self.c1.p.val_SI, 14e5, msg)

        # upper boundary
        self.c1.set_attr(v=1.5)
        self.nw.solve('design')
        msg = ('Value of power must be ' + str(0) + ', is ' +
               str(self.c2.p.val_SI - self.c1.p.val_SI) + '.')
        eq_(self.c2.p.val_SI - self.c1.p.val_SI, 0, msg)
        shutil.rmtree('./tmp', ignore_errors=True)

    def test_turbine(self):
        """
        Test component properties of turbines.
        """
        instance = turbine('turbine')
        self.setup_network(instance)
        fl = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129, 'INCOMP::DowQ': 0,
              'NH3': 0}
        self.c1.set_attr(fluid=fl, m=15, p=10)
        self.c2.set_attr(p=1, T=20)
        instance.set_attr(eta_s=0.85)
        self.nw.solve('design')
        self.nw.save('tmp')

        # design value of isentropic efficiency
        eta_s_d = ((self.c2.h.val_SI - self.c1.h.val_SI) /
                   (instance.h_os('') - self.c1.h.val_SI))
        msg = ('Value of isentropic efficiency must be ' +
               str(round(eta_s_d, 3)) + ', is ' + str(instance.eta_s.val) +
               '.')
        eq_(round(eta_s_d, 3), round(instance.eta_s.val, 3), msg)

        # trigger invalid value for isentropic efficiency
        instance.set_attr(eta_s=1.1)
        self.nw.solve('design')
        eta_s = round((self.c2.h.val_SI - self.c1.h.val_SI) /
                      (instance.h_os('') - self.c1.h.val_SI), 3)
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
               ', is ' + str(instance.eta_s.val) + '.')
        eq_(eta_s, round(instance.eta_s.val, 3), msg)

        # unset isentropic efficiency and inlet pressure,
        # use characteristcs and cone law instead, parameters have to be in
        # design state
        self.c1.set_attr(p=np.nan)
        instance.cone.is_set = True
        instance.eta_s_char.is_set = True
        instance.eta_s.is_set = False
        self.nw.solve('offdesign', design_path='tmp')
        # check efficiency
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
               ') must be identical to design case (' + str(eta_s_d) + ').')
        eq_(round(eta_s_d, 2), round(instance.eta_s.val, 2), msg)
        # check pressure
        msg = ('Value of inlet pressure (' + str(round(self.c1.p.val_SI)) +
               ') must be identical to design case (' +
               str(round(self.c1.p.design)) + ').')
        eq_(round(self.c1.p.design), round(self.c1.p.val_SI), msg)

        # lowering mass flow, inlet pressure must sink according to cone law
        self.c1.set_attr(m=self.c1.m.val * 0.8)
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of pressure ratio (' + str(instance.pr.val) +
               ') must be at (' + str(0.128) + ').')
        eq_(0.128, round(instance.pr.val, 3), msg)

        # testing more parameters for eta_s_char
        # test parameter specification v
        self.c1.set_attr(m=10)
        instance.eta_s_char.param = 'v'
        self.nw.solve('offdesign', design_path='tmp')
        expr = self.c1.v.val_SI / self.c1.v.design
        eta_s = round(eta_s_d * instance.eta_s_char.func.evaluate(expr), 3)
        msg = ('Value of isentropic efficiency (' +
               str(round(instance.eta_s.val, 3)) +
               ') must be (' + str(eta_s) + ').')
        eq_(eta_s, round(instance.eta_s.val, 3), msg)

        # test parameter specification pr
        instance.eta_s_char.param = 'pr'
        self.nw.solve('offdesign', design_path='tmp')
        expr = (self.c2.p.val_SI * self.c1.p.design /
                (self.c2.p.design * self.c1.p.val_SI))
        eta_s = round(eta_s_d * instance.eta_s_char.func.evaluate(expr), 3)
        msg = ('Value of isentropic efficiency (' +
               str(round(instance.eta_s.val, 3)) +
               ') must be (' + str(eta_s) + ').')
        eq_(eta_s, round(instance.eta_s.val, 3), msg)

        # test parameter specification dh_s
        instance.eta_s_char.param = 'dh_s'
        self.nw.solve('offdesign', design_path='tmp')
        expr = (instance.h_os('') - self.c1.h.val_SI) / instance.dh_s_ref
        eta_s = round(eta_s_d * instance.eta_s_char.func.evaluate(expr), 3)
        msg = ('Value of isentropic efficiency (' +
               str(round(instance.eta_s.val, 3)) +
               ') must be (' + str(eta_s) + ').')
        eq_(eta_s, round(instance.eta_s.val, 3), msg)
        shutil.rmtree('./tmp', ignore_errors=True)

    def test_turbomachine(self):
        """
        Test component properties of turbomachines.
        """
        instance = turbomachine('turbomachine')
        self.setup_network(instance)
        fl = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129, 'INCOMP::DowQ': 0,
              'H2O': 0, 'NH3': 0, 'CO2': 0, 'CH4': 0}
        self.c1.set_attr(fluid=fl, m=10, p=1, h=1e5)
        self.c2.set_attr(p=3, h=2e5)

        # pressure ratio and power are the basic functions for turbomachines,
        # these are inherited by all children, thus only tested here
        self.nw.solve('design')
        power = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        pr = self.c2.p.val_SI / self.c1.p.val_SI
        msg = ('Value of power must be ' + str(power) + ', is ' +
               str(instance.P.val) + '.')
        eq_(power, instance.P.val, msg)
        msg = ('Value of power must be ' + str(pr) + ', is ' +
               str(instance.pr.val) + '.')
        eq_(pr, instance.pr.val, )

        # test pressure ratio specification
        self.c2.set_attr(p=np.nan)
        instance.set_attr(pr=5)
        self.nw.solve('design')
        pr = self.c2.p.val_SI / self.c1.p.val_SI
        msg = ('Value of power must be ' + str(pr) + ', is ' +
               str(instance.pr.val) + '.')
        eq_(pr, instance.pr.val, msg)

        # test power specification
        self.c2.set_attr(h=np.nan)
        instance.set_attr(P=1e5)
        self.nw.solve('design')
        power = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        msg = ('Value of power must be ' + str(pr) + ', is ' +
               str(instance.pr.val) + '.')
        eq_(pr, instance.pr.val, msg)
        instance.set_attr(eta_s=0.8)
        self.c2.set_attr(h=np.nan)
