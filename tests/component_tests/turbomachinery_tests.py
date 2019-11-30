# -*- coding: utf-8

from nose.tools import eq_

from tespy.components.basics import sink, source
from tespy.components.turbomachinery import (turbomachine, compressor, pump,
                                             turbine)
from tespy.connections import connection, bus
from tespy.networks.networks import network
from tespy.tools.data_containers import dc_cc, dc_cm
from tespy.tools.fluid_properties import s_mix_ph
import numpy as np
import shutil


class turbomachinery_tests:

    def setup_network(self, instance):
        self.nw = network(['INCOMP::DowQ', 'H2O', 'NH3', 'N2',
                           'O2', 'Ar', 'CO2', 'CH4'],
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
        fl = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 0,
              'H2O': 0, 'NH3': 1, 'CO2': 0, 'CH4': 0}
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
        instance.set_attr(char_map=dc_cm(method='GENERIC', is_set=True),
                          eta_s=np.nan)

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
        eta_s = eta_s_d * instance.char_map.z2[6, 0]
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
                ') must be at (' + str(round(eta_s, 4)) + ').')
        eq_(round(eta_s, 4), round(instance.eta_s.val, 4), msg)

        # going below lowest available speedline, above highest mass flow at
        # that line
        self.c1.set_attr(T=300)
        self.nw.solve('offdesign', design_path='tmp')
        # should be value
        eta_s = eta_s_d * instance.char_map.z2[0, 9]
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
                ') must be at (' + str(round(eta_s, 4)) + ').')
        eq_(round(eta_s, 4), round(instance.eta_s.val, 4), msg)

        # back to design properties, test eta_s_char
        self.c2.set_attr(p=7)
        self.c1.set_attr(v=1, T=100, m=np.nan)

        # test parameter specification for eta_s_char, unset char map
        instance.set_attr(eta_s_char=dc_cc(method='GENERIC',
                                           is_set=True, param='m'))
        instance.char_map.is_set = False
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
                ', is ' + str(instance.eta_s.val) + '.')
        eq_(round(eta_s, 3), round(instance.eta_s.val, 3), msg)

        # move up in volumetric flow
        self.c1.set_attr(v=1.5)
        self.nw.solve('offdesign', design_path='tmp')
        eta_s = 0.88
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
                ', is ' + str(instance.eta_s.val) + '.')
        eq_(eta_s, round(instance.eta_s.val, 3), msg)

        # test parameter specification for pr
        instance.eta_s_char.set_attr(param='pr')
        self.c1.set_attr(v=1)
        self.c2.set_attr(p=7.5)
        self.nw.solve('offdesign', design_path='tmp')
        eta_s = 0.829
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
                ', is ' + str(instance.eta_s.val) + '.')
        eq_(eta_s, round(instance.eta_s.val, 3), msg)

        shutil.rmtree('./tmp', ignore_errors=True)

    # def test_pump(self):
    #     """
    #     Test component properties of pumps.
    #     """
    #     instance = pump('pump')
    #     c1, c2 = self.setup_network_11(instance)
    #     fl = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 1, 'H2O': 0,
    #           'NH3': 0, 'CO2': 0, 'CH4': 0}
    #     self.c1.set_attr(fluid=fl, v=1, p=5, T=50)
    #     self.c2.set_attr(p=7)
    #     instance.set_attr(eta_s=1)
    #     self.nw.solve('design')
    #     # calculate isentropic efficiency the old fashioned way
    #     eta_s = (instance.h_os('') - self.c1.h.val_SI) / (self.c2.h.val_SI - self.c1.h.val_SI)
    #     msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
    #             ', is ' + str(instance.eta_s.val) + '.')
    #     eq_(eta_s, instance.eta_s.val, msg)
    #     s1 = round(s_mix_ph(self.c1.to_flow()), 4)
    #     s2 = round(s_mix_ph(self.c2.to_flow()), 4)
    #     msg = ('Value of entropy must be identical for inlet (' + str(s1) +
    #             ') and outlet (' + str(s2) +
    #             ') at 100 % isentropic efficiency.')
    #     eq_(s1, s2, msg)
    #     instance.set_attr(eta_s=0.7)
    #     self.nw.solve('design')
    #     self.nw.save('tmp')
    #     self.c2.set_attr(p=np.nan)
    #     # flow char (pressure rise vs. volumetric flow)
    #     x = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4]
    #     y = np.array([14, 13.5, 12.5, 11, 9, 6.5, 3.5, 0]) * 1e5
    #     char = dc_cc(x=x, y=y, is_set=True)
    #     # apply flow char and eta_s char
    #     instance.set_attr(flow_char=char, eta_s=np.nan,
    #                       eta_s_char=dc_cc(method='GENERIC', is_set=True))
    #     self.nw.solve('offdesign', design_path='tmp')
    #     msg = ('Value of pressure rise must be ' + str(650000) + ', is ' +
    #             str(self.c2.p.val_SI - self.c1.p.val_SI) + '.')
    #     eq_(round(self.c2.p.val_SI - self.c1.p.val_SI, 0), 650000, msg)
    #     self.c1.set_attr(v=0.9)
    #     self.nw.solve('offdesign', design_path='tmp')
    #     msg = ('Value of pressure rise must be ' + str(775000.0) + ', is ' +
    #             str(self.c2.p.val_SI - self.c1.p.val_SI) + '.')
    #     eq_(self.c2.p.val_SI - self.c1.p.val_SI, 775000.0, msg)
    #     msg = ('Value of isentropic efficiency must be ' + str(0.694) +
    #             ', is ' + str(instance.eta_s.val) + '.')
    #     eq_(0.694, round(instance.eta_s.val, 3), msg)
    #     instance.eta_s_char.is_set = False
    #     # test boundaries of characteristic line
    #     self.c2.set_attr(T=ref(c1, 0, 20))
    #     self.c1.set_attr(v=-0.1)
    #     self.nw.solve('design')
    #     msg = ('Value of power must be ' + str(14e5) + ', is ' +
    #             str(self.c2.p.val_SI - self.c1.p.val_SI) + '.')
    #     eq_(self.c2.p.val_SI - self.c1.p.val_SI, 14e5, msg)
    #     self.c1.set_attr(v=1.5)
    #     self.nw.solve('design')
    #     msg = ('Value of power must be ' + str(0) + ', is ' +
    #             str(self.c2.p.val_SI - self.c1.p.val_SI) + '.')
    #     eq_(self.c2.p.val_SI - self.c1.p.val_SI, 0, msg)
    #     shutil.rmtree('./tmp', ignore_errors=True)

    # def test_turbine(self):
    #     """
    #     Test component properties of turbines.
    #     """
    #     instance = turbine('turbine')
    #     c1, c2 = self.setup_network_11(instance)
    #     fl = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129, 'INCOMP::DowQ': 0,
    #           'H2O': 0, 'NH3': 0, 'CO2': 0, 'CH4': 0}
    #     self.c1.set_attr(fluid=fl, m=15, p=10)
    #     self.c2.set_attr(p=1, T=20)
    #     instance.set_attr(eta_s=0.8)
    #     self.nw.solve('design')
    #     self.nw.save('tmp')
    #     # calculate isentropic efficiency the old fashioned way
    #     eta_s_d = ((self.c2.h.val_SI - self.c1.h.val_SI) /
    #                 (instance.h_os('') - self.c1.h.val_SI))
    #     msg = ('Value of isentropic efficiency must be ' + str(eta_s_d) +
    #             ', is ' + str(instance.eta_s.val) + '.')
    #     eq_(round(eta_s_d, 3), round(instance.eta_s.val, 3), msg)
    #     # trigger invalid isentropic efficiency
    #     instance.set_attr(eta_s=1.1)
    #     self.nw.solve('design')
    #     # calculate isentropic efficiency the old fashioned way
    #     eta_s = (self.c2.h.val_SI - self.c1.h.val_SI) / (instance.h_os('') - self.c1.h.val_SI)
    #     msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
    #             ', is ' + str(instance.eta_s.val) + '.')
    #     eq_(round(eta_s, 3), round(instance.eta_s.val, 3), msg)
    #     self.c1.set_attr(p=np.nan)
    #     instance.cone.is_set = True
    #     instance.eta_s_char.is_set = True
    #     instance.eta_s.is_set = False
    #     self.nw.solve('offdesign', design_path='tmp')
    #     msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
    #             ') must be identical to design case (' + str(eta_s_d) + ').')
    #     eq_(round(eta_s_d, 2), round(instance.eta_s.val, 2), msg)
    #     # lowering mass flow, inlet pressure must sink according to cone law
    #     self.c1.set_attr(m=self.c1.m.val * 0.8)
    #     self.nw.solve('offdesign', design_path='tmp')
    #     msg = ('Value of pressure ratio (' + str(instance.pr.val) +
    #             ') must be at (' + str(0.128) + ').')
    #     eq_(0.128, round(instance.pr.val, 3), msg)
    #     self.nw.print_results()
    #     # testing more parameters for eta_s_char
    #     self.c1.set_attr(m=10)
    #     # test param specification v
    #     instance.eta_s_char.param = 'v'
    #     self.nw.solve('offdesign', design_path='tmp')
    #     msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
    #             ') must be (' + str(0.8) + ').')
    #     eq_(0.8, round(instance.eta_s.val, 3), msg)
    #     # test param specification pr
    #     instance.eta_s_char.param = 'pr'
    #     self.nw.solve('offdesign', design_path='tmp')
    #     msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
    #             ') must be (' + str(0.768) + ').')
    #     eq_(0.768, round(instance.eta_s.val, 3), msg)
    #     # test param specification dh_s
    #     instance.eta_s_char.param = 'dh_s'
    #     self.nw.solve('offdesign', design_path='tmp')
    #     msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
    #             ') must be (' + str(0.798) + ').')
    #     eq_(0.798, round(instance.eta_s.val, 3), msg)
    #     instance.eta_s_char.param = None
    #     # test for missing parameter declaration
    #     try:
    #         self.nw.solve('offdesign', design_path='tmp')
    #     except ValueError:
    #         pass
    #     shutil.rmtree('./tmp', ignore_errors=True)

    # def test_turbomachine(self):
    #     """
    #     Test component properties of turbomachines.
    #     """
    #     instance = turbomachine('turbomachine')
    #     c1, c2 = self.setup_network_11(instance)
    #     fl = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129, 'INCOMP::DowQ': 0,
    #           'H2O': 0, 'NH3': 0, 'CO2': 0, 'CH4': 0}
    #     self.c1.set_attr(fluid=fl, m=10, p=1, h=1e5)
    #     self.c2.set_attr(p=1, h=2e5)
    #     self.nw.solve('design')
    #     power = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
    #     pr = self.c2.p.val_SI / self.c1.p.val_SI
    #     # pressure ratio and power are the basic functions for turbomachines,
    #     # these are inherited by all children, thus only tested here
    #     msg = ()
    #     eq_(power, instance.P.val, 'Value of power must be ' + str(power) +
    #         ', is ' + str(instance.P.val) + '.')
    #     msg = ()
    #     eq_(pr, instance.pr.val, 'Value of power must be ' + str(pr) +
    #         ', is ' + str(instance.pr.val) + '.')
    #     self.c2.set_attr(p=np.nan)
    #     instance.set_attr(pr=5)
    #     self.nw.solve('design')
    #     pr = self.c2.p.val_SI / self.c1.p.val_SI
    #     msg = ('Value of power must be ' + str(pr) + ', is ' +
    #             str(instance.pr.val) + '.')
    #     eq_(pr, instance.pr.val, msg)
    #     self.c2.set_attr(h=np.nan)
    #     instance.set_attr(P=1e5)
    #     self.nw.solve('design')
    #     power = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
    #     msg = ('Value of power must be ' + str(pr) + ', is ' +
    #             str(instance.pr.val) + '.')
    #     eq_(pr, instance.pr.val, msg)
    #     instance.set_attr(eta_s=0.8)
    #     self.c2.set_attr(h=np.nan)
