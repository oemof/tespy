# -*- coding: utf-8
# -*- coding: utf-8

from nose.tools import ok_, eq_

from tespy import nwk, cmp, con, hlp
from CoolProp.CoolProp import PropsSI as CP
import numpy as np
import shutil


class component_tests:

    def setup(self):
        self.nw = nwk.network(['INCOMP::DowQ', 'H2O', 'NH3', 'N2', 'O2', 'Ar'],
                              T_unit='C', p_unit='bar', v_unit='m3 / s')
        self.source = cmp.source('source')
        self.sink = cmp.sink('sink')

    def setup_network_11(self, instance):
        c1 = con.connection(self.source, 'out1', instance, 'in1')
        c2 = con.connection(instance, 'out1', self.sink, 'in1')
        self.nw.add_conns(c1, c2)
        return c1, c2

    def test_turbomachine(self):
        instance = cmp.turbomachine('turbomachine')
        c1, c2 = self.setup_network_11(instance)
        fl = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129, 'INCOMP::DowQ': 0, 'H2O': 0, 'NH3': 0}
        c1.set_attr(fluid=fl, m=10, p=1, h=1e5)
        c2.set_attr(p=1, h=2e5)
        self.nw.solve('design')
        power = c1.m.val_SI * (c2.h.val_SI - c1.h.val_SI)
        pr = c2.p.val_SI / c1.p.val_SI
        # pressure ratio and power are the basic functions for turbomachines, these are inherited by all children, thus only tested here
        eq_(power, instance.P.val, 'Value of power must be ' + str(power) + ', is ' + str(instance.P.val) + '.')
        eq_(pr, instance.pr.val, 'Value of power must be ' + str(pr) + ', is ' + str(instance.pr.val) + '.')
        c2.set_attr(p=np.nan)
        instance.set_attr(pr=5)
        self.nw.solve('design')
        pr = c2.p.val_SI / c1.p.val_SI
        eq_(pr, instance.pr.val, 'Value of power must be ' + str(pr) + ', is ' + str(instance.pr.val) + '.')
        c2.set_attr(h=np.nan)
        instance.set_attr(P=1e5)
        self.nw.solve('design')
        power = c1.m.val_SI * (c2.h.val_SI - c1.h.val_SI)
        eq_(pr, instance.pr.val, 'Value of power must be ' + str(pr) + ', is ' + str(instance.pr.val) + '.')
        instance.set_attr(eta_s=0.8)
        c2.set_attr(h=np.nan)
        try:
            self.nw.solve('design')
        except hlp.TESPyComponentError:
            pass
        try:
            instance.eta_s_deriv()
        except hlp.TESPyComponentError:
            pass

    def test_pump(self):
        instance = cmp.pump('pump')
        c1, c2 = self.setup_network_11(instance)
        fl = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 1, 'H2O': 0, 'NH3': 0}
        c1.set_attr(fluid=fl, v=1, p=5,T=50)
        c2.set_attr(p=7)
        instance.set_attr(eta_s=1)
        self.nw.solve('design')
        eta_s = (instance.h_os('') - c1.h.val_SI) / (c2.h.val_SI - c1.h.val_SI)
        eq_(eta_s, instance.eta_s.val, 'Value of isentropic efficiency must be ' + str(eta_s) + ', is ' + str(instance.eta_s.val) + '.')
        s1 = round(hlp.s_mix_ph(c1.to_flow()), 4)
        s2 = round(hlp.s_mix_ph(c2.to_flow()), 4)
        eq_(s1, s2, 'Value of entropy must be identical for inlet (' + str(s1) + ') and outlet (' + str(s2) + ') at 100 % isentropic efficiency.')
        instance.set_attr(eta_s=0.7)
        self.nw.solve('design')
        self.nw.save('tmp')
        c2.set_attr(p=np.nan)
        x = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4]
        y = np.array([14, 13.5, 12.5, 11, 9, 6.5, 3.5, 0]) * 1e5
        char = hlp.dc_cc(x=x, y=y, is_set=True)
        instance.set_attr(flow_char=char, eta_s=np.nan, eta_s_char=hlp.dc_cc(method='GENERIC', is_set=True))
        self.nw.solve('offdesign', design_path='tmp')
        eq_(c2.p.val_SI - c1.p.val_SI, 6.5e5, 'Value of pressure rise must be ' + str(6.5e5) + ', is ' + str(c2.p.val_SI - c1.p.val_SI) + '.')
        c1.set_attr(v=0.9)
        self.nw.solve('offdesign', design_path='tmp')
        eq_(c2.p.val_SI - c1.p.val_SI, 781250.0, 'Value of pressure rise must be ' + str(781250.0) + ', is ' + str(c2.p.val_SI - c1.p.val_SI) + '.')
        eq_(0.695, round(instance.eta_s.val, 3), 'Value of isentropic efficiency must be ' + str(0.695) + ', is ' + str(instance.eta_s.val) + '.')
        instance.eta_s_char.is_set = False
        c2.set_attr(T=con.ref(c1, 0, 20))
        c1.set_attr(v=-0.1)
        self.nw.solve('design')
        eq_(c2.p.val_SI - c1.p.val_SI, 14e5, 'Value of power must be ' + str(14e5) + ', is ' + str(c2.p.val_SI - c1.p.val_SI) + '.')
        c1.set_attr(v=1.5)
        self.nw.solve('design')
        eq_(c2.p.val_SI - c1.p.val_SI, 0, 'Value of power must be ' + str(0) + ', is ' + str(c2.p.val_SI - c1.p.val_SI) + '.')
        shutil.rmtree('./tmp', ignore_errors=True)

    def test_compressor(self):
        instance = cmp.compressor('compressor')
        c1, c2 = self.setup_network_11(instance)
        fl = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 0, 'H2O': 0, 'NH3': 1}
        c1.set_attr(fluid=fl, v=1, p=5, T=100)
        c2.set_attr(p=7)
        instance.set_attr(eta_s=0.8)
        self.nw.solve('design')
        self.nw.save('tmp')
        eta_s = (instance.h_os('') - c1.h.val_SI) / (c2.h.val_SI - c1.h.val_SI)
        eq_(round(eta_s, 3), round(instance.eta_s.val, 3), 'Value of isentropic efficiency must be ' + str(eta_s) + ', is ' + str(instance.eta_s.val) + '.')
        c2.set_attr(p=np.nan)
        instance.set_attr(char_map=hlp.dc_cm(method='GENERIC', is_set=True), eta_s=np.nan)
        self.nw.solve('offdesign', design_path='tmp')
        eq_(round(eta_s, 2), round(instance.eta_s.val, 2), 'Value of isentropic efficiency (' + str(instance.eta_s.val) + ') must be identical to design case (' + str(eta_s) + ').')
        # going above highes available speedline, beneath lowest mass flow at that line
        c1.set_attr(v=np.nan, m=c1.m.val*0.8, T=30)
        self.nw.solve('offdesign', design_path='tmp')
        eq_(round(eta_s * instance.char_map.z2[6, 0], 4), round(instance.eta_s.val, 4), 'Value of isentropic efficiency (' + str(instance.eta_s.val) + ') must be at (' + str(round(eta_s * instance.char_map.z2[6, 0], 4)) + ').')
        # going below lowest available speedline, above highest mass flow at that line
        c1.set_attr(T=300)
        self.nw.solve('offdesign', design_path='tmp')
        eq_(round(eta_s * instance.char_map.z2[0, 9], 4), round(instance.eta_s.val, 4), 'Value of isentropic efficiency (' + str(instance.eta_s.val) + ') must be at (' + str(round(eta_s * instance.char_map.z2[0, 9], 4)) + ').')
        c2.set_attr(p=7)
        c1.set_attr(v=1, T=100, m=np.nan)
        instance.set_attr(eta_s_char=hlp.dc_cc(method='GENERIC', is_set=True, param='m'))
        instance.char_map.is_set = False
        self.nw.solve('offdesign', design_path='tmp')
        eq_(round(eta_s, 3), round(instance.eta_s.val, 3), 'Value of isentropic efficiency must be ' + str(eta_s) + ', is ' + str(instance.eta_s.val) + '.')
        c1.set_attr(v=1.5)
        self.nw.solve('offdesign', design_path='tmp')
        eq_(0.82, round(instance.eta_s.val, 3), 'Value of isentropic efficiency must be ' + str(0.82) + ', is ' + str(instance.eta_s.val) + '.')
        instance.eta_s_char.set_attr(param='pr')
        c1.set_attr(v=1)
        c2.set_attr(p=7.5)
        self.nw.solve('offdesign', design_path='tmp')
        eq_(0.797, round(instance.eta_s.val, 3), 'Value of isentropic efficiency must be ' + str(0.797) + ', is ' + str(instance.eta_s.val) + '.')
        instance.eta_s_char.set_attr(param=None)
        try:
            self.nw.solve('offdesign', design_path='tmp')
        except ValueError:
            pass
        shutil.rmtree('./tmp', ignore_errors=True)

    def test_turbine(self):
        instance = cmp.turbine('turbine')
        c1, c2 = self.setup_network_11(instance)
        fl = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129, 'INCOMP::DowQ': 0, 'H2O': 0, 'NH3': 0}
        c1.set_attr(fluid=fl, m=15, p=10, T=120)
        c2.set_attr(p=1)
        instance.set_attr(eta_s=0.8)
        self.nw.solve('design')
        self.nw.save('tmp')
        eta_s = (c2.h.val_SI - c1.h.val_SI) / (instance.h_os('') - c1.h.val_SI)
        eq_(eta_s, instance.eta_s.val, 'Value of isentropic efficiency must be ' + str(eta_s) + ', is ' + str(instance.eta_s.val) + '.')
        c1.set_attr(p=np.nan)
        instance.cone.is_set = True
        instance.eta_s_char.is_set = True
        instance.eta_s.is_set = False
        self.nw.solve('offdesign', design_path='tmp')
        eq_(round(eta_s, 2), round(instance.eta_s.val, 2), 'Value of isentropic efficiency (' + str(instance.eta_s.val) + ') must be identical to design case (' + str(eta_s) + ').')
        # lowering mass flow, inlet pressure must sink according to cone law
        c1.set_attr(m=c1.m.val*0.8)
        self.nw.solve('offdesign', design_path='tmp')
        eq_(0.125, round(instance.pr.val, 3), 'Value of pressure ratio (' + str(instance.pr.val) + ') must be at (' + str(0.125) + ').')
        self.nw.print_results()
        shutil.rmtree('./tmp', ignore_errors=True)
