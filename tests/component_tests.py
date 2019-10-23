# -*- coding: utf-8

from nose.tools import eq_

from tespy import nwk, cmp, con, hlp
import numpy as np
import shutil


class component_tests:

    def setup(self):
        self.nw = nwk.network(['INCOMP::DowQ', 'H2O', 'NH3', 'N2',
                               'O2', 'Ar', 'CO2', 'CH4'],
                              T_unit='C', p_unit='bar', v_unit='m3 / s')
        self.source = cmp.source('source')
        self.sink = cmp.sink('sink')

    def setup_network_11(self, instance):
        """
        Set up network for components with one inlet and one outlet.
        """
        c1 = con.connection(self.source, 'out1', instance, 'in1')
        c2 = con.connection(instance, 'out1', self.sink, 'in1')
        self.nw.add_conns(c1, c2)
        return c1, c2

    def setup_network_21(self, instance):
        """
        Set up network for components with two inlets and one outlet.
        """
        c1 = con.connection(self.source, 'out1', instance, 'in1')
        c2 = con.connection(cmp.source('fuel'), 'out1', instance, 'in2')
        c3 = con.connection(instance, 'out1', self.sink, 'in1')
        self.nw.add_conns(c1, c2, c3)
        return c1, c2, c3

    def setup_network_electrolyzer(self, instance):
        """
        Set up network for electrolyzer tests.
        """
        fw = cmp.source('feed water')
        cw_in = cmp.source('cooling water')
        o2 = cmp.sink('oxygen sink')
        h2 = cmp.sink('hydrogen sink')
        cw_out = cmp.sink('cooling water sink')

        instance.set_attr(pr_c=0.99)

        cw_el = con.connection(cw_in, 'out1', instance, 'in1',
                               fluid={'H2O': 1, 'H2': 0, 'O2': 0}, T=20, p=1)
        el_cw = con.connection(instance, 'out1', cw_out, 'in1', T=45)

        self.nw.add_conns(cw_el, el_cw)

        fw_el = con.connection(fw, 'out1', instance, 'in2', m=0.1, T=20, p=10)
        el_o2 = con.connection(instance, 'out2', o2, 'in1')
        el_h2 = con.connection(instance, 'out3', h2, 'in1', T=50)

        self.nw.add_conns(fw_el, el_o2, el_h2)

    def setup_network_individual_offdesign(self):
        """
        Set up network for individual offdesign tests.
        """
        so = cmp.source('source')
        sp = cmp.splitter('splitter', num_out=2)
        self.pump1 = cmp.pump('pump 1')
        self.sc1 = cmp.solar_collector('collector field 1')
        v1 = cmp.valve('valve1')
        self.pump2 = cmp.pump('pump 2')
        self.sc2 = cmp.solar_collector('collector field 2')
        v2 = cmp.valve('valve2')
        me = cmp.merge('merge', num_in=2)
        si = cmp.sink('sink')

        self.pump1.set_attr(eta_s=0.8, design=['eta_s'],
                            offdesign=['eta_s_char'])
        self.pump2.set_attr(eta_s=0.8, design=['eta_s'],
                            offdesign=['eta_s_char'])
        self.sc1.set_attr(pr=0.95, lkf_lin=3.33, lkf_quad=0.011, A=1252, E=700,
                          Tamb=20, design=['pr'], offdesign=['zeta'])
        self.sc2.set_attr(pr=0.95, lkf_lin=3.5, lkf_quad=0.011, A=700, E=800,
                          Tamb=20, design=['pr'], offdesign=['zeta'])

        fl = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 0,
              'H2O': 1, 'NH3': 0, 'CO2': 0, 'CH4': 0}

        inlet = con.connection(so, 'out1', sp, 'in1', T=50, p=3, fluid=fl)
        outlet = con.connection(me, 'out1', si, 'in1', p=3)

        sp_p1 = con.connection(sp, 'out1', self.pump1, 'in1')
        p1_sc1 = con.connection(self.pump1, 'out1', self.sc1, 'in1')
        self.sc1_v1 = con.connection(self.sc1, 'out1', v1, 'in1', p=3.1, T=90)
        v1_me = con.connection(v1, 'out1', me, 'in1')

        self.sp_p2 = con.connection(sp, 'out2', self.pump2, 'in1')
        self.p2_sc2 = con.connection(self.pump2, 'out1', self.sc2, 'in1')
        self.sc2_v2 = con.connection(self.sc2, 'out1', v2, 'in1', p=3.1, m=0.1)
        v2_me = con.connection(v2, 'out1', me, 'in2')

        self.nw.add_conns(inlet, outlet, sp_p1, p1_sc1, self.sc1_v1, v1_me,
                          self.sp_p2, self.p2_sc2, self.sc2_v2, v2_me)

    def test_turbomachine(self):
        """
        Test component properties of turbomachines.
        """
        instance = cmp.turbomachine('turbomachine')
        c1, c2 = self.setup_network_11(instance)
        fl = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129, 'INCOMP::DowQ': 0,
              'H2O': 0, 'NH3': 0, 'CO2': 0, 'CH4': 0}
        c1.set_attr(fluid=fl, m=10, p=1, h=1e5)
        c2.set_attr(p=1, h=2e5)
        self.nw.solve('design')
        power = c1.m.val_SI * (c2.h.val_SI - c1.h.val_SI)
        pr = c2.p.val_SI / c1.p.val_SI
        # pressure ratio and power are the basic functions for turbomachines,
        # these are inherited by all children, thus only tested here
        msg = ()
        eq_(power, instance.P.val, 'Value of power must be ' + str(power) +
            ', is ' + str(instance.P.val) + '.')
        msg = ()
        eq_(pr, instance.pr.val, 'Value of power must be ' + str(pr) +
            ', is ' + str(instance.pr.val) + '.')
        c2.set_attr(p=np.nan)
        instance.set_attr(pr=5)
        self.nw.solve('design')
        pr = c2.p.val_SI / c1.p.val_SI
        msg = ('Value of power must be ' + str(pr) + ', is ' +
               str(instance.pr.val) + '.')
        eq_(pr, instance.pr.val, msg)
        c2.set_attr(h=np.nan)
        instance.set_attr(P=1e5)
        self.nw.solve('design')
        power = c1.m.val_SI * (c2.h.val_SI - c1.h.val_SI)
        msg = ('Value of power must be ' + str(pr) + ', is ' +
               str(instance.pr.val) + '.')
        eq_(pr, instance.pr.val, msg)
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
        """
        Test component properties of pumps.
        """
        instance = cmp.pump('pump')
        c1, c2 = self.setup_network_11(instance)
        fl = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 1, 'H2O': 0,
              'NH3': 0, 'CO2': 0, 'CH4': 0}
        c1.set_attr(fluid=fl, v=1, p=5, T=50)
        c2.set_attr(p=7)
        instance.set_attr(eta_s=1)
        self.nw.solve('design')
        # calculate isentropic efficiency the old fashioned way
        eta_s = (instance.h_os('') - c1.h.val_SI) / (c2.h.val_SI - c1.h.val_SI)
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
               ', is ' + str(instance.eta_s.val) + '.')
        eq_(eta_s, instance.eta_s.val, msg)
        s1 = round(hlp.s_mix_ph(c1.to_flow()), 4)
        s2 = round(hlp.s_mix_ph(c2.to_flow()), 4)
        msg = ('Value of entropy must be identical for inlet (' + str(s1) +
               ') and outlet (' + str(s2) +
               ') at 100 % isentropic efficiency.')
        eq_(s1, s2, msg)
        instance.set_attr(eta_s=0.7)
        self.nw.solve('design')
        self.nw.save('tmp')
        c2.set_attr(p=np.nan)
        # flow char (pressure rise vs. volumetric flow)
        x = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4]
        y = np.array([14, 13.5, 12.5, 11, 9, 6.5, 3.5, 0]) * 1e5
        char = hlp.dc_cc(x=x, y=y, is_set=True)
        # apply flow char and eta_s char
        instance.set_attr(flow_char=char, eta_s=np.nan,
                          eta_s_char=hlp.dc_cc(method='GENERIC', is_set=True))
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of pressure rise must be ' + str(650000) + ', is ' +
               str(c2.p.val_SI - c1.p.val_SI) + '.')
        eq_(round(c2.p.val_SI - c1.p.val_SI, 0), 650000, msg)
        c1.set_attr(v=0.9)
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of pressure rise must be ' + str(775000.0) + ', is ' +
               str(c2.p.val_SI - c1.p.val_SI) + '.')
        eq_(c2.p.val_SI - c1.p.val_SI, 775000.0, msg)
        msg = ('Value of isentropic efficiency must be ' + str(0.694) +
               ', is ' + str(instance.eta_s.val) + '.')
        eq_(0.694, round(instance.eta_s.val, 3), msg)
        instance.eta_s_char.is_set = False
        # test boundaries of characteristic line
        c2.set_attr(T=con.ref(c1, 0, 20))
        c1.set_attr(v=-0.1)
        self.nw.solve('design')
        msg = ('Value of power must be ' + str(14e5) + ', is ' +
               str(c2.p.val_SI - c1.p.val_SI) + '.')
        eq_(c2.p.val_SI - c1.p.val_SI, 14e5, msg)
        c1.set_attr(v=1.5)
        self.nw.solve('design')
        msg = ('Value of power must be ' + str(0) + ', is ' +
               str(c2.p.val_SI - c1.p.val_SI) + '.')
        eq_(c2.p.val_SI - c1.p.val_SI, 0, msg)
        shutil.rmtree('./tmp', ignore_errors=True)

    def test_compressor(self):
        """
        Test component properties of compressors.
        """
        instance = cmp.compressor('compressor')
        c1, c2 = self.setup_network_11(instance)
        fl = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 0,
              'H2O': 0, 'NH3': 1, 'CO2': 0, 'CH4': 0}
        c1.set_attr(fluid=fl, v=1, p=5, T=100)
        c2.set_attr(p=7)
        instance.set_attr(eta_s=0.8)
        self.nw.solve('design')
        self.nw.save('tmp')
        # calculate isentropic efficiency the old fashioned way
        eta_s_d = ((instance.h_os('') - c1.h.val_SI) /
                   (c2.h.val_SI - c1.h.val_SI))
        msg = ('Value of isentropic efficiency must be ' + str(eta_s_d) +
               ', is ' + str(instance.eta_s.val) + '.')
        eq_(round(eta_s_d, 3), round(instance.eta_s.val, 3), msg)
        # trigger invalid isentropic efficiency
        instance.set_attr(eta_s=1.1)
        self.nw.solve('design')
        # calculate isentropic efficiency the old fashioned way
        eta_s = (instance.h_os('') - c1.h.val_SI) / (c2.h.val_SI - c1.h.val_SI)
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
               ', is ' + str(instance.eta_s.val) + '.')
        eq_(round(eta_s, 3), round(instance.eta_s.val, 3), msg)
        c2.set_attr(p=np.nan)
        instance.set_attr(char_map=hlp.dc_cm(method='GENERIC', is_set=True),
                          eta_s=np.nan)
        self.nw.solve('offdesign', design_path='tmp')
        eta_s = (instance.h_os('') - c1.h.val_SI) / (c2.h.val_SI - c1.h.val_SI)
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
               ') must be identical to design case (' + str(eta_s) + ').')
        eq_(round(eta_s, 2), round(instance.eta_s.val, 2), msg)
        # going above highes available speedline, beneath lowest mass flow at
        # that line
        c1.set_attr(v=np.nan, m=c1.m.val*0.8, T=30)
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
               ') must be at (' +
               str(round(eta_s * instance.char_map.z2[6, 0], 4)) + ').')
        eq_(round(eta_s * instance.char_map.z2[6, 0], 4),
            round(instance.eta_s.val, 4), msg)
        # going below lowest available speedline, above highest mass flow at
        # that line
        c1.set_attr(T=300)
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
               ') must be at (' +
               str(round(eta_s * instance.char_map.z2[0, 9], 4)) + ').')
        eq_(round(eta_s * instance.char_map.z2[0, 9], 4),
            round(instance.eta_s.val, 4), msg)
        # back to design properties, test eta_s_char
        c2.set_attr(p=7)
        c1.set_attr(v=1, T=100, m=np.nan)
        # test param specification m
        instance.set_attr(eta_s_char=hlp.dc_cc(method='GENERIC',
                                               is_set=True, param='m'))
        instance.char_map.is_set = False
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
               ', is ' + str(instance.eta_s.val) + '.')
        eq_(round(eta_s, 3), round(instance.eta_s.val, 3), msg)
        c1.set_attr(v=1.5)
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of isentropic efficiency must be ' + str(0.88) +
               ', is ' + str(instance.eta_s.val) + '.')
        eq_(0.88, round(instance.eta_s.val, 3), msg)
        # test param specification pr
        instance.eta_s_char.set_attr(param='pr')
        c1.set_attr(v=1)
        c2.set_attr(p=7.5)
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of isentropic efficiency must be ' + str(0.829) +
               ', is ' + str(instance.eta_s.val) + '.')
        eq_(0.829, round(instance.eta_s.val, 3), msg)
        instance.eta_s_char.set_attr(param=None)
        # test for missing parameter declaration
        try:
            self.nw.solve('offdesign', design_path='tmp')
        except ValueError:
            pass
        shutil.rmtree('./tmp', ignore_errors=True)

    def test_turbine(self):
        """
        Test component properties of turbines.
        """
        instance = cmp.turbine('turbine')
        c1, c2 = self.setup_network_11(instance)
        fl = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129, 'INCOMP::DowQ': 0,
              'H2O': 0, 'NH3': 0, 'CO2': 0, 'CH4': 0}
        c1.set_attr(fluid=fl, m=15, p=10)
        c2.set_attr(p=1, T=20)
        instance.set_attr(eta_s=0.8)
        self.nw.solve('design')
        self.nw.save('tmp')
        # calculate isentropic efficiency the old fashioned way
        eta_s_d = ((c2.h.val_SI - c1.h.val_SI) /
                   (instance.h_os('') - c1.h.val_SI))
        msg = ('Value of isentropic efficiency must be ' + str(eta_s_d) +
               ', is ' + str(instance.eta_s.val) + '.')
        eq_(round(eta_s_d, 3), round(instance.eta_s.val, 3), msg)
        # trigger invalid isentropic efficiency
        instance.set_attr(eta_s=1.1)
        self.nw.solve('design')
        # calculate isentropic efficiency the old fashioned way
        eta_s = (c2.h.val_SI - c1.h.val_SI) / (instance.h_os('') - c1.h.val_SI)
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
               ', is ' + str(instance.eta_s.val) + '.')
        eq_(round(eta_s, 3), round(instance.eta_s.val, 3), msg)
        c1.set_attr(p=np.nan)
        instance.cone.is_set = True
        instance.eta_s_char.is_set = True
        instance.eta_s.is_set = False
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
               ') must be identical to design case (' + str(eta_s_d) + ').')
        eq_(round(eta_s_d, 2), round(instance.eta_s.val, 2), msg)
        # lowering mass flow, inlet pressure must sink according to cone law
        c1.set_attr(m=c1.m.val * 0.8)
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of pressure ratio (' + str(instance.pr.val) +
               ') must be at (' + str(0.128) + ').')
        eq_(0.128, round(instance.pr.val, 3), msg)
        self.nw.print_results()
        # testing more parameters for eta_s_char
        c1.set_attr(m=10)
        # test param specification v
        instance.eta_s_char.param = 'v'
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
               ') must be (' + str(0.8) + ').')
        eq_(0.8, round(instance.eta_s.val, 3), msg)
        # test param specification pr
        instance.eta_s_char.param = 'pr'
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
               ') must be (' + str(0.768) + ').')
        eq_(0.768, round(instance.eta_s.val, 3), msg)
        # test param specification dh_s
        instance.eta_s_char.param = 'dh_s'
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
               ') must be (' + str(0.798) + ').')
        eq_(0.798, round(instance.eta_s.val, 3), msg)
        instance.eta_s_char.param = None
        # test for missing parameter declaration
        try:
            self.nw.solve('offdesign', design_path='tmp')
        except ValueError:
            pass
        shutil.rmtree('./tmp', ignore_errors=True)

    def test_combustion_chamber(self):
        """
        Test component properties of combustion chambers.
        """
        instance = cmp.combustion_chamber('combustion chamber', fuel='CH4')
        c1, c2, c3 = self.setup_network_21(instance)
        air = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129, 'INCOMP::DowQ': 0,
               'H2O': 0, 'NH3': 0, 'CO2': 0, 'CH4': 0}
        fuel = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 0, 'H2O': 0,
                'NH3': 0, 'CO2': 0.04, 'CH4': 0.96}
        c1.set_attr(fluid=air, p=1, T=30)
        c2.set_attr(fluid=fuel, T=30)
        c3.set_attr(T=1200)
        b = con.bus('thermal input', P=1e6)
        b.add_comps({'c': instance})
        self.nw.add_busses(b)
        self.nw.solve('design')
        self.nw.print_results()
        msg = ('Value of thermal input must be ' + str(b.P.val) + ', is ' +
               str(instance.ti.val) + '.')
        eq_(round(b.P.val, 1), round(instance.ti.val, 1), msg)

        # test unspecified fuel
        instance.set_attr(fuel=np.nan)
        try:
            self.nw.solve('design')
        except hlp.TESPyComponentError:
            pass

        # test wrongly specified fuel
        instance.set_attr(fuel='Ar')
        try:
            self.nw.solve('design')
        except hlp.TESPyComponentError:
            pass

    def test_valve(self):
        """
        Test component properties of valves.
        """
        instance = cmp.valve('valve')
        c1, c2 = self.setup_network_11(instance)
        fl = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129, 'INCOMP::DowQ': 0,
              'H2O': 0, 'NH3': 0, 'CO2': 0, 'CH4': 0}
        c1.set_attr(fluid=fl, m=10, p=10, T=120)
        c2.set_attr(p=1)
        instance.set_attr(pr='var')
        self.nw.solve('design')
        msg = ('Value of pressure ratio must be ' +
               str(c2.p.val_SI / c1.p.val_SI) + ', is ' +
               str(instance.pr.val) + '.')
        eq_(round(c2.p.val_SI / c1.p.val_SI, 1),
            round(instance.pr.val, 1), msg)
        zeta = instance.zeta.val
        instance.set_attr(zeta='var', pr=np.nan)
        self.nw.solve('design')
        msg = ('Value of pressure ratio must be ' + str(zeta) + ', is ' +
               str(instance.zeta.val) + '.')
        eq_(round(zeta, 1), round(instance.zeta.val, 1), msg)

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
        amb_comb = con.connection(amb, 'out1', chp, 'in3')
        sf_comb = con.connection(sf, 'out1', chp, 'in4')
        comb_fg = con.connection(chp, 'out3', fg, 'in1')
        self.nw.add_conns(sf_comb, amb_comb, comb_fg)
        cw1_chp1 = con.connection(cw_in1, 'out1', chp, 'in1')
        cw2_chp2 = con.connection(cw_in2, 'out1', chp, 'in2')
        self.nw.add_conns(cw1_chp1, cw2_chp2)
        chp1_cw = con.connection(chp, 'out1', cw_out1, 'in1')
        chp2_cw = con.connection(chp, 'out2', cw_out2, 'in1')
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

        TI = con.bus('thermal input')
        Q1 = con.bus('heat output 1')
        Q2 = con.bus('heat output 2')
        Q = con.bus('heat output')
        Qloss = con.bus('thermal heat loss')

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

    def test_heat_ex_simple(self):
        """
        Test component properties of simple heat exchanger.
        """
        instance = cmp.heat_exchanger_simple('heat exchanger')
        c1, c2 = self.setup_network_11(instance)
        fl = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 0, 'H2O': 1,
              'NH3': 0, 'CO2': 0, 'CH4': 0}
        c1.set_attr(fluid=fl, m=1, p=10, T=100)
        # trigger heat exchanger parameter groups
        instance.set_attr(hydro_group='HW', L=100, ks=100, pr=0.99, Tamb=20)
        instance.hydro_group.is_set = True
        instance.kA_group.is_set = True
        self.nw.solve('design', init_only=True)
        msg = ('Hydro group must no be set, if one parameter is missing!')
        eq_(instance.hydro_group.is_set, False, msg)
        msg = ('kA group must no be set, if one parameter is missing!')
        eq_(instance.kA_group.is_set, False, msg)
        instance.set_attr(hydro_group='HW', D='var', L=100,
                          ks=100, pr=0.99, Tamb=20)
        b = con.bus('heat', P=-1e5)
        b.add_comps({'c': instance})
        self.nw.add_busses(b)
        self.nw.solve('design')
        msg = ('Value of pressure ratio must be ' +
               str(c2.p.val_SI / c1.p.val_SI) + ', is ' +
               str(instance.pr.val) + '.')
        eq_(round(c2.p.val_SI / c1.p.val_SI, 2),
            round(instance.pr.val, 2), msg)
        zeta = instance.zeta.val
        instance.set_attr(D=instance.D.val, zeta='var', pr=np.nan)
        instance.D.is_var = False
        self.nw.solve('design')
        msg = ('Value of zeta must be ' + str(zeta) + ', is ' +
               str(instance.zeta.val) + '.')
        eq_(round(zeta, 1), round(instance.zeta.val, 1), msg)
        pr = instance.pr.val
        instance.set_attr(zeta=np.nan, pr='var')
        self.nw.solve('design')
        msg = ('Value of pressure ratio must be ' + str(pr) +
               ', is ' + str(instance.pr.val) + '.')
        eq_(round(pr, 2), round(instance.pr.val, 2), msg)
        instance.set_attr(kA='var', pr=np.nan)
        b.set_attr(P=-5e4)
        self.nw.solve('design')
        # due to heat output being half of reference (for Tamb) kA should be
        # somewhere near to that (actual value is 677)
        msg = ('Value of heat transfer coefficient must be ' + str(677) +
               ', is ' + str(instance.kA.val) + '.')
        eq_(677, round(instance.kA.val, 0), msg)

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
        instance = cmp.solar_collector('solar collector')
        c1, c2 = self.setup_network_11(instance)
        fl = {'N2': 0, 'O2': 0, 'Ar': 0, 'INCOMP::DowQ': 0,
              'H2O': 1, 'NH3': 0, 'CO2': 0, 'CH4': 0}
        c1.set_attr(fluid=fl, m=1, p=10, T=100)
        # trigger solar collector parameter groups (see heat exchanger simple)
        instance.set_attr(hydro_group='HW', L=100, ks=100, pr=0.99, Tamb=20)
        instance.hydro_group.is_set = True
        instance.energy_group.is_set = True
        self.nw.solve('design', init_only=True)
        msg = ('Hydro group must no be set, if one parameter is missing!')
        eq_(instance.hydro_group.is_set, False, msg)
        msg = ('Energy group must no be set, if one parameter is missing!')
        eq_(instance.energy_group.is_set, False, msg)

    def test_heat_ex(self):
        """
        Test component properties of heat exchanger.
        """
        tesin = cmp.sink('TES in')
        tesout = cmp.source('TES out')
        hsin = cmp.sink('HS in')
        hsout = cmp.source('HS out')
        he = cmp.heat_exchanger('heat exchanger')
        tes_he = con.connection(tesout, 'out1', he, 'in2')
        he_tes = con.connection(he, 'out2', tesin, 'in1')
        hs_he = con.connection(hsout, 'out1', he, 'in1')
        he_hs = con.connection(he, 'out1', hsin, 'in1')
        self.nw.add_conns(tes_he, he_tes, hs_he, he_hs)
        # design specification
        he.set_attr(pr1=0.98, pr2=0.98, ttd_u=5,
                    design=['pr1', 'pr2', 'ttd_u'],
                    offdesign=['zeta1', 'zeta2', 'kA'])
        hs_he.set_attr(T=120, p=3, fluid={'N2': 0, 'O2': 0, 'Ar': 0,
                                          'INCOMP::DowQ': 0, 'H2O': 1,
                                          'NH3': 0, 'CO2': 0, 'CH4': 0})
        he_hs.set_attr(T=70)
        tes_he.set_attr(T=40, p=5, fluid={'N2': 0, 'O2': 0, 'Ar': 1,
                                          'INCOMP::DowQ': 0, 'H2O': 0,
                                          'NH3': 0, 'CO2': 0, 'CH4': 0})
        b = con.bus('heat transfer', P=-80e3)
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
        tesin = cmp.sink('TES in')
        tesout = cmp.source('TES out')
        hsin = cmp.sink('Cond in')
        hsout = cmp.source('Cond out')
        he = cmp.condenser('condenser')
        tes_he = con.connection(tesout, 'out1', he, 'in2')
        he_tes = con.connection(he, 'out2', tesin, 'in1')
        hs_he = con.connection(hsout, 'out1', he, 'in1')
        he_hs = con.connection(he, 'out1', hsin, 'in1')
        self.nw.add_conns(tes_he, he_tes, hs_he, he_hs)
        # design specification
        he.set_attr(pr1=0.98, pr2=0.98, ttd_u=5, design=['pr2', 'ttd_u'],
                    offdesign=['zeta2', 'kA'])
        hs_he.set_attr(T=100, p0=0.5, fluid={'N2': 0, 'O2': 0, 'Ar': 0,
                                             'INCOMP::DowQ': 0, 'H2O': 1,
                                             'NH3': 0, 'CO2': 0, 'CH4': 0})
        tes_he.set_attr(T=30, p=5, fluid={'N2': 0, 'O2': 0, 'Ar': 0,
                                          'INCOMP::DowQ': 0, 'H2O': 1,
                                          'NH3': 0, 'CO2': 0, 'CH4': 0})
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

    def test_water_electrolyzer(self):
        """
        Test component properties of water electrolyzer.
        """
        self.nw = nwk.network(['H2O', 'O2'])
        instance = cmp.water_electrolyzer('electrolyzer')
        self.setup_network_electrolyzer(instance)

        try:
            self.nw.solve('design')
        except ValueError:
            pass

        self.nw = nwk.network(['H2O', 'H2'])
        instance = cmp.water_electrolyzer('electrolyzer')
        self.setup_network_electrolyzer(instance)

        try:
            self.nw.solve('design')
        except ValueError:
            pass

        self.nw = nwk.network(['O2', 'H2'])
        instance = cmp.water_electrolyzer('electrolyzer')
        self.setup_network_electrolyzer(instance)

        try:
            self.nw.solve('design')
        except ValueError:
            pass

        self.nw = nwk.network(['O2', 'H2', 'H2O'], T_unit='C', p_unit='bar')
        instance = cmp.water_electrolyzer('electrolyzer')
        self.setup_network_electrolyzer(instance)

        # check bus functions
        power = con.bus('power')
        power.add_comps({'c': instance, 'p': 'P'})
        power.set_attr(P=2.5e6)
        self.nw.add_busses(power)

        self.nw.solve('design')
        msg = ('Value of power must be ' + str(power.P.val) + ', is ' +
               str(instance.P.val) + '.')
        eq_(round(power.P.val, 1), round(instance.P.val), msg)

        # check bus functions
        power.set_attr(P=np.nan)
        heat = con.bus('heat')
        heat.add_comps({'c': instance, 'p': 'Q'})
        heat.set_attr(P=-8e5)
        self.nw.add_busses(heat)

        self.nw.solve('design')
        msg = ('Value of heat flow must be ' + str(heat.P.val) +
               ', is ' + str(instance.Q.val) + '.')
        eq_(round(heat.P.val, 1), round(instance.Q.val), msg)
        self.nw.save('tmp')

        Q = heat.P.val * 0.9
        heat.set_attr(P=Q)
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of heat flow must be ' + str(Q) +
               ', is ' + str(instance.Q.val) + '.')
        eq_(round(Q, 1), round(instance.Q.val), msg)

        self.nw.del_busses(heat, power)

        # check invalid bus parameter error
        some_bus = con.bus('some_bus')
        some_bus.add_comps({'c': instance, 'p': 'G'})
        some_bus.set_attr(P=-8e5)
        self.nw.add_busses(some_bus)
        try:
            self.nw.solve('design')
        except ValueError:
            pass

        # check invald bus parameter error in derivatives
        try:
            instance.bus_deriv(some_bus.comps.loc[instance])
        except ValueError:
            pass

        self.nw.del_busses(some_bus)

        instance.set_attr(eta=0.9, e='var')
        self.nw.solve('design')
        msg = ('Value of efficiency must be ' + str(instance.eta.val) +
               ', is ' + str(instance.e0 / instance.e.val) + '.')
        eq_(round(instance.eta.val, 2), round(instance.e0 / instance.e.val, 2),
            msg)

        # isentropic efficiency value > 1
        e = 130e6
        instance.set_attr(e=np.nan, eta=np.nan)
        instance.set_attr(e=e)
        self.nw.solve('design')
        msg = ('Value of efficiency must be ' + str(instance.e0 / e) +
               ', is ' + str(instance.eta.val) + '.')
        eq_(round(instance.e0 / e, 2), round(instance.eta.val, 2), msg)

        e = 150e6
        instance.set_attr(e=np.nan, eta=np.nan)
        instance.set_attr(e=e)
        self.nw.solve('design')
        msg = ('Value of efficiency must be ' + str(e) + ', is ' +
               str(instance.e.val) + '.')
        eq_(round(e, 1), round(instance.e.val, 1), msg)

        pr = 0.95
        instance.set_attr(pr_c=pr, e=np.nan, zeta='var',
                          P=2.5e6, design=['pr_c'])
        self.nw.solve('design')
        self.nw.save('tmp')
        msg = ('Value of pressure ratio must be ' + str(pr) + ', is ' +
               str(instance.pr_c.val) + '.')
        eq_(round(pr, 2), round(instance.pr_c.val, 2), msg)

        instance.set_attr(zeta=np.nan, offdesign=['zeta'])
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of pressure ratio must be ' + str(pr) + ', is ' +
               str(instance.pr_c.val) + '.')
        eq_(round(pr, 2), round(instance.pr_c.val, 2), msg)

        Q = instance.Q.val * 0.9
        instance.set_attr(Q=Q, P=np.nan)
        self.nw.solve('offdesign', design_path='tmp')
        msg = ('Value of heat must be ' + str(Q) + ', is ' +
               str(instance.Q.val) + '.')
        eq_(round(Q, 0), round(instance.Q.val, 0), msg)
        shutil.rmtree('./tmp', ignore_errors=True)

    def test_individual_design_path_on_connections_and_components(self):

        self.setup_network_individual_offdesign()
        self.nw.solve('design')
        self.sc2_v2.set_attr(m=0)
        self.nw.solve('design')
        self.nw.save('design1')
        v1_design = self.sc1_v1.v.val_SI
        zeta_sc1_design = self.sc1.zeta.val

        self.sc2_v2.set_attr(T=95, m=np.nan)
        self.sc1_v1.set_attr(T=np.nan, m=0.001)
        self.nw.solve('design')
        self.nw.save('design2')
        v2_design = self.sc2_v2.v.val_SI
        zeta_sc2_design = self.sc2.zeta.val

        self.sc1_v1.set_attr(m=np.nan)
        self.sc1_v1.set_attr(design=['T'], offdesign=['v'], state='l')
        self.sc2_v2.set_attr(design=['T'], offdesign=['v'], state='l')

        self.sc2.set_attr(design_path='design2')
        self.pump2.set_attr(design_path='design2')
        self.sp_p2.set_attr(design_path='design2')
        self.p2_sc2.set_attr(design_path='design2')
        self.sc2_v2.set_attr(design_path='design2')
        self.nw.solve('offdesign', design_path='design1')

        self.sc1.set_attr(E=500)
        self.sc2.set_attr(E=950)

        self.nw.solve('offdesign', design_path='design1')
        self.sc2_v2.set_attr(design_path=np.nan)

        # volumetric flow comparison
        msg = ('Design path was set to None, is ' +
               str(self.sc2_v2.design_path) + '.')
        eq_(None, self.sc2_v2.design_path, msg)

        # volumetric flow comparison
        msg = ('Value of volumetric flow must be ' + str(v1_design) + ', is ' +
               str(self.sc1_v1.v.val_SI) + '.')
        eq_(round(v1_design, 5), round(self.sc1_v1.v.val_SI, 5), msg)

        msg = ('Value of volumetric flow must be ' + str(v2_design) + ', is ' +
               str(self.sc2_v2.v.val_SI) + '.')
        eq_(round(v2_design, 5), round(self.sc2_v2.v.val_SI, 5), msg)

        # zeta value of solar collector comparison
        msg = ('Value of zeta must be ' + str(zeta_sc1_design) + ', is ' +
               str(self.sc1.zeta.val) + '.')
        eq_(round(zeta_sc1_design, 0), round(self.sc1.zeta.val, 0), msg)

        msg = ('Value of zeta must be ' + str(zeta_sc2_design) + ', is ' +
               str(self.sc2.zeta.val) + '.')
        eq_(round(zeta_sc2_design, 0), round(self.sc2.zeta.val, 0), msg)

        shutil.rmtree('./design1', ignore_errors=True)
        shutil.rmtree('./design2', ignore_errors=True)
