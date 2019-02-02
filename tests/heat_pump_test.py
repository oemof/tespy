# -*- coding: utf-8

from nose.tools import eq_

from tespy import cmp, con, nwk, hlp, cmp_char

import numpy as np


class test_heat_pump_ebsilon:

    def setup(self):
        # %% network

        self.nw = nwk.network(fluids=['water', 'NH3'], T_unit='C', p_unit='bar', h_unit='kJ / kg', m_unit='kg / s')

        # %% components

        # sources & sinks

        c_in = cmp.source('coolant in')
        cb = cmp.source('consumer back flow')
        cf = cmp.sink('consumer feed flow')
        amb_in = cmp.source('source ambient')
        amb_out = cmp.sink('sink ambient')
        ic_in = cmp.source('source intercool')
        ic_out = cmp.sink('sink intercool')

        c_out = cmp.sink('coolant out')

        # consumer system

        cd = cmp.heat_exchanger('condenser')
        rp = cmp.pump('recirculation pump')
        cons = cmp.heat_exchanger_simple('consumer')

        # evaporator system

        va = cmp.valve('valve')
        dr = cmp.drum('drum')
        ev = cmp.heat_exchanger('evaporator')
        su = cmp.heat_exchanger('superheater')
        pu = cmp.pump('pump evaporator')

        # compressor-system

        cp1 = cmp.compressor('compressor 1')
        cp2 = cmp.compressor('compressor 2')
        he = cmp.heat_exchanger('intercooler')

        # busses

        x = np.array([0, 0.7, 1, 1.3])
        y = 1 / np.array([0.5, 0.95, 1, 0.98]) / 0.9583794
        motor = cmp_char.characteristics(x=x, y=y)

        self.power = con.bus('total compressor power')
        self.power.add_comps({'c': cp1, 'char': motor}, {'c': cp2, 'char': motor})
        self.heat = con.bus('total delivered heat')
        self.heat.add_comps({'c': cd, 'char': -1})
        self.nw.add_busses(self.power, self.heat)

        # %% connections

        # consumer system

        c_in_cd = con.connection(c_in, 'out1', cd, 'in1')

        cb_rp = con.connection(cb, 'out1', rp, 'in1')
        rp_cd = con.connection(rp, 'out1', cd, 'in2')
        self.cd_cons = con.connection(cd, 'out2', cons, 'in1')
        cons_cf = con.connection(cons, 'out1', cf, 'in1')

        self.nw.add_conns(c_in_cd, cb_rp, rp_cd, self.cd_cons, cons_cf)

        # connection condenser - evaporator system

        cd_va = con.connection(cd, 'out1', va, 'in1')

        self.nw.add_conns(cd_va)

        # evaporator system

        va_dr = con.connection(va, 'out1', dr, 'in1')
        dr_pu = con.connection(dr, 'out1', pu, 'in1')
        pu_ev = con.connection(pu, 'out1', ev, 'in2')
        ev_dr = con.connection(ev, 'out2', dr, 'in2')
        dr_su = con.connection(dr, 'out2', su, 'in2')

        self.nw.add_conns(va_dr, dr_pu, pu_ev, ev_dr, dr_su)

        self.amb_in_su = con.connection(amb_in, 'out1', su, 'in1')
        su_ev = con.connection(su, 'out1', ev, 'in1')
        ev_amb_out = con.connection(ev, 'out1', amb_out, 'in1')

        self.nw.add_conns(self.amb_in_su, su_ev, ev_amb_out)

        # connection evaporator system - compressor system

        su_cp1 = con.connection(su, 'out2', cp1, 'in1')

        self.nw.add_conns(su_cp1)

        # compressor-system

        cp1_he = con.connection(cp1, 'out1', he, 'in1')
        he_cp2 = con.connection(he, 'out1', cp2, 'in1')
        cp2_c_out = con.connection(cp2, 'out1', c_out, 'in1')

        ic_in_he = con.connection(ic_in, 'out1', he, 'in2')
        he_ic_out = con.connection(he, 'out2', ic_out, 'in1')

        self.nw.add_conns(cp1_he, he_cp2, ic_in_he, he_ic_out, cp2_c_out)

        # helping connection to reference boiling temperature

        ref_Ts_su = con.connection(cmp.source('so'), 'out1', cmp.sink('si'), 'in1')
        ref_Ts_cd = con.connection(cmp.source('socd'), 'out1', cmp.sink('sicd'), 'in1')
        self.nw.add_conns(ref_Ts_su, ref_Ts_cd)
        # %% component parametrization

        # condenser system

        rp.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'])
        cons.set_attr(pr=1, design=['pr'], offdesign=['zeta'])

        # evaporator system

        ev.set_attr(pr1=1, pr2=.999, ttd_l=5, design=['ttd_l'], offdesign=['kA'],
                    kA_char1='EVA_HOT', kA_char2='EVA_COLD')

        # characteristic line for superheater kA
        x = np.array([0, 0.045, 0.136, 0.244, 0.43, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2])
        y = np.array([0, 0.037, 0.112, 0.207, 0.5, 0.8, 0.85, 0.9, 0.95, 1, 1.04, 1.07])
        su_char = hlp.dc_cc(x=x, y=y, param='m')
        su.set_attr(kA_char1='default', kA_char2=su_char, offdesign=['zeta1', 'zeta2', 'kA'])
        pu.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'])

        # compressor system

        cp1.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'])
        cp2.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'])

        # characteristic line for intercooler kA
        x = np.linspace(0, 2.5, 26)
        y = np.array([
                0.000, 0.164, 0.283, 0.389, 0.488, 0.581, 0.670, 0.756, 0.840, 0.921,
                1.000, 1.078, 1.154, 1.228, 1.302, 1.374, 1.446, 1.516, 1.585, 1.654,
                1.722, 1.789, 1.855, 1.921, 1.986, 2.051])
        he_char_cold = hlp.dc_cc(x=x, y=y, param='m')

        he.set_attr(kA_char1='default', kA_char2=he_char_cold, offdesign=['zeta1', 'zeta2', 'kA'])
        cd.set_attr(pr2=0.998, design=['pr2'], offdesign=['zeta2', 'kA'])

        # %% connection parametrization

        # condenser system

        c_in_cd.set_attr(fluid={'water': 0, 'NH3': 1}, p=60)
        cb_rp.set_attr(T=60, p=10, fluid={'water': 1, 'NH3': 0})
        self.cd_cons.set_attr(T=105)
        cons_cf.set_attr(h=con.ref(cb_rp, 1, 0), p=con.ref(cb_rp, 1, 0))
        cd_va.set_attr(p=con.ref(c_in_cd, 1, -1000), T=con.ref(ref_Ts_cd, 1, -5), h0=500, design=['T'])

        # evaporator system cold side

        pu_ev.set_attr(m=con.ref(va_dr, 10, 0), p0=5)
        dr_su.set_attr(p0=5, T=5)
        su_cp1.set_attr(p=con.ref(dr_su, 1, -5000), T=con.ref(ref_Ts_su , 1, 5), h0=1700, design=['T', 'p'])

        # evaporator system hot side

        self.amb_in_su.set_attr(m=20, T=12, p=1, fluid={'water': 1, 'NH3': 0})
        su_ev.set_attr(p=con.ref(self.amb_in_su, 1, -100), design=['p'])
        ev_amb_out.set_attr()

        ref_Ts_su.set_attr(m=5, fluid={'water': 0, 'NH3': 1}, p=con.ref(su_cp1, 1, 0), x=1)
        ref_Ts_cd.set_attr(m=5, fluid={'water': 0, 'NH3': 1}, p=con.ref(cd_va, 1, 0), x=0)

        # compressor-system

        cp1_he.set_attr(p=15)
        he_cp2.set_attr(T=40, p=con.ref(cp1_he, 1, -1000), design=['T', 'p'])
        ic_in_he.set_attr(p=1, T=20, m=5, fluid={'water': 1, 'NH3': 0})
        he_ic_out.set_attr(p=con.ref(ic_in_he, 1, -200), design=['p'])
        cp2_c_out.set_attr(p=con.ref(c_in_cd, 1, 0), h=con.ref(c_in_cd, 1, 0))


    def test_model(self):
        """
        Tests the operating points of the heat pump against an ebsilon model.
        By now, not all characteristic functions are available in detail,
        thus perfect matching is not possible!
        """

        self.nw.solve('design')
        self.nw.save('tmp')

        # input values from ebsilon
        T = [105, 100, 90, 80]
        m_source = np.array([[23, 22, 20, 18, 16],
                             [27, 24, 20, 16, 12],
                             [31, 25, 20, 15, 10],
                             [33, 25, 20, 15, 10]])
        COP = np.array([[2.323, 2.306, 2.269, 2.221, 2.153],
                        [2.464, 2.410, 2.329, 2.202, 2.031],
                        [2.641, 2.524, 2.409, 2.250, 1.94],
                        [2.726, 2.577, 2.462, 2.310, 1.988]])
        P = np.array([[0.1064, 0.1041, 0.0991, 0.0937, 0.0883],
                      [0.1098, 0.1039, 0.0948, 0.0842, 0.0725],
                      [0.1099, 0.1000, 0.0888, 0.0754, 0.0620],
                      [0.1099, 0.0968, 0.0852, 0.0717, 0.0585]]) * 1e6

        i = 0
        for T in T:
            j = 0
            self.cd_cons.set_attr(T=T)
            for m in m_source[i]:
                self.amb_in_su.set_attr(m=m)
                self.nw.solve('offdesign', design_path='tmp', init_path='tmp')
                # relative deviation should not exceed 10 %
                # this should be much less, unfortunately not all ebsilon characteristics are available,
                # thus it is difficult/impossible to match the models perfectly!
                d_rel_COP = abs(self.heat.P.val / self.power.P.val - COP[i, j]) / COP[i, j]
                d_rel_P = abs(self.power.P.val - P[i, j]) / P[i, j]
                print(P[i, j], COP[i, j], d_rel_COP, d_rel_P)
                eq_(d_rel_COP < 0.1, True, 'The deviation in COP should be less than 0.1, is ' + str(d_rel_COP) + ' at mass flow ' + str(m) + '.')
                eq_(d_rel_P < 0.1, True, 'The deviation in power should be less than 0.1, is ' + str(d_rel_P) + ' at mass flow ' + str(m) + '.')
                j += 1
            i += 1
