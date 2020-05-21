# -*- coding: utf-8

"""Module for testing a tespy simulation vs results from a different simulator.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_models/test_heat_pump_model.py

SPDX-License-Identifier: MIT
"""
import shutil

import numpy as np

from tespy.components.basics import sink
from tespy.components.basics import source
from tespy.components.heat_exchangers import heat_exchanger
from tespy.components.heat_exchangers import heat_exchanger_simple
from tespy.components.nodes import drum
from tespy.components.piping import valve
from tespy.components.turbomachinery import compressor
from tespy.components.turbomachinery import pump
from tespy.connections import bus
from tespy.connections import connection
from tespy.connections import ref
from tespy.networks.networks import network
from tespy.tools.characteristics import char_line
from tespy.tools.data_containers import dc_cc


class TestHeatPump:

    def setup(self):
        # %% network setup
        self.nw = network(fluids=['water', 'NH3'], T_unit='C', p_unit='bar',
                          h_unit='kJ / kg', m_unit='kg / s')

        # %% components
        # sources & sinks
        c_in = source('coolant in')
        cb = source('consumer back flow')
        cf = sink('consumer feed flow')
        amb_in = source('source ambient')
        amb_out = sink('sink ambient')
        ic_in = source('source intercool')
        ic_out = sink('sink intercool')

        c_out = sink('coolant out')

        # consumer system
        cd = heat_exchanger('condenser')
        rp = pump('recirculation pump')
        cons = heat_exchanger_simple('consumer')

        # evaporator system
        va = valve('valve')
        dr = drum('drum')
        ev = heat_exchanger('evaporator')
        su = heat_exchanger('superheater')
        pu = pump('pump evaporator')

        # compressor-system
        cp1 = compressor('compressor 1')
        cp2 = compressor('compressor 2')
        he = heat_exchanger('intercooler')

        # busses
        self.power = bus('total compressor power')
        self.power.add_comps({'comp': cp1}, {'comp': cp2})
        self.heat = bus('total delivered heat')
        self.heat.add_comps({'comp': cd, 'char': -1})
        self.nw.add_busses(self.power, self.heat)

        # %% connections
        # consumer system
        c_in_cd = connection(c_in, 'out1', cd, 'in1')

        cb_rp = connection(cb, 'out1', rp, 'in1')
        rp_cd = connection(rp, 'out1', cd, 'in2')
        self.cd_cons = connection(cd, 'out2', cons, 'in1')
        cons_cf = connection(cons, 'out1', cf, 'in1')

        self.nw.add_conns(c_in_cd, cb_rp, rp_cd, self.cd_cons, cons_cf)

        # connection condenser - evaporator system
        cd_va = connection(cd, 'out1', va, 'in1')

        self.nw.add_conns(cd_va)

        # evaporator system
        va_dr = connection(va, 'out1', dr, 'in1')
        dr_pu = connection(dr, 'out1', pu, 'in1')
        pu_ev = connection(pu, 'out1', ev, 'in2')
        ev_dr = connection(ev, 'out2', dr, 'in2')
        dr_su = connection(dr, 'out2', su, 'in2')

        self.nw.add_conns(va_dr, dr_pu, pu_ev, ev_dr, dr_su)

        self.amb_in_su = connection(amb_in, 'out1', su, 'in1')
        su_ev = connection(su, 'out1', ev, 'in1')
        ev_amb_out = connection(ev, 'out1', amb_out, 'in1')

        self.nw.add_conns(self.amb_in_su, su_ev, ev_amb_out)

        # connection evaporator system - compressor system
        su_cp1 = connection(su, 'out2', cp1, 'in1')

        self.nw.add_conns(su_cp1)

        # compressor-system
        cp1_he = connection(cp1, 'out1', he, 'in1')
        he_cp2 = connection(he, 'out1', cp2, 'in1')
        cp2_c_out = connection(cp2, 'out1', c_out, 'in1')

        ic_in_he = connection(ic_in, 'out1', he, 'in2')
        he_ic_out = connection(he, 'out2', ic_out, 'in1')

        self.nw.add_conns(cp1_he, he_cp2, ic_in_he, he_ic_out, cp2_c_out)

        # %% component parametrization
        # condenser system
        x = np.array(
            [0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5,
             0.5625, 0.6375, 0.7125, 0.7875, 0.9, 0.9875, 1, 1.0625, 1.125,
             1.175, 1.2125, 1.2375, 1.25])
        y = np.array(
            [0.0076, 0.1390, 0.2731, 0.4003, 0.5185, 0.6263, 0.7224, 0.8056,
             0.8754, 0.9312, 0.9729, 1.0006, 1.0203, 1.0158, 1.0051, 1.0000,
             0.9746, 0.9289, 0.8832, 0.8376, 0.7843, 0.7614])
        rp.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'],
                    eta_s_char=dc_cc(func=char_line(x, y), param='m'))
        cons.set_attr(pr=1, design=['pr'], offdesign=['zeta'])

        # evaporator system
        x = np.linspace(0, 2.5, 26)
        y = np.array(
            [0.000, 0.164, 0.283, 0.389, 0.488, 0.581, 0.670, 0.756, 0.840,
             0.921, 1.000, 1.078, 1.154, 1.228, 1.302, 1.374, 1.446, 1.516,
             1.585, 1.654, 1.722, 1.789, 1.855, 1.921, 1.986, 2.051])
        kA_char1 = dc_cc(func=char_line(x, y), param='m')

        x = np.array(
            [0.0100, 0.0400, 0.0700, 0.1100, 0.1500, 0.2000, 0.2500, 0.3000,
             0.3500, 0.4000, 0.4500, 0.5000, 0.5500, 0.6000, 0.6500, 0.7000,
             0.7500, 0.8000, 0.8500, 0.9000, 0.9500, 1.0000, 1.5000, 2.0000])

        y = np.array(
            [0.0185, 0.0751, 0.1336, 0.2147, 0.2997, 0.4118, 0.5310, 0.6582,
             0.7942, 0.9400, 0.9883, 0.9913, 0.9936, 0.9953, 0.9966, 0.9975,
             0.9983, 0.9988, 0.9992, 0.9996, 0.9998, 1.0000, 1.0008, 1.0014])
        kA_char2 = dc_cc(func=char_line(x, y), param='m')
        ev.set_attr(
            pr1=1, pr2=.999, ttd_l=5, design=['ttd_l'], offdesign=['kA_char'],
            kA_char1=kA_char1, kA_char2=kA_char2)

        # no kA modification for hot side!
        x = np.array([0, 1])
        y = np.array([1, 1])
        kA_char1 = dc_cc(func=char_line(x, y), param='m')

        # characteristic line for superheater kA
        x = np.array(
            [0, 0.045, 0.136, 0.244, 0.43, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2])
        y = np.array(
            [0, 0.037, 0.112, 0.207, 0.5, 0.8, 0.85, 0.9, 0.95, 1, 1.04, 1.07])
        kA_char2 = dc_cc(func=char_line(x, y), param='m')
        su.set_attr(kA_char1=kA_char1, kA_char2=kA_char2,
                    offdesign=['zeta1', 'zeta2', 'kA_char'])

        x = np.array(
            [0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5,
             0.5625, 0.6375, 0.7125, 0.7875, 0.9, 0.9875, 1, 1.0625, 1.125,
             1.175, 1.2125, 1.2375, 1.25])
        y = np.array(
            [0.0076, 0.1390, 0.2731, 0.4003, 0.5185, 0.6263, 0.7224, 0.8056,
             0.8754, 0.9312, 0.9729, 1.0006, 1.0203, 1.0158, 1.0051, 1.0000,
             0.9746, 0.9289, 0.8832, 0.8376, 0.7843, 0.7614])
        pu.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'],
                    eta_s_char=dc_cc(func=char_line(x, y), param='m'))

        # compressor system
        x = np.array([0, 0.4, 1, 1.2])
        y = np.array([0.5, 0.9, 1, 1.1])

        cp1.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'],
                     eta_s_char=dc_cc(func=char_line(x, y), param='m'))
        cp2.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'],
                     eta_s_char=dc_cc(func=char_line(x, y), param='m'))

        # characteristic line for intercooler kA
        x = np.linspace(0, 2.5, 26)
        y = np.array(
            [0.0000, 0.2455, 0.3747, 0.4798, 0.5718, 0.6552, 0.7323, 0.8045,
             0.8727, 0.9378, 1.0000, 1.0599, 1.1176, 1.1736, 1.2278, 1.2806,
             1.3320, 1.3822, 1.4313, 1.4792, 1.5263, 1.5724, 1.6176, 1.6621,
             1.7058, 1.7488])
        kA_char1 = dc_cc(func=char_line(x, y), param='m')

        x = np.linspace(0, 2.5, 26)
        y = np.array(
            [0.000, 0.164, 0.283, 0.389, 0.488, 0.581, 0.670, 0.756, 0.840,
             0.921, 1.000, 1.078, 1.154, 1.228, 1.302, 1.374, 1.446, 1.516,
             1.585, 1.654, 1.722, 1.789, 1.855, 1.921, 1.986, 2.051])
        kA_char2 = dc_cc(func=char_line(x, y), param='m')

        he.set_attr(kA_char1=kA_char1, kA_char2=kA_char2,
                    offdesign=['zeta1', 'zeta2', 'kA_char'])

        # characteristic line for condenser kA
        x = np.linspace(0, 2.5, 26)
        y = np.array(
            [0.0000, 0.2455, 0.3747, 0.4798, 0.5718, 0.6552, 0.7323, 0.8045,
             0.8727, 0.9378, 1.0000, 1.0599, 1.1176, 1.1736, 1.2278, 1.2806,
             1.3320, 1.3822, 1.4313, 1.4792, 1.5263, 1.5724, 1.6176, 1.6621,
             1.7058, 1.7488])
        kA_char1 = dc_cc(func=char_line(x, y), param='m')

        x = np.linspace(0, 2.5, 26)
        y = np.array(
            [0.000, 0.164, 0.283, 0.389, 0.488, 0.581, 0.670, 0.756, 0.840,
             0.921, 1.000, 1.078, 1.154, 1.228, 1.302, 1.374, 1.446, 1.516,
             1.585, 1.654, 1.722, 1.789, 1.855, 1.921, 1.986, 2.051])
        kA_char2 = dc_cc(func=char_line(x, y), param='m')

        cd.set_attr(kA_char1=kA_char1, kA_char2=kA_char2, pr2=0.9998,
                    design=['pr2'], offdesign=['zeta2', 'kA_char'])

        # %% connection parametrization
        # condenser system
        c_in_cd.set_attr(fluid={'water': 0, 'NH3': 1}, p=60)
        rp_cd.set_attr(T=60, fluid={'water': 1, 'NH3': 0}, p=10)
        self.cd_cons.set_attr(T=105)
        cons_cf.set_attr(h=ref(cb_rp, 1, 0), p=ref(cb_rp, 1, 0))
        cd_va.set_attr(
            p=ref(c_in_cd, 1, -1000), Td_bp=-5, h0=500, design=['Td_bp'])

        # evaporator system cold side
        pu_ev.set_attr(m=ref(va_dr, 10, 0), p0=5)
        dr_su.set_attr(p0=5, T=5)
        su_cp1.set_attr(
            p=ref(dr_su, 1, -5000), Td_bp=5, h0=1700, design=['Td_bp', 'p'])

        # evaporator system hot side
        self.amb_in_su.set_attr(m=20, T=12, p=1, fluid={'water': 1, 'NH3': 0})
        su_ev.set_attr(p=ref(self.amb_in_su, 1, -100), design=['p'])
        ev_amb_out.set_attr()

        # compressor-system
        cp1_he.set_attr(p=15)
        he_cp2.set_attr(T=40, p=ref(cp1_he, 1, -1000), design=['T', 'p'])
        ic_in_he.set_attr(p=1, T=20, m=5, fluid={'water': 1, 'NH3': 0})
        he_ic_out.set_attr(p=ref(ic_in_he, 1, -200), design=['p'])
        cp2_c_out.set_attr(p=ref(c_in_cd, 1, 0), h=ref(c_in_cd, 1, 0))

    def test_model(self):
        """
        Test the operating points of the heat pump against a different model.

        By now, not all characteristic functions of the original model are
        available in detail, thus perfect matching is not possible!
        """
        self.nw.solve('design')
        self.nw.save('tmp')
        self.nw.print_results()

        self.nw.print_results()

        # input values from ebsilon
        T = [105, 100, 90, 80]
        m_source = np.array([[23, 22, 20, 18, 16],
                             [27, 24, 20, 16, 12],
                             [31, 25, 20, 15, 10],
                             [33, 26, 20, 15, 10]])
        COP = np.array([[2.436, 2.414, 2.368, 2.338, 2.287],
                        [2.591, 2.523, 2.448, 2.355, 2.216],
                        [2.777, 2.635, 2.557, 2.442, 2.243],
                        [2.866, 2.711, 2.629, 2.528, 2.351]])

        i = 0
        for T in T:
            j = 0
            self.cd_cons.set_attr(T=T)
            for m in m_source[i]:
                self.amb_in_su.set_attr(m=m)
                if j == 0:
                    self.nw.solve(
                        'offdesign', design_path='tmp', init_path='tmp')

                else:
                    self.nw.solve('offdesign', design_path='tmp')

                # relative deviation should not exceed 20 %
                # this should be much less, unfortunately not all ebsilon
                # characteristics are available, thus it is
                # difficult/impossible to match the models perfectly! Another
                # issue is, that the component characteristics for generators
                # and motors do not work on the same basis (tespy: mechanical,
                # ebsilon: electrical)!
                d_rel_COP = abs(
                    self.heat.P.val / self.power.P.val - COP[i, j]) / COP[i, j]
                msg = ('The deviation in COP should be less than 0.065, is ' +
                       str(d_rel_COP) + ' at mass flow ' + str(m) +
                       ' and temperature ' + str(T) + '.')
                assert d_rel_COP < 0.065, msg
                j += 1
            i += 1
        shutil.rmtree('./tmp', ignore_errors=True)
