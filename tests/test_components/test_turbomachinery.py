# -*- coding: utf-8

"""Module for testing components of type turbomachinery.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_turbomachinery.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components import Compressor
from tespy.components import Pump
from tespy.components import Sink
from tespy.components import Source
from tespy.components import SteamTurbine
from tespy.components import Turbine
from tespy.components.turbomachinery.base import Turbomachine
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap
from tespy.tools.characteristics import load_default_char as ldc
from tespy.tools.fluid_properties import isentropic
from tespy.tools.fluid_properties import s_mix_ph


class TestTurbomachinery:

    def setup_network(self, instance):
        self.nw = Network(T_unit='C', p_unit='bar', v_unit='m3 / s')
        self.source = Source('source')
        self.sink = Sink('sink')
        self.c1 = Connection(self.source, 'out1', instance, 'in1')
        self.c2 = Connection(instance, 'out1', self.sink, 'in1')
        self.nw.add_conns(self.c1, self.c2)

    def test_Compressor(self, tmp_path):
        """Test component properties of compressors."""
        instance = Compressor('compressor')
        self.setup_network(instance)

        # compress NH3, other fluids in network are for turbine, pump, ...
        fl = {'N2': 1}
        self.c1.set_attr(fluid=fl, v=1, p=1, T=5)
        self.c2.set_attr(p=6)
        instance.set_attr(eta_s=0.8)
        self.nw.solve('design')
        self.nw._convergence_check()
        self.nw.save(tmp_path)

        # test isentropic efficiency value
        eta_s_d = (
            (isentropic(self.c1.p.val_SI, self.c1.h.val_SI, self.c2.p.val_SI, self.c1.fluid_data, self.c1.mixing_rule) -
             self.c1.h.val_SI) / (self.c2.h.val_SI - self.c1.h.val_SI))
        msg = ('Value of isentropic efficiency must be ' + str(eta_s_d) +
               ', is ' + str(instance.eta_s.val) + '.')
        assert round(eta_s_d, 3) == round(instance.eta_s.val, 3), msg

        # trigger invalid value for isentropic efficiency
        instance.set_attr(eta_s=1.1)
        self.nw.solve('design')
        self.nw._convergence_check()

        # test calculated value
        eta_s = (
            (isentropic(self.c1.p.val_SI, self.c1.h.val_SI, self.c2.p.val_SI, self.c1.fluid_data, self.c1.mixing_rule) -
             self.c1.h.val_SI) / (self.c2.h.val_SI - self.c1.h.val_SI))
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
               ', is ' + str(instance.eta_s.val) + '.')
        assert round(eta_s, 3) == round(instance.eta_s.val, 3), msg

        # remove pressure at outlet, use characteristic map for pressure
        # rise calculation
        self.c2.set_attr(p=None)
        instance.set_attr(
            char_map_pr={'char_func': ldc(
                'compressor', 'char_map_pr', 'DEFAULT', CharMap),
                'is_set': True},
            char_map_eta_s={'char_func': ldc(
                'compressor', 'char_map_eta_s', 'DEFAULT', CharMap),
                'is_set': True},
            eta_s=None, igva=0)

        # offdesign test, efficiency value should be at design value
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw._convergence_check()
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
               ') must be identical to design case (' + str(eta_s) + ').')
        assert round(eta_s_d, 2) == round(instance.eta_s.val, 2), msg

        # move to highest available speedline, mass flow below lowest value
        # at that line
        self.c1.set_attr(v=None, m=self.c1.m.val * 0.8, T=-30)
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw._convergence_check()

        # should be value
        eta_s = eta_s_d * instance.char_map_eta_s.char_func.z[6, 0]
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
               ') must be at (' + str(round(eta_s, 4)) + ').')
        assert round(eta_s, 4) == round(instance.eta_s.val, 4), msg

        # going below lowest available speedline, above highest mass flow at
        # that line
        self.c1.set_attr(T=175)
        self.nw.solve('offdesign', design_path=tmp_path, init_path=tmp_path)
        self.nw._convergence_check()
        # should be value
        eta_s = eta_s_d * instance.char_map_eta_s.char_func.z[0, 9]
        msg = ('Value of isentropic efficiency (' + str(instance.eta_s.val) +
               ') must be at (' + str(round(eta_s, 4)) + ').')
        assert round(eta_s, 4) == round(instance.eta_s.val, 4), msg

        # back to design properties, test eta_s_char
        self.c2.set_attr(p=6)
        self.c1.set_attr(v=1, T=5, m=None)

        # test parameter specification for eta_s_char with unset char map
        instance.set_attr(eta_s_char={'char_func': ldc(
            'compressor', 'eta_s_char', 'DEFAULT', CharLine),
            'is_set': True, 'param': 'm'})
        instance.char_map_eta_s.is_set = False
        instance.char_map_pr.is_set = False
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw._convergence_check()
        msg = ('Value of isentropic efficiency must be ' + str(eta_s_d) +
               ', is ' + str(instance.eta_s.val) + '.')
        assert round(eta_s_d, 3) == round(instance.eta_s.val, 3), msg

        # move up in volumetric flow
        self.c1.set_attr(v=1.5)
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw._convergence_check()
        eta_s = round(eta_s_d * instance.eta_s_char.char_func.evaluate(
            self.c1.m.val_SI / self.c1.m.design), 3)
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
               ', is ' + str(round(instance.eta_s.val, 3)) + '.')
        assert eta_s == round(instance.eta_s.val, 3), msg

        # test parameter specification for pr
        instance.eta_s_char.set_attr(param='pr')
        self.c1.set_attr(v=1)
        self.c2.set_attr(p=6)
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw._convergence_check()
        expr = (self.c2.p.val_SI * self.c1.p.design /
                (self.c2.p.design * self.c1.p.val_SI))
        eta_s = round(
            eta_s_d * instance.eta_s_char.char_func.evaluate(expr), 3)
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
               ', is ' + str(round(instance.eta_s.val, 3)) + '.')
        assert eta_s == round(instance.eta_s.val, 3), msg

    def test_Pump(self, tmp_path):
        """Test component properties of pumps."""
        instance = Pump('pump')
        self.setup_network(instance)
        self.c1.set_attr(fluid={'INCOMP::DowQ': 1}, v=1, p=5, T=50)
        self.c2.set_attr(p=7)
        instance.set_attr(eta_s=1)
        self.nw.solve('design')
        self.nw._convergence_check()

        # test calculated value for efficiency
        eta_s = (
            (isentropic(self.c1.p.val_SI, self.c1.h.val_SI, self.c2.p.val_SI, self.c1.fluid_data, self.c1.mixing_rule) -
             self.c1.h.val_SI) / (self.c2.h.val_SI - self.c1.h.val_SI))
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
               ', is ' + str(instance.eta_s.val) + '.')
        assert eta_s == instance.eta_s.val, msg

        # isentropic efficiency of 1 means inlet and outlet entropy are
        # identical
        s1 = round(self.c1.calc_s(), 4)
        s2 = round(self.c2.calc_s(), 4)
        msg = ('Value of entropy must be identical for inlet (' + str(s1) +
               ') and outlet (' + str(s2) +
               ') at 100 % isentropic efficiency.')
        assert s1 == s2, msg

        # specify realistic value for efficiency, outlet pressure from flow
        # char
        eta_s_d = 0.8
        instance.set_attr(eta_s=eta_s_d)
        self.nw.solve('design')
        self.nw._convergence_check()
        self.nw.save(tmp_path)
        self.c2.set_attr(p=None)

        # flow char (pressure rise vs. volumetric flow)
        x = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4]
        y = np.array([14, 13.5, 12.5, 11, 9, 6.5, 3.5, 0]) * 1e5
        char = {'char_func': CharLine(x, y), 'is_set': True}
        # apply flow char and eta_s char
        instance.set_attr(
            flow_char=char, eta_s=None, eta_s_char={
                'char_func': ldc('pump', 'eta_s_char', 'DEFAULT', CharLine),
                'is_set': True})
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw._convergence_check()

        # value for difference pressure
        dp = 650000.0
        msg = ('Value of pressure rise must be ' + str(dp) + ', is ' +
               str(self.c2.p.val_SI - self.c1.p.val_SI) + '.')
        assert round(self.c2.p.val_SI - self.c1.p.val_SI, 0) == dp, msg

        # test ohter volumetric flow on flow char
        self.c1.set_attr(v=0.9)
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw._convergence_check()
        dp = 775000.0
        msg = ('Value of pressure rise must be ' + str(dp) + ', is ' +
               str(round(self.c2.p.val_SI - self.c1.p.val_SI, 0)) + '.')
        assert round(self.c2.p.val_SI - self.c1.p.val_SI, 0) == dp, msg

        # test value of isentropic efficiency
        eta_s = round(eta_s_d * instance.eta_s_char.char_func.evaluate(
            self.c1.v.val_SI / self.c1.v.design), 3)
        msg = ('Value of isentropic efficiency must be ' + str(eta_s) +
               ', is ' + str(instance.eta_s.val) + '.')
        assert eta_s == round(instance.eta_s.val, 3), msg
        instance.eta_s_char.is_set = False

        # test boundaries of characteristic line:
        # lower boundary
        instance.set_attr(eta_s=0.8)
        self.c1.set_attr(m=0, v=None)
        self.nw.solve('design')
        self.nw._convergence_check()
        msg = ('Value of power must be ' + str(14e5) + ', is ' +
               str(round(self.c2.p.val_SI - self.c1.p.val_SI, 0)) + '.')
        assert round(self.c2.p.val_SI - self.c1.p.val_SI, 0) == 14e5, msg

        # upper boundary
        self.c1.set_attr(v=1.5, m=None)
        self.nw.solve('design')
        self.nw._convergence_check()
        msg = ('Value of power must be 0, is ' +
               str(round(self.c2.p.val_SI - self.c1.p.val_SI, 0)) + '.')
        assert round(self.c2.p.val_SI - self.c1.p.val_SI, 0) == 0, msg

    def test_Turbine(self, tmp_path):
        """Test component properties of turbines."""
        instance = Turbine('turbine')
        self.setup_network(instance)
        fl = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129}
        self.c1.set_attr(fluid=fl, m=15, p=10)
        self.c2.set_attr(p=1, T=25)
        instance.set_attr(eta_s=0.85)
        self.nw.solve('design')
        self.nw._convergence_check()
        self.nw.save(tmp_path)

        # design value of isentropic efficiency
        eta_s_d = (
            (self.c2.h.val_SI - self.c1.h.val_SI)
            / (
                isentropic(
                    self.c1.p.val_SI,
                    self.c1.h.val_SI,
                    self.c2.p.val_SI,
                    self.c1.fluid_data,
                    self.c1.mixing_rule
                )
               - self.c1.h.val_SI
            )
        )
        msg = (
            f'Value of isentropic efficiency must be {round(eta_s_d, 3)}, is '
            f'{instance.eta_s.val}.'
        )
        assert round(eta_s_d, 3) == round(instance.eta_s.val, 3), msg

        # trigger invalid value for isentropic efficiency
        instance.set_attr(eta_s=1.1)
        self.nw.solve('design')
        self.nw._convergence_check()
        eta_s = (
            (self.c2.h.val_SI - self.c1.h.val_SI)
            / (
                isentropic(
                    self.c1.p.val_SI,
                    self.c1.h.val_SI,
                    self.c2.p.val_SI,
                    self.c1.fluid_data,
                    self.c1.mixing_rule
                )
               - self.c1.h.val_SI
            )
        )
        msg = (
            f'Value of isentropic efficiency must be {round(eta_s, 3)}, is '
            f'{instance.eta_s.val}.'
        )
        assert round(eta_s, 3) == round(instance.eta_s.val, 3), msg

        # unset isentropic efficiency and inlet pressure,
        # use characteristcs and cone law instead, parameters have to be in
        # design state
        self.c1.set_attr(p=None)
        instance.cone.is_set = True
        instance.eta_s_char.is_set = True
        instance.eta_s.is_set = False
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw._convergence_check()
        # check efficiency
        msg = (
            f'Value of isentropic efficiency ({instance.eta_s.val}) must be '
            f'identical to design case ({eta_s_d}).'
        )
        assert round(eta_s_d, 2) == round(instance.eta_s.val, 2), msg
        # check pressure
        msg = (
            f'Value of inlet pressure ({round(self.c1.p.val_SI)}) must be '
            f'identical to design case ({round(self.c1.p.design)}).'
        )
        assert round(self.c1.p.design) == round(self.c1.p.val_SI), msg

        # lowering mass flow, inlet pressure must sink according to cone law
        self.c1.set_attr(m=self.c1.m.val * 0.8)
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw._convergence_check()
        msg = f'Value of pressure ratio ({instance.pr.val}) must be at 0.128.'
        assert 0.128 == round(instance.pr.val, 3), msg

        # testing more parameters for eta_s_char
        # test parameter specification v
        self.c1.set_attr(m=10)
        instance.eta_s_char.param = 'v'
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw._convergence_check()
        expr = self.c1.v.val_SI / self.c1.v.design
        eta_s = round(eta_s_d * instance.eta_s_char.char_func.evaluate(expr), 3)
        msg = (
            f'Value of isentropic efficiency ({round(instance.eta_s.val, 3)}) '
            f'must be {eta_s}.'
        )
        assert eta_s == round(instance.eta_s.val, 3), msg

        # test parameter specification pr
        instance.eta_s_char.param = 'pr'
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw._convergence_check()
        expr = (
            self.c2.p.val_SI * self.c1.p.design
            / (self.c2.p.design * self.c1.p.val_SI)
        )
        eta_s = round(eta_s_d * instance.eta_s_char.char_func.evaluate(expr), 3)
        msg = (
            f'Value of isentropic efficiency ({round(instance.eta_s.val, 3)}) '
            f'must be {eta_s}.'
        )
        assert eta_s == round(instance.eta_s.val, 3), msg

    def test_SteamTurbine(self, tmp_path):
        instance = SteamTurbine('turbine')
        self.setup_network(instance)
        fl = {'H2O': 1}
        # start in gas, end in gas
        self.c1.set_attr(fluid=fl, m=15, p=100, T=500)
        self.c2.set_attr(Td_bp=1)

        eta_s_dry = 0.9
        instance.set_attr(eta_s_dry=0.9, alpha=1)

        self.nw.solve('design')
        self.nw._convergence_check()

        eta_s = round(instance.eta_s.val, 4)
        msg = (
            f"In gas expansion isentropic ({eta_s}) and dry "
            f"efficiency ({eta_s_dry}) have to be identical."
        )
        assert eta_s_dry == eta_s, msg

        # end in two phase
        self.c2.set_attr(Td_bp=None, x=0.9)
        self.nw.solve('design')
        self.nw._convergence_check()

        eta_s = round(instance.eta_s.val, 4)
        msg = (
            f"In gas expansion isentropic ({eta_s}) and dry "
            f"efficiency ({eta_s_dry}) cannot be identical."
        )
        assert eta_s_dry != eta_s, msg

        # start in two phase
        self.c1.set_attr(T=None, x=0.99)
        self.nw.solve('design')
        self.nw._convergence_check()

        eta_s = round(instance.eta_s.val, 4)
        msg = (
            f"In gas expansion isentropic ({eta_s}) and dry "
            f"efficiency ({eta_s_dry}) cannot be identical."
        )
        assert eta_s_dry != eta_s, msg

    def test_Turbomachine(self):
        """Test component properties of turbomachines."""
        instance = Turbomachine('turbomachine')
        msg = f'Component name must be turbomachine, is {instance.component()}.'
        assert 'turbomachine' == instance.component(), msg
        self.setup_network(instance)
        fl = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129}
        self.c1.set_attr(fluid=fl, m=10, p=1, h=1e5)
        self.c2.set_attr(p=3, h=2e5)

        # pressure ratio and power are the basic functions for turbomachines,
        # these are inherited by all children, thus only tested here
        self.nw.solve('design')
        self.nw._convergence_check()
        power = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        pr = self.c2.p.val_SI / self.c1.p.val_SI
        msg = ('Value of power must be ' + str(power) + ', is ' +
               str(instance.P.val) + '.')
        assert power == instance.P.val, msg
        msg = ('Value of power must be ' + str(pr) + ', is ' +
               str(instance.pr.val) + '.')
        assert pr == instance.pr.val, msg

        # test pressure ratio specification
        self.c2.set_attr(p=None)
        instance.set_attr(pr=5)
        self.nw.solve('design')
        self.nw._convergence_check()
        pr = self.c2.p.val_SI / self.c1.p.val_SI
        msg = ('Value of power must be ' + str(pr) + ', is ' +
               str(instance.pr.val) + '.')
        assert pr == instance.pr.val, msg

        # test power specification
        self.c2.set_attr(h=None)
        instance.set_attr(P=1e5)
        self.nw.solve('design')
        self.nw._convergence_check()
        power = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        msg = ('Value of power must be ' + str(power) + ', is ' +
               str(instance.P.val) + '.')
        assert power == instance.P.val, msg
