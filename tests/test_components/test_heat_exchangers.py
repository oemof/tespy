# -*- coding: utf-8

"""Module for testing components of type heat exchanger.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_heat_exchangers.py

SPDX-License-Identifier: MIT
"""
import math

import numpy as np

from tespy.components import Condenser
from tespy.components import HeatExchanger
from tespy.components import ParabolicTrough
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import SolarCollector
from tespy.components import Source
from tespy.connections import Bus
from tespy.connections import Connection
from tespy.networks import Network


class TestHeatExchangers:

    def setup_method(self):

        self.nw = Network(T_unit='C', p_unit='bar', v_unit='m3 / s')
        self.inl1 = Source('inlet 1')
        self.outl1 = Sink('outlet 1')

    def setup_SimpleHeatExchanger_network(self, instance):

        self.c1 = Connection(self.inl1, 'out1', instance, 'in1')
        self.c2 = Connection(instance, 'out1', self.outl1, 'in1')

        self.nw.add_conns(self.c1, self.c2)

    def setup_HeatExchanger_network(self, instance):

        self.inl2 = Source('inlet 2')
        self.outl2 = Sink('outlet 2')

        self.c1 = Connection(self.inl1, 'out1', instance, 'in1')
        self.c2 = Connection(instance, 'out1', self.outl1, 'in1')
        self.c3 = Connection(self.inl2, 'out1', instance, 'in2')
        self.c4 = Connection(instance, 'out2', self.outl2, 'in1')

        self.nw.add_conns(self.c1, self.c2, self.c3, self.c4)

    def test_SimpleHeatExchanger(self):
        """Test component properties of simple heat exchanger."""
        instance = SimpleHeatExchanger('heat exchanger')
        self.setup_SimpleHeatExchanger_network(instance)
        fl = {'H2O': 1}
        self.c1.set_attr(fluid=fl, m=1, p=10, T=100)
        # trigger heat exchanger parameter groups
        instance.set_attr(L=100, ks_HW=100, pr=0.99, Tamb=20)

        # test grouped parameter settings with missing parameters
        instance.darcy_group.is_set = True
        instance.kA_group.is_set = True
        instance.kA_char_group.is_set = True
        self.nw.solve('design', init_only=True)
        msg = ('Darcy group must no be set, if one parameter is missing!')
        assert not instance.darcy_group.is_set, msg
        msg = ('kA group must no be set, if one parameter is missing!')
        assert not instance.kA_group.is_set, msg
        msg = ('kA char group must no be set, if one parameter is missing!')
        assert not instance.kA_char_group.is_set, msg

        # test diameter calculation from specified dimensions (as pipe)
        # with Hazen-Williams method
        instance.set_attr(D='var', L=100, ks_HW=100, pr=0.99, Tamb=20)
        b = Bus('heat', P=-1e5)
        b.add_comps({'comp': instance})
        self.nw.add_busses(b)
        self.nw.solve('design')
        self.nw._convergence_check()
        pr = round(self.c2.p.val_SI / self.c1.p.val_SI, 3)
        msg = f"Value of pressure ratio must be {pr}, is {instance.pr.val}."
        assert pr == round(instance.pr.val, 3), msg

        # make zeta system variable and use previously calculated diameter
        # to calculate zeta. The value for zeta must not change
        zeta = round(instance.zeta.val, 0)
        instance.set_attr(D=instance.D.val, zeta='var', pr=None)
        instance.D.is_var = False
        self.nw.solve('design')
        self.nw._convergence_check()
        msg = f"Value of pressure ratio must be {zeta}, is {instance.zeta.val}."
        assert zeta == round(instance.zeta.val, 0), msg

        # same test with pressure ratio as sytem variable
        pr = round(instance.pr.val, 3)
        instance.set_attr(zeta=None, pr='var')
        self.nw.solve('design')
        self.nw._convergence_check()
        msg = f"Value of pressure ratio must be {pr}, is {instance.pr.val}."
        assert pr == round(instance.pr.val, 3), msg

        # test heat transfer coefficient as variable of the system (ambient
        # temperature required)
        instance.set_attr(kA='var', pr=None)
        b.set_attr(P=-5e4)
        self.nw.solve('design')
        self.nw._convergence_check()

        # due to heat output being half of reference (for Tamb) kA should be
        # somewhere near to that (actual value is 677)
        msg = (
            "Value of heat transfer coefficient must be 677, is "
            f"{instance.kA.val}."
        )
        assert 677 == round(instance.kA.val, 0), msg

        # test heat transfer as variable of the system
        instance.set_attr(Q='var', kA=None)
        Q = -5e4
        b.set_attr(P=Q)
        self.nw.solve('design')
        self.nw._convergence_check()
        msg = f"Value of heat transfer must be {Q}, is {instance.Q.val}."
        assert Q == round(instance.Q.val, 0), msg

        # test kA as network results parameter
        instance.set_attr(Q=-5e4, Tamb=None)
        b.set_attr(P=None)
        self.nw.solve('design')
        self.nw._convergence_check()
        kA_network = self.nw.results['SimpleHeatExchanger'].loc[
            instance.label, 'kA']
        print(kA_network)
        msg = 'kA value must not be included in network results.'
        expr = not instance.kA.is_result and np.isnan(kA_network)
        assert expr, msg

        # test kA as network results parameter
        instance.set_attr(Tamb=20)
        self.nw.solve('design')
        kA_network = self.nw.results['SimpleHeatExchanger'].loc[
            instance.label, 'kA']
        kA_comp = instance.kA.val
        msg = 'kA value needs to be identical on network and component level.'
        assert kA_network == kA_comp, msg

    def test_ParabolicTrough(self):
        """Test component properties of parabolic trough."""
        instance = ParabolicTrough('parabolic trough')
        self.setup_SimpleHeatExchanger_network(instance)
        self.c1.set_attr(fluid={'INCOMP::S800': 1}, p=2, T=200)
        self.c2.set_attr(T=350)

        # test grouped parameter settings with missing parameters
        instance.darcy_group.is_set = True
        instance.energy_group.is_set = True
        self.nw.solve('design', init_only=True)
        msg = ('Darcy group must no be set, if one parameter is missing!')
        assert not instance.darcy_group.is_set, msg
        msg = ('Energy group must no be set, if one parameter is missing!')
        assert not instance.energy_group.is_set, msg

        # test solar collector params as system variables
        instance.set_attr(
            pr=1, aoi=10, doc=0.95, Q=1e6, Tamb=25, A='var', eta_opt=0.816,
            c_1=0.0622, c_2=0.00023, E=8e2, iam_1=-1.59e-3, iam_2=9.77e-5)
        self.nw.solve('design')
        self.nw._convergence_check()
        # heat loss must be identical to E * A - Q (internal heat loss
        # calculation)
        T_diff = (self.c2.T.val + self.c1.T.val) / 2 - instance.Tamb.val
        iam = (
            1 - instance.iam_1.val * abs(instance.aoi.val) -
            instance.iam_2.val * instance.aoi.val ** 2)

        Q_loss = -round(instance.A.val * (
            instance.E.val * (
                1 - instance.eta_opt.val * instance.doc.val ** 1.5 * iam
            ) + T_diff * instance.c_1.val + T_diff ** 2 * instance.c_2.val), 0)
        msg = (
            'Value for heat loss of parabolic trough must be ' + str(Q_loss) +
            ', is ' + str(round(instance.Q_loss.val, 0)) + '.')
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: E
        # going to a different operating point first
        area = instance.A.val
        instance.set_attr(A=area * 1.2, E='var')
        self.nw.solve('design')
        instance.set_attr(A=area)
        self.nw.solve('design')
        self.nw._convergence_check()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: eta_opt
        instance.set_attr(E=5e2, eta_opt='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        self.nw._convergence_check()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: c_1
        instance.set_attr(E=5e2, eta_opt=instance.eta_opt.val, c_1='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        self.nw._convergence_check()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: c_2
        instance.set_attr(E=5e2, c_1=instance.c_1.val, c_2='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        self.nw._convergence_check()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: iam_1
        instance.set_attr(E=5e2, c_2=instance.c_2.val, iam_1='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        self.nw._convergence_check()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: iam_2
        instance.set_attr(E=5e2, iam_1=instance.iam_1.val, iam_2='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        self.nw._convergence_check()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: aoi
        instance.set_attr(E=5e2, iam_2=instance.iam_2.val, aoi='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        self.nw._convergence_check()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: doc
        instance.set_attr(E=5e2, aoi=instance.aoi.val, doc='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        self.nw._convergence_check()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

    def test_SolarCollector(self, tmp_path):
        """Test component properties of solar collector."""
        instance = SolarCollector('solar collector')
        self.setup_SimpleHeatExchanger_network(instance)
        fl = {'H2O': 1}
        self.c1.set_attr(fluid=fl, p=10, T=30)
        self.c2.set_attr(T=70)

        # test grouped parameter settings with missing parameters
        instance.darcy_group.is_set = True
        instance.energy_group.is_set = True
        self.nw.solve('design', init_only=True)
        msg = ('Darcy group must no be set, if one parameter is missing!')
        assert not instance.darcy_group.is_set, msg
        msg = ('Energy group must no be set, if one parameter is missing!')
        assert not instance.energy_group.is_set, msg

        # test solar collector params as system variables
        instance.set_attr(E=1e3, lkf_lin=1.0, lkf_quad=0.005, A='var',
                          eta_opt=0.9, Q=1e5, Tamb=20, pr=0.99)
        self.nw.solve('design')
        self.nw._convergence_check()
        # heat loss must be identical to E * A - Q (internal heat loss
        # calculation)
        T_diff = (self.c2.T.val + self.c1.T.val) / 2 - instance.Tamb.val
        Q_loss = -round(instance.A.val * (
            instance.E.val * (1 - instance.eta_opt.val) +
            T_diff * instance.lkf_lin.val +
            T_diff ** 2 * instance.lkf_quad.val), 0)
        msg = ('Value for heat loss of solar collector must be '
               + str(Q_loss) + ', is ' + str(round(instance.Q_loss.val, 0)) +
               '.')
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: E
        area = instance.A.val
        instance.set_attr(A=area * 1.2, E='var')
        self.nw.solve('design')
        instance.set_attr(A=area)
        self.nw.solve('design')
        self.nw._convergence_check()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: eta_opt
        instance.set_attr(E=8e2, eta_opt='var')
        self.nw.solve('design')
        instance.set_attr(E=1e3)
        self.nw.solve('design')
        self.nw._convergence_check()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: lkf_lin
        instance.set_attr(E=8e2, eta_opt=instance.eta_opt.val, lkf_lin='var')
        self.nw.solve('design')
        instance.set_attr(E=1e3)
        self.nw.solve('design')
        self.nw._convergence_check()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: lkf_quad
        instance.set_attr(E=8e2, lkf_lin=instance.lkf_lin.val, lkf_quad='var')
        self.nw.solve('design')
        instance.set_attr(E=1e3)
        self.nw.solve('design')
        self.nw._convergence_check()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: Tamb
        instance.set_attr(E=8e2, lkf_lin=instance.lkf_lin.val, lkf_quad='var')
        self.nw.solve('design')
        instance.set_attr(E=1e3)
        self.nw.solve('design')
        self.nw._convergence_check()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

    def test_HeatExchanger(self, tmp_path):
        """Test component properties of heat exchanger."""
        instance = HeatExchanger('heat exchanger')
        self.setup_HeatExchanger_network(instance)

        # design specification
        instance.set_attr(
            pr1=0.98, pr2=0.98, ttd_u=5,
            design=['pr1', 'pr2', 'ttd_u'],
            offdesign=['zeta1', 'zeta2', 'kA_char']
        )
        self.c1.set_attr(T=120, p=3, fluid={'H2O': 1})
        self.c2.set_attr(T=70)
        self.c3.set_attr(T=40, p=5, fluid={'Ar': 1})
        b = Bus('heat transfer', P=-80e3)
        b.add_comps({'comp': instance})
        self.nw.add_busses(b)
        self.nw.solve('design')
        self.nw._convergence_check()
        self.nw.save(tmp_path)
        Q_design = instance.Q.val

        # test specified kA value
        instance.set_attr(kA=instance.kA.val * 2 / 3)
        b.set_attr(P=None)
        self.nw.solve('design')
        self.nw._convergence_check()

        # test heat transfer
        Q = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        msg = (
            'Value of heat flow must be ' + str(round(Q_design * 2 / 3, 0)) +
            ', is ' + str(round(Q, 0)) + '.')
        assert round(Q, 1) == round(Q_design * 2 / 3, 1), msg

        # back to design case
        instance.set_attr(kA=None)
        b.set_attr(P=Q_design)
        self.nw.solve('design')
        self.nw._convergence_check()

        # check heat transfer
        Q = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        td_log = (
            (self.c2.T.val - self.c3.T.val - self.c1.T.val + self.c4.T.val)
            / math.log(
                (self.c2.T.val - self.c3.T.val)
                / (self.c1.T.val - self.c4.T.val)
            )
        )
        kA = round(-Q / td_log, 0)
        msg = (
            f"Value of heat transfer must be {round(Q, 0)}, is "
            f"{round(instance.Q.val, 0)}."
        )
        assert round(Q, 0) == round(instance.Q.val, 0), msg

        # check upper terminal temperature difference
        msg = ('Value of terminal temperature difference must be ' +
               str(round(instance.ttd_u.val, 1)) + ', is ' +
               str(round(self.c1.T.val - self.c4.T.val, 1)) + '.')
        ttd_u_calc = round(self.c1.T.val - self.c4.T.val, 1)
        ttd_u = round(instance.ttd_u.val, 1)
        assert ttd_u_calc == ttd_u, msg

        # check lower terminal temperature difference
        self.c2.set_attr(T=None)
        instance.set_attr(ttd_l=20)
        self.nw.solve('design')
        self.nw._convergence_check()
        msg = ('Value of terminal temperature difference must be ' +
               str(instance.ttd_l.val) + ', is ' +
               str(self.c2.T.val - self.c3.T.val) + '.')
        ttd_l_calc = round(self.c2.T.val - self.c3.T.val, 1)
        ttd_l = round(instance.ttd_l.val, 1)
        assert ttd_l_calc == ttd_l, msg

        # check specified kA value (by offdesign parameter), reset temperatures
        # to design state
        self.c2.set_attr(T=70)
        instance.set_attr(ttd_l=None)
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw._convergence_check()
        msg = ('Value of heat flow must be ' + str(instance.Q.val) + ', is ' +
               str(round(Q, 0)) + '.')
        assert round(Q, 0) == round(instance.Q.val, 0), msg
        msg = ('Value of heat transfer coefficient must be ' + str(kA) +
               ', is ' + str(round(instance.kA.val, 0)) + '.')
        assert kA == round(instance.kA.val, 0), msg

        # trigger negative lower terminal temperature difference as result
        self.c4.set_attr(T=None)
        self.c2.set_attr(T=30)
        self.nw.solve('design')
        self.nw._convergence_check()
        msg = ('Value of upper terminal temperature differences must be '
               'smaller than zero, is ' + str(round(instance.ttd_l.val, 1)) +
               '.')
        assert instance.ttd_l.val < 0, msg

        # trigger negative upper terminal temperature difference as result
        self.c4.set_attr(T=100)
        self.c2.set_attr(h=200e3, T=None)
        instance.set_attr(
            pr1=0.98, pr2=0.98, ttd_u=None, design=['pr1', 'pr2']
        )
        self.c1.set_attr(h=150e3, T=None)
        self.c3.set_attr(T=40)
        self.nw.solve('design')
        self.nw._convergence_check()
        msg = (
            'Value of upper terminal temperature differences must be '
            f'smaller than zero, is {round(instance.ttd_u.val, 1)}.'
        )
        assert instance.ttd_u.val < 0, msg

        # test heat exchanger effectiveness
        b.set_attr(P=None)
        self.c1.set_attr(m=1, T=80, h=None)
        self.c2.set_attr(T=None, h=None)

        self.c3.set_attr(m=4, T=25)
        self.c4.set_attr(T=None)

        instance.set_attr(eff_max=0.9, ttd_u=None)
        self.nw.solve("design")
        self.nw._convergence_check()
        msg = (
            'Value of cold side effectiveness must be equal to 0.9 but is '
            f'{round(instance.eff_cold.val, 1)}.'
        )
        assert round(instance.eff_cold.val, 1) == 0.9, msg

        self.c3.set_attr(m=None)
        instance.set_attr(eff_max=None, eff_hot=0.9, eff_cold=0.9)
        self.nw.solve("design")
        self.nw._convergence_check()
        msg = (
            'Value of max effectiveness must be equal to 0.9 but is '
            f'{round(instance.eff_max.val, 1)}.'
        )
        assert round(instance.eff_max.val, 1) == 0.9, msg

        msg = (
            'Value of cold side mass flow must be equal to 7.96 but is '
            f'{round(self.c3.m.val, 2)}.'
        )
        assert round(self.c3.m.val, 2) == 7.96

    def test_HeatExchanger_effectiveness_invalid(self):

        instance = HeatExchanger('heat exchanger')
        self.setup_HeatExchanger_network(instance)

        # remove fluid specifications
        self.c1.set_attr(fluid={f: None for f in self.c1.fluid.val})
        self.c3.set_attr(fluid={f: None for f in self.c3.fluid.val})

        # add new fluids
        # temperature range > 300 °C
        self.c1.set_attr(fluid={"INCOMP::NaK": 1}, m=10, T=400, p=1)
        self.c2.set_attr(T=350, p=1)
        # temperature range < 100 °C at 1 bar
        self.c3.set_attr(fluid={"INCOMP::Water": 1}, T=25, p=1)
        self.c4.set_attr(T=50, p=1)
        instance.set_attr(eff_cold=None, eff_hot=None, pr1=None, pr2=None)

        self.nw.solve("design")
        self.nw._convergence_check()
        msg = (
            'Value of cold effectiveness must be nan but is '
            f'{round(instance.eff_cold.val, 1)}.'
        )
        assert np.isnan(instance.eff_cold.val), msg
        msg = (
            'Value of hot effectiveness must be nan but is '
            f'{round(instance.eff_hot.val, 1)}.'
        )
        assert np.isnan(instance.eff_hot.val), msg

    def test_Condenser(self, tmp_path):
        """Test component properties of Condenser."""
        instance = Condenser('condenser')
        self.setup_HeatExchanger_network(instance)

        # design specification
        instance.set_attr(
            pr1=0.98, pr2=0.98, ttd_u=5, offdesign=['zeta2', 'kA_char']
        )
        self.c1.set_attr(T=100, p0=0.5, fluid={'H2O': 1})
        self.c3.set_attr(T=30, p=5, fluid={'H2O': 1})
        self.c4.set_attr(T=40)
        instance.set_attr(Q=-80e3)
        self.nw.solve('design')
        self.nw._convergence_check()
        self.nw.save(tmp_path)
        Q_design = instance.Q.val

        # test specified kA value
        instance.set_attr(kA=instance.kA.val * 2 / 3, Q=None)
        self.nw.solve('design')
        self.nw._convergence_check()

        # test heat transfer
        Q = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        msg = (
            'Value of heat flow must be ' + str(round(Q_design * 2 / 3, 0)) +
            ', is ' + str(round(Q, 0)) + '.')
        assert round(Q, 1) == round(Q_design * 2 / 3, 1), msg

        # back to design case
        instance.set_attr(kA=None, Q=Q_design)
        self.nw.solve('design')
        self.nw._convergence_check()

        # test heat transfer
        Q = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        msg = ('Value of heat flow must be ' + str(round(instance.Q.val, 0)) +
               ', is ' + str(round(Q, 0)) + '.')
        assert round(Q, 1) == round(instance.Q.val, 1), msg

        # test upper terminal temperature difference. For the component
        # condenser the temperature of the condensing fluid is relevant.
        ttd_u = round(self.c1.calc_T_sat() - self.c4.T.val_SI, 1)
        p = round(self.c1.p.val_SI, 5)
        msg = ('Value of terminal temperature difference must be ' +
               str(round(instance.ttd_u.val, 1)) + ', is ' +
               str(ttd_u) + '.')
        assert ttd_u == round(instance.ttd_u.val, 1), msg

        # test lower terminal temperature difference
        instance.set_attr(ttd_l=20, ttd_u=None, design=['pr2', 'ttd_l'])
        self.nw.solve('design')
        self.nw._convergence_check()
        msg = ('Value of terminal temperature difference must be ' +
               str(instance.ttd_l.val) + ', is ' +
               str(self.c2.T.val - self.c3.T.val) + '.')
        ttd_l_calc = round(self.c2.T.val - self.c3.T.val, 1)
        ttd_l = round(instance.ttd_l.val, 1)
        assert ttd_l_calc == ttd_l, msg

        # check kA value with condensing pressure in offdesign mode:
        # no changes to design point means: identical pressure
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw._convergence_check()
        msg = ('Value of condensing pressure be ' + str(p) + ', is ' +
               str(round(self.c1.p.val_SI, 5)) + '.')
        assert p == round(self.c1.p.val_SI, 5), msg

    def test_CondenserWithEvaporation(self):
        """Test a Condenser that evaporates a fluid."""
        instance = Condenser('condenser')
        self.setup_HeatExchanger_network(instance)

        # design specification
        instance.set_attr(pr1=1, pr2=1, offdesign=["kA"])
        self.c1.set_attr(x=1, p=1, fluid={'H2O': 1}, m=1)
        self.c3.set_attr(x=0, p=0.7, fluid={'H2O': 1}, m=2, design=["m"])
        self.nw.solve('design')
        self.nw._convergence_check()
        ttd_l = round(instance.ttd_l.val, 3)
        ttd_u = round(instance.ttd_u.val, 3)
        td_log = round(instance.td_log.val, 3)
        msg = (
            "Value of upper and lower terminal temperature differences must be "
            f"identical, but they are not: ttd_u={ttd_u}, ttd_l={ttd_l}."
        )
        assert instance.ttd_l.val == instance.ttd_u.val, msg

        msg = (
            "Value of logarithmic and lower terminal temperature differences "
            f"must be identical, but they are not: td_log={td_log}, "
            f"ttd_l={ttd_l}."
        )
        assert instance.td_log.val == instance.ttd_l.val, msg

        # self.nw.save(tmp_path)
        # self.c1.set_attr(m=1)
        # self.nw.solve("offdesign", design_path="tmp")
        # self.nw._convergence_check()
        # msg = (
        #     "Value of logarithmic and lower terminal temperature differences "
        #     f"must be identical, but they are not: td_log={td_log}, "
        #     f"ttd_l={self.c3.m.val}."
        # )
        # assert instance.td_log.val == instance.ttd_l.val, msg
