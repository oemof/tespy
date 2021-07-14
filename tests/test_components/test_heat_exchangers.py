# -*- coding: utf-8

"""Module for testing components of type heat exchanger.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_heat_exchangers.py

SPDX-License-Identifier: MIT
"""
import shutil

import numpy as np

from tespy.components import Condenser
from tespy.components import HeatExchanger
from tespy.components import HeatExchangerSimple
from tespy.components import ParabolicTrough
from tespy.components import Sink
from tespy.components import SolarCollector
from tespy.components import Source
from tespy.connections import Bus
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.fluid_properties import T_bp_p


def convergence_check(lin_dep):
    """Check convergence status of a simulation."""
    msg = 'Calculation did not converge!'
    assert lin_dep is False, msg


class TestHeatExchangers:

    def setup(self):

        self.nw = Network(
            ['H2O', 'Ar', 'INCOMP::S800'], T_unit='C', p_unit='bar',
            v_unit='m3 / s')
        self.inl1 = Source('inlet 1')
        self.outl1 = Sink('outlet 1')

    def setup_HeatExchangerSimple_network(self, instance):

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

    def test_HeatExchangerSimple(self):
        """Test component properties of simple heat exchanger."""
        instance = HeatExchangerSimple('heat exchanger')
        self.setup_HeatExchangerSimple_network(instance)
        fl = {'Ar': 0, 'H2O': 1, 'S800': 0}
        self.c1.set_attr(fluid=fl, m=1, p=10, T=100)
        # trigger heat exchanger parameter groups
        instance.set_attr(hydro_group='HW', L=100, ks=100, pr=0.99, Tamb=20)

        # test grouped parameter settings with missing parameters
        instance.hydro_group.is_set = True
        instance.kA_group.is_set = True
        instance.kA_char_group.is_set = True
        self.nw.solve('design', init_only=True)
        msg = ('Hydro group must no be set, if one parameter is missing!')
        assert instance.hydro_group.is_set is False, msg
        msg = ('kA group must no be set, if one parameter is missing!')
        assert instance.kA_group.is_set is False, msg
        msg = ('kA char group must no be set, if one parameter is missing!')
        assert instance.kA_char_group.is_set is False, msg

        # test diameter calculation from specified dimensions (as pipe)
        # with Hazen-Williams method
        instance.set_attr(hydro_group='HW', D='var', L=100,
                          ks=100, pr=0.99, Tamb=20)
        b = Bus('heat', P=-1e5)
        b.add_comps({'comp': instance})
        self.nw.add_busses(b)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        pr = round(self.c2.p.val_SI / self.c1.p.val_SI, 3)
        msg = ('Value of pressure ratio must be ' + str(pr) + ', is ' +
               str(instance.pr.val) + '.')
        assert pr == round(instance.pr.val, 3), msg

        # make zeta system variable and use previously calculated diameter
        # to calculate zeta. The value for zeta must not change
        zeta = round(instance.zeta.val, 0)
        instance.set_attr(D=instance.D.val, zeta='var', pr=np.nan)
        instance.D.is_var = False
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of zeta must be ' + str(zeta) + ', is ' +
               str(round(instance.zeta.val, 0)) + '.')
        assert zeta == round(instance.zeta.val, 0), msg

        # same test with pressure ratio as sytem variable
        pr = round(instance.pr.val, 3)
        instance.set_attr(zeta=np.nan, pr='var')
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of pressure ratio must be ' + str(pr) +
               ', is ' + str(round(instance.pr.val, 3)) + '.')
        assert pr == round(instance.pr.val, 3), msg

        # test heat transfer coefficient as variable of the system (ambient
        # temperature required)
        instance.set_attr(kA='var', pr=np.nan)
        b.set_attr(P=-5e4)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)

        # due to heat output being half of reference (for Tamb) kA should be
        # somewhere near to that (actual value is 677)
        msg = ('Value of heat transfer coefficient must be 677, is ' +
               str(instance.kA.val) + '.')
        assert 677 == round(instance.kA.val, 0), msg

        # test heat transfer as variable of the system
        instance.set_attr(Q='var', kA=np.nan)
        Q = -5e4
        b.set_attr(P=Q)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of heat transfer must be ' + str(Q) +
               ', is ' + str(instance.Q.val) + '.')
        assert Q == round(instance.Q.val, 0), msg

        # test kA as network results parameter
        instance.set_attr(Q=-5e4, Tamb=None)
        b.set_attr(P=None)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        kA_network = self.nw.results['HeatExchangerSimple'].loc[
            instance.label, 'kA']
        print(kA_network)
        msg = 'kA value must not be included in network results.'
        expr = not instance.kA.is_result and np.isnan(kA_network)
        assert expr, msg

        # test kA as network results parameter
        instance.set_attr(Tamb=20)
        self.nw.solve('design')
        kA_network = self.nw.results['HeatExchangerSimple'].loc[
            instance.label, 'kA']
        kA_comp = instance.kA.val
        msg = 'kA value needs to be identical on network and component level.'
        assert kA_network == kA_comp, msg

    def test_ParabolicTrough(self):
        """Test component properties of parabolic trough."""
        instance = ParabolicTrough('parabolic trough')
        self.setup_HeatExchangerSimple_network(instance)
        fl = {'Ar': 0, 'H2O': 0, 'S800': 1}
        self.c1.set_attr(fluid=fl, p=2, T=200)
        self.c2.set_attr(T=350)

        # test grouped parameter settings with missing parameters
        instance.hydro_group.is_set = True
        instance.energy_group.is_set = True
        self.nw.solve('design', init_only=True)
        msg = ('Hydro group must no be set, if one parameter is missing!')
        assert instance.hydro_group.is_set is False, msg
        msg = ('Energy group must no be set, if one parameter is missing!')
        assert instance.energy_group.is_set is False, msg

        # test solar collector params as system variables
        instance.set_attr(
            pr=1, aoi=10, doc=0.95, Q=1e6, Tamb=25, A='var', eta_opt=0.816,
            c_1=0.0622, c_2=0.00023, E=8e2, iam_1=-1.59e-3, iam_2=9.77e-5)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
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
        convergence_check(self.nw.lin_dep)
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: eta_opt
        instance.set_attr(E=5e2, eta_opt='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: c_1
        instance.set_attr(E=5e2, eta_opt=instance.eta_opt.val, c_1='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: c_2
        instance.set_attr(E=5e2, c_1=instance.c_1.val, c_2='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: iam_1
        instance.set_attr(E=5e2, c_2=instance.c_2.val, iam_1='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: iam_2
        instance.set_attr(E=5e2, iam_1=instance.iam_1.val, iam_2='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: aoi
        instance.set_attr(E=5e2, iam_2=instance.iam_2.val, aoi='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: doc
        instance.set_attr(E=5e2, aoi=instance.aoi.val, doc='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        assert Q_loss == round(instance.Q_loss.val, 0), msg

    def test_SolarCollector(self):
        """Test component properties of solar collector."""
        instance = SolarCollector('solar collector')
        self.setup_HeatExchangerSimple_network(instance)
        fl = {'Ar': 0, 'H2O': 1, 'S800': 0}
        self.c1.set_attr(fluid=fl, p=10, T=30)
        self.c2.set_attr(T=70)

        # test grouped parameter settings with missing parameters
        instance.hydro_group.is_set = True
        instance.energy_group.is_set = True
        self.nw.solve('design', init_only=True)
        msg = ('Hydro group must no be set, if one parameter is missing!')
        assert instance.hydro_group.is_set is False, msg
        msg = ('Energy group must no be set, if one parameter is missing!')
        assert instance.energy_group.is_set is False, msg

        # test solar collector params as system variables
        instance.set_attr(E=1e3, lkf_lin=1.0, lkf_quad=0.005, A='var',
                          eta_opt=0.9, Q=1e5, Tamb=20, pr=0.99)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
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
        convergence_check(self.nw.lin_dep)
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: eta_opt
        instance.set_attr(E=8e2, eta_opt='var')
        self.nw.solve('design')
        instance.set_attr(E=1e3)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: lkf_lin
        instance.set_attr(E=8e2, eta_opt=instance.eta_opt.val, lkf_lin='var')
        self.nw.solve('design')
        instance.set_attr(E=1e3)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: lkf_quad
        instance.set_attr(E=8e2, lkf_lin=instance.lkf_lin.val, lkf_quad='var')
        self.nw.solve('design')
        instance.set_attr(E=1e3)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: Tamb
        instance.set_attr(E=8e2, lkf_lin=instance.lkf_lin.val, lkf_quad='var')
        self.nw.solve('design')
        instance.set_attr(E=1e3)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        assert Q_loss == round(instance.Q_loss.val, 0), msg

    def test_HeatExchanger(self):
        """Test component properties of heat exchanger."""
        instance = HeatExchanger('heat exchanger')
        self.setup_HeatExchanger_network(instance)

        # design specification
        instance.set_attr(pr1=0.98, pr2=0.98, ttd_u=5,
                          design=['pr1', 'pr2', 'ttd_u'],
                          offdesign=['zeta1', 'zeta2', 'kA_char'])
        self.c1.set_attr(T=120, p=3, fluid={'Ar': 0, 'H2O': 1, 'S800': 0})
        self.c2.set_attr(T=70)
        self.c3.set_attr(T=40, p=5, fluid={'Ar': 1, 'H2O': 0, 'S800': 0})
        b = Bus('heat transfer', P=-80e3)
        b.add_comps({'comp': instance})
        self.nw.add_busses(b)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        self.nw.save('tmp')
        Q_design = instance.Q.val

        # test specified kA value
        instance.set_attr(kA=instance.kA.val * 2 / 3)
        b.set_attr(P=None)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)

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
        convergence_check(self.nw.lin_dep)

        # check heat transfer
        Q = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        td_log = ((self.c2.T.val - self.c3.T.val -
                   self.c1.T.val + self.c4.T.val) /
                  np.log((self.c2.T.val - self.c3.T.val) /
                         (self.c1.T.val - self.c4.T.val)))
        kA = round(-Q / td_log, 0)
        msg = ('Value of heat transfer must be ' + str(round(Q, 0)) + ', is ' +
               str(round(instance.Q.val, 0)) + '.')
        assert round(Q, 0) == round(instance.Q.val, 0), msg

        # check upper terminal temperature difference
        msg = ('Value of terminal temperature difference must be ' +
               str(round(instance.ttd_u.val, 1)) + ', is ' +
               str(round(self.c1.T.val - self.c4.T.val, 1)) + '.')
        ttd_u_calc = round(self.c1.T.val - self.c4.T.val, 1)
        ttd_u = round(instance.ttd_u.val, 1)
        assert ttd_u_calc == ttd_u, msg

        # check lower terminal temperature difference
        self.c2.set_attr(T=np.nan)
        instance.set_attr(ttd_l=20)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of terminal temperature difference must be ' +
               str(instance.ttd_l.val) + ', is ' +
               str(self.c2.T.val - self.c3.T.val) + '.')
        ttd_l_calc = round(self.c2.T.val - self.c3.T.val, 1)
        ttd_l = round(instance.ttd_l.val, 1)
        assert ttd_l_calc == ttd_l, msg

        # check specified kA value (by offdesign parameter), reset temperatures
        # to design state
        self.c2.set_attr(T=70)
        instance.set_attr(ttd_l=np.nan)
        self.nw.solve('offdesign', design_path='tmp')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of heat flow must be ' + str(instance.Q.val) + ', is ' +
               str(round(Q, 0)) + '.')
        assert round(Q, 0) == round(instance.Q.val, 0), msg
        msg = ('Value of heat transfer coefficient must be ' + str(kA) +
               ', is ' + str(round(instance.kA.val, 0)) + '.')
        assert kA == round(instance.kA.val, 0), msg

        # trigger negative lower terminal temperature difference as result
        self.c4.set_attr(T=np.nan)
        self.c2.set_attr(T=30)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of upper terminal temperature differences must be '
               'smaller than zero, is ' + str(round(instance.ttd_l.val, 1)) +
               '.')
        assert instance.ttd_l.val < 0, msg

        # trigger negative upper terminal temperature difference as result
        self.c4.set_attr(T=100)
        self.c2.set_attr(h=200e3, T=np.nan)
        instance.set_attr(pr1=0.98, pr2=0.98, ttd_u=np.nan,
                          design=['pr1', 'pr2'])
        self.c1.set_attr(h=150e3, T=np.nan)
        self.c3.set_attr(T=40)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of upper terminal temperature differences must be '
               'smaller than zero, is ' + str(round(instance.ttd_u.val, 1)) +
               '.')
        assert instance.ttd_u.val < 0, msg

        shutil.rmtree('./tmp', ignore_errors=True)

    def test_Condenser(self):
        """Test component properties of Condenser."""
        instance = Condenser('condenser')
        self.setup_HeatExchanger_network(instance)

        # design specification
        instance.set_attr(pr1=0.98, pr2=0.98, ttd_u=5,
                          offdesign=['zeta2', 'kA_char'])
        self.c1.set_attr(T=100, p0=0.5, fluid={'Ar': 0, 'H2O': 1, 'S800': 0})
        self.c3.set_attr(T=30, p=5, fluid={'Ar': 0, 'H2O': 1, 'S800': 0})
        self.c4.set_attr(T=40)
        instance.set_attr(Q=-80e3)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        self.nw.save('tmp')
        Q_design = instance.Q.val

        # test specified kA value
        instance.set_attr(kA=instance.kA.val * 2 / 3, Q=None)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)

        # test heat transfer
        Q = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        msg = (
            'Value of heat flow must be ' + str(round(Q_design * 2 / 3, 0)) +
            ', is ' + str(round(Q, 0)) + '.')
        assert round(Q, 1) == round(Q_design * 2 / 3, 1), msg

        # back to design case
        instance.set_attr(kA=None, Q=Q_design)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)

        # test heat transfer
        Q = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        msg = ('Value of heat flow must be ' + str(round(instance.Q.val, 0)) +
               ', is ' + str(round(Q, 0)) + '.')
        assert round(Q, 1) == round(instance.Q.val, 1), msg

        # test upper terminal temperature difference. For the component
        # condenser the temperature of the condensing fluid is relevant.
        ttd_u = round(T_bp_p(self.c1.get_flow()) - self.c4.T.val_SI, 1)
        p = round(self.c1.p.val_SI, 5)
        msg = ('Value of terminal temperature difference must be ' +
               str(round(instance.ttd_u.val, 1)) + ', is ' +
               str(ttd_u) + '.')
        assert ttd_u == round(instance.ttd_u.val, 1), msg

        # test lower terminal temperature difference
        instance.set_attr(ttd_l=20, ttd_u=np.nan, design=['pr2', 'ttd_l'])
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of terminal temperature difference must be ' +
               str(instance.ttd_l.val) + ', is ' +
               str(self.c2.T.val - self.c3.T.val) + '.')
        ttd_l_calc = round(self.c2.T.val - self.c3.T.val, 1)
        ttd_l = round(instance.ttd_l.val, 1)
        assert ttd_l_calc == ttd_l, msg

        # check kA value with condensing pressure in offdesign mode:
        # no changes to design point means: identical pressure
        self.nw.solve('offdesign', design_path='tmp')
        convergence_check(self.nw.lin_dep)
        msg = ('Value of condensing pressure be ' + str(p) + ', is ' +
               str(round(self.c1.p.val_SI, 5)) + '.')
        assert p == round(self.c1.p.val_SI, 5), msg
        shutil.rmtree('./tmp', ignore_errors=True)
