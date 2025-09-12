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
from pytest import approx
from pytest import fixture
from pytest import mark

from tespy.components import Condenser
from tespy.components import Desuperheater
from tespy.components import HeatExchanger
from tespy.components import ParallelFlowHeatExchanger
from tespy.components import ParabolicTrough
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import SolarCollector
from tespy.components import Source
from tespy.connections import Bus
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.fluid_properties import h_mix_pT


def _calc_Q(c1, c2):
    return c1.m.val_SI * (c2.h.val_SI - c1.h.val_SI)


def _calc_zeta(c1, c2):
    return (
        (c1.p.val_SI - c2.p.val_SI) * math.pi ** 2
        / (4 * c1.m.val_SI ** 2 * (c1.vol.val_SI + c2.vol.val_SI))
    )


def _calc_pr(c1, c2):
    return c2.p.val_SI / c1.p.val_SI


def _calc_dp(c1, c2):
    return c1.p.val_SI - c2.p.val_SI


def _calc_ttd_l(c2, c3):
    return c2.T.val_SI - c3.T.val_SI


def _calc_ttd_u(c1, c4):
    return c1.T.val_SI - c4.T.val_SI


def _calc_ttd_u_condenser(c1, c4):
    return c1.calc_T_sat() - c4.T.val_SI


def _calc_td_log(ttd_u, ttd_l):
    return (ttd_l - ttd_u) / math.log(ttd_l / ttd_u)


def _calc_kA(Q, td_log):
    return - Q / td_log


def _calc_eff_hot(c1, c2, c3):
    h_at_T_in_cold = h_mix_pT(
        c2.p.val_SI, c3.T.val_SI, c2.fluid_data, c2.mixing_rule
    )
    return (c2.h.val_SI - c1.h.val_SI) / (h_at_T_in_cold - c1.h.val_SI)


def _calc_eff_cold(c3, c4, c1):
    h_at_T_in_hot = h_mix_pT(
        c4.p.val_SI, c1.T.val_SI, c4.fluid_data, c4.mixing_rule
    )
    return (c4.h.val_SI - c3.h.val_SI) / (h_at_T_in_hot - c3.h.val_SI)


@fixture
def heatexchanger_network(request):

    nw = Network()
    nw.units.set_defaults(**{
        "pressure": "bar", "temperature": "degC", "volumetric_flow": "m3/s"
    })

    inl1 = Source('inlet 1')
    outl1 = Sink('outlet 1')
    inl2 = Source('inlet 2')
    outl2 = Sink('outlet 2')
    instance = request.param('heat exchanger')

    c1 = Connection(inl1, 'out1', instance, 'in1', label='1')
    c2 = Connection(instance, 'out1', outl1, 'in1', label='2')
    c3 = Connection(inl2, 'out1', instance, 'in2', label='3')
    c4 = Connection(instance, 'out2', outl2, 'in1', label='4')

    nw.add_conns(c1, c2, c3, c4)

    return nw


@fixture
def heatexchanger_converged_network(heatexchanger_network):
    # rename for shorter variable name
    nw = heatexchanger_network

    c1, c2, c3, c4 = nw.get_conn(['1', '2', '3', '4'])
    instance = nw.get_comp('heat exchanger')

    c1.set_attr(T=120, p=1, fluid={'H2O': 1}, m=1)
    c2.set_attr(p=0.95)
    c3.set_attr(T=40, p=5, fluid={'Ar': 1})
    c4.set_attr(T=90, p=4.9)

    if type(instance) == HeatExchanger:
        c2.set_attr(T=105)

    nw.solve('design')
    nw.assert_convergence()

    return nw


class TestHeatExchanger:

    @mark.parametrize(
        'heatexchanger_network',
        [HeatExchanger, Condenser, Desuperheater],
        indirect=True
    )
    def test_postprocessing(self, heatexchanger_converged_network):
        # rename for shorter variable name
        nw = heatexchanger_converged_network

        c1, c2, c3, c4 = nw.get_conn(['1', '2', '3', '4'])
        instance = nw.get_comp('heat exchanger')

        assert approx(instance.Q.val) == _calc_Q(c1, c2)

        assert approx(instance.pr1.val) == _calc_pr(c1, c2)
        assert approx(instance.pr2.val) == _calc_pr(c3, c4)

        # this already utilizes SI values, it should be changed!
        assert approx(instance.dp1.val_SI) == _calc_dp(c1, c2)
        assert approx(instance.dp2.val_SI) == _calc_dp(c3, c4)

        assert approx(instance.zeta1.val) == _calc_zeta(c1, c2)
        assert approx(instance.zeta2.val) == _calc_zeta(c3, c4)

        assert approx(instance.ttd_l.val) == _calc_ttd_l(c2, c3)
        # this is not really beautiful
        if type(instance) == Condenser:
            assert approx(instance.ttd_u.val) == _calc_ttd_u_condenser(c1, c4)
        else:
            assert approx(instance.ttd_u.val) == _calc_ttd_u(c1, c4)

        assert instance.ttd_min.val == min(instance.ttd_l.val, instance.ttd_u.val)

        assert approx(instance.td_log.val) == _calc_td_log(instance.ttd_u.val, instance.ttd_l.val)
        assert approx(instance.kA.val) == _calc_kA(instance.Q.val, instance.td_log.val)

        assert approx(instance.eff_hot.val) == _calc_eff_hot(c1, c2, c3)
        assert approx(instance.eff_cold.val) == _calc_eff_cold(c3, c4, c1)

        assert instance.eff_max.val == max(instance.eff_cold.val, instance.eff_hot.val)

    @mark.parametrize(
        'heatexchanger_network',
        [HeatExchanger, Condenser, Desuperheater],
        indirect=True
    )
    def test_pr1(self, heatexchanger_converged_network):
        nw = heatexchanger_converged_network

        c1, c2 = nw.get_conn(['1', '2'])
        instance = nw.get_comp('heat exchanger')

        pr = 0.95
        c2.set_attr(p=None)
        instance.set_attr(pr1=pr)

        nw.solve('design')
        nw.assert_convergence()

        assert approx(pr) == _calc_pr(c1, c2)

    @mark.parametrize(
        'heatexchanger_network',
        [HeatExchanger, Condenser, Desuperheater],
        indirect=True
    )
    def test_zeta1(self, heatexchanger_converged_network):
        nw = heatexchanger_converged_network

        c1, c2 = nw.get_conn(['1', '2'])
        instance = nw.get_comp('heat exchanger')

        zeta = 100
        c2.set_attr(p=None)
        instance.set_attr(zeta1=zeta)

        nw.solve('design')
        nw.assert_convergence()

        assert approx(zeta) == _calc_zeta(c1, c2)

    @mark.parametrize(
        'heatexchanger_network',
        [HeatExchanger, Condenser, Desuperheater],
        indirect=True
    )
    def test_pr2(self, heatexchanger_converged_network):
        nw = heatexchanger_converged_network

        c3, c4 = nw.get_conn(['3', '4'])
        instance = nw.get_comp('heat exchanger')

        pr = 0.95
        c4.set_attr(p=None)
        instance.set_attr(pr2=pr)

        nw.solve('design')
        nw.assert_convergence()

        assert approx(pr) == _calc_pr(c3, c4)

    @mark.parametrize(
        'heatexchanger_network',
        [HeatExchanger, Condenser, Desuperheater],
        indirect=True
    )
    def test_zeta2(self, heatexchanger_converged_network):
        nw = heatexchanger_converged_network

        c3, c4 = nw.get_conn(['3', '4'])
        instance = nw.get_comp('heat exchanger')

        zeta = 100
        c4.set_attr(p=None)
        instance.set_attr(zeta2=zeta)

        nw.solve('design')
        nw.assert_convergence()

        assert approx(zeta) == _calc_zeta(c3, c4)

class TestHeatExchangers:

    def setup_method(self):

        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "temperature": "degC",
            "volumetric_flow": "m3/s"
        })
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
        instance.set_attr(D='var')
        b = Bus('heat', P=-1e5)
        b.add_comps({'comp': instance})
        self.nw.add_busses(b)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert self.nw.status == 0
        pr = round(self.c2.p.val_SI / self.c1.p.val_SI, 3)
        msg = f"Value of pressure ratio must be {pr}, is {instance.pr.val}."
        assert pr == round(instance.pr.val, 3), msg

        # make zeta system variable and use previously calculated diameter
        # to calculate zeta. The value for zeta must not change
        zeta = round(instance.zeta.val, 0)
        diameter = instance.D.val
        instance.set_attr(D=None, zeta='var', pr=None)
        instance.set_attr(D=diameter)
        self.nw.solve('design')
        self.nw.assert_convergence()
        msg = f"Value of pressure ratio must be {zeta}, is {instance.zeta.val}."
        assert zeta == round(instance.zeta.val, 0), msg
        assert round(diameter, 3) == round(instance.D.val, 3)

        # test heat transfer coefficient as variable of the system (ambient
        # temperature required)
        instance.set_attr(kA='var', zeta=None)
        b.set_attr(P=-5e4)
        self.nw.solve('design')
        self.nw.assert_convergence()

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
        self.nw.assert_convergence()
        msg = f"Value of heat transfer must be {Q}, is {instance.Q.val}."
        assert Q == round(instance.Q.val, 0), msg

        # test kA as network results parameter
        instance.set_attr(Q=-5e4, Tamb=None)
        b.set_attr(P=None)
        self.nw.solve('design')
        self.nw.assert_convergence()
        kA_network = self.nw.results['SimpleHeatExchanger'].loc[
            instance.label, 'kA'
        ]
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

    def test_SimpleHeatExchanger_kA_convergence_cooling(self):
        """Test kA group equation convergence for near ambient temperature outflow."""
        instance = SimpleHeatExchanger("heatexchanger")
        self.setup_SimpleHeatExchanger_network(instance)
        self.c1.set_attr(fluid={"H2O": 1}, x=0, p=28, m=0.005)
        instance.D.min_val = 0.001
        instance.set_attr(
            pr=0.95,
            Tamb=20,
            kA=250,
            L=1000,
            D='var',
            ks=4.57e-5
        )
        self.nw.solve("design")
        self.nw.assert_convergence()
        assert round(self.c2.T.val - instance.Tamb.val, 3) == 0.002

    def test_SimpleHeatExchanger_kA_convergence_heating(self):
        instance = SimpleHeatExchanger("heatexchanger")
        self.setup_SimpleHeatExchanger_network(instance)

        instance.set_attr(Tamb=10, D=0.0215, L=50, ks=0.00001, kA=3000.0)
        self.c1.set_attr(p=4, m=0.1, fluid={"Water": 1}, T=0)
        self.nw.solve("design")
        self.nw.assert_convergence()
        assert round(instance.Tamb.val - self.c2.T.val, 3) == 0.008

    def test_ParabolicTrough(self):
        """Test component properties of parabolic trough."""
        instance = ParabolicTrough('parabolic trough')
        self.setup_SimpleHeatExchanger_network(instance)
        self.c1.set_attr(fluid={'INCOMP::S800': 1}, p=10, T=200)
        self.c2.set_attr(T=350)

        # test grouped parameter settings with missing parameters
        instance.darcy_group.is_set = True
        instance.energy_group.is_set = True
        self.nw.solve('design', init_only=True)
        msg = 'Darcy group must no be set, if one parameter is missing!'
        assert not instance.darcy_group.is_set, msg
        msg = 'Energy group must no be set, if one parameter is missing!'
        assert not instance.energy_group.is_set, msg

        # test solar collector params as system variables
        instance.set_attr(
            pr=1, aoi=10, doc=0.95, Q=1e6, Tamb=25, A='var', eta_opt=0.816,
            c_1=0.0622, c_2=0.00023, E=8e2, iam_1=-1.59e-3, iam_2=9.77e-5
        )
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert self.nw.status == 0
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
            f'Value for heat loss of parabolic trough must be {Q_loss}, is '
            f'{round(instance.Q_loss.val, 0)}.'
        )
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: E
        # going to a different operating point first
        area = instance.A.val
        instance.set_attr(A=area * 1.2, E='var')
        self.nw.solve('design')
        instance.set_attr(A=area)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: eta_opt
        instance.set_attr(E=5e2, eta_opt='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: c_1
        instance.set_attr(E=5e2, eta_opt=instance.eta_opt.val, c_1='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: c_2
        instance.set_attr(E=5e2, c_1=instance.c_1.val, c_2='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: iam_1
        instance.set_attr(E=5e2, c_2=instance.c_2.val, iam_1='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: iam_2
        instance.set_attr(E=5e2, iam_1=instance.iam_1.val, iam_2='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: aoi
        instance.set_attr(E=5e2, iam_2=instance.iam_2.val, aoi='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert Q_loss == round(instance.Q_loss.val, 0), msg

        # test all parameters of the energy group: doc
        instance.set_attr(E=5e2, aoi=instance.aoi.val, doc='var')
        self.nw.solve('design')
        instance.set_attr(E=8e2)
        self.nw.solve('design')
        self.nw.assert_convergence()
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
        instance.set_attr(
            E=1e3, lkf_lin=1.0, lkf_quad=0.005, A='var',
            eta_opt=0.9, Q=1e5, Tamb=20, pr=0.99
        )
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert self.nw.status == 0
        # heat loss must be identical to E * A - Q (internal heat loss
        # calculation)
        T_diff = (self.c2.T.val_SI + self.c1.T.val_SI) / 2 - instance.Tamb.val_SI
        Q_loss = -round(instance.A.val_SI * (
            instance.E.val_SI * (1 - instance.eta_opt.val_SI) +
            T_diff * instance.lkf_lin.val_SI +
            T_diff ** 2 * instance.lkf_quad.val_SI), 0)
        msg = (
            f"Value for heat loss of solar collector must be {Q_loss}, is "
            f"{round(instance.Q_loss.val_SI, 0)}."
        )
        assert Q_loss == round(instance.Q_loss.val_SI, 0), msg

        # test all parameters of the energy group: E
        area = instance.A.val_SI
        instance.set_attr(A=area * 1.2, E='var')
        self.nw.solve('design')
        instance.set_attr(A=area)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert Q_loss == round(instance.Q_loss.val_SI, 0), msg

        # test all parameters of the energy group: eta_opt
        instance.set_attr(E=8e2, eta_opt='var')
        self.nw.solve('design')
        instance.set_attr(E=1e3)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert Q_loss == round(instance.Q_loss.val_SI, 0), msg

        # test all parameters of the energy group: lkf_lin
        instance.set_attr(E=8e2, eta_opt=instance.eta_opt.val, lkf_lin='var')
        self.nw.solve('design')
        instance.set_attr(E=1e3)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert Q_loss == round(instance.Q_loss.val_SI, 0), msg

        # test all parameters of the energy group: lkf_quad
        instance.set_attr(E=8e2, lkf_lin=instance.lkf_lin.val, lkf_quad='var')
        self.nw.solve('design')
        instance.set_attr(E=1e3)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert Q_loss == round(instance.Q_loss.val_SI, 0), msg

        # test all parameters of the energy group: Tamb
        instance.set_attr(E=8e2, lkf_lin=instance.lkf_lin.val, lkf_quad='var')
        self.nw.solve('design')
        instance.set_attr(E=1e3)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert Q_loss == round(instance.Q_loss.val_SI, 0), msg

    def test_HeatExchanger(self, tmp_path):
        """Test component properties of heat exchanger."""
        tmp_path = f'{tmp_path}.json'
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
        self.nw.assert_convergence()
        assert self.nw.status == 0
        self.nw.save(tmp_path)
        Q_design = instance.Q.val

        # test specified kA value
        instance.set_attr(kA=instance.kA.val * 2 / 3)
        b.set_attr(P=None)
        self.nw.solve('design')
        self.nw.assert_convergence()

        # test heat transfer
        Q = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        msg = (
            f'Value of heat flow must be {round(Q_design * 2 / 3, 0)}, is '
            f'{round(Q, 0)}.'
        )
        assert round(Q, 1) == round(Q_design * 2 / 3, 1), msg

        # back to design case
        instance.set_attr(kA=None)
        b.set_attr(P=Q_design)
        self.nw.solve('design')
        self.nw.assert_convergence()

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
        msg = (
            'Value of terminal temperature difference must be '
            f'{round(instance.ttd_u.val, 1)}, is '
            f'{round(self.c1.T.val - self.c4.T.val, 1)}.'
        )
        ttd_u_calc = round(self.c1.T.val - self.c4.T.val, 1)
        ttd_u = round(instance.ttd_u.val, 1)
        assert ttd_u_calc == ttd_u, msg

        # check lower terminal temperature difference
        self.c2.set_attr(T=None)
        instance.set_attr(ttd_l=20)
        self.nw.solve('design')
        self.nw.assert_convergence()
        msg = (
            'Value of terminal temperature difference must be '
            f'{instance.ttd_l.val}, is {self.c2.T.val - self.c3.T.val}.'
        )
        ttd_l_calc = round(self.c2.T.val - self.c3.T.val, 1)
        ttd_l = round(instance.ttd_l.val, 1)
        assert ttd_l_calc == ttd_l, msg

        # check specified kA value (by offdesign parameter), reset temperatures
        # to design state
        self.c2.set_attr(T=70)
        instance.set_attr(ttd_l=None)
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw.assert_convergence()
        msg = f'Value of heat flow must be {instance.Q.val}, is {round(Q, 0)}.'
        assert round(Q, 0) == round(instance.Q.val, 0), msg
        msg = (
            f'Value of heat transfer coefficient must be {kA}, is '
            f'{round(instance.kA.val, 0)}.'
        )
        assert kA == round(instance.kA.val, 0), msg

        # trigger negative lower terminal temperature difference as result
        self.c4.set_attr(T=None)
        self.c2.set_attr(T=30)
        self.nw.solve('design')
        self.nw.assert_convergence()
        msg = (
            'Value of upper terminal temperature differences must be '
            f'smaller than zero, is {round(instance.ttd_l.val, 1)}.'
        )
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
        self.nw.assert_convergence()
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
        self.nw.assert_convergence()
        msg = (
            'Value of cold side effectiveness must be equal to 0.9 but is '
            f'{round(instance.eff_cold.val, 1)}.'
        )
        assert round(instance.eff_cold.val, 1) == 0.9, msg

        self.c1.set_attr(p=None)
        self.c2.set_attr(p=3)
        self.c3.set_attr(m=None, p=None)
        self.c4.set_attr(p=5)
        instance.set_attr(eff_max=None, eff_hot=0.9, eff_cold=0.9)
        self.nw.solve("design")
        self.nw.assert_convergence()
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
        self.nw.assert_convergence()
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

    def test_HeatExchanger_ttd_zero(self):

        instance = HeatExchanger('heat exchanger')
        self.setup_HeatExchanger_network(instance)

        # add new fluids
        # temperature range > 300 °C
        self.c1.set_attr(fluid={"H2O": 1}, m=10, T=100, p=1)
        self.c2.set_attr(T=50)
        # temperature range < 100 °C at 1 bar
        self.c3.set_attr(fluid={"air": 1}, T=50, p=1)
        self.c4.set_attr(T=75)
        instance.set_attr(dp1=0, dp2=0)

        self.nw.solve("design")
        self.nw.assert_convergence()
        msg = (
            'Value of heat transfer coefficient must be nan but is '
            f'{round(instance.kA.val, 1)}.'
        )
        assert np.isnan(instance.kA.val), msg

    def test_Condenser(self, tmp_path):
        """Test component properties of Condenser."""
        tmp_path = f'{tmp_path}.json'
        instance = Condenser('condenser')
        self.setup_HeatExchanger_network(instance)

        # design specification
        instance.set_attr(
            pr1=0.98, pr2=0.98, ttd_u=5, offdesign=['zeta2', 'kA_char']
        )
        self.c1.set_attr(T=100, p0=0.5, fluid={'H2O': 1})
        self.c2.set_attr(p0=0.5)
        self.c3.set_attr(T=30, p=5, fluid={'H2O': 1})
        self.c4.set_attr(T=40)
        instance.set_attr(Q=-80e3)
        self.nw.solve('design')
        self.nw.assert_convergence()
        assert self.nw.status == 0
        self.nw.save(tmp_path)
        Q_design = instance.Q.val

        # test specified kA value
        instance.set_attr(kA=instance.kA.val * 2 / 3, Q=None)
        self.nw.solve('design')
        self.nw.assert_convergence()

        # test heat transfer
        Q = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        msg = (
            f'Value of heat flow must be {round(Q_design * 2 / 3, 0)}, is '
            f'{round(Q, 0)}.'
        )
        assert round(Q, 1) == round(Q_design * 2 / 3, 1), msg

        # back to design case
        instance.set_attr(kA=None, Q=Q_design)
        self.nw.solve('design')
        self.nw.assert_convergence()

        # test heat transfer
        Q = self.c1.m.val_SI * (self.c2.h.val_SI - self.c1.h.val_SI)
        msg = (
            f'Value of heat flow must be {round(instance.Q.val, 0)}, is '
            f'{round(Q, 0)}.'
        )
        assert round(Q, 1) == round(instance.Q.val, 1), msg

        # test upper terminal temperature difference. For the component
        # condenser the temperature of the condensing fluid is relevant.
        ttd_u = round(self.c1.calc_T_sat() - self.c4.T.val_SI, 1)
        p = round(self.c1.p.val_SI, 5)
        msg = (
            'Value of terminal temperature difference must be '
            f'{round(instance.ttd_u.val, 1)}, is {ttd_u}.'
        )
        assert ttd_u == round(instance.ttd_u.val, 1), msg

        # test lower terminal temperature difference
        instance.set_attr(ttd_l=20, ttd_u=None, design=['pr2', 'ttd_l'])
        self.nw.solve('design')
        self.nw.assert_convergence()
        msg = (
            'Value of terminal temperature difference must be '
            f'{instance.ttd_l.val}, is {self.c2.T.val - self.c3.T.val}.'
        )
        ttd_l_calc = round(self.c2.T.val - self.c3.T.val, 1)
        ttd_l = round(instance.ttd_l.val, 1)
        assert ttd_l_calc == ttd_l, msg

        # check kA value with condensing pressure in offdesign mode:
        # no changes to design point means: identical pressure
        self.nw.solve('offdesign', design_path=tmp_path)
        self.nw.assert_convergence()
        msg = (
            f'Value of condensing pressure be {p}, is '
            f'{round(self.c1.p.val_SI, 5)}.'
        )
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
        self.nw.assert_convergence()
        ttd_l = round(instance.ttd_l.val, 3)
        ttd_u = round(instance.ttd_u.val, 3)
        td_log = round(instance.td_log.val, 3)
        msg = (
            "Value of upper and lower terminal temperature differences must be "
            f"identical, but they are not: ttd_u={ttd_u}, ttd_l={ttd_l}."
        )
        assert approx(instance.ttd_l.val) == instance.ttd_u.val, msg

        msg = (
            "Value of logarithmic and lower terminal temperature differences "
            f"must be identical, but they are not: td_log={td_log}, "
            f"ttd_l={ttd_l}."
        )
        assert instance.td_log.val == instance.ttd_l.val, msg

    def test_Condenser_ttd_zero(self):

        instance = Condenser('condenser')
        self.setup_HeatExchanger_network(instance)

        # add new fluids
        # temperature range > 300 °C
        self.c1.set_attr(fluid={"H2O": 1}, m=1, T=100)
        self.c2.set_attr(p=0.5)
        # temperature range < 100 °C at 1 bar
        self.c3.set_attr(fluid={"air": 1}, p=1)
        self.c4.set_attr(T=75)
        instance.set_attr(dp1=0, dp2=0, ttd_l=0)

        self.nw.solve("design")
        self.nw.print_results()
        self.nw.assert_convergence()
        msg = (
            'Value of heat transfer coefficient must be nan but is '
            f'{round(instance.kA.val, 1)}.'
        )
        assert np.isnan(instance.kA.val), msg


    def test_ParallelFlowHeatExchanger(self):
        instance = ParallelFlowHeatExchanger("heat exchanger")
        self.setup_HeatExchanger_network(instance)

        self.c1.set_attr(fluid={"air": 1}, m=1, T=85, p=1)
        self.c3.set_attr(fluid={"water": 1}, m=3, T=25, p=1)
        instance.set_attr(dp1=0.1, dp2=0.1, ttd_u=15)

        self.nw.solve("design")
        self.nw.assert_convergence()

        ttd_l = self.c1.T.val - self.c3.T.val
        ttd_u = self.c2.T.val - self.c4.T.val
        assert approx(instance.ttd_l.val) == ttd_l
        assert approx(instance.ttd_u.val) == ttd_u

        assert approx(instance.dp1.val_SI) == _calc_dp(self.c1, self.c2)
        assert approx(instance.dp2.val_SI) == _calc_dp(self.c3, self.c4)

        assert approx(instance.pr1.val_SI) == _calc_pr(self.c1, self.c2)
        assert approx(instance.pr2.val_SI) == _calc_pr(self.c3, self.c4)

        assert approx(instance.td_log.val_SI) == _calc_td_log(ttd_u, ttd_l)
        assert approx(instance.kA.val) == _calc_kA(
            instance.Q.val, instance.td_log.val
        )
