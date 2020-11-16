# -*- coding: utf-8

"""Module for testing network properties.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_exergy_analysis.py

SPDX-License-Identifier: MIT
"""

import os
import shutil

import numpy as np
from pytest import raises

from CoolProp.CoolProp import PropsSI as PSI

from tespy.components.basics import cycle_closer
from tespy.components.basics import sink
from tespy.components.basics import source
from tespy.components.heat_exchangers import heat_exchanger_simple
from tespy.components.nodes import merge
from tespy.components.nodes import splitter
from tespy.components.piping import valve
from tespy.components.turbomachinery import compressor
from tespy.components.turbomachinery import pump
from tespy.components.turbomachinery import turbine
from tespy.connections import bus
from tespy.connections import connection
from tespy.networks.networks import network
from tespy.tools.fluid_properties import calc_physical_exergy
from tespy.tools.global_vars import err
from tespy.tools.helpers import TESPyComponentError
from tespy.tools.helpers import TESPyNetworkError


def convergence_check(lin_dep):
    """Check convergence status of a simulation."""
    msg = 'Calculation did not converge!'
    assert lin_dep is False, msg


class TestClausiusRankine:

    def setup(self):
        """Set up clausis rankine cycle with turbine driven feed water pump."""
        self.Tamb = 20
        self.pamb = 1
        fluids = ['water']
        self.nw = network(fluids=fluids)
        self.nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')

        # create components
        splitter1 = splitter('splitter 1')
        merge1 = merge('merge 1')
        turb = turbine('turbine')
        fwp_turb = turbine('feed water pump turbine')
        condenser = heat_exchanger_simple('condenser')
        fwp = pump('pump')
        steam_generator = heat_exchanger_simple('steam generator')
        cycle_close = cycle_closer('cycle closer')

        # create busses
        # power output bus
        self.power = bus('power_output')
        self.power.add_comps({'comp': turb, 'char': 1})
        # turbine driven feed water pump internal bus
        self.fwp_power = bus('feed water pump power', P=0)
        self.fwp_power.add_comps(
            {'comp': fwp_turb, 'char': 1},
            {'comp': fwp, 'char': 1, 'base': 'bus'})
        # heat input bus
        self.heat = bus('heat_input')
        self.heat.add_comps({'comp': steam_generator, 'base': 'bus'})
        self.nw.add_busses(self.power, self.fwp_power, self.heat)

        # create connections
        fs_in = connection(cycle_close, 'out1', splitter1, 'in1', label='fs')
        fs_fwpt = connection(splitter1, 'out1', fwp_turb, 'in1')
        fs_t = connection(splitter1, 'out2', turb, 'in1')
        fwpt_ws = connection(fwp_turb, 'out1', merge1, 'in1')
        t_ws = connection(turb, 'out1', merge1, 'in2')
        ws = connection(merge1, 'out1', condenser, 'in1')
        cond = connection(condenser, 'out1', fwp, 'in1', label='cond')
        fw = connection(fwp, 'out1', steam_generator, 'in1', label='fw')
        fs_out = connection(steam_generator, 'out1', cycle_close, 'in1')
        self.nw.add_conns(fs_in, fs_fwpt, fs_t, fwpt_ws, t_ws, ws, cond, fw,
                          fs_out)

        # component parameters
        turb.set_attr(eta_s=1)
        fwp_turb.set_attr(eta_s=1)
        condenser.set_attr(pr=1)
        fwp.set_attr(eta_s=1)
        steam_generator.set_attr(pr=1)

        # connection parameters
        fs_in.set_attr(m=10, p=120, T=600, fluid={'water': 1})
        cond.set_attr(T=self.Tamb, x=0)

        # solve network
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)

    def test_exergy_analysis_perfect_cycle(self):
        """Test exergy analysis in the perfect clausius rankine cycle."""
        self.nw.exergy_analysis(
            self.pamb, self.Tamb,
            E_P=[self.power], E_F=[self.heat], internal_busses=[self.fwp_power]
        )
        msg = (
            'Exergy destruction of this network must be 0 (smaller than ' +
            str(err ** 0.5) + ') for this test but is ' +
            str(round(abs(self.nw.E_D), 4)) + ' .')
        assert abs(self.nw.E_D) <= err ** 0.5, msg

        msg = (
            'Exergy efficiency of this network must be 1 for this test but '
            'is ' + str(round(self.nw.epsilon, 4)) + ' .')
        assert round(self.nw.epsilon, 4) == 1, msg

        exergy_balance = self.nw.E_F - self.nw.E_P - self.nw.E_L - self.nw.E_D
        msg = (
            'Exergy balance must be closed (residual value smaller than ' +
            str(err ** 0.5) + ') for this test but is ' +
            str(round(abs(exergy_balance), 4)) + ' .')
        assert abs(exergy_balance) <= err ** 0.5, msg

        msg = (
            'Fuel exergy and product exergy must be identical for this test. '
            'Fuel exergy value: ' + str(round(self.nw.E_F, 4)) +
            '. Product exergy value: ' + str(round(self.nw.E_P, 4)) + '.')
        assert round(abs(self.nw.E_F - self.nw.E_P), 4) < err ** 0.5, msg

    def test_exergy_analysis_violated_balance(self):
        """Test exergy analysis with violated balance."""
        # specify efficiency values for the internal bus
        self.nw.del_busses(self.fwp_power)
        self.fwp_power = bus('feed water pump power', P=0)
        self.fwp_power.add_comps(
            {'comp': self.nw.components['feed water pump turbine'],
             'char': 0.99},
            {'comp': self.nw.components['pump'], 'char': 0.98, 'base': 'bus'})
        self.nw.add_busses(self.fwp_power)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        # miss out on internal bus in exergy_analysis
        self.nw.exergy_analysis(
            self.pamb, self.Tamb,
            E_P=[self.power], E_F=[self.heat]
        )

        exergy_balance = self.nw.E_F - self.nw.E_P - self.nw.E_L - self.nw.E_D
        msg = (
            'Exergy balance must be violated for this test (larger than ' +
            str(err ** 0.5) + ') but is ' +
            str(round(abs(exergy_balance), 4)) + ' .')
        assert abs(exergy_balance) > err ** 0.5, msg

    def test_exergy_analysis_bus_conversion(self):
        """Test exergy analysis bus conversion factors."""
        # specify efficiency values for the internal bus
        self.nw.del_busses(self.fwp_power)
        self.fwp_power = bus('feed water pump power', P=0)
        self.fwp_power.add_comps(
            {'comp': self.nw.components['feed water pump turbine'],
             'char': 0.99},
            {'comp': self.nw.components['pump'], 'char': 0.98, 'base': 'bus'})
        self.nw.add_busses(self.fwp_power)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)
        self.nw.exergy_analysis(
            self.pamb, self.Tamb,
            E_P=[self.power], E_F=[self.heat], internal_busses=[self.fwp_power]
        )
        self.nw.print_exergy_analysis()
        eps = self.nw.component_exergy_data.loc['pump', 'epsilon']
        msg = (
            'Pump exergy efficiency must be 0.98 but is ' +
            str(round(eps, 4)) + ' .')
        assert round(eps, 4) == 0.98, msg

        eps = self.nw.component_exergy_data.loc[
            'feed water pump turbine', 'epsilon']
        msg = (
            'Feed water pump turbine exergy efficiency must be 0.99 but is ' +
            str(round(eps, 4)) + ' .')
        assert round(eps, 4) == 0.99, msg

    def test_exergy_analysis_missing_E_F_E_P_information(self):
        """Test exergy analysis errors with missing information."""
        with raises(TESPyNetworkError):
            self.nw.exergy_analysis(
                self.pamb, self.Tamb, E_P=[self.power], E_F=[])

        with raises(TESPyNetworkError):
            self.nw.exergy_analysis(
                self.pamb, self.Tamb, E_P=[], E_F=[self.heat])

    def test_exergy_analysis_component_on_two_busses(self):
        """Test exergy analysis errors with components on more than one bus."""
        with raises(TESPyNetworkError):
            self.nw.exergy_analysis(
                self.pamb, self.Tamb, E_P=[self.power],
                E_F=[self.heat, self.power])


class TestRefrigerator:

    def setup(self):
        """Set up simple refrigerator."""
        self.Tamb = 20
        self.pamb = 1
        fluids = ['R134a']
        self.nw = network(fluids=fluids)
        self.nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')

        # create components
        va = valve('expansion valve')
        cp = compressor('compressor')
        cond = heat_exchanger_simple('condenser')
        eva = heat_exchanger_simple('evaporator')
        cc = cycle_closer('cycle closer')

        # create busses
        # power output bus
        self.power = bus('power input')
        self.power.add_comps({'comp': cp, 'char': 0.97, 'base': 'bus'})
        # cooling bus
        self.cool = bus('heat from fridge')
        self.cool.add_comps({'comp': eva})
        # heat input bus
        self.heat = bus('heat to ambient')
        self.heat.add_comps({'comp': cond})
        self.nw.add_busses(self.power, self.cool, self.heat)

        # create connections
        cc_cp = connection(cc, 'out1', cp, 'in1', label='from eva')
        cp_cond = connection(cp, 'out1', cond, 'in1', label='to cond')
        cond_va = connection(cond, 'out1', va, 'in1', label='from cond')
        va_eva = connection(va, 'out1', eva, 'in1', label='to eva')
        eva_cc = connection(eva, 'out1', cc, 'in1')
        self.nw.add_conns(cc_cp, cp_cond, cond_va, va_eva, eva_cc)

        # component parameters
        cp.set_attr(eta_s=0.9)
        cond.set_attr(pr=0.97)
        eva.set_attr(pr=0.96)

        # connection parameters
        cc_cp.set_attr(m=1, x=1, T=-25, fluid={'R134a': 1})
        cond_va.set_attr(x=0, T=self.Tamb)

        # solve network
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)

    def test_exergy_analysis_bus_conversion(self):
        """Test exergy analysis at product exergy with T < Tamb."""
        self.nw.exergy_analysis(
            self.pamb, self.Tamb, E_P=[self.cool], E_F=[self.power])

        exergy_balance = self.nw.E_F - self.nw.E_P - self.nw.E_L - self.nw.E_D
        msg = (
            'Exergy balance must be closed (residual value smaller than ' +
            str(err ** 0.5) + ') for this test but is ' +
            str(round(abs(exergy_balance), 4)) + ' .')
        assert abs(exergy_balance) <= err ** 0.5, msg

    def test_exergy_analysis_heat_transfer_to_higher_temperature(self):
        """Test error with heat transfer to higher temperature."""
        # we do not need to recalculate when changing only ambient temperature
        with raises(TESPyComponentError):
            self.nw.exergy_analysis(
                self.pamb, self.Tamb + 100, E_P=[self.cool], E_F=[self.power])


class TestCompressedAirIn:

    def setup(self):
        """Set up air compressor."""
        self.Tamb = 20
        self.pamb = 1
        fluids = ['Air']

        # compressor part
        self.nw = network(fluids=fluids)
        self.nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')

        # components
        amb = source('air intake', exergy='fuel')
        cp = compressor('compressor')
        cooler = heat_exchanger_simple('cooling')
        cas = sink('compressed air storage', exergy='product')

        # power input bus
        self.power_in = bus('power input')
        self.power_in.add_comps({'comp': cp, 'char': 1, 'base': 'bus'})
        # compressed air bus (not sure about this!)
        self.cas_in = bus('massflow into storage')
        self.cas_in.add_comps({'comp': cas}, {'comp': amb, 'base': 'bus'})
        self.nw.add_busses(self.power_in, self.cas_in)

        # create connections
        amb_cp = connection(amb, 'out1', cp, 'in1')
        cp_cool = connection(cp, 'out1', cooler, 'in1')
        cool_cas = connection(cooler, 'out1', cas, 'in1')
        self.nw.add_conns(amb_cp, cp_cool, cool_cas)

        # component parameters
        cp.set_attr(eta_s=1)
        cooler.set_attr(pr=1)

        # connection parameters
        amb_cp.set_attr(m=2, T=self.Tamb, p=self.pamb, fluid={'Air': 1})
        cool_cas.set_attr(T=self.Tamb, p=10)

        # solve network
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)

    def test_exergy_analysis_bus_conversion(self):
        """Test exergy analysis at product exergy with T < Tamb."""
        self.nw.exergy_analysis(
            self.pamb, self.Tamb, E_P=[self.cas_in], E_F=[self.power_in])

        exergy_balance = (
            self.nw.E_F - self.nw.E_P - self.nw.E_L - self.nw.E_D)
        msg = (
            'Exergy balance must be closed (residual value smaller than ' +
            str(err ** 0.5) + ') for this test but is ' +
            str(round(abs(exergy_balance), 4)) + ' .')
        assert abs(exergy_balance) <= err ** 0.5, msg


class TestCompressedAirOut:

    def setup(self):
        """Set up air compressed air turbine."""
        self.Tamb = 20
        self.pamb = 1
        fluids = ['Air']

        # turbine part
        self.nw = network(fluids=fluids)
        self.nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')

        # components
        cas = source('compressed air storage', exergy='fuel')
        reheater = heat_exchanger_simple('reheating')
        turb = turbine('turbine')
        amb = sink('air outlet', exergy='loss')

        # power ouput bus
        self.power_out = bus('power output')
        self.power_out.add_comps({'comp': turb, 'char': 1})
        # compressed air bus
        self.cas_out = bus('exergy in')
        self.cas_out.add_comps(
            {'comp': cas, 'base': 'bus'},
            {'comp': reheater, 'base': 'bus'})
        self.nw.add_busses(self.power_out, self.cas_out)

        # create connections
        cas_reheater = connection(cas, 'out1', reheater, 'in1')
        reheater_turb = connection(reheater, 'out1', turb, 'in1')
        turb_amb = connection(turb, 'out1', amb, 'in1', label='outlet')
        self.nw.add_conns(cas_reheater, reheater_turb, turb_amb)

        # component parameters
        turb.set_attr(eta_s=1)
        reheater.set_attr(pr=1)

        # connection parameters
        cas_reheater.set_attr(m=2, T=self.Tamb, p=10, fluid={'Air': 1})
        reheater_turb.set_attr()
        turb_amb.set_attr(p=self.pamb, T=self.Tamb)

        # solve network
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)

    def test_exergy_analysis_bus_conversion(self):
        """Test exergy analysis at product exergy with T < Tamb."""
        self.nw.exergy_analysis(
            self.pamb, self.Tamb, E_P=[self.power_out], E_F=[self.cas_out])

        exergy_balance = (
            self.nw.E_F - self.nw.E_P - self.nw.E_L - self.nw.E_D)
        msg = (
            'Exergy balance must be closed (residual value smaller than ' +
            str(err ** 0.5) + ') for this test but is ' +
            str(round(abs(exergy_balance), 4)) + ' .')
        assert abs(exergy_balance) <= err ** 0.5, msg

        msg = (
            'Exergy efficiency must be equal to 1.0 for this test but is ' +
            str(round(self.nw.epsilon, 4)) + ' .')
        assert round(self.nw.epsilon, 4) == 1, msg

        c = self.nw.connections['outlet']
        c.set_attr(T=self.Tamb - 20)
        self.nw.solve('design')
        convergence_check(self.nw.lin_dep)

        self.nw.exergy_analysis(
            self.pamb, self.Tamb, E_P=[self.power_out], E_F=[self.cas_out])

        msg = (
            'Exergy destruction must be equal to 0.0 for this test but is ' +
            str(round(self.nw.E_D, 4)) + ' .')
        assert round(self.nw.E_D, 4) == 0, msg

        msg = (
            'Exergy loss must be equal to ' + str(round(c.Ex_physical, 4)) +
            ' for this test but is ' + str(round(self.nw.E_L, 4)) + ' .')
        assert round(self.nw.E_L, 4) == round(c.Ex_physical, 4), msg
