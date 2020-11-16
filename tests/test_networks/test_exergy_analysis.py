# -*- coding: utf-8

"""Module for testing network properties.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_network.py

SPDX-License-Identifier: MIT
"""

import os
import shutil

import numpy as np
from pytest import raises

from tespy.components.basics import cycle_closer
from tespy.components.heat_exchangers import heat_exchanger_simple
from tespy.components.nodes import merge
from tespy.components.nodes import splitter
from tespy.components.turbomachinery import pump
from tespy.components.turbomachinery import turbine
from tespy.connections import bus
from tespy.connections import connection
from tespy.networks.networks import network
from tespy.tools.global_vars import err


def convergence_check(lin_dep):
    """Check convergence status of a simulation."""
    msg = 'Calculation did not converge!'
    assert lin_dep is False, msg


class TestClausiusRankine:

    def setup(self):
        """Setup clausis rankine cycle with turbine driven feed water pump."""
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
        condenser.set_attr(pr=1, Tamb=self.Tamb)
        fwp.set_attr(eta_s=1)
        steam_generator.set_attr(pr=1, Tamb=self.Tamb)

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
        # miss out on internal bus in exergy_analysis
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
