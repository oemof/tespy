# -*- coding: utf-8

"""Module for testing busses.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_busses.py

SPDX-License-Identifier: MIT
"""
import shutil

import numpy as np

from tespy.components import CombustionChamber
from tespy.components import Compressor
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Turbine
from tespy.connections import Bus
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools.characteristics import CharLine


class TestBusses:

    def setup(self):
        """Set up the model."""
        # %% network setup
        fluid_list = ['Ar', 'N2', 'O2', 'CO2', 'CH4', 'H2O']
        self.nw = Network(
            fluids=fluid_list, p_unit='bar', T_unit='C', p_range=[0.5, 20])

        # %% components
        amb = Source('ambient')
        sf = Source('fuel')
        cc = CombustionChamber('combustion')
        cp = Compressor('compressor')
        gt = Turbine('turbine')
        fg = Sink('flue gas outlet')

        # %% connections
        amb_cp = Connection(amb, 'out1', cp, 'in1', label='ambient air flow')
        cp_cc = Connection(cp, 'out1', cc, 'in1')
        sf_cc = Connection(sf, 'out1', cc, 'in2')
        cc_gt = Connection(cc, 'out1', gt, 'in1')
        gt_fg = Connection(gt, 'out1', fg, 'in1')

        self.nw.add_conns(amb_cp, cp_cc, sf_cc, cc_gt, gt_fg)

        # %% component parameters
        cc.set_attr(lamb=3)
        cp.set_attr(eta_s=0.9, pr=15)
        gt.set_attr(eta_s=0.9)

        # %% connection parameters
        amb_cp.set_attr(
            T=20, p=1, m=100,
            fluid={'Ar': 0.0129, 'N2': 0.7553, 'H2O': 0, 'CH4': 0,
                   'CO2': 0.0004, 'O2': 0.2314})
        sf_cc.set_attr(
            T=20, fluid={'CO2': 0.04, 'Ar': 0, 'N2': 0,
                         'O2': 0, 'H2O': 0, 'CH4': 0.96})
        gt_fg.set_attr(p=1)

        # motor efficiency
        x = np.array(
            [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
             0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15,
             1.2, 10])
        y = np.array(
            [0.01, 0.3148, 0.5346, 0.6843, 0.7835, 0.8477, 0.8885, 0.9145,
             0.9318, 0.9443, 0.9546, 0.9638, 0.9724, 0.9806, 0.9878, 0.9938,
             0.9982, 0.999, 0.9995, 0.9999, 1, 0.9977, 0.9947, 0.9909, 0.9853,
             0.9644]) * 0.975
        self.motor_bus_based = CharLine(x=x, y=y)
        self.motor_comp_based = CharLine(x=x, y=1 / y)

        # generator efficiency
        x = np.array(
            [0.100, 0.345, 0.359, 0.383, 0.410, 0.432, 0.451, 0.504, 0.541,
             0.600, 0.684, 0.805, 1.000, 1.700, 10])
        y = np.array(
            [0.976, 0.989, 0.990, 0.991, 0.992, 0.993, 0.994, 0.995, 0.996,
             0.997, 0.998, 0.999, 1.000, 0.999, 0.99]) * 0.975
        self.generator = CharLine(x=x, y=y)

        power_bus_total = Bus('total power output')
        power_bus_total.add_comps(
            {'comp': cp, 'char': self.motor_bus_based, 'base': 'bus'},
            {'comp': gt, 'char': self.generator})

        thermal_input = Bus('thermal input')
        thermal_input.add_comps({'comp': cc})

        compressor_power_comp = Bus('compressor power input')
        compressor_power_comp.add_comps(
            {'comp': cp, 'char': self.motor_comp_based})

        compressor_power_bus = Bus('compressor power input bus based')
        compressor_power_bus.add_comps(
            {'comp': cp, 'char': self.motor_bus_based, 'base': 'bus'})

        self.nw.add_busses(
            power_bus_total, thermal_input, compressor_power_comp,
            compressor_power_bus)

        # %% solving
        self.nw.solve('design')
        self.nw.save('tmp')

    def test_model(self):
        """Test the bus functionalities in a gas turbine model."""
        tpo = self.nw.busses['total power output']
        ti = self.nw.busses['thermal input']
        cpi = self.nw.busses['compressor power input']
        cpibb = self.nw.busses['compressor power input bus based']

        cp = self.nw.get_comp('compressor')
        gt = self.nw.get_comp('turbine')
        cc = self.nw.get_comp('combustion')

        # test results of design case

        eta_cpi = round(1 / cp.calc_bus_efficiency(cpi), 6)
        eta_cp_tpo = round(cp.calc_bus_efficiency(tpo), 6)
        msg = (
            'The efficiency value of the compressor on the bus ' + tpo.label +
            ' (' + str(eta_cp_tpo) + ') must be identical to the efficiency '
            'on the bus ' + cpi.label + ' (' + str(eta_cpi) + ').')
        assert eta_cp_tpo == eta_cpi, msg

        P_cp_tpo = cp.calc_bus_value(tpo)
        eta_cp_tpo = cp.calc_bus_efficiency(tpo)
        P_cp = round(P_cp_tpo * eta_cp_tpo, 0)
        msg = (
            'The compressor power must be ' + str(round(cp.P.val, 0)) + ' on '
            'the bus ' + tpo.label + ' but is ' + str(P_cp) + ').')
        assert round(cp.P.val, 0) == P_cp, msg

        P_cp_tpo = round(
            cp.calc_bus_value(tpo) * cp.calc_bus_efficiency(tpo), 0)
        P_cp_cpi = round(
            cp.calc_bus_value(cpi) / cp.calc_bus_efficiency(cpi), 0)
        P_cp_cpibb = round(
            cp.calc_bus_value(cpibb) * cp.calc_bus_efficiency(cpibb), 0)
        msg = (
            'The busses\' component power value for the compressor on bus ' +
            tpo.label + ' (' + str(P_cp_tpo) + ') must be equal to the '
            'component power on all other busses. Bus ' + cpi.label + ' (' +
            str(P_cp_cpi) + ') and bus ' + cpibb.label + ' (' +
            str(P_cp_cpibb) + ').')
        assert P_cp_tpo == P_cp_cpi and P_cp_tpo == P_cp_cpibb, msg

        eta_gt_tpo = gt.calc_bus_efficiency(tpo)
        msg = (
            'The efficiency value of the turbine on the bus ' + tpo.label +
            ' (' + str(eta_gt_tpo) + ') must be equal to 0.975.')
        assert eta_gt_tpo == 0.975, msg

        eta_ti = cc.calc_bus_efficiency(ti)
        msg = (
            'The efficiency value of the combustion chamber on the bus ' +
            ti.label + ' (' + str(eta_ti) + ') must be equal to 1.0.')
        assert eta_ti == 1.0, msg

        # test partload for bus functions
        # first test in identical conditions

        self.nw.get_conn('ambient air flow').set_attr(m=None)
        P_design = cpibb.P.val
        cpibb.set_attr(P=P_design)
        self.nw.solve('offdesign', design_path='tmp')

        eta_cpi = round(1 / cp.calc_bus_efficiency(cpi), 6)
        eta_cp_tpo = round(cp.calc_bus_efficiency(tpo), 6)
        msg = (
            'The efficiency value of the compressor on the bus ' + tpo.label +
            ' (' + str(eta_cp_tpo) + ') must be identical to the efficiency '
            'on the bus ' + cpi.label + ' (' + str(eta_cpi) + ').')
        assert eta_cp_tpo == eta_cpi, msg

        eta_gt_tpo = gt.calc_bus_efficiency(tpo)
        msg = (
            'The efficiency value of the turbine on the bus ' + tpo.label +
            ' (' + str(eta_gt_tpo) + ') must be equal to 0.975.')
        assert eta_gt_tpo == 0.975, msg

        P_cp_tpo = round(
            cp.calc_bus_value(tpo) * cp.calc_bus_efficiency(tpo), 0)
        P_cp_cpi = round(
            cp.calc_bus_value(cpi) / cp.calc_bus_efficiency(cpi), 0)
        P_cp_cpibb = round(
            cp.calc_bus_value(cpibb) * cp.calc_bus_efficiency(cpibb), 0)
        msg = (
            'The busses\' component power value for the compressor on bus ' +
            tpo.label + ' (' + str(P_cp_tpo) + ') must be equal to the '
            'component power on all other busses. Bus ' + cpi.label + ' (' +
            str(P_cp_cpi) + ') and bus ' + cpibb.label + ' (' +
            str(P_cp_cpibb) + ').')
        assert P_cp_tpo == P_cp_cpi and P_cp_tpo == P_cp_cpibb, msg

        # 60 % load
        load = 0.6
        cpibb.set_attr(P=P_design * load)
        self.nw.solve('offdesign', design_path='tmp')

        eta_cp_tpo = round(cp.calc_bus_efficiency(tpo), 6)
        eta_cp_char = self.motor_bus_based.evaluate(load)
        msg = (
            'The efficiency value of the compressor on the bus ' + tpo.label +
            ' (' + str(eta_cp_tpo) + ') must be identical to the efficiency '
            'on the characteristic line (' + str(eta_cp_char) + ').')
        assert eta_cp_tpo == eta_cp_char, msg

        load_frac = round(
            cp.calc_bus_value(tpo) / tpo.comps.loc[cp, 'P_ref'], 6)
        msg = (
            'The load fraction value of the compressor on the bus ' +
            tpo.label + ' (' + str(load_frac) + ') must be identical to the '
            'load fraction value on the bus ' + cpibb.label + ' (' +
            str(load) + ').')
        assert load == load_frac, msg

        eta_cpi = round(1 / cp.calc_bus_efficiency(cpi), 6)
        eta_cp_tpo = round(cp.calc_bus_efficiency(tpo), 6)
        msg = (
            'The efficiency value of the compressor on the bus ' + tpo.label +
            ' (' + str(eta_cp_tpo) + ') must be higher than the efficiency '
            'on the bus ' + cpi.label + ' (' + str(eta_cpi) + ').')
        assert eta_cp_tpo > eta_cpi, msg

        P_cp_tpo = round(
            cp.calc_bus_value(tpo) * cp.calc_bus_efficiency(tpo), 0)
        P_cp_cpi = round(
            cp.calc_bus_value(cpi) / cp.calc_bus_efficiency(cpi), 0)
        P_cp_cpibb = round(
            cp.calc_bus_value(cpibb) * cp.calc_bus_efficiency(cpibb), 0)
        msg = (
            'The busses\' component power value for the compressor on bus ' +
            tpo.label + ' (' + str(P_cp_tpo) + ') must be equal to the '
            'component power on all other busses. Bus ' + cpi.label + ' (' +
            str(P_cp_cpi) + ') and bus ' + cpibb.label + ' (' +
            str(P_cp_cpibb) + ').')
        assert P_cp_tpo == P_cp_cpi and P_cp_tpo == P_cp_cpibb, msg

        shutil.rmtree('tmp', ignore_errors=True)
