# -*- coding: utf-8 -*-

"""Module for testing two tespy simulation against each other.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_models/test_gasturbine_model.py

SPDX-License-Identifier: MIT
"""
import shutil

from tespy.components import CombustionChamber
from tespy.components import CombustionChamberStoich
from tespy.components import Compressor
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Turbine
from tespy.connections import Connection
from tespy.networks import Network


class TestGasturbine:

    def setup_CombustionChamber_model(self):
        """Set up the model using the combustion chamber."""
        # %% network setup
        fluid_list = ['Ar', 'N2', 'O2', 'CO2', 'CH4', 'H2O']
        self.nw1 = Network(
            fluids=fluid_list, p_unit='bar', T_unit='C', p_range=[0.5, 20])

        # %% components
        amb = Source('ambient')
        sf = Source('fuel')
        cc = CombustionChamber('combustion')
        cp = Compressor('compressor')
        gt = Turbine('turbine')
        fg = Sink('flue gas outlet')

        # %% connections
        amb_cp = Connection(amb, 'out1', cp, 'in1')
        cp_cc = Connection(cp, 'out1', cc, 'in1')
        sf_cc = Connection(sf, 'out1', cc, 'in2')
        cc_gt = Connection(cc, 'out1', gt, 'in1', label='flue gas after cc')
        gt_fg = Connection(gt, 'out1', fg, 'in1', label='flue gas after gt')

        self.nw1.add_conns(amb_cp, cp_cc, sf_cc, cc_gt, gt_fg)

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

        # %% solving
        mode = 'design'
        self.nw1.solve(mode=mode)

    def setup_CombustionChamberStoich_model(self):
        """Set up the model using the stoichimetric combustion chamber."""
        # %% network setup
        fluid_list = ['myAir', 'myFuel', 'myFuel_fg']
        self.nw2 = Network(
            fluids=fluid_list, p_unit='bar', T_unit='C',
            p_range=[0.5, 20], T_range=[10, 2000])

        # %% components
        amb = Source('ambient')
        sf = Source('fuel')
        cc = CombustionChamberStoich('combustion')
        cp = Compressor('compressor')
        gt = Turbine('turbine')
        fg = Sink('flue gas outlet')

        # %% connections
        amb_cp = Connection(amb, 'out1', cp, 'in1')
        cp_cc = Connection(cp, 'out1', cc, 'in1')
        sf_cc = Connection(sf, 'out1', cc, 'in2')
        cc_gt = Connection(cc, 'out1', gt, 'in1', label='flue gas after cc')
        gt_fg = Connection(gt, 'out1', fg, 'in1', label='flue gas after gt')

        self.nw2.add_conns(amb_cp, cp_cc, sf_cc, cc_gt, gt_fg)

        # %% component parameters
        cc.set_attr(
            fuel={'CH4': 0.96, 'CO2': 0.04},
            air={'Ar': 0.0129, 'N2': 0.7553, 'CO2': 0.0004, 'O2': 0.2314},
            fuel_alias='myFuel', air_alias='myAir', lamb=3)
        cp.set_attr(eta_s=0.9, pr=15)
        gt.set_attr(eta_s=0.9)

        # %% connection parameters
        amb_cp.set_attr(
            T=20, p=1, m=100, fluid={'myAir': 1, 'myFuel': 0, 'myFuel_fg': 0})
        sf_cc.set_attr(
            T=20, fluid={'myAir': 0, 'myFuel': 1, 'myFuel_fg': 0})
        gt_fg.set_attr(p=1)

        # %% solving
        self.nw2.solve(mode='design')

    def test_models(self):
        """Tests the results of both gas turbine models."""
        self.setup_CombustionChamber_model()
        self.setup_CombustionChamberStoich_model()
        m1 = round(self.nw1.get_conn('flue gas after cc').m.val, 6)
        m2 = round(self.nw2.get_conn('flue gas after cc').m.val, 6)
        msg = (
            'The outlet mass flow of the combustion chamber model is ' +
            str(m1) + ' while the outlet mass flow of the combustion chamber '
            'stoich model is ' + str(m2) + '. Both values should match.')

        assert m1 == m2, msg

        T1 = self.nw1.get_conn('flue gas after cc').T.val_SI
        T2 = self.nw2.get_conn('flue gas after cc').T.val_SI
        d_rel = abs(T2 - T1) / T1
        msg = (
            'The relative deviation in temperature after combustion is ' +
            str(d_rel) + ' with a maximum allowed value of 1e-3.')

        assert d_rel <= 1e-3, msg

        T1 = self.nw1.get_conn('flue gas after gt').T.val_SI
        T2 = self.nw2.get_conn('flue gas after gt').T.val_SI
        d_rel = abs(T2 - T1) / T1
        msg = (
            'The relative deviation in temperature after the turbine is ' +
            str(d_rel) + ' with a maximum allowed value of 1e-3.')

        assert d_rel <= 1e-3, msg

        shutil.rmtree('LUT', ignore_errors=True)
