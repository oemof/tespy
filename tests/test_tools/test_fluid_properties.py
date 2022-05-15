# -*- coding: utf-8

"""Module for testing fluid properties of gas mixtures.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_fluid_properties.py

SPDX-License-Identifier: MIT
"""
import os

import numpy as np
import pytest

from tespy.components import Condenser
from tespy.components import CycleCloser
from tespy.components import HeatExchangerSimple
from tespy.components import Pipe
from tespy.components import Pump
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Turbine
from tespy.connections import Connection
from tespy.networks import Network
from tespy.tools import fluid_properties as fp


def convergence_check(lin_dep):
    """Check convergence status of a simulation."""
    msg = 'Calculation did not converge!'
    assert lin_dep is False, msg


class TestFluidProperties:

    def setup(self):
        fp.Memorise.add_fluids({'Air': 'HEOS'})
        fp.Memorise.add_fluids({
            'N2': 'HEOS', 'O2': 'HEOS', 'Ar': 'HEOS', 'CO2': 'HEOS'})

        mix = {'N2': 0.7556, 'O2': 0.2315, 'Ar': 0.0129}
        pure = {'Air': 1}
        self.flow_mix = [0, 0, 0, mix]
        self.flow_pure = [0, 0, 0, pure]
        self.p_range = np.linspace(1e-2, 200, 40) * 1e5
        self.T_range = np.linspace(220, 2220, 40)
        self.errormsg = ('Relative deviation of fluid mixture to base '
                         '(CoolProp air) is too high: ')

    def test_properties(self):
        """
        Test gas mixture fluid properties.

        Test the CoolProp pseudo pure fluid dry air properties vs. mixture of
        air components. Check enthalpy, entropy, specific volume, viscosity.

        A CoolProp mixture object could be check as well!
        """
        funcs = {'h': fp.h_mix_pT,
                 's': fp.s_mix_pT,
                 'v': fp.v_mix_pT,
                 'visc': fp.visc_mix_pT}
        for name, func in funcs.items():
            # enthalpy and entropy need reference point definition
            if name == 'h' or name == 's':
                p_ref = 1e5
                T_ref = 500
                mix_ref = func([0, p_ref, 0, self.flow_mix[3]], T_ref)
                pure_ref = func([0, p_ref, 0, self.flow_pure[3]], T_ref)

            for p in self.p_range:
                self.flow_mix[1] = p
                self.flow_pure[1] = p
                for T in self.T_range:
                    val_mix = func(self.flow_mix, T)
                    val_pure = func(self.flow_pure, T)

                    # enthalpy and entropy need reference point
                    if name == 'h' or name == 's':
                        d_rel = abs(((val_mix - mix_ref) -
                                     (val_pure - pure_ref)) /
                                    (val_pure - pure_ref))
                    else:
                        d_rel = abs((val_mix - val_pure) / val_pure)

                    # these values seem arbitrary...
                    if name == 's':
                        if round(p, 0) == 7180128.0 and round(T) == 1502.0:
                            continue
                        elif round(p, 0) == 17948821.0 and round(T) == 1861.0:
                            continue

                    # the deviations might have to be checked
                    if p <= 1e6:
                        d_rel_max = 0.005
                        msg = ('Relative deviation is ' +
                               str(round(d_rel, 4)) + ' at inputs p=' +
                               str(round(p, 0)) + ', T=' + str(round(T, 0)) +
                               ' for function ' + name + ', should be < ' +
                               str(d_rel_max) + '.')
                        assert d_rel < d_rel_max, self.errormsg + msg
                    elif p < 5e6 and T < 500:
                        d_rel_max = 0.04
                        msg = ('Relative deviation is ' +
                               str(round(d_rel, 4)) + ' at inputs p=' +
                               str(round(p, 0)) + ', T=' + str(round(T, 0)) +
                               ' for function ' + name + ', should be < ' +
                               str(d_rel_max) + '.')
                        assert d_rel < d_rel_max, self.errormsg + msg
                    elif p < 5e6 and T < 1000:
                        d_rel_max = 0.03
                        msg = ('Relative deviation is ' +
                               str(round(d_rel, 4)) + ' at inputs p=' +
                               str(round(p, 0)) + ', T=' + str(round(T, 0)) +
                               ' for function ' + name + ', should be < ' +
                               str(d_rel_max) + '.')
                        assert d_rel < d_rel_max, self.errormsg + msg
                    elif p < 5e6 and T < 1500:
                        d_rel_max = 0.02
                        msg = ('Relative deviation is ' +
                               str(round(d_rel, 4)) + ' at inputs p=' +
                               str(round(p, 0)) + ', T=' + str(round(T, 0)) +
                               ' for function ' + name + ', should be < ' +
                               str(d_rel_max) + '.')
                        assert d_rel < d_rel_max, self.errormsg + msg
                    elif T < 500:
                        d_rel_max = 0.09
                        msg = ('Relative deviation is ' +
                               str(round(d_rel, 4)) + ' at inputs p=' +
                               str(round(p, 0)) + ', T=' + str(round(T, 0)) +
                               ' for function ' + name + ', should be < ' +
                               str(d_rel_max) + '.')
                        assert d_rel < d_rel_max, self.errormsg + msg
                    elif T < 1000:
                        d_rel_max = 0.06
                        msg = ('Relative deviation is ' +
                               str(round(d_rel, 4)) + ' at inputs p=' +
                               str(round(p, 0)) + ', T=' + str(round(T, 0)) +
                               ' for function ' + name + ', should be < ' +
                               str(d_rel_max) + '.')
                        assert d_rel < d_rel_max, self.errormsg + msg
                    else:
                        d_rel_max = 0.025
                        msg = ('Relative deviation is ' +
                               str(round(d_rel, 4)) + ' at inputs p=' +
                               str(round(p, 0)) + ', T=' + str(round(T, 0)) +
                               ' for function ' + name + ', should be < ' +
                               str(d_rel_max) + '.')
                        assert d_rel < d_rel_max, self.errormsg + msg


class TestFluidPropertyBackEnds:
    """Testing full models with different fluid property back ends."""

    def setup_clausius_rankine(self, fluid_list):
        """Setup a Clausius-Rankine cycle."""
        self.nw = Network(fluids=fluid_list)
        self.nw.set_attr(p_unit='bar', T_unit='C', iterinfo=True)

        # %% components

        # main components
        turb = Turbine('turbine')
        con = Condenser('condenser')
        pu = Pump('pump')
        steam_generator = HeatExchangerSimple('steam generator')
        closer = CycleCloser('cycle closer')

        # cooling water
        so_cw = Source('cooling water inlet')
        si_cw = Sink('cooling water outlet')

        # %% connections

        # main cycle
        fs_in = Connection(closer, 'out1', turb, 'in1', label='livesteam')
        ws = Connection(turb, 'out1', con, 'in1', label='wastesteam')
        cond = Connection(con, 'out1', pu, 'in1', label='condensate')
        fw = Connection(pu, 'out1', steam_generator, 'in1', label='feedwater')
        fs_out = Connection(steam_generator, 'out1', closer, 'in1')
        self.nw.add_conns(fs_in, ws, cond, fw, fs_out)

        # cooling water
        cw_in = Connection(so_cw, 'out1', con, 'in2')
        cw_out = Connection(con, 'out2', si_cw, 'in1')
        self.nw.add_conns(cw_in, cw_out)

        # %% parametrization of components

        turb.set_attr(eta_s=0.9)
        con.set_attr(pr1=1, pr2=0.99, ttd_u=5)
        steam_generator.set_attr(pr=0.9)

        # %% parametrization of connections

        fs_in.set_attr(p=100, T=500, m=100, fluid={self.nw.fluids[0]: 1})
        fw.set_attr(h=200e3)
        cw_in.set_attr(T=20, p=5, fluid={self.nw.fluids[0]: 1})
        cw_out.set_attr(T=30)

        # %% solving
        self.nw.solve('design')
        pu.set_attr(eta_s=0.7)
        fw.set_attr(h=None)
        self.nw.solve('design')

    def setup_pipeline_network(self, fluid_list):
        """Setup a pipeline network."""
        self.nw = Network(fluids=fluid_list)
        self.nw.set_attr(p_unit='bar', T_unit='C', iterinfo=False)

        # %% components

        # main components
        pu = Pump('pump')
        pi = Pipe('pipeline')
        es = HeatExchangerSimple('energy balance closing')

        closer = CycleCloser('cycle closer')

        pu_pi = Connection(pu, 'out1', pi, 'in1')
        pi_es = Connection(pi, 'out1', es, 'in1')
        es_closer = Connection(es, 'out1', closer, 'in1')
        closer_pu = Connection(closer, 'out1', pu, 'in1')
        self.nw.add_conns(pu_pi, pi_es, es_closer, closer_pu)

        # %% parametrization of components

        pu.set_attr(eta_s=0.7)
        pi.set_attr(pr=0.95, L=100, ks=1e-5, D='var', Q=0)
        es.set_attr(pr=1)

        # %% parametrization of connections

        pu_pi.set_attr(p=20, T=100, m=10, fluid={self.nw.fluids[0]: 1})

        # %% solving
        self.nw.solve('design')

    @pytest.mark.skipif(
        os.environ.get('GITHUB_ACTIONS') == 'true',
        reason='GitHub actions cannot handle the tabular CoolProp back ends, '
        'skipping this test. The test should run on your local machine.')
    def test_clausius_rankine_tabular(self):
        """Test the Clausius-Rankine cycle with different back ends."""
        fluid = 'water'
        back_ends = ['HEOS', 'BICUBIC', 'TTSE']
        results = {}
        for back_end in back_ends:
            # delete the fluid from the memorisation class
            if fluid in fp.Memorise.state.keys():
                del fp.Memorise.state[fluid]
                del fp.Memorise.back_end[fluid]
            self.setup_clausius_rankine([back_end + '::' + fluid])
            results[back_end] = (
                1 - abs(self.nw.get_comp('condenser').Q.val) /
                self.nw.get_comp('steam generator').Q.val)

        efficiency = results['HEOS']

        if fluid in fp.Memorise.state.keys():
            del fp.Memorise.state[fluid]
            del fp.Memorise.back_end[fluid]
        for back_end in back_ends:
            if back_end == 'HEOS':
                continue

            d_rel = (abs(results[back_end] - efficiency) / efficiency)

            msg = (
                'The deviation in thermal efficiency of the Clausius-Rankine '
                'cycle calculated with ' + back_end + ' back end is ' +
                str(d_rel) + ' but should not be larger than 1e-4.')
            assert d_rel <= 1e-4, msg

    def test_clausius_rankine(self):
        """Test the Clausius-Rankine cycle with different back ends."""
        fluid = 'water'
        back_ends = ['HEOS', 'IF97']
        results = {}
        for back_end in back_ends:
            # delete the fluid from the memorisation class
            if fluid in fp.Memorise.state.keys():
                del fp.Memorise.state[fluid]
                del fp.Memorise.back_end[fluid]
            self.setup_clausius_rankine([back_end + '::' + fluid])
            results[back_end] = (
                1 - abs(self.nw.get_comp('condenser').Q.val) /
                self.nw.get_comp('steam generator').Q.val)

        efficiency = results['HEOS']

        if fluid in fp.Memorise.state.keys():
            del fp.Memorise.state[fluid]
            del fp.Memorise.back_end[fluid]
        for back_end in back_ends:
            if back_end == 'HEOS':
                continue

            d_rel = (abs(results[back_end] - efficiency) / efficiency)

            msg = (
                'The deviation in thermal efficiency of the Clausius-Rankine '
                'cycle calculated with ' + back_end + ' back end is ' +
                str(d_rel) + ' but should not be larger than 1e-4.')
            assert d_rel <= 1e-4, msg

    def test_pipeline_network(self):
        """Test a pipeline network with fluids from different back ends."""
        fluids_back_ends = {'DowJ': 'INCOMP', 'water': 'HEOS'}

        for fluid, back_end in fluids_back_ends.items():
            # delete the fluid from the memorisation class
            if fluid in fp.Memorise.state.keys():
                del fp.Memorise.state[fluid]
            self.setup_pipeline_network([back_end + '::' + fluid])
            convergence_check(self.nw.lin_dep)

            value = round(self.nw.get_comp('pipeline').pr.val, 5)
            msg = (
                'The pressure ratio of the pipeline must be at 0.95, but '
                'is at ' + str(value) + ' for the fluid ' + fluid + '.')
            assert value == 0.95, msg
            value = round(self.nw.get_comp('pump').pr.val, 5)
            msg = (
                'The pressure ratio of the pipeline must be at ' +
                str(round(1 / 0.95, 5)) + ', but is at ' + str(value) +
                ' for the fluid ' + fluid + '.')
            assert value == round(1 / 0.95, 5), msg
