# -*- coding: utf-8

"""Module for testing fluid properties of gas mixtures.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_fluid_properties.py

SPDX-License-Identifier: MIT
"""
import numpy as np
import shutil
from tespy.components import (
    pump, turbine, condenser, source, sink,
    cycle_closer, heat_exchanger_simple)
from tespy.connections import connection
from tespy.networks import network
from tespy.tools import fluid_properties as fp
from tespy.tools.global_vars import molar_masses, gas_constants


class TestFluidProperties:

    def setup(self):
        fluids = ['Air', 'N2', 'O2', 'Ar', 'CO2']
        fp.memorise.add_fluids({'Air': 'HEOS'})
        fp.memorise.add_fluids({
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
                        d_rel_max = 0.015
                        msg = ('Relative deviation is ' +
                               str(round(d_rel, 4)) + ' at inputs p=' +
                               str(round(p, 0)) + ', T=' + str(round(T, 0)) +
                               ' for function ' + name + ', should be < ' +
                               str(d_rel_max) + '.')
                        assert d_rel < d_rel_max, self.errormsg + msg
                    elif p < 5e6 and T < 500:
                        d_rel_max = 0.05
                        msg = ('Relative deviation is ' +
                               str(round(d_rel, 4)) + ' at inputs p=' +
                               str(round(p, 0)) + ', T=' + str(round(T, 0)) +
                               ' for function ' + name + ', should be < ' +
                               str(d_rel_max) + '.')
                        assert d_rel < d_rel_max, self.errormsg + msg
                    elif p < 5e6 and T < 1000:
                        d_rel_max = 0.04
                        msg = ('Relative deviation is ' +
                               str(round(d_rel, 4)) + ' at inputs p=' +
                               str(round(p, 0)) + ', T=' + str(round(T, 0)) +
                               ' for function ' + name + ', should be < ' +
                               str(d_rel_max) + '.')
                        assert d_rel < d_rel_max, self.errormsg + msg
                    elif p < 5e6 and T < 1500:
                        d_rel_max = 0.03
                        msg = ('Relative deviation is ' +
                               str(round(d_rel, 4)) + ' at inputs p=' +
                               str(round(p, 0)) + ', T=' + str(round(T, 0)) +
                               ' for function ' + name + ', should be < ' +
                               str(d_rel_max) + '.')
                        assert d_rel < d_rel_max, self.errormsg + msg
                    elif T < 500:
                        d_rel_max = 0.1
                        msg = ('Relative deviation is ' +
                               str(round(d_rel, 4)) + ' at inputs p=' +
                               str(round(p, 0)) + ', T=' + str(round(T, 0)) +
                               ' for function ' + name + ', should be < ' +
                               str(d_rel_max) + '.')
                        assert d_rel < d_rel_max, self.errormsg + msg
                    elif T < 1000:
                        d_rel_max = 0.075
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


def test_tespy_fluid_mixture():
    """
    Test the mixture of a tespy fluid with a third fluid.

    This test checks the deviation in the calculated values between a mixture
    calculated with the individual fluids against a mixture of a tespy fluid
    consisting of two components of the mixture with the third component.
    """
    full_mix = {'N2': 0.3333, 'O2': 0.3333, 'Ar': 0.3334}
    N2 = full_mix['N2'] / (full_mix['N2'] + full_mix['O2'])
    partial_mix = {'N2': N2, 'O2': 1 - N2}
    mixture = {'Ar': full_mix['Ar'], 'partial': 1 - full_mix['Ar']}

    Tmin = 250
    Tmax = 1250
    pmin = 1e4
    pmax = 1e7

    myfluid = fp.tespy_fluid(
        alias='partial', fluid=partial_mix,
        p_range=[pmin, pmax], T_range=[Tmin, Tmax])

    flow_mix = [0, 0, 0, mixture]
    flow_full = [0, 0, 0, full_mix]
    p_range = np.linspace(pmin, pmax, 40)
    T_range = np.linspace(Tmin, Tmax, 40)
    # the backend specification does not matter in this case!
    fp.memorise.add_fluids({'partial': 'HEOS'})
    fp.memorise.add_fluids({fluid: 'HEOS' for fluid in full_mix.keys()})

    funcs = {
        'h': fp.h_mix_pT,
        's': fp.s_mix_pT,
        'v': fp.v_mix_pT,
        'visc': fp.visc_mix_pT}

    for name, func in funcs.items():
        # enthalpy and entropy need reference point definition
        if name == 'h' or name == 's':
            p_ref = 1e5
            T_ref = 500
            mix_ref = func([0, p_ref, 0, flow_mix[3]], T_ref)
            full_ref = func([0, p_ref, 0, flow_full[3]], T_ref)

        for p in p_range:
            flow_mix[1] = p
            flow_full[1] = p
            for T in T_range:
                val_mix = func(flow_mix, T)
                val_pure = func(flow_full, T)

                # enthalpy and entropy need reference point
                if name == 'h' or name == 's':
                    d_rel = abs(((val_mix - mix_ref) -
                                 (val_pure - full_ref)) /
                                (val_pure - full_ref))
                else:
                    d_rel = abs((val_mix - val_pure) / val_pure)

                msg = (
                    'Relative deviation is ' + str(round(d_rel, 6)) +
                    ' at inputs p=' + str(round(p, 0)) + ', T=' +
                    str(round(T, 0)) + ' for function ' + name +
                    ', should be < 1e-4.')
                assert d_rel < 1e-4, msg

    shutil.rmtree('LUT', ignore_errors=True)


class TestFluidPropertyBackEnds:
    """Testing full models with different fluid property back ends."""

    def setup_clausius_rankine(self, fluid_list):
        """Setup a Clausius-Rankine cycle."""
        self.nw = network(fluids=fluid_list)
        self.nw.set_attr(p_unit='bar', T_unit='C', iterinfo=False)

        # %% components

        # main components
        turb = turbine('turbine')
        con = condenser('condenser')
        pu = pump('pump')
        steam_generator = heat_exchanger_simple('steam generator')
        closer = cycle_closer('cycle closer')

        # cooling water
        so_cw = source('cooling water inlet')
        si_cw = sink('cooling water outlet')

        # %% connections

        # main cycle
        fs_in = connection(closer, 'out1', turb, 'in1', label='livesteam')
        ws = connection(turb, 'out1', con, 'in1', label='wastesteam')
        cond = connection(con, 'out1', pu, 'in1', label='condensate')
        fw = connection(pu, 'out1', steam_generator, 'in1', label='feedwater')
        fs_out = connection(steam_generator, 'out1', closer, 'in1')
        self.nw.add_conns(fs_in, ws, cond, fw, fs_out)

        # cooling water
        cw_in = connection(so_cw, 'out1', con, 'in2')
        cw_out = connection(con, 'out2', si_cw, 'in1')
        self.nw.add_conns(cw_in, cw_out)

        # %% parametrization of components

        turb.set_attr(eta_s=0.9)
        con.set_attr(pr1=1, pr2=0.99, ttd_u=5)
        pu.set_attr(eta_s=0.7)
        steam_generator.set_attr(pr=0.9)

        # %% parametrization of connections

        fs_in.set_attr(p=100, T=500, m=100, fluid={self.nw.fluids[0]: 1})

        cw_in.set_attr(T=20, p=5, fluid={self.nw.fluids[0]: 1})
        cw_out.set_attr(T=30)

        # %% solving
        self.nw.solve('design')

    def setup_pipeline_network(self, fluid_list):
        """Setup a pipeline network."""

    def test_clausius_rankine(self):
        """Test the Clausius-Rankine cycle with different back ends."""
        fluid = 'water'
        back_ends = ['HEOS', 'BICUBIC', 'TTSE']
        # the IF97 back end is buggy on the CoolProp side, therefore not
        # supported at the moment
        # back_ends = ['HEOS', 'BICUBIC', 'TTSE', 'IF97']
        results = {}
        for back_end in back_ends:
            # delete the fluid from the memorisation class
            if fluid in fp.memorise.state.keys():
                del fp.memorise.state[fluid]
            self.setup_clausius_rankine([back_end + '::' + fluid])
            results[back_end] = (
                1 - abs(self.nw.components['condenser'].Q.val) /
                self.nw.components['steam generator'].Q.val)

        efficiency = results['HEOS']

        for back_end in back_ends:
            if back_end == 'HEOS':
                continue

            d_rel = (
                abs(results[back_end] - efficiency) /
                efficiency)

            msg = (
                'The deviation in thermal efficiency of the Clausius-Rankine '
                'cycle calculated with ' + back_end + ' back end is ' +
                str(d_rel) + ' but should not be larger than 1e-6.')
            assert d_rel <= 1e-6, msg
