# -*- coding: utf-8

"""Module for testing network properties.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_networks/test_exergy_and_entropy_analysis.py

SPDX-License-Identifier: MIT
"""
from tespy.components import CycleCloser
from tespy.components import HeatExchangerSimple
from tespy.components import Merge
from tespy.components import Pump
from tespy.components import Splitter
from tespy.components import Turbine
from tespy.connections import Bus
from tespy.connections import Connection
from tespy.networks import Network


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
        self.nw = Network(fluids=fluids)
        self.nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')

        # create components
        splitter1 = Splitter('splitter 1')
        merge1 = Merge('merge 1')
        turb = Turbine('turbine')
        fwp_turb = Turbine('feed water pump turbine')
        condenser = HeatExchangerSimple('condenser')
        fwp = Pump('pump')
        steam_generator = HeatExchangerSimple('steam generator')
        cycle_close = CycleCloser('cycle closer')

        # create busses
        # power output bus
        self.power = Bus('power_output')
        self.power.add_comps({'comp': turb, 'char': 1})
        # turbine driven feed water pump internal bus
        self.fwp_power = Bus('feed water pump power', P=0)
        self.fwp_power.add_comps(
            {'comp': fwp_turb, 'char': 1},
            {'comp': fwp, 'char': 1, 'base': 'bus'})
        # heat input bus
        self.heat = Bus('heat_input')
        self.heat.add_comps({'comp': steam_generator, 'base': 'bus'})
        self.nw.add_busses(self.power, self.fwp_power, self.heat)

        # create connections
        fs_in = Connection(cycle_close, 'out1', splitter1, 'in1', label='fs')
        fs_fwpt = Connection(splitter1, 'out1', fwp_turb, 'in1')
        fs_t = Connection(splitter1, 'out2', turb, 'in1')
        fwpt_ws = Connection(fwp_turb, 'out1', merge1, 'in1')
        t_ws = Connection(turb, 'out1', merge1, 'in2')
        ws = Connection(merge1, 'out1', condenser, 'in1')
        cond = Connection(condenser, 'out1', fwp, 'in1', label='cond')
        fw = Connection(fwp, 'out1', steam_generator, 'in1', label='fw')
        fs_out = Connection(steam_generator, 'out1', cycle_close, 'in1')
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
        for cp in self.nw.comps['object']:
            cp.entropy_balance()
        convergence_check(self.nw.lin_dep)

    def test_entropy_perfect_cycle(self):
        """Test entropy values in the perfect clausius rankine cycle."""
        labels = [
            'turbine', 'feed water pump turbine', 'condenser',
            'steam generator', 'pump'
        ]
        for label in labels:
            cp = self.nw.get_comp(label)
            msg = (
                'Entropy production due to irreversibility must be 0 for all '
                'components in this test but is ' + str(round(cp.S_irr, 4)) +
                ' at component ' + label + ' of type ' + cp.component() + '.')
            assert round(cp.S_irr, 4) == 0, msg
        sg = self.nw.get_comp('steam generator')
        cd = self.nw.get_comp('condenser')
        msg = (
            'Value of entropy production due to heat input at steam generator '
            '(S_Q=' + str(round(sg.S_Q, 4)) + ') must equal the negative '
            'value of entropy reduction in condenser (S_Q=' +
            str(round(cd.S_Q, 4)) + ').')
        assert round(sg.S_Q, 4) == -round(cd.S_Q, 4), msg
