# -*- coding: utf-8

"""Module for testing busses.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_busses.py

SPDX-License-Identifier: MIT
"""

from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.connections import Ref
from tespy.networks import Network


class TestConnections:

    def setup_method(self):
        """Set up the model."""
        # %% network setup
        fluid_list = ['Air']
        self.nw = Network(
            fluids=fluid_list, p_unit='bar', T_unit='C', p_range=[0.5, 20])

        # %% components
        so1 = Source('source 1')
        so2 = Source('source 2')
        si1 = Sink('sink 1')
        si2 = Sink('sink 2')

        # %% connections
        c1 = Connection(so1, 'out1', si1, 'in1', label='Some example label')
        c2 = Connection(so2, 'out1', si2, 'in1')

        self.nw.add_conns(c1, c2)

        # %% component parameters
        c1.set_attr(m=1, p=1, T=25, fluid={'Air': 1})
        c2.set_attr(m=0.5, p=10, T=25, fluid={'Air': 1})

        # %% solving
        self.nw.solve('design')

    def test_volumetric_flow_reference(self):
        """Test the referenced volumetric flow."""
        c1, c2 = self.nw.get_conn(
            ['Some example label', 'source 2:out1_sink 2:in1']
        )
        c2.set_attr(m=None, v=Ref(c1, 1, 0))

        self.nw.solve('design')

        # expected result: mass flow of second connection is lower by
        # the fraction of the specific volumes

        m_expected = round(c1.m.val * c1.vol.val / c2.vol.val, 4)
        m_is = round(c2.m.val, 4)
        msg = (
            'The mass flow of the connection 2 should be equal to '
            f'{m_expected} kg/s, but is {m_is} kg/s'
        )
        assert m_is == m_expected, msg
