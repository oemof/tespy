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
        self.nw = Network()
        self.nw.units.set_defaults(**{
            "pressure": "bar", "temperature": "C", "volumetric_flow": "l/s",
            "mass_flow": "t/h"
        })

        so1 = Source('source 1')
        so2 = Source('source 2')
        si1 = Sink('sink 1')
        si2 = Sink('sink 2')

        c1 = Connection(so1, 'out1', si1, 'in1', label='Some example label')
        c2 = Connection(so2, 'out1', si2, 'in1')

        self.nw.add_conns(c1, c2)

        c1.set_attr(m=1, p=1, T=25, fluid={'Air': 1})
        c2.set_attr(m=0.5, p=10, T=25, fluid={'Air': 1})

        self.nw.solve('design')

    def test_volumetric_flow_reference(self):
        """Test the referenced volumetric flow."""
        c1, c2 = self.nw.get_conn(
            ['Some example label', 'source 2:out1_sink 2:in1']
        )
        c2.set_attr(m=None, v=Ref(c1, 1, 0))
        self.nw.solve('design')

        m_expected = round(c1.m.val * c1.vol.val / c2.vol.val, 4)
        m_is = round(c2.m.val, 4)
        msg = (
            'The mass flow of the connection 2 should be equal to '
            f'{m_expected} kg/s, but is {m_is} kg/s'
        )
        assert m_is == m_expected, msg

        c2.set_attr(v=Ref(c1, 2, 10))
        self.nw.solve('design')

        v_expected = round(c1.v.val * 2 + 10, 4)
        v_is = round(c2.v.val, 4)
        msg = (
            'The mass flow of the connection 2 should be equal to '
            f'{v_expected} l/s, but is {v_is} l/s'
        )
        assert v_is == v_expected, msg

    def test_temperature_reference(self):
        """Test the referenced temperature."""
        c1, c2 = self.nw.get_conn(
            ['Some example label', 'source 2:out1_sink 2:in1']
        )
        c2.set_attr(T=None)
        c2.set_attr(T=Ref(c1, 1, 0))

        self.nw.solve('design')

        T_expected = round(c1.T.val, 4)
        T_is = round(c2.T.val, 4)
        msg = (
            'The temperature of the connection 2 should be equal to '
            f'{T_expected} C, but is {T_is} C'
        )
        assert T_is == T_expected, msg

        delta = -75
        c2.set_attr(T=Ref(c1, 1.5, delta))
        self.nw.solve('design')

        delta_SI = self.nw.units.ureg.Quantity(
            delta, self.nw.units.default["temperature_difference"]
        ).to_base_units().magnitude
        assert round(delta_SI, 4) == round(c2.T_ref.ref.delta_SI, 4)

        T_expected = round(c1.T.val_SI * 1.5 + delta_SI, 4)
        T_is = round(c2.T.val_SI, 4)
        msg = (
            'The temperature of the connection 2 should be equal to '
            f'{T_expected} C, but is {T_is} C'
        )
        assert T_is == T_expected, msg

    def test_primary_reference(self):
        """Test referenced primary variable."""
        c1, c2 = self.nw.get_conn(
            ['Some example label', 'source 2:out1_sink 2:in1']
        )
        c2.set_attr(m=None)
        c2.set_attr(m=Ref(c1, 1, 0))

        self.nw.solve('design')

        m_expected = round(c1.m.val, 4)
        m_is = round(c2.m.val, 4)
        msg = (
            'The mass flow of the connection 2 should be equal to '
            f'{m_expected} kg/s, but is {m_is} kg/s'
        )
        assert m_is == m_expected, msg

        delta = -0.5
        c2.set_attr(m=Ref(c1, 2, delta))
        self.nw.solve('design')
        delta_SI = self.nw.units.ureg.Quantity(
            delta, self.nw.units.default["mass_flow"]
        ).to_base_units().magnitude
        assert round(delta_SI, 4) == round(c2.m_ref.ref.delta_SI, 4)

        m_expected = round(c1.m.val_SI * 2 + delta_SI, 4)
        m_is = round(c2.m.val_SI, 4)
        msg = (
            'The mass flow of the connection 2 should be equal to '
            f'{m_expected} kg/s, but is {m_is} kg/s'
        )
        assert m_is == m_expected, msg
