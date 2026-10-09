# -*- coding: utf-8

"""Module for testing components of type separator.
This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_components/test_separator.py
SPDX-License-Identifier: MIT
"""

from pytest import approx

from tespy.components import Separator
from tespy.components import Sink
from tespy.components import Source
from tespy.connections import Connection
from tespy.connections import Ref
from tespy.networks import Network


class TestSeparator:

    def setup_method(self):
        self.nwk = Network()
        self.nwk.units.set_defaults(**{
            "pressure": "bar", "pressure_difference": "bar",
            "temperature": "degC", "enthalpy": "kJ/kg"
        })
        self.nwk.iterinfo = False

        so = Source("Source")
        sep = Separator("Separator")
        si1 = Sink("Sink 1")
        si2 = Sink("Sink 2")

        c1 = Connection(so, "out1", sep, "in1", label="1")
        c2 = Connection(sep, "out1", si1, "in1", label="2")
        c3 = Connection(sep, "out2", si2, "in1", label="3")

        self.nwk.add_conns(c1, c2, c3)

    def test_coupled_enthalpies(self):
        """This tests if the energy balance equation handles linear dependent
        enthalpy at one inlet and one outlet"""
        c1, c2, c3 = self.nwk.get_conn(["1", "2", "3"])
        c1.set_attr(fluid={"N2": 0.5, "O2": 0.5}, m=5, p=10, T=50)
        c2.set_attr(fluid={"N2": 1, "O2": 0})
        c3.set_attr(fluid={"O2": 1, "N2": 0})

        self.nwk.solve("design")
        self.nwk.assert_convergence()
        assert self.nwk.status == 0
        assert c2.T.val == approx(c1.T.val)

        c1.set_attr(T=None)
        c2.set_attr(h=Ref(c1, 1, 23))
        self.nwk.solve("design", max_iter=500)
        self.nwk.assert_convergence()
        assert c2.T.val_SI == approx(c1.T.val_SI, abs=1e-3)
        assert c2.T.val == approx(102.007, abs=1e-3)

    def test_outlet_composition_omits_fluid(self):
        """An outlet composition that sums to one without naming every fluid
        of the branch must work: the presolver fixes the omitted fluid at
        zero and removes it from the fluid vector of that outlet"""
        c1, c2, c3 = self.nwk.get_conn(["1", "2", "3"])
        c1.set_attr(fluid={"N2": 0.7, "O2": 0.2, "Ar": 0.1}, m=10, p=1, T=20)
        c2.set_attr(fluid={"N2": 0.5, "O2": 0.5}, m=1)

        self.nwk.solve("design")
        self.nwk.assert_convergence()
        assert c2.fluid.val.get("Ar", 0) == 0
        assert c3.m.val_SI == approx(9)
        assert c3.fluid.val["N2"] == approx(6.5 / 9)
        assert c3.fluid.val["O2"] == approx(1.5 / 9)
        assert c3.fluid.val["Ar"] == approx(1 / 9)
