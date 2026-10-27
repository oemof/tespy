# -*- coding: utf-8

"""Backwards compatibility of the data container names before version 0.12.

A custom component and a network are written in the old style only: the
former container classes, their aliases, the direct access to the
characteristic function object and the former reference attributes of
connections.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_tools/test_legacy_containers.py

SPDX-License-Identifier: MIT
"""
import warnings

import pytest

from tespy.components import Sink
from tespy.components import Source
from tespy.components.component import Component
from tespy.connections import Connection
from tespy.connections import Ref
from tespy.networks import Network
from tespy.tools.characteristics import CharLine
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import FluidProperties as dc_prop
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple

# the module exercises the deprecated names on purpose, the warnings of the
# constraint shims are raised on every solve and are not the subject here
pytestmark = pytest.mark.filterwarnings(
    "ignore:ComponentMandatoryConstraints is deprecated:FutureWarning"
)


class LegacyHeater(Component):
    """Heater written against the container API of tespy 0.11."""

    @staticmethod
    def inlets():
        return ["in1"]

    @staticmethod
    def outlets():
        return ["out1"]

    def get_mandatory_constraints(self):
        return {
            "mass_flow_constraints": dc_cmc(**{
                "structure_matrix": self.variable_equality_structure_matrix,
                "num_eq_sets": self.num_i,
                "func_params": {"variable": "m"}
            }),
            "fluid_constraints": dc_cmc(**{
                "structure_matrix": self.variable_equality_structure_matrix,
                "num_eq_sets": self.num_i,
                "func_params": {"variable": "fluid"}
            })
        }

    def get_parameters(self):
        return {
            "Q": dc_cp(
                num_eq_sets=1, func=self.Q_func,
                dependents=self.Q_dependents, quantity="heat",
                calc=self._calc_Q
            ),
            "pr": dc_cp(
                num_eq_sets=1, structure_matrix=self.pr_structure_matrix,
                func_params={"pr": "pr"}, quantity="ratio",
                calc=self._calc_pr
            ),
            "Q_design": dc_cp(quantity="heat"),
            "Q_char": dc_cc(param="m", char_params={"type": "abs"}),
            "Q_char_group": dc_gcp(
                elements=["Q_design", "Q_char"], num_eq_sets=1,
                func=self.Q_char_group_func, dependents=self.Q_dependents
            ),
            "num_sections": dc_simple(dtype="int"),
            "Q_loss": dc_cp(
                is_result=True, quantity="heat", calc=self._calc_Q_loss
            ),
        }

    def _calc_Q(self):
        i, o = self.inl[0], self.outl[0]
        return i.m.val_SI * (o.h.val_SI - i.h.val_SI)

    def _calc_Q_loss(self):
        return -0.01 * self._calc_Q()

    def Q_func(self):
        return self._calc_Q() - self.Q.val_SI

    def Q_dependents(self):
        i, o = self.inl[0], self.outl[0]
        return [i.m, i.h, o.h]

    def Q_char_group_func(self):
        expr = self.get_char_expr(self.Q_char.param, **self.Q_char.char_params)
        heat = self.Q_design.val_SI * self.Q_char.char_func.evaluate(expr)
        return self._calc_Q() - heat


@pytest.fixture
def legacy_network():
    with pytest.warns(FutureWarning):
        heater = LegacyHeater("heater")
    nw = Network(iterinfo=False)
    nw.units.set_defaults(
        pressure="bar", pressure_difference="bar", temperature="degC"
    )
    source = Source("source")
    sink = Sink("sink")
    c1 = Connection(source, "out1", heater, "in1", label="c1")
    c2 = Connection(heater, "out1", sink, "in1", label="c2")
    nw.add_conns(c1, c2)
    c1.set_attr(fluid={"water": 1}, m=1, p=1, T=20)
    return nw, heater, c1, c2


def test_legacy_component_solves(legacy_network):
    nw, heater, c1, c2 = legacy_network
    heater.set_attr(Q=1e4, num_sections=3)
    c2.set_attr(p=Ref(c1, 1, 0))
    nw.solve("design")
    nw.assert_convergence()
    assert c2.T.val > c1.T.val
    assert c2.p.val == pytest.approx(c1.p.val)
    assert heater.Q_loss.val == pytest.approx(-100)
    assert heater.num_sections.val == 3


def test_legacy_characteristic_group(legacy_network):
    nw, heater, c1, c2 = legacy_network
    heater.set_attr(
        Q_design=1e4, pr=1,
        Q_char={
            "char_func": CharLine(x=[0, 1, 2], y=[0.5, 1, 1.5]),
            "is_set": True
        }
    )
    nw.solve("design")
    nw.assert_convergence()
    assert heater.Q_char_group.is_set
    assert heater.Q.val == pytest.approx(1e4)
    c1.set_attr(m=2)
    nw.solve("design")
    nw.assert_convergence()
    assert heater.Q.val == pytest.approx(1.5e4)


def test_legacy_type_checks(legacy_network):
    nw, heater, c1, c2 = legacy_network
    assert isinstance(heater.Q, dc_cp)
    assert isinstance(heater.Q_char, dc_cc)
    assert isinstance(heater.Q_char_group, dc_gcp)
    assert isinstance(heater.num_sections, dc_simple)
    assert isinstance(c1.m, dc_prop)
    assert not isinstance(heater.Q, dc_prop)
    assert heater.Q_loss.result_only
    assert not heater.Q.result_only


def test_legacy_reference_attributes(legacy_network):
    nw, heater, c1, c2 = legacy_network
    c2.set_attr(m=Ref(c1, 1, 0))
    with pytest.warns(FutureWarning):
        assert c2.m_ref.is_set
    with pytest.warns(FutureWarning):
        assert c2.get_attr("m_ref").ref.obj is c1
    c2.set_attr(m=None)
    with pytest.warns(FutureWarning):
        assert not c2.m_ref.is_set


def test_legacy_simple_container_with_equation_is_optional_equation():
    """Components of 0.11 declared user switched equations such as cone or
    fluid_balance as SimpleDataContainer with an equation."""
    with pytest.warns(FutureWarning, match="OptionalEquation"):
        switch = dc_simple(
            dtype="bool", func=lambda: 0, _val=False, num_eq_sets=1,
            dependents=lambda: []
        )
    assert type(switch).__name__ == "OptionalEquation"
    assert not switch.is_set
    switch.set_attr(is_set=True)
    assert switch.is_set
    with pytest.warns(FutureWarning):
        config = dc_simple(val=5)
    assert config.val == 5


def test_legacy_container_dispatch_matches_new_instances():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        heater = LegacyHeater("heater")
    for key, container in heater.parameters.items():
        assert isinstance(container, (dc_cp, dc_cc, dc_gcp, dc_simple)), key
