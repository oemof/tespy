# -*- coding: utf-8

"""Module of the base classes of the sorption components.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/sorption/base.py

SPDX-License-Identifier: MIT
"""

from tespy.components.nodes.base import NodeBase
from tespy.tools.data_containers import ComponentParameter
from tespy.tools.data_containers import MandatoryConstraint
from tespy.tools.fluid_properties.solutions import get_solution_pair
from tespy.tools.global_vars import FLUID_ALIASES


class SorptionBase(NodeBase):
    r"""
    Shared equations of the absorption cycle sorption components.

    The components are written against the working pair of the solution
    connections (:class:`.SolutionConnection`), e.g. LiBr-water or
    NH3-water. The mass fraction :math:`w` is the mass fraction of the key
    component of the pair (LiBr, NH3). The equations of all sorption
    components are

    - total mass balance over the sorption ports,
    - key component balance over the sorption ports, a vapor without the
      key component contributes zero,
    - pressure equality of the sorption ports,
    - the solution outlet at the bubble point, approached by the optional
      temperature difference :math:`\Delta T_\text{eq}` (:code:`td_eq`),

      .. math::

          0 = h_\text{out} - h\left(p, T_\text{bubble}\left(p,
          w_\text{out}\right) - \Delta T_\text{eq}, w_\text{out}\right)

    - the fluid closure of the solution outlet if its composition is
      variable.

    Subclasses define the ports through :code:`_sorption_inlets`,
    :code:`_sorption_outlets`, :code:`_solution_inlet` and
    :code:`_solution_outlet`.
    """

    def get_parameters(self):
        return {
            "Q": ComponentParameter(
                func=self.heat_func,
                dependents=self.heat_dependents,
                num_eq_sets=1,
                quantity="heat",
                description="heat transferred into the component",
            ),
            "td_eq": ComponentParameter(
                min_val=0,
                func=self.saturation_func,
                dependents=self.saturation_dependents,
                num_eq_sets=1,
                quantity="temperature_difference",
                description=(
                    "temperature difference of the solution outlet to its "
                    "bubble point"
                ),
            ),
        }

    def _sorption_constraints(self):
        return {
            "key_balance_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.key_balance_func,
                dependents=self.key_balance_dependents,
                description="key component mass balance",
            ),
            "saturation_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.saturation_func,
                dependents=self.saturation_dependents,
                description="solution outlet at the bubble point",
            ),
            "fluid_balance_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.fluid_balance_func,
                dependents=self.fluid_balance_dependents,
                description="fluid closure of the solution outlet",
            ),
        }

    def _sorption_inlets(self):
        return self.inl

    def _sorption_outlets(self):
        return self.outl

    @property
    def pair(self):
        """Working pair of the solution outlet."""
        return get_solution_pair(self._solution_outlet().mixing_rule)

    def _key_alias(self, c):
        aliases = FLUID_ALIASES.get_fluid(self.pair.key) & set(c.fluid.val)
        return next(iter(aliases), None)

    def _key_fraction(self, c):
        alias = self._key_alias(c)
        if alias is None:
            return 0.0
        return c.fluid.val[alias]

    def _key_dependents(self, c):
        alias = self._key_alias(c)
        if alias is None:
            return {}
        return {c.fluid: {alias} & c.fluid.is_var}

    def mass_flow_func(self):
        r"""
        Equation for the mass balance of the sorption ports.

        .. math::

            0 = \sum \dot m_\text{in} - \sum \dot m_\text{out}
        """
        return (
            sum(c.m.val_SI for c in self._sorption_inlets())
            - sum(c.m.val_SI for c in self._sorption_outlets())
        )

    def mass_flow_dependents(self):
        return [c.m for c in self._sorption_inlets() + self._sorption_outlets()]

    def pressure_structure_matrix(self, k):
        r"""
        Equations for pressure equality of the sorption ports.

        .. math::

            0 = p_\text{first} - p_i \; \forall i \in \text{other ports}
        """
        conns = self._sorption_inlets() + self._sorption_outlets()
        first = conns[0]
        for eq, conn in enumerate(conns[1:]):
            self._structure_matrix[k + eq, first.p.sm_col] = 1
            self._structure_matrix[k + eq, conn.p.sm_col] = -1

    def key_balance_func(self):
        r"""
        Equation for the mass balance of the key component.

        .. math::

            0 = \sum \dot m_\text{in} w_\text{in}
            - \sum \dot m_\text{out} w_\text{out}
        """
        return (
            sum(
                c.m.val_SI * self._key_fraction(c)
                for c in self._sorption_inlets()
            )
            - sum(
                c.m.val_SI * self._key_fraction(c)
                for c in self._sorption_outlets()
            )
        )

    def key_balance_dependents(self):
        # a stream without the key component, e.g. water vapor of a
        # non-volatile absorbent, does not enter the balance
        conns = [
            c for c in self._sorption_inlets() + self._sorption_outlets()
            if self._key_alias(c) is not None
        ]
        vectors = {}
        for c in conns:
            vectors.update(self._key_dependents(c))
        return {
            "scalars": [[c.m for c in conns]],
            "vectors": [vectors],
        }

    def saturation_func(self):
        r"""
        Equation for the solution outlet at its bubble point.

        .. math::

            0 = h_\text{out} - h\left(p, T_\text{bubble}\left(p,
            w_\text{out}\right) - \Delta T_\text{eq}, w_\text{out}\right)
        """
        c = self._solution_outlet()
        pair = self.pair
        w = self._key_fraction(c)
        T = pair.T_bubble(c.p.val_SI, w)
        if self.td_eq.is_set:
            T -= self.td_eq.val_SI
        return c.h.val_SI - pair.h_pT(c.p.val_SI, T, w)

    def saturation_dependents(self):
        c = self._solution_outlet()
        return {
            "scalars": [[c.p, c.h]],
            "vectors": [self._key_dependents(c)],
        }

    def fluid_balance_func(self):
        r"""
        Equation for the fluid closure of the solution outlet.

        .. math::

            0 = 1 - \sum_i x_{i,\text{out}}
        """
        c = self._solution_outlet()
        return sum(c.fluid.val.values()) - 1

    def fluid_balance_dependents(self):
        c = self._solution_outlet()
        return {"scalars": [[]], "vectors": [{c.fluid: c.fluid.is_var}]}

    def heat_func(self):
        return self.Q.val_SI - self.calc_Q()

    def heat_dependents(self):
        conns = self._sorption_inlets() + self._sorption_outlets()
        return [v for c in conns for v in (c.m, c.h)]

    def calc_Q(self):
        return (
            sum(c.m.val_SI * c.h.val_SI for c in self._sorption_outlets())
            - sum(c.m.val_SI * c.h.val_SI for c in self._sorption_inlets())
        )

    def _update_num_eq(self):
        c = self._solution_outlet()
        self.constraints["fluid_balance_constraints"].num_eq = (
            1 if c.fluid.is_var else 0
        )
        # the approach temperature replaces the equilibrium equation
        self.constraints["saturation_constraints"].num_eq = (
            0 if self.td_eq.is_set else 1
        )

    def _initial_affine_edges(self):
        conns = self._sorption_inlets() + self._sorption_outlets()
        first = conns[0]
        edges = [(first.p, c.p, 1.0, 0.0) for c in conns[1:]]
        edges += [(
            self._solution_inlet().h, self._solution_outlet().h, 1.0, 0.0, 0.1
        )]
        return edges

    def _initial_temperature_edges(self):
        return [(self._solution_inlet(), self._solution_outlet(), 0.0, 0.3)]

    def calc_parameters(self):
        self.Q.val_SI = self.calc_Q()
        c = self._solution_outlet()
        w = self._key_fraction(c)
        self.td_eq.val_SI = self.pair.T_bubble(c.p.val_SI, w) - c.T.val_SI


class TwoSidedSorptionBase(SorptionBase):
    r"""
    Shared equations of sorption components with a heat transfer fluid side.

    On top of the sorption side equations of :class:`SorptionBase`, the heat
    transfer fluid side adds

    - mass flow and fluid composition equality between its inlet
      :code:`in2` and outlet :code:`out2`,
    - an overall energy balance, all heat is exchanged between both sides.

    Optional parameters are the heat flow :code:`Q` (transferred from the hot
    side, thus always negative in analogy to the two-sided heat exchangers),
    the terminal temperature differences :code:`ttd_u` and :code:`ttd_l`
    and the heat transfer fluid side pressure specifications :code:`pr2` and
    :code:`dp2`.
    """

    def get_parameters(self):
        params = super().get_parameters()
        params.update({
            "Q": ComponentParameter(
                max_val=0,
                func=self.heat_func,
                dependents=self.heat_dependents,
                num_eq_sets=1,
                quantity="heat",
                description="heat transferred from the hot side",
            ),
            "ttd_u": ComponentParameter(
                min_val=0,
                func=self.ttd_u_func,
                dependents=self.ttd_u_dependents,
                num_eq_sets=1,
                quantity="temperature_difference",
                description="upper terminal temperature difference",
            ),
            "ttd_l": ComponentParameter(
                min_val=0,
                func=self.ttd_l_func,
                dependents=self.ttd_l_dependents,
                num_eq_sets=1,
                quantity="temperature_difference",
                description="lower terminal temperature difference",
            ),
            "pr2": ComponentParameter(
                min_val=1e-4, max_val=1, num_eq_sets=1,
                structure_matrix=self.pr_structure_matrix,
                func_params={"pr": "pr2", "inconn": 1, "outconn": 1},
                quantity="ratio",
                description="heat transfer fluid side outlet to inlet pressure ratio",
            ),
            "dp2": ComponentParameter(
                min_val=0, max_val=1e15, num_eq_sets=1,
                structure_matrix=self.dp_structure_matrix,
                func_params={"dp": "dp2", "inconn": 1, "outconn": 1},
                quantity="pressure_difference",
                description="heat transfer fluid side inlet to outlet absolute pressure change",
            ),
        })
        return params

    def _two_sided_constraints(self):
        return {
            "htf_mass_flow_constraints": MandatoryConstraint(
                num_eq_sets=1,
                structure_matrix=self.htf_variable_equality_structure_matrix,
                func_params={"variable": "m"},
                description="heat transfer fluid side mass flow equality",
            ),
            "htf_fluid_constraints": MandatoryConstraint(
                num_eq_sets=1,
                structure_matrix=self.htf_variable_equality_structure_matrix,
                func_params={"variable": "fluid"},
                description="heat transfer fluid side fluid composition equality",
            ),
            "energy_balance_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.energy_balance_func,
                dependents=self.energy_balance_dependents,
                description="energy balance between both sides",
            ),
        }

    def htf_variable_equality_structure_matrix(self, k, **kwargs):
        variable = kwargs.get("variable")
        self._structure_matrix[k, self.inl[1].get_attr(variable).sm_col] = 1
        self._structure_matrix[k, self.outl[1].get_attr(variable).sm_col] = -1

    def energy_balance_func(self):
        r"""
        Equation for the energy balance of both sides.

        .. math::

            0 = \sum \dot m_\text{out} h_\text{out}
            - \sum \dot m_\text{in} h_\text{in}
        """
        return (
            sum(c.h.val_SI * c.m.val_SI for c in self.outl)
            - sum(c.h.val_SI * c.m.val_SI for c in self.inl)
        )

    def energy_balance_dependents(self):
        return [v for c in self.inl + self.outl for v in (c.m, c.h)]

    def heat_dependents(self):
        return [self.inl[1].m, self.inl[1].h, self.outl[1].h]

    def _htf_heat(self):
        i = self.inl[1]
        o = self.outl[1]
        return i.m.val_SI * (o.h.val_SI - i.h.val_SI)

    def _separate_flat_enthalpy_starts(self, seeded):
        i = self.inl[1]
        o = self.outl[1]
        if i.h._reference_container is o.h._reference_container:
            return 0
        delta = 1e3
        if abs(o.h.val_SI - i.h.val_SI) >= delta:
            return 0
        if o.h.is_var and o.h._reference_container not in seeded:
            o.h.set_reference_val_SI(i.h.val_SI + self._htf_dh_sign * delta)
        elif i.h.is_var and i.h._reference_container not in seeded:
            i.h.set_reference_val_SI(o.h.val_SI - self._htf_dh_sign * delta)
        else:
            return 0
        return 1

    def _initial_affine_edges(self):
        edges = super()._initial_affine_edges()
        edges += [
            (self.inl[1].p, self.outl[1].p, 1.0, 0.0),
            (self.inl[1].h, self.outl[1].h, 1.0, 0.0),
        ]
        return edges

    def _ttd_temperature_edges(self, upper, lower):
        edges = super()._initial_temperature_edges()
        edges += [(self.inl[1], self.outl[1], 0.0, 1.0)]
        ttd_upper = self.ttd_u.val_SI if self.ttd_u.is_set else 10.0
        weight_upper = 5.0 if self.ttd_u.is_set else 0.3
        ttd_lower = self.ttd_l.val_SI if self.ttd_l.is_set else 10.0
        weight_lower = 5.0 if self.ttd_l.is_set else 0.3
        edges += [
            (*upper, ttd_upper, weight_upper),
            (*lower, ttd_lower, weight_lower),
        ]
        return edges

    def calc_parameters(self):
        super().calc_parameters()
        self.pr2.val_SI = self.outl[1].p.val_SI / self.inl[1].p.val_SI
        self.dp2.val_SI = self.inl[1].p.val_SI - self.outl[1].p.val_SI
