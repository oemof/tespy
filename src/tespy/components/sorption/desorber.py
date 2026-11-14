# -*- coding: utf-8

"""Module of the desorber components.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/sorption/desorber.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.sorption.base import SorptionBase
from tespy.components.sorption.base import TwoSidedSorptionBase
from tespy.tools import logger
from tespy.tools.data_containers import ConfigParameter
from tespy.tools.data_containers import MandatoryConstraint
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.helpers import TESPyComponentError

_VAPOR_STATES = ["inlet-equilibrium", "outlet-equilibrium", "free"]


class _DesorberMixin:
    r"""
    Vapor outlet equations of the desorbers.

    The option :code:`vapor_state` selects the state of the vapor leaving
    the desorber:

    - :code:`"inlet-equilibrium"` (default): in equilibrium with the
      refrigerant-rich solution entering the desorber, i.e. at the bubble
      point temperature of the inlet composition, as in a counter-current
      generator.
    - :code:`"outlet-equilibrium"`: in equilibrium with the
      refrigerant-poor solution leaving the desorber.
    - :code:`"free"`: no equation, the vapor state is specified on the
      connection.

    The vapor temperature equals the equilibrium temperature
    :math:`T_\text{eq} = T_\text{bubble}\left(p, w_\text{eq}\right)`. If the
    absorbent is volatile, e.g. water in NH3-water, the vapor composition
    :math:`y` is the equilibrium composition and the vapor is saturated:

    .. math::

        0 = y - y_\text{eq}\left(p, w_\text{eq}\right)\\
        0 = h_\text{vapor} - h\left(p, Q=1, y\right)

    Otherwise the vapor is the pure refrigerant at the equilibrium
    temperature:

    .. math::

        0 = h_\text{vapor} - h\left(p, T_\text{eq}\right)
    """

    def _vapor_state_parameter(self):
        return {
            "vapor_state": ConfigParameter(
                _val="inlet-equilibrium", dtype="str", options=_VAPOR_STATES,
                description="state of the vapor leaving the desorber",
            )
        }

    def _vapor_constraints(self):
        return {
            "vapor_state_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.vapor_state_func,
                dependents=self.vapor_state_dependents,
                description="vapor outlet in equilibrium with the solution",
            ),
            "vapor_composition_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.vapor_composition_func,
                dependents=self.vapor_composition_dependents,
                description="vapor outlet composition in equilibrium",
            ),
            "vapor_fluid_balance_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.vapor_fluid_balance_func,
                dependents=self.vapor_fluid_balance_dependents,
                description="fluid closure of the vapor outlet",
            ),
        }

    def _equilibrium_solution(self):
        if self.vapor_state.val == "outlet-equilibrium":
            return self._solution_outlet()
        return self._solution_inlet()

    def vapor_state_func(self):
        r"""
        Equation for the vapor outlet in equilibrium with the solution.

        Returns
        -------
        residual : float
            Residual value of the equation, see :class:`_DesorberMixin`.
        """
        vapor = self._vapor_outlet()
        pair = self.pair
        p = vapor.p.val_SI
        if pair.volatile_absorbent:
            return vapor.h.val_SI - pair.h_pQ(p, 1, self._key_fraction(vapor))
        w = self._key_fraction(self._equilibrium_solution())
        T = pair.T_bubble(p, w)
        return vapor.h.val_SI - h_mix_pT(
            p, T, vapor.fluid_data, vapor.mixing_rule
        )

    def vapor_state_dependents(self):
        vapor = self._vapor_outlet()
        if self.pair.volatile_absorbent:
            vectors = self._key_dependents(vapor)
        else:
            vectors = self._key_dependents(self._equilibrium_solution())
        return {"scalars": [[vapor.p, vapor.h]], "vectors": [vectors]}

    def vapor_composition_func(self):
        r"""
        Equation for the vapor outlet composition in equilibrium.

        Returns
        -------
        residual : float
            Residual value of the equation, see :class:`_DesorberMixin`.
        """
        vapor = self._vapor_outlet()
        w = self._key_fraction(self._equilibrium_solution())
        return (
            self._key_fraction(vapor) - self.pair.y_eq(vapor.p.val_SI, w)
        )

    def vapor_composition_dependents(self):
        vapor = self._vapor_outlet()
        vectors = {}
        vectors.update(self._key_dependents(self._equilibrium_solution()))
        vectors.update(self._key_dependents(vapor))
        return {"scalars": [[vapor.p]], "vectors": [vectors]}

    def vapor_fluid_balance_func(self):
        r"""
        Equation for the fluid closure of the vapor outlet.

        .. math::

            0 = 1 - \sum_i y_{i}
        """
        return sum(self._vapor_outlet().fluid.val.values()) - 1

    def vapor_fluid_balance_dependents(self):
        vapor = self._vapor_outlet()
        return {"scalars": [[]], "vectors": [{vapor.fluid: vapor.fluid.is_var}]}

    def _update_num_eq(self):
        super()._update_num_eq()
        vapor = self._vapor_outlet()
        equilibrium = self.vapor_state.val != "free"
        variable_composition = equilibrium and self.pair.volatile_absorbent
        if variable_composition and not vapor.fluid.is_var:
            msg = (
                f"The vapor composition at {vapor.label} leaving the "
                f"desorber {self.label} is the equilibrium composition of "
                f"its vapor_state '{self.vapor_state.val}'. Remove the "
                "composition specification of the vapor or set the "
                "vapor_state to 'free'."
            )
            logger.error(msg)
            raise TESPyComponentError(msg)
        self.constraints["vapor_state_constraints"].num_eq = int(equilibrium)
        self.constraints["vapor_composition_constraints"].num_eq = int(
            variable_composition
        )
        self.constraints["vapor_fluid_balance_constraints"].num_eq = int(
            variable_composition
        )

    def _desorber_temperature_edges(self):
        return [(self._solution_outlet(), self._vapor_outlet(), 0.0, 0.3)]

    def start_fluid_wrapper_branch(self):
        branches = {}
        for outconn in [self._vapor_outlet(), self._solution_outlet()]:
            branch = {"connections": [outconn], "components": [self]}
            outconn.target.propagate_wrapper_to_target(branch)
            branches[outconn.label] = branch
        return branches


@component_registry
class Desorber(_DesorberMixin, SorptionBase):
    r"""
    Ideal desorber (generator) of absorption cycles.

    The desorber separates a refrigerant-rich solution (:code:`in1`) into a
    refrigerant-poor solution (:code:`out1`) at its bubble point and a
    refrigerant vapor (:code:`out2`). The working pair is the mixing rule of
    the solution connections, see :class:`.SolutionConnection`. For
    LiBr-water the vapor is pure water on a plain :class:`.Connection`, for
    NH3-water it is a :class:`.SolutionConnection` of the same pair.

    **Mandatory Equations**

    - mass balance :py:meth:`mass_flow_func`
    - key component balance :py:meth:`key_balance_func`
    - pressure equality :py:meth:`pressure_structure_matrix`
    - solution outlet at the bubble point :py:meth:`saturation_func`
    - fluid closure of the solution outlet :py:meth:`fluid_balance_func`
    - vapor outlet in equilibrium :py:meth:`vapor_state_func`, unless
      :code:`vapor_state` is :code:`"free"`
    - vapor composition in equilibrium and its fluid closure for pairs with
      a volatile absorbent :py:meth:`vapor_composition_func`,
      :py:meth:`vapor_fluid_balance_func`

    **Optional Equations**

    - heat flow :py:meth:`heat_func`
    - approach to the bubble point :code:`td_eq` replacing the bubble point
      equation, :py:meth:`saturation_func`

    Inlets/Outlets

    - in1: refrigerant-rich solution
    - out1: refrigerant-poor solution
    - out2: refrigerant vapor

    Parameters
    ----------
    label : str
        The label of the component.

    Q : float, dict
        Heat supplied to the desorber :math:`\dot Q`/W.

    td_eq : float, dict
        Temperature difference of the solution outlet below its bubble point
        :math:`\Delta T_\text{eq}`/K.

    vapor_state : str
        State of the vapor leaving the desorber,
        :code:`"inlet-equilibrium"` (default), :code:`"outlet-equilibrium"`
        or :code:`"free"`, see :class:`_DesorberMixin`.

    Example
    -------
    A LiBr solution of 55 % LiBr is heated to 87 degC at the saturation
    pressure of water at 45 degC. The vapor leaves in equilibrium with the
    entering solution.

    >>> from CoolProp.CoolProp import PropsSI
    >>> from tespy.components import Desorber, Sink, Source
    >>> from tespy.connections import Connection, SolutionConnection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> desorber = Desorber("desorber")
    >>> rich = SolutionConnection(
    ...     Source("rich"), "out1", desorber, "in1", mixing_rule="libr-water"
    ... )
    >>> poor = SolutionConnection(desorber, "out1", Sink("poor"), "in1")
    >>> vapor = Connection(desorber, "out2", Sink("vapor"), "in1")
    >>> nw.add_conns(rich, poor, vapor)
    >>> p_cd = PropsSI("P", "T", 318.15, "Q", 0, "water")
    >>> rich.set_attr(xi=0.55, p=p_cd, m=10, T=315)
    >>> poor.set_attr(T=360)
    >>> vapor.set_attr(fluid={"H2O": 1})
    >>> nw.solve("design")
    >>> round(poor.xi.val, 3)
    0.581
    >>> round(vapor.T.val - 273.15, 1)
    80.2
    """

    _is_wrapper_branch_source = True

    def get_parameters(self):
        return {**super().get_parameters(), **self._vapor_state_parameter()}

    @staticmethod
    def inlets():
        return ["in1"]

    @staticmethod
    def outlets():
        return ["out1", "out2"]

    def _solution_inlet(self):
        return self.inl[0]

    def _solution_outlet(self):
        return self.outl[0]

    def _vapor_outlet(self):
        return self.outl[1]

    def get_mandatory_constraints(self):
        return {
            "mass_flow_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.mass_flow_func,
                dependents=self.mass_flow_dependents,
                description="mass balance",
            ),
            "pressure_constraints": MandatoryConstraint(
                num_eq_sets=2,
                structure_matrix=self.pressure_structure_matrix,
                description="pressure equality",
            ),
            **self._sorption_constraints(),
            **self._vapor_constraints(),
        }

    def _initial_temperature_edges(self):
        return (
            super()._initial_temperature_edges()
            + self._desorber_temperature_edges()
        )

    def propagate_wrapper_to_target(self, branch):
        if self in branch["components"]:
            return
        outconn = self.outl[0]
        branch["connections"] += [outconn]
        branch["components"] += [self]
        outconn.target.propagate_wrapper_to_target(branch)


@component_registry
class HeatedDesorber(_DesorberMixin, TwoSidedSorptionBase):
    r"""
    Desorber (generator) with an integrated heating fluid side.

    Like :class:`Desorber`, the component separates a refrigerant-rich
    solution (:code:`in1`) into a refrigerant-poor solution (:code:`out1`)
    at its bubble point and a refrigerant vapor (:code:`out3`). The heat of
    desorption is supplied by a heating fluid (:code:`in2` to :code:`out2`)
    through the energy balance of both sides.

    **Mandatory Equations**

    - sorption side mass balance :py:meth:`mass_flow_func`
    - heating fluid mass flow and composition equality
    - key component balance :py:meth:`key_balance_func`
    - sorption side pressure equality :py:meth:`pressure_structure_matrix`
    - solution outlet at the bubble point :py:meth:`saturation_func`
    - fluid closure of the solution outlet :py:meth:`fluid_balance_func`
    - energy balance :py:meth:`energy_balance_func`
    - vapor outlet equations as for :class:`Desorber`

    **Optional Equations**

    - heat flow :py:meth:`heat_func`
    - terminal temperature differences :py:meth:`ttd_u_func`,
      :py:meth:`ttd_l_func`
    - heating fluid side pressure ratio or pressure drop
    - approach to the bubble point :code:`td_eq`

    Inlets/Outlets

    - in1: refrigerant-rich solution
    - in2: heating fluid
    - out1: refrigerant-poor solution
    - out2: heating fluid
    - out3: refrigerant vapor

    .. note::

        The heating fluid occupies :code:`in2` and :code:`out2` in analogy to
        the two-sided heat exchangers (matching :code:`pr2` and :code:`dp2`),
        therefore the vapor outlet is :code:`out3`.

    Parameters
    ----------
    label : str
        The label of the component.

    Q : float, dict
        Heat transferred from the heating fluid (hot) side to the sorption
        process :math:`\dot Q`/W, always negative.

    ttd_u : float, dict
        Upper terminal temperature difference
        :math:`T_\text{in,2} - T_\text{out,1}`/K.

    ttd_l : float, dict
        Lower terminal temperature difference
        :math:`T_\text{out,2} - T_\text{in,1}`/K.

    pr2 : float, dict
        Heating fluid side outlet to inlet pressure ratio.

    dp2 : float, dict
        Heating fluid side inlet to outlet pressure drop/Pa.

    td_eq : float, dict
        Temperature difference of the solution outlet below its bubble point
        :math:`\Delta T_\text{eq}`/K.

    vapor_state : str
        State of the vapor leaving the desorber, see :class:`Desorber`.
    """

    _is_wrapper_branch_source = True
    _htf_dh_sign = -1

    def get_parameters(self):
        return {**super().get_parameters(), **self._vapor_state_parameter()}

    @staticmethod
    def inlets():
        return ["in1", "in2"]

    @staticmethod
    def outlets():
        return ["out1", "out2", "out3"]

    def _sorption_inlets(self):
        return [self.inl[0]]

    def _sorption_outlets(self):
        return [self.outl[0], self.outl[2]]

    def _solution_inlet(self):
        return self.inl[0]

    def _solution_outlet(self):
        return self.outl[0]

    def _vapor_outlet(self):
        return self.outl[2]

    def get_mandatory_constraints(self):
        return {
            "mass_flow_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.mass_flow_func,
                dependents=self.mass_flow_dependents,
                description="sorption side mass balance",
            ),
            "pressure_constraints": MandatoryConstraint(
                num_eq_sets=2,
                structure_matrix=self.pressure_structure_matrix,
                description="sorption side pressure equality",
            ),
            **self._sorption_constraints(),
            **self._vapor_constraints(),
            **self._two_sided_constraints(),
        }

    def calc_Q(self):
        return self._htf_heat()

    def ttd_u_func(self):
        r"""
        Equation for the upper terminal temperature difference.

        .. math::

            0 = \Delta T_\text{u} - T_\text{in,2} + T_\text{out,1}
        """
        return self.ttd_u.val_SI - self.inl[1].calc_T() + self.outl[0].calc_T()

    def ttd_u_dependents(self):
        solution = self.outl[0]
        return {
            "scalars": [[
                self.inl[1].p, self.inl[1].h, solution.p, solution.h
            ]],
            "vectors": [self._key_dependents(solution)],
        }

    def ttd_l_func(self):
        r"""
        Equation for the lower terminal temperature difference.

        .. math::

            0 = \Delta T_\text{l} - T_\text{out,2} + T_\text{in,1}
        """
        return self.ttd_l.val_SI - self.outl[1].calc_T() + self.inl[0].calc_T()

    def ttd_l_dependents(self):
        solution = self.inl[0]
        return {
            "scalars": [[
                self.outl[1].p, self.outl[1].h, solution.p, solution.h
            ]],
            "vectors": [self._key_dependents(solution)],
        }

    def _initial_temperature_edges(self):
        return self._ttd_temperature_edges(
            upper=(self.outl[0], self.inl[1]),
            lower=(self.inl[0], self.outl[1]),
        ) + self._desorber_temperature_edges()

    def calc_parameters(self):
        super().calc_parameters()
        self.ttd_u.val_SI = self.inl[1].T.val_SI - self.outl[0].T.val_SI
        self.ttd_l.val_SI = self.outl[1].T.val_SI - self.inl[0].T.val_SI

    def propagate_wrapper_to_target(self, branch):
        inconn = branch["connections"][-1]
        if inconn == self.inl[1]:
            outconn = self.outl[1]
            branch["connections"] += [outconn]
            branch["components"] += [self]
            outconn.target.propagate_wrapper_to_target(branch)
            return
        if self in branch["components"]:
            return
        outconn = self.outl[0]
        branch["connections"] += [outconn]
        branch["components"] += [self]
        outconn.target.propagate_wrapper_to_target(branch)
