# -*- coding: utf-8

"""Module of the absorber components.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/sorption/absorber.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.sorption.base import SorptionBase
from tespy.components.sorption.base import TwoSidedSorptionBase
from tespy.tools.data_containers import MandatoryConstraint


@component_registry
class Absorber(SorptionBase):
    r"""
    Ideal absorber of absorption cycles.

    The absorber merges a refrigerant vapor (:code:`in2`) into a
    refrigerant-poor solution (:code:`in1`) and produces a refrigerant-rich
    solution (:code:`out1`) at its bubble point. The working pair is the
    mixing rule of the solution connections, see
    :class:`.SolutionConnection`. For LiBr-water the vapor is pure water on a
    plain :class:`.Connection`, for NH3-water it is a
    :class:`.SolutionConnection` of the same pair.

    **Mandatory Equations**

    - mass balance :py:meth:`mass_flow_func`
    - key component balance :py:meth:`key_balance_func`
    - pressure equality :py:meth:`pressure_structure_matrix`
    - solution outlet at the bubble point :py:meth:`saturation_func`
    - fluid closure of the solution outlet :py:meth:`fluid_balance_func`

    **Optional Equations**

    - heat flow :py:meth:`heat_func`
    - approach to the bubble point :code:`td_eq` replacing the bubble point
      equation, :py:meth:`saturation_func`

    Inlets/Outlets

    - in1: refrigerant-poor solution
    - in2: refrigerant vapor
    - out1: refrigerant-rich solution

    Parameters
    ----------
    label : str
        The label of the component.

    Q : float, dict
        Heat transferred into the absorber :math:`\dot Q`/W, negative for
        the heat of absorption that is removed.

    td_eq : float, dict
        Temperature difference of the solution outlet below its bubble point
        :math:`\Delta T_\text{eq}`/K.

    Example
    -------
    Water vapor at 5 degC saturation pressure is absorbed in a LiBr solution
    of 62 % LiBr, the refrigerant-rich solution leaves at 34 degC.

    >>> from CoolProp.CoolProp import PropsSI
    >>> from tespy.components import Absorber, Sink, Source
    >>> from tespy.connections import Connection, SolutionConnection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> absorber = Absorber("absorber")
    >>> poor = SolutionConnection(
    ...     Source("poor"), "out1", absorber, "in1", mixing_rule="libr-water"
    ... )
    >>> vapor = Connection(Source("vapor"), "out1", absorber, "in2")
    >>> rich = SolutionConnection(absorber, "out1", Sink("rich"), "in1")
    >>> nw.add_conns(poor, vapor, rich)
    >>> p_ev = PropsSI("P", "T", 278.15, "Q", 1, "water")
    >>> poor.set_attr(xi=0.62, T=320)
    >>> vapor.set_attr(fluid={"H2O": 1}, p=p_ev, m=1, x=1)
    >>> rich.set_attr(T=307.15)
    >>> nw.solve("design")
    >>> round(rich.xi.val, 3)
    0.548
    >>> round(absorber.Q.val / 1e6, 2)
    -2.9
    """

    @staticmethod
    def inlets():
        return ["in1", "in2"]

    @staticmethod
    def outlets():
        return ["out1"]

    def _solution_inlet(self):
        return self.inl[0]

    def _solution_outlet(self):
        return self.outl[0]

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
        }

    def propagate_wrapper_to_target(self, branch):
        inconn = branch["connections"][-1]
        if inconn == self.inl[1]:
            branch["components"] += [self]
            return
        if self in branch["components"]:
            return
        outconn = self.outl[0]
        branch["connections"] += [outconn]
        branch["components"] += [self]
        outconn.target.propagate_wrapper_to_target(branch)


@component_registry
class CooledAbsorber(TwoSidedSorptionBase):
    r"""
    Absorber with an integrated cooling fluid side.

    Like :class:`Absorber`, the component merges a refrigerant vapor
    (:code:`in3`) into a refrigerant-poor solution (:code:`in1`) and
    produces a refrigerant-rich solution (:code:`out1`) at its bubble point.
    The heat of absorption is transferred to a cooling fluid (:code:`in2` to
    :code:`out2`) through the energy balance of both sides.

    **Mandatory Equations**

    - sorption side mass balance :py:meth:`mass_flow_func`
    - cooling fluid mass flow and composition equality
    - key component balance :py:meth:`key_balance_func`
    - sorption side pressure equality :py:meth:`pressure_structure_matrix`
    - solution outlet at the bubble point :py:meth:`saturation_func`
    - fluid closure of the solution outlet :py:meth:`fluid_balance_func`
    - energy balance :py:meth:`energy_balance_func`

    **Optional Equations**

    - heat flow :py:meth:`heat_func`
    - terminal temperature differences :py:meth:`ttd_u_func`,
      :py:meth:`ttd_l_func`
    - cooling fluid side pressure ratio or pressure drop
    - approach to the bubble point :code:`td_eq`

    Inlets/Outlets

    - in1: refrigerant-poor solution
    - in2: cooling fluid
    - in3: refrigerant vapor
    - out1: refrigerant-rich solution
    - out2: cooling fluid

    .. note::

        The cooling fluid occupies :code:`in2` and :code:`out2` in analogy to
        the two-sided heat exchangers (matching :code:`pr2` and :code:`dp2`),
        therefore the vapor inlet is :code:`in3`.

    Parameters
    ----------
    label : str
        The label of the component.

    Q : float, dict
        Heat transferred from the sorption (hot) side to the cooling fluid
        :math:`\dot Q`/W, always negative.

    ttd_u : float, dict
        Upper terminal temperature difference
        :math:`T_\text{in,1} - T_\text{out,2}`/K.

    ttd_l : float, dict
        Lower terminal temperature difference
        :math:`T_\text{out,1} - T_\text{in,2}`/K.

    pr2 : float, dict
        Cooling fluid side outlet to inlet pressure ratio.

    dp2 : float, dict
        Cooling fluid side inlet to outlet pressure drop/Pa.

    td_eq : float, dict
        Temperature difference of the solution outlet below its bubble point
        :math:`\Delta T_\text{eq}`/K.
    """

    _htf_dh_sign = 1

    @staticmethod
    def inlets():
        return ["in1", "in2", "in3"]

    @staticmethod
    def outlets():
        return ["out1", "out2"]

    def _sorption_inlets(self):
        return [self.inl[0], self.inl[2]]

    def _sorption_outlets(self):
        return [self.outl[0]]

    def _solution_inlet(self):
        return self.inl[0]

    def _solution_outlet(self):
        return self.outl[0]

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
            **self._two_sided_constraints(),
        }

    def calc_Q(self):
        return -self._htf_heat()

    def ttd_u_func(self):
        r"""
        Equation for the upper terminal temperature difference.

        .. math::

            0 = \Delta T_\text{u} - T_\text{in,1} + T_\text{out,2}
        """
        return self.ttd_u.val_SI - self.inl[0].calc_T() + self.outl[1].calc_T()

    def ttd_u_dependents(self):
        solution = self.inl[0]
        return {
            "scalars": [[
                solution.p, solution.h, self.outl[1].p, self.outl[1].h
            ]],
            "vectors": [self._key_dependents(solution)],
        }

    def ttd_l_func(self):
        r"""
        Equation for the lower terminal temperature difference.

        .. math::

            0 = \Delta T_\text{l} - T_\text{out,1} + T_\text{in,2}
        """
        return self.ttd_l.val_SI - self.outl[0].calc_T() + self.inl[1].calc_T()

    def ttd_l_dependents(self):
        solution = self.outl[0]
        return {
            "scalars": [[
                solution.p, solution.h, self.inl[1].p, self.inl[1].h
            ]],
            "vectors": [self._key_dependents(solution)],
        }

    def _initial_temperature_edges(self):
        return self._ttd_temperature_edges(
            upper=(self.outl[1], self.inl[0]),
            lower=(self.inl[1], self.outl[0]),
        )

    def calc_parameters(self):
        super().calc_parameters()
        self.ttd_u.val_SI = self.inl[0].T.val_SI - self.outl[1].T.val_SI
        self.ttd_l.val_SI = self.outl[0].T.val_SI - self.inl[1].T.val_SI

    def propagate_wrapper_to_target(self, branch):
        inconn = branch["connections"][-1]
        if inconn == self.inl[1]:
            outconn = self.outl[1]
            branch["connections"] += [outconn]
            branch["components"] += [self]
            outconn.target.propagate_wrapper_to_target(branch)
            return
        if inconn == self.inl[2]:
            branch["components"] += [self]
            return
        if self in branch["components"]:
            return
        outconn = self.outl[0]
        branch["connections"] += [outconn]
        branch["components"] += [self]
        outconn.target.propagate_wrapper_to_target(branch)
