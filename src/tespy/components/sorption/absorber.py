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

    Ports
    -----

    - Fluid inlets: in1, in2
    - Fluid outlets: out1

    Mandatory Equations
    -------------------

    - mass balance: :py:meth:`mass_flow_func <tespy.components.sorption.base.SorptionBase.mass_flow_func>`
    - pressure equality: :py:meth:`pressure_structure_matrix <tespy.components.sorption.base.SorptionBase.pressure_structure_matrix>`
    - key component mass balance: :py:meth:`key_balance_func <tespy.components.sorption.base.SorptionBase.key_balance_func>`
    - solution outlet at the bubble point: :py:meth:`saturation_func <tespy.components.sorption.base.SorptionBase.saturation_func>`
    - fluid closure of the solution outlet: :py:meth:`fluid_balance_func <tespy.components.sorption.base.SorptionBase.fluid_balance_func>`

    Parameters
    ----------

    char_warnings : bool
        Ignore warnings on default characteristics usage for this component.

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the components design case.

    label : str
        The label of the component.

    local_design : bool
        Treat this component in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this component in offdesign mode in a design calculation.

    offdesign : list
        List containing offdesign parameters (stated as String).

    printout : bool
        Include this component in the network's results printout.

    Q : float, dict
        Heat transferred into the component. Quantity: :code:`heat`.
        Equation: :py:meth:`heat_func <tespy.components.sorption.base.SorptionBase.heat_func>`.

    td_eq : float, dict
        Temperature difference of the solution outlet to its bubble point.
        Quantity: :code:`temperature_difference`.
        Equation: :py:meth:`saturation_func <tespy.components.sorption.base.SorptionBase.saturation_func>`.

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

    Ports
    -----

    - Fluid inlets: in1, in2, in3
    - Fluid outlets: out1, out2

    Mandatory Equations
    -------------------

    - sorption side mass balance: :py:meth:`mass_flow_func <tespy.components.sorption.base.SorptionBase.mass_flow_func>`
    - sorption side pressure equality: :py:meth:`pressure_structure_matrix <tespy.components.sorption.base.SorptionBase.pressure_structure_matrix>`
    - key component mass balance: :py:meth:`key_balance_func <tespy.components.sorption.base.SorptionBase.key_balance_func>`
    - solution outlet at the bubble point: :py:meth:`saturation_func <tespy.components.sorption.base.SorptionBase.saturation_func>`
    - fluid closure of the solution outlet: :py:meth:`fluid_balance_func <tespy.components.sorption.base.SorptionBase.fluid_balance_func>`
    - heat transfer fluid side mass flow equality: :py:meth:`htf_variable_equality_structure_matrix <tespy.components.sorption.base.TwoSidedSorptionBase.htf_variable_equality_structure_matrix>`
    - heat transfer fluid side fluid composition equality: :py:meth:`htf_variable_equality_structure_matrix <tespy.components.sorption.base.TwoSidedSorptionBase.htf_variable_equality_structure_matrix>`
    - energy balance between both sides: :py:meth:`energy_balance_func <tespy.components.sorption.base.TwoSidedSorptionBase.energy_balance_func>`

    Parameters
    ----------

    char_warnings : bool
        Ignore warnings on default characteristics usage for this component.

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the components design case.

    dp2 : float, dict
        Heat transfer fluid side inlet to outlet absolute pressure change.
        Quantity: :code:`pressure_difference`.
        Equation: :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`.

    label : str
        The label of the component.

    local_design : bool
        Treat this component in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this component in offdesign mode in a design calculation.

    offdesign : list
        List containing offdesign parameters (stated as String).

    pr2 : float, dict
        Heat transfer fluid side outlet to inlet pressure ratio. Quantity:
        :code:`ratio`.
        Equation: :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`.

    printout : bool
        Include this component in the network's results printout.

    Q : float, dict
        Heat transferred from the hot side. Quantity: :code:`heat`.
        Equation: :py:meth:`heat_func <tespy.components.sorption.base.SorptionBase.heat_func>`.

    td_eq : float, dict
        Temperature difference of the solution outlet to its bubble point.
        Quantity: :code:`temperature_difference`.
        Equation: :py:meth:`saturation_func <tespy.components.sorption.base.SorptionBase.saturation_func>`.

    ttd_l : float, dict
        Lower terminal temperature difference. Quantity:
        :code:`temperature_difference`.
        Equation: :py:meth:`ttd_l_func <tespy.components.sorption.absorber.CooledAbsorber.ttd_l_func>`.

    ttd_u : float, dict
        Upper terminal temperature difference. Quantity:
        :code:`temperature_difference`.
        Equation: :py:meth:`ttd_u_func <tespy.components.sorption.absorber.CooledAbsorber.ttd_u_func>`.

    Notes
    -----

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
