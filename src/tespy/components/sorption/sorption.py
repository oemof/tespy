# -*- coding: utf-8

"""Module for Absorber and Desorber sorption components.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/sorption/sorption.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.nodes.base import NodeBase
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.fluid_properties.mixtures import _xi_sat_libr

_LIBR_KEY = "LiBr"


class _SorptionBase(NodeBase):
    r"""
    Shared equations for absorption-cycle sorption components.

    Both :class:`Absorber` and :class:`Desorber` enforce:

    - total mass balance (from :class:`.NodeBase`)
    - LiBr mass balance between solution connections
    - pressure equality across all ports (from :class:`.NodeBase`)
    - saturation condition at the solution outlet

    The optional heat-flow parameter :code:`Q` adds an energy balance
    equation when set.
    """

    def get_parameters(self):
        return {
            "Q": dc_cp(
                func=self.heat_func,
                deriv=self.heat_deriv,
                num_eq=1,
                description="heat transferred into the component",
            )
        }

    # ------------------------------------------------------------------
    # LiBr mass balance: m_sol_in * xi_in = m_sol_out * xi_out
    # ------------------------------------------------------------------

    def libr_balance_func(self):
        sol_in = self._solution_inlet()
        sol_out = self._solution_outlet()
        return (
            sol_in.m.val_SI * sol_in.fluid.val[_LIBR_KEY]
            - sol_out.m.val_SI * sol_out.fluid.val[_LIBR_KEY]
        )

    def libr_balance_dependents(self):
        sol_in = self._solution_inlet()
        sol_out = self._solution_outlet()
        return {
            "scalars": [[sol_in.m, sol_out.m]],
            "vectors": [{
                sol_in.fluid: {_LIBR_KEY} & sol_in.fluid.is_var,
                sol_out.fluid: {_LIBR_KEY} & sol_out.fluid.is_var,
            }],
        }

    # ------------------------------------------------------------------
    # Saturation condition at solution outlet: xi_out = xi_sat(p_out, T_out)
    # ------------------------------------------------------------------

    def saturation_func(self):
        sol_out = self._solution_outlet()
        T = sol_out.calc_T()
        return _xi_sat_libr(sol_out.p.val_SI, T, sol_out.fluid_data) - sol_out.fluid.val[_LIBR_KEY]

    def saturation_dependents(self):
        sol_out = self._solution_outlet()
        return {
            "scalars": [[sol_out.p, sol_out.h]],
            "vectors": [{sol_out.fluid: {_LIBR_KEY} & sol_out.fluid.is_var}],
        }

    # ------------------------------------------------------------------
    # Binary fluid balance at solution outlet: xi_LiBr + x_H2O = 1
    # ------------------------------------------------------------------

    def fluid_balance_func(self):
        sol_out = self._solution_outlet()
        return sol_out.fluid.val["H2O"] + sol_out.fluid.val[_LIBR_KEY] - 1.0

    def fluid_balance_dependents(self):
        sol_out = self._solution_outlet()
        return {
            "scalars": [[]],
            "vectors": [{sol_out.fluid: {"H2O", _LIBR_KEY} & sol_out.fluid.is_var}],
        }

    # ------------------------------------------------------------------
    # Energy balance (optional)
    # ------------------------------------------------------------------

    def heat_func(self):
        return self.Q.val_SI - self.calc_Q()

    def heat_deriv(self, increment_filter, k, dependents=None):
        for c in self.inl:
            if c.m.is_var:
                self.jacobian[k, c.m.J_col] = -c.h.val_SI
            if c.h.is_var:
                self.jacobian[k, c.h.J_col] = -c.m.val_SI
        for c in self.outl:
            if c.m.is_var:
                self.jacobian[k, c.m.J_col] = c.h.val_SI
            if c.h.is_var:
                self.jacobian[k, c.h.J_col] = c.m.val_SI

    def calc_Q(self):
        return (
            sum(c.h.val_SI * c.m.val_SI for c in self.outl)
            - sum(c.h.val_SI * c.m.val_SI for c in self.inl)
        )

    def _update_num_eq(self):
        sol_out = self._solution_outlet()
        num_eq = 0 if len(sol_out.fluid.is_var) == 0 else 1
        self.constraints["fluid_balance_constraints"].num_eq = num_eq

    def calc_parameters(self):
        self.Q.val = self.calc_Q()

    def convergence_check(self):
        from tespy.tools.fluid_properties.functions import h_mix_pT
        sol_out = self._solution_outlet()
        if not sol_out.h.is_var:
            return
        try:
            T = sol_out.calc_T()
        except Exception:
            T = 250
        if T < 274 or T > 499:
            try:
                h = h_mix_pT(
                    sol_out.p.val_SI, 320, sol_out.fluid_data, sol_out.mixing_rule
                )
                sol_out.h.set_reference_val_SI(h)
            except Exception:
                pass

    def propagate_to_target(self, branch):
        return


@component_registry
class Absorber(_SorptionBase):
    r"""
    Ideal single-effect absorber for LiBr-water absorption cycles.

    The absorber merges a water-vapour stream (:code:`in2`) into a poor
    LiBr solution (:code:`in1`) and produces a rich solution
    (:code:`out1`).  The outlet solution is assumed to be at thermodynamic
    equilibrium (saturation condition).

    Ports
    -----
    - :code:`in1` - poor LiBr solution (SolutionConnection)
    - :code:`in2` - water vapour (Connection)
    - :code:`out1` - rich LiBr solution (SolutionConnection)

    Mandatory Equations
    -------------------
    - mass balance: :math:`\dot m_{in,1} + \dot m_{in,2} = \dot m_{out,1}`
    - LiBr balance:
      :math:`\dot m_{in,1} \xi_{in,1} = \dot m_{out,1} \xi_{out,1}`
    - pressure equality: all ports at the same pressure
    - saturation:
      :math:`\xi_{out,1} = \xi_\text{sat}(p_{out,1},\, T_{out,1})`

    Parameters
    ----------
    Q : float
        Heat removed from the absorber (W); negative by convention.
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
            "mass_flow_constraints": dc_cmc(
                num_eq_sets=1,
                func=self.mass_flow_func,
                dependents=self.mass_flow_dependents,
                description="mass balance",
            ),
            "libr_balance_constraints": dc_cmc(
                num_eq_sets=1,
                func=self.libr_balance_func,
                dependents=self.libr_balance_dependents,
                description="LiBr mass balance",
            ),
            "pressure_constraints": dc_cmc(
                num_eq_sets=2,
                structure_matrix=self.pressure_structure_matrix,
                description="pressure equality",
            ),
            "saturation_constraints": dc_cmc(
                num_eq_sets=1,
                func=self.saturation_func,
                dependents=self.saturation_dependents,
                description="saturation at solution outlet",
            ),
            "fluid_balance_constraints": dc_cmc(
                num_eq_sets=1,
                func=self.fluid_balance_func,
                dependents=self.fluid_balance_dependents,
                description="binary fluid balance at solution outlet",
            ),
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
class Desorber(_SorptionBase):

    _is_wrapper_branch_source = True
    r"""
    Ideal single-effect desorber (generator) for LiBr-water absorption cycles.

    The desorber heats a rich LiBr solution (:code:`in1`) and separates it
    into a poor solution (:code:`out1`) and water vapour (:code:`out2`).
    The poor solution outlet is assumed to be at thermodynamic equilibrium
    (saturation condition).

    Ports
    -----
    - :code:`in1` - rich LiBr solution (SolutionConnection)
    - :code:`out1` - poor LiBr solution (SolutionConnection)
    - :code:`out2` - water vapour (Connection)

    Mandatory Equations
    -------------------
    - mass balance: :math:`\dot m_{in,1} = \dot m_{out,1} + \dot m_{out,2}`
    - LiBr balance:
      :math:`\dot m_{in,1} \xi_{in,1} = \dot m_{out,1} \xi_{out,1}`
    - pressure equality: all ports at the same pressure
    - saturation:
      :math:`\xi_{out,1} = \xi_\text{sat}(p_{out,1},\, T_{out,1})`

    Parameters
    ----------
    Q : float
        Heat supplied to the desorber (W); positive by convention.
    """

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

    def get_mandatory_constraints(self):
        return {
            "mass_flow_constraints": dc_cmc(
                num_eq_sets=1,
                func=self.mass_flow_func,
                dependents=self.mass_flow_dependents,
                description="mass balance",
            ),
            "libr_balance_constraints": dc_cmc(
                num_eq_sets=1,
                func=self.libr_balance_func,
                dependents=self.libr_balance_dependents,
                description="LiBr mass balance",
            ),
            "pressure_constraints": dc_cmc(
                num_eq_sets=2,
                structure_matrix=self.pressure_structure_matrix,
                description="pressure equality",
            ),
            "saturation_constraints": dc_cmc(
                num_eq_sets=1,
                func=self.saturation_func,
                dependents=self.saturation_dependents,
                description="saturation at solution outlet",
            ),
            "fluid_balance_constraints": dc_cmc(
                num_eq_sets=1,
                func=self.fluid_balance_func,
                dependents=self.fluid_balance_dependents,
                description="binary fluid balance at solution outlet",
            ),
        }

    def convergence_check(self):
        super().convergence_check()
        vap_out = self.outl[1]
        if not vap_out.h.is_var:
            return
        try:
            from tespy.tools.fluid_properties.functions import h_mix_pQ
            from tespy.tools.fluid_properties.functions import phase_mix_ph
            phase = phase_mix_ph(vap_out.p.val_SI, vap_out.h.val_SI, vap_out.fluid_data)
            if phase != "g":
                h_g = h_mix_pQ(vap_out.p.val_SI, 1.0, vap_out.fluid_data)
                vap_out.h.set_reference_val_SI(h_g)
        except Exception:
            pass

    def start_fluid_wrapper_branch(self):
        vap_conn = self.outl[1]
        vap_branch = {
            "connections": [vap_conn],
            "components": [self],
        }
        vap_conn.target.propagate_wrapper_to_target(vap_branch)

        sol_conn = self.outl[0]
        sol_branch = {
            "connections": [sol_conn],
            "components": [self],
        }
        sol_conn.target.propagate_wrapper_to_target(sol_branch)

        return {vap_conn.label: vap_branch, sol_conn.label: sol_branch}

    def propagate_wrapper_to_target(self, branch):
        if self in branch["components"]:
            return
        outconn = self.outl[0]
        branch["connections"] += [outconn]
        branch["components"] += [self]
        outconn.target.propagate_wrapper_to_target(branch)


class _TwoSidedSorptionBase(_SorptionBase):
    r"""
    Shared equations for sorption components with an integrated heat
    transfer fluid side.

    The sorption side equations of :class:`_SorptionBase` are restricted to
    the sorption ports. On top of those, the heat transfer fluid side adds:

    - mass flow and fluid composition equality between the heat transfer
      fluid inlet and outlet
    - an overall energy balance (the component itself is adiabatic, all heat
      is exchanged between the two sides)

    By convention the heat transfer fluid enters at :code:`in2` and leaves
    at :code:`out2` on both components. Optional parameters are the heat
    flow :code:`Q` (transferred from the hot side, thus always negative in
    analogy to the two-sided heat exchangers), the terminal temperature
    differences :code:`ttd_u` and :code:`ttd_l` and the heat transfer fluid
    side pressure specifications :code:`pr2` and :code:`dp2`.
    """

    def get_parameters(self):
        return {
            "Q": dc_cp(
                max_val=0,
                func=self.heat_func,
                dependents=self.heat_dependents,
                num_eq_sets=1,
                quantity="heat",
                description="heat transferred from the hot side",
            ),
            "ttd_u": dc_cp(
                min_val=0,
                func=self.ttd_u_func,
                dependents=self.ttd_u_dependents,
                num_eq_sets=1,
                quantity="temperature_difference",
                description="upper terminal temperature difference",
            ),
            "ttd_l": dc_cp(
                min_val=0,
                func=self.ttd_l_func,
                dependents=self.ttd_l_dependents,
                num_eq_sets=1,
                quantity="temperature_difference",
                description="lower terminal temperature difference",
            ),
            "pr2": dc_cp(
                min_val=1e-4, max_val=1, num_eq_sets=1,
                structure_matrix=self.pr_structure_matrix,
                func_params={"pr": "pr2", "inconn": 1, "outconn": 1},
                quantity="ratio",
                description="heat transfer fluid side outlet to inlet pressure ratio",
            ),
            "dp2": dc_cp(
                min_val=0, max_val=1e15, num_eq_sets=1,
                structure_matrix=self.dp_structure_matrix,
                func_params={"dp": "dp2", "inconn": 1, "outconn": 1},
                quantity="pressure_difference",
                description="heat transfer fluid side inlet to outlet absolute pressure change",
            ),
        }

    def get_mandatory_constraints(self):
        return {
            "mass_flow_constraints": dc_cmc(
                num_eq_sets=1,
                func=self.mass_flow_func,
                dependents=self.mass_flow_dependents,
                description="sorption side mass balance",
            ),
            "htf_mass_flow_constraints": dc_cmc(
                num_eq_sets=1,
                structure_matrix=self.htf_variable_equality_structure_matrix,
                func_params={"variable": "m"},
                description="heat transfer fluid side mass flow equality",
            ),
            "htf_fluid_constraints": dc_cmc(
                num_eq_sets=1,
                structure_matrix=self.htf_variable_equality_structure_matrix,
                func_params={"variable": "fluid"},
                description="heat transfer fluid side fluid composition equality",
            ),
            "libr_balance_constraints": dc_cmc(
                num_eq_sets=1,
                func=self.libr_balance_func,
                dependents=self.libr_balance_dependents,
                description="LiBr mass balance",
            ),
            "pressure_constraints": dc_cmc(
                num_eq_sets=2,
                structure_matrix=self.pressure_structure_matrix,
                description="sorption side pressure equality",
            ),
            "saturation_constraints": dc_cmc(
                num_eq_sets=1,
                func=self.saturation_func,
                dependents=self.saturation_dependents,
                description="saturation at solution outlet",
            ),
            "fluid_balance_constraints": dc_cmc(
                num_eq_sets=1,
                func=self.fluid_balance_func,
                dependents=self.fluid_balance_dependents,
                description="binary fluid balance at solution outlet",
            ),
            "energy_balance_constraints": dc_cmc(
                num_eq_sets=1,
                func=self.energy_balance_func,
                dependents=self.energy_balance_dependents,
                description="energy balance between both sides",
            ),
        }

    def mass_flow_func(self):
        res = 0
        for c in self._sorption_inlets():
            res += c.m.val_SI
        for c in self._sorption_outlets():
            res -= c.m.val_SI
        return res

    def mass_flow_dependents(self):
        return [c.m for c in self._sorption_inlets() + self._sorption_outlets()]

    def htf_variable_equality_structure_matrix(self, k, **kwargs):
        variable = kwargs.get("variable")
        self._structure_matrix[k, self.inl[1].get_attr(variable).sm_col] = 1
        self._structure_matrix[k, self.outl[1].get_attr(variable).sm_col] = -1

    def pressure_structure_matrix(self, k):
        conns = self._sorption_inlets() + self._sorption_outlets()
        first = conns[0]
        for eq, conn in enumerate(conns[1:]):
            self._structure_matrix[k + eq, first.p.sm_col] = 1
            self._structure_matrix[k + eq, conn.p.sm_col] = -1

    def energy_balance_func(self):
        return (
            sum(c.h.val_SI * c.m.val_SI for c in self.outl)
            - sum(c.h.val_SI * c.m.val_SI for c in self.inl)
        )

    def energy_balance_dependents(self):
        return [v for c in self.inl + self.outl for v in (c.m, c.h)]

    def heat_dependents(self):
        return [self.inl[1].m, self.inl[1].h, self.outl[1].h]

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
        sorption_ports = self._sorption_inlets() + self._sorption_outlets()
        first = sorption_ports[0]
        edges = []
        for c in sorption_ports[1:]:
            edges += [
                (first.p, c.p, 1.0, 0.0),
                (first.h, c.h, 1.0, 0.0),
            ]
        edges += [
            (self.inl[1].p, self.outl[1].p, 1.0, 0.0),
            (self.inl[1].h, self.outl[1].h, 1.0, 0.0),
        ]
        return edges

    def calc_parameters(self):
        self.Q.val_SI = self.calc_Q()
        self.pr2.val_SI = self.outl[1].p.val_SI / self.inl[1].p.val_SI
        self.dp2.val_SI = self.inl[1].p.val_SI - self.outl[1].p.val_SI


@component_registry
class CooledAbsorber(_TwoSidedSorptionBase):
    r"""
    Absorber with integrated cooling fluid side for LiBr-water absorption
    cycles.

    Like :class:`Absorber`, the component merges a water-vapour stream
    (:code:`in3`) into a poor LiBr solution (:code:`in1`) and produces a
    saturated rich solution (:code:`out1`). The heat of absorption is
    transferred to a cooling fluid (:code:`in2` to :code:`out2`) via a
    mandatory energy balance, so no external heat flow specification is
    required.

    Ports
    -----
    - :code:`in1` - poor LiBr solution (SolutionConnection)
    - :code:`in2` - cooling fluid (Connection)
    - :code:`in3` - water vapour (Connection)
    - :code:`out1` - rich LiBr solution (SolutionConnection)
    - :code:`out2` - cooling fluid (Connection)

    .. note::

        The port numbering deviates from :class:`Absorber`: the cooling
        fluid occupies :code:`in2` and :code:`out2` in analogy to the cold
        side of the two-sided heat exchangers (matching :code:`pr2` and
        :code:`dp2`), therefore the water vapour inlet moves to
        :code:`in3`.

    Mandatory Equations
    -------------------
    - sorption side mass balance:
      :math:`\dot m_{in,1} + \dot m_{in,3} = \dot m_{out,1}`
    - cooling fluid mass flow and fluid composition equality
    - LiBr balance:
      :math:`\dot m_{in,1} \xi_{in,1} = \dot m_{out,1} \xi_{out,1}`
    - pressure equality of the sorption side ports
    - saturation:
      :math:`\xi_{out,1} = \xi_\text{sat}(p_{out,1},\, T_{out,1})`
    - energy balance:
      :math:`0 = \sum \dot m_{out} h_{out} - \sum \dot m_{in} h_{in}`

    Parameters
    ----------
    Q : float
        Heat transferred from the sorption (hot) side to the cooling fluid
        (W); always negative.

    ttd_u : float
        Upper terminal temperature difference
        :math:`T_{in,1} - T_{out,2}` (K).

    ttd_l : float
        Lower terminal temperature difference
        :math:`T_{out,1} - T_{in,2}` (K).

    pr2 : float
        Cooling fluid side outlet to inlet pressure ratio.

    dp2 : float
        Cooling fluid side inlet to outlet absolute pressure change.
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

    def heat_func(self):
        i = self.inl[1]
        o = self.outl[1]
        return self.Q.val_SI + i.m.val_SI * (o.h.val_SI - i.h.val_SI)

    def calc_Q(self):
        i = self.inl[1]
        o = self.outl[1]
        return -i.m.val_SI * (o.h.val_SI - i.h.val_SI)

    def ttd_u_func(self):
        T_sol_in = self.inl[0].calc_T()
        T_cool_out = self.outl[1].calc_T()
        return self.ttd_u.val_SI - T_sol_in + T_cool_out

    def ttd_u_dependents(self):
        sol_in = self.inl[0]
        return {
            "scalars": [[
                sol_in.p, sol_in.h, self.outl[1].p, self.outl[1].h
            ]],
            "vectors": [{sol_in.fluid: {_LIBR_KEY} & sol_in.fluid.is_var}],
        }

    def ttd_l_func(self):
        T_sol_out = self.outl[0].calc_T()
        T_cool_in = self.inl[1].calc_T()
        return self.ttd_l.val_SI - T_sol_out + T_cool_in

    def ttd_l_dependents(self):
        sol_out = self.outl[0]
        return {
            "scalars": [[
                sol_out.p, sol_out.h, self.inl[1].p, self.inl[1].h
            ]],
            "vectors": [{sol_out.fluid: {_LIBR_KEY} & sol_out.fluid.is_var}],
        }

    def _initial_temperature_edges(self):
        sorption_ports = self._sorption_inlets() + self._sorption_outlets()
        first = sorption_ports[0]
        edges = [(first, c, 0.0, 1.0) for c in sorption_ports[1:]]
        edges += [(self.inl[1], self.outl[1], 0.0, 1.0)]
        ttd_upper = self.ttd_u.val_SI if self.ttd_u.is_set else 10.0
        weight_upper = 5.0 if self.ttd_u.is_set else 0.3
        ttd_lower = self.ttd_l.val_SI if self.ttd_l.is_set else 10.0
        weight_lower = 5.0 if self.ttd_l.is_set else 0.3
        edges += [
            (self.outl[1], self.inl[0], ttd_upper, weight_upper),
            (self.inl[1], self.outl[0], ttd_lower, weight_lower),
        ]
        return edges

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


@component_registry
class HeatedDesorber(_TwoSidedSorptionBase):
    r"""
    Desorber (generator) with integrated heating fluid side for LiBr-water
    absorption cycles.

    Like :class:`Desorber`, the component separates a rich LiBr solution
    (:code:`in1`) into a saturated poor solution (:code:`out1`) and water
    vapour (:code:`out3`). The heat of desorption is supplied by a heating
    fluid (:code:`in2` to :code:`out2`) via a mandatory energy balance, so
    no external heat flow specification is required.

    Ports
    -----
    - :code:`in1` - rich LiBr solution (SolutionConnection)
    - :code:`in2` - heating fluid (Connection)
    - :code:`out1` - poor LiBr solution (SolutionConnection)
    - :code:`out2` - heating fluid (Connection)
    - :code:`out3` - water vapour (Connection)

    .. note::

        The port numbering deviates from :class:`Desorber`: the heating
        fluid occupies :code:`in2` and :code:`out2` in analogy to the
        two-sided heat exchangers (matching :code:`pr2` and :code:`dp2`),
        therefore the water vapour outlet moves to :code:`out3`.

    Mandatory Equations
    -------------------
    - sorption side mass balance:
      :math:`\dot m_{in,1} = \dot m_{out,1} + \dot m_{out,3}`
    - heating fluid mass flow and fluid composition equality
    - LiBr balance:
      :math:`\dot m_{in,1} \xi_{in,1} = \dot m_{out,1} \xi_{out,1}`
    - pressure equality of the sorption side ports
    - saturation:
      :math:`\xi_{out,1} = \xi_\text{sat}(p_{out,1},\, T_{out,1})`
    - energy balance:
      :math:`0 = \sum \dot m_{out} h_{out} - \sum \dot m_{in} h_{in}`

    Parameters
    ----------
    Q : float
        Heat transferred from the heating fluid (hot) side to the sorption
        process (W); always negative.

    ttd_u : float
        Upper terminal temperature difference
        :math:`T_{in,2} - T_{out,1}` (K).

    ttd_l : float
        Lower terminal temperature difference
        :math:`T_{out,2} - T_{in,1}` (K).

    pr2 : float
        Heating fluid side outlet to inlet pressure ratio.

    dp2 : float
        Heating fluid side inlet to outlet absolute pressure change.
    """

    _is_wrapper_branch_source = True
    _htf_dh_sign = -1

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

    def heat_func(self):
        i = self.inl[1]
        o = self.outl[1]
        return self.Q.val_SI - i.m.val_SI * (o.h.val_SI - i.h.val_SI)

    def calc_Q(self):
        i = self.inl[1]
        o = self.outl[1]
        return i.m.val_SI * (o.h.val_SI - i.h.val_SI)

    def ttd_u_func(self):
        T_heat_in = self.inl[1].calc_T()
        T_sol_out = self.outl[0].calc_T()
        return self.ttd_u.val_SI - T_heat_in + T_sol_out

    def ttd_u_dependents(self):
        sol_out = self.outl[0]
        return {
            "scalars": [[
                self.inl[1].p, self.inl[1].h, sol_out.p, sol_out.h
            ]],
            "vectors": [{sol_out.fluid: {_LIBR_KEY} & sol_out.fluid.is_var}],
        }

    def ttd_l_func(self):
        T_heat_out = self.outl[1].calc_T()
        T_sol_in = self.inl[0].calc_T()
        return self.ttd_l.val_SI - T_heat_out + T_sol_in

    def ttd_l_dependents(self):
        sol_in = self.inl[0]
        return {
            "scalars": [[
                self.outl[1].p, self.outl[1].h, sol_in.p, sol_in.h
            ]],
            "vectors": [{sol_in.fluid: {_LIBR_KEY} & sol_in.fluid.is_var}],
        }

    def _initial_temperature_edges(self):
        sorption_ports = self._sorption_inlets() + self._sorption_outlets()
        first = sorption_ports[0]
        edges = [(first, c, 0.0, 1.0) for c in sorption_ports[1:]]
        edges += [(self.inl[1], self.outl[1], 0.0, 1.0)]
        ttd_upper = self.ttd_u.val_SI if self.ttd_u.is_set else 10.0
        weight_upper = 5.0 if self.ttd_u.is_set else 0.3
        ttd_lower = self.ttd_l.val_SI if self.ttd_l.is_set else 10.0
        weight_lower = 5.0 if self.ttd_l.is_set else 0.3
        edges += [
            (self.outl[0], self.inl[1], ttd_upper, weight_upper),
            (self.inl[0], self.outl[1], ttd_lower, weight_lower),
        ]
        return edges

    def calc_parameters(self):
        super().calc_parameters()
        self.ttd_u.val_SI = self.inl[1].T.val_SI - self.outl[0].T.val_SI
        self.ttd_l.val_SI = self.outl[1].T.val_SI - self.inl[0].T.val_SI

    def convergence_check(self):
        super().convergence_check()
        vap_out = self.outl[2]
        if not vap_out.h.is_var:
            return
        try:
            from tespy.tools.fluid_properties.functions import h_mix_pQ
            from tespy.tools.fluid_properties.functions import phase_mix_ph
            phase = phase_mix_ph(vap_out.p.val_SI, vap_out.h.val_SI, vap_out.fluid_data)
            if phase != "g":
                h_g = h_mix_pQ(vap_out.p.val_SI, 1.0, vap_out.fluid_data)
                vap_out.h.set_reference_val_SI(h_g)
        except Exception:
            pass

    def start_fluid_wrapper_branch(self):
        vap_conn = self.outl[2]
        vap_branch = {
            "connections": [vap_conn],
            "components": [self],
        }
        vap_conn.target.propagate_wrapper_to_target(vap_branch)

        sol_conn = self.outl[0]
        sol_branch = {
            "connections": [sol_conn],
            "components": [self],
        }
        sol_conn.target.propagate_wrapper_to_target(sol_branch)

        return {vap_conn.label: vap_branch, sol_conn.label: sol_branch}

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
