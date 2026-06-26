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
