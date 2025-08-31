# -*- coding: utf-8

"""Module of class Merge.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/nodes/merge.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import component_registry
from tespy.components.nodes.base import NodeBase
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import s_mix_pT


@component_registry
class Merge(NodeBase):
    r"""
    Class for merge points with multiple inflows and one outflow.

    **Mandatory Equations**

    - :py:meth:`tespy.components.nodes.base.NodeBase.mass_flow_func`
    - :py:meth:`tespy.components.nodes.base.NodeBase.pressure_structure_matrix`
    - :py:meth:`tespy.components.nodes.merge.Merge.fluid_func`
    - :py:meth:`tespy.components.nodes.merge.Merge.energy_balance_func`

    Inlets/Outlets

    - specify number of inlets with :code:`num_in` (default value: 2)
    - out1

    Image

    .. image:: /api/_images/Merge.svg
       :alt: flowsheet of the merge
       :align: center
       :class: only-light

    .. image:: /api/_images/Merge_darkmode.svg
       :alt: flowsheet of the merge
       :align: center
       :class: only-dark

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    design_path : str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings : boolean
        Ignore warnings on default characteristics usage for this component.

    printout : boolean
        Include this component in the network's results printout.

    num_in : float, dict
        Number of inlets for this component, default value: 2.

    Example
    -------
    The merge mixes a specified number of mass flows and has a single outlet.
    At the outlet, fluid composition and enthalpy are calculated by mass
    weighted fluid composition and enthalpy of the inlets.

    >>> from tespy.components import Sink, Source, Merge
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar"
    ... })
    >>> so1 = Source('source1')
    >>> so2 = Source('source2')
    >>> so3 = Source('source3')
    >>> si1 = Sink('sink')
    >>> m = Merge('merge', num_in=3)
    >>> inc1 = Connection(so1, 'out1', m, 'in1')
    >>> inc2 = Connection(so2, 'out1', m, 'in2')
    >>> inc3 = Connection(so3, 'out1', m, 'in3')
    >>> outg = Connection(m, 'out1', si1, 'in1')
    >>> nw.add_conns(inc1, inc2, inc3, outg)

    Consider a merge with three inlets which mixes three mass flows of the
    same fluid. In this case, the outlet mass flow will be the sum of both
    inlet mass flows and the outlet enthalpy will be the weighted sum of the
    inlet enthalpies. The pressure is equal for all connections of the merge
    imposed by the component's mandatory constraints.

    >>> T = 293.15
    >>> inc1.set_attr(fluid={'O2': 1}, p=1, T=300, m=5)
    >>> inc2.set_attr(fluid={'O2': 1}, T=450, m=5)
    >>> inc3.set_attr(fluid={'O2': 1}, T=350, m=5)
    >>> nw.solve('design')
    >>> round(outg.m.val_SI, 1)
    15.0
    >>> round(outg.h.val_SI, 0)
    334919.0
    >>> round(outg.T.val_SI, 0)
    367.0

    We could also fix the outlet temperature and by that determine a missing
    mass flow, e.g. the hottest incoming stream.

    >>> outg.set_attr(T=360)
    >>> inc2.set_attr(m=None)
    >>> nw.solve("design")
    >>> round(inc2.m.val_SI, 1)
    3.8

    More interesting things can happen, if we want to take the fluid
    composition into account. For example, air (O2 + N2) is mixed with pure
    nitrogen and pure oxygen flows. At the outlet we want to have a new mixture
    with a fixed amount of nitrogen, e.g. 40 %. All gases enter the component
    at the same temperature. When changing the fluids, we have to rerun the
    network fluid detection, which is part of the topological setup. This
    usually only happens if you run a network with changed topology.

    >>> T = 293.15
    >>> inc1.reset_fluid_vector()
    >>> inc2.reset_fluid_vector()
    >>> inc3.reset_fluid_vector()
    >>> outg.reset_fluid_vector()
    >>> inc1.set_attr(fluid={'O2': 0.23, 'N2': 0.77}, p=1, T=T, m=5)
    >>> inc2.set_attr(fluid={'O2': 1}, T=T, m=5)
    >>> inc3.set_attr(fluid={'N2': 1}, T=T, m=None)
    >>> outg.set_attr(fluid={'N2': 0.4}, T=None)
    >>> nw.solve('design')
    >>> m_expected = (
    ...     (inc1.fluid.val["O2"] * inc1.m.val_SI + inc2.m.val_SI)
    ...     / (1 - outg.fluid.val["N2"])
    ... )
    >>> round(outg.m.val_SI, 2) == round(m_expected, 2)
    True
    >>> abs((outg.T.val_SI - T) / T) < 0.01
    True

    >>> T = 173.15
    >>> inc1.set_attr(T=T)
    >>> inc2.set_attr(T=T)
    >>> inc3.set_attr(T=T)
    >>> nw.solve('design')
    >>> abs((outg.T.val_SI - T) / T) < 0.01
    True
    """

    @staticmethod
    def get_parameters():
        return {'num_in': dc_simple()}

    def _update_num_eq(self):
        self.variable_fluids = set(
            [fluid for c in self.inl + self.outl for fluid in c.fluid.is_var]
        )
        set_fluids = set(
            [fluid for c in self.inl + self.outl for fluid in c.fluid.is_set]
        )
        self.all_fluids = self.variable_fluids | set_fluids
        if len(self.variable_fluids) == 0 and len(set_fluids) == 0:
            fluid_eq = 0
            self.constraints["mass_flow_constraints"].num_eq = 1
        elif len(self.variable_fluids) == 0:
            fluid_eq = len(self.all_fluids)
            self.constraints["mass_flow_constraints"].num_eq = 0
        else:
            fluid_eq = len(self.variable_fluids)
        self.constraints["fluid_constraints"].num_eq = fluid_eq

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': dc_cmc(**{
                'num_eq_sets': 1,
                'func': self.mass_flow_func,
                'dependents': self.mass_flow_dependents,
            }),
            'fluid_constraints': dc_cmc(**{
                'num_eq_sets': 1,
                'func': self.fluid_func,
                'deriv': self.fluid_deriv,
                'dependents': self.fluid_dependents
            }),
            'energy_balance_constraints': dc_cmc(**{
                'num_eq_sets': 1,
                'func': self.energy_balance_func,
                'dependents': self.energy_balance_dependents,
            }),
            'pressure_constraints': dc_cmc(**{
                'structure_matrix': self.pressure_structure_matrix,
                'num_eq_sets': self.num_i + self.num_o - 1
            })
        }

    def inlets(self):
        if self.num_in.is_set:
            return [f'in{i + 1}' for i in range(self.num_in.val)]
        else:
            self.set_attr(num_in=2)
            return self.inlets()

    @staticmethod
    def outlets():
        return ['out1']

    def fluid_func(self):
        r"""
        Calculate the vector of residual values for fluid balance equations.

        Returns
        -------
        residual : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = \sum_i \dot{m}_{in,i} \cdot x_{fl,in,i} -
                \dot {m}_{out} \cdot x_{fl,out}\\
                \forall fl \in \text{network fluids},
                \; \forall i \in \text{inlets}
        """
        residual = []
        # we take the total mass flow to handle more than one outlet if necessary
        total_mass_flow = sum([c.m.val_SI for c in self.outl])
        for fluid in self.all_fluids:
            res = -self.outl[0].fluid.val.get(fluid, 0) * total_mass_flow
            for i in self.inl:
                res += i.fluid.val.get(fluid, 0) * i.m.val_SI
            residual += [res]
        return residual

    def fluid_deriv(self, increment_filter, k, dependents=None):
        r"""
        Calculate partial derivatives of fluid balance.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        # we take the total mass flow to handle more than one outlet if necessary
        total_mass_flow = sum([c.m.val_SI for c in self.outl])
        for fluid in self.all_fluids:
            for i in self.inl:
                if i.m.is_var:
                    self.jacobian[k, i.m.J_col] = i.fluid.val.get(fluid, 0)
                if fluid in i.fluid.is_var:
                    self.jacobian[k, i.fluid.J_col[fluid]] = i.m.val_SI
            for o in self.outl:
                if o.m.is_var:
                    self.jacobian[k, o.m.J_col] = -self.outl[0].fluid.val.get(fluid, 0)
            if fluid in self.outl[0].fluid.is_var:
                self.jacobian[k, self.outl[0].fluid.J_col[fluid]] = -total_mass_flow
            k += 1

    def fluid_dependents(self):
        return {
            "scalars": [
                [c.m for c in self.inl + self.outl] for f in self.all_fluids
            ],
            "vectors": [{
                # only depends on first outlet (there is only one in merge)
                # but there may be more in inheriting components
                c.fluid: c.fluid.is_var for c in self.inl + self.outl[:1]
            } for f in self.all_fluids]
        }

    def energy_balance_func(self):
        r"""
        Calculate energy balance.

        Returns
        -------
        residual : float
            Residual value of energy balance.

            .. math::

                0 = \sum_i \left(\dot{m}_{in,i} \cdot h_{in,i} \right) -
                \dot{m}_{out} \cdot h_{out}\\
                \forall i \in \text{inlets}
        """
        # we take the total mass flow to handle more than one outlet if necessary
        total_mass_flow = sum([c.m.val_SI for c in self.outl])
        res = -total_mass_flow * self.outl[0].h.val_SI
        for i in self.inl:
            res += i.m.val_SI * i.h.val_SI
        return res

    def energy_balance_dependents(self):
        dependents = []
        for c in self.inl + self.outl:
            dependents += [c.m, c.h]
        return dependents

    def propagate_wrapper_to_target(self, branch):
        if self in branch["components"]:
            return

        branch["components"] += [self]
        for outconn in self.outl:
            branch["connections"] += [outconn]
            outconn.target.propagate_wrapper_to_target(branch)

    def entropy_balance(self):
        r"""
        Calculate entropy balance of a merge.

        Note
        ----
        A definition of reference points is included for compensation of
        differences in zero point definitions of different fluid compositions.

        - Reference temperature: 298.15 K.
        - Reference pressure: 1 bar.

        .. math::

            \dot{S}_\mathrm{irr}= \dot{m}_\mathrm{out} \cdot
            \left( s_\mathrm{out} - s_\mathrm{out,ref} \right)
            - \sum_{i} \dot{m}_{\mathrm{in,}i} \cdot
            \left( s_{\mathrm{in,}i} - s_{\mathrm{in,ref,}i} \right)\\
        """
        T_ref = 298.15
        p_ref = 1e5
        o = self.outl[0]
        self.S_irr = o.m.val_SI * (
            o.s.val_SI - s_mix_pT(p_ref, T_ref, o.fluid_data, o.mixing_rule)
        )
        for i in self.inl:
            self.S_irr -= i.m.val_SI * (
                i.s.val_SI - s_mix_pT(p_ref, T_ref, i.fluid_data, i.mixing_rule)
            )

    def exergy_balance(self, T0):
        r"""
        Calculate exergy balance of a merge.

        Parameters
        ----------
        T0 : float
            Ambient temperature T0 / K.

        Note
        ----
        Please note, that the exergy balance accounts for physical exergy only.

        .. math ::

            \dot{E}_\mathrm{P} =
            \begin{cases}
            \begin{cases}
            \sum_i \dot{m}_i \cdot \left(e_\mathrm{out}^\mathrm{PH} -
            e_{\mathrm{in,}i}^\mathrm{PH}\right)
            & T_{\mathrm{in,}i} < T_\mathrm{out} \text{ \& }
            T_{\mathrm{in,}i} \geq T_0 \\
            \sum_i \dot{m}_i \cdot e_\mathrm{out}^\mathrm{PH}
            & T_{\mathrm{in,}i} < T_\mathrm{out} \text{ \& }
            T_{\mathrm{in,}i} < T_0 \\
            \end{cases} & T_\mathrm{out} > T_0\\

            \text{not defined (nan)} & T_\mathrm{out} = T_0\\

            \begin{cases}
            \sum_i \dot{m}_i \cdot e_\mathrm{out}^\mathrm{PH}
            & T_{\mathrm{in,}i} > T_\mathrm{out} \text{ \& }
            T_{\mathrm{in,}i} \geq T_0 \\
            \sum_i \dot{m}_i \cdot \left(e_\mathrm{out}^\mathrm{PH} -
            e_{\mathrm{in,}i}^\mathrm{PH}\right)
            & T_{\mathrm{in,}i} > T_\mathrm{out} \text{ \& }
            T_{\mathrm{in,}i} < T_0 \\
            \end{cases} & T_\mathrm{out} < T_0\\
            \end{cases}

            \dot{E}_\mathrm{F} =
            \begin{cases}
            \begin{cases}
            \sum_i \dot{m}_i \cdot \left(e_{\mathrm{in,}i}^\mathrm{PH} -
            e_\mathrm{out}^\mathrm{PH}\right)
            & T_{\mathrm{in,}i} > T_\mathrm{out} \\
            \sum_i \dot{E}_{\mathrm{in,}i}^\mathrm{PH}
            & T_{\mathrm{in,}i} < T_\mathrm{out} \text{ \& }
            T_{\mathrm{in,}i} < T_0 \\
            \end{cases} & T_\mathrm{out} > T_0\\

            \sum_i \dot{E}_{\mathrm{in,}i}^\mathrm{PH} & T_\mathrm{out} = T_0\\

            \begin{cases}
            \sum_i \dot{E}_{\mathrm{in,}i}^\mathrm{PH}
            & T_{\mathrm{in,}i} > T_\mathrm{out} \text{ \& }
            T_{\mathrm{in,}i} \geq T_0 \\
            \sum_i \dot{m}_i \cdot \left(e_{\mathrm{in,}i}^\mathrm{PH} -
            e_\mathrm{out}^\mathrm{PH}\right)
            & T_{\mathrm{in,}i} < T_\mathrm{out} \\
            \end{cases} & T_\mathrm{out} < T_0\\
            \end{cases}

            \forall i \in \text{merge inlets}

            \dot{E}_\mathrm{bus} = \text{not defined (nan)}
        """
        self.E_P = 0
        self.E_F = 0
        if self.outl[0].T.val_SI > T0:
            for i in self.inl:
                if i.T.val_SI < self.outl[0].T.val_SI:
                    if i.T.val_SI >= T0:
                        self.E_P += i.m.val_SI * (
                            self.outl[0].ex_physical - i.ex_physical)
                    else:
                        self.E_P += i.m.val_SI * self.outl[0].ex_physical
                        self.E_F += i.Ex_physical
                else:
                    self.E_F += i.m.val_SI * (
                        i.ex_physical - self.outl[0].ex_physical)
        elif self.outl[0].T.val_SI == T0:
            self.E_P = np.nan
            for i in self.inl:
                self.E_F += i.Ex_physical
        else:
            for i in self.inl:
                if i.T.val_SI > self.outl[0].T.val_SI:
                    if i.T.val_SI >= T0:
                        self.E_P += i.m.val_SI * self.outl[0].ex_physical
                        self.E_F += i.Ex_physical
                    else:
                        self.E_P += i.m.val_SI * (
                            self.outl[0].ex_physical - i.ex_physical)
                else:
                    self.E_F += i.m.val_SI * (
                        i.ex_physical - self.outl[0].ex_physical)

        self.E_bus = {
            "chemical": np.nan, "physical": np.nan, "massless": np.nan
        }
        self.E_D = self.E_F - self.E_P
        self.epsilon = self._calc_epsilon()

    def get_plotting_data(self):
        """Generate a dictionary containing FluProDia plotting information.

        Returns
        -------
        data : dict
            A nested dictionary containing the keywords required by the
            :code:`calc_individual_isoline` method of the
            :code:`FluidPropertyDiagram` class. First level keys are the
            connection index ('in1' -> 'out1', therefore :code:`1` etc.).
        """
        return {
            i + 1: {
                'isoline_property': 'p',
                'isoline_value': self.inl[i].p.val,
                'isoline_value_end': self.outl[0].p.val,
                'starting_point_property': 'v',
                'starting_point_value': self.inl[i].vol.val,
                'ending_point_property': 'v',
                'ending_point_value': self.outl[0].vol.val
            } for i in range(self.num_i)}
