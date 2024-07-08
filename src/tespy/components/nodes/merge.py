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
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import s_mix_pT


@component_registry
class Merge(NodeBase):
    r"""
    Class for merge points with multiple inflows and one outflow.

    **Mandatory Equations**

    - :py:meth:`tespy.components.nodes.base.NodeBase.mass_flow_func`
    - :py:meth:`tespy.components.nodes.base.NodeBase.pressure_equality_func`
    - :py:meth:`tespy.components.nodes.merge.Merge.fluid_func`
    - :py:meth:`tespy.components.nodes.merge.Merge.energy_balance_func`

    Inlets/Outlets

    - specify number of outlets with :code:`num_in` (default value: 2)
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
    >>> import shutil
    >>> nw = Network(p_unit='bar', iterinfo=False)
    >>> so1 = Source('source1')
    >>> so2 = Source('source2')
    >>> so3 = Source('source3')
    >>> si1 = Sink('sink')
    >>> m = Merge('merge', num_in=3)
    >>> m.component()
    'merge'
    >>> inc1 = Connection(so1, 'out1', m, 'in1')
    >>> inc2 = Connection(so2, 'out1', m, 'in2')
    >>> inc3 = Connection(so3, 'out1', m, 'in3')
    >>> outg = Connection(m, 'out1', si1, 'in1')
    >>> nw.add_conns(inc1, inc2, inc3, outg)

    A merge with three inlets mixes air (simplified) with pure nitrogen and
    pure oxygen. All gases enter the component at the same temperature. As
    mixing effects are not considered, the outlet temperature should thus be
    similar to the three inlet temperatures (difference might occur due to
    rounding in fluid property functions, let's check it for two different
    temperatures). It is e.g. possible to find the required mass flow of pure
    nitrogen given the nitrogen mass fraction in the outlet.

    >>> T = 293.15
    >>> inc1.set_attr(fluid={'O2': 0.23, 'N2': 0.77}, p=1, T=T, m=5)
    >>> inc2.set_attr(fluid={'O2': 1}, T=T, m=5)
    >>> inc3.set_attr(fluid={'N2': 1}, T=T)
    >>> outg.set_attr(fluid={'N2': 0.4})
    >>> nw.solve('design')
    >>> round(inc3.m.val_SI, 2)
    0.25
    >>> abs(round((outg.T.val_SI - T) / T, 5)) < 0.01
    True
    >>> T = 173.15
    >>> inc1.set_attr(T=T)
    >>> inc2.set_attr(T=T)
    >>> inc3.set_attr(T=T)
    >>> nw.solve('design')
    >>> abs(round((outg.T.val_SI - T) / T, 5)) < 0.01
    True
    """

    @staticmethod
    def component():
        return 'merge'

    @staticmethod
    def get_parameters():
        return {'num_in': dc_simple()}

    def get_mandatory_constraints(self):
        variable_fluids = set(
            [fluid for c in self.inl + self.outl for fluid in c.fluid.is_var]
        )
        num_fluid_eq = len(variable_fluids)

        if num_fluid_eq == 0:
            num_fluid_eq = len(self.inl[0].fluid.val)
            num_m_eq = 0
        else:
            num_m_eq = 1
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': num_m_eq},
            'fluid_constraints': {
                'func': self.fluid_func, 'deriv': self.fluid_deriv,
                'constant_deriv': False, 'latex': self.fluid_func_doc,
                'num_eq': num_fluid_eq},
            'energy_balance_constraints': {
                'func': self.energy_balance_func,
                'deriv': self.energy_balance_deriv,
                'constant_deriv': False, 'latex': self.energy_balance_func_doc,
                'num_eq': 1},
            'pressure_constraints': {
                'func': self.pressure_equality_func,
                'deriv': self.pressure_equality_deriv,
                'constant_deriv': True,
                'latex': self.pressure_equality_func_doc,
                'num_eq': self.num_i + self.num_o - 1}
        }

    def inlets(self):
        if self.num_in.is_set:
            return ['in' + str(i + 1) for i in range(self.num_in.val)]
        else:
            self.set_attr(num_in=2)
            return self.inlets()

    @staticmethod
    def outlets():
        return ['out1']

    @staticmethod
    def is_branch_source():
        return True

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
        for fluid, x in self.outl[0].fluid.val.items():
            res = -x * self.outl[0].m.val_SI
            for i in self.inl:
                res += i.fluid.val[fluid] * i.m.val_SI
            residual += [res]
        return residual

    def fluid_func_doc(self, label):
        r"""
        Calculate the vector of residual values for fluid balance equations.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = (
            r'0=\sum_i \dot{m}_{\mathrm{in,}i} \cdot x_{fl\mathrm{,in,}i}'
            r'- \dot {m}_\mathrm{out} \cdot x_{fl\mathrm{,out}}'
            r'\; \forall fl \in \text{network fluids,} \; \forall i \in'
            r'\text{inlets}'
        )
        return generate_latex_eq(self, latex, label)

    def fluid_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of fluid balance.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        o = self.outl[0]
        for fluid, x in self.outl[0].fluid.val.items():
            for i in self.inl:
                if i.m.is_var:
                    self.jacobian[k, i.m.J_col] = i.fluid.val[fluid]
                if fluid in i.fluid.is_var:
                    self.jacobian[k, i.fluid.J_col[fluid]] = i.m.val_SI
            if o.m.is_var:
                self.jacobian[k, o.m.J_col] = -x
            if fluid in o.fluid.is_var:
                self.jacobian[k, o.fluid.J_col[fluid]] = -o.m.val_SI
            k += 1

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
        res = -self.outl[0].m.val_SI * self.outl[0].h.val_SI
        for i in self.inl:
            res += i.m.val_SI * i.h.val_SI
        return res

    def energy_balance_func_doc(self, label):
        r"""
        Calculate energy balance.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = (
            r'0=\sum_i\left(\dot{m}_{\mathrm{in,}i}\cdot h_{\mathrm{in,}i}'
            r'\right) - \dot{m}_\mathrm{out} \cdot h_\mathrm{out} '
            r'\; \forall i \in \text{inlets}'
        )
        return generate_latex_eq(self, latex, label)

    def energy_balance_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of energy balance.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        for i in self.inl:
            if i.m.is_var:
                self.jacobian[k, i.m.J_col] = i.h.val_SI
            if i.h.is_var:
                self.jacobian[k, i.h.J_col] = i.m.val_SI
        o = self.outl[0]
        if o.m.is_var:
            self.jacobian[k, o.m.J_col] = -o.h.val_SI
        if o.h.is_var:
            self.jacobian[k, o.h.J_col] = -o.m.val_SI

    @staticmethod
    def is_branch_source():
        return True

    def start_branch(self):
        outconn = self.outl[0]
        branch = {
            "connections": [outconn],
            "components": [self, outconn.target],
            "subbranches": {}
        }
        outconn.target.propagate_to_target(branch)

        return {outconn.label: branch}

    def propagate_to_target(self, branch):
        return

    def propagate_wrapper_to_target(self, branch):
        if self in branch["components"]:
            return

        outconn = self.outl[0]
        branch["connections"] += [outconn]
        branch["components"] += [self]
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
