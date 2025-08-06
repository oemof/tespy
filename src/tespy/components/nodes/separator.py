# -*- coding: utf-8

"""Module of class Separator.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/nodes/separator.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.nodes.base import NodeBase
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import dT_mix_dph
from tespy.tools.fluid_properties import dT_mix_pdh

# from tespy.tools.fluid_properties import dT_mix_ph_dfluid


@component_registry
class Separator(NodeBase):
    r"""
    A separator separates fluid components from a mass flow.

    **Mandatory Equations**

    - :py:meth:`tespy.components.nodes.base.NodeBase.mass_flow_func`
    - :py:meth:`tespy.components.nodes.base.NodeBase.pressure_structure_matrix`
    - :py:meth:`tespy.components.nodes.separator.Separator.fluid_func`
    - :py:meth:`tespy.components.nodes.separator.Separator.energy_balance_func`

    Inlets/Outlets

    - in1
    - specify number of outlets with :code:`num_out` (default value: 2)

    Image

    .. image:: /api/_images/Splitter.svg
       :alt: flowsheet of the splitter
       :align: center
       :class: only-light

    .. image:: /api/_images/Splitter_darkmode.svg
       :alt: flowsheet of the splitter
       :align: center
       :class: only-dark

    Note
    ----
    Fluid separation requires power and cooling, equations have not been
    implemented, yet!

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

    num_out : float, dict
        Number of outlets for this component, default value: 2.

    Example
    -------
    The separator is used to split up a single mass flow into a specified
    number of different parts at identical pressure and temperature but
    different fluid composition. Fluids can be separated from each other.

    >>> from tespy.components import Sink, Source, Separator
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network(p_unit='bar', T_unit='C', iterinfo=False)
    >>> so = Source('source')
    >>> si1 = Sink('sink1')
    >>> si2 = Sink('sink2')
    >>> s = Separator('separator', num_out=2)
    >>> inc = Connection(so, 'out1', s, 'in1')
    >>> outg1 = Connection(s, 'out1', si1, 'in1')
    >>> outg2 = Connection(s, 'out2', si2, 'in1')
    >>> nw.add_conns(inc, outg1, outg2)

    An Air (simplified) mass flow of 5 kg/s is split up into two mass flows.
    One mass flow of 1 kg/s containing 10 % oxygen and 90 % nitrogen leaves the
    separator. It is possible to calculate the fluid composition of the second
    mass flow. Specify starting values for the second mass flow fluid
    composition for calculation stability.

    >>> inc.set_attr(fluid={'O2': 0.23, 'N2': 0.77}, p=1, T=20, m=5)
    >>> outg1.set_attr(fluid={'O2': 0.1, 'N2': 0.9}, m=1)
    >>> outg2.set_attr(fluid0={'O2': 0.5, 'N2': 0.5})
    >>> nw.solve('design')
    >>> outg2.fluid.val['O2']
    0.2625

    In the same way, it is possible to specify one of the fluid components in
    the second mass flow instead of the first mass flow. The solver will find
    the mass flows matching the desired composition. 65 % of the mass flow
    will leave the separator at the second outlet the case of 30 % oxygen
    mass fraction for this outlet.

    >>> outg1.set_attr(m=None)
    >>> outg2.set_attr(fluid={'O2': 0.3})
    >>> nw.solve('design')
    >>> outg2.fluid.val['O2']
    0.3
    >>> round(outg2.m.val_SI / inc.m.val_SI, 2)
    0.65
    """

    @staticmethod
    def get_parameters():
        return {'num_out': dc_simple()}

    def _update_num_eq(self):
        self.variable_fluids = set(
            [fluid for c in self.inl + self.outl for fluid in c.fluid.is_var]
        )
        num_fluid_eq = len(self.variable_fluids)
        if num_fluid_eq == 0:
            num_fluid_eq = 1
            self.variable_fluids = [list(self.inl[0].fluid.is_set)[0]]

        self.constraints["fluid_constraints"].num_eq = num_fluid_eq

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': dc_cmc(**{
                'num_eq_sets': 1,
                'func': self.mass_flow_func,
                'dependents': self.mass_flow_dependents,
            }),
            'fluid_constraints': dc_cmc(**{
                'num_eq_sets': self.num_o,
                'func': self.fluid_func,
                'deriv': self.fluid_deriv,
                'dependents': self.fluid_dependents
            }),
            'energy_balance_constraints': dc_cmc(**{
                'num_eq_sets': self.num_o,
                'func': self.energy_balance_func,
                'deriv': self.energy_balance_deriv,
                'dependents': self.energy_balance_dependents
            }),
            'pressure_constraints': dc_cmc(**{
                'num_eq_sets': self.num_o,
                'structure_matrix': self.pressure_structure_matrix,
            })
        }

    @staticmethod
    def inlets():
        return ['in1']

    def outlets(self):
        if self.num_out.is_set:
            return [f'out{i + 1}' for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=2)
            return self.outlets()

    def propagate_wrapper_to_target(self, branch):
        branch["components"] += [self]
        for outconn in self.outl:
            branch["connections"] += [outconn]
            outconn.target.propagate_wrapper_to_target(branch)

    def fluid_func(self):
        r"""
        Calculate the vector of residual values for fluid balance equations.

        Returns
        -------
        residual : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = \dot{m}_{in} \cdot x_{fl,in} - \dot {m}_{out,j}
                \cdot x_{fl,out,j}\\
                \forall fl \in \text{network fluids,}
                \; \forall j \in \text{outlets}
        """
        i = self.inl[0]
        residual = []
        for fluid in self.variable_fluids:
            res = i.fluid.val[fluid] * i.m.val_SI
            for o in self.outl:
                res -= o.fluid.val[fluid] * o.m.val_SI
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
        i = self.inl[0]
        for fluid in self.variable_fluids:
            for o in self.outl:
                self._partial_derivative(o.m, k, -o.fluid.val[fluid], increment_filter)
                if fluid in o.fluid.is_var:
                    self.jacobian[k, o.fluid.J_col[fluid]] = -o.m.val_SI

            self._partial_derivative(i.m, k, i.fluid.val[fluid], increment_filter)
            if fluid in i.fluid.is_var:
                self.jacobian[k, i.fluid.J_col[fluid]] = i.m.val_SI

            k += 1

    def fluid_dependents(self):
        return {
            "scalars": [
                [c.m for c in self.inl + self.outl]
                for f in self.variable_fluids
            ],
            "vectors": [{
                c.fluid: set(f) & c.fluid.is_var for c in self.inl + self.outl
            } for f in self.variable_fluids]
        }

    def energy_balance_func(self):
        r"""
        Calculate energy balance.

        Returns
        -------
        residual : list
            Residual value of energy balance.

            .. math::

                0 = T_{in} - T_{out,j}\\
                \forall j \in \text{outlets}
        """
        residual = []
        T_in = self.inl[0].calc_T()
        for o in self.outl:
            residual += [T_in - o.calc_T()]
        return residual

    def energy_balance_deriv(self, increment_filter, k, dependents=None):
        r"""
        Calculate partial derivatives of energy balance.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        i = self.inl[0]
        dT_dp_in = 0
        dT_dh_in = 0
        if i.p.is_var:
            # outlet pressure must be variable as well in this case!
            dT_dp_in = dT_mix_dph(i.p.val_SI, i.h.val_SI, i.fluid_data, i.mixing_rule)
        if i.h.is_var:
            dT_dh_in = dT_mix_pdh(i.p.val_SI, i.h.val_SI, i.fluid_data, i.mixing_rule)

        for o in self.outl:
            args = (o.p.val_SI, o.h.val_SI, o.fluid_data, o.mixing_rule)

            dT_dp_out = 0
            if o.p.is_var:
                dT_dp_out = -dT_mix_dph(*args)
            # pressure is always coupled
            self._partial_derivative(i.p, k, dT_dp_in - dT_dp_out)

            if o.h.is_var:
                dT_dh_out = -dT_mix_pdh(*args)

            # enthalpy is not necessarily coupled
            if i.h._reference_container == o.h._reference_container:
                self._partial_derivative(i.h, k, dT_dh_in - dT_dh_out)
            else:
                self._partial_derivative(i.h, k, dT_dh_in)
                self._partial_derivative(o.h, k, dT_dh_out)

            k += 1

    def energy_balance_dependents(self):
        return [
            [self.inl[0].p, self.inl[0].h, o.p, o.h] for o in self.outl
        ]
