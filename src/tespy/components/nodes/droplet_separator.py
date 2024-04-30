# -*- coding: utf-8

"""Module of class DropletSeparator.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/nodes/droplet_separator.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.nodes.base import NodeBase
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import dh_mix_dpQ
from tespy.tools.fluid_properties import h_mix_pQ


@component_registry
class DropletSeparator(NodeBase):
    r"""
    Separate liquid phase from gas phase of a single fluid.

    This component is the parent component of the Drum.

    **Mandatory Equations**

    - :py:meth:`tespy.components.nodes.base.NodeBase.mass_flow_func`
    - :py:meth:`tespy.components.nodes.base.NodeBase.pressure_equality_func`
    - :py:meth:`tespy.components.nodes.droplet_separator.DropletSeparator.fluid_func`
    - :py:meth:`tespy.components.nodes.droplet_separator.DropletSeparator.energy_balance_func`
    - :py:meth:`tespy.components.nodes.droplet_separator.DropletSeparator.outlet_states_func`

    Inlets/Outlets

    - in1
    - out1, out2 (index 1: saturated liquid, index 2: saturated gas)

    Image

    .. image:: /api/_images/DropletSeparator.svg
       :alt: flowsheet of the droplet separator
       :align: center
       :class: only-light

    .. image:: /api/_images/DropletSeparator_darkmode.svg
       :alt: flowsheet of the droplet separator
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

    Example
    -------
    The droplet separator separates gas from liquid phase. From a stream of
    water the liquid phase will be separated.

    >>> from tespy.components import Sink, Source, DropletSeparator
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import shutil
    >>> nw = Network(T_unit='C', p_unit='bar', h_unit='kJ / kg', iterinfo=False)
    >>> so = Source('two phase inflow')
    >>> sig = Sink('gas outflow')
    >>> sil = Sink('liquid outflow')
    >>> ds = DropletSeparator('droplet separator')
    >>> ds.component()
    'droplet separator'
    >>> so_ds = Connection(so, 'out1', ds, 'in1')
    >>> ds_sig = Connection(ds, 'out2', sig, 'in1')
    >>> ds_sil = Connection(ds, 'out1', sil, 'in1')
    >>> nw.add_conns(so_ds, ds_sig, ds_sil)

    We specify the fluid's state at the inlet. At the gas outflow saturated
    gas enthalpy is expected, at the liquid gas outflow saturated liquid
    enthalpy. The mass flow at the outlets is expected to split according to
    the vapor mass fraction:

    .. math::

        \dot{m}_\mathrm{out,1} = \left(1 - \frac{h_\mathrm{in} - h'}{h'' - h'}
        \right) \cdot \dot{m}_\mathrm{in}

        \dot{m}_\mathrm{out,2} = \frac{h_\mathrm{in} - h'}{h'' - h'} \cdot
        \dot{m}_\mathrm{in}

    >>> so_ds.set_attr(fluid={'water': 1}, p=1, h=1500, m=10)
    >>> nw.solve('design')
    >>> Q_in = so_ds.calc_Q()
    >>> round(Q_in * so_ds.m.val_SI, 6) == round(ds_sig.m.val_SI, 6)
    True
    >>> round((1 - Q_in) * so_ds.m.val_SI, 6) == round(ds_sil.m.val_SI, 6)
    True
    >>> ds_sig.calc_Q()
    1.0
    >>> ds_sil.calc_Q()
    0.0

    In a different setup, we unset pressure and enthalpy and specify gas
    temperature and mass flow instead. The temperature specification must yield
    the corresponding boiling point pressure and the mass flow must yield the
    inlet enthalpy. The inlet vapor mass fraction must be equal to fraction of
    gas mass flow to inlet mass flow (0.95 in this example).

    >>> so_ds.set_attr(fluid={'water': 1}, p=None, h=None, T=150, m=10)
    >>> ds_sig.set_attr(m=9.5)
    >>> nw.solve('design')
    >>> round(so_ds.calc_Q(), 6)
    0.95
    >>> T_boil = so_ds.calc_T_sat()
    >>> round(T_boil, 6) == round(so_ds.T.val_SI, 6)
    True
    """

    @staticmethod
    def component():
        return 'droplet separator'

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': 1},
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
                'num_eq': self.num_i + self.num_o - 1},
            'outlet_constraints': {
                'func': self.outlet_states_func,
                'deriv': self.outlet_states_deriv,
                'constant_deriv': False,
                'latex': self.outlet_states_func_doc,
                'num_eq': 2}
        }

    @staticmethod
    def inlets():
        return ['in1']

    @staticmethod
    def outlets():
        return ['out1', 'out2']

    def energy_balance_func(self):
        r"""
        Calculate energy balance.

        Returns
        -------
        residual : float
            Residual value of energy balance.

            .. math::

                0 = \sum_i \left(\dot{m}_{in,i} \cdot h_{in,i} \right) -
                \sum_j \left(\dot{m}_{out,j} \cdot h_{out,j} \right)\\
                \forall i \in \text{inlets} \; \forall j \in \text{outlets}
        """
        res = 0
        for i in self.inl:
            res += i.m.val_SI * i.h.val_SI
        for o in self.outl:
            res -= o.m.val_SI * o.h.val_SI

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
            r'\right) - \sum_j \left(\dot{m}_{\mathrm{out,}j} \cdot '
            r'h_{\mathrm{out,}j} \right) \; \forall i \in \text{inlets} \;'
            r'\forall j \in \text{outlets}'
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

        for o in self.outl:
            if o.m.is_var:
                self.jacobian[k, o.m.J_col] = -o.h.val_SI
            if o.h.is_var:
                self.jacobian[k, o.h.J_col] = -o.m.val_SI


    def outlet_states_func(self):
        r"""
        Calculate energy balance.

        Returns
        -------
        residual : list
            Residual values of outlet state equations.

            .. math::

                0 = h_{out,1} - h\left(p, x=0 \right)\\
                0 = h_{out,2} - h\left(p, x=1 \right)
        """
        o0 = self.outl[0]
        o1 = self.outl[1]
        return [
            h_mix_pQ(o0.p.val_SI, 0, o0.fluid_data) - o0.h.val_SI,
            h_mix_pQ(o1.p.val_SI, 1, o1.fluid_data) - o1.h.val_SI
        ]

    def outlet_states_func_doc(self, label):
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
            r'\begin{split}' + '\n'
            r'0 =&h_\mathrm{out,1} -h\left(p_\mathrm{out,1}, x=0\right)\\'
            r'0 =&h_\mathrm{out,2} -h\left(p_\mathrm{out,2}, x=1\right)\\'
            r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def outlet_states_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of outlet states.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        o0 = self.outl[0]
        o1 = self.outl[1]
        if o0.p.is_var:
            self.jacobian[k, o0.p.J_col] = (
                dh_mix_dpQ(o0.p.val_SI, 0, o0.fluid_data)
            )
        if o0.h.is_var and self.it == 0:
            self.jacobian[k, o0.h.J_col] = -1

        if o1.p.is_var:
            self.jacobian[k + 1, o1.p.J_col] = (
                dh_mix_dpQ(o1.p.val_SI, 1, o1.fluid_data)
            )
        if o1.h.is_var and self.it == 0:
            self.jacobian[k + 1, o1.h.J_col] = -1

    def propagate_wrapper_to_target(self, branch):
        if self in branch["components"]:
            return

        for outconn in self.outl:
            branch["connections"] += [outconn]
            branch["components"] += [self]
            outconn.target.propagate_wrapper_to_target(branch)

    @staticmethod
    def initialise_source(c, key):
        r"""
        Return a starting value for pressure and enthalpy at outlet.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                10^6 & \text{key = 'p'}\\
                h\left(p, x=1 \right) & \text{key = 'h' at outlet 1}\\
                h\left(p, x=0 \right) & \text{key = 'h' at outlet 2}
                \end{cases}
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            if c.source_id == 'out1':
                return h_mix_pQ(c.p.val_SI, 0, c.fluid_data)
            else:
                return h_mix_pQ(c.p.val_SI, 1, c.fluid_data)

    @staticmethod
    def initialise_target(c, key):
        r"""
        Return a starting value for pressure and enthalpy at inlet.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                10^6 & \text{key = 'p'}\\
                h\left(p, x=0.5 \right) & \text{key = 'h' at inlet 1}
                \end{cases}
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            return h_mix_pQ(c.p.val_SI, 0.5, c.fluid_data)

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
                'isoline_value': self.inl[0].p.val,
                'isoline_value_end': self.outl[i].p.val,
                'starting_point_property': 'v',
                'starting_point_value': self.inl[0].vol.val,
                'ending_point_property': 'v',
                'ending_point_value': self.outl[i].vol.val
            } for i in range(2)}
