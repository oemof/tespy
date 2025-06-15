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
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
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
            'mass_flow_constraints': dc_cmc(**{
                'func': self.mass_flow_func,
                'dependents': self.mass_flow_dependents,
                'num_eq_sets': 1
            }),
            'energy_balance_constraints': dc_cmc(**{
                'func': self.energy_balance_func,
                'dependents': self.energy_balance_dependents,
                'num_eq_sets': 1
            }),
            'pressure_constraints': dc_cmc(**{
                'structure_matrix': self.pressure_structure_matrix,
                'num_eq_sets': self.num_i + self.num_o - 1
            }),
            'outlet_constraint_liquid': dc_cmc(**{
                'func': self.saturated_outlet_func,
                'deriv': self.saturated_outlet_deriv,
                'dependents': self.saturated_outlet_dependents,
                'num_eq_sets': 1,
                'func_params': {'outconn': 0, 'quality': 0}
            }),
            'outlet_constraint_gas': dc_cmc(**{
                'func': self.saturated_outlet_func,
                'deriv': self.saturated_outlet_deriv,
                'dependents': self.saturated_outlet_dependents,
                'num_eq_sets': 1,
                'func_params': {'outconn': 1, 'quality': 1}
            }),
            'fluid_constraints': dc_cmc(**{
                'structure_matrix': self.fluid_structure_matrix,
                'num_eq_sets': self.num_o
            })
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

    def energy_balance_dependents(self):
        dependents = []
        for c in self.inl + self.outl:
            dependents += [c.m, c.h]
        return dependents

    def saturated_outlet_func(self, outconn=None, quality=None):
        r"""
        Set the outlet state.

        Returns
        -------
        residual : list
            Residual values of outlet state equations.

            .. math::

                0 = h_{out,1} - h\left(p, x=0 \right)\
        """
        o = self.outl[outconn]
        return h_mix_pQ(o.p.val_SI, quality, o.fluid_data) - o.h.val_SI

    def saturated_outlet_deriv(self, increment_filter, k, dependents=None, outconn=None, quality=None):

        o = self.outl[outconn]
        if o.p.is_var:
            self._partial_derivative(
                o.p, k, dh_mix_dpQ(o.p.val_SI, quality, o.fluid_data), increment_filter
            )
        self._partial_derivative(o.h, k, -1)

    def saturated_outlet_dependents(self, outconn=None, quality=None):
        return [
            self.outl[outconn].p,
            self.outl[outconn].h
        ]

    def fluid_structure_matrix(self, k):
        r"""
        Calculate partial derivatives for all pressure equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        for eq, conn in enumerate(self.outl):
            self._structure_matrix[k + eq, self.inl[0].fluid.sm_col] = 1
            self._structure_matrix[k + eq, conn.fluid.sm_col] = -1

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
