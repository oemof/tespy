# -*- coding: utf-8

"""Module of class EquilibriumSeparator.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/nodes/equilibrium_separator.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.nodes.droplet_separator import DropletSeparator
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.fluid_properties import h_mix_pQ


@component_registry
class EquilibriumSeparator(DropletSeparator):
    r"""
    Separate the phases of a single fluid in any inlet state.

    In contrast to the
    :class:`DropletSeparator <tespy.components.nodes.droplet_separator.DropletSeparator>`
    this component also handles single phase inlet states: The full mass flow
    then leaves through the outlet matching the inlet phase without a change of
    state, while the other outlet sees zero mass flow at saturated state. For
    two-phase inlets the behavior is identical to the parent component.

    .. image:: /api/_images/components/DropletSeparator.svg
       :alt: flowsheet of the equilibriumseparator
       :align: center
       :class: only-light

    .. image:: /api/_images/components/DropletSeparator_darkmode.svg
       :alt: flowsheet of the equilibriumseparator
       :align: center
       :class: only-dark

    Ports
    -----

    - Fluid inlets: in1
    - Fluid outlets: out1, out2

    Mandatory Equations
    -------------------

    - mass balance constraint: :py:meth:`mass_flow_func <tespy.components.nodes.base.NodeBase.mass_flow_func>`
    - energy balance constraint: :py:meth:`energy_balance_func <tespy.components.nodes.equilibrium_separator.EquilibriumSeparator.energy_balance_func>`
    - pressure equality constraints: :py:meth:`pressure_structure_matrix <tespy.components.nodes.base.NodeBase.pressure_structure_matrix>`
    - outlet 0 liquid state constraint: :py:meth:`outlet_state_func <tespy.components.nodes.equilibrium_separator.EquilibriumSeparator.outlet_state_func>`
    - outlet 1 gas state constraint: :py:meth:`outlet_state_func <tespy.components.nodes.equilibrium_separator.EquilibriumSeparator.outlet_state_func>`
    - fluid equality constraints: :py:meth:`fluid_structure_matrix <tespy.components.nodes.droplet_separator.DropletSeparator.fluid_structure_matrix>`

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

    Example
    -------
    With a two-phase inlet the equilibrium separator splits the mass flow
    according to the vapor mass fraction, just like the droplet separator.

    >>> from tespy.components import Sink, Source, EquilibriumSeparator
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC"
    ... })
    >>> so = Source('inflow')
    >>> sig = Sink('gas outflow')
    >>> sil = Sink('liquid outflow')
    >>> es = EquilibriumSeparator('separator')
    >>> c1 = Connection(so, 'out1', es, 'in1', label='1')
    >>> c2 = Connection(es, 'out1', sil, 'in1', label='2')
    >>> c3 = Connection(es, 'out2', sig, 'in1', label='3')
    >>> nw.add_conns(c1, c2, c3)
    >>> c1.set_attr(fluid={'water': 1}, p=1, x=0.6, m=10)
    >>> nw.solve('design')
    >>> round(c3.m.val_SI, 6)
    6.0

    With a subcooled liquid inlet all mass leaves through the liquid outlet
    at inlet enthalpy, the gas outlet sees zero mass flow at saturated state.

    >>> c1.set_attr(x=None, td_bubble=10)
    >>> nw.solve('design')
    >>> round(c3.m.val_SI, 6)
    0.0
    >>> round(c2.h.val_SI - c1.h.val_SI, 6)
    0.0
    >>> round(c3.calc_Q(), 6)
    1.0
    """

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        constraints['outlet_constraint_liquid'] = dc_cmc(**{
            'func': self.outlet_state_func,
            'dependents': self.outlet_state_dependents,
            'num_eq_sets': 1,
            'func_params': {'outconn': 0, 'phase': 'l'},
            'description': 'outlet 0 liquid state constraint'
        })
        constraints['outlet_constraint_gas'] = dc_cmc(**{
            'func': self.outlet_state_func,
            'dependents': self.outlet_state_dependents,
            'num_eq_sets': 1,
            'func_params': {'outconn': 1, 'phase': 'g'},
            'description': 'outlet 1 gas state constraint'
        })
        return constraints

    def energy_balance_func(self):
        r"""
        Calculate energy balance.

        For a two-phase inlet the energy balance is applied. For a single
        phase inlet it is replaced by a zero mass flow condition on the outlet
        opposite of the inlet phase, the energy balance is then implied by the
        outlet state equations.

        Returns
        -------
        residual : float
            Residual value of energy balance.

            .. math::

                0 = \begin{cases}
                \dot{m}_{in} \cdot h_{in} -
                \dot{m}_{out,1} \cdot h_{out,1} -
                \dot{m}_{out,2} \cdot h_{out,2} & \text{two-phase inlet}\\
                \dot{m}_{out,2} & \text{liquid inlet}\\
                \dot{m}_{out,1} & \text{gas inlet}
                \end{cases}
        """
        phase = self.inl[0].calc_phase()
        if phase == "l":
            return self.outl[1].m.val_SI
        elif phase == "g":
            return self.outl[0].m.val_SI
        elif phase == "tp":
            return super().energy_balance_func()
        else:
            msg = (
                "The equilibrium separator inlet must be liquid, gas or "
                f"two-phase, but the phase is '{phase}'."
            )
            raise ValueError(msg)

    def outlet_state_func(self, outconn=None, phase=None):
        r"""
        Calculate the state of the specified outlet.

        If the inlet phase matches the outlet's phase, the state passes
        through unchanged. Otherwise the outlet is at saturated state, i.e.
        for a two-phase inlet both outlets are saturated, for a single phase
        inlet the opposite outlet is saturated at zero mass flow.

        Returns
        -------
        residual : float
            Residual value of outlet state equation.

            .. math::

                0 = \begin{cases}
                h_{in} - h_{out} & \text{inlet phase matches outlet}\\
                h\left(p_{out}, x \right) - h_{out} & \text{otherwise}
                \end{cases}
        """
        o = self.outl[outconn]
        if self.inl[0].calc_phase() == phase:
            return self.inl[0].h.val_SI - o.h.val_SI
        else:
            quality = 0 if phase == "l" else 1
            return h_mix_pQ(o.p.val_SI, quality, o.fluid_data) - o.h.val_SI

    def outlet_state_dependents(self, outconn=None, phase=None):
        return [
            self.inl[0].h,
            self.outl[outconn].p,
            self.outl[outconn].h
        ]
