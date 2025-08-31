# -*- coding: utf-8

"""Module of class PowerBus.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/power/bus.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import Component
from tespy.components.component import component_registry
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import SimpleDataContainer as dc_simple


@component_registry
class PowerBus(Component):
    """
    A PowerBus can hold any number incoming and outgoing power flows.

    For example, it can be used to model single shaft gas turbine systems or to
    calculate the net power generation of a rankine cycle plant

    **Mandatory Equations**

    - :py:meth:`tespy.components.power.bus.PowerBus.energy_balance_func`

    Inlets/Outlets

    - None

    PowerConnection inlets/outlets

    - specify number of inlets with :code:`num_in`: 'power_in1', ...
    - specify number of outlets with :code:`num_out` 'power_out1', ...

    Image

    .. image:: /api/_images/PowerBus.svg
       :alt: flowsheet of the power bus
       :align: center
       :class: only-light

    .. image:: /api/_images/PowerBus_darkmode.svg
       :alt: flowsheet of the power bus
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

    num_in : float
        Number of inlets

    num_out : float
        Number of outlets

    Example
    -------
    In a very simple example, a PowerBus is utilized to distribute power from
    the grid to 3 different consumers.

    >>> from tespy.components import PowerSource, PowerSink, PowerBus
    >>> from tespy.connections import PowerConnection
    >>> from tespy.networks import Network
    >>> import os
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC"
    ... })

    We can add a PowerSource representing the grid and three PowerSink
    components representing different power demands.

    >>> grid = PowerSource('grid')
    >>> bus = PowerBus('power bus', num_in=1, num_out=3)
    >>> demand1 = PowerSink('power demand 1')
    >>> demand2 = PowerSink('power demand 2')
    >>> demand3 = PowerSink('power demand 3')
    >>> e1 = PowerConnection(grid, 'power', bus, 'power_in1')
    >>> e2 = PowerConnection(bus, 'power_out1', demand1, 'power')
    >>> e3 = PowerConnection(bus, 'power_out2', demand2, 'power')
    >>> e4 = PowerConnection(bus, 'power_out3', demand3, 'power')
    >>> nw.add_conns(e1, e2, e3, e4)

    We have 4 variables (4 energy flows) and one equation (bus energy balance)
    in our system. That means, we have to fix three values of the variables,
    e.g. we can fix the three demand values:

    >>> e2.set_attr(E=10e3)
    >>> e3.set_attr(E=20e3)
    >>> e4.set_attr(E=30e3)
    >>> nw.solve('design')
    >>> nw.assert_convergence()
    >>> round(e1.E.val_SI) == 60000
    True
    """

    def powerinlets(self):
        return [f"power_in{i + 1}" for i in range(self.num_in.val)]

    def poweroutlets(self):
        return [f"power_out{i + 1}" for i in range(self.num_out.val)]

    def get_parameters(self):
        return {
            "num_in": dc_simple(val=0),
            "num_out": dc_simple(val=0)
        }

    def get_mandatory_constraints(self):
        return {
            "energy_balance_constraint": dc_cmc(**{
                "func": self.energy_balance_func,
                "dependents": self.energy_balance_dependents,
                "num_eq_sets": 1
            })
        }

    def energy_balance_func(self):
        r"""
        Equation for energy balance of the component

        Returns
        -------
        residual : float
            Residual value of equation

            .. math::

                0=\sum_{i} \dot E_\text{i} - \sum_{o} \dot E_\text{o}\\
                \forall i \in \text{inlets}, o \in \text{outlets}
        """
        residual = 0
        for i in self.power_inl:
            residual += i.E.val_SI
        for o in self.power_outl:
            residual -= o.E.val_SI
        return residual

    def energy_balance_dependents(self):
        return [c.E for c in self.power_inl + self.power_outl]
