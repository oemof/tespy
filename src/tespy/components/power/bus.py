# -*- coding: utf-8

"""Module of class PowerBus.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/power/bus.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.energy._bus import _EnergyBus


@component_registry
class PowerBus(_EnergyBus):
    """
    A PowerBus can hold any number incoming and outgoing power flows.

    For example, it can be used to model single shaft gas turbine systems or to
    calculate the net power generation of a rankine cycle plant

    .. image:: /api/_images/components/PowerBus.svg
       :alt: flowsheet of the powerbus
       :align: center
       :class: only-light

    .. image:: /api/_images/components/PowerBus_darkmode.svg
       :alt: flowsheet of the powerbus
       :align: center
       :class: only-dark

    Ports
    -----

    Power inlets: power_in1, power_in2, ... (variable, count set by :code:`num_in`)

    Power outlets: power_out1, power_out2, ... (variable, count set by :code:`num_out`)

    Mandatory Equations
    -------------------

    - energy balance over all inflows and outflows: :py:meth:`energy_balance_func <tespy.components.energy._bus._EnergyBus.energy_balance_func>`

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

    num_in : int
        Number of inlets.

    num_out : int
        Number of outlets.

    offdesign : list
        List containing offdesign parameters (stated as String).

    printout : bool
        Include this component in the network's results printout.

    Example
    -------
    In a very simple example, a PowerBus is utilized to distribute power from
    the grid to 3 different consumers.

    >>> from tespy.components import PowerSource, PowerSink, PowerBus
    >>> from tespy.connections import PowerConnection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC"
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
    @classmethod
    def port_schema(cls):
        return {
            "inlets": {"type": "fixed", "ports": []},
            "outlets": {"type": "fixed", "ports": []},
            "powerinlets": {"type": "variable", "parameter": "num_in", "pattern": "power_in{n}", "min": 1},
            "poweroutlets": {"type": "variable", "parameter": "num_out", "pattern": "power_out{n}", "min": 1},
            "heatinlets": {"type": "fixed", "ports": []},
            "heatoutlets": {"type": "fixed", "ports": []},
        }

    _energy_port = "power"
