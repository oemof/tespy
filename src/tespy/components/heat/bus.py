# -*- coding: utf-8

"""Module of class HeatBus.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/heat/bus.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.energy._bus import _EnergyBus


@component_registry
class HeatBus(_EnergyBus):
    """
    A HeatBus can hold any number of incoming and outgoing heat flows.

    It is the heat-domain equivalent of :class:`PowerBus` and follows the
    same usage pattern - replace :class:`~tespy.connections.PowerConnection`
    with :class:`~tespy.connections.HeatConnection` and
    :class:`~tespy.components.PowerSource` / :class:`~tespy.components.PowerSink`
    with :class:`~tespy.components.HeatSource` / :class:`~tespy.components.HeatSink`.

    **Mandatory Equations**

    - :py:meth:`tespy.components.energy._bus._EnergyBus.energy_balance_func`

    HeatConnection inlets/outlets

    - specify number of inlets with :code:`num_in`: 'heat_in1', ...
    - specify number of outlets with :code:`num_out`: 'heat_out1', ...

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
        Number of inlets.

    num_out : float
        Number of outlets.

    Example
    -------
    A HeatBus collects heat from two sources and distributes it to one sink.

    >>> from tespy.components import HeatSource, HeatSink, HeatBus
    >>> from tespy.connections import HeatConnection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> src1 = HeatSource('source 1')
    >>> src2 = HeatSource('source 2')
    >>> bus = HeatBus('heat bus', num_in=2, num_out=1)
    >>> sink = HeatSink('sink')
    >>> h1 = HeatConnection(src1, 'heat', bus, 'heat_in1')
    >>> h2 = HeatConnection(src2, 'heat', bus, 'heat_in2')
    >>> h3 = HeatConnection(bus, 'heat_out1', sink, 'heat')
    >>> nw.add_conns(h1, h2, h3)
    >>> h1.set_attr(E=10e3)
    >>> h2.set_attr(E=20e3)
    >>> nw.solve('design')
    >>> nw.assert_convergence()
    >>> round(h3.E.val_SI) == 30000
    True
    """

    _energy_port = "heat"
