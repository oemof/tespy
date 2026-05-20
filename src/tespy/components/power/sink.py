# -*- coding: utf-8

"""Module of class PowerSink.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/power/sink.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.energy._sink import _EnergySink


@component_registry
class PowerSink(_EnergySink):
    r"""
    A power flow drains in a PowerSink.

    Ports
    -----

    Power inlets: power

    Mandatory Equations
    -------------------

    None

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
    Create a PowerSink and specify a label.

    >>> from tespy.components import PowerSink
    >>> si = PowerSink('a labeled sink')
    >>> si.label
    'a labeled sink'
    """

    _energy_port = "power"
