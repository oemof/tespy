# -*- coding: utf-8

"""Module of class HeatSink.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/heat/sink.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.energy._sink import _EnergySink


@component_registry
class HeatSink(_EnergySink):
    r"""
    A heat flow drains in a HeatSink.

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
    Create a HeatSink and specify a label.

    >>> from tespy.components import HeatSink
    >>> si = HeatSink('a labeled heat sink')
    >>> si.label
    'a labeled heat sink'
    """

    _energy_port = "heat"
