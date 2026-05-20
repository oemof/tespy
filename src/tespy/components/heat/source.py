# -*- coding: utf-8

"""Module of class HeatSource.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/heat/source.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.energy._source import _EnergySource


@component_registry
class HeatSource(_EnergySource):
    r"""
    A heat flow emerges from a HeatSource.

    Ports
    -----

    Heat outlets: heat

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
    Create a HeatSource and specify a label.

    >>> from tespy.components import HeatSource
    >>> so = HeatSource('a labeled heat source')
    >>> so.label
    'a labeled heat source'
    """

    _energy_port = "heat"
