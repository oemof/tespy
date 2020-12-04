# -*- coding: utf-8

"""Module for class Source.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/basics/source.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import Component


class Source(Component):
    r"""
    A flow originates from a Source.

    Equations
        This component is unconstrained.

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    design_path: str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings: boolean
        Ignore warnings on default characteristics usage for this component.

    printout: boolean
        Include this component in the network's results printout.

    Example
    -------
    Create a source and specify a label.

    >>> from tespy.components import Source
    >>> so = Source('a labeled source')
    >>> so.component()
    'source'
    >>> so.label
    'a labeled source'
    """

    @staticmethod
    def component():
        return 'source'

    @staticmethod
    def outlets():
        return ['out1']

    def propagate_fluid_to_source(self, outconn, start):
        r"""
        Fluid propagation to source stops here.

        Parameters
        ----------
        outconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        return
