# -*- coding: utf-8

"""Module for custom component groups.

It is possible to create subsystems of component groups in tespy. The subsystem
class is the base class for custom subsystems.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/subsystems.py

SPDX-License-Identifier: MIT
"""

from tespy.components import SubsystemInterface
from tespy.tools import logger


class Subsystem:
    r"""
    Class Subsystem is the base class of all TESPy subsystems.

    Parameters
    ----------
    label : str
        The label of the subsystem.

    Example
    -------
    Basic example for a setting up a Subsystem object. This example does not
    run a TESPy calculation!

    >>> from tespy.components import Subsystem
    >>> class MySubsystem(Subsystem):
    ...     def create_network(self):
    ...         pass
    >>> mysub = MySubsystem('mySubsystem')
    >>> type(mysub)
    <class 'tespy.components.subsystem.MySubsystem'>
    >>> mysub.label
    'mySubsystem'
    >>> type(mysub.inlet)
    <class 'tespy.components.basics.subsystem_interface.SubsystemInterface'>
    >>> type(mysub.outlet)
    <class 'tespy.components.basics.subsystem_interface.SubsystemInterface'>

    If you want to connect to the subsystem from outside of it in a Network,
    then you have to pass the respective number of inlet and outlet connections.
    The number is to your choice, but  for the `Subsystem` to be functional,
    all of the available interfaces must be wired properly internally in the
    :code:`create_network` method. For example, consider a subsystem which is
    just passing its inlet to the outlet:

    >>> from tespy.components import Source, Sink
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> class MySubsystem(Subsystem):
    ...     def __init__(self, label):
    ...         self.num_in = 1
    ...         self.num_out = 1
    ...         super().__init__(label)
    ...
    ...     def create_network(self):
    ...         c1 = Connection(self.inlet, "out1", self.outlet, "in1", label="1")
    ...         self.add_conns(c1)
    >>> mysub = MySubsystem('mySubsystem')
    >>> nw = Network()
    >>> so = Source("source")
    >>> si = Sink("sink")
    >>> c1 = Connection(so, "out1", mysub, "in1", label="1")
    >>> c2 = Connection(mysub, "out1", si, "in1", label="2")
    >>> nw.add_conns(c1, c2)
    >>> nw.add_subsystems(mysub)

    We can run the :code:`check_network` method to check if everything is
    properly connected and a valid topology was created, without needing to
    parametrize the system (for the sake of simplicity in this example).

    >>> nw.check_network()

    You can retrieve components and connections from inside the subsystem with
    their label, which is used inside the :code:`create_network` method of the
    subsystem.

    >>> type(mysub.get_conn("1"))
    <class 'tespy.connections.connection.Connection'>
    >>> type(mysub.get_comp("inlet"))
    <class 'tespy.components.basics.subsystem_interface.SubsystemInterface'>

    Their actual label is prefixed with the subsystem's label, and therefore to
    get it from the network level, you must use that label:

    >>> mysub.get_conn("1").label
    'mySubsystem_1'
    >>> type(nw.get_conn('mySubsystem_1'))
    <class 'tespy.connections.connection.Connection'>

    The same is true for components:

    >>> mysub.get_comp("inlet").label
    'mySubsystem_inlet'
    >>> type(nw.get_comp("mySubsystem_inlet"))
    <class 'tespy.components.basics.subsystem_interface.SubsystemInterface'>
    """

    def __init__(self, label):

        forbidden = [';', ', ', '.']
        if not isinstance(label, str):
            msg = 'Subsystem label must be of type str!'
            logger.error(msg)
            raise ValueError(msg)

        elif len([x for x in forbidden if x in label]) > 0:
            msg = f'Can\'t use {", ".join(forbidden)} in label.'
            logger.error(msg)
            raise ValueError(msg)
        else:
            self.label = label

        self.comps = {}
        self.conns = {}

        if not hasattr(self, "num_in"):
            msg = (
                "When creating your own Subsystem class you need to define "
                "the number of inlets (interfaces to external parts) before "
                "calling super.__init__() in the 'self.num_in' attribute."
            )
            logger.warning(msg)
            self.num_in = 0

        if not hasattr(self, "num_out"):
            msg = (
                "When creating your own Subsystem class you need to define "
                "the number of outlets (interfaces to external parts) before "
                "calling super.__init__() in the 'self.num_out' attribute."
            )
            logger.warning(msg)
            self.num_out = 0

        if self.num_in == 0 and self.num_out == 0:
            msg = (
                "Your subsystem has no interfaces at all. To make interfaces "
                "available to connect to outside of the subsystem components "
                "you have to provide a number of inlets and a number of "
                "outlets in the same style as they are provided for "
                "component classes, i.e. by defining the 'inlets' and "
                "'outlets' methods."
            )
            logger.warning(msg)

        self.inlet = SubsystemInterface("inlet", num_inter=self.num_in)
        self.outlet = SubsystemInterface("outlet", num_inter=self.num_out)

        self.create_network()

    def add_conns(self, *args):

        for conn in args:
            self.conns[conn.label] = conn
            self.conns[conn.label].label = f"{self.label}_{conn.label}"

        self._add_comps(*args)

    def _add_comps(self, *args):
        # get unique components in new connections
        comps = list({cp for c in args for cp in [c.source, c.target]})
        # add to the dict of components

        for comp in comps:
            if comp.label in self.comps.keys():
                msg = "Component name in subsystem is not unique"
                raise ValueError(msg)

            self.comps[comp.label] = comp
            self.comps[comp.label].label = f"{self.label}_{comp.label}"

    def get_conn(self, label):

        try:
            return self.conns[label]
        except KeyError:
            msg = (
                f"Connection with label {label} not found. Note: The label "
                "should not include the Subsystem label when retrieving the "
                "connection from the subsystem."
            )
            logger.warning(msg)
            return None

    def get_comp(self, label):

        try:
            return self.comps[label]
        except KeyError:
            msg = (
                f"Component with label {label} not found. Note: The label "
                "should not include the Subsystem label when retrieving the "
                "component from the subsystem."
            )
            logger.warning(msg)
            return None

    def create_network(self):
        """Create the subsystem's network."""
        msg = (
            "Your subsystem's network has to be set up through the "
            "'create_network' method. To do this, inherit from the "
            "Subsystem class and overwrite this method."
        )
        raise NotImplementedError(msg)
