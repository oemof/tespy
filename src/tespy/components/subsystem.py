# -*- coding: utf-8

"""Module for custom component groups.

It is possible to create subsystems of component groups in tespy. The subsystem
class is the base class for custom subsystems.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/subsystems.py

SPDX-License-Identifier: MIT
"""

import copy

from tespy.components import SubsystemInterface
from tespy.tools import logger
from tespy.tools.helpers import TESPyComponentError
from tespy.tools.helpers import TESPyNetworkError

# connection classes with the boundary component name and port prefix used to
# rewire them in Subsystem.from_network
_BOUNDARY_REWIRING = [
    (["Connection", "HAConnection"], "", ""),
    (["HeatConnection"], "Heat", "heat_"),
    (["PowerConnection"], "Power", "power_"),
]


class Subsystem:
    r"""
    Class Subsystem is the base class of all TESPy subsystems.

    Parameters
    ----------
    label : str
        The label of the subsystem.

    Notes
    -----
    Subclasses must set the following attributes *before* calling
    :code:`super().__init__()` to declare the subsystem's external interfaces:

    - :code:`num_in` / :code:`num_out` - number of fluid inlet/outlet ports
      (exposed as :code:`in{n}` / :code:`out{n}` on :code:`self.inlet` /
      :code:`self.outlet`). Default: 0 with a warning.
    - :code:`num_power_in` / :code:`num_power_out` - number of
      :py:class:`~tespy.connections.powerconnection.PowerConnection` inlet/
      outlet ports (exposed as :code:`power_in{n}` / :code:`power_out{n}`).
      Default: 0, no warning.
    - :code:`num_heat_in` / :code:`num_heat_out` - number of
      :py:class:`~tespy.connections.heatconnection.HeatConnection` inlet/
      outlet ports (exposed as :code:`heat_in{n}` / :code:`heat_out{n}`).
      Default: 0, no warning.

    For every power or heat port pair the underlying
    :py:class:`~tespy.components.basics.subsystem_interface.SubsystemInterface`
    enforces :math:`\dot E_\text{in} = \dot E_\text{out}`, so energy passes
    through the boundary unchanged.

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

    We can run the :code:`check_topology` method to check if everything is
    properly connected and a valid topology was created, without needing to
    parametrize the system (for the sake of simplicity in this example).

    >>> nw.check_topology()

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

        if not isinstance(label, str):
            msg = 'Subsystem label must be of type str!'
            logger.error(msg)
            raise ValueError(msg)

        self.label = label

        self.comps = {}
        self.conns = {}
        self.interface_map = {}

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

        if not hasattr(self, "num_power_in"):
            self.num_power_in = 0
        if not hasattr(self, "num_power_out"):
            self.num_power_out = 0
        if not hasattr(self, "num_heat_in"):
            self.num_heat_in = 0
        if not hasattr(self, "num_heat_out"):
            self.num_heat_out = 0

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

        self.inlet = SubsystemInterface(
            "inlet",
            num_inter=self.num_in,
            num_power_inter=self.num_power_in,
            num_heat_inter=self.num_heat_in,
        )
        self.outlet = SubsystemInterface(
            "outlet",
            num_inter=self.num_out,
            num_power_inter=self.num_power_out,
            num_heat_inter=self.num_heat_out,
        )

        self.create_network()

    def add_conns(self, *args):

        for conn in args:
            if conn.label in self.conns:
                msg = (
                    f"A connection with the label {conn.label} has already "
                    "been added to this Subsystem. All connections must have "
                    "unique labels."
                )
                raise TESPyComponentError(msg)
            self.conns[conn.label] = conn
            self.conns[conn.label].label = f"{self.label}_{conn.label}"

        self._add_comps(*args)

    def _add_comps(self, *args):
        # get unique components in new connections and remove existing ones
        comps = (
            {cp for c in args for cp in [c.source, c.target]}
            - set(self.comps.values())
        )
        for comp in comps:
            if comp.label in self.comps.keys():
                msg = "Component name in subsystem is not unique"
                raise TESPyComponentError(msg)

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

    @classmethod
    def from_dict(cls, subsystem_data, label=None):
        r"""
        Create a :code:`Subsystem` instance from serialized state.

        Parameters
        ----------
        subsystem_data : dict
            Serialized structure and parametrization of the subsystem, e.g.
            created through :py:meth:`from_network`. The data are not modified,
            so they can be reused to create multiple instances.

        label : str, optional
            Label of the subsystem. If not provided, the label stored in the
            :code:`Subsystem` section of the data is used.

        Returns
        -------
        tespy.components.subsystem.Subsystem
            The subsystem instance.
        """
        from tespy.components.component import component_registry
        from tespy.connections.connection import connection_registry
        from tespy.networks.network import _construct_components
        from tespy.networks.network import _construct_connections

        subsystem_data = copy.deepcopy(subsystem_data)

        if label is None:
            label = subsystem_data["Subsystem"].get("label")
        if label is None:
            msg = (
                "A label has to be provided, either through the label keyword "
                "or through the 'label' key in the 'Subsystem'  section of "
                "the data."
            )
            logger.error(msg)
            raise ValueError(msg)
        elif not isinstance(label, str):
            msg = "Subsystem label must be of type str!"
            logger.error(msg)
            raise ValueError(msg)

        sub = cls.__new__(cls)
        sub.label = label
        for attr in [
                "num_in", "num_out", "num_heat_in", "num_heat_out",
                "num_power_in", "num_power_out"]:
            setattr(sub, attr, subsystem_data["Subsystem"][attr])
        sub.interface_map = subsystem_data["Subsystem"].get(
            "interface_map", {}
        )
        sub.conns = {}
        sub.comps = {}

        comps = {}
        for component, data in subsystem_data["Component"].items():
            if component not in component_registry.items:
                msg = (
                    f"A class {component} is not available through the "
                    "tespy.components.component.component_registry decorator. "
                    "If you are using a custom component make sure to "
                    "decorate the class."
                )
                logger.error(msg)
                raise TESPyNetworkError(msg)

            target_class = component_registry.items[component]
            # no unit context here: values stay plain with their unit string
            # and are bound to the registry of the network the subsystem
            # becomes part of
            comps.update(_construct_components(target_class, data))

        try:
            sub.inlet = comps["inlet"]
            sub.outlet = comps["outlet"]
        except KeyError:
            msg = (
                "The data must contain SubsystemInterface components with "
                "the labels 'inlet' and 'outlet'."
            )
            logger.error(msg)
            raise TESPyNetworkError(msg)

        conns = {}
        for connection, data in subsystem_data["Connection"].items():
            if connection not in connection_registry.items:
                msg = (
                    f"A class {connection} is not available through the "
                    "tespy.connections.connection.connection_registry "
                    "decorator. If you are using a custom connection make "
                    "sure to decorate the class."
                )
                logger.error(msg)
                raise TESPyNetworkError(msg)

            target_class = connection_registry.items[connection]
            conns.update(_construct_connections(target_class, data, comps))

        for c in conns.values():
            sub.add_conns(c)

        msg = f"Created subsystem {label} from data."
        logger.info(msg)

        return sub

    @classmethod
    def from_network(cls, label, nw, interface_exceptions=None):
        r"""
        Convert a network into a subsystem.

        The system boundary components of the network are removed and the
        connections previously attached to them are rewired to ports of the
        subsystem's interfaces, i.e.:

        - :py:class:`~tespy.components.basics.source.Source`,
        - :py:class:`~tespy.components.basics.sink.Sink`,
        - :py:class:`~tespy.components.heat.source.HeatSource`,
        - :py:class:`~tespy.components.heat.sink.HeatSink`,
        - :py:class:`~tespy.components.power.source.PowerSource` and
        - :py:class:`~tespy.components.power.sink.PowerSink`.

        The subsystem's :code:`interface_map` attribute indicates which source
        or sink was transformed to which port.

        Parameters
        ----------
        label : str
            Label of the subsystem.

        nw : tespy.networks.network.Network
            Network to convert into a subsystem.

        interface_exceptions : list, optional
            Labels of boundary components to keep inside the subsystem instead
            of converting them to interface ports.

        Returns
        -------
        tespy.components.subsystem.Subsystem
            The subsystem instance.

        Example
        -------
        Set up a simple network with a heater and convert it to a subsystem.

        >>> from tespy.components import SimpleHeatExchanger, Sink, Source
        >>> from tespy.components import Subsystem
        >>> from tespy.connections import Connection
        >>> from tespy.networks import Network
        >>> nw = Network(iterinfo=False)
        >>> so = Source("water inlet")
        >>> hx = SimpleHeatExchanger("heater")
        >>> si = Sink("water outlet")
        >>> c1 = Connection(so, "out1", hx, "in1", label="1")
        >>> c2 = Connection(hx, "out1", si, "in1", label="2")
        >>> nw.add_conns(c1, c2)
        >>> c1.set_attr(fluid={"water": 1}, T=20, p=1, m=1)
        >>> sub = Subsystem.from_network("heating block", nw)

        The :code:`interface_map` attribute shows which port of the subsystem
        replaces which boundary component of the network.

        >>> sub.interface_map
        {'in1': 'water inlet', 'out1': 'water outlet'}
        >>> sub.get_comp("heater").label
        'heating block_heater'

        All parameter specifications are taken over into the subsystem,
        including those on the rewired boundary connections.

        >>> sub.get_conn("1").T.is_set
        True
        """
        if interface_exceptions is None:
            interface_exceptions = []

        nw_dict = nw.export()

        for comps in nw_dict["Component"].values():
            for reserved in ("inlet", "outlet"):
                if reserved in comps:
                    msg = (
                        f"The component label '{reserved}' is reserved for "
                        "the interfaces of the subsystem. Please rename the "
                        "component in your network."
                    )
                    logger.error(msg)
                    raise TESPyComponentError(msg)

        interfaces = nw_dict["Component"].setdefault("SubsystemInterface", {})
        interfaces.update(SubsystemInterface("inlet")._serialize())
        interfaces.update(SubsystemInterface("outlet")._serialize())

        nw_dict["Subsystem"] = {"label": label}
        for attr in [
                "num_in", "num_out", "num_heat_in", "num_heat_out",
                "num_power_in", "num_power_out"]:
            nw_dict["Subsystem"][attr] = 0

        interface_map = {}

        for conn_keys, comp_name, prefix in _BOUNDARY_REWIRING:
            for conn_key in conn_keys:
                for con in nw_dict["Connection"].get(conn_key, {}).values():
                    for end, interface, comp_class, count_key, port in [
                            ("source", "inlet", f"{comp_name}Source",
                             f"num_{prefix}in", f"{prefix}out"),
                            ("target", "outlet", f"{comp_name}Sink",
                             f"num_{prefix}out", f"{prefix}in")]:
                        comp = con[end]
                        if (comp not in nw_dict["Component"].get(comp_class, {})
                                or comp in interface_exceptions):
                            continue
                        del nw_dict["Component"][comp_class][comp]
                        nw_dict["Subsystem"][count_key] += 1
                        num = nw_dict["Subsystem"][count_key]
                        con[end] = interface
                        con[f"{end}_id"] = f"{port}{num}"
                        interface_map[
                            f"{count_key.removeprefix('num_')}{num}"
                        ] = comp

        for attr, port_count in [
                ("num_inter", "num_in"),
                ("num_heat_inter", "num_heat_in"),
                ("num_power_inter", "num_power_in")]:
            interfaces["inlet"][attr] = {
                "val": nw_dict["Subsystem"][port_count], "is_set": True
            }
        for attr, port_count in [
                ("num_inter", "num_out"),
                ("num_heat_inter", "num_heat_out"),
                ("num_power_inter", "num_power_out")]:
            interfaces["outlet"][attr] = {
                "val": nw_dict["Subsystem"][port_count], "is_set": True
            }

        nw_dict["Subsystem"]["interface_map"] = interface_map

        return cls.from_dict(nw_dict)
