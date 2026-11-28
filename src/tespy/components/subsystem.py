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
import weakref

from tespy.components import SubsystemInterface
from tespy.tools import logger
from tespy.tools.helpers import TESPyComponentError
from tespy.tools.helpers import TESPyNetworkError

PATH_SEPARATOR = "/"

# the subsystem a component is a direct member of, kept outside of the
# components like the network ownership
_subsystem_membership = weakref.WeakKeyDictionary()


def _owner(component):
    reference = _subsystem_membership.get(component)
    return reference() if reference is not None else None


# connection classes with the boundary component name and port prefix used to
# rewire them in Subsystem.from_network
_BOUNDARY_REWIRING = [
    (["Connection", "HAConnection", "SolutionConnection"], "", ""),
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

    # members of a group keep their labels, the members of other subsystems
    # are prefixed with the labels of the subsystem path
    _prefix_members = True

    def __init__(self, label):

        self._init_structure(label)

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

        self._add_interface_ports()

        self.create_network()

    def _init_structure(self, label):
        self.parent = None
        self.subsystems = {}
        self.comps = {}
        self.conns = {}
        self.ports = {}
        self.interface_map = {}
        self.label = label

    @classmethod
    def _create_group(cls, label):
        group = cls.__new__(cls)
        group._prefix_members = False
        group._init_structure(label)
        return group

    @property
    def is_group(self):
        """The subsystem is a group of existing components without interfaces."""
        return not self._prefix_members

    def _add_interface_ports(self):
        inlet, outlet = self.inlet, self.outlet
        for interface, ports in [
                (inlet, inlet.inlets() + inlet.powerinlets() + inlet.heatinlets()),
                (outlet, outlet.outlets() + outlet.poweroutlets() + outlet.heatoutlets())]:
            for port in ports:
                self.ports[port] = (interface, port)

    def add_port(self, key, component, port):
        r"""
        Make a port of a component available as port of the group.

        Connections to the group through this port are connected to the
        component, e.g. :code:`Connection(source, "out1", group, key)`.

        Parameters
        ----------
        key : str
            Name of the port of the group.

        component : tespy.components.component.Component
            A component of the group or of one of its nested subsystems.

        port : str
            Port id of the component, e.g. :code:`"in1"`.
        """
        if not self.is_group:
            msg = (
                f"The ports of the subsystem {self.path} are given by its "
                "interfaces. Ports can only be added to groups."
            )
        elif component not in self.get_components(recursive=True):
            msg = (
                f"The component {component.label} is not part of the group "
                f"{self.path}."
            )
        elif port not in [
                p for ports in [
                    component.inlets(), component.outlets(),
                    component.powerinlets(), component.poweroutlets(),
                    component.heatinlets(), component.heatoutlets()]
                for p in ports]:
            msg = (
                f"The component {component.label} has no port {port}."
            )
        elif key in self.ports:
            msg = f"The group {self.path} already has a port {key}."
        else:
            self.ports[key] = (component, port)
            return
        logger.error(msg)
        raise TESPyComponentError(msg)

    def _resolve_port(self, key):
        try:
            return self.ports[key]
        except KeyError:
            msg = (
                f"Error creating connection. The subsystem {self.path} has no "
                f"port {key}. Select one of the following ports "
                f"{', '.join(self.ports)}."
            )
            logger.error(msg)
            raise ValueError(msg)

    @property
    def label(self):
        """Label of the subsystem, unique among its siblings."""
        return self._label

    @label.setter
    def label(self, label):
        if not isinstance(label, str):
            msg = "Subsystem label must be of type str!"
            logger.error(msg)
            raise ValueError(msg)
        if PATH_SEPARATOR in label:
            msg = (
                f"The label {label!r} of the subsystem must not contain "
                f"'{PATH_SEPARATOR}', it separates the levels of a subsystem "
                "path."
            )
            logger.error(msg)
            raise ValueError(msg)
        if self.parent is not None or self._network() is not None:
            msg = (
                f"The subsystem {self.path} cannot be renamed after it has "
                "been added to another subsystem or to a network."
            )
            logger.error(msg)
            raise TESPyComponentError(msg)
        self._label = label
        self._relabel()

    @property
    def path(self):
        """Labels of the subsystem and of its parents, separated by '/'."""
        return PATH_SEPARATOR.join(sub.label for sub in self._lineage())

    def _lineage(self):
        lineage = [self]
        while lineage[0].parent is not None:
            lineage.insert(0, lineage[0].parent)
        return lineage

    def _walk(self):
        yield self
        for child in self.subsystems.values():
            yield from child._walk()

    def _contains(self, subsystem):
        return any(sub is subsystem for sub in self._walk())

    def _network(self):
        from tespy.networks.network import _network_ownership

        reference = _network_ownership.get(self._lineage()[0])
        return reference() if reference is not None else None

    def _interfaces(self):
        return [
            comp for comp in self.comps.values()
            if isinstance(comp, SubsystemInterface)
        ]

    def _global_label(self, local):
        if not self._prefix_members:
            return local
        prefix = "".join(
            f"{sub.label}_" for sub in self._lineage() if sub._prefix_members
        )
        return f"{prefix}{local}"

    def _relabel(self):
        # the global labels are derived from the local labels, never from the
        # current labels, so the prefix is applied once per path
        for local, comp in self.comps.items():
            comp.label = self._global_label(local)
        for local, conn in self.conns.items():
            conn.label = self._global_label(local)
        for child in self.subsystems.values():
            child._relabel()

    def _check_unique_labels(self):
        for kind, members in [
                ("component", lambda sub: sub.comps),
                ("connection", lambda sub: sub.conns)]:
            seen = {}
            for sub in self._walk():
                for local in members(sub):
                    label = sub._global_label(local)
                    if label in seen:
                        other_path, other_local = seen[label]
                        msg = (
                            f"The {kind} {local!r} of the subsystem {sub.path} "
                            f"and the {kind} {other_local!r} of the subsystem "
                            f"{other_path} both get the label {label!r}. "
                            "Rename one of them."
                        )
                        logger.error(msg)
                        raise TESPyComponentError(msg)
                    seen[label] = (sub.path, local)

    def add_conns(self, *args):
        r"""
        Add connections to the subsystem.

        The components at both ends of the connections become members of the
        subsystem, unless they are the interfaces of another subsystem: that
        subsystem is added as a child, see
        :py:meth:`~tespy.components.subsystem.Subsystem.add_subsystems`.

        Parameters
        ----------
        *args : tespy.connections.connection.ConnectionBase
            The connections to add, labelled with their local label.
        """
        for conn in args:
            if conn.label in self.conns:
                msg = (
                    f"A connection with the label {conn.label} has already "
                    "been added to this Subsystem. All connections must have "
                    "unique labels."
                )
                raise TESPyComponentError(msg)
            self.conns[conn.label] = conn

        self._add_comps(*args)

        for conn in args:
            conn.label = self._global_label(conn.label)

    def _add_comps(self, *args):
        comps = dict.fromkeys(cp for c in args for cp in [c.source, c.target])
        for comp in comps:
            owner = _owner(comp)
            if owner is None:
                self._add_member(comp)
            elif self._contains(owner):
                continue
            elif comp in owner._interfaces() and owner.parent is None:
                self.add_subsystems(owner)
            else:
                msg = (
                    f"The component {comp.label} is part of the subsystem "
                    f"{owner.path}. Connect to the interfaces of that "
                    "subsystem instead."
                )
                logger.error(msg)
                raise TESPyComponentError(msg)

    def _add_member(self, comp):
        local = comp.label
        if local in self.comps:
            msg = (
                f"A component with the label {local} has already been added "
                f"to the subsystem {self.path}. All components must have "
                "unique labels."
            )
            logger.error(msg)
            raise TESPyComponentError(msg)
        self.comps[local] = comp
        _subsystem_membership[comp] = weakref.ref(self)
        comp.label = self._global_label(local)

    def add_comps(self, *args):
        r"""
        Add components to the subsystem.

        Components at the ends of the connections added with
        :py:meth:`~tespy.components.subsystem.Subsystem.add_conns` are added
        automatically. Use this method for components without connections
        inside of the subsystem.

        Parameters
        ----------
        *args : tespy.components.component.Component
            The components to add, labelled with their local label.
        """
        for comp in args:
            owner = _owner(comp)
            if owner is self:
                continue
            if owner is not None:
                msg = (
                    f"The component {comp.label} is already part of the "
                    f"subsystem {owner.path}."
                )
                logger.error(msg)
                raise TESPyComponentError(msg)
            self._add_member(comp)

    def add_subsystems(self, *args):
        r"""
        Add subsystems as children of this subsystem.

        The labels of the members of the children are prefixed with the label
        of this subsystem. Connecting to the interfaces of a subsystem inside
        of :code:`create_network` adds it automatically.

        Parameters
        ----------
        *args : tespy.components.subsystem.Subsystem
            The subsystems to add.
        """
        if self._network() is not None:
            msg = (
                f"Add the subsystems to {self.path} before adding it to a "
                "network."
            )
            logger.error(msg)
            raise TESPyComponentError(msg)
        for child in args:
            if child._contains(self):
                msg = (
                    f"The subsystem {child.label} cannot be added to "
                    f"{self.path}, because it contains it."
                )
            elif child.parent is not None:
                msg = (
                    f"The subsystem {child.label} is already part of the "
                    f"subsystem {child.parent.path}."
                )
            elif child._network() is not None:
                msg = (
                    f"The subsystem {child.label} is already part of a "
                    "network. Add it to the subsystem before adding the "
                    "subsystem to the network."
                )
            elif child.label in self.subsystems:
                msg = (
                    f"The subsystem {self.path} already contains a subsystem "
                    f"with the label {child.label}."
                )
            else:
                msg = None
            if msg is not None:
                logger.error(msg)
                raise TESPyComponentError(msg)

            self.subsystems[child.label] = child
            child.parent = self
            try:
                self._lineage()[0]._check_unique_labels()
            except TESPyComponentError:
                del self.subsystems[child.label]
                child.parent = None
                raise
            child._relabel()

    def get_subsystem(self, label):
        r"""
        Get a child subsystem by its label or its path relative to this one.

        Parameters
        ----------
        label : str
            Label of the child, or path of a deeper subsystem, e.g.
            :code:`"drive/pump group"`.

        Returns
        -------
        tespy.components.subsystem.Subsystem
        """
        sub = self
        for part in label.split(PATH_SEPARATOR):
            try:
                sub = sub.subsystems[part]
            except KeyError:
                msg = (
                    f"Subsystem with label {label} not found in subsystem "
                    f"{self.path}."
                )
                logger.error(msg)
                raise KeyError(msg)
        return sub

    def get_components(self, recursive=False):
        r"""
        Get the components of the subsystem.

        Parameters
        ----------
        recursive : boolean
            Include the components of the nested subsystems, default:
            :code:`False`.

        Returns
        -------
        list
        """
        subsystems = self._walk() if recursive else [self]
        return [comp for sub in subsystems for comp in sub.comps.values()]

    def get_connections(self, recursive=False):
        r"""
        Get the connections of the subsystem.

        A connection belongs to the deepest subsystem that contains the
        components at both of its ends. Once the subsystem is part of a
        network, this rule is applied to all connections of the network,
        before that the connections added to the subsystems are returned.

        Parameters
        ----------
        recursive : boolean
            Include the connections of the nested subsystems, default:
            :code:`False`.

        Returns
        -------
        list
        """
        nw = self._network()
        if nw is None:
            subsystems = self._walk() if recursive else [self]
            return [conn for sub in subsystems for conn in sub.conns.values()]
        connections = []
        for conn in nw.conns["object"]:
            owner = nw.subsystem_of(conn)
            if owner is self or (
                    recursive and owner is not None and self._contains(owner)):
                connections.append(conn)
        return connections

    @property
    def boundary_connections(self):
        """Connections of the network with exactly one end in the subsystem."""
        nw = self._network()
        if nw is None:
            msg = (
                f"The boundary connections of the subsystem {self.path} are "
                "known once it is part of a network."
            )
            logger.error(msg)
            raise TESPyComponentError(msg)
        members = set(self.get_components(recursive=True))
        return [
            conn for conn in nw.conns["object"]
            if (conn.source in members) != (conn.target in members)
        ]

    def get_conn(self, label):

        try:
            return self.conns[label]
        except KeyError:
            msg = (
                f"Connection with label {label} not found in subsystem "
                f"{self.label}. Use the label without the prefix of the "
                "subsystem."
            )
            logger.error(msg)
            raise KeyError(msg)

    def get_comp(self, label):

        try:
            return self.comps[label]
        except KeyError:
            msg = (
                f"Component with label {label} not found in subsystem "
                f"{self.label}. Use the label without the prefix of the "
                "subsystem."
            )
            logger.error(msg)
            raise KeyError(msg)

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
        from tespy.networks.network import _check_meta
        from tespy.networks.network import _construct_components
        from tespy.networks.network import _construct_connections

        _check_meta(subsystem_data, "subsystem")
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

        sub = cls.__new__(cls)
        sub._init_structure(label)
        for attr in [
                "num_in", "num_out", "num_heat_in", "num_heat_out",
                "num_power_in", "num_power_out"]:
            setattr(sub, attr, subsystem_data["Subsystem"][attr])
        sub.interface_map = subsystem_data["Subsystem"].get(
            "interface_map", {}
        )

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
            sub._add_interface_ports()
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
    def from_network(cls, label, nw, interface_exceptions=None, keep="specifications"):
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

        keep : str, optional
            Which state of the network to carry into the subsystem, default
            :code:`"specifications"`.

            - :code:`"specifications"`: all specifications and values are
              carried over unchanged.
            - :code:`"starting_values"`: all specifications are deactivated,
              the current values of the variables become starting values for
              the next solve.
            - :code:`"nothing"`: all specifications are deactivated and no
              values are carried over.

            Model data, e.g. characteristic lines, polynomial coefficients or
            reference definitions, as well as the design/offdesign setup stay
            on the objects in all modes, only their activation state is cleared

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

        With the :code:`keep` parameter the specifications can be deactivated
        instead, turning the network's values into starting values for the
        next solve, or dropping them altogether.

        >>> sub2 = Subsystem.from_network(
        ...     "heating block 2", nw, keep="starting_values"
        ... )
        >>> sub2.get_conn("1").T.is_set
        False
        >>> sub2.get_conn("1").p.val0
        1.0
        >>> sub3 = Subsystem.from_network("heating block 3", nw, keep="nothing")
        >>> sub3.get_conn("1").p.val0
        nan
        """
        if interface_exceptions is None:
            interface_exceptions = []
        if keep not in ("specifications", "starting_values", "nothing"):
            msg = (
                "The keep parameter must be one of 'specifications', "
                "'starting_values' or 'nothing'."
            )
            logger.error(msg)
            raise ValueError(msg)

        nw_dict = nw.export()

        if keep != "specifications":
            for section in ("Connection", "Component"):
                for class_data in nw_dict[section].values():
                    for obj_data in class_data.values():
                        for param_data in obj_data.values():
                            if not isinstance(param_data, dict):
                                continue
                            if "engine" in param_data:
                                param_data["is_set"] = []
                                if keep == "starting_values":
                                    param_data["val0"] = dict(param_data["val"])
                                continue
                            if "is_set" in param_data:
                                param_data["is_set"] = False
                            if "is_var" in param_data:
                                param_data["is_var"] = False
                            if "val" in param_data and "unit" in param_data:
                                if keep == "starting_values":
                                    param_data["val0"] = param_data.pop("val")
                                    param_data["val0_unit"] = param_data.pop("unit")
                                else:
                                    param_data.pop("val")
                                param_data.pop("val_SI", None)

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

        nw_dict["meta"]["kind"] = "subsystem"
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
