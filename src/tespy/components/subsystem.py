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
from tespy.tools.helpers import TESPyComponentError
from tespy.tools.units import Units
import tespy
import importlib
#from tespy.components.component import component_registry
#from tespy.connections.connection import connection_registry
#from tespy.tools import helpers as hlp


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

#def from_network(nw: tespy.networks.Network):
#    nw_dict=nw.export()
#    nw_import=tespy.networks.Network().from_dict(nw_dict)
    @classmethod
    def from_dict(cls, subsystem_data):
        # create network
        # get method to ensure compatibility with old style export
        units = Units.from_json(subsystem_data["Network"].get("units", {}))
        subsystem_data["Network"]["units"] = units
        nw = tespy.networks.Network(**subsystem_data["Network"])
        sub = cls.__new__(cls)
        sub.label = subsystem_data['Subsystem']['label']
        # load components
        comps = {}
        module_name = "tespy.components"
        _ = importlib.import_module(module_name)
        sub.num_in = subsystem_data['Subsystem']['num_in']
        sub.num_out = subsystem_data['Subsystem']['num_out']
        sub.num_heat_in = subsystem_data['Subsystem']['num_heat_in']
        sub.num_heat_out = subsystem_data['Subsystem']['num_heat_out']
        sub.num_power_in = subsystem_data['Subsystem']['num_power_in']
        sub.num_power_out = subsystem_data['Subsystem']['num_power_out']

        sub.conns={}
        sub.comps={}

        for component, data in subsystem_data["Component"].items():
            if component not in tespy.components.component.component_registry.items:
                msg = (
                    f"A class {component} is not available through the "
                    "tespy.components.component.component_registry decorator. "
                    "If you are using a custom component make sure to "
                    "decorate the class."
                )
                ##logger.error(msg)
                raise hlp.TESPyNetworkError(msg)

            target_class = tespy.components.component.component_registry.items[component]
            comps.update(tespy.networks.network._construct_components(target_class, data, nw)) 
        sub.inlet= comps[f'{sub.label}_inlet']
        sub.outlet= comps[f'{sub.label}_outlet']
        msg = 'Created network components.'
        ##logger.info(msg)

        conns = {}
        # load connections
        if "Connection" not in subsystem_data["Connection"]:
            # v0.8 compatibility
            target_class = tespy.connections.connection.connection_registry.items["Connection"]
            conns.update(tespy.networks.network._construct_connections(
                target_class, subsystem_data["Connection"], comps)
            )
        else:
            for connection, data in subsystem_data["Connection"].items():
                if connection not in tespy.connections.connection.connection_registry.items:
                    msg = (
                        f"A class {connection} is not available through the "
                        "tespy.connections.connection.connection_registry "
                        "decorator. If you are using a custom connection make "
                        "sure to decorate the class."
                    )
                    ##logger.error(msg)
                    raise hlp.TESPyNetworkError(msg)

                target_class = tespy.connections.connection.connection_registry.items[connection]
                conns.update(tespy.networks.network._construct_connections(target_class, data, comps))

        # add connections to network
        for c in conns.values():
            sub.add_conns(c)

        msg = 'Created connections.'
        logger.info(msg)

        msg = 'Created Subsystem.'
        logger.info(msg)


        return sub

    @classmethod
    def from_network(cls, label, nw, interface_exceptions=[]):


        nw_dict=nw.export()
        nw_dict['Component']['SubsystemInterface']={}
        nw_dict['Component']['SubsystemInterface'].update(
            SubsystemInterface(f'{label}_inlet')._serialize())
        nw_dict['Component']['SubsystemInterface'].update(
            SubsystemInterface(f'{label}_outlet')._serialize())
        

        nw_dict['Subsystem']={}
        for attr in ['num_in', 'num_out', 'num_heat_in', 'num_heat_out', 'num_power_in','num_power_out']:
            nw_dict['Subsystem'][attr]=0

        for typ, con in nw_dict['Connection'].get('Connection',{}).items():
            if con['source'] in nw_dict['Component'].get('Source', {}).keys() and con['source'] not in interface_exceptions:
                del nw_dict['Component']['Source'][con['source']]
                con['source']= f'{label}_inlet'
                nw_dict['Subsystem']['num_in'] +=1
                con['source_id']=f'out{nw_dict['Subsystem']['num_in']}'
                
            if con['target'] in nw_dict['Component'].get('Sink', {}).keys() and con['target'] not in interface_exceptions:
                del nw_dict['Component']['Sink'][con['target']]
                con['target']= f'{label}_outlet'
                nw_dict['Subsystem']['num_out']  +=1
                con['target_id']=f'in{nw_dict['Subsystem']['num_out'] }'

        for typ, con in nw_dict['Connection'].get('HeatConnection',{}).items():
            if con['source'] in nw_dict['Component'].get('HeatSource',{}).keys() and con['source'] not in interface_exceptions:
                del nw_dict['Component']['HeatSource'][con['source']]
                con['source']= f'{label}_inlet'
                nw_dict['Subsystem']['num_heat_in'] +=1
                con['source_id']=f'heat_out{nw_dict['Subsystem']['num_heat_in']}'
                
            if con['target'] in nw_dict['Component'].get('HeatSink',{}).keys() and con['target'] not in interface_exceptions:
                del nw_dict['Component']['HeatSink'][con['target']]
                con['target']= f'{label}_outlet'
                nw_dict['Subsystem']['num_heat_out']  +=1
                con['target_id']=f'heat_in{nw_dict['Subsystem']['num_heat_out'] }'
                
        for typ, con in nw_dict['Connection'].get('PowerConnection', {}).items():
            if con['source'] in nw_dict['Component'].get('PowerSource',{}).keys() and con['source'] not in interface_exceptions:
                del nw_dict['Component']['PowerSource'][con['source']]
                con['source']= f'{label}_inlet'
                nw_dict['Subsystem']['num_power_in'] +=1
                con['source_id']=f'power_out{nw_dict['Subsystem']['num_power_in']}'
            
            if con['target'] in nw_dict['Component'].get('PowerSink',{}).keys() and con['target'] not in interface_exceptions:
                del nw_dict['Component']['PowerSink'][con['target']]
                con['target']= f'{label}_outlet'
                nw_dict['Subsystem']['num_power_out']  +=1
                con['target_id']=f'power_in{nw_dict['Subsystem']['num_power_out'] }'
        
        nw_dict['Component']['SubsystemInterface'][f'{label}_outlet']['num_inter']={'val':nw_dict['Subsystem']['num_out'], 'is_set':True}
        nw_dict['Component']['SubsystemInterface'][f'{label}_inlet']['num_inter']={'val':nw_dict['Subsystem']['num_in'], 'is_set':True}
        nw_dict['Component']['SubsystemInterface'][f'{label}_outlet']['num_heat_inter']={'val':nw_dict['Subsystem']['num_heat_out'], 'is_set':True}
        nw_dict['Component']['SubsystemInterface'][f'{label}_inlet']['num_heat_inter']={'val':nw_dict['Subsystem']['num_heat_in'], 'is_set':True}
        nw_dict['Component']['SubsystemInterface'][f'{label}_outlet']['num_power_inter']={'val':nw_dict['Subsystem']['num_power_out'], 'is_set':True}
        nw_dict['Component']['SubsystemInterface'][f'{label}_inlet']['num_power_inter']={'val':nw_dict['Subsystem']['num_power_in'], 'is_set':True}
        nw_dict['Subsystem']['label']= label
        sub_network= cls.from_dict(nw_dict)

        
        return sub_network  