# -*- coding: utf-8

"""Module of class Connection and class Ref.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/connections/connection.py

SPDX-License-Identifier: MIT
"""

import logging

import numpy as np

from tespy.components.component import Component
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.data_containers import FluidComposition as dc_flu
from tespy.tools.data_containers import FluidProperties as dc_prop
from tespy.tools.helpers import TESPyConnectionError

# pass the warning messages to the logger
logging.captureWarnings(True)


class Connection:
    r"""
    Class connection is the container for fluid properties between components.

    Parameters
    ----------
    m : float, tespy.connections.connection.Ref, tespy.tools.data_containers.FluidProperties
        Mass flow specification.

    m0 : float
        Starting value specification for mass flow.

    p : float, tespy.connections.connection.Ref, tespy.tools.data_containers.FluidProperties
        Pressure specification.

    p0 : float
        Starting value specification for pressure.

    h : float, tespy.connections.connection.Ref, tespy.tools.data_containers.FluidProperties
        Enthalpy specification.

    h0 : float
        Starting value specification for enthalpy.

    fluid : dict, tespy.tools.data_containers.FluidComposition
        Fluid compostition specification.

    fluid0 : dict
        Starting value specification for fluid compostition.

    fluid_balance : boolean
        Fluid balance equation specification.

    x : float, tespy.tools.data_containers.FluidProperties
        Gas phase mass fraction specification.

    T : float, tespy.connections.connection.Ref, tespy.tools.data_containers.FluidProperties
        Temperature specification.

    Td_bp : float, tespy.tools.data_containers.FluidProperties
        Temperature difference to boiling point at pressure corresponding
        pressure of this connection in K.

    v : float, tespy.tools.data_containers.FluidProperties
        Volumetric flow specification.

    state : str
        State of the pure fluid on this connection: liquid ('l') or gaseous
        ('g').

    design : list
        List containing design parameters (stated as string).

    offdesign : list
        List containing offdesign parameters (stated as string).

    design_path : str
        Path to individual design case for this connection.

    local_offdesign : boolean
        Treat this connection in offdesign mode in a design calculation.

    local_design : boolean
        Treat this connection in design mode in an offdesign calculation.

    printout: boolean
        Include this connection in the network's results printout.

    label : str
        Label of the connection. The default value is:
        :code:`'source:source_id_target:target_id'`.

    Note
    ----
    - The fluid balance parameter applies a balancing of the fluid vector on
      the specified conntion to 100 %. For example, you have four fluid
      components (a, b, c and d) in your vector, you set two of them
      (a and b) and want the other two (components c and d) to be a result of
      your calculation. If you set this parameter to True, the equation
      (0 = 1 - a - b - c - d) will be applied.

    - The specification of values for design and/or offdesign is used for
      automatic switch from design to offdesign calculation: All parameters
      given in 'design', e.g. :code:`design=['T', 'p']`, are unset in any
      offdesign calculation, parameters given in 'offdesign' are set for
      offdesign calculation.

    Example
    -------
    This example shows how to create connections and specify parameters. First
    create the required components and connect them in the next step. After
    that, it is possible specify parameters with the :code:`set_attr` method.

    >>> from tespy.components import Sink, Source
    >>> from tespy.connections import Connection, Ref
    >>> from tespy.tools import FluidComposition as dc_flu
    >>> from tespy.tools import FluidProperties as dc_prop
    >>> import numpy as np
    >>> so1 = Source('source1')
    >>> so2 = Source('source2')
    >>> si1 = Sink('sink1')
    >>> si2 = Sink('sink2')
    >>> so_si1 = Connection(so1, 'out1', si1, 'in1')
    >>> so_si2 = Connection(so2, 'out1', si2, 'in1')

    There are different ways of setting parameters on connections: Specify

    - a numeric value  (for attributes mass flow, pressure and enthalpy)
    - a numeric starting value (for attributes mass flow, pressure and
      enthalpy)
    - a dictionary (for attributes fluid and fluid0)
    - a boolean value (for attributes fluid_balance, local_design,
      local_offdesign).
    - a referenced value (mass flow, pressure, temperature, enthalpy).
    - a data_container (for attributes fluid and fluid0 dc_flu, for other
      fluid propertie attributes dc_prop).
    - numpy.nan or None (unsetting a value).
    - a string (for attributes design_paht and state).
    - a list (for attributes design and offdesign).

    >>> so_si1.set_attr(v=0.012, m0=10, p=5, h=400, fluid={'H2O': 1, 'N2': 0})

    Individual specification of pressure with a data container. Setting the
    unit individually will overwrite the network's unit specification for this
    connection.

    >>> p = dc_prop(val=50, val_set=True, unit='bar', unit_set=True)
    >>> so_si2.set_attr(m=Ref(so_si1, 2, -5), p=p, h0=700, T=200,
    ... fluid={'N2': 1}, fluid_balance=True)

    The set_attr method automatically converts your input in data_container
    information.

    >>> type(so_si1.v)
    <class 'tespy.tools.data_containers.FluidProperties'>
    >>> type(so_si1.fluid)
    <class 'tespy.tools.data_containers.FluidComposition'>

    If you want get a spcific value use the logic: connection.property.*.
    Aditionally, it is possible to use the :code:`get_attr` method.

    >>> so_si1.m.val0
    10
    >>> so_si1.m.get_attr('val_set')
    False
    >>> type(so_si2.m.ref)
    <class 'tespy.connections.connection.Ref'>
    >>> so_si2.fluid.get_attr('balance')
    True
    >>> so_si2.m.ref.get_attr('d')
    -5
    >>> so_si2.m.ref_set
    True
    >>> type(so_si2.m.ref.get_attr('obj'))
    <class 'tespy.connections.connection.Connection'>

    Unset the specified temperature and specify temperature difference to
    boiling point instead.

    >>> so_si2.T.val_set
    True
    >>> so_si2.set_attr(Td_bp=5, T=np.nan)
    >>> so_si2.T.val_set
    False
    >>> so_si2.Td_bp.val
    5
    >>> so_si2.set_attr(Td_bp=None)
    >>> so_si2.Td_bp.val_set
    False

    Specify the state keyword: The fluid will be forced to liquid or gaseous
    state in this case.

    >>> so_si2.set_attr(state='l')
    >>> so_si2.state.is_set
    True
    >>> so_si2.set_attr(state=np.nan)
    >>> so_si2.state.is_set
    False
    >>> so_si2.set_attr(state='g')
    >>> so_si2.state.is_set
    True
    >>> so_si2.set_attr(state=None)
    >>> so_si2.state.is_set
    False
    >>> so_si2.set_attr(label='myconnection')
    >>> so_si2.label
    'myconnection'
    """

    def __init__(self, comp1, outlet_id, comp2, inlet_id, **kwargs):

        # check input parameters
        if not (isinstance(comp1, Component) and
                isinstance(comp2, Component)):
            msg = ('Error creating connection. Check if comp1, comp2 are of '
                   'type component.')
            logging.error(msg)
            raise TypeError(msg)

        if comp1 == comp2:
            msg = ('Error creating connection. Cannot connect component ' +
                   comp1.label + ' to itself.')
            logging.error(msg)
            raise TESPyConnectionError(msg)

        if outlet_id not in comp1.outlets():
            msg = ('Error creating connection. Specified oulet_id (' +
                   outlet_id + ') is not valid for component ' +
                   comp1.component() + '. Valid ids are: ' +
                   str(comp1.outlets()) + '.')
            logging.error(msg)
            raise ValueError(msg)

        if inlet_id not in comp2.inlets():
            msg = (
                'Error creating connection. Specified inlet_id (' + inlet_id +
                ') is not valid for component ' + comp2.component() +
                '. Valid ids are: ' + str(comp2.inlets()) + '.')
            logging.error(msg)
            raise ValueError(msg)

        self.label = (
            comp1.label + ':' + outlet_id + '_' + comp2.label + ':' + inlet_id)

        # set specified values
        self.source = comp1
        self.source_id = outlet_id
        self.target = comp2
        self.target_id = inlet_id

        # defaults
        self.new_design = True
        self.design_path = None
        self.design = []
        self.offdesign = []
        self.local_design = False
        self.local_offdesign = False
        self.printout = True

        # set default values for kwargs
        self.variables = self.attr()
        self.variables0 = [x + '0' for x in self.variables.keys()]
        self.__dict__.update(self.variables)
        self.set_attr(**kwargs)

        msg = (
            'Created connection ' + self.source.label + ' (' + self.source_id +
            ') -> ' + self.target.label + ' (' + self.target_id + ').')
        logging.debug(msg)

    def set_attr(self, **kwargs):
        r"""
        Set, reset or unset attributes of a connection.

        Parameters
        ----------
        m : float, tespy.connections.connection.Ref, tespy.tools.data_containers.FluidProperties
            Mass flow specification.

        m0 : float
            Starting value specification for mass flow.

        p : float, tespy.connections.connection.Ref, tespy.tools.data_containers.FluidProperties
            Pressure specification.

        p0 : float
            Starting value specification for pressure.

        h : float, tespy.connections.connection.Ref, tespy.tools.data_containers.FluidProperties
            Enthalpy specification.

        h0 : float
            Starting value specification for enthalpy.

        fluid : dict, tespy.tools.data_containers.FluidComposition
            Fluid compostition specification.

        fluid0 : dict
            Starting value specification for fluid compostition.

        fluid_balance : boolean
            Fluid balance equation specification.

        x : float, tespy.tools.data_containers.FluidProperties
            Gas phase mass fraction specification.

        T : float, tespy.connections.connection.Ref, tespy.tools.data_containers.FluidProperties
            Temperature specification.

        Td_bp : float, tespy.tools.data_containers.FluidProperties
            Temperature difference to boiling point at pressure corresponding
            pressure of this connection in K.

        v : float, tespy.tools.data_containers.FluidProperties
            Volumetric flow specification.

        state : str
            State of the pure fluid on this connection: liquid ('l') or gaseous
            ('g').

        design : list
            List containing design parameters (stated as string).

        offdesign : list
            List containing offdesign parameters (stated as string).

        design_path : str
            Path to individual design case for this connection.

        local_offdesign : boolean
            Treat this connection in offdesign mode in a design calculation.

        local_design : boolean
            Treat this connection in design mode in an offdesign calculation.

        printout: boolean
            Include this connection in the network's results printout.

        label : str
            Label of the connection. The default value is:
            :code:`'source:source_id_target:target_id'`.

        Note
        ----
        - The fluid balance parameter applies a balancing of the fluid vector
          on the specified conntion to 100 %. For example, you have four fluid
          components (a, b, c and d) in your vector, you set two of them
          (a and b) and want the other two (components c and d) to be a result
          of your calculation. If you set this parameter to True, the equation
          (0 = 1 - a - b - c - d) will be applied.
        - The specification of values for design and/or offdesign is used for
          automatic switch from design to offdesign calculation: All parameters
          given in 'design', e.g. :code:`design=['T', 'p']`, are unset in any
          offdesign calculation, parameters given in 'offdesign' are set for
          offdesign calculation.
        - The property state is applied on pure fluids only. If you specify the
          desired state of the fluid at a connection the convergence check will
          adjust the enthalpy values of that connection for the first
          iterations in order to meet the state requirement.
        """
        # set specified values
        for key in kwargs:
            if key in self.variables.keys() or key in self.variables0:
                # fluid specification
                try:
                    float(kwargs[key])
                    is_numeric = True
                except (TypeError, ValueError):
                    is_numeric = False
                if 'fluid' in key:
                    if isinstance(kwargs[key], dict):
                        # starting values
                        if key in self.variables0:
                            self.fluid.set_attr(val0=kwargs[key].copy())
                        # specified parameters
                        else:
                            self.fluid.set_attr(val=kwargs[key].copy())
                            self.fluid.set_attr(
                                val_set={f: True for f in kwargs[key].keys()})

                    elif isinstance(kwargs[key], dc_flu):
                        # data container for fluids
                        self.__dict__.update({key: kwargs[key]})

                    else:
                        # bad datatype
                        msg = (
                            'Datatype for fluid vector specification must be '
                            'tespy.tools.data_containers.FluidComposition or dict.')
                        logging.error(msg)
                        raise TypeError(msg)

                elif key == 'state':
                    if kwargs[key] in ['l', 'g']:
                        self.state.set_attr(val=kwargs[key], is_set=True)
                    elif isinstance(kwargs[key], dc_simple):
                        self.state = kwargs[key]
                    elif kwargs[key] is None:
                        self.state.set_attr(is_set=False)
                    elif is_numeric:
                        if np.isnan(kwargs[key]):
                            self.get_attr(key).set_attr(is_set=False)
                        else:
                            msg = (
                                'To unset the state specification either use '
                                'np.nan or None.')
                            logging.error(msg)
                            raise ValueError(msg)
                    else:
                        msg = (
                            'Keyword argument "state" must either be '
                            '"l" or "g" or be None or np.nan.')
                        logging.error(msg)
                        raise TypeError(msg)

                elif kwargs[key] is None:
                    self.get_attr(key).set_attr(val_set=False)
                    self.get_attr(key).set_attr(ref_set=False)

                elif is_numeric:
                    if np.isnan(kwargs[key]):
                        self.get_attr(key).set_attr(val_set=False)
                        self.get_attr(key).set_attr(ref_set=False)
                    else:
                        # value specification
                        if key in self.variables:
                            self.get_attr(key).set_attr(
                                val_set=True,
                                val=kwargs[key])
                        # starting value specification
                        else:
                            self.get_attr(key.replace('0', '')).set_attr(
                                val0=kwargs[key])

                # reference object
                elif isinstance(kwargs[key], Ref):
                    if key in ['x', 'v', 'Td_bp']:
                        msg = (
                            'References for volumetric flow, vapor mass '
                            'fraction and subcooling/superheating are not '
                            'implemented.')
                        logging.error(msg)
                        raise NotImplementedError(msg)
                    else:
                        self.get_attr(key).set_attr(ref=kwargs[key])
                        self.get_attr(key).set_attr(ref_set=True)

                # data container specification
                elif isinstance(kwargs[key], dc_prop):
                    self.__dict__.update({key: kwargs[key]})

                # invalid datatype for keyword
                else:
                    msg = 'Bad datatype for keyword argument ' + key + '.'
                    logging.error(msg)
                    raise TypeError(msg)

            # fluid balance
            elif key == 'fluid_balance':
                if isinstance(kwargs[key], bool):
                    self.get_attr('fluid').set_attr(balance=kwargs[key])
                else:
                    msg = (
                        'Datatype for keyword argument fluid_balance must be '
                        'boolean.')
                    logging.error(msg)
                    raise TypeError(msg)

            # design/offdesign parameter list
            elif key == 'design' or key == 'offdesign':
                if not isinstance(kwargs[key], list):
                    msg = 'Please provide the ' + key + ' parameters as list!'
                    logging.error(msg)
                    raise TypeError(msg)
                elif set(kwargs[key]).issubset(self.variables.keys()):
                    self.__dict__.update({key: kwargs[key]})
                else:
                    msg = (
                        'Available parameters for (off-)design specification '
                        'are: ' + str(self.variables.keys()) + '.')
                    logging.error(msg)
                    raise ValueError(msg)

            # design path
            elif key == 'design_path':
                if isinstance(kwargs[key], str):
                    self.__dict__.update({key: kwargs[key]})
                elif np.isnan(kwargs[key]):
                    self.design_path = None
                else:
                    msg = (
                        'Please provide the design_path parameter as string '
                        'or as nan.')
                    logging.error(msg)
                    raise TypeError(msg)

                self.new_design = True

            # other boolean keywords
            elif key in ['printout', 'local_design', 'local_offdesign']:
                if not isinstance(kwargs[key], bool):
                    msg = ('Please provide the ' + key + ' as boolean.')
                    logging.error(msg)
                    raise TypeError(msg)
                else:
                    self.__dict__.update({key: kwargs[key]})

            # label
            elif key == 'label':
                if isinstance(kwargs[key], str):
                    self.__dict__.update({key: kwargs[key]})
                else:
                    msg = 'Please provide the label as string.'
                    logging.error(msg)
                    raise TypeError(msg)

            # invalid keyword
            else:
                msg = 'Connection has no attribute ' + key + '.'
                logging.error(msg)
                raise KeyError(msg)

    def get_attr(self, key):
        r"""
        Get the value of a connection's attribute.

        Parameters
        ----------
        key : str
            The attribute you want to retrieve.

        Returns
        -------
        out :
            Specified attribute.
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            msg = 'Connection has no attribute \"' + key + '\".'
            logging.error(msg)
            raise KeyError(msg)

    @staticmethod
    def attr():
        r"""
        Return available attributes of a connection.

        Returns
        -------
        out : list
            List of available attributes of a connection.
        """
        return {'m': dc_prop(), 'p': dc_prop(), 'h': dc_prop(), 'T': dc_prop(),
                'x': dc_prop(), 'v': dc_prop(), 'vol': dc_prop(),
                's': dc_prop(),
                'fluid': dc_flu(), 'Td_bp': dc_prop(), 'state': dc_simple()}

    def to_flow(self):
        r"""
        Return the SI-values for the network variables.

        Returns
        -------
        out : list
            List of mass flow and fluid property information.
        """
        return [self.m.val_SI, self.p.val_SI, self.h.val_SI, self.fluid.val]

    def to_flow_design(self):
        r"""
        Return the SI-values for the network variables at design point.

        Returns
        -------
        out : list
            List of mass flow and fluid property information.
        """
        return [self.m.design, self.p.design, self.h.design, self.fluid.design]


class Ref:
    r"""
    Reference fluid properties from one connection to another connection.

    For example, reference the mass flow of one connection :math:`\dot{m}` to
    another mass flow :math:`\dot{m}_{ref}`:

    .. math::

        \dot{m} = \dot{m}_\mathrm{ref} \cdot f + d

    Parameters
    ----------
    obj : tespy.connections.connection.Connection
        Connection to be referenced.

    f : float
        Factor to multiply specified property with.

    d : float
        Delta to add after multiplication.
    """

    def __init__(self, ref_obj, factor, delta):

        if not isinstance(ref_obj, Connection):
            msg = 'First parameter must be object of type connection.'
            logging.error(msg)
            raise TypeError(msg)

        if not (isinstance(factor, int) or isinstance(factor, float)):
            msg = 'Second parameter must be of type int or float.'
            logging.error(msg)
            raise TypeError(msg)

        if not (isinstance(delta, int) or isinstance(delta, float)):
            msg = 'Thrid parameter must be of type int or float.'
            logging.error(msg)
            raise TypeError(msg)

        self.obj = ref_obj
        self.f = factor
        self.d = delta

        msg = ('Created reference object with factor ' + str(self.f) +
               ' and delta ' + str(self.d) + ' referring to connection ' +
               ref_obj.source.label + ' (' + ref_obj.source_id + ') -> ' +
               ref_obj.target.label + ' (' + ref_obj.target_id + ').')
        logging.debug(msg)

    def get_attr(self, key):
        r"""
        Get the value of a reference attribute.

        Parameters
        ----------
        key : str
            The attribute you want to retrieve.

        Returns
        -------
        out :
            Specified attribute.
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            msg = 'Reference has no attribute \"' + key + '\".'
            logging.error(msg)
            raise KeyError(msg)
