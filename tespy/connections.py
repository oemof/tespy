"""
.. module:: connections
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import numpy as np
import pandas as pd

from tespy.helpers import MyConnectionError, data_container, dc_prop, dc_flu
from tespy.components import components as cmp


class connection:
    """
    creates connection between two components

    - check argument consistency
    - set attributes to specified values

    :param comp1: connections source
    :type comp1: tespy.components.components.component
    :param outlet_id: outlet id at the connections source
    :type outlet_id: str
    :param comp2: connections target
    :type comp2: tespy.components.components.component
    :param inlet_id: inlet id at the connections target
    :type inlet_id: str
    :returns: no return value
    :raises: - :code:`TypeError`, if comp1 and comp2 are not of type
               components
             - :code:`ValueError`, if comp1 and comp2 are the same object
             - :code:`ValueError`, if outlet_id or inlet_id are not allowed
               for ids for comp1 or comp2

    **allowed keywords** in kwargs (also see connections.attr().keys()):

    - m (*numeric*, *ref object*), m0 (*numeric*)
    - p (*numeric*, *ref object*), p0 (*numeric*)
    - h (*numeric*, *ref object*), h0 (*numeric*)
    - T (*numeric*, *ref object*)
    - x (*numeric*)
    - fluid (*dict*), fluid_balance (*bool*)
    - design (*list*), offdesign (*list*)

    **(off-)design parameters**

    The specification of values for design and/or offdesign is used for
    automatic switch from design to offdesign calculation: All parameters given
    in 'design', e. g. :code:`design=['T', 'p']`, are unset in any offdesign
    calculation, parameters given in 'offdesign' are set for offdesign
    calculation.

    .. note::
        The fluid balance parameter applies a balancing of the fluid vector on
        the specified conntion to 100 %. For example, you have four fluid
        components (a, b, c and d) in your vector, you set two of them
        (a and b) and want the other two (components c and d) to be a result of
        your calculation. If you set this parameter to True, the equation
        (0 = 1 - a - b - c - d) will be applied.

    **example**

    .. code-block:: python

        from tespy import con
        conn = con.connection(turbine, 'out1', condenser, 'in1', m=10, p=0.05)

    creates connection from turbine to condenser (hot side inlet) and sets
    values for mass flow and pressure

    **improvements**

    """
    def __init__(self, comp1, outlet_id, comp2, inlet_id, **kwargs):

        # check input parameters
        if not (isinstance(comp1, cmp.component) and
                isinstance(comp2, cmp.component)):
            msg = ('Error creating connection.'
                   'Check if comp1, comp2 are of type component.')
            raise TypeError(msg)

        if comp1 == comp2:
            msg = ('Error creating connection. '
                   'Can\'t connect component to itself.')
            raise ValueError(msg)

        if outlet_id not in comp1.outlets():
            msg = ('Error creating connection. '
                   'Specified oulet_id is not valid for component ' +
                   comp1.component() + '. '
                   'Valid ids are: ' + str(comp1.outlets()) + '.')
            raise ValueError(msg)

        if inlet_id not in comp2.inlets():
            msg = ('Error creating connection. '
                   'Specified inlet_id is not valid for component ' +
                   comp2.component() + '. '
                   'Valid ids are: ' + str(comp2.inlets()) + '.')
            raise ValueError(msg)

        # set specified values
        self.s = comp1
        self.s_id = outlet_id
        self.t = comp2
        self.t_id = inlet_id

        self.design = []
        self.offdesign = []

        # set default values for kwargs
        var = self.attr()

        for key in self.attr().keys():
            self.__dict__.update({key: var[key]})

        self.set_attr(**kwargs)

    def set_attr(self, **kwargs):
        """
        sets, resets or unsets attributes of a connection, for the keyword
        arguments, return values and errors see object initialisation
        """
        var = self.attr()
        var0 = [x + '0' for x in var.keys()]

        # set specified values
        for key in kwargs:
            if key in var.keys() or key in var0:
                if (isinstance(kwargs[key], float) or
                        isinstance(kwargs[key], int)):
                    # unset
                    if np.isnan(kwargs[key]) and key not in var0:
                        self.get_attr(key).set_attr(val_set=False)
                        self.get_attr(key).set_attr(ref_set=False)
                    # starting value
                    elif key in var0:
                        self.get_attr(key.replace('0', '')).set_attr(
                                val0=kwargs[key])
                    # set/reset
                    else:
                        self.get_attr(key).set_attr(val_set=True,
                                                    val=kwargs[key],
                                                    val0=kwargs[key])

                # reference object
                elif isinstance(kwargs[key], ref):
                    if key == 'fluid' or key == 'x':
                        print('References for fluid vector and vapour mass '
                              'fraction not implemented.')
                    else:
                        self.get_attr(key).set_attr(ref=kwargs[key])
                        self.get_attr(key).set_attr(ref_set=True)

                # fluid specification
                elif isinstance(kwargs[key], dict):
                    # starting values
                    if key in var0:
                        self.get_attr(key).set_attr(val0=kwargs[key])
                    # specified parameters
                    else:
                        self.get_attr(key).set_attr(val=kwargs[key].copy())
                        self.get_attr(key).set_attr(val0=kwargs[key].copy())
                        for f in kwargs[key]:
                            kwargs[key][f] = True
                        self.get_attr(key).set_attr(val_set=kwargs[key])

                # data container specification
                elif isinstance(kwargs[key], data_container):
                    self.__dict__.update({key: kwargs[key]})

                # invalid datatype for keyword
                else:
                    msg = 'Bad datatype for keyword argument ' + str(key)
                    raise TypeError(msg)

            # fluid balance
            elif key == 'fluid_balance':
                if isinstance(kwargs[key], bool):
                    self.get_attr('fluid').set_attr(balance=kwargs[key])
                else:
                    msg = ('Datatype for keyword argument ' + str(key) +
                           ' must be bool.')
                    raise TypeError(msg)

            elif key == 'design' or key == 'offdesign':
                if not isinstance(kwargs[key], list):
                    msg = 'Please provide the design parameters as list!'
                    raise TypeError(msg)
                if set(kwargs[key]).issubset(var.keys()):
                    self.__dict__.update({key: kwargs[key]})
                else:
                    msg = ('Available parameters for (off-)design'
                           'specification are: ' + str(var.keys()) + '.')
                    raise ValueError(msg)

            # invalid keyword
            else:
                msg = 'Connection has no attribute ' + str(key)
                raise ValueError(msg)

    def get_attr(self, key):
        """
        get the value of a connection attribute

        :param key: attribute to return its value
        :type key: str
        :returns:
            - :code:`self.__dict__[key]` if object has attribute key
            - :code:`None` if object has no attribute key
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            print(self.connection(), ' has no attribute \"', key, '\"')
            return None

    def attr(self):
        """
        get the list of attributes allowed for a connection object

        :returns: list object
        """
        return {'m': dc_prop(), 'p': dc_prop(), 'h': dc_prop(), 'T': dc_prop(),
                'x': dc_prop(),
                'fluid': dc_flu()}

    def to_flow(self):
        """
        create a list containing the connections fluid information

        :returns: :code:`[mass flow, pressure, enthalpy, fluid vector]`
        """
        return [self.m.val_SI, self.p.val_SI, self.h.val_SI, self.fluid.val]


class bus:
    """
    establish power connection between turbines, pumps, heat exchanger

    **TODO**

    - add documentation

    **Improvements**

    - improve architecture (e. g. make it similar to connections)
    """

    def __init__(self, label, **kwargs):

        self.comps = pd.DataFrame(columns=['factor'])

        self.label = label
        self.P_set = False
        self.P = kwargs.get('P', np.nan)

        if not np.isnan(self.P):
            self.P_set = True

    def set_attr(self, **kwargs):

        self.label = kwargs.get('label')
        self.P = kwargs.get('P')

        if not np.isnan(self.P):
            self.P_set = True
        else:
            self.P_set = False

    def get_attr(self, key):
        """
        get the value of a bus attribute

        :param key: attribute to return its value
        :type key: str
        :returns:
            - :code:`self.__dict__[key]` if object has attribute key
            - :code:`None` if object has no attribute key
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            print(self.bus(), ' has no attribute \"', key, '\"')
            return None

    def add_comps(self, *args):
        """
        add component to bus
        """
        for c in args:
            if isinstance(c, list):
                if len(c) == 2:
                    self.comps.loc[c[0]] = [c[1]]
                    if not self.check_comp(c[0]):
                        self.comps = self.comps[:-1]
                else:
                    msg = 'List must have two elements: [component, factor].'
                    raise MyConnectionError(msg)
            else:
                msg = 'Provide argument as list: [component, factor].'
                raise MyConnectionError(msg)

    def check_comp(self, comp):
        """
        check component

        **TODO**

        - add documentation
        """
        for c in self.comps.index:
            if type(comp) != type(c):
                if (type(comp).__bases__[0] == type(c).__bases__[0] and
                        type(comp).__bases__[0] == cmp.component):

                    if type(c).__bases__[0] == cmp.component:
                        msg = ('Error adding component to power bus. '
                               'This bus accepts components of type ' +
                               str(type(c)) + '.')
                        raise TypeError(msg)
                    else:
                        msg = ('Error adding component to power bus. '
                               'This bus accepts components of type ' +
                               str(type(c).__bases__[0]) + '.')
                        raise TypeError(msg)
                return True
        return True


class ref:
    r"""
    creates reference object for network parametetrisation

    :param ref_obj: connection to be referenced
    :type ref_obj: tespy.connections.connection
    :param factor: factor for the reference
    :type factor: numeric
    :param delta: delta for the reference
    :type delta: numeric
    :returns: no return value

    **example**

    .. math::
        \dot{m} = \dot{m}_\mathrm{ref} \cdot factor + delta\\
        m_{ref}: \text{mass flow of referenced connection}

    """
    def __init__(self, ref_obj, factor, delta):
        if not isinstance(ref_obj, connection):
            msg = 'First parameter must be object of type connection.'
            raise TypeError(msg)

        if not (isinstance(factor, int) or isinstance(factor, float)):
            msg = 'Second parameter must be of type int or float.'
            raise TypeError(msg)

        if not (isinstance(delta, int) or isinstance(delta, float)):
            msg = 'Thrid parameter must be of type int or float.'
            raise TypeError(msg)

        self.obj = ref_obj
        self.f = factor
        self.d = delta

    def get_attr(self, key):
        """
        get the value of a ref attribute

        :param key: attribute to return its value
        :type key: str
        :returns:
            - :code:`self.__dict__[key]` if object has attribute key
            - :code:`None` if object has no attribute key
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            print(self.bus(), ' has no attribute \"', key, '\"')
            return None
