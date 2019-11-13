# -*- coding: utf-8

"""This module contains the data_container class and its children. The
different ata_containers are used to store component or connection properties.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/helpers.py

SPDX-License-Identifier: MIT
"""

import collections
import logging

# %%


class data_container:
    r"""
    Class data_container is the base class for dc_cc, dc_cp, dc_flu, dc_prop,
    dc_simple.

    Parameters
    ----------
    **kwargs :
        See the class documentation of desired data_container for available
        keywords.

    Note
    ----
    The initialisation method (__init__), setter method (set_attr) and getter
    method (get_attr) are used for instances of class data_container and its
    children. TESPy uses different data_containers for specific tasks:
    Component characteristics (dc_cc), component maps (dc_cm), component
    properties (dc_cp), grouped component properites (dc_gcp), fluid
    composition (dc_flu), fluid properties (dc_prop).

    Grouped component properties are used, if more than one component property
    has to be specified in order to apply one equation, e. g. pressure drop in
    pipes by specified length, diameter and roughness. If you specify all three
    of these properties, the data_container for the group will be created
    automatically!

    For the full list of available parameters for each data container, see its
    documentation.

    Example
    -------
    >>> from tespy import hlp, cmp
    >>> type(hlp.dc_cm(is_set=True))
    <class 'tespy.tools.helpers.dc_cm'>
    >>> type(hlp.dc_cc(x=[1, 2, 3, 4], y=[1, 4, 9, 16], is_set=True))
    <class 'tespy.tools.helpers.dc_cc'>
    >>> type(hlp.dc_cp(val=100, is_set=True, is_var=True, printout=True,
    ...     max_val=1000, min_val=1))
    <class 'tespy.tools.helpers.dc_cp'>
    >>> pipe = cmp.pipe('testpipe', L=100, D=0.5, ks=5e-5)
    >>> type(hlp.dc_gcp(is_set=True, elements=[pipe.L, pipe.D, pipe.ks],
    ...     method='default'))
    <class 'tespy.tools.helpers.dc_gcp'>
    >>> type(hlp.dc_flu(val={'CO2': 0.1, 'H2O': 0.11, 'N2': 0.75, 'O2': 0.03},
    ...     val_set={'CO2': False, 'H2O': False, 'N2': False, 'O2': True},
    ...     balance=False))
    <class 'tespy.tools.helpers.dc_flu'>
    >>> type(hlp.dc_prop(val=5, val_SI=500000, val_set=True, unit='bar',
    ...     unit_set=False, ref=None, ref_set=False))
    <class 'tespy.tools.helpers.dc_prop'>
    >>> type(hlp.dc_simple(val=5, val_set=False))
    <class 'tespy.tools.helpers.dc_simple'>
    """

    def __init__(self, **kwargs):

        var = self.attr()

        # default values
        for key in var.keys():
            self.__dict__.update({key: var[key]})

        self.set_attr(**kwargs)

    def set_attr(self, **kwargs):
        r"""
        Sets, resets or unsets attributes of a data_container type object.

        Parameters
        ----------
        **kwargs :
            See the class documentation of desired data_container for available
            keywords.
        """
        var = self.attr()
        # specify values
        for key in kwargs:
            if key in var.keys():
                self.__dict__.update({key: kwargs[key]})

    def get_attr(self, key):
        r"""
        Get the value of a data_container's attribute.

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
            msg = ('Datacontainer of type ' + self.__class__.__name__ +
                   ' has no attribute \"' + str(key) + '\".')
            logging.error(msg)
            raise KeyError(msg)

    def attr(self):
        r"""
        Return the available attributes for a data_container type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {}

# %%


class dc_prop(data_container):
    r"""
    Parameters
    ----------
    val : float
        Value in user specified unit (or network unit) if unit is unspecified,
        default: val=np.nan.

    val0 : float
        Starting value in user specified unit (or network unit) if unit is
        unspecified, default: val0=np.nan.

    val_SI : float
        Value in SI_unit, default: val_SI=0.

    val_set : boolean
        Has the value for this property been set? default: val_set=False.

    ref : tespy.connections.ref
        Reference object, default: ref=None.

    ref_set : boolean
        Has a value for this property been referenced to another connection?
        default: ref_set=False.

    unit : str
        Unit for this property, default: ref=None.

    unit : boolean
        Has the unit for this property been specified manually by the user?
        default: unit_set=False.
    """
    def attr(self):
        return {'val': np.nan, 'val0': np.nan, 'val_SI': 0, 'val_set': False,
                'ref': None, 'ref_set': False,
                'unit': None, 'unit_set': False, 'design': np.nan}

# %%


class dc_simple(data_container):
    r"""
    Simple data container without data type restrictions to val field.

    Parameters
    ----------
    val : no specific datatype
        Value for the property, no predefined datatype. Unset this property by
        stating val=np.nan.

    val_set : boolean
        Has the value for this property been set? default: val_set=False.
    """
    def attr(self):
        return {'val': np.nan, 'val_set': False}

# %%


class dc_flu(data_container):
    r"""
    Parameters
    ----------

    val : dict
        Mass fractions of the fluids in a mixture, default: val={}.
        Pattern for dictionary: keys are fluid name, values are mass fractions.

    val0 : dict
        Starting values for mass fractions of the fluids in a mixture,
        default: val0={}. Pattern for dictionary: keys are fluid name, values
        are mass fractions.

    val_set : dict
        Which fluid mass fractions have been set, default val_set={}.
        Pattern for dictionary: keys are fluid name, values are True or False.

    balance : boolean
        Should the fluid balance equation be applied for this mixture?
        default: False.
    """
    def attr(self):
        r"""
        Return the available attributes for a data_container type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {'val': {}, 'val0': {}, 'val_set': {},
                'design': collections.OrderedDict(), 'balance': False}

# %%


class dc_cp(data_container):
    r"""
    Parameters
    ----------

    val : float
        Value for this component attribute, default: val=1.

    val_SI : float
        Value in SI_unit (available for temperatures only, unit transformation
        according to network's temperature unit), default: val_SI=0.

    is_set : boolean
        Has the value for this attribute been set?, default: is_set=False.

    is_var : boolean
        Is this attribute part of the system variables?, default: is_var=False.

    d : float
        Interval width for numerical calculation of partial derivative towards
        this attribute, it is part of the system variables, default d=1e-4.

    min_val : float
        Minimum value for this attribute, used if attribute is part of the
        system variables, default: min_val=1.1e-4.

    max_val : float
        Maximum value for this attribute, used if attribute is part of the
        system variables, default: max_val=1e12.

    printout : boolean
        Should the value of this attribute be printed in the results overview?
    """
    def attr(self):
        r"""
        Return the available attributes for a data_container type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {'val': 1, 'val_SI': 0, 'is_set': False, 'printout': True,
                'd': 1e-4, 'min_val': -1e12, 'max_val': 1e12, 'is_var': False,
                'val_ref': 1, 'design': np.nan}

# %%


class dc_cc(data_container):
    r"""
    Parameters
    ----------

    func : tespy.components.characteristics.characteristics
        Function to be applied for this characteristics, default: None.

    is_set : boolean
        Should this equation be applied?, default: is_set=False.

    method : str
        Which default method for this characteristic function should be used?
        default: method='default'.

    param : str
        Which parameter should be applied as the x value?
        default: method='default'.

    x : numpy.array
        Array for the x-values of the characteristic line, default x=None.

    y : numpy.array
        Array for the y-values of the characteristic line, default y=None.

    Note
    ----
    If you do not specify x-values or y-values, default values according to the
    specified method will be used. If you specify a method as well as x-values
    and/or y-values, these will override the defaults values of the chosen
    method.
    """
    def attr(self):
        r"""
        Return the available attributes for a data_container type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {'func': None, 'is_set': False,
                'method': 'default', 'param': None,
                'x': None, 'y': None}

# %%


class dc_cm(data_container):
    r"""
    Parameters
    ----------

    func : tespy.components.characteristics.characteristics
        Function to be applied for this characteristic map, default: None.

    is_set : boolean
        Should this equation be applied?, default: is_set=False.

    method : str
        Which default method for this characteristic function should be used?
        default: method='default'.

    param : str
        Which parameter should be applied as the x value?
        default: method='default'.

    x : numpy.array
        Array for the x-values of the characteristic line, default x=None.

    y : numpy.array
        Array for the y-values of the characteristic line, default y=None.

    z1 : numpy.array
        Array for the y-values of the characteristic line, default y=None.

    z2 : numpy.array
        Array for the y-values of the characteristic line, default y=None.

    Note
    ----
    If you do not specify any interpolation points (x, y, z1, z2), default
    values according to the specified method will be used. If you specify a
    method as well as interpolation points, these will override the defaults
    values of the chosen method.
    """
    def attr(self):
        r"""
        Return the available attributes for a data_container type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {'func': None, 'is_set': False,
                'method': 'default', 'param': None,
                'x': None, 'y': None, 'z1': None, 'z2': None}

# %%


class dc_gcp(data_container):
    r"""
    Parameters
    ----------
    is_set : boolean
        Should the equation for this parameter group be applied?
        default: is_set=False.

    method : str
        Which calculation method for this parameter group should be used?
        default: method='default'.

    elements : list
        Which component properties are part of this component group?
        default elements=[].
    """
    def attr(self):
        r"""
        Return the available attributes for a data_container type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {'is_set': False, 'method': 'default', 'elements': []}
