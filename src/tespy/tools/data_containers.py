# -*- coding: utf-8

"""Module for data container classes.

The DataContainer class and its subclasses are used to store component or
connection properties.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/data_containers.py

SPDX-License-Identifier: MIT
"""

import collections
import logging

import numpy as np

# %%


class DataContainer:
    """
    The DataContainer is parent class for all data containers.

    Parameters
    ----------
    **kwargs :
        See the class documentation of desired DataContainer for available
        keywords.

    Note
    ----
    The initialisation method (:code:`__init__`), setter method
    (:code:`set_attr`) and getter method (:code:`get_attr`) are used for
    instances of class DataContainer and its children. TESPy uses different
    :code:`DataContainer` classes for specific objectives:

    - component characteristics
      :py:class:`tespy.tools.data_containers.ComponentCharacteristics`
    - component characteristic maps
      :py:class:`tespy.tools.data_containers.ComponentCharacteristicMaps`
    - component properties
      :py:class:`tespy.tools.data_containers.ComponentProperties`
    - grouped component properites
      :py:class:`tespy.tools.data_containers.GroupedComponentProperties`
    - fluid composition
      :py:class:`tespy.tools.data_containers.FluidComposition`
    - fluid properties
      :py:class:`tespy.tools.data_containers.FluidProperties`

    Grouped component properties are used, if more than one component property
    has to be specified in order to apply one equation, e.g. pressure drop in
    pipes by specified length, diameter and roughness. If you specify all three
    of these properties, the DataContainer for the group will be created
    automatically!

    For the full list of available parameters for each data container, see its
    documentation.

    Example
    -------
    The examples below show the different (sub-)classes of DataContainers
    available.

    >>> from tespy.tools.data_containers import (
    ... ComponentCharacteristics, ComponentCharacteristicMaps,
    ... ComponentProperties, FluidComposition, GroupedComponentProperties,
    ... FluidProperties, DataContainerSimple)
    >>> from tespy.components import Pipe
    >>> type(ComponentCharacteristicMaps(is_set=True))
    <class 'tespy.tools.data_containers.ComponentCharacteristicMaps'>
    >>> type(ComponentCharacteristics(is_set=True, param='m'))
    <class 'tespy.tools.data_containers.ComponentCharacteristics'>
    >>> type(ComponentProperties(val=100, is_set=True, is_var=True,
    ...      max_val=1000, min_val=1))
    <class 'tespy.tools.data_containers.ComponentProperties'>
    >>> pi = Pipe('testpipe', L=100, D=0.5, ks=5e-5)
    >>> type(GroupedComponentProperties(is_set=True,
    ...      elements=[pi.L, pi.D, pi.ks], method='default'))
    <class 'tespy.tools.data_containers.GroupedComponentProperties'>
    >>> type(FluidComposition(
    ... val={'CO2': 0.1, 'H2O': 0.11, 'N2': 0.75, 'O2': 0.03},
    ... val_set={'CO2': False, 'H2O': False, 'N2': False, 'O2': True},
    ... balance=False))
    <class 'tespy.tools.data_containers.FluidComposition'>
    >>> type(FluidProperties(val=5, val_SI=500000, val_set=True, unit='bar',
    ...      ref=None, ref_set=False))
    <class 'tespy.tools.data_containers.FluidProperties'>
    >>> type(DataContainerSimple(val=5, is_set=False))
    <class 'tespy.tools.data_containers.DataContainerSimple'>
    """

    def __init__(self, **kwargs):

        var = self.attr()

        # default values
        for key in var.keys():
            self.__dict__.update({key: var[key]})

        self.set_attr(**kwargs)

    def set_attr(self, **kwargs):
        """
        Sets, resets or unsets attributes of a DataContainer type object.

        Parameters
        ----------
        **kwargs :
            See the class documentation of desired DataContainer for available
            keywords.
        """
        var = self.attr()
        # specify values
        for key in kwargs:
            if key in var:
                self.__dict__.update({key: kwargs[key]})

            else:
                msg = ('Data container of type ' + self.__class__.__name__ +
                       ' has no attribute ' + key + '.')
                logging.error(msg)
                raise KeyError(msg)

    def get_attr(self, key):
        """
        Get the value of a DataContainer's attribute.

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

    @staticmethod
    def attr():
        """
        Return the available attributes for a DataContainer type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {}


class ComponentCharacteristics(DataContainer):
    """
    Data container for component characteristics.

    Parameters
    ----------
    func : tespy.components.characteristics.characteristics
        Function to be applied for this characteristics, default: None.

    is_set : boolean
        Should this equation be applied?, default: is_set=False.

    param : str
        Which parameter should be applied as the x value?
        default: method='default'.
    """

    @staticmethod
    def attr():
        """
        Return the available attributes for a ComponentCharacteristics
        type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {
            'char_func': None, 'is_set': False, 'param': None,
            'func_params': {}, 'func': None, 'deriv': None, 'latex': None,
            'char_params': {'type': 'rel', 'inconn': 0, 'outconn': 0},
            'num_eq': 0}


class ComponentCharacteristicMaps(DataContainer):
    """
    Data container for characteristic maps.

    Parameters
    ----------
    func : tespy.components.characteristics.characteristics
        Function to be applied for this characteristic map, default: None.

    is_set : boolean
        Should this equation be applied?, default: is_set=False.

    param : str
        Which parameter should be applied as the x value?
        default: method='default'.
    """

    @staticmethod
    def attr():
        """
        Return the available attributes for a ComponentCharacteristicMaps type
        object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {
            'char_func': None, 'is_set': False, 'param': None, 'latex': None,
            'func_params': {}, 'func': None, 'deriv': None,
            'num_eq': 0}


class ComponentProperties(DataContainer):
    """
    Data container for component properties.

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
    """

    @staticmethod
    def attr():
        """
        Return the available attributes for a ComponentProperties type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {
            'val': 1, 'val_SI': 0, 'is_set': False, 'd': 1e-4,
            'min_val': -1e12, 'max_val': 1e12, 'is_var': False,
            'val_ref': 1, 'design': np.nan, 'is_result': False,
            'num_eq': 0, 'func_params': {}, 'func': None, 'deriv': None,
            'latex': None}


class FluidComposition(DataContainer):
    """
    Data container for fluid composition.

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

    @staticmethod
    def attr():
        """
        Return the available attributes for a FluidComposition type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {'val': {}, 'val0': {}, 'val_set': {},
                'design': collections.OrderedDict(), 'balance': False}


class GroupedComponentProperties(DataContainer):
    """
    Data container for grouped component parameters.

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

    @staticmethod
    def attr():
        """
        Return the available attributes for a GroupedComponentProperties type
        object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {'is_set': False, 'method': 'default', 'elements': [],
                'func': None, 'deriv': None, 'num_eq': 0, 'latex': None,
                'func_params': {}}


class GroupedComponentCharacteristics(DataContainer):
    """
    Data container for grouped component characteristics.

    Parameters
    ----------
    is_set : boolean
        Should the equation for this parameter group be applied?
        default: is_set=False.

    elements : list
        Which component properties are part of this component group?
        default elements=[].
    """

    @staticmethod
    def attr():
        """
        Return the available attributes for a GroupedComponentCharacteristics
        type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {'is_set': False, 'elements': [], 'func': None, 'deriv': None,
                'num_eq': 0, 'latex': None, 'func_params': {}}


class FluidProperties(DataContainer):
    """
    Data container for fluid properties.

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

    @staticmethod
    def attr():
        r"""
        Return the available attributes for a FluidProperties type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {'val': np.nan, 'val0': np.nan, 'val_SI': 0, 'val_set': False,
                'ref': None, 'ref_set': False, 'unit': None, 'design': np.nan}


class DataContainerSimple(DataContainer):
    """
    Simple data container without data type restrictions to val field.

    Parameters
    ----------
    val : no specific datatype
        Value for the property, no predefined datatype.

    is_set : boolean
        Has the value for this property been set? default: val_set=False.
    """

    @staticmethod
    def attr():
        r"""
        Return the available attributes for a DataContainerSimple type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {
            'val': np.nan, 'is_set': False,
            'func_params': {}, 'func': None, 'deriv': None, 'latex': None,
            'num_eq': 0}
