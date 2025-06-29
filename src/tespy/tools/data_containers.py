# -*- coding: utf-8

"""Module for data container classes.

The DataContainer class and its subclasses are used to store component or
connection properties.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/data_containers.py

SPDX-License-Identifier: MIT
"""
import numpy as np

from tespy.tools import logger


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
    ... FluidProperties, SimpleDataContainer)
    >>> from tespy.components import Pipe
    >>> type(ComponentCharacteristicMaps(is_set=True))
    <class 'tespy.tools.data_containers.ComponentCharacteristicMaps'>
    >>> type(ComponentCharacteristics(is_set=True, param='m'))
    <class 'tespy.tools.data_containers.ComponentCharacteristics'>
    >>> type(ComponentProperties(_val=100, is_set=True, _is_var=True,
    ...      max_val=1000, min_val=1))
    <class 'tespy.tools.data_containers.ComponentProperties'>
    >>> pi = Pipe('testpipe', L=100, D=0.5, ks=5e-5)
    >>> type(GroupedComponentProperties(
    ... is_set=True, elements=["L", "D", "ks"]
    ... ))
    <class 'tespy.tools.data_containers.GroupedComponentProperties'>
    >>> type(FluidComposition(
    ... _val={'CO2': 0.1, 'H2O': 0.11, 'N2': 0.75, 'O2': 0.03}, _is_set={'O2'}
    ... ))
    <class 'tespy.tools.data_containers.FluidComposition'>
    >>> type(FluidProperties(val=5, _val_SI=500000, is_set=True, unit='bar'))
    <class 'tespy.tools.data_containers.FluidProperties'>
    >>> type(SimpleDataContainer(_val=5, is_set=False))
    <class 'tespy.tools.data_containers.SimpleDataContainer'>
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
            elif f"_{key}" in var:
                self.__dict__.update({f"_{key}": kwargs[key]})
            else:
                msg = (
                    f"Datacontainer of type {self.__class__.__name__} has no "
                    f"attribute \"{key}\"."
                )
                logger.error(msg)
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
            msg = (
                f"Datacontainer of type {self.__class__.__name__} has no "
                f"attribute \"{key}\"."
            )
            logger.error(msg)
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

    def _serialize(self):
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
            'char_func': None,
            'is_set': False,
            'param': None,
            'func_params': {},
            'func': None,
            'deriv': None,
            'char_params': {'type': 'rel', 'inconn': 0, 'outconn': 0},
            'num_eq_sets': 0,
            '_num_eq': None,
            'structure_matrix': None,
            'dependents': None,
            'constant_deriv': False
        }

    def _serialize(self):
        export = {}
        if self.char_func is not None:
            export.update({"char_func": self.char_func._serialize()})

        for k in ["is_set", "param", "char_params"]:
            export.update({k: self.get_attr(k)})
        return export

    def get_num_eq(self):
        if self._num_eq is None:
            return self.num_eq_sets
        else:
            return self._num_eq

    def set_num_eq(self, value):
        self._num_eq = value

    num_eq = property(get_num_eq, set_num_eq)


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
            'char_func': None,
            'is_set': False,
            'param': None,
            'func_params': {},
            'func': None,
            'deriv': None,
            'num_eq_sets': 0,
            'structure_matrix': None,
            'constant_deriv': False,
            'dependents': None
        }

    def _serialize(self):
        export = {}
        if self.char_func is not None:
            export.update({"char_func": self.char_func._serialize()})

        for k in ["is_set", "param"]:
            export.update({k: self.get_attr(k)})
        return export

    def get_num_eq(self):
        if self._num_eq is None:
            return self.num_eq_sets
        else:
            return self._num_eq

    def set_num_eq(self, value):
        self._num_eq = value

    num_eq = property(get_num_eq, set_num_eq)


class ComponentMandatoryConstraints(DataContainer):
    """
    Data container for component mandatory constraints.
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
            'num_eq_sets': 0,
            '_num_eq': None,
            'func_params': {},
            'func': None,
            'deriv': None,
            'constant_deriv': False,
            'structure_matrix': None,
            'dependents': None
        }

    def _serialize(self):
        keys = self._serializable_keys()
        return {k: self.get_attr(k) for k in keys}

    @staticmethod
    def _serializable_keys():
        return [
            "val", "val_SI", "is_set", "d", "min_val", "max_val", "is_var",
        ]

    def get_num_eq(self):
        if self._num_eq is None:
            return self.num_eq_sets
        else:
            return self._num_eq

    def set_num_eq(self, value):
        self._num_eq = value

    num_eq = property(get_num_eq, set_num_eq)


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
            '_val': 1,
            'val_SI': 0,
            'is_set': False,
            'd': 1e-4,
            'min_val': -1e12,
            'max_val': 1e12,
            '_is_var': False,
            'design': np.nan,
            'is_result': False,
            'num_eq_sets': 0,
            '_num_eq': None,
            'func_params': {},
            'func': None,
            'deriv': None,
            'structure_matrix': None,
            'constant_deriv': False,
            '_reference_container': None,
            '_factor': None,
            '_offset': None,
            'dependents': None
        }

    def _serialize(self):
        keys = self._serializable_keys()
        return {k: getattr(self, k) for k in keys}

    @staticmethod
    def _serializable_keys():
        return [
            "val", "val_SI", "is_set", "d", "min_val", "max_val", "is_var",
        ]

    def get_num_eq(self):
        if self._num_eq is None:
            return self.num_eq_sets
        else:
            return self._num_eq

    def set_num_eq(self, value):
        self._num_eq = value

    def get_J_col(self):
        reference = self._reference_container
        if reference:
            return reference.J_col
        else:
            raise ValueError("")

    def get_reference_val_SI(self):
        """Get value of the reference corresponding to own value

        Returns
        -------
        float
            Value of reference container corresponding to this data container's
            value.
        """
        return (self._val - self._offset) / self._factor

    def set_reference_val_SI(self, value):
        if self._reference_container is not None:
            self._reference_container.val_SI = (value - self._offset) / self._factor
        else:
            raise ValueError()

    def get_J_col(self):
        if self._reference_container is not None:
            return self._reference_container.J_col
        else:
            raise ValueError()

    def get_val_SI(self):
        if self._reference_container is not None:
            return self._reference_container.val_SI * self._factor + self._offset
        else:
            return float(self._val)

    def set_val_SI(self, value):
        self._val = value

    def get_is_var(self):
        if self._reference_container is not None:
            return self._reference_container.is_var
        else:
            return self._is_var

    def set_is_var(self, value):
        if self._reference_container is not None:
            self._reference_container.is_var = value
        else:
            raise ValueError()

    num_eq = property(get_num_eq, set_num_eq)
    val = property(get_val_SI, set_val_SI)
    J_col = property(get_J_col)
    is_var = property(get_is_var, set_is_var)

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
        return {
            'is_set': False,
            'elements': [],
            '_num_eq': None,
            'func': None,
            'deriv': None,
            'num_eq_sets': 0,
            'func_params': {},
            'structure_matrix': None,
            'constant_deriv': False,
            'dependents': None
        }

    def get_num_eq(self):
        if self._num_eq is None:
            return self.num_eq_sets
        else:
            return self._num_eq

    def set_num_eq(self, value):
        self._num_eq = value

    num_eq = property(get_num_eq, set_num_eq)


class GroupedComponentCharacteristics(GroupedComponentProperties):
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
    pass

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

    is_set : boolean
        Has the value for this property been set? default: is_set=False.

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
        return {
            "design": np.nan,
            "val": np.nan,
            "val0": np.nan,
            "_val_SI": 0,
            "d": 1e-1,
            "unit": None,
            "is_set": False,
            "_potential_var": False,
            "func": None,
            "deriv": None,
            "structure_matrix": None,
            "constant_deriv": False,
            "num_eq": 0,
            "func_params": {},
            "_reference_container": None,
            "_offset": None,
            "_factor": None,
            'dependents': None
        }

    def _serialize(self):
        keys = ["val", "val0", "val_SI", "is_set", "unit"]
        return {k: getattr(self, k) for k in keys}

    def get_reference_val_SI(self):
        """Get value of the reference corresponding to own value

        Returns
        -------
        float
            Value of reference container corresponding to this data container's
            value.
        """
        return (self._val_SI - self._offset) / self._factor

    def set_reference_val_SI(self, value):
        if self._reference_container is not None:
            self._reference_container.val_SI = (value - self._offset) / self._factor
        else:
            raise ValueError()

    def get_J_col(self):
        if self._reference_container is not None:
            return self._reference_container.J_col
        else:
            raise ValueError()

    def get_val_SI(self):
        if self._reference_container is not None:
            return self._reference_container.val_SI * self._factor + self._offset
        else:
            return float(self._val_SI)

    def set_val_SI(self, value):
        self._val_SI = value

    def get_is_var(self):
        if self._reference_container is not None:
            return self._reference_container.is_var
        else:
            raise ValueError()

    def set_is_var(self, value):
        if self._reference_container is not None:
            self._reference_container.is_var = value
        else:
            raise ValueError()

    val_SI = property(get_val_SI, set_val_SI)
    J_col = property(get_J_col)
    is_var = property(get_is_var, set_is_var)


class ScalarVariable(DataContainer):

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
        return {
            "_val_SI": 0,
            "_is_var": True,
            "_J_col": None,
            "_d": 1e-4,
            "min_val": None,
            "max_val": None
        }

    def get_val_SI(self):
        return float(self._val_SI)

    def set_val_SI(self, value):
        self._val_SI = value

    def get_is_var(self):
        return self._is_var

    def set_is_var(self, value):
        if type(value) != bool:
            raise TypeError()

        self._is_var = value

    def get_J_col(self):
        if self.is_var:
            return self._J_col
        else:
            raise ValueError()

    def set_J_col(self, value):
        if self.is_var:
            self._J_col = value
        else:
            raise ValueError()

    def get_d(self):
        return self._d

    J_col = property(get_J_col, set_J_col)
    is_var = property(get_is_var, set_is_var)
    val_SI = property(get_val_SI, set_val_SI)
    d = property(get_d)


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

    is_set : dict
        Which fluid mass fractions have been set, default is_set={}.
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
        return {
            '_val': dict(),
            'val0': dict(),
            'd': 1e-5,
            '_is_set': set(),
            'design': dict(),
            'wrapper': dict(),
            'back_end': dict(),
            'engine': dict(),
            '_is_var': set(),
            '_J_col': dict(),
            '_reference_container': None,
            '_offset': None,
            '_factor': None
        }

    def _serialize(self):
        export = {"val": self.val}
        export["is_set"] = list(self.is_set)
        export["engine"] = {k: e.__name__ for k, e in self.engine.items()}
        export["back_end"] = {k: b for k, b in self.back_end.items()}
        return export

    def get_is_var(self):
        reference = self._reference_container
        if reference:
            return reference.is_var
        else:
            return self._is_var

    def get_J_col(self):
        reference = self._reference_container
        if reference:
            return reference.J_col
        else:
            raise ValueError("")

    def get_is_set(self):
        return self._is_set

    def get_val(self):
        reference = self._reference_container
        if reference:
            return {
                f: val * self._factor + self._offset
                for f, val in reference.val.items()
            }
        else:
            return self._val

    def set_val(self, value):
        self._val = value

    def get_reference_val(self):
        reference = self._reference_container
        if reference:
            return {
                f: val * self._factor + self._offset
                for f, val in reference.val.items()
            }
        else:
            return self._val

    val = property(get_val, set_val)
    is_set = property(get_is_set)
    is_var = property(get_is_var)
    J_col = property(get_J_col)

class VectorVariable(DataContainer):
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

    is_set : dict
        Which fluid mass fractions have been set, default is_set={}.
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
        return {
            "_val": dict(),
            "_is_var": set(),
            "_J_col": dict(),
            "_d": 1e-4
        }

    def get_val_SI(self):
        return self._val

    def set_val_SI(self, value):
        self._val = value

    def get_is_var(self):
        return self._is_var

    def set_is_var(self, value):
        if type(value) != set:
            raise TypeError()

        self._is_var = value

    def get_J_col(self):
        if self.is_var:
            return self._J_col
        else:
            raise ValueError()

    def set_J_col(self, value):
        if self.is_var:
            self._J_col = value
        else:
            raise ValueError()

    def get_d(self):
        return self._d

    J_col = property(get_J_col, set_J_col)
    is_var = property(get_is_var, set_is_var)
    val = property(get_val_SI, set_val_SI)
    d = property(get_d)


class ReferencedFluidProperties(DataContainer):

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
        return {
            "ref": None,
            "is_set": False,
            "unit": None,
            "func": None,
            "deriv": None,
            "structure_matrix": None,
            'constant_deriv': False,
            "num_eq": 0,
            "func_params": {},
            "_solved": False,
            "dependents": None,
        }

    def _serialize(self):
        if self.ref is not None:
            keys = ["is_set", "unit"]
            export = {k: self.get_attr(k) for k in keys}
            export["conn"] = self.ref.obj.label
            export["factor"] = self.ref.factor
            export["delta"] = self.ref.delta
            return export
        else:
            return {}


class SimpleDataContainer(DataContainer):
    """
    Simple data container without data type restrictions to val field.

    Parameters
    ----------
    val : no specific datatype
        Value for the property, no predefined datatype.

    is_set : boolean
        Has the value for this property been set? default: is_set=False.
    """

    @staticmethod
    def attr():
        r"""
        Return the available attributes for a SimpleDataContainer type object.

        Returns
        -------
        out : dict
            Dictionary of available attributes (dictionary keys) with default
            values.
        """
        return {
            "_val": np.nan,
            "is_set": False,
            "func_params": {},
            "func": None,
            "deriv": None,
            "num_eq_sets": 0,
            'constant_deriv': False,
            "_num_eq": None,
            "structure_matrix": None,
            "_solved": False,
            'dependents': None
        }

    def _serialize(self):
        return {"val": self.val, "is_set": self.is_set}

    def get_num_eq(self):
        if self._num_eq is None:
            return self.num_eq_sets
        else:
            return self._num_eq

    def set_num_eq(self, value):
        self._num_eq = value

    def get_val(self):
        return self._val

    def set_val(self, value):
        self._val = value

    num_eq = property(get_num_eq, set_num_eq)
    val = property(get_val, set_val)
