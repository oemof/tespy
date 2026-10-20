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
import pint
from tabulate import tabulate

from tespy.tools import logger
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap
from tespy.tools.global_vars import get_display_mode
from tespy.tools.units import SI_UNITS
from tespy.tools.units import is_compatible


def _is_numeric(potentially_a_number):
    """Checks if the value provided is a number by trying to convert it to
    float

    Parameters
    ----------
    potentially_a_number : any
        Value to check

    Returns
    -------
    bool
        True if the value is a number

    Example
    -------
    >>> from tespy.tools.data_containers import _is_numeric
    >>> _is_numeric(5)
    True
    >>> _is_numeric("var")
    False
    """
    try:
        float(potentially_a_number)
        return True
    except (TypeError, ValueError):
        return False


def _display_repr(obj):
    mode = get_display_mode()
    if mode != "none":
        if mode == "extensive" and hasattr(obj, "_repr_extensive"):
            return obj._repr_extensive()
        if hasattr(obj, "_repr_compact"):
            return obj._repr_compact()
    return object.__repr__(obj)


def _format_value(value):
    if isinstance(value, str):
        return value
    if np.isnan(value):
        return "nan"
    if value == 0 or 1e-4 <= abs(value) < 1e5:
        return f"{value:.4f}"
    return f"{value:.4e}"



class _NumEqMixin:
    """Mixin that provides the num_eq property for data containers."""

    def get_num_eq(self):
        if self._num_eq is None:
            return self.num_eq_sets
        else:
            return self._num_eq

    def set_num_eq(self, value):
        self._num_eq = value

    num_eq = property(get_num_eq, set_num_eq)


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
    - grouped component properties
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
    >>> type(ComponentProperties(
    ...     val=100, is_set=True, is_var=True, max_val=1000, min_val=1
    ... ))
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
    >>> type(FluidProperties(
    ...     val=5, val_SI=500000, is_set=True, quantity="pressure"
    ... ))
    <class 'tespy.tools.data_containers.FluidProperties'>
    >>> type(SimpleDataContainer(val=5, is_set=False))
    <class 'tespy.tools.data_containers.SimpleDataContainer'>
    """

    def __init__(self, **kwargs):

        var = self.attr()

        # default values
        for key in var.keys():
            self.__dict__.update({key: var[key]})

        self.set_attr(**kwargs)

    def __repr__(self):
        return _display_repr(self)

    __str__ = __repr__

    def set_attr(self, **kwargs):
        """
        Sets, resets or unsets attributes of a DataContainer type object.

        Parameters
        ----------
        **kwargs :
            See the class documentation of desired DataContainer for available
            keywords.
        """
        for key, value in kwargs.items():
            if key not in self.__dict__ and not hasattr(type(self), key):
                raise AttributeError(f"Unknown parameter '{key}'.")
            setattr(self, key, value)

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
        return {"result_only": False}

    def _serialize(self, units=None):
        return {}

    def accept(self, value):
        """
        Interpret a user-supplied value and update this container's state.

        Each subclass implements this to encapsulate the validation and
        assignment logic that formerly lived in Component.set_attr /
        Connection.set_attr.

        Parameters
        ----------
        value :
            The value supplied by the user (numeric, string, dict, None, ...).
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not implement accept()."
        )


class ComponentCharacteristics(_NumEqMixin, DataContainer):
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
            "char_func": None,
            "is_set": False,
            "param": None,
            "func_params": {},
            "func": None,
            "deriv": None,
            "char_params": {"type": "rel", "inconn": 0, "outconn": 0},
            "num_eq_sets": 0,
            "_num_eq": None,
            "structure_matrix": None,
            "dependents": None,
            "constant_deriv": False,
            "description": None,
            "result_only": False
        }

    def _serialize(self, units=None):
        export = {}
        if self.char_func is not None:
            export.update({"char_func": self.char_func._serialize()})

        for k in ["is_set", "param", "char_params"]:
            export.update({k: self.get_attr(k)})
        return export

    def _repr_compact(self):
        char = None
        if self.char_func is not None:
            char = type(self.char_func).__name__
        return (
            f"{type(self).__name__}(param={self.param!r}, "
            f"char_func={char}, is_set={self.is_set})"
        )

    def accept(self, value):
        if value is None:
            self.is_set = False
        elif isinstance(value, dict):
            self.set_attr(**value)
        elif isinstance(value, CharLine):
            self.char_func = value
        else:
            raise TypeError(
                f"Expected a CharLine, dict, or None for a "
                f"ComponentCharacteristics parameter, got {type(value).__name__}."
            )


class ComponentCharacteristicMaps(_NumEqMixin, DataContainer):
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
            "char_func": None,
            "is_set": False,
            "param": None,
            "func_params": {},
            "func": None,
            "deriv": None,
            "num_eq_sets": 0,
            "_num_eq": None,
            "structure_matrix": None,
            "constant_deriv": False,
            "dependents": None,
            "description": None,
            "result_only": False
        }

    def _serialize(self, units=None):
        export = {}
        if self.char_func is not None:
            export.update({"char_func": self.char_func._serialize()})

        for k in ["is_set", "param"]:
            export.update({k: self.get_attr(k)})
        return export

    def _repr_compact(self):
        char = None
        if self.char_func is not None:
            char = type(self.char_func).__name__
        return (
            f"{type(self).__name__}(param={self.param!r}, "
            f"char_func={char}, is_set={self.is_set})"
        )

    def accept(self, value):
        if value is None:
            self.is_set = False
        elif isinstance(value, dict):
            self.set_attr(**value)
        elif isinstance(value, CharMap):
            self.char_func = value
        else:
            raise TypeError(
                f"Expected a CharMap, dict, or None for a "
                f"ComponentCharacteristicMaps parameter, got {type(value).__name__}."
            )


class ComponentMandatoryConstraints(_NumEqMixin, DataContainer):
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
            "num_eq_sets": 0,
            "_num_eq": None,
            "func_params": {},
            "func": None,
            "deriv": None,
            "constant_deriv": False,
            "structure_matrix": None,
            "dependents": None,
            "description": None,
            "result_only": False
        }

    def _repr_compact(self):
        return (
            f"{type(self).__name__}(description={self.description!r}, "
            f"num_eq_sets={self.num_eq_sets})"
        )


class GroupedComponentProperties(_NumEqMixin, DataContainer):
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
            "is_set": False,
            "elements": [],
            "_num_eq": None,
            "func": None,
            "deriv": None,
            "num_eq_sets": 0,
            "func_params": {},
            "structure_matrix": None,
            "constant_deriv": False,
            "dependents": None,
            "description": None,
            "result_only": False
        }

    def accept(self, value):
        if value is None:
            self.is_set = False
        elif isinstance(value, bool):
            self.is_set = value
        elif isinstance(value, dict):
            self.set_attr(**value)
        else:
            raise TypeError(
                f"Expected a bool, dict, or None for a grouped parameter, "
                f"got {type(value).__name__}."
            )

    def _repr_compact(self):
        return (
            f"{type(self).__name__}(elements={self.elements}, "
            f"is_set={self.is_set})"
        )


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

class FluidProperties(_NumEqMixin, DataContainer):
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
            "_val": np.nan,
            "_val0": np.nan,
            "_val_SI": np.nan,
            "_unit_pinned": False,
            "_val0_unit": None,
            "_units": None,
            "_user_quantity": None,
            "_is_var": False,
            "is_result": False,
            "min_val": -1e12,
            "max_val": 1e12,
            "limit_scale": None,
            "d": 1e-4,
            "_unit": None,
            "is_set": False,
            "_potential_var": False,
            "func": None,
            "deriv": None,
            "structure_matrix": None,
            "constant_deriv": False,
            "num_eq_sets": 0,
            "_num_eq": None,
            "func_params": {},
            "calc": None,
            "calc_params": {},
            "calc_deps": [],
            "_reference_container": None,
            "_offset": None,
            "_factor": None,
            'dependents': None,
            "quantity": None,
            "description": None,
            "result_only": False
        }

    def _serialize(self, units=None):
        return {
            "val": self.val, "val_SI": self.val_SI, "is_set": self.is_set,
            "unit": self._unit_for_export(units)
        }

    def _unit_for_export(self, units):
        """Unit of the exported value: the unit the container is bound to or
        a specification is pinned to. A specification not yet bound to a
        network takes the default unit of the given unit system. Values that
        are neither specified nor calculated have no unit."""
        if self._unit is not None or units is None or not self.is_set:
            return self._unit
        return units.canonical(units.default[self.quantity])

    def _display_unit(self):
        if self._unit is None:
            return ""
        if self._units is not None:
            return f"{self._units.parsed_unit(self._unit):~}"
        return self._unit

    def _repr_compact(self):
        return (
            f"{type(self).__name__}(val={self.val:.6g}, "
            f"unit={self._display_unit()!r}, "
            f"val_SI={self.val_SI:.6g}, is_set={self.is_set})"
        )

    def _repr_extensive(self):
        spec = "var" if self.is_var else ("set" if self.is_set else "")
        rows = [[
            _format_value(self.val),
            self._display_unit(),
            f"{self.val_SI:.4e}",
            spec
        ]]
        return tabulate(
            rows, headers=["value", "unit", "SI value", "spec"],
            tablefmt="simple", disable_numparse=True,
            colalign=("right", "left", "right", "left")
        )

    def accept(self, value):
        if value is None:
            self.is_set = False
            self._unit = None
            self._unit_pinned = False
            self._user_quantity = None
        elif isinstance(value, dict):
            self.set_attr(**value)
        elif isinstance(value, pint.Quantity) or _is_numeric(value):
            self.val = value
            self.is_set = True
        else:
            raise TypeError(
                f"Expected a numeric value, pint.Quantity, dict, or None, "
                f"got {type(value).__name__}."
            )

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
            return self._is_var

    def set_is_var(self, value):
        if self._reference_container is not None:
            self._reference_container.is_var = value
        else:
            self._is_var = value

    def get_val_with_unit(self):
        if self._units is not None and self._unit is not None:
            return self._units.ureg.Quantity(self._val, self._unit)
        if self._user_quantity is not None:
            return self._user_quantity
        return self._val

    def get_val0_with_unit(self):
        if self._units is not None and self._val0_unit is not None:
            return self._units.ureg.Quantity(self._val0, self._val0_unit)
        return self._val0

    def _bind_unit(self, units):
        """Attach the network's unit system: a plain number takes the default
        unit, a number pinned to the earlier default is converted to the
        current one, explicitly pinned units are kept."""
        self._units = units
        default = units.canonical(units.default[self.quantity])
        if self._unit is None:
            self._unit = default
        elif not self._unit_pinned and self._unit != default:
            self._val = units.convert(self._val, self._unit, default)
            self._unit = default

    def _bind_val0_unit(self, units):
        self._units = units
        if self._val0_unit is None:
            self._val0_unit = units.canonical(units.default[self.quantity])

    def set_SI_from_val(self, units):
        self._bind_unit(units)
        self.val_SI = units.value_to_SI(self._val, self._unit, self.quantity)

    def set_SI_from_val0(self, units):
        self._bind_val0_unit(units)
        self.val_SI = units.value_to_SI(
            self._val0, self._val0_unit, self.quantity
        )

    def value_from_SI(self, units):
        """SI value expressed in the unit of this container without
        modifying the container."""
        if self._unit is None:
            unit = units.canonical(units.default[self.quantity])
        else:
            unit = self._unit
        return units.value_from_SI(self.val_SI, SI_UNITS[self.quantity], unit)

    def set_val_from_SI(self, units):
        self._bind_unit(units)
        self._val = units.value_from_SI(
            self.val_SI, SI_UNITS[self.quantity], self._unit
        )

    def set_val0_from_SI(self, units):
        self._bind_val0_unit(units)
        self._val0 = units.value_from_SI(
            self.val_SI, SI_UNITS[self.quantity], self._val0_unit
        )

    def detach(self):
        if self._reference_container is not None:
            self._val_SI = self.val_SI
        self._reference_container = None

    def get_val(self):
        return float(self._val)

    def set_val(self, value):
        if isinstance(value, pint.Quantity):
            self._check_compatible(str(value.units))
            self._val = value.magnitude
            self._unit = str(value.units)
            self._unit_pinned = True
            self._user_quantity = value
        else:
            self._val = value
            self._unit = None
            self._unit_pinned = False
            self._user_quantity = None

    def get_val0(self):
        return self._val0

    def set_val0(self, value):
        if isinstance(value, pint.Quantity):
            self._check_compatible(str(value.units))
            self._val0 = value.magnitude
            self._val0_unit = str(value.units)
        else:
            self._val0 = value
            self._val0_unit = None

    def _check_compatible(self, unit):
        if not is_compatible(self.quantity, unit):
            msg = f"Unit '{unit}' is not compatible with {self.quantity}."
            raise ValueError(msg)

    def get_unit(self):
        return self._unit

    def _replace_unit_for_compatibility(self, value):
        if value == "C":
            if self.quantity == "temperature":
                value = "degC"
            elif self.quantity == "temperature_difference":
                value = "delta_degC"
        elif "kgK" in value:
            value = value.replace("kgK", "kg/K")
        elif value == "-":
            value = "1"
        return value

    def set_unit(self, value):
        # a unit string describes the stored value, the container still
        # follows the default units of the network. Only quantities pin
        # the unit of a container
        self._unit_pinned = False
        if value is None:
            self._unit = None
        else:
            value = self._replace_unit_for_compatibility(value)
            self._check_compatible(value)
            self._unit = value

    def get_val0_unit(self):
        return self._val0_unit

    def set_val0_unit(self, value):
        if value is None:
            self._val0_unit = None
        else:
            value = self._replace_unit_for_compatibility(value)
            self._check_compatible(value)
            self._val0_unit = value

    unit = property(get_unit, set_unit)
    val0_unit = property(get_val0_unit, set_val0_unit)
    val = property(get_val, set_val)
    val0 = property(get_val0, set_val0)
    val_SI = property(get_val_SI, set_val_SI)
    val_with_unit = property(get_val_with_unit)
    val0_with_unit = property(get_val0_with_unit)
    J_col = property(get_J_col)
    is_var = property(get_is_var, set_is_var)


class ComponentArrayProperties(DataContainer):
    """Data container for array-valued component result properties.

    .. note::

        This class is currently intended for result-only use. Instances are
        populated automatically by the network after each solve and cannot
        meaningfully be set by the user. The API of this class may change in
        future versions.

    Attributes
    ----------
    val_SI : numpy.ndarray
        Values in SI units.

    val : numpy.ndarray
        Values in the network's user-specified unit. Populated by the network
        after each solve via :py:meth:`set_val_from_SI`.

    quantity : str
        Physical quantity type (e.g. :code:`"heat"`, :code:`"temperature"`).
    """

    @staticmethod
    def attr():
        return {
            "val": None,
            "val_SI": None,
            "unit": None,
            "_units": None,
            "quantity": None,
            "structure_matrix": None,
            "func": None,
            "is_set": False,
            "num_eq_sets": 0,
            "description": None,
            "result_only": True,
        }

    def accept(self, value):
        if self.result_only:
            raise TypeError("This parameter is result-only and cannot be set by the user.")

    def set_val_from_SI(self, units):
        if self.val_SI is None:
            return
        self._units = units
        self.unit = units.canonical(units.default[self.quantity])
        self.val = units.value_from_SI(
            self.val_SI, SI_UNITS[self.quantity], self.unit
        )

    def _display_unit(self):
        if self.unit is None:
            return ""
        if self._units is not None:
            return f"{self._units.parsed_unit(self.unit):~}"
        return self.unit

    def _si_display_unit(self):
        si_unit = SI_UNITS[self.quantity]
        if self._units is not None:
            return f"{self._units.parsed_unit(si_unit):~}"
        return si_unit

    def _repr_compact(self):
        val = self.val if self.val is not None else self.val_SI
        if val is not None:
            val = np.array2string(np.asarray(val), precision=4, threshold=6)
        return (
            f"{type(self).__name__}(val={val}, "
            f"unit={self._display_unit()!r}, quantity={self.quantity!r})"
        )

    def _repr_extensive(self):
        val = self.val if self.val is not None else self.val_SI
        if val is None:
            return f"{type(self).__name__}: no values"
        val = np.atleast_1d(val)
        val_SI = None
        if self.val_SI is not None:
            val_SI = np.atleast_1d(self.val_SI)
        if val.ndim > 1:
            return self._repr_compact()
        opts = np.get_printoptions()
        summarize = len(val) > opts["threshold"]
        if summarize:
            edge = opts["edgeitems"]
            indices = [*range(edge), None, *range(len(val) - edge, len(val))]
        else:
            indices = range(len(val))
        rows = []
        for i in indices:
            if i is None:
                rows.append(["...", "...", "..."])
            else:
                rows.append([
                    str(i),
                    _format_value(val[i]),
                    _format_value(val_SI[i]) if val_SI is not None else ""
                ])
        lines = []
        if self.quantity is not None:
            title = self.quantity
            unit = self._display_unit()
            if unit:
                title += f" in {unit}"
            si_unit = self._si_display_unit()
            if si_unit:
                title += f" (SI value in {si_unit})"
            lines += [title, ""]
        lines.append(tabulate(
            rows, headers=["index", "value", "SI value"], tablefmt="simple",
            disable_numparse=True, colalign=("left", "right", "right")
        ))
        if summarize:
            lines.append(
                f"({len(val)} values, adjust numpy.set_printoptions or "
                "access val/val_SI directly to see all)"
            )
        return "\n".join(lines)


class ComponentProperties(FluidProperties):

    @staticmethod
    def attr():
        attrs = FluidProperties.attr()
        attrs["_allows_var"] = False
        return attrs

    def _serialize(self, units=None):
        export = super()._serialize(units)
        export["is_var"] = self.is_var
        return export

    def _repr_compact(self):
        return (
            f"{type(self).__name__}(val={self.val:.6g}, "
            f"unit={self._display_unit()!r}, "
            f"val_SI={self.val_SI:.6g}, is_set={self.is_set}, "
            f"is_var={self.is_var})"
        )

    def accept(self, value):
        if value == "var":
            if not self._allows_var:
                raise ValueError(
                    "This parameter cannot be made a system variable."
                )
            self.is_set = True
            self.is_var = True
        else:
            super().accept(value)
            self.is_var = False


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
            "max_val": None,
            "result_only": False
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
            "_val": dict(),
            "val0": dict(),
            "d": 1e-5,
            "_is_set": set(),
            "design": dict(),
            "wrapper": dict(),
            "back_end": dict(),
            "engine": dict(),
            "wrapper_kwargs": dict(),
            "description": None,
            "quantity": None,
            "_is_var": set(),
            "_J_col": dict(),
            "_reference_container": None,
            "_offset": None,
            "_factor": None,
            "result_only": False
        }

    def _serialize(self, units=None):
        export = {"val": self.val}
        export["is_set"] = list(self.is_set)
        export["engine"] = {k: e.__name__ for k, e in self.engine.items()}
        export["back_end"] = {k: b for k, b in self.back_end.items()}
        return export

    def _repr_compact(self):
        return (
            f"{type(self).__name__}(val={self.val}, "
            f"is_set={sorted(self.is_set)})"
        )

    def _repr_extensive(self):
        if not self.val:
            return "no fluid composition"
        rows = [
            [f, _format_value(x), "set" if f in self.is_set else ""]
            for f, x in self.val.items()
        ]
        return tabulate(
            rows, headers=["fluid", "mass fraction", "spec"],
            tablefmt="simple", disable_numparse=True,
            colalign=("left", "right", "left")
        )

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

    def set_is_set(self, value):
        self._is_set = value

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

    def set_reference_val(self, key, value):
        if self._reference_container is not None:
            self._reference_container.val[key] = value
        else:
            raise ValueError()

    def detach(self):
        if self._reference_container is not None:
            self._val = self.val
        self._reference_container = None

    val = property(get_val, set_val)
    is_set = property(get_is_set, set_is_set)
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
            "_d": 1e-4,
            "result_only": False
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
            "constant_deriv": False,
            "num_eq": 0,
            "func_params": {},
            "_solved": False,
            "dependents": None,
            "quantity": None,
            "description": None,
            "result_only": False
        }

    def _serialize(self, units=None):
        if self.ref is not None:
            keys = ["is_set", "unit"]
            export = {k: self.get_attr(k) for k in keys}
            export["conn"] = self.ref.obj.label
            export["factor"] = self.ref.factor
            if isinstance(self.ref.delta, pint.Quantity):
                export["delta"] = self.ref.delta.magnitude
                export["delta_unit"] = str(self.ref.delta.units)
            else:
                export["delta"] = self.ref.delta
            return export
        else:
            return {}

    def set_SI_from_val(self, units):
        ref = self.ref
        si_unit = SI_UNITS[self.quantity]
        if isinstance(ref.delta, pint.Quantity):
            ref.delta_SI = units.convert(
                ref.delta.magnitude, str(ref.delta.units), si_unit
            )
            return
        default = units.canonical(units.default[self.quantity])
        if ref._delta_unit is None:
            ref._delta_unit = default
        elif ref.delta_SI is not None and ref._delta_unit != default:
            ref.delta = units.convert(ref.delta, ref._delta_unit, default)
            ref._delta_unit = default
        ref.delta_SI = units.convert(ref.delta, ref._delta_unit, si_unit)

    def _repr_compact(self):
        return f"{type(self).__name__}(ref={self.ref!r}, is_set={self.is_set})"


class SimpleDataContainer(_NumEqMixin, DataContainer):
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
            'dependents': None,
            "description": None,
            "dtype": None,
            "result_only": False,
        }

    def _serialize(self, units=None):
        return {"val": self.val, "is_set": self.is_set}

    def _repr_compact(self):
        return f"{type(self).__name__}(val={self.val!r}, is_set={self.is_set})"

    def accept(self, value):
        if value is None:
            self.is_set = False
        else:
            self.val = value
            self.is_set = True

    def get_val(self):
        return self._val

    def set_val(self, value):
        self._val = value

    val = property(get_val, set_val)
