# -*- coding: utf-8

"""Module for data container classes.

The data containers hold the parameters of components and connections. Every
container answers two independent questions: what kind of value it stores and
whether it can put an equation into the system of equations.

- :py:class:`ConnectionParameter` and :py:class:`ComponentParameter` hold a
  scalar with a physical quantity. They can be system variables and can carry
  an equation.
- :py:class:`ArrayParameter` holds an array with a physical quantity, as input
  or as result, and can carry an equation.
- :py:class:`ConfigParameter` holds a typed configuration value without unit
  and without equation.
- :py:class:`ParameterGroup` activates one equation, either when all of its
  element parameters are specified or when switched on directly.
- :py:class:`CharacteristicData` holds a characteristic function with a
  physical quantity per axis, :py:class:`CharacteristicEquation` additionally
  activates an equation.
- :py:class:`MandatoryConstraint` is an always active equation of a component.
- :py:class:`ReferenceEquation` ties a connection parameter to the parameter
  of another connection.
- :py:class:`FluidComposition` holds the mass fractions of a mixture.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/data_containers.py

SPDX-License-Identifier: MIT
"""
import warnings

import numpy as np
import pint
from tabulate import tabulate

from tespy.tools import logger
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap
from tespy.tools.characteristics import characteristic_from_dict
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


def _equation_attrs():
    """Attributes of a container that can place an equation in the system.

    :code:`func`, :code:`deriv`, :code:`dependents` and
    :code:`structure_matrix` are declared by the component or connection,
    :code:`num_eq_sets` is the number of equations per activation. The
    underscored attributes are written by the solver preparation.
    """
    return {
        "func": None,
        "deriv": None,
        "dependents": None,
        "structure_matrix": None,
        "constant_deriv": False,
        "func_params": {},
        "num_eq_sets": 0,
        "_num_eq": None,
        "_first_eq_index": None,
        "_scalar_dependents": None,
        "_vector_dependents": None,
    }


class _EquationMixin:
    """Provides the :code:`num_eq` property for equation carrying containers.

    :code:`num_eq` defaults to :code:`num_eq_sets` and can be overridden,
    e.g. by connections that count their equations individually.
    """

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

    - scalar connection parameters
      :py:class:`tespy.tools.data_containers.ConnectionParameter`
    - scalar component parameters
      :py:class:`tespy.tools.data_containers.ComponentParameter`
    - array valued component parameters
      :py:class:`tespy.tools.data_containers.ArrayParameter`
    - configuration values
      :py:class:`tespy.tools.data_containers.ConfigParameter`
    - parameter groups and equation switches
      :py:class:`tespy.tools.data_containers.ParameterGroup`
    - characteristic functions as data
      :py:class:`tespy.tools.data_containers.CharacteristicData`
    - characteristic functions with an equation
      :py:class:`tespy.tools.data_containers.CharacteristicEquation`
    - mandatory equations of components
      :py:class:`tespy.tools.data_containers.MandatoryConstraint`
    - references between connections
      :py:class:`tespy.tools.data_containers.ReferenceEquation`
    - fluid composition
      :py:class:`tespy.tools.data_containers.FluidComposition`

    A :py:class:`ParameterGroup` is used when more than one component
    parameter has to be specified in order to apply one equation, e.g.
    pressure drop in pipes by specified length, diameter and roughness. If you
    specify all three of these parameters, the group is activated
    automatically.

    For the full list of available parameters for each data container, see its
    documentation.

    Example
    -------
    The examples below show the different (sub-)classes of DataContainers
    available.

    >>> from tespy.tools.data_containers import (
    ... CharacteristicData, CharacteristicEquation, ComponentParameter,
    ... ConfigParameter, ConnectionParameter, FluidComposition,
    ... ParameterGroup)
    >>> type(CharacteristicData(input_quantities=[None, None]))
    <class 'tespy.tools.data_containers.CharacteristicData'>
    >>> type(CharacteristicEquation(is_set=True, param='m'))
    <class 'tespy.tools.data_containers.CharacteristicEquation'>
    >>> type(ComponentParameter(
    ...     val=100, is_set=True, is_var=True, max_val=1000, min_val=1
    ... ))
    <class 'tespy.tools.data_containers.ComponentParameter'>
    >>> type(ParameterGroup(is_set=True, elements=["L", "D", "ks"]))
    <class 'tespy.tools.data_containers.ParameterGroup'>
    >>> type(FluidComposition(
    ... _val={'CO2': 0.1, 'H2O': 0.11, 'N2': 0.75, 'O2': 0.03}, _is_set={'O2'}
    ... ))
    <class 'tespy.tools.data_containers.FluidComposition'>
    >>> type(ConnectionParameter(
    ...     val=5, val_SI=500000, is_set=True, quantity="pressure"
    ... ))
    <class 'tespy.tools.data_containers.ConnectionParameter'>
    >>> type(ConfigParameter(val=5, dtype="int", is_set=False))
    <class 'tespy.tools.data_containers.ConfigParameter'>
    """

    # consumers probe every parameter for an equation, containers without
    # equation answer with None and do not accept a value for these
    func = None
    structure_matrix = None

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
            is_property = isinstance(getattr(type(self), key, None), property)
            if key not in self.__dict__ and not is_property:
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
        return {}

    def _serialize(self, units=None):
        return {}

    @property
    def result_only(self):
        """Can the user not specify this parameter?"""
        return False

    def accept(self, value):
        """
        Interpret a user-supplied value and update this container's state.

        Each subclass implements this to encapsulate the validation and
        assignment logic of :code:`Component.set_attr` and
        :code:`Connection.set_attr`.

        Parameters
        ----------
        value :
            The value supplied by the user (numeric, string, dict, None, ...).
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not implement accept()."
        )


class ScalarParameter(_EquationMixin, DataContainer):
    """
    Base class for scalar parameters with a physical quantity.

    The value is stored in the unit of the parameter (:code:`val`) and in SI
    units (:code:`val_SI`). A parameter can be a system variable of the
    solver and can carry an equation. The solver groups variables with an
    affine relation and lets all of them read their value through a shared
    reference container.

    Parameters
    ----------
    val : float, pint.Quantity
        Value in the unit of the parameter. A plain number takes the default
        unit of the network for the quantity, a :code:`pint.Quantity` pins
        its own unit, default: :code:`np.nan`.

    val_SI : float
        Value in SI units, default: :code:`np.nan`.

    unit : str
        Unit of :code:`val`, default: :code:`None` (network default).

    is_set : boolean
        Has the value been specified? default: :code:`False`.

    quantity : str
        Physical quantity of the parameter, e.g. :code:`"pressure"`.

    min_val, max_val : float
        Bounds of the value when it is a system variable, default:
        :code:`-1e12` and :code:`1e12`.

    d : float
        Increment for numerical derivatives, default: :code:`1e-4`.

    val0 : float, pint.Quantity
        Starting value for the solver in the unit of the parameter, default:
        :code:`np.nan`.

    is_result : boolean
        Report the value in the results even if it is neither specified nor
        bound to an equation, default: :code:`False`.

    description : str
        Description of the parameter for documentation and schema.
    """

    @staticmethod
    def attr():
        return {
            "_val": np.nan,
            "_val0": np.nan,
            "_val0_unit": None,
            "_val_SI": np.nan,
            "_unit": None,
            "_unit_pinned": False,
            "_units": None,
            "_user_quantity": None,
            "quantity": None,
            "description": None,
            "is_set": False,
            "is_result": False,
            "design": np.nan,
            "_is_var": False,
            "_potential_var": False,
            "_reference_container": None,
            "_factor": None,
            "_offset": None,
            "sm_col": None,
            "d": 1e-4,
            "min_val": -1e12,
            "max_val": 1e12,
            **_equation_attrs(),
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

    def set_SI_from_val(self, units):
        self._bind_unit(units)
        self.val_SI = units.value_to_SI(self._val, self._unit, self.quantity)

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

    def get_val0_with_unit(self):
        if self._units is not None and self._val0_unit is not None:
            return self._units.ureg.Quantity(self._val0, self._val0_unit)
        return self._val0

    def _bind_val0_unit(self, units):
        self._units = units
        if self._val0_unit is None:
            self._val0_unit = units.canonical(units.default[self.quantity])

    def set_SI_from_val0(self, units):
        self._bind_val0_unit(units)
        self.val_SI = units.value_to_SI(
            self._val0, self._val0_unit, self.quantity
        )

    def set_val0_from_SI(self, units):
        self._bind_val0_unit(units)
        self._val0 = units.value_from_SI(
            self.val_SI, SI_UNITS[self.quantity], self._val0_unit
        )

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

    def get_val0_unit(self):
        return self._val0_unit

    def set_val0_unit(self, value):
        if value is None:
            self._val0_unit = None
        else:
            value = self._replace_unit_for_compatibility(value)
            self._check_compatible(value)
            self._val0_unit = value

    val0 = property(get_val0, set_val0)
    val0_unit = property(get_val0_unit, set_val0_unit)
    val0_with_unit = property(get_val0_with_unit)

    unit = property(get_unit, set_unit)
    val = property(get_val, set_val)
    val_SI = property(get_val_SI, set_val_SI)
    val_with_unit = property(get_val_with_unit)
    J_col = property(get_J_col)
    is_var = property(get_is_var, set_is_var)


class ConnectionParameter(ScalarParameter):
    """
    Data container for scalar connection parameters.

    Takes the parameters of :py:class:`ScalarParameter`.

    Parameters
    ----------
    reference : ReferenceEquation
        The equation applied when the parameter is specified as reference to
        the same parameter of another connection, default: :code:`None`,
        the parameter cannot be referenced.
    """

    @staticmethod
    def attr():
        return {
            **ScalarParameter.attr(),
            "reference": None,
        }


class ComponentParameter(ScalarParameter):
    """
    Data container for scalar component parameters.

    Extends :py:class:`ScalarParameter` with the specification as system
    variable through :code:`"var"`, with result post-processing and with
    bound scaling. A parameter is an input when it carries an equation, is an
    element of a parameter group or can be a system variable. A result
    parameter without any of these is result-only and refuses user values.

    Parameters
    ----------
    is_var : boolean
        Is the parameter a system variable? default: :code:`False`.

    calc : callable
        Method calculating the SI value after convergence.

    calc_params : dict
        Keyword arguments for :code:`calc`.

    calc_deps : list
        Names of parameters :code:`calc` depends on, evaluated first.

    limit_scale : str
        Name of a parameter whose magnitude scales the bound tolerance.
    """

    @staticmethod
    def attr():
        return {
            **ScalarParameter.attr(),
            "_allows_var": False,
            "_group_element": False,
            "limit_scale": None,
            "calc": None,
            "calc_params": {},
            "calc_deps": [],
        }

    def _is_input(self):
        return (
            not self.is_result
            or self.func is not None
            or self.structure_matrix is not None
            or self._allows_var
            or self._group_element
        )

    @property
    def result_only(self):
        return not self._is_input()

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
        if isinstance(value, str) and value == "var":
            if not self._allows_var:
                raise ValueError(
                    "This parameter cannot be made a system variable."
                )
            self.is_set = True
            self.is_var = True
        elif value is not None and self.result_only:
            raise TypeError(
                "This parameter is result-only and cannot be set by the user."
            )
        else:
            super().accept(value)
            self.is_var = False


class ArrayParameter(_EquationMixin, DataContainer):
    """
    Data container for array valued component parameters.

    All entries of the array share one physical quantity. The container is a
    result by default, an input array declares :code:`is_result=False` or
    carries an equation.

    Parameters
    ----------
    val : numpy.ndarray
        Values in the unit of the parameter.

    val_SI : numpy.ndarray
        Values in SI units.

    quantity : str
        Physical quantity of the entries, e.g. :code:`"heat"`.

    is_result : boolean
        Is the array a pure result? default: :code:`True`.

    calc : callable
        Method calculating the SI values after convergence.

    calc_params : dict
        Keyword arguments for :code:`calc`.

    calc_deps : list
        Names of parameters :code:`calc` depends on, evaluated first.

    description : str
        Description of the parameter for documentation and schema.
    """

    @staticmethod
    def attr():
        return {
            "val": None,
            "val_SI": None,
            "unit": None,
            "_units": None,
            "quantity": None,
            "description": None,
            "is_set": False,
            "is_result": True,
            "_group_element": False,
            "design": None,
            "calc": None,
            "calc_params": {},
            "calc_deps": [],
            **_equation_attrs(),
        }

    def accept(self, value):
        if value is None:
            self.is_set = False
        elif isinstance(value, dict):
            self.set_attr(**value)
        elif self._is_input():
            self.val = np.asarray(value, dtype=float)
            self.is_set = True
        else:
            raise TypeError(
                "This parameter is result-only and cannot be set by the user."
            )

    def _is_input(self):
        return (
            not self.is_result
            or self.func is not None
            or self.structure_matrix is not None
            or self._group_element
        )

    @property
    def result_only(self):
        return not self._is_input()

    def _default_unit(self, units):
        return units.canonical(units.default[self.quantity])

    def set_SI_from_val(self, units):
        if self.val is None:
            return
        self._units = units
        if self.unit is None:
            self.unit = self._default_unit(units)
        self.val_SI = units.value_to_SI(
            np.asarray(self.val, dtype=float), self.unit, self.quantity
        )

    def set_val_from_SI(self, units):
        if self.val_SI is None:
            return
        self._units = units
        self.unit = self._default_unit(units)
        self.val = units.value_from_SI(
            self.val_SI, SI_UNITS[self.quantity], self.unit
        )

    def _serialize(self, units=None):
        if not self.is_set:
            return {}
        unit = self.unit
        if unit is None and units is not None:
            unit = self._default_unit(units)
        return {
            "val": np.asarray(self.val).tolist(), "unit": unit,
            "is_set": self.is_set
        }

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


_CONFIG_DTYPES = {
    "int": (int, np.integer),
    "float": (int, float, np.number),
    "bool": (bool, np.bool_),
    "str": (str,),
    "dict": (dict,),
    "list": (list, tuple, np.ndarray),
}


class ConfigParameter(DataContainer):
    """
    Data container for configuration values.

    A configuration value has no physical unit and no equation, e.g. the
    number of outlets of a splitter or the flow arrangement of a heat
    exchanger.

    Parameters
    ----------
    val : int, float, bool, str, dict, list
        The value, default: :code:`np.nan`.

    dtype : str
        Expected type of the value: :code:`"int"`, :code:`"float"`,
        :code:`"bool"`, :code:`"str"`, :code:`"dict"` or :code:`"list"`.
        Default: :code:`None`, any type is accepted.

    is_set : boolean
        Has the value been set? default: :code:`False`.

    description : str
        Description of the parameter for documentation and schema.
    """

    @staticmethod
    def attr():
        return {
            "_val": np.nan,
            "is_set": False,
            "dtype": None,
            "description": None,
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
        empty = value is None or (isinstance(value, float) and np.isnan(value))
        if self.dtype is not None and not empty:
            if not isinstance(value, _CONFIG_DTYPES[self.dtype]):
                raise TypeError(
                    f"Expected a value of type {self.dtype}, got "
                    f"{type(value).__name__}."
                )
        self._val = value

    val = property(get_val, set_val)


class ParameterGroup(_EquationMixin, DataContainer):
    """
    Data container for a group of parameters that jointly activate one
    equation.

    With :code:`auto_activate=True` the group is active when all of its
    :code:`elements` are specified. With :code:`auto_activate=False` the
    group is switched on and off by the user, e.g. with
    :code:`component.set_attr(cone=True)`.

    Parameters
    ----------
    is_set : boolean
        Should the equation of this group be applied? default: :code:`False`.

    elements : list
        Names of the parameters that belong to this group, default:
        :code:`[]`.

    auto_activate : boolean
        Derive :code:`is_set` from the elements, default: :code:`True`.

    description : str
        Description of the equation for documentation and schema.
    """

    @staticmethod
    def attr():
        return {
            "is_set": False,
            "elements": [],
            "auto_activate": True,
            "description": None,
            **_equation_attrs(),
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
                f"Expected a bool, dict, or None for a parameter group, "
                f"got {type(value).__name__}."
            )

    def _serialize(self, units=None):
        return {"is_set": self.is_set}

    def _repr_compact(self):
        return (
            f"{type(self).__name__}(elements={self.elements}, "
            f"is_set={self.is_set})"
        )

    def get_val(self):
        return self.is_set

    def set_val(self, value):
        self.is_set = bool(value)

    val = property(get_val, set_val)


class MandatoryConstraint(_EquationMixin, DataContainer):
    """
    Data container for the mandatory equations of a component.

    Mandatory constraints are always active and never specified by the user.

    Parameters
    ----------
    description : str
        Description of the equation for documentation and schema.
    """

    @staticmethod
    def attr():
        return {
            "description": None,
            **_equation_attrs(),
        }

    def _repr_compact(self):
        return (
            f"{type(self).__name__}(description={self.description!r}, "
            f"num_eq_sets={self.num_eq_sets})"
        )


_PARAM_QUANTITIES = {
    "m": "mass_flow",
    "m_out": "mass_flow",
    "v": "volumetric_flow",
    "pr": None,
}


class CharacteristicData(DataContainer):
    """
    Data container for a characteristic function.

    The container holds a characteristic function object, e.g. a
    :py:class:`tespy.tools.characteristics.CharLine`,
    :py:class:`tespy.tools.characteristics.CharMap` or
    :py:class:`tespy.tools.characteristics.Polynomial`, and knows the
    physical quantity of every axis. Equations evaluate the characteristic
    through :py:meth:`evaluate` with SI inputs and receive SI outputs, the
    conversion to the units of the characteristic data is handled here. The
    unit of an axis is the one declared with :code:`input_units` or
    :code:`output_unit`, otherwise the one the characteristic function
    carries itself (e.g. the °C of an
    :py:class:`tespy.tools.characteristics.EN12900Polynomial`), otherwise
    the default unit of the network, like a plain number of a scalar
    parameter. The default characteristics shipped with tespy are SI data.

    A characteristic function matches the container when its number of inputs
    equals the number of input quantities declared, the representation is
    not restricted otherwise.

    Parameters
    ----------
    char_func : CharLine, CharMap, Polynomial
        The characteristic function, default: :code:`None`.

    input_quantities : list
        Physical quantity per input axis, :code:`None` for a dimensionless
        axis. Default: :code:`None`, which means one dimensionless input.

    output_quantity : str
        Physical quantity of the output, :code:`None` for dimensionless.

    input_units : list
        Unit per input axis the data is given in. Default: :code:`None`, the
        default unit of the network for the quantity.

    output_unit : str
        Unit the output data is given in. Default: :code:`None`, the default
        unit of the network for the quantity.

    param : str
        Connection parameter used as first input, e.g. :code:`"m"`,
        default: :code:`None`.

    char_params : dict
        Evaluation of the first input: :code:`type` is :code:`"rel"` for
        the ratio to the design point and :code:`"abs"` for the absolute
        value, :code:`inconn` and :code:`outconn` select the connections.

    is_set : boolean
        Has the characteristic been specified? default: :code:`False`.
        Providing a characteristic function sets it.

    description : str
        Description of the parameter for documentation and schema.
    """

    @staticmethod
    def attr():
        return {
            "char_func": None,
            "input_quantities": None,
            "output_quantity": None,
            "input_units": None,
            "output_unit": None,
            "param": None,
            "char_params": {"type": "rel", "inconn": 0, "outconn": 0},
            "is_set": False,
            "description": None,
            "_units": None,
            "_axes": None,
            "_label": None,
            "_is_default": False,
        }

    def set_attr(self, **kwargs):
        if isinstance(kwargs.get("char_func"), dict):
            kwargs["char_func"] = characteristic_from_dict(kwargs["char_func"])
        super().set_attr(**kwargs)
        self._axes = None

    def get_val(self):
        return self.char_func

    val = property(get_val)

    def get_num_inputs(self):
        if self.input_quantities is not None:
            return len(self.input_quantities)
        if self.char_func is not None:
            return self.char_func.num_inputs
        return 1

    num_inputs = property(get_num_inputs)

    def _resolved_input_quantities(self):
        """Quantities of the input axes: declared ones, otherwise
        dimensionless except for an absolute first input that follows the
        connection parameter it is taken from."""
        if self.input_quantities is not None:
            return list(self.input_quantities)
        quantities = [None] * self.num_inputs
        if self.char_params.get("type") == "abs" and self.param is not None:
            quantities[0] = _PARAM_QUANTITIES.get(self.param)
        return quantities

    def _check_signature(self, char_func):
        if char_func.num_inputs != self.num_inputs:
            raise ValueError(
                f"The characteristic {type(char_func).__name__} takes "
                f"{char_func.num_inputs} input(s), this parameter expects "
                f"{self.num_inputs}."
            )

    def accept(self, value):
        if value is None:
            self.is_set = False
        elif isinstance(value, dict):
            if "char_func" in value:
                self._check_signature(value["char_func"])
            self.set_attr(**value)
        elif hasattr(value, "evaluate") and hasattr(value, "num_inputs"):
            self._check_signature(value)
            self.set_attr(char_func=value, is_set=True, _is_default=False)
        else:
            raise TypeError(
                "Expected a characteristic function (CharLine, CharMap, "
                f"Polynomial), dict, or None, got {type(value).__name__}."
            )

    def set_default(self, char_func):
        """Assign a default characteristic shipped with tespy, its data are
        SI regardless of the network's default units."""
        self.set_attr(char_func=char_func, _is_default=True)

    def _bind_units(self, units, label=None):
        self._units = units
        self._label = label
        self._axes = None

    def _axis_units(self):
        """Resolve the unit of every axis once per binding: tuples of the
        SI unit and the data unit, :code:`None` for dimensionless axes."""
        if self._axes is None:
            quantities = self._resolved_input_quantities()
            input_units = self._declared_units()
            axes = []
            for quantity, unit in zip(quantities, input_units):
                axes.append(self._axis(quantity, unit))
            output_unit = self.output_unit
            if output_unit is None:
                output_unit = getattr(self.char_func, "output_unit", None)
            self._axes = (axes, self._axis(self.output_quantity, output_unit))
            self._warn_undeclared_units()
        return self._axes

    def _declared_units(self):
        """Units of the input axes declared on the container, otherwise
        carried by the characteristic function, :code:`None` per axis
        without either."""
        if self.input_units is not None:
            return list(self.input_units)
        carried = getattr(self.char_func, "input_units", None)
        if carried is not None:
            return list(carried)
        return [None] * self.num_inputs

    def _axis(self, quantity, unit):
        if quantity is None or self._units is None:
            return None
        if unit is None:
            if self._is_default:
                unit = SI_UNITS[quantity]
            else:
                unit = self._units.canonical(self._units.default[quantity])
        return (SI_UNITS[quantity], unit)

    def _warn_undeclared_units(self):
        """Characteristic data were SI before version 0.12 and follow the
        network's default units now, affected users must check their data."""
        if self._is_default or self.char_func is None:
            return
        axes, output = self._axis_units()
        quantities = self._resolved_input_quantities() + [self.output_quantity]
        declared = self._declared_units() + [
            self.output_unit
            if self.output_unit is not None
            else getattr(self.char_func, "output_unit", None)
        ]
        changed = [
            f"{quantity} in {axis[1]}"
            for quantity, unit, axis in zip(quantities, declared, axes + [output])
            if axis is not None and unit is None
            and self._units.canonical(axis[0]) != axis[1]
        ]
        if changed:
            label = f" of {self._label}" if self._label else ""
            warnings.warn(
                f"The characteristic data{label} are interpreted in the "
                f"default units of the network ({', '.join(changed)}). "
                "Before version 0.12 these data were SI data regardless of "
                "the network's units. Check your data or declare their units "
                "with input_units and output_unit.",
                UserWarning, stacklevel=2
            )

    def _inputs_to_data(self, inputs):
        axes, _ = self._axis_units()
        return [
            self._units.convert(x, *axis) if axis is not None else x
            for x, axis in zip(inputs, axes)
        ]

    def _input_to_data(self, index, value):
        axes, _ = self._axis_units()
        axis = axes[index]
        if axis is None:
            return value
        return self._units.convert(value, *axis)

    def _input_to_SI(self, index, value):
        axes, _ = self._axis_units()
        axis = axes[index]
        if axis is None:
            return value
        return self._units.convert(value, axis[1], axis[0])

    def _output_to_SI(self, value):
        _, axis = self._axis_units()
        if axis is None:
            return value
        return self._units.convert(value, axis[1], axis[0])

    def _output_to_data(self, value):
        _, axis = self._axis_units()
        if axis is None:
            return value
        return self._units.convert(value, *axis)

    def evaluate(self, *inputs):
        """Evaluate the characteristic at SI inputs, returns the SI output."""
        return self._output_to_SI(
            self.char_func.evaluate(*self._inputs_to_data(inputs))
        )

    def check_domain(self, *inputs, label=None):
        """Log warnings if the SI inputs leave the range of the data."""
        self.char_func.get_domain_errors(*self._inputs_to_data(inputs), label)

    def evaluate_x(self, x):
        """First step of a two step map evaluation at SI input, returns the
        second axis and the output arrays in SI units."""
        yarr, zarr = self.char_func.evaluate_x(self._input_to_data(0, x))
        return self._input_to_SI(1, yarr), self._output_to_SI(zarr)

    def evaluate_y(self, y, yarr, zarr):
        """Second step of a two step map evaluation with SI arrays from
        :py:meth:`evaluate_x`, returns the SI output."""
        return self._output_to_SI(self.char_func.evaluate_y(
            self._input_to_data(1, y),
            self._input_to_data(1, yarr),
            self._output_to_data(zarr)
        ))

    def check_domain_x(self, x, label=None):
        yarr = self.char_func.get_domain_errors_x(
            self._input_to_data(0, x), label
        )
        return self._input_to_SI(1, yarr)

    def check_domain_y(self, y, yarr, label=None):
        self.char_func.get_domain_errors_y(
            self._input_to_data(1, y), self._input_to_data(1, yarr), label
        )

    def input_bounds(self, index=0):
        """Lower and upper bound of the data on an input axis in SI units,
        :code:`None` if the characteristic has no bounded data."""
        bounds = self.char_func.input_bounds(index)
        if bounds is None:
            return None
        return tuple(self._input_to_SI(index, b) for b in bounds)

    def _serialize(self, units=None):
        export = {
            "is_set": self.is_set, "param": self.param,
            "char_params": self.char_params, "input_units": self.input_units,
            "output_unit": self.output_unit
        }
        if self.char_func is not None:
            export["char_func"] = self.char_func._serialize()
        return export

    def _repr_compact(self):
        char = None
        if self.char_func is not None:
            char = type(self.char_func).__name__
        return (
            f"{type(self).__name__}(param={self.param!r}, "
            f"char_func={char}, is_set={self.is_set})"
        )


class CharacteristicEquation(_EquationMixin, CharacteristicData):
    """
    Data container for a characteristic function that activates an equation.

    Takes the parameters of :py:class:`CharacteristicData`. :code:`is_set`
    activates the equation, therefore providing a characteristic function
    does not change it: a custom characteristic can be supplied for the
    offdesign calculation while the equation stays inactive in design mode.
    """

    @staticmethod
    def attr():
        return {
            **CharacteristicData.attr(),
            **_equation_attrs(),
        }

    def accept(self, value):
        is_set = self.is_set
        super().accept(value)
        if not isinstance(value, dict) and value is not None:
            self.is_set = is_set


class ReferenceEquation(_EquationMixin, DataContainer):
    """
    Data container for a reference of a connection parameter to the same
    parameter of another connection.

    Parameters
    ----------
    ref : tespy.connections.connection.Ref
        The reference object with target connection, factor and delta.

    is_set : boolean
        Should the equation be applied? default: :code:`False`.

    quantity : str
        Physical quantity of :code:`delta`.
    """

    @staticmethod
    def attr():
        return {
            "ref": None,
            "is_set": False,
            "unit": None,
            "quantity": None,
            "description": None,
            **_equation_attrs(),
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

    def accept(self, value):
        if value is None:
            self.is_set = False
        else:
            self.set_attr(ref=value, is_set=True)

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

    is_set : set
        Which fluid mass fractions have been set, default is_set=set().
    """

    @staticmethod
    def attr():
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
            "sm_col": None,
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


# Solver internal containers: the reference containers of affine variable
# groups, created by the solver, never declared by components or connections.

class ScalarVariable(DataContainer):
    """Reference container of a group of affinely related scalar variables."""

    @staticmethod
    def attr():
        return {
            "_val_SI": 0,
            "_is_var": True,
            "_J_col": None,
            "_d": 1e-4,
            "min_val": None,
            "max_val": None,
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


class VectorVariable(DataContainer):
    """Reference container of a group of affinely related fluid
    compositions."""

    @staticmethod
    def attr():
        return {
            "_val": dict(),
            "_is_var": set(),
            "_J_col": dict(),
            "_d": 1e-4,
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


# Deprecated names of the containers before version 0.12. Code written
# against the old names dispatches with isinstance on them, so the deprecated
# classes recognise instances of their replacement as their own.

class _DeprecatedContainerMeta(type):

    def __instancecheck__(cls, instance):
        return cls._instancecheck(instance)


def _deprecated_container(old_name, new_cls, instancecheck=None, **defaults):
    if instancecheck is None:
        instancecheck = lambda obj: isinstance(obj, new_cls)

    def __init__(self, **kwargs):
        warnings.warn(
            f"{old_name} is deprecated, use {new_cls.__name__} instead.",
            FutureWarning, stacklevel=2
        )
        for key, value in defaults.items():
            kwargs.setdefault(key, value)
        new_cls.__init__(self, **kwargs)

    return _DeprecatedContainerMeta(old_name, (new_cls,), {
        "__init__": __init__,
        "__doc__": f"Deprecated, use :py:class:`{new_cls.__name__}` instead.",
        "__module__": __name__,
        "_instancecheck": staticmethod(instancecheck),
    })


FluidProperties = _deprecated_container(
    "FluidProperties", ConnectionParameter
)
ComponentProperties = _deprecated_container(
    "ComponentProperties", ComponentParameter
)
ComponentArrayProperties = _deprecated_container(
    "ComponentArrayProperties", ArrayParameter
)
SimpleDataContainer = _deprecated_container(
    "SimpleDataContainer", ConfigParameter
)
GroupedComponentProperties = _deprecated_container(
    "GroupedComponentProperties", ParameterGroup,
    lambda obj: isinstance(obj, ParameterGroup) and obj.auto_activate
)
GroupedComponentCharacteristics = _deprecated_container(
    "GroupedComponentCharacteristics", ParameterGroup,
    lambda obj: isinstance(obj, ParameterGroup) and not obj.auto_activate,
    auto_activate=False
)
ComponentCharacteristics = _deprecated_container(
    "ComponentCharacteristics", CharacteristicEquation,
    lambda obj: isinstance(obj, CharacteristicData) and obj.num_inputs == 1
)
ComponentCharacteristicMaps = _deprecated_container(
    "ComponentCharacteristicMaps", CharacteristicEquation,
    lambda obj: isinstance(obj, CharacteristicData) and obj.num_inputs == 2,
    input_quantities=[None, None]
)
ComponentMandatoryConstraints = _deprecated_container(
    "ComponentMandatoryConstraints", MandatoryConstraint
)
ReferencedFluidProperties = _deprecated_container(
    "ReferencedFluidProperties", ReferenceEquation
)
