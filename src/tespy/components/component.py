# -*- coding: utf-8

"""Module class component.

All tespy components inherit from this class.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/components.py

SPDX-License-Identifier: MIT
"""

import math

import numpy as np
import pandas as pd
import pint

from tespy.tools import logger
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap
from tespy.tools.characteristics import load_default_char as ldc
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentCharacteristics as dc_gcc
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.global_vars import ERR
from tespy.tools.helpers import _get_dependents
from tespy.tools.helpers import _get_vector_dependents
from tespy.tools.helpers import _partial_derivative
from tespy.tools.helpers import _partial_derivative_vecvar
from tespy.tools.helpers import bus_char_derivative
from tespy.tools.helpers import bus_char_evaluation
from tespy.tools.helpers import newton_with_kwargs
from tespy.tools.units import _UNITS
from tespy.tools.units import SI_UNITS


def component_registry(type):
    component_registry.items[type.__name__] = type
    return type


component_registry.items = {}


@component_registry
class Component:
    r"""
    Class Component is the base class of all TESPy components.

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    design_path : str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings : boolean
        Ignore warnings on default characteristics usage for this component.

    printout : boolean
        Include this component in the network's results printout.

    **kwargs :
        See the class documentation of desired component for available
        keywords.

    Note
    ----
    The initialisation method (__init__), setter method (set_attr) and getter
    method (get_attr) are used for instances of class component and its
    children.

    Allowed keywords in kwargs are 'design_path', 'design' and 'offdesign'.
    Additional keywords depend on the type of component you want to create.

    Example
    -------
    Basic example for a setting up a
    :py:class:`tespy.components.component.Component` object. This example does
    not run a tespy calculation.

    >>> from tespy.components.component import Component
    >>> comp = Component('myComponent')
    >>> type(comp)
    <class 'tespy.components.component.Component'>
    """

    def __init__(self, label, **kwargs):

        # check if components label is of type str and for prohibited chars
        _forbidden = [';', ',', '.']
        if not isinstance(label, str):
            msg = 'Component label must be of type str!'
            logger.error(msg)
            raise ValueError(msg)

        elif any([True for x in _forbidden if x in label]):
            msg = (
                f"You cannot use any of {', '.join(_forbidden)} in a "
                f"component label ({self.__class__.__name__})"
            )
            logger.error(msg)
            raise ValueError(msg)

        else:
            self.label = label

        # defaults
        self.new_design = True
        self.design_path = None
        self.design = []
        self.offdesign = []
        self.local_design = False
        self.local_offdesign = False
        self.char_warnings = True
        self.printout = True
        self.bypass = False
        self.fkt_group = self.label

        # add container for components attributes
        self.parameters = self.get_parameters().copy()
        self.__dict__.update(self.parameters)
        self.set_attr(**kwargs)

    def set_attr(self, **kwargs):
        r"""
        Set, reset or unset attributes of a component for provided arguments.

        Parameters
        ----------
        design : list
            List containing design parameters (stated as String).

        offdesign : list
            List containing offdesign parameters (stated as String).

        design_path: str
            Path to the components design case.

        **kwargs :
            See the class documentation of desired component for available
            keywords.

        Note
        ----
        Allowed keywords in kwargs are obtained from class documentation as all
        components share the
        :py:meth:`tespy.components.component.Component.set_attr` method.
        """
        # set specified values
        for key in kwargs:
            if key in self.parameters:
                data = self.get_attr(key)
                if kwargs[key] is None:
                    data.set_attr(is_set=False)
                    try:
                        data.set_attr(is_var=False)
                    except KeyError:
                        pass
                    continue


                is_numeric = False
                is_quantity = False

                if isinstance(kwargs[key], pint.Quantity):
                    is_quantity = True
                else:
                    try:
                        float(kwargs[key])
                        is_numeric = True
                    except (TypeError, ValueError):
                        pass

                # dict specification
                if (isinstance(kwargs[key], dict) and
                        not isinstance(data, dc_simple)):
                    data.set_attr(**kwargs[key])

                # value specification for component properties
                elif isinstance(data, dc_cp) or isinstance(data, dc_simple):
                    if is_numeric or is_quantity:
                        data.set_attr(val=kwargs[key], is_set=True)
                        if isinstance(data, dc_cp):
                            data.set_attr(is_var=False)

                    elif isinstance(data, dc_simple):
                        data.set_attr(val=kwargs[key], is_set=True)

                    elif kwargs[key] == 'var' and isinstance(data, dc_cp):
                        data.set_attr(is_set=True, is_var=True)

                    # invalid datatype for keyword
                    else:
                        msg = (
                            f"Bad datatype for keyword argument {key} for "
                            f"component {self.label}."
                        )
                        logger.error(msg)
                        raise TypeError(msg)

                elif isinstance(data, dc_cc) or isinstance(data, dc_cm):
                    # value specification for characteristics
                    if (isinstance(kwargs[key], CharLine) or
                            isinstance(kwargs[key], CharMap)):
                        data.char_func = kwargs[key]

                    # invalid datatype for keyword
                    else:
                        msg = (
                            f"Bad datatype for keyword argument {key} for "
                            f"component {self.label}."
                        )
                        logger.error(msg)
                        raise TypeError(msg)

            elif key in ['design', 'offdesign']:
                if not isinstance(kwargs[key], list):
                    msg = (
                        f"Please provide the {key} parameters as list for "
                        f"component {self.label}."
                    )
                    logger.error(msg)
                    raise TypeError(msg)
                if set(kwargs[key]).issubset(list(self.parameters.keys())):
                    self.__dict__.update({key: kwargs[key]})

                else:
                    keys = ", ".join(self.parameters.keys())
                    msg = (
                        "Available parameters for (off-)design specification "
                        f"of component {self.label} are: {keys}."
                    )
                    logger.error(msg)
                    raise ValueError(msg)

            elif key in ['local_design', 'local_offdesign',
                         'printout', 'char_warnings', 'bypass']:
                if not isinstance(kwargs[key], bool):
                    msg = (
                        f"Please provide the {key} parameters as bool for "
                        f"component {self.label}."
                    )
                    logger.error(msg)
                    raise TypeError(msg)

                else:
                    self.__dict__.update({key: kwargs[key]})

            elif key == 'design_path' or key == 'fkt_group':
                self.__dict__.update({key: kwargs[key]})

                self.new_design = True

            # invalid keyword
            else:
                msg = f"Component {self.label} has no attribute {key}."
                logger.error(msg)
                raise KeyError(msg)

    def get_attr(self, key):
        r"""
        Get the value of a component's attribute.

        Parameters
        ----------
        key : str
            The attribute you want to retrieve.

        Returns
        -------
        out :
            Value of specified attribute.
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            msg = f"Component {self.label} has no attribute {key}."
            logger.error(msg)
            raise KeyError(msg)

    def _serialize(self):
        export = {}
        for k in self._serializable():
            export.update({k: self.get_attr(k)})
        for k in self.parameters:
            data = self.get_attr(k)
            export.update({k: data._serialize()})

        return {self.label: export}

    @staticmethod
    def _serializable():
        return [
            "design", "offdesign", "local_design", "local_offdesign",
            "design_path", "printout", "fkt_group", "char_warnings", "bypass"
        ]

    def _get_result_attributes(self):
        return [
            key for key, p in self.parameters.items() if isinstance(p, dc_cp)
        ]

    def collect_results(self):
        result = {}
        for key in self._get_result_attributes():
            p = self.get_attr(key)
            if (p.func is not None or (p.func is None and p.is_set) or
                    p.is_result or p.structure_matrix is not None):
                result[key] = p.val
            else:
                result[key] = np.nan

            result[f"{key}_unit"] = p.unit
        return pd.Series(result)

    def propagate_wrapper_to_target(self, branch):
        inconn = branch["connections"][-1]
        conn_idx = self.inl.index(inconn)
        outconn = self.outl[conn_idx]

        branch["connections"] += [outconn]
        branch["components"] += [self]

        outconn.target.propagate_wrapper_to_target(branch)

    def get_variables(self):
        variables = {}
        for key, data in self.parameters.items():
            if isinstance(data, dc_cp):
                if data.is_var:
                    variables.update({key: data})

        return variables

    def _preprocess(self, row_idx):
        r"""
        Perform component initialization in network preprocessing.

        Parameters
        ----------
        nw : tespy.networks.network.Network
            Network this component is integrated in.
        """
        self.num_vars = 0

        if self.bypass:
            self.constraints = self.get_bypass_constraints().copy()
        else:
            self.constraints = self.get_mandatory_constraints().copy()

        self.parameters

        self.prop_specifications = {}
        self.var_specifications = {}
        self.group_specifications = {}
        self.char_specifications = {}
        # ???
        self.__dict__.update(self.constraints)

        self._structure_matrix = {}
        self._rhs = {}
        self._equation_set_lookup = {}

        sum_eq = 0

        for name, constraint in self.constraints.items():
            for i in range(sum_eq, sum_eq + constraint.num_eq_sets):
                self._rhs[i + row_idx] = 0
                self._equation_set_lookup[i + row_idx] = name

            if constraint.structure_matrix is not None:
                constraint.structure_matrix(row_idx + sum_eq, **constraint.func_params)

            sum_eq += constraint.num_eq_sets

        if not self.bypass:
            sum_eq = self._setup_user_imposed_constraints(row_idx, sum_eq)

        self.num_eq = sum_eq

    def _setup_user_imposed_constraints(self, row_idx, sum_eq):
        for key, data in self.parameters.items():
            if isinstance(data, dc_cp):
                if data.is_var:
                    self.num_vars += 1
                    data._potential_var = True
                else:
                    data._potential_var = False
                self.prop_specifications[key] = data.is_set
                self.var_specifications[key] = data.is_var

            # component characteristics
            elif isinstance(data, dc_cc):
                if data.func is not None:
                    self.char_specifications[key] = data.is_set
                if data.char_func is None:
                    try:
                        data.char_func = ldc(
                            self.__class__.__name__, key, 'DEFAULT', CharLine
                        )
                    except KeyError:
                        data.char_func = CharLine(x=[0, 1], y=[1, 1])

            # component characteristics
            elif isinstance(data, dc_cm):
                if data.func is not None:
                    self.char_specifications[key] = data.is_set
                if data.char_func is None:
                    try:
                        data.char_func = ldc(
                            self.__class__.__name__, key, 'DEFAULT', CharMap
                        )
                    except KeyError:
                        data.char_func = CharLine(x=[0, 1], y=[1, 1])

            # grouped component properties
            elif type(data) == dc_gcp:
                is_set = True
                for e in data.elements:
                    if not self.get_attr(e).is_set:
                        is_set = False

                if is_set:
                    data.set_attr(is_set=True)
                elif data.is_set:
                    msg = (
                        'All parameters of the component group have to be '
                        'specified! This component group uses the following '
                        f'parameters: {", ".join(data.elements)} at '
                        f'{self.label}. Group will be set to False.'
                    )
                    logger.warning(msg)
                    data.set_attr(is_set=False)
                else:
                    data.set_attr(is_set=False)
                self.group_specifications[key] = data.is_set

            # grouped component characteristics
            elif type(data) == dc_gcc:
                self.group_specifications[key] = data.is_set
            # add equations to structure matrix
            structure_matrix = data.structure_matrix is not None
            func = data.func is not None
            if data.is_set and (structure_matrix or func):
                for i in range(sum_eq, sum_eq + data.num_eq_sets):
                    self._rhs[i + row_idx] = 0
                    self._equation_set_lookup[i + row_idx] = key

                if structure_matrix:
                    data.structure_matrix(row_idx + sum_eq, **data.func_params)

                sum_eq += data.num_eq_sets

        return sum_eq

    def _update_num_eq(self):
        pass

    def _check_dependents_implemented(self, deriv, dependents):
        if deriv is None and len(dependents) > 1:
            msg = (
                "Retrieving the derivatives of component parameters "
                "associated with more than one equation is not yet "
                "supported. For these equations, you have to implement "
                "a separate derivate calculation method yourself and "
                "specify it in the component's parameter dictionaries."
            )
            raise NotImplementedError(msg)

    def _assign_dependents_and_eq_mapping(self, value, data, eq_dict, eq_counter):
        if data.dependents is None:
            scalar_dependents = [[] for _ in range(data.num_eq)]
            vector_dependents = [{} for _ in range(data.num_eq)]
        else:
            dependents = data.dependents(**data.func_params)
            if type(dependents) == list:
                scalar_dependents = _get_dependents(dependents)
                vector_dependents = [{} for _ in range(data.num_eq)]
            else:
                scalar_dependents = _get_dependents(dependents["scalars"])
                vector_dependents = _get_vector_dependents(dependents["vectors"])

                # this is a temporary fix
                if len(vector_dependents) < data.num_eq:
                    vector_dependents = [{} for _ in range(data.num_eq)]

            self._check_dependents_implemented(data.deriv, scalar_dependents)

        eq_dict[value]._scalar_dependents = scalar_dependents
        eq_dict[value]._vector_dependents = vector_dependents
        eq_dict[value]._first_eq_index = eq_counter

        for i in range(data.num_eq):
            self._equation_lookup[eq_counter + i] = (value, i)
            self._equation_scalar_dependents_lookup[eq_counter + i] = scalar_dependents[i]
            self._equation_vector_dependents_lookup[eq_counter + i] = vector_dependents[i]

    def _prepare_for_solver(self, system_dependencies, eq_counter):
        r"""
        Perform component initialization in network preprocessing.

        Parameters
        ----------
        nw : tespy.networks.network.Network
            Network this component is integrated in.
        """
        self.num_eq = 0
        self.it = 0
        self.equations = {}
        self._equation_lookup = {}
        self._equation_scalar_dependents_lookup = {}
        self._equation_vector_dependents_lookup = {}

        self._update_num_eq()

        for eq_num, value in self._equation_set_lookup.items():
            if eq_num in system_dependencies:
                continue

            eq_dict = self.equations
            if value in self.constraints:
                data = self.constraints[value]
            elif value in self.parameters:
                data = self.parameters[value]

            if data.num_eq == 0:
                continue

            if value not in eq_dict:
                eq_dict.update({value: data})
                self._assign_dependents_and_eq_mapping(
                    value, data, eq_dict, eq_counter
                )
                self.num_eq += data.num_eq
                eq_counter += data.num_eq

        self.jacobian = {}
        self.residual = {}

        # this could in principle apply for all equations!
        for constraint in self.equations.values():
            eq_num = constraint._first_eq_index
            if constraint.constant_deriv:
                constraint.deriv(None, eq_num)

        return eq_counter

    def get_parameters(self):
        return {}

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': self.num_i,
                'func_params': {'variable': 'm'}
            }),
            'fluid_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': self.num_i,
                'func_params': {'variable': 'fluid'}
            })
        }

    def get_bypass_constraints(self):
        msg = (
            f"The component {self.label} of type {self.__class__.__name__} "
            "does not have bypassing functionality yet."
        )
        logger.exception(msg)
        raise NotImplementedError(msg)

    @staticmethod
    def inlets():
        return []

    @staticmethod
    def outlets():
        return []

    @staticmethod
    def powerinlets():
        return []

    @staticmethod
    def poweroutlets():
        return []

    def _partial_derivative(self, var, eq_num, value, increment_filter=None, **kwargs):
        result = _partial_derivative(var, value, increment_filter, **kwargs)
        if result is not None:
            self.jacobian[eq_num, var.J_col] = result

    def _partial_derivative_fluid(self, var, eq_num, value, dx, increment_filter=None, **kwargs):
        result = _partial_derivative_vecvar(var, value, dx, increment_filter, **kwargs)
        if result is not None:
            self.jacobian[eq_num, var.J_col[dx]] = result

    def _add_missing_fluids(self, connections):
        return []

    def get_char_expr(self, param, type='rel', inconn=0, outconn=0):
        r"""
        Generic method to access characteristic function parameters.

        Parameters
        ----------
        param : str
            Parameter for characteristic function evaluation.

        type : str
            Type of expression:

            - :code:`rel`: relative to design value
            - :code:`abs`: absolute value

        inconn : int
            Index of inlet connection.

        outconn : int
            Index of outlet connection.

        Returns
        -------
        expr : float
            Value of expression
        """
        if type == 'rel':
            if param == 'm':
                return self.inl[inconn].m.val_SI / self.inl[inconn].m.design
            elif param == 'm_out':
                return self.outl[outconn].m.val_SI / self.outl[outconn].m.design
            elif param == 'v':
                v = self.inl[inconn].m.val_SI * self.inl[inconn].calc_vol()
                return v / self.inl[inconn].v.design
            elif param == 'pr':
                return (
                    (self.outl[outconn].p.val_SI * self.inl[inconn].p.design)
                    / (self.inl[inconn].p.val_SI * self.outl[outconn].p.design)
                )
            else:
                msg = (
                    f"The parameter {param} is not available for "
                    "characteristic function evaluation."
                )
                logger.error(msg)
                raise ValueError(msg)
        else:
            if param == 'm':
                return self.inl[inconn].m.val_SI
            elif param == 'm_out':
                return self.outl[outconn].m.val_SI
            elif param == 'v':
                return self.inl[inconn].m.val_SI * self.inl[inconn].calc_vol()
            elif param == 'pr':
                return self.outl[outconn].p.val_SI / self.inl[inconn].p.val_SI
            else:
                return False

    def bus_func(self, bus):
        r"""
        Base method for calculation of the value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        residual : float
            Residual value of bus equation.
        """
        return 0

    def calc_bus_expr(self, bus):
        r"""
        Return the busses' characteristic line input expression.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            Bus to calculate the characteristic function expression for.

        Returns
        -------
        expr : float
            Ratio of power to power design depending on the bus base
            specification.
        """
        b = bus.comps.loc[self]
        if np.isnan(b['P_ref']) or b['P_ref'] == 0:
            return 1
        else:
            comp_val = self.bus_func(b)
            if b['base'] == 'component':
                return abs(comp_val / b['P_ref'])
            else:
                kwargs = {
                    "function": bus_char_evaluation,
                    "parameter": "bus_value",
                    "component_value": comp_val,
                    "reference_value": b["P_ref"],
                    "char_func": b["char"]
                }
                bus_value = newton_with_kwargs(
                    derivative=bus_char_derivative,
                    target_value=0,
                    val0=b['P_ref'],
                    valmin=-1e15,
                    valmax=1e15,
                    **kwargs
                )
                return bus_value / b['P_ref']

    def calc_bus_efficiency(self, bus):
        r"""
        Return the busses' efficiency.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            Bus to calculate the efficiency value on.

        Returns
        -------
        efficiency : float
            Efficiency value of the bus.

            .. math::

                \eta_\mathrm{bus} = \begin{cases}
                \eta\left(
                \frac{\dot{E}_\mathrm{bus}}{\dot{E}_\mathrm{bus,ref}}\right) &
                \text{bus base = 'bus'}\\
                \eta\left(
                \frac{\dot{E}_\mathrm{component}}
                {\dot{E}_\mathrm{component,ref}}\right) &
                \text{bus base = 'component'}
                \end{cases}

        Note
        ----
        If the base value of the bus is the bus value itself, a newton
        iteration is used to find the bus value satisfying the corresponding
        equation (case 1).
        """
        return bus.comps.loc[self, 'char'].evaluate(self.calc_bus_expr(bus))

    def calc_bus_value(self, bus):
        r"""
        Return the busses' value of the component's energy transfer.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            Bus to calculate energy transfer on.

        Returns
        -------
        bus_value : float
            Value of the energy transfer on the specified bus.

            .. math::

                \dot{E}_\mathrm{bus} = \begin{cases}
                \frac{\dot{E}_\mathrm{component}}{f\left(
                \frac{\dot{E}_\mathrm{bus}}{\dot{E}_\mathrm{bus,ref}}\right)} &
                \text{bus base = 'bus'}\\
                \dot{E}_\mathrm{component} \cdot f\left(
                \frac{\dot{E}_\mathrm{component}}
                {\dot{E}_\mathrm{component,ref}}\right) &
                \text{bus base = 'component'}
                \end{cases}

        Note
        ----
        If the base value of the bus is the bus value itself, a newton
        iteration is used to find the bus value satisfying the corresponding
        equation (case 1).
        """
        b = bus.comps.loc[self]
        comp_val = self.bus_func(b)
        expr = self.calc_bus_expr(bus)
        if b['base'] == 'component':
            return comp_val * b['char'].evaluate(expr)
        else:
            return comp_val / b['char'].evaluate(expr)

    def initialise_source(self, c, key):
        r"""
        Return a starting value for pressure and enthalpy at outlet.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                0 & \text{key = 'p'}\\
                0 & \text{key = 'h'}
                \end{cases}
        """
        return 0

    def initialise_target(self, c, key):
        r"""
        Return a starting value for pressure and enthalpy at inlet.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                0 & \text{key = 'p'}\\
                0 & \text{key = 'h'}
                \end{cases}
        """
        return 0

    def _set_design_parameters(self, mode, data):
        r"""
        Set or unset design values of component parameters.

        Parameters
        ----------
        mode : str
            Setting component design values for :code:`mode='offdesign'`
            and unsetting them for :code:`mode='design'`.

        df : pandas.core.series.Series
            Series containing the component parameters.
        """
        if mode == 'design' or self.local_design:
            self.new_design = True

        for key, dc in self.parameters.items():
            if isinstance(dc, dc_cp):
                if (
                        ((mode == 'offdesign' and not self.local_design) or
                        (mode == 'design' and self.local_offdesign)) and
                        (data[key] is not None)
                    ):
                    if f"{key}_unit" in data:
                        value = _UNITS.ureg.Quantity(
                            data[key], data[f"{key}_unit"]
                        ).to(SI_UNITS[dc.quantity]).magnitude
                    else:
                        value = data[key]
                    self.get_attr(key).design = float(value)

                else:
                    self.get_attr(key).design = np.nan

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        return

    def check_parameter_bounds(self):
        r"""Check parameter value limits."""
        _no_limit_violated = True
        for p in self.parameters.keys():
            data = self.get_attr(p)
            if isinstance(data, dc_cp):
                if data.val_SI > data.max_val + ERR:
                    msg = (
                        f"Invalid value for {p}: {p} = {data.val_SI} above "
                        f"maximum value ({data.max_val}) at component "
                        f"{self.label}."
                    )
                    logger.warning(msg)
                    _no_limit_violated = False

                elif data.val_SI < data.min_val - ERR:
                    msg = (
                        f"Invalid value for {p}: {p} = {data.val_SI} below "
                        f"minimum value ({data.min_val}) at component "
                        f"{self.label}."
                    )
                    logger.warning(msg)
                    _no_limit_violated = False

            elif isinstance(data, dc_cc) and data.is_set:
                if data.param is not None:
                    expr = self.get_char_expr(data.param, **data.char_params)
                    data.char_func.get_domain_errors(expr, self.label)

            elif isinstance(data, dc_gcc) and data.is_set:
                for char in data.elements:
                    char_data = self.get_attr(char)
                    expr = self.get_char_expr(
                        char_data.param, **char_data.char_params)
                    char_data.char_func.get_domain_errors(expr, self.label)

        return _no_limit_violated

    def convergence_check(self):
        return

    def entropy_balance(self):
        r"""Entropy balance calculation method."""
        return

    def exergy_balance(self, T0):
        r"""
        Exergy balance calculation method.

        Parameters
        ----------
        T0 : float
            Ambient temperature T0 / K.
        """
        self.E_P = np.nan
        self.E_F = np.nan
        self.E_bus = {
            "chemical": np.nan, "physical": np.nan, "massless": np.nan
        }
        self.E_D = np.nan
        self.epsilon = self._calc_epsilon()

    def _calc_epsilon(self):
        if self.E_F == 0:
            return np.nan
        else:
            return self.E_P / self.E_F

    def get_plotting_data(self):
        return

    def pr_structure_matrix(self, k, pr=None, inconn=0, outconn=0):
        r"""
        Create linear relationship between inflow and outflow pressure

        .. math::

            p_{in} \cdot pr = p_{out}

        Parameters
        ----------
        k : int
            equation number in systems of equations

        pr : str
            Component parameter, e.g. :code:`pr1`.

        inconn : int
            Connection index of inlet.

        outconn : int
            Connection index of outlet.
        """
        pr = self.get_attr(pr)
        i = self.inl[inconn]
        o = self.outl[outconn]

        self._structure_matrix[k, i.p.sm_col] = pr.val_SI
        self._structure_matrix[k, o.p.sm_col] = -1

    def variable_equality_structure_matrix(self, k, **kwargs):
        r"""
        Create pairwise linear relationship between two variables :code:`var`
        for all inlets and the respective outlets. This usually is applied to
        mass flow, pressure, enthalpy and fluid composition.

        .. math::

            var_\text{in,i} = var_\text{out,i}

        Parameters
        ----------
        k : int
            equation number in systems of equations

        variable : str
            Connection variable name, e.g. :code:`h`.
        """
        variable = kwargs.get("variable")
        for count, (i, o) in enumerate(zip(self.inl, self.outl)):
            self._structure_matrix[k + count, i.get_attr(variable).sm_col] = 1
            self._structure_matrix[k + count, o.get_attr(variable).sm_col] = -1

    def calc_zeta(self, i, o):
        if abs(i.m.val_SI) <= 1e-4:
            return 0
        else:
            return (
                (i.p.val_SI - o.p.val_SI) * math.pi ** 2
                / (4 * i.m.val_SI ** 2 * (i.vol.val_SI + o.vol.val_SI))
            )

    def zeta_func(self, zeta=None, inconn=0, outconn=0):
        r"""
        Calculate residual value of :math:`\zeta`-function.

        Parameters
        ----------
        zeta : str
            Component parameter to evaluate the zeta_func on, e.g.
            :code:`zeta1`.

        inconn : int
            Connection index of inlet.

        outconn : int
            Connection index of outlet.

        Returns
        -------
        residual : float
            Residual value of function.

            .. math::

                0 = \begin{cases}
                p_{in} - p_{out} & |\dot{m}| < \epsilon \\
                \frac{\zeta}{D^4} - \frac{(p_{in} - p_{out}) \cdot \pi^2}
                {8 \cdot \dot{m}_{in} \cdot |\dot{m}_{in}| \cdot \frac{v_{in} +
                v_{out}}{2}} &
                |\dot{m}| > \epsilon
                \end{cases}

        Note
        ----
        The zeta value is caluclated on the basis of a given pressure loss at
        a given flow rate in the design case. As the cross sectional area A
        will not change, it is possible to handle the equation in this way:

        .. math::

            \frac{\zeta}{D^4} = \frac{\Delta p \cdot \pi^2}
            {8 \cdot \dot{m}^2 \cdot v}
        """
        data = self.get_attr(zeta)
        i = self.inl[inconn]
        o = self.outl[outconn]

        if abs(i.m.val_SI) < 1e-4:
            return i.p.val_SI - o.p.val_SI

        else:
            v_i = i.calc_vol(T0=i.T.val_SI)
            v_o = o.calc_vol(T0=o.T.val_SI)
            return (
                data.val_SI - (i.p.val_SI - o.p.val_SI) * math.pi ** 2
                / (8 * abs(i.m.val_SI) * i.m.val_SI * (v_i + v_o) / 2)
            )

    def zeta_dependents(self, zeta=None, inconn=0, outconn=0):
        return [
            self.inl[inconn].m,
            self.inl[inconn].p,
            self.inl[inconn].h,
            self.outl[outconn].p,
            self.outl[outconn].h,
            self.get_attr(zeta)
        ]

    def dp_structure_matrix(self, k, dp=None, inconn=0, outconn=0):
        r"""
        Create linear relationship between inflow and outflow pressure

        .. math::

            p_{in} - dp = p_{out}

        Parameters
        ----------
        k : int
            equation number in systems of equations

        dp : str
            Component parameter, e.g. :code:`dp1`.

        inconn : int
            Connection index of inlet.

        outconn : int
            Connection index of outlet.
        """
        inlet_conn = self.inl[inconn]
        outlet_conn = self.outl[outconn]
        self._structure_matrix[k, inlet_conn.p.sm_col] = 1
        self._structure_matrix[k, outlet_conn.p.sm_col] = -1
        self._rhs[k] = self.get_attr(dp).val_SI
