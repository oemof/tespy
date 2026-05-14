# -*- coding: utf-8

"""Module class component.

All tespy components inherit from this class.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/components.py

SPDX-License-Identifier: MIT
"""

import math
from collections import deque

import numpy as np
import pandas as pd

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
from tespy.tools.global_vars import ERR
from tespy.tools.helpers import TESPyNetworkError
from tespy.tools.helpers import _get_dependents
from tespy.tools.helpers import _get_vector_dependents
from tespy.tools.helpers import _partial_derivative
from tespy.tools.helpers import _partial_derivative_vecvar
from tespy.tools.units import _UNITS
from tespy.tools.units import SI_UNITS


def component_registry(type):
    component_registry.items[type.__name__] = type
    return type


component_registry.items = {}


def _topological_sort(calc_items):
    """Return keys of calc_items in topological order based on calc_deps.

    Raises ValueError if a dependency cycle is detected.
    """
    successors = {k: [] for k in calc_items}
    in_degree = {k: 0 for k in calc_items}
    for k, dc in calc_items.items():
        for dep in dc.calc_deps:
            if dep in calc_items:
                successors[dep].append(k)
                in_degree[k] += 1

    queue = deque(k for k in calc_items if in_degree[k] == 0)
    ordered = []
    while queue:
        k = queue.popleft()
        ordered.append(k)
        for succ in successors[k]:
            in_degree[succ] -= 1
            if in_degree[succ] == 0:
                queue.append(succ)

    if len(ordered) != len(calc_items):
        cycle = set(calc_items) - set(ordered)
        raise ValueError(f"Cycle detected in calc_deps: {cycle}")
    return ordered


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

        if not isinstance(label, str):
            msg = 'Component label must be of type str!'
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
        self._local_connection_design_state = {}

        # add container for components attributes
        self.parameters = self.get_parameters().copy()
        self.__dict__.update(self.parameters)
        self.set_attr(**kwargs)

        self.num_i = len(self.inlets())
        self.num_o = len(self.outlets())
        self.num_power_i = len(self.powerinlets())
        self.num_power_o = len(self.poweroutlets())
        self.num_heat_i = len(self.heatinlets())
        self.num_heat_o = len(self.heatoutlets())

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
        for key, value in kwargs.items():
            if key in self.parameters:
                self._set_parameter(key, value)
            elif key in ('design', 'offdesign'):
                self._set_design_list(key, value)
            elif key in ('local_design', 'local_offdesign',
                         'printout', 'char_warnings', 'bypass'):
                self._set_bool_attr(key, value)
            elif key in ('design_path', 'fkt_group'):
                self._set_path_attr(key, value)
            else:
                msg = f"Component {self.label} has no attribute {key}."
                logger.error(msg)
                raise KeyError(msg)

    def _set_parameter(self, key, value):
        try:
            self.parameters[key].accept(value)
        except TypeError as e:
            msg = (
                f"Bad datatype for keyword argument '{key}' on "
                f"component {self.label}: {e}"
            )
            logger.error(msg)
            raise TypeError(msg) from e

    def _set_design_list(self, key, value):
        if not isinstance(value, list):
            msg = (
                f"Please provide the {key} parameters as list for "
                f"component {self.label}."
            )
            logger.error(msg)
            raise TypeError(msg)
        if not set(value).issubset(self.parameters.keys()):
            keys = ", ".join(self.parameters.keys())
            msg = (
                "Available parameters for (off-)design specification "
                f"of component {self.label} are: {keys}."
            )
            logger.error(msg)
            raise ValueError(msg)
        self.__dict__[key] = value

    def _set_bool_attr(self, key, value):
        if not isinstance(value, bool):
            msg = (
                f"Please provide the {key} parameter as bool for "
                f"component {self.label}."
            )
            logger.error(msg)
            raise TypeError(msg)
        self.__dict__[key] = value
        if key == 'local_offdesign' and not value:
            self._local_connection_design_state = {}

    def _set_path_attr(self, key, value):
        self.__dict__[key] = value
        self.new_design = True
        if key == 'design_path' and value is None:
            self._local_connection_design_state = {}

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
                        data.char_func = CharLine()

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
                        data.char_func = CharMap()

            # grouped component properties
            elif type(data) == dc_gcp:
                is_set = True
                for e in data.elements:
                    if not self.get_attr(e).is_set:
                        is_set = False

                if is_set:
                    if self._mode == "design" and key not in self.offdesign:
                        data.set_attr(is_set=True)
                    elif self._mode == "offdesign" and key not in self.design:
                        data.set_attr(is_set=True)

                elif data.is_set:
                    msg = (
                        f"Not all parameters of the component group {key} "
                        f"of component {self.label} have to been specified, "
                        "the group equation will not be applied. This "
                        "component group uses the following "
                        f"parameters: {', '.join(data.elements)}."
                    )
                    logger.debug(msg)
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
                'func_params': {'variable': 'm'},
                'description': "mass flow equality constraint(s)"
            }),
            'fluid_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': self.num_i,
                'func_params': {'variable': 'fluid'},
                'description': "fluid composition equality constraint(s)"
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

    @staticmethod
    def heatinlets():
        return []

    @staticmethod
    def heatoutlets():
        return []

    @property
    def all_connections(self):
        return self.all_inlets + self.all_outlets

    @property
    def all_inlets(self):
        return self.inl + self.power_inl + self.heat_inl

    @property
    def all_outlets(self):
        return self.outl + self.power_outl + self.heat_outl

    def _validate_connections(self):
        if len(self.outl) != self.num_o:
            msg = (
                f"The component {self.label} is missing "
                f"{self.num_o - len(self.outl)} outgoing connections. "
                "Make sure all outlets are connected and all connections "
                "have been added to the network."
            )
            logger.error(msg)
            raise TESPyNetworkError(msg)
        if len(self.inl) != self.num_i:
            msg = (
                f"The component {self.label} is missing "
                f"{self.num_i - len(self.inl)} incoming connections. "
                "Make sure all inlets are connected and all connections "
                "have been added to the network."
            )
            logger.error(msg)
            raise TESPyNetworkError(msg)

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
                return self.inl[inconn].m.val_SI / self._conn_design(self.inl[inconn], 'm')
            elif param == 'm_out':
                return self.outl[outconn].m.val_SI / self._conn_design(self.outl[outconn], 'm')
            elif param == 'v':
                v = self.inl[inconn].m.val_SI * self.inl[inconn].calc_vol()
                return v / self._conn_design(self.inl[inconn], 'v')
            elif param == 'pr':
                return (
                    (self.outl[outconn].p.val_SI * self._conn_design(self.inl[inconn], 'p'))
                    / (self.inl[inconn].p.val_SI * self._conn_design(self.outl[outconn], 'p'))
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

    def _conn_design(self, conn, param):
        r"""
        Return the design point value of a connection parameter.

        When a component has an individual :code:`design_path` (either because
        it has :code:`local_offdesign=True` in a design-mode solve, or because
        it carries its own :code:`design_path` in an offdesign solve), the
        adjacent connection design values are stored in
        :code:`_local_connection_design_state` during preprocessing.  This
        method returns those local values when available and falls back to the
        connection's own :code:`.design` attribute otherwise.

        Parameters
        ----------
        conn : tespy.connections.connection.Connection
            Adjacent connection object.

        param : str
            Connection parameter name, e.g. :code:`'m'`, :code:`'p'`, :code:`'h'`,
            :code:`'T'`, :code:`'v'`, :code:`'vol'`.

        Returns
        -------
        float
            Design point value in SI units.
        """
        if self._local_connection_design_state:
            local_state = self._local_connection_design_state.get(conn.label)
            if local_state is not None and param in local_state:
                return local_state[param]
        return getattr(conn, param).design

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
        self._mode = mode
        if mode == 'design' or self.local_design:
            self.new_design = True

        for key, dc in self.parameters.items():
            if isinstance(dc, dc_cp):
                if (
                        ((mode == 'offdesign' and not self.local_design) or
                        (mode == 'design' and self.local_offdesign)) and
                        (data.get(key) is not None)
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

    def _calc_pr(self, inconn=0, outconn=0):
        return self.outl[outconn].p.val_SI / self.inl[inconn].p.val_SI

    def _calc_dp(self, inconn=0, outconn=0):
        return self.inl[inconn].p.val_SI - self.outl[outconn].p.val_SI

    def calc_parameters(self):
        r"""Postprocessing parameter calculation.

        Each :py:class:`~tespy.tools.data_containers.ComponentProperties`
        whose :code:`calc` attribute is set is called here in topological
        order (respecting :code:`calc_deps` dependencies).

        .. note::

            Two patterns exist for :code:`calc` methods, and it is important
            to keep them distinct:

            **Pattern A - solver variables only** (p, h, m, fluid composition,
            connection energies E): methods like :py:meth:`_calc_P` read only
            quantities that are unknowns of the solver, therefore these methods
            can also be used in the residual calculations during iterations.

            **Pattern B - derived properties** (T, v, x, ...): methods like
            :code:`_calc_ttd_u` or :code:`_calc_td_log` call helpers such as
            :code:`calc_T()` which rely on values that are computed during
            connection postprocessing (e.g. :code:`connection.T.val_SI`). These
            methods must **only** be called in postprocessing, never during
            iteration, because the derived values are not yet available.

            When adding a new :code:`calc` method, choose Pattern A if the
            result depends solely on solver variables; choose Pattern B
            otherwise, and make sure no caller invokes it during iteration.
        """
        calc_items = {
            k: dc for k, dc in self.parameters.items()
            if isinstance(dc, dc_cp) and dc.calc is not None
        }
        for k in _topological_sort(calc_items):
            dc = calc_items[k]
            dc.val_SI = dc.calc(**dc.calc_params)

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

    def _calc_zeta(self, inconn=0, outconn=0):
        i, o = self.inl[inconn], self.outl[outconn]

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
        The zeta value is calculated on the basis of a given pressure loss at
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
