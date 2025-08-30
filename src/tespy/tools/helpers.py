# -*- coding: utf-8

"""Module for helper functions used by several other modules.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/helpers.py

SPDX-License-Identifier: MIT
"""

import json
import os
from collections.abc import Mapping
from copy import deepcopy

import CoolProp.CoolProp as CP
import pandas as pd

from tespy import __datapath__
from tespy.tools import logger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.global_vars import ERR
from tespy.tools.global_vars import FLUID_ALIASES
from tespy.tools.global_vars import fluid_property_data


def get_all_subdictionaries(data):
    subdictionaries = []
    for value in data.values():
        if len(value["subbranches"]) == 0:
            subdictionaries.append(
                {k: v for k, v in value.items() if k != "subbranches"}
            )
        else:
            subdictionaries.append(
                {k: v for k, v in value.items() if k != "subbranches"}
            )
            subdictionaries.extend(get_all_subdictionaries(value["subbranches"]))

    return subdictionaries


def get_chem_ex_lib(name):
    """Return a new dictionary by merging two dictionaries recursively."""
    path = os.path.join(__datapath__, "ChemEx", f"{name}.json")
    with open(path, "r") as f:
        return json.load(f)


def fluidalias_in_list(fluid, fluid_list):
    aliases = FLUID_ALIASES.get_fluid(fluid)
    return set(fluid_list) & aliases


def merge_dicts(dict1, dict2):
    """Return a new dictionary by merging two dictionaries recursively."""

    result = deepcopy(dict1)

    for key, value in dict2.items():
        if isinstance(value, Mapping):
            result[key] = merge_dicts(result.get(key, {}), value)
        else:
            result[key] = deepcopy(dict2[key])

    return result


class TESPyNetworkError(Exception):
    """Custom message for network related errors."""

    pass


class TESPyConnectionError(Exception):
    """Custom message for connection related errors."""

    pass


class TESPyComponentError(Exception):
    """Custom message for component related errors."""

    pass


def convert_to_SI(property, value, unit):
    r"""
    Convert a value to its SI value.

    Parameters
    ----------
    property : str
        Fluid property to convert.

    value : float
        Value to convert.

    unit : str
        Unit of the value.

    Returns
    -------
    SI_value : float
        Specified fluid property in SI value.
    """
    if property == 'T':
        converters = fluid_property_data['T']['units'][unit]
        return (value + converters[0]) * converters[1]

    else:
        return value * fluid_property_data[property]['units'][unit]


def convert_from_SI(property, SI_value, unit):
    r"""
    Get a value in the network's unit system from SI value.

    Parameters
    ----------
    property : str
        Fluid property to convert.

    SI_value : float
        SI value to convert.

    unit : str
        Unit of the value.

    Returns
    -------
    value : float
        Specified fluid property value in network's unit system.
    """
    if property == 'T':
        converters = fluid_property_data['T']['units'][unit]
        return SI_value / converters[1] - converters[0]
    else:
        return SI_value / fluid_property_data[property]['units'][unit]


class UserDefinedEquation:

    def __init__(self, label: str, func: callable, dependents:callable, deriv: callable=None, conns: list=[], comps:list=[], params:dict={}):
        r"""
        A UserDefinedEquation allows use of generic user specified equations.

        Parameters
        ----------
        label : str
            Label of the user defined function

        func : function
            Method returning residual value of equation to evaluate

        dependents : function
            Function, that returns the variables the equation depends on

        deriv : function, optional
            Method calculating partial derivatives of the equation, be default
            None

        conns : list
            List of connections used

        comps : list
            List of components used

        params : dict
            Dictionary containing keyword arguments required by the function
            and/or derivative.

        Note
        ----
        The derivative function specification is optional.

        - In case you provide the dependents and no derivative, the equation
          will numerically derived to all of the specified dependents.
        - In case you provide the dependents and the derivative, the derivative
          method will be used to calculate the partial derivatives.

        Example
        -------
        Consider a pipeline transporting hot water with measurement data on
        temperature reduction in the pipeline as function of volumetric flow.
        First, we set up the TESPy model. Additionally, we will import the
        :py:class:`tespy.tools.helpers.UserDefinedEquation` class as well as
        some fluid property functions. We specify fluid property information
        for the inflow and assume that no pressure losses occur in the
        pipeline.

        >>> from tespy.components import Source, Sink, Pipe
        >>> from tespy.networks import Network
        >>> from tespy.connections import Connection
        >>> from tespy.tools import UserDefinedEquation
        >>> from tespy.tools import CharLine
        >>> from tespy.tools.fluid_properties import T_mix_ph, v_mix_ph
        >>> nw = Network()
        >>> nw.set_attr(iterinfo=False)
        >>> nw.units.set_defaults(**{"pressure": "bar", "temperature": "degC"})
        >>> so = Source('source')
        >>> si = Sink('sink')
        >>> pipeline = Pipe('pipeline')
        >>> inflow = Connection(so, 'out1', pipeline, 'in1')
        >>> outflow = Connection(pipeline, 'out1', si, 'in1')
        >>> nw.add_conns(inflow, outflow)
        >>> inflow.set_attr(T=120, p=10, v=1, fluid={'water': 1})
        >>> pipeline.set_attr(pr=1)

        Let's assume, the temperature reduction is measured from inflow and
        outflow temperature. The mathematical description of the relation
        we want the model to follow therefore is:

        .. math::

            0 = T_{in} - T_{out} + f \left( \dot{m}_{in} \cdot v_{in} \right)

        We can define a function, describing exactly that relation using a
        :py:class:`tespy.tools.characteristics.CharLine` object with volumetric
        flow as input values and temperature drop as output values. The
        function should look like this:

        >>> def myfunc(ude):
        ...    char = ude.params['char']
        ...    return (
        ...        ude.conns[0].calc_T() - ude.conns[1].calc_T()
        ...        - char.evaluate(
        ...            ude.conns[0].m.val_SI *
        ...            ude.conns[0].calc_vol()
        ...        )
        ...    )

        The function does only take one parameter, we name it :code:`ude` in
        this case. This parameter will hold all relevant information you pass
        to your UserDefinedEquation later, i.e. a list of the connections
        (:code:`.conns`) required by the UserDefinedEquation as well as
        a dictionary of arbitrary parameters required for your function
        (:code:`.params`). The index of the :code:`.conns` indicates the
        position of the connection in the list of connections required for the
        UserDefinedEquation (see below).

        On top of the equation the solver requires its derivatives with respect
        to all relevant primary variables of the network, which are mass flow
        pressure, enthalpy and fluid composition, potentially also component
        variables. In this case, the derivatives to the mass flow, pressure and
        enthalpy of the inflow as well as the derivatives to the pressure and
        enthalpy of the outflow will be required.

        We can now apply two different approaches:

        - Either we specify the dependents and let the solver automatically
          derive the method, or
        - we specify the derivative function directly ourselves.

        In the first example, we will only specify the dependents. With that,
        we only need to specify the characteristic line we want our temperature
        drop to follow as well as create the UserDefinedEquation instance and
        add it to the network. We apply extrapolation for the characteristic
        line as it helps with convergence.

        >>> def mydependents(ude):
        ...    c0, c1 = ude.conns
        ...    return [c0.m, c0.p, c0.h, c1.p, c1.h]
        >>> char = CharLine(
        ...    x=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0],
        ...    y=[17, 12, 9, 6.5, 4.5, 3, 2, 1.5, 1.25, 1.125, 1.1, 1.05],
        ...    extrapolate=True)
        >>> my_ude = UserDefinedEquation(
        ...    'myudelabel', myfunc,
        ...    conns=[inflow, outflow],
        ...    dependents=mydependents,
        ...    params={'char': char}
        ... )
        >>> nw.add_ude(my_ude)
        >>> nw.solve('design')

        Clearly the result is obvious here as the volumetric flow is exactly
        at one of the supporting points.

        >>> round(inflow.T.val - outflow.T.val, 3)
        1.125

        So now, let's say, we want to calculate the volumetric flow necessary
        to at least maintain a specific temperature at the outflow.

        >>> inflow.set_attr(v=None)
        >>> outflow.set_attr(T=110)
        >>> nw.solve('design')
        >>> round(inflow.v.val, 3)
        0.267

        Or calculate volumetric flow and/or temperature drop corresponding to a
        specified heat loss.

        >>> outflow.set_attr(T=None)
        >>> pipeline.set_attr(Q=-5e6)
        >>> nw.solve('design')
        >>> round(inflow.v.val, 3)
        0.067

        To make the example complete, we will quickly have a look at, how you
        can specify the derivative method. In this case, you have to define a
        function placing the derivatives in the correct locations of the
        Jacobian matrix. The highlevel method :code:`partial_derivative` of the
        class can handle this for your. You have to pass the variable to this
        method, for which you want to calculate the partial derivative. In case
        the value can be determined analytically, you can additionally pass the
        value as second argument to the function. If you do not pass a second
        argument to the function as in the example below, a numerical
        derivative will be calculated.

        >>> def myjacobian(increment_filter, k, dependents=None, ude=None):
        ...     c0 = ude.conns[0]
        ...     c1 = ude.conns[1]
        ...     ude.partial_derivative(c0.m)
        ...     ude.partial_derivative(c0.p)
        ...     ude.partial_derivative(c0.h)
        ...     ude.partial_derivative(c1.p)
        ...     ude.partial_derivative(c1.h)

        >>> my_ude.deriv = myjacobian
        >>> nw.solve('design')

        """
        if isinstance(label, str):
            self.label = label
        else:
            msg = 'Label of UserDefinedEquation object must be of type String.'
            logger.exception(msg)
            raise TypeError(msg)

        self.func = func
        self.deriv = deriv
        self.dependents = dependents
        self.conns = conns
        self.comps = comps
        self.params = params

    def _preprocess(self, row_idx):
        self.num_eq = 0

        self._structure_matrix = {}
        self._rhs = {}
        self._equation_set_lookup = {}

        self._equation_set_lookup[row_idx] = "equation"

        self.num_eq += 1

    def _prepare_for_solver(self, system_dependencies, eq_counter):
        self.num_eq = 0
        self.it = 0
        self.equations = {}
        self._equation_lookup = {}
        self._equation_scalar_dependents_lookup = {}
        self._equation_vector_dependents_lookup = {}

        for eq_num, value in self._equation_set_lookup.items():
            if eq_num in system_dependencies:
                continue

            if value not in self.equations:
                self.equations[value] = dc_cp(
                    func_params={"ude": self},
                    dependents=self.dependents,
                    num_eq_sets=1,
                    func=self.func,
                    deriv=self.deriv
                )
                self._assign_dependents_and_eq_mapping(
                    value, self.equations[value], self.equations, eq_counter
                )
                self.num_eq += 1
                eq_counter += 1

        self.residual = {}
        self.jacobian = {}

        return eq_counter

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

        eq_dict[value]._scalar_dependents = scalar_dependents
        eq_dict[value]._vector_dependents = vector_dependents
        eq_dict[value]._first_eq_index = eq_counter

        for i in range(data.num_eq):
            self._equation_lookup[eq_counter + i] = (value, i)
            self._equation_scalar_dependents_lookup[eq_counter + i] = scalar_dependents[i]
            self._equation_vector_dependents_lookup[eq_counter + i] = vector_dependents[i]

    def get_structure_matrix(self):
        return

    def partial_derivative(self, var, value=None):
        eq_num = self.equations["equation"]._first_eq_index
        if value is None:
            value = self.func
        self._partial_derivative(var, eq_num, value, **{"ude": self})

    def _partial_derivative(self, var, eq_num, value, increment_filter=None, **kwargs):
        result = _partial_derivative(var, value, increment_filter, **kwargs)
        if result is not None:
            self.jacobian[eq_num, var.J_col] = result

    def _partial_derivative_fluid(self, var, eq_num, value, dx, increment_filter=None, **kwargs):
        result = _partial_derivative_vecvar(var, value, dx, increment_filter, **kwargs)
        if result is not None:
            self.jacobian[eq_num, var.J_col[dx]] = result


def solve(obj, increment_filter):
    """
    Solve equations and calculate partial derivatives of a component.

    Parameters
    ----------
    increment_filter : ndarray
        Matrix for filtering non-changing variables.
    """
    for label, data in obj.equations.items():
        eq_num = data._first_eq_index
        _solve_residual(obj, data, eq_num)
        _solve_jacobian(obj, data, increment_filter, eq_num)


def _solve_residual(obj, data, eq_num):
    result = data.func(**data.func_params)
    if isinstance(result, list):
        result = {eq_num + k: value for k, value in enumerate(result)}
    else:
        result = {eq_num: result}

    obj.residual.update(result)


def _solve_jacobian(obj, data, increment_filter, eq_num):
    if data.constant_deriv:
        return
    elif data.deriv is not None:
        data.deriv(
            increment_filter,
            eq_num,
            dependents={
                "scalars": data._scalar_dependents,
                "vectors": data._vector_dependents
            },
            **data.func_params
        )
    else:
        # these can only be parameters with a single equation for now
        for dependent in data._scalar_dependents[0]:
            f = data.func
            obj._partial_derivative(
                dependent, eq_num, f, increment_filter, **data.func_params
            )

        for dependent, dx in data._vector_dependents[0].items():
            f = data.func
            obj._partial_derivative_fluid(
                dependent, eq_num, f, dx, increment_filter, **data.func_params
            )


def _is_variable(var, increment_filter=None):
    if var.is_var:
        if increment_filter is None or not increment_filter[var.J_col]:
            return True
    return False


def _is_variable_vecvar(var, dx, increment_filter=None):
    if dx in var.is_var:
        if increment_filter is None or not increment_filter[var.J_col[dx]]:
            return True
    return False

from tespy.tools.data_containers import ScalarVariable as dc_scavar
from tespy.tools.data_containers import VectorVariable as dc_vecvar


def _partial_derivative(var, value, increment_filter, **kwargs):
    if _is_variable(var, increment_filter):
        if callable(value):
            if type(var) != dc_scavar:
                var = var._reference_container
            return _numeric_deriv(var, value, **kwargs)
        else:
            return value
    else:
        return None


def _partial_derivative_vecvar(var, value, dx, increment_filter, **kwargs):
    if _is_variable_vecvar(var, dx, increment_filter):
        if callable(value):
            if type(var) != dc_vecvar:
                var = var._reference_container
            return _numeric_deriv_vecvar(var, value, dx, **kwargs)
        else:
            return value
    else:
        return None


def _numeric_deriv(variable, func, **kwargs):
    r"""
    Calculate partial derivative of the function func to dx.

    Parameters
    ----------
    variable : object
        Variable container.

    func : function
        Function :math:`f` to calculate the partial derivative for.

    Returns
    -------
    deriv : float/list
        Partial derivative(s) of the function :math:`f` to variable(s)
        :math:`x`.

        .. math::

            \frac{\partial f}{\partial x} = \frac{f(x + d) + f(x - d)}{2 d}
    """
    d = variable.d
    variable.val_SI += d
    exp = func(**kwargs)

    variable.val_SI -= 2 * d
    exp -= func(**kwargs)
    deriv = exp / (2 * d)

    variable.val_SI += d

    return deriv


def _numeric_deriv_vecvar(variable, func, dx, **kwargs):
    r"""
    Calculate partial derivative of the function func to dx.

    Parameters
    ----------
    obj : object
        Instance, which provides the equation to calculate the derivative for.

    func : function
        Function :math:`f` to calculate the partial derivative for.

    dx : str
        Partial derivative.

    conn : tespy.connections.connection.Connection
        Connection to calculate the numeric derivative for.

    Returns
    -------
    deriv : float/list
        Partial derivative(s) of the function :math:`f` to variable(s)
        :math:`x`.

        .. math::

            \frac{\partial f}{\partial x} = \frac{f(x + d) + f(x - d)}{2 d}
    """
    original_vector = variable.val.copy()
    # this is specific to fluids right now (upper limit of 1, lower limit of 0)
    d1 = min(variable.d, 1 - variable.val[dx])
    variable.val[dx] += d1
    exp = func(**kwargs)
    d2 = min(variable.d * 2, variable.val[dx])
    variable.val[dx] -= d2
    exp -= func(**kwargs)

    variable.val = original_vector
    # d2 is the complete delta of the central difference no matter how big
    # d1 is actually
    deriv = exp / d2
    return deriv


def bus_char_evaluation(component_value, char_func, reference_value, bus_value, **kwargs):
    r"""
    Calculate the value of a bus.

    Parameters
    ----------
    comp_value : float
        Value of the energy transfer at the component.

    reference_value : float
        Value of the bus in reference state.

    char_func : tespy.tools.characteristics.char_line
        Characteristic function of the bus.

    Returns
    -------
    residual : float
        Residual of the equation.

        .. math::

            residual = \dot{E}_\mathrm{bus} - \frac{\dot{E}_\mathrm{component}}
            {f\left(\frac{\dot{E}_\mathrm{bus}}
            {\dot{E}_\mathrm{bus,ref}}\right)}
    """
    return bus_value - component_value / char_func.evaluate(
        bus_value / reference_value
    )


def bus_char_derivative(component_value, char_func, reference_value, bus_value, **kwargs):
    """Calculate derivative for bus char evaluation."""
    d = 1e-3
    return (1 - (
        1 / char_func.evaluate((bus_value + d) / reference_value) -
        1 / char_func.evaluate((bus_value - d) / reference_value)
    ) / (2 * d))


def newton_with_kwargs(
        derivative, target_value, val0=300, valmin=70, valmax=3000, max_iter=10,
        tol_rel=ERR, tol_abs=ERR ** 0.5, tol_mode="rel", **function_kwargs
    ):

    # start newton loop
    iteration = 0
    expr = True
    x = val0
    parameter = function_kwargs["parameter"]
    function = function_kwargs["function"]
    relax = 1

    if tol_mode == "rel" and abs(target_value) <= 2 * tol_rel:
        tol_mode = "abs"

    while expr:
        # calculate function residual and new value
        function_kwargs[parameter] = x
        residual = target_value - function(**function_kwargs)
        x += residual / derivative(**function_kwargs) * relax

        # check for value ranges
        if x < valmin:
            x = valmin
        if x > valmax:
            x = valmax

        iteration += 1
        # relaxation to help convergence in case of jumping
        if iteration == 5:
            relax = 0.75
            max_iter = 12

        if iteration > max_iter:
            msg = (
                'The Newton algorithm was not able to find a feasible value '
                f'for function {function}. Current value with x={x} is '
                f'{function(**function_kwargs)}, target value is '
                f'{target_value}, residual is {residual} after {iteration} '
                'iterations.'
            )
            logger.debug(msg)

            break
        if tol_mode == 'abs':
            expr = abs(residual) >= tol_abs
        elif tol_mode == 'rel':
            expr = abs(residual / target_value) >= tol_rel

    return x


def central_difference(function=None, parameter=None, delta=None, **kwargs):
    upper = kwargs.copy()
    upper[parameter] += delta
    lower = kwargs
    lower[parameter] -= delta
    return (function(**upper) - function(**lower)) / (2 * delta)


def get_basic_path():
    """
    Return the basic tespy path and creates it if necessary.

    The basic path is the '.tespy' folder in the $HOME directory.
    """
    basicpath = os.path.join(os.path.expanduser('~'), '.tespy')
    if not os.path.isdir(basicpath):
        os.mkdir(basicpath)
    return basicpath


def extend_basic_path(subfolder):
    """
    Return a path based on the basic tespy path and creates it if necessary.

    The subfolder is the name of the path extension.
    """
    extended_path = os.path.join(get_basic_path(), subfolder)
    if not os.path.isdir(extended_path):
        os.mkdir(extended_path)
    return extended_path


def _get_vector_dependents(variable_list):
    if len(variable_list) == 0:
        return []
    if isinstance(variable_list, list):
        variable_list_new = []
        for _eq_set in variable_list:
            variable_list_new.append({})
            for key, value in _eq_set.items():
                if not value:
                    continue
                k = key._reference_container
                if k not in variable_list_new[-1]:
                    variable_list_new[-1][k] = set()
                variable_list_new[-1][k] |= value
        return variable_list_new
    else:
        return [
            set(var for var in variable_list)
        ]


def _get_dependents(variable_list):
    if isinstance(variable_list[0], list):
        return [set(
            var._reference_container
            for var in sublist if var.is_var
        ) for sublist in variable_list
    ]
    else:
        return [set(
            var._reference_container
            for var in variable_list if var.is_var
        )]


def _nested_dict_of_dataframes_to_dict(dictionary):
    """Transpose a nested dict with dataframes or series in a json style dict

    Parameters
    ----------
    dictionary : dict
        Dictionary of dataframes

    Returns
    -------
    dict
        json style dictionary containing all data from the dataframes
    """
    for key, value in dictionary.items():
        if isinstance(value, dict):
            dictionary[key] = _nested_dict_of_dataframes_to_dict(value)
        else:
            # Otherwise this is assumed a series, then orient is not available
            kwargs = {}
            if isinstance(value, pd.DataFrame):
                kwargs = {"orient": "index"}
            dictionary[key] = value.to_dict(**kwargs)

    return dictionary


def _nested_dict_of_dataframes_to_filetree(dictionary, basepath):
    """Dump a nested dict with dataframes into a folder structrue

    The upper level keys with subdictionaries are folder names, the lower
    level keys (where a dataframe is the value) will be the names of the
    csv files.

    Parameters
    ----------
    dictionary : dict
        Nested dictionary to write to filesystem.
    basepath : str
        path to dump data to
    """
    os.makedirs(basepath, exist_ok=True)
    for key, value in dictionary.items():
        if isinstance(value, dict):
            basepath = os.path.join(basepath, key)
            _nested_dict_of_dataframes_to_filetree(value, basepath)
        else:
            value.to_csv(os.path.join(basepath, f"{key}.csv"))
