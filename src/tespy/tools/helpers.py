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
from tespy.tools.global_vars import ERR
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
    aliases = [alias.replace(' ', '') for alias in CP.get_aliases(fluid)]
    return any(alias in fluid_list for alias in aliases)


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


def latex_unit(unit):
    r"""
    Convert unit to LaTeX.

    Parameters
    ----------
    unit : str
        Value of unit for input, e.g. :code:`m3 / kg`.

    Returns
    -------
    unit : str
        Value of unit for output, e.g. :code:`$\unitfrac{m3}{kg}$`.
    """
    if '/' in unit:
        numerator = unit.split('/')[0].replace(' ', '')
        denominator = unit.split('/')[1].replace(' ', '')
        return r'$\unitfrac[]{' + numerator + '}{' + denominator + '}$'
    else:
        if unit == 'C' or unit == 'F':
            unit = r'^\circ ' + unit
        return r'$\unit[]{' + unit + '}$'


class UserDefinedEquation:

    def __init__(self, label, func, deriv, conns, params={},
                 latex={}, structure_matrix=None):
        r"""
        A UserDefinedEquation allows use of generic user specified equations.

        Parameters
        ----------
        label : str
            Label of the user defined function.

        func : function
            Equation to evaluate.

        deriv : function
            Partial derivatives of the equation.

        conns : list
            List of connections used by the function.

        params : dict
            Dictionary containing keyword arguments required by the function
            and/or derivative.

        latex : dict
            Dictionary holding LaTeX string of the equation as well as
            CharLine and CharMap instances applied in the equation for the
            automatic model documentation module.

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
        >>> nw = Network(p_unit='bar', T_unit='C')
        >>> nw.set_attr(iterinfo=False)
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
        pressure, enthalpy and fluid composition. In this case, the derivatives
        to the mass flow, pressure and enthalpy of the inflow as well as the
        derivatives to the pressure and enthalpy of the outflow will be
        required. You have to define a function placing the derivatives in the
        Jacobian matrix. The Jacobian is a dictionary containing tuples as keys
        with the derivative as their value. The tuples indicate the equation
        number (always 0 for user defined equations, since there is only a
        single equation) and the position of the variable in the system matrix.
        The position of the variables is stored in the :code:`J_col` attribute.
        Before calculating and placing a result in the Jacobian, you have to
        make sure, that the variable you want to calculate the partial
        derivative for is actually a variable. For example, in case you
        specified a value for the mass flow, it will not be part of the
        variables' space, since it has a constant value, and thus, no derivate
        needs to be calculated. You can use the :code:`is_var` keyword to check,
        whether a mass flow, pressure or enthalpy is actually variable.

        We can calculate the derivatives numerically, if an easy analytical
        solution is not available. Simply use the :code:`numeric_deriv` method
        passing the variable ('m', 'p', 'h', 'fluid') as well as the
        connection.

        >>> def myjacobian(ude):
        ...     c0 = ude.conns[0]
        ...     c1 = ude.conns[1]
        ...     ude._partial_derivative(c0.m)
        ...     ude._partial_derivative(c0.p)
        ...     ude._partial_derivative(c0.h)
        ...     ude._partial_derivative(c1.p)
        ...     ude._partial_derivative(c1.h)

        After that, we only need to th specify the characteristic line we want
        out temperature drop to follow as well as create the
        UserDefinedEquation instance and add it to the network. Its equation
        is automatically applied. We apply extrapolation for the characteristic
        line as it helps with convergence, in case a paramter

        >>> char = CharLine(
        ...    x=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0],
        ...    y=[17, 12, 9, 6.5, 4.5, 3, 2, 1.5, 1.25, 1.125, 1.1, 1.05],
        ...    extrapolate=True)
        >>> my_ude = UserDefinedEquation(
        ...    'myudelabel', myfunc, myjacobian, [inflow, outflow],
        ...    params={'char': char})
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
        """
        if isinstance(label, str):
            self.label = label
        else:
            msg = 'Label of UserDefinedEquation object must be of type String.'
            logger.error(msg)
            raise TypeError(msg)

        if isinstance(conns, list):
            self.conns = conns
        else:
            msg = (
                'Parameter conns must be a list of '
                'tespy.connections.connection.Connection objects.'
            )
            logger.error(msg)
            raise TypeError(msg)

        self.func = func
        self.deriv = deriv
        if structure_matrix is None:
            self.structure_matrix = {}
        else:
            self.structure_matrix = structure_matrix

        self.right_hand_side = {}

        if isinstance(params, dict):
            self.params = params
        else:
            msg = 'The parameter params must be passed as dictionary.'
            logger.error(msg)
            raise TypeError(msg)

        self.latex = {
            'equation': r'\text{equation string not available}',
            'lines': [],
            'maps': []
        }
        if isinstance(latex, dict):
            self.latex.update(latex)
        else:
            msg = 'The parameter latex must be passed as dictionary.'
            logger.error(msg)
            raise TypeError(msg)

    def solve(self):
        self.residual = self.func(self)
        self.deriv(self)
        logger.error(self.jacobian)

    def get_structure_matrix(self):
        return

    def _partial_derivative(self, var, value=None):
        if value is None:
            value = self.func
        result = _partial_derivative(var, value, None, ude=self)
        if result is not None:
            self.jacobian[var.J_col] = result

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
        tol_rel=ERR, tol_abs=ERR, tol_mode="rel", **function_kwargs
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
