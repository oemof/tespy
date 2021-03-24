# -*- coding: utf-8

"""Module for helper functions used by several other modules.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/helpers.py

SPDX-License-Identifier: MIT
"""
import logging
import os

import CoolProp as CP
import numpy as np

from tespy.tools.global_vars import err
from tespy.tools.global_vars import fluid_property_data
from tespy.tools.global_vars import molar_masses

# %%


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
                 latex={}):
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
        >>> from tespy.tools.helpers import UserDefinedEquation
        >>> from tespy.tools import CharLine
        >>> from tespy.tools.fluid_properties import T_mix_ph, v_mix_ph
        >>> nw = Network(fluids=['water'], p_unit='bar', T_unit='C')
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
        ...        T_mix_ph(ude.conns[0].get_flow()) -
        ...        T_mix_ph(ude.conns[1].get_flow()) - char.evaluate(
        ...            ude.conns[0].m.val_SI *
        ...            v_mix_ph(ude.conns[0].get_flow()))
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
        required. Similar to the equation definition, define a function
        returning the corresponding jacobian matrix. The jacobian is a
        dictionary containing numpy arrays for every connection. Therefore
        the first key is the connection you want to calculate the derivative
        for and the second key is the index of the variable in the jacobian.
        The indices correspond to

        - 0: mass flow
        - 1: pressure
        - 2: enthalpy
        - 3 until end (:code:`3:`): fluid composition

        We can calculate the derivatives numerically, if an easy analytical
        solution is not available. Simply use the :code:`numeric_deriv` method
        passing the variable ('m', 'p', 'h', 'fluid') as well as the
        connection's index.

        >>> def myjacobian(ude):
        ...    ude.jacobian[ude.conns[0]][0] = ude.numeric_deriv('m', 0)
        ...    ude.jacobian[ude.conns[0]][1] = ude.numeric_deriv('p', 0)
        ...    ude.jacobian[ude.conns[0]][2] = ude.numeric_deriv('h', 0)
        ...    ude.jacobian[ude.conns[1]][1] = ude.numeric_deriv('p', 1)
        ...    ude.jacobian[ude.conns[1]][2] = ude.numeric_deriv('h', 1)
        ...    return ude.jacobian

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
            logging.error(msg)
            raise TypeError(msg)

        if isinstance(conns, list):
            self.conns = conns
        else:
            msg = (
                'Parameter conns must be a list of '
                'tespy.connections.connection.Connection objects.')
            logging.error(msg)
            raise TypeError(msg)

        self.func = func
        self.deriv = deriv

        if isinstance(params, dict):
            self.params = params
        else:
            msg = 'The parameter params must be passed as dictionary.'
            logging.error(msg)
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
            logging.error(msg)
            raise TypeError(msg)

    def numeric_deriv(self, param, idx):
        r"""
        Calculate partial derivative of the function func to dx numerically.

        Parameters
        ----------
        param : str
            Parameter to calculate partial derivative for.

        idx : int
            Position of the connection to calculate the partial derivative for
            within the list of the connections :code:`conns`.

        Returns
        -------
        deriv : float/list
            Partial derivative(s) of the function :math:`f` to variable(s)
            :math:`x`.

            .. math::

                \frac{\partial f}{\partial x}=\frac{f(x+d)+f(x-d)}{2\cdot d}
        """
        if param == 'fluid':
            d = 1e-5
            deriv = []
            for f in self.conns[0].fluid.val.keys():
                val = self.conns[idx].fluid.val[f]
                if self.conns[idx].fluid.val[f] + d <= 1:
                    self.conns[idx].fluid.val[f] += d
                else:
                    self.conns[idx].fluid.val[f] = 1
                exp = self.func(self)
                if self.conns[idx].fluid.val[f] - 2 * d >= 0:
                    self.conns[idx].fluid.val[f] -= 2 * d
                else:
                    self.conns[idx].fluid.val[f] = 0
                exp -= self.func(self)
                self.conns[idx].fluid.val[f] = val

                deriv += [exp / (2 * d)]

        elif param in ['m', 'p', 'h']:

            if param == 'm':
                d = 1e-4
            else:
                d = 1e-1

            self.conns[idx].get_attr(param).val_SI += d
            exp = self.func(self)
            self.conns[idx].get_attr(param).val_SI -= 2 * d
            exp -= self.func(self)
            self.conns[idx].get_attr(param).val_SI += d

            deriv = exp / (2 * d)

        else:
            msg = (
                'Can only calculate numerical derivative to primary variables.'
                'Please specify "m", "p", "h" or "fluid" as param.')
            logging.error(msg)
            raise ValueError(msg)

        return deriv


def newton(func, deriv, params, y, **kwargs):
    r"""
    Find zero crossings with 1-D newton algorithm.

    Parameters
    ----------
    func : function
        Function to find zero crossing in,
        :math:`0=y-func\left(x,\text{params}\right)`.

    deriv : function
        First derivative of the function.

    params : list
        Additional parameters for function, optional.

    y : float
        Target function value.

    val0 : float
        Starting value, default: val0=300.

    valmin : float
        Lower value boundary, default: valmin=70.

    valmax : float
        Upper value boundary, default: valmax=3000.

    max_iter : int
        Maximum number of iterations, default: max_iter=10.

    tol_rel : float
        Maximum relative tolerance :math:`|\frac{y - f(x)}{f(x)}|`, default
        value: 1e-6.

    tol_abs : float
        Maximum absolute tolerance :math:`|y - f(x)|`, default value: 1e-6.

    tol_mode : str
        Check for relative, absolute or both tolerances:

        - :code:`tol_mode='abs'` (default)
        - :code:`tol_mode='rel'`
        - :code:`tol_mode='both'`

    Returns
    -------
    val : float
        x-value of zero crossing.

    Note
    ----
    Algorithm

    .. math::

        x_{i+1} = x_{i} - \frac{f(x_{i})}{\frac{df}{dx}(x_{i})}\\
        f(x_{i}) \leq \epsilon
    """
    # default valaues
    x = kwargs.get('val0', 300)
    valmin = kwargs.get('valmin', 70)
    valmax = kwargs.get('valmax', 3000)
    max_iter = kwargs.get('max_iter', 10)
    tol_rel = kwargs.get('tol_rel', err)
    tol_abs = kwargs.get('tol_abs', err)
    tol_mode = kwargs.get('tol_mode', 'abs')

    # start newton loop
    expr = True
    i = 0
    while expr:
        # calculate function residual and new value
        res = y - func(params, x)
        x += res / deriv(params, x)

        # check for value ranges
        if x < valmin:
            x = valmin
        if x > valmax:
            x = valmax
        i += 1

        if i > max_iter:
            msg = ('Newton algorithm was not able to find a feasible value '
                   'for function ' + str(func) + '. Current value with x=' +
                   str(x) + ' is ' + str(func(params, x)) +
                   ', target value is ' + str(y) + '.')
            logging.debug(msg)

            break
        if tol_mode == 'abs':
            expr = abs(res) >= tol_abs
        elif tol_mode == 'rel':
            expr = abs(res / y) >= tol_rel
        else:
            expr = abs(res / y) >= tol_rel or abs(res) >= tol_abs

    return x

# %%


def reverse_2d(params, y):
    r"""
    Calculate the residual value of an inverse function.

    Parameters
    ----------
    params : list
        Variable function parameters.

    y : float
        Function value of function :math:`y = f \left( x_1, x_2 \right)`.

    Returns
    -------
    deriv : float
        Residual value of inverse function :math:`x_2 - f\left(x_1, y \right)`.
    """
    func, x1, x2 = params[0], params[1], params[2]
    return x2 - func.ev(x1, y)


def reverse_2d_deriv(params, y):
    r"""
    Calculate derivative of an inverse function.

    Parameters
    ----------
    params : list
        Variable function parameters.

    y : float
        Function value of function :math:`y = f \left( x_1, x_2 \right)`,
        so that :math:`x_2 - f\left(x_1, y \right) = 0`

    Returns
    -------
    deriv : float
        Partial derivative :math:`\frac{\partial f}{\partial y}`.
    """
    func, x1 = params[0], params[1]
    return - func.ev(x1, y, dy=1)


def bus_char_evaluation(params, bus_value):
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
    comp_value = params[0]
    reference_value = params[1]
    char_func = params[2]
    return bus_value - comp_value / char_func.evaluate(
        bus_value / reference_value)


def bus_char_derivative(params, bus_value):
    """Calculate derivative for bus char evaluation."""
    reference_value = params[1]
    char_func = params[2]
    d = 1e-3
    return (1 - (
        1 / char_func.evaluate((bus_value + d) / reference_value) -
        1 / char_func.evaluate((bus_value - d) / reference_value)) / (2 * d))


def molar_mass_flow(flow):
    r"""
    Calculate molar mass flow.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

    Returns
    -------
    m_m : float
        Molar mass flow m_m / (mol/s).

        .. math::

            \dot{m}_\mathrm{m} = \sum_{i} \left( \frac{x_{i}}{M_{i}} \right)
    """
    mm = 0
    for fluid, x in flow.items():
        if x > err:
            mm += x / molar_masses[fluid]
    return mm

# %%


def num_fluids(fluids):
    r"""
    Return number of fluids in fluid mixture.

    Parameters
    ----------
    fluids : dict
        Fluid mass fractions.

    Returns
    -------
    n : int
        Number of fluids in fluid mixture n / 1.

        .. math::

            n = \sum_{i} \left( \begin{cases}
            0 & x_{i} < \epsilon \\
            1 & x_{i} \geq \epsilon
            \end{cases} \right)\;
            \forall i \in \text{network fluids}
    """
    n = 0
    for fluid, x in fluids.items():
        if x > err:
            n += 1

    return n

# %%


def single_fluid(fluids):
    r"""
    Return the name of the pure fluid in a fluid vector.

    Parameters
    ----------
    fluids : dict
        Fluid mass fractions.

    Returns
    -------
    fluid : str
        Name of the single fluid or None in case of mixtures.
    """
    if num_fluids(fluids) == 1:
        for fluid, x in fluids.items():
            if x > err:
                return fluid
    else:
        return None

# %%


def fluid_structure(fluid):
    r"""
    Return the checmical formula of fluid.

    Parameters
    ----------
    fluid : str
        Name of the fluid.

    Returns
    -------
    parts : dict
        Dictionary of the chemical base elements as keys and the number of
        atoms in a molecule as values.

    Example
    -------
    Get the chemical formula of methane.

    >>> from tespy.tools.helpers import fluid_structure
    >>> elements = fluid_structure('methane')
    >>> elements['C'], elements['H']
    (1, 4)
    """
    parts = {}
    for element in CP.CoolProp.get_fluid_param_string(
            fluid, 'formula').split('}'):
        if element != '':
            el = element.split('_{')
            parts[el[0]] = int(el[1])

    return parts

# %%


def darcy_friction_factor(re, ks, d):
    r"""
    Calculate the Darcy friction factor.

    Parameters
    ----------
    re : float
        Reynolds number re / 1.

    ks : float
        Pipe roughness ks / m.

    d : float
        Pipe diameter/characteristic lenght d / m.

    Returns
    -------
    darcy_friction_factor : float
        Darcy friction factor :math:`\lambda` / 1

    Note
    ----
    **Laminar flow** (:math:`re \leq 2320`)

    .. math::

        \lambda = \frac{64}{re}

    **turbulent flow** (:math:`re > 2320`)

    *hydraulically smooth:* :math:`\frac{re \cdot k_{s}}{d} < 65`

    .. math::

        \lambda = \begin{cases}
        0.03164 \cdot re^{-0.25} & re \leq 10^4\\
        \left(1.8 \cdot \log \left(re\right) -1.5 \right)^{-2} &
        10^4 < re < 10^6\\
        solve \left(0 = 2 \cdot \log\left(re \cdot \sqrt{\lambda} \right) -0.8
        - \frac{1}{\sqrt{\lambda}}\right) & re \geq 10^6\\
        \end{cases}

    *transition zone and hydraulically rough:*

    .. math::

        \lambda = solve \left( 0 = 2 \cdot \log \left( \frac{2.51}{re \cdot
        \sqrt{\lambda}} + \frac{k_{s}}{d \cdot 3.71} \right) -
        \frac{1}{\sqrt{\lambda}} \right)

    Reference: :cite:`Nirschl2018`.

    Example
    -------
    Calculate the Darcy friction factor at different hydraulic states.

    >>> from tespy.tools.helpers import darcy_friction_factor
    >>> ks = 5e-5
    >>> d = 0.05
    >>> re_laminar = 2000
    >>> re_turb_smooth = 5000
    >>> re_turb_trans = 70000
    >>> re_high = 1000000
    >>> d_high = 0.8
    >>> re_very_high = 6000000
    >>> d_very_high = 1
    >>> ks_low = 1e-5
    >>> ks_rough = 1e-3
    >>> darcy_friction_factor(re_laminar, ks, d)
    0.032
    >>> round(darcy_friction_factor(re_turb_smooth, ks, d), 3)
    0.038
    >>> round(darcy_friction_factor(re_turb_trans, ks, d), 3)
    0.023
    >>> round(darcy_friction_factor(re_turb_trans, ks_rough, d), 3)
    0.049
    >>> round(darcy_friction_factor(re_high, ks, d_high), 3)
    0.012
    >>> round(darcy_friction_factor(re_very_high, ks_low, d_very_high), 3)
    0.009
    """
    if re <= 2320:
        return 64 / re
    else:
        if re * ks / d < 65:
            if re <= 1e4:
                return blasius(re)
            elif re < 1e6:
                return hanakov(re)
            else:
                l0 = 0.02
                return newton(
                    prandtl_karman, prandtl_karman_derivative, [re],
                    0, val0=l0, valmin=0.00001, valmax=0.2)

        else:
            l0 = 0.002
            return newton(
                colebrook, colebrook_derivative, [re, ks, d], 0,
                val0=l0, valmin=0.0001, valmax=0.2)


def blasius(re):
    """
    Calculate friction coefficient according to Blasius.

    Parameters
    ----------
    re : float
        Reynolds number.

    Returns
    -------
    darcy_friction_factor : float
        Darcy friction factor.
    """
    return 0.3164 * re ** (-0.25)


def hanakov(re):
    """
    Calculate friction coefficient according to Hanakov.

    Parameters
    ----------
    re : float
        Reynolds number.

    Returns
    -------
    darcy_friction_factor : float
        Darcy friction factor.
    """
    return (1.8 * np.log10(re) - 1.5) ** (-2)


def prandtl_karman(params, darcy_friction_factor):
    """
    Calculate friction coefficient according to Prandtl and v. Kármán.

    Applied in smooth conditions.

    Parameters
    ----------
    re : float
        Reynolds number.

    darcy_friction_factor : float
        Darcy friction factor.

    Returns
    -------
    darcy_friction_factor : float
        Darcy friction factor.
    """
    re = params[0]
    return (
        2 * np.log10(re * darcy_friction_factor ** 0.5) - 0.8 -
        1 / darcy_friction_factor ** 0.5)


def prandtl_karman_derivative(params, darcy_friction_factor):
    """Calculate derivative for Prandtl and v. Kármán equation."""
    return (
        1 / (darcy_friction_factor * np.log(10)) +
        1 / 2 * darcy_friction_factor ** (-1.5))


def colebrook(params, darcy_friction_factor):
    """
    Calculate friction coefficient accroding to Colebrook-White equation.

    Applied in transition zone and rough conditions.

    Parameters
    ----------
    re : float
        Reynolds number.

    ks : float
        Equivalent sand roughness.

    d : float
        Pipe's diameter.

    darcy_friction_factor : float
        Darcy friction factor.

    Returns
    -------
    darcy_friction_factor : float
        Darcy friction factor.
    """
    re, ks, d = params[0], params[1], params[2]
    return (
        2 * np.log10(
            2.51 / (re * darcy_friction_factor ** 0.5) + ks / (3.71 * d)) +
        1 / darcy_friction_factor ** 0.5)


def colebrook_derivative(params, darcy_friction_factor):
    """Calculate derivative for Colebrook-White equation."""
    d = 0.001
    return (colebrook(params, darcy_friction_factor + d) -
            colebrook(params, darcy_friction_factor - d)) / (2 * d)

# %%


def modify_path_os(path):
    """
    Modify a path according the os.

    Also detects weather the path specification is absolute or relative and
    adjusts the path respectively.

    Parameters
    ----------
    path : str
        Path to modify.

    Returns
    -------
    path : str
        Modified path.
    """
    if os.name == 'nt':
        # windows
        path = path.replace('/', '\\')
        if path[0] != '\\' and path[1:2] != ':' and path[0] != '.':
            # relative path
            path = '.\\' + path
    else:
        # linux, mac
        path = path.replace('\\', '/')
        if path[0] != '/' and path[0] != '.':
            # relative path
            path = './' + path

    return path

# %%


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
