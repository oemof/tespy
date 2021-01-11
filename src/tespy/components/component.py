# -*- coding: utf-8

"""Module class component.

All tespy components inherit from this class.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/components.py

SPDX-License-Identifier: MIT
"""

import logging
from collections import OrderedDict

import numpy as np

from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap
from tespy.tools.characteristics import CompressorMap
from tespy.tools.characteristics import load_default_char as ldc
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.fluid_properties import v_mix_ph
from tespy.tools.global_vars import err
from tespy.tools.helpers import bus_char_derivative
from tespy.tools.helpers import bus_char_evaluation
from tespy.tools.helpers import newton

# %%


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
        if not isinstance(label, str):
            msg = 'Component label must be of type str!'
            logging.error(msg)
            raise ValueError(msg)

        elif len([x for x in [';', ',', '.'] if x in label]) > 0:
            msg = (
                'You must not use ' + str([';', ',', '.']) + ' in label (' +
                str(self.component()) + ').')
            logging.error(msg)
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

        # add container for components attributes
        self.variables = OrderedDict(self.attr().copy())
        self.__dict__.update(self.variables)
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
            if key in self.variables.keys():
                data = self.get_attr(key)
                if kwargs[key] is None:
                    data.set_attr(is_set=False)
                    try:
                        data.set_attr(is_var=False)
                    except KeyError:
                        pass
                    continue

                try:
                    float(kwargs[key])
                    is_numeric = True
                except (TypeError, ValueError):
                    is_numeric = False

                # dict specification
                if (isinstance(kwargs[key], dict) and
                        not isinstance(data, dc_simple)):
                    data.set_attr(**kwargs[key])

                # value specification for component properties
                elif isinstance(data, dc_cp) or isinstance(data, dc_simple):
                    if is_numeric:
                        if np.isnan(kwargs[key]):
                            data.set_attr(is_set=False)
                            if isinstance(data, dc_cp):
                                data.set_attr(is_var=False)

                        else:
                            data.set_attr(val=kwargs[key], is_set=True)
                            if isinstance(data, dc_cp):
                                data.set_attr(is_var=False)

                    elif (kwargs[key] == 'var' and
                          isinstance(data, dc_cp)):
                        data.set_attr(is_set=True, is_var=True)

                    elif isinstance(data, dc_simple):
                        data.set_attr(val=kwargs[key], is_set=True)

                    # invalid datatype for keyword
                    else:
                        msg = (
                            'Bad datatype for keyword argument ' + key +
                            ' at ' + self.label + '.')
                        logging.error(msg)
                        raise TypeError(msg)

                elif isinstance(data, dc_cc) or isinstance(data, dc_cm):
                    # value specification for characteristics
                    if (isinstance(kwargs[key], CharLine) or
                            isinstance(kwargs[key], CharMap) or
                            isinstance(kwargs[key], CompressorMap)):
                        data.char_func = kwargs[key]

                    # invalid datatype for keyword
                    else:
                        msg = (
                            'Bad datatype for keyword argument ' + key +
                            ' at ' + self.label + '.')
                        logging.error(msg)
                        raise TypeError(msg)

                elif isinstance(data, dc_gcp):
                    # value specification of grouped component parameter method
                    if isinstance(kwargs[key], str):
                        data.method = kwargs[key]

                    # invalid datatype for keyword
                    else:
                        msg = (
                            'Bad datatype for keyword argument ' + key +
                            ' at ' + self.label + '.')
                        logging.error(msg)
                        raise TypeError(msg)

            elif key in ['design', 'offdesign']:
                if not isinstance(kwargs[key], list):
                    msg = (
                        'Please provide the ' + key + ' parameters as list '
                        'at ' + self.label + '.')
                    logging.error(msg)
                    raise TypeError(msg)
                if set(kwargs[key]).issubset(list(self.variables.keys())):
                    self.__dict__.update({key: kwargs[key]})

                else:
                    msg = (
                        'Available parameters for (off-)design specification '
                        'are: ' + str(list(self.variables.keys())) + ' at ' +
                        self.label + '.')
                    logging.error(msg)
                    raise ValueError(msg)

            elif key in ['local_design', 'local_offdesign',
                         'printout', 'char_warnings']:
                if not isinstance(kwargs[key], bool):
                    msg = (
                        'Please provide the parameter ' + key + ' as boolean '
                        'at component ' + self.label + '.')
                    logging.error(msg)
                    raise TypeError(msg)

                else:
                    self.__dict__.update({key: kwargs[key]})

            elif key == 'design_path':
                if isinstance(kwargs[key], str):
                    self.__dict__.update({key: kwargs[key]})
                elif kwargs[key] is None:
                    self.design_path = None
                elif np.isnan(kwargs[key]):
                    self.design_path = None
                else:
                    msg = (
                        'Please provide the design_path parameter as string. '
                        'For unsetting use np.nan or None.')
                    logging.error(msg)
                    raise TypeError(msg)

                self.new_design = True

            # invalid keyword
            else:
                msg = (
                    'Component ' + self.label + ' has no attribute ' +
                    str(key) + '.')
                logging.error(msg)
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
            msg = ('Component ' + self.label + ' has no attribute \"' +
                   key + '\".')
            logging.error(msg)
            raise KeyError(msg)

    def comp_init(self, nw, num_eq=0):
        r"""
        Perform component initialization in network preprocessing.

        Parameters
        ----------
        nw : tespy.networks.network.Network
            Network this component is integrated in.
        """
        self.num_nw_fluids = len(nw.fluids)
        self.nw_fluids = nw.fluids
        self.always_all_equations = nw.always_all_equations
        self.num_nw_vars = self.num_nw_fluids + 3
        self.it = 0
        self.residual = []
        self.jacobian = None
        self.num_eq = num_eq
        self.vars = {}
        self.num_vars = 0

        for key, val in self.variables.items():
            data = self.get_attr(key)
            # component properties
            if isinstance(val, dc_cp):
                if data.is_var:
                    data.var_pos = self.num_vars
                    self.num_vars += 1
                    self.vars[data] = key
                if data.is_set and data.func is not None:
                    self.num_eq += 1

            # simple component properties
            elif isinstance(val, dc_simple):
                if data.is_set and data.func is not None:
                    self.num_eq += 1

            # component characteristics
            elif isinstance(val, dc_cc):
                if data.char_func is None:
                    try:
                        data.char_func = ldc(
                            self.component(), key, 'DEFAULT', CharLine)
                    except KeyError:
                        data.char_func = CharLine(x=[0, 1], y=[1, 1])

                if data.is_set and data.func is not None:
                    self.num_eq += 1

            # grouped component properties
            elif isinstance(val, dc_gcp):
                is_set = True
                for e in data.elements:
                    if not self.get_attr(e).is_set:
                        is_set = False

                if is_set:
                    data.set_attr(is_set=True)
                    self.num_eq += 1
                elif data.is_set:
                    start = (
                        'All parameters of the component group have to be '
                        'specified! This component group uses the following '
                        'parameters: ')
                    end = ' at ' + self.label + '. Group will be set to False.'
                    logging.warning(start + ', '.join(val.elements) + end)
                    val.set_attr(is_set=False)
                else:
                    val.set_attr(is_set=False)

        # set up Jacobian matrix and residual vector
        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))
        self.residual = np.zeros(self.num_eq)
        # done
        msg = (
            'The component ' + self.label + ' has ' + str(self.num_vars) +
            ' custom variables.')
        logging.debug(msg)

    @staticmethod
    def attr():
        return {}

    @staticmethod
    def inlets():
        return []

    @staticmethod
    def outlets():
        return []

    def get_char_expr(self, param, type='rel', inconn=0, outconn=0, doc=False):
        r"""
        Generic method to access characteristic function parameters.

        Parameters
        ----------
        param : str
            Parameter for characteristic function evaluation.

        type : str ('rel' or 'abs')
            Type of expression:

            - :code:`rel`: relative to design value
            - :code:`abs`: absolute value

        inconn : int
            Index of inlet connection.

        outconn : int
            Index of outlet connection.

        doc : boolean
            Return equation in LaTeX format instead of value.

        Returns
        -------
        expr : float, str
            Value of expression or LaTeX code for documentation if doc is True
        """
        if not doc:
            if type == 'rel':
                if param == 'm':
                    return (
                        self.inl[inconn].m.val_SI / self.inl[inconn].m.design)
                elif param == 'm_out':
                    return (
                        self.outl[outconn].m.val_SI /
                        self.outl[outconn].m.design)
                elif param == 'v':
                    v = self.inl[inconn].m.val_SI * v_mix_ph(
                        self.inl[inconn].to_flow(),
                        T0=self.inl[inconn].T.val_SI)
                    return v / self.inl[inconn].v.design
                elif param == 'pr':
                    return (
                        (self.outl[outconn].p.val_SI *
                         self.inl[inconn].p.design) /
                        (self.inl[inconn].p.val_SI *
                         self.outl[outconn].p.design))
                else:
                    msg = (
                        'The parameter ' + str(param) + ' is not available '
                        'for characteristic function evaluation.')
                    logging.error(msg)
                    raise ValueError(msg)
            else:
                if param == 'm':
                    return self.inl[inconn].m.val_SI
                elif param == 'm_out':
                    return self.outl[outconn].m.val_SI
                elif param == 'v':
                    return self.inl[inconn].m.val_SI * v_mix_ph(
                        self.inl[inconn].to_flow(),
                        T0=self.inl[inconn].T.val_SI)
                elif param == 'pr':
                    return (
                        self.outl[outconn].p.val_SI /
                        self.inl[inconn].p.val_SI)
                else:
                    return False
        else:
            if type == 'rel':
                if param == 'm':
                    return (
                        r'\frac{\dot{m}_\mathrm{in,' + str(inconn + 1) + r'}}'
                        r'{\dot{m}_\mathrm{in,' + str(inconn + 1) +
                        r',design}}')
                elif param == 'm_out':
                    return (
                        r'\frac{\dot{m}_\mathrm{out,' + str(outconn + 1) +
                        r'}}{\dot{m}_\mathrm{out,' + str(outconn + 1) +
                        r',design}}')
                elif param == 'v':
                    return (
                        r'\frac{\dot{V}_\mathrm{in,' + str(inconn + 1) + r'}}'
                        r'{\dot{V}_\mathrm{in,' + str(inconn + 1) +
                        r',design}}')
                elif param == 'pr':
                    return (
                        r'\frac{p_\mathrm{out,' + str(outconn + 1) +
                        r'}\cdot p_\mathrm{in,' + str(inconn + 1) +
                        r',design}}{p_\mathrm{out,' + str(outconn + 1) +
                        r',design}\cdot p_\mathrm{in,' + str(inconn + 1) +
                        r'}}')
            else:
                if param == 'm':
                    return r'\dot{m}_\mathrm{in,' + str(inconn + 1) + r'}'
                elif param == 'm_out':
                    return r'\dot{m}_\mathrm{out,' + str(outconn + 1) + r'}'
                elif param == 'v':
                    return r'\dot{V}_\mathrm{in,' + str(inconn + 1) + r'}'
                elif param == 'pr':
                    return (
                        r'\frac{p_\mathrm{out,' + str(outconn + 1) +
                        r'}}{p_\mathrm{in,' + str(inconn + 1) + r'}}')

    def solve(self, increment_filter, doc=False):
        """
        Solve equations and calculate partial derivatives of a component.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        doc : boolean
            Return equation in LaTeX format instead of value.
        """
        if doc:
            self.equation_docs = ['' for i in range(self.num_eq)]

        k = self.mandatory_equations(doc=doc)
        if not doc:
            self.mandatory_derivatives(increment_filter)

        for parameter, data in self.variables.items():
            if data.is_set and data.func is not None:
                residuals = data.func(doc=doc, **data.func_params)
                try:
                    num_eq = len(residuals)
                except TypeError:
                    num_eq = 1

                if not doc:
                    self.residual[k:k + num_eq] = residuals
                    data.deriv(increment_filter, k, **data.func_params)
                else:
                    self.equation_docs[k:k + num_eq] = residuals

                k += num_eq

    def mandatory_equations(self, doc=False):
        r"""
        Calculate residual vector of mandatory equations.

        Parameters
        ----------
        doc : boolean
            Return equation in LaTeX format instead of value.

        Returns
        -------
        k : int
            Position of last equation in residual value vector (k-th equation).
        """
        return 0

    def mandatory_derivatives(self, increment_filter):
        r"""
        Calculate partial derivatives for mandatory equations.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        Returns
        -------
        k : int
            Position of last equation in residual value vector (k-th equation).
        """
        # return value is required in some class of parent methods
        # e.g. desuperheater
        return 0

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

    def bus_deriv(self, bus):
        r"""
        Base method for partial derivatives of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        deriv : ndarray
            Matrix of partial derivatives.
        """
        return np.zeros((1, self.num_i + self.num_o, self.num_nw_vars))

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
        b = bus.comps.loc[self]
        comp_val = self.bus_func(b)
        if np.isnan(b['P_ref']) or b['P_ref'] == 0:
            expr = 1
        else:
            if b['base'] == 'component':
                expr = abs(comp_val / b['P_ref'])
            else:
                bus_value = newton(
                    bus_char_evaluation,
                    bus_char_derivative,
                    [comp_val, b['P_ref'], b['char']], 0,
                    val0=b['P_ref'], valmin=-1e15, valmax=1e15)
                expr = bus_value / b['P_ref']

        return b['char'].evaluate(expr)

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
        if np.isnan(b['P_ref']):
            expr = 1
        else:
            if b['base'] == 'component':
                expr = abs(comp_val / b['P_ref'])
            else:
                bus_value = newton(
                    bus_char_evaluation,
                    bus_char_derivative,
                    [comp_val, b['P_ref'], b['char']], 0,
                    val0=b['P_ref'], valmin=-1e15, valmax=1e15)
                return bus_value

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

    def propagate_fluid_to_target(self, inconn, start):
        r"""
        Propagate the fluids towards connection's target in recursion.

        Parameters
        ----------
        inconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        conn_idx = self.inl.index(inconn)
        outconn = self.outl[conn_idx]

        for fluid, x in inconn.fluid.val.items():
            if (outconn.fluid.val_set[fluid] is False and
                    outconn.good_starting_values is False):
                outconn.fluid.val[fluid] = x

        outconn.target.propagate_fluid_to_target(outconn, start)

    def propagate_fluid_to_source(self, outconn, start):
        r"""
        Propagate the fluids towards connection's source in recursion.

        Parameters
        ----------
        outconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        conn_idx = self.outl.index(outconn)
        inconn = self.inl[conn_idx]

        for fluid, x in outconn.fluid.val.items():
            if (inconn.fluid.val_set[fluid] is False and
                    inconn.good_starting_values is False):
                inconn.fluid.val[fluid] = x

        inconn.source.propagate_fluid_to_source(inconn, start)

    def set_parameters(self, mode, data):
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

        for key, dc in self.variables.items():
            if isinstance(dc, dc_cp):
                if ((mode == 'offdesign' and not self.local_design) or
                        (mode == 'design' and self.local_offdesign)):
                    self.get_attr(key).design = data[key]

                else:
                    self.get_attr(key).design = np.nan

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        return

    def check_parameter_bounds(self):
        for p in self.variables.keys():
            data = self.get_attr(p)
            if isinstance(data, dc_cp):
                if data.val > data.max_val + err:
                    msg = (
                        'Invalid value for ' + p + ': ' + p + ' = ' +
                        str(data.val) + ' above maximum value (' +
                        str(data.max_val) + ') at component ' + self.label +
                        '.')
                    logging.warning(msg)

                elif data.val < data.min_val - err:
                    msg = (
                        'Invalid value for ' + p + ': ' + p + ' = ' +
                        str(data.val) + ' below minimum value (' +
                        str(data.min_val) + ') at component ' + self.label +
                        '.')
                    logging.warning(msg)

            elif isinstance(data, dc_cc) and data.is_set:
                expr = self.get_char_expr(data.param, **data.char_params)
                data.char_func.get_bound_errors(expr, self.label)

            elif isinstance(data, dc_cm) and data.is_set:
                if isinstance(data.char_func, CompressorMap):
                    x = np.sqrt(self.inl[0].T.design / self.inl[0].T.val_SI)
                    y = (self.inl[0].m.val_SI * self.inl[0].p.design) / (
                        self.inl[0].m.design * self.inl[0].p.val_SI * x)
                    self.char_map.char_func.get_bound_errors(
                        x, y, self.igva.val, self.label)

    def initialise_fluids(self):
        return

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
        self.E_bus = np.nan
        self.E_D = np.nan
        self.epsilon = np.nan

    def get_plotting_data(self):
        return

    def generate_latex(self, eqn, func):
        latex = (
            r'\begin{equation}' + '\n' + r'\label{eq:' +
            self.__class__.__name__ + '_' + func + r'}' + '\n'
        )
        latex += eqn + '\n'
        latex += r'\end{equation}'
        return latex

    def fluid_func(self, doc=False):
        r"""
        Calculate the vector of residual values for fluid balance equations.

        Returns
        -------
        residual : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = x_{fl,in,i} - x_{fl,out,i} \; \forall fl \in
                \text{network fluids,} \; \forall i \in \text{inlets}
        """
        if not doc:
            residual = []
            for i in range(self.num_i):
                for fluid, x in self.inl[0].fluid.val.items():
                    residual += [x - self.outl[0].fluid.val[fluid]]
            return residual
        else:
            indices = list(range(1, self.num_i + 1))
            if len(indices) > 1:
                indices = ', '.join(str(idx) for idx in indices)
            else:
                indices = str(indices[0])
            latex = (
                r'0=x_{fl\mathrm{,in,}i}-x_{fl\mathrm{,out,}i}\;'
                r'\forall fl \in\text{network fluids,}'
                r'\; \forall i \in [' + indices + r']')
            return (
                [self.generate_latex(latex, 'fluid_func')] +
                (self.num_i - 1) * [''])

    def fluid_deriv(self):
        r"""
        Calculate partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((self.num_nw_fluids * self.num_i,
                          2 * self.num_i + self.num_vars,
                          self.num_nw_vars))
        for i in range(self.num_i):
            for j in range(self.num_nw_fluids):
                deriv[i * self.num_nw_fluids + j, i, j + 3] = 1
                deriv[i * self.num_nw_fluids + j, self.num_i + i, j + 3] = -1
        return deriv

    def mass_flow_func(self, doc=False):
        r"""
        Calculate the residual value for mass flow balance equation.

        Returns
        -------
        residual : list
            Vector with residual value for component's mass flow balance.

            .. math::

                0 = \dot{m}_{in,i} -\dot{m}_{out,i} \;\forall i\in\text{inlets}
        """
        if not doc:
            residual = []
            for i in range(self.num_i):
                residual += [self.inl[i].m.val_SI - self.outl[i].m.val_SI]
            return residual
        else:
            indices = list(range(1, self.num_i + 1))
            if len(indices) > 1:
                indices = ', '.join(str(idx) for idx in indices)
            else:
                indices = str(indices[0])
            latex = (
                r'0=\dot{m}_{\mathrm{in,}i}-\dot{m}_{\mathrm{out,}i}'
                r'\; \forall i \in [' + indices + r']')
            return (
                [self.generate_latex(latex, 'mass_flow_func')] +
                (self.num_i - 1) * [''])

    def mass_flow_deriv(self):
        r"""
        Calculate partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the mass flow balance
            equations.
        """
        deriv = np.zeros((
            self.num_i,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))
        for i in range(self.num_i):
            deriv[i, i, 0] = 1
        for j in range(self.num_o):
            deriv[j, j + i + 1, 0] = -1
        return deriv

    def numeric_deriv(self, func, dx, pos, **kwargs):
        r"""
        Calculate partial derivative of the function func to dx.

        Parameters
        ----------
        func : function
            Function :math:`f` to calculate the partial derivative for.

        dx : str
            Partial derivative.

        pos : int
            Position of connection regarding to inlets and outlet of the
            component, logic: ['in1', 'in2', ..., 'out1', ...] ->
            0, 1, ..., n, n + 1, ..., n + m

        Returns
        -------
        deriv : float/list
            Partial derivative(s) of the function :math:`f` to variable(s)
            :math:`x`.

            .. math::

                \frac{\partial f}{\partial x} = \frac{f(x + d) + f(x - d)}{2 d}
        """
        dm, dp, dh, df = 0, 0, 0, 0
        if dx == 'm':
            dm = 1e-4
        elif dx == 'p':
            dp = 1e-1
        elif dx == 'h':
            dh = 1e-1
        else:
            df = 1e-5

        if dx == 'fluid':
            deriv = []
            for f in self.inl[0].fluid.val.keys():
                val = (self.inl + self.outl)[pos].fluid.val[f]
                exp = 0
                if (self.inl + self.outl)[pos].fluid.val[f] + df <= 1:
                    (self.inl + self.outl)[pos].fluid.val[f] += df
                else:
                    (self.inl + self.outl)[pos].fluid.val[f] = 1
                exp += func(**kwargs)
                if (self.inl + self.outl)[pos].fluid.val[f] - 2 * df >= 0:
                    (self.inl + self.outl)[pos].fluid.val[f] -= 2 * df
                else:
                    (self.inl + self.outl)[pos].fluid.val[f] = 0
                exp -= func(**kwargs)
                (self.inl + self.outl)[pos].fluid.val[f] = val

                deriv += [exp / (2 * (dm + dp + dh + df))]

        elif dx in ['m', 'p', 'h']:
            exp = 0
            (self.inl + self.outl)[pos].m.val_SI += dm
            (self.inl + self.outl)[pos].p.val_SI += dp
            (self.inl + self.outl)[pos].h.val_SI += dh
            exp += func(**kwargs)

            (self.inl + self.outl)[pos].m.val_SI -= 2 * dm
            (self.inl + self.outl)[pos].p.val_SI -= 2 * dp
            (self.inl + self.outl)[pos].h.val_SI -= 2 * dh
            exp -= func(**kwargs)
            deriv = exp / (2 * (dm + dp + dh + df))

            (self.inl + self.outl)[pos].m.val_SI += dm
            (self.inl + self.outl)[pos].p.val_SI += dp
            (self.inl + self.outl)[pos].h.val_SI += dh

        else:
            d = self.get_attr(dx).d
            exp = 0
            self.get_attr(dx).val += d
            exp += func(**kwargs)

            self.get_attr(dx).val -= 2 * d
            exp -= func(**kwargs)
            deriv = exp / (2 * d)

            self.get_attr(dx).val += d

        return deriv

    def pr_func(self, pr='', inconn=0, outconn=0, doc=False):
        r"""
        Calculate residual value of pressure ratio function.

        Parameters
        ----------
        pr : str
            Component parameter to evaluate the pr_func on, e.g.
            :code:`pr1`.

        inconn : int
            Connection index of inlet.

        outconn : int
            Connection index of outlet.

        doc : boolean
            Return equation in LaTeX format instead of value.

        Returns
        -------
        residual : float
            Residual value of function.

            .. math::

                0 = p_{in} \cdot pr - p_{out}
        """
        if not doc:
            pr = self.get_attr(pr)
            return (self.inl[inconn].p.val_SI * pr.val -
                    self.outl[outconn].p.val_SI)
        else:
            latex = (
                r'0=p_\mathrm{in,' + str(inconn + 1) + r'}\cdot ' + pr +
                r' - p_\mathrm{out,' + str(outconn + 1) + r'}'
            )
            return [self.generate_latex(latex, 'pr_func_' + pr)]

    def pr_deriv(self, increment_filter, k, pr='', inconn=0, outconn=0):
        r"""
        Calculate residual value of pressure ratio function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.

        pr : str
            Component parameter to evaluate the pr_func on, e.g.
            :code:`pr1`.

        inconn : int
            Connection index of inlet.

        outconn : int
            Connection index of outlet.
        """
        pr = self.get_attr(pr)
        self.jacobian[k, inconn, 1] = pr.val
        self.jacobian[k, self.num_i + outconn, 1] = -1
        if pr.is_var:
            pos = self.num_i + self.num_o + pr.var_pos
            self.jacobian[k, pos, 0] = self.inl[inconn].p.val_SI

    def zeta_func(self, zeta='', inconn=0, outconn=0, doc=False):
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

        doc : boolean
            Return equation in LaTeX format instead of value.

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
        if not doc:
            data = self.get_attr(zeta)
            i = self.inl[inconn].to_flow()
            o = self.outl[outconn].to_flow()

            if abs(i[0]) < 1e-4:
                return i[1] - o[1]

            else:
                v_i = v_mix_ph(i, T0=self.inl[inconn].T.val_SI)
                v_o = v_mix_ph(o, T0=self.outl[outconn].T.val_SI)
                return (data.val - (i[1] - o[1]) * np.pi ** 2 /
                        (8 * abs(i[0]) * i[0] * (v_i + v_o) / 2))
        else:
            inl = r'_\mathrm{in,' + str(inconn + 1) + r'}'
            outl = r'_\mathrm{out,' + str(outconn + 1) + r'}'
            latex = (
                r'0 = \begin{cases}' + '\n' +
                r'p' + inl + r'- p' + outl + r' & |\dot{m}' + inl +
                r'| < \unitfrac[0.0001]{kg}{s} \\' + '\n' +
                r'\frac{\zeta}{D^4}-\frac{(p' + inl + r'-p' + outl + r')'
                r'\cdot\pi^2}{8\cdot\dot{m}' + inl + r'\cdot|\dot{m}' + inl +
                r'|\cdot\frac{v' + inl + r'v' + outl + r'}{2}}' +
                r'& |\dot{m}' + inl + r'| \geq \unitfrac[0.0001]{kg}{s}' + '\n'
                r'\end{cases}'
            )
            return [self.generate_latex(latex, 'zeta_func_' + zeta)]

    def zeta_deriv(self, increment_filter, k, zeta='', inconn=0, outconn=0):
        r"""
        Calculate partial derivatives of zeta function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.

        zeta : str
            Component parameter to evaluate the zeta_func on, e.g.
            :code:`zeta1`.

        inconn : int
            Connection index of inlet.

        outconn : int
            Connection index of outlet.
        """
        data = self.get_attr(zeta)
        f = self.zeta_func
        outpos = self.num_i + outconn
        if not increment_filter[inconn, 0]:
            self.jacobian[k, inconn, 0] = self.numeric_deriv(
                f, 'm', inconn, zeta=zeta, inconn=inconn, outconn=outconn)
        if not increment_filter[inconn, 2]:
            self.jacobian[k, inconn, 1] = self.numeric_deriv(
                f, 'p', inconn, zeta=zeta, inconn=inconn, outconn=outconn)
        if not increment_filter[inconn, 2]:
            self.jacobian[k, inconn, 2] = self.numeric_deriv(
                f, 'h', inconn, zeta=zeta, inconn=inconn, outconn=outconn)
        if not increment_filter[outpos, 1]:
            self.jacobian[k, outpos, 1] = self.numeric_deriv(
                f, 'p', outpos, zeta=zeta, inconn=inconn, outconn=outconn)
        if not increment_filter[outpos, 2]:
            self.jacobian[k, outpos, 2] = self.numeric_deriv(
                f, 'h', outpos, zeta=zeta, inconn=inconn, outconn=outconn)
        # custom variable zeta
        if data.is_var:
            pos = self.num_i + self.num_o + data.var_pos
            self.jacobian[k, pos, 0] = self.numeric_deriv(
                f, zeta, 2, zeta=zeta, inconn=inconn, outconn=outconn)
