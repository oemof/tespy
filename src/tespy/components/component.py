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

from tespy.tools import logger
from tespy.tools.characteristics import CharLine
from tespy.tools.characteristics import CharMap
from tespy.tools.characteristics import load_default_char as ldc
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentCharacteristics as dc_gcc
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import v_mix_ph
from tespy.tools.global_vars import ERR
from tespy.tools.helpers import _numeric_deriv
from tespy.tools.helpers import bus_char_derivative
from tespy.tools.helpers import bus_char_evaluation
from tespy.tools.helpers import newton_with_kwargs


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
                f"You cannot use any of " + ", ".join(_forbidden) + " in a "
                f"component label ({self.component()}"
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
                        data.set_attr(val=kwargs[key], is_set=True)
                        if isinstance(data, dc_cp):
                            data.set_attr(is_var=False)

                    elif kwargs[key] == 'var' and isinstance(data, dc_cp):
                        data.set_attr(is_set=True, is_var=True)

                    elif isinstance(data, dc_simple):
                        data.set_attr(val=kwargs[key], is_set=True)

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

    @staticmethod
    def is_branch_source():
        return False

    def propagate_to_target(self, branch):
        inconn = branch["connections"][-1]
        conn_idx = self.inl.index(inconn)
        outconn = self.outl[conn_idx]

        branch["connections"] += [outconn]
        branch["components"] += [outconn.target]

        outconn.target.propagate_to_target(branch)

    def propagate_wrapper_to_target(self, branch):
        inconn = branch["connections"][-1]
        conn_idx = self.inl.index(inconn)
        outconn = self.outl[conn_idx]

        branch["connections"] += [outconn]
        branch["components"] += [self]

        outconn.target.propagate_wrapper_to_target(branch)

    def preprocess(self, num_nw_vars):
        r"""
        Perform component initialization in network preprocessing.

        Parameters
        ----------
        nw : tespy.networks.network.Network
            Network this component is integrated in.
        """
        self.it = 0
        self.num_eq = 0
        self.vars = {}
        self.num_vars = 0

        if self.bypass:
            self.constraints = self.get_bypass_constraints().copy()
        else:
            self.constraints = self.get_mandatory_constraints().copy()

        self.prop_specifications = {}
        self.var_specifications = {}
        self.group_specifications = {}
        self.char_specifications = {}
        self.__dict__.update(self.constraints)

        for constraint in self.constraints.values():
            self.num_eq += constraint['num_eq']

        if not self.bypass:
            self._setup_user_imposed_constraints(num_nw_vars)

        self.jacobian = {}
        self.residual = np.zeros(self.num_eq)

        sum_eq = 0
        for constraint in self.constraints.values():
            num_eq = constraint['num_eq']
            if constraint['constant_deriv']:
                constraint["deriv"](None, sum_eq)
            sum_eq += num_eq

        # done
        msg = f"The component {self.label} has {self.num_vars} variables."
        logger.debug(msg)

    def _setup_user_imposed_constraints(self, num_nw_vars):
        for key, val in self.parameters.items():
            data = self.get_attr(key)
            if isinstance(val, dc_cp):
                if data.is_var:
                    data.J_col = num_nw_vars + self.num_vars
                    self.num_vars += 1
                    self.vars[data] = key

                self.prop_specifications[key] = val.is_set
                self.var_specifications[key] = val.is_var

            # component characteristics
            elif isinstance(val, dc_cc):
                if data.func is not None:
                    self.char_specifications[key] = val.is_set
                if data.char_func is None:
                    try:
                        data.char_func = ldc(
                            self.component(), key, 'DEFAULT', CharLine)
                    except KeyError:
                        data.char_func = CharLine(x=[0, 1], y=[1, 1])

            # component characteristics
            elif isinstance(val, dc_cm):
                if data.func is not None:
                    self.char_specifications[key] = val.is_set
                if data.char_func is None:
                    try:
                        data.char_func = ldc(
                            self.component(), key, 'DEFAULT', CharMap)
                    except KeyError:
                        data.char_func = CharLine(x=[0, 1], y=[1, 1])

            # grouped component properties
            elif isinstance(val, dc_gcp):
                is_set = True
                for e in data.elements:
                    if not self.get_attr(e).is_set:
                        is_set = False

                if is_set:
                    data.set_attr(is_set=True)
                elif data.is_set:
                    start = (
                        'All parameters of the component group have to be '
                        'specified! This component group uses the following '
                        'parameters: '
                    )
                    end = f" at {self.label}. Group will be set to False."
                    logger.warning(start + ', '.join(val.elements) + end)
                    val.set_attr(is_set=False)
                else:
                    val.set_attr(is_set=False)
                self.group_specifications[key] = val.is_set

            # grouped component characteristics
            elif isinstance(val, dc_gcc):
                self.group_specifications[key] = val.is_set

            # component properties
            if data.is_set and data.func is not None:
                self.num_eq += data.num_eq

    def get_parameters(self):
        return {}

    def get_mandatory_constraints(self):
        return {}

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
    def is_variable(var, increment_filter=None):
        if var.is_var:
            if increment_filter is None or not increment_filter[var.J_col]:
                return True
        return False

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
                v = self.inl[inconn].m.val_SI * v_mix_ph(
                    self.inl[inconn].p.val_SI, self.inl[inconn].h.val_SI,
                    self.inl[inconn].fluid_data, self.inl[inconn].mixing_rule,
                    T0=self.inl[inconn].T.val_SI
                )
                return v / self.inl[inconn].v.design
            elif param == 'pr':
                return (
                    (self.outl[outconn].p.val_SI * self.inl[inconn].p.design)
                    / (self.inl[inconn].p.val_SI * self.outl[outconn].p.design)
                )
            else:
                msg = (
                    f"The parameter {param}) is not available for "
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
                return self.inl[inconn].m.val_SI * v_mix_ph(
                    self.inl[inconn].p.val_SI, self.inl[inconn].h.val_SI,
                    self.inl[inconn].fluid_data, self.inl[inconn].mixing_rule,
                    T0=self.inl[inconn].T.val_SI
                )
            elif param == 'pr':
                return self.outl[outconn].p.val_SI / self.inl[inconn].p.val_SI
            else:
                return False

    def get_char_expr_doc(self, param, type='rel', inconn=0, outconn=0):
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
        expr : str
            LaTeX code for documentation
        """
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

    def solve(self, increment_filter):
        """
        Solve equations and calculate partial derivatives of a component.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.
        """
        sum_eq = 0
        for constraint in self.constraints.values():
            num_eq = constraint['num_eq']
            if num_eq > 0:
                self.residual[sum_eq:sum_eq + num_eq] = constraint['func']()
            if not constraint['constant_deriv']:
                constraint['deriv'](increment_filter, sum_eq)
            sum_eq += num_eq

        if self.bypass:
            return

        for data in self.parameters.values():
            if data.is_set and data.func is not None:
                self.residual[sum_eq:sum_eq + data.num_eq] = data.func(
                    **data.func_params
                )
                data.deriv(increment_filter, sum_eq, **data.func_params)

                sum_eq += data.num_eq

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

        for key, dc in self.parameters.items():
            if isinstance(dc, dc_cp):
                if (
                        ((mode == 'offdesign' and not self.local_design) or
                        (mode == 'design' and self.local_offdesign)) and
                        (data[key] is not None)
                    ):
                    self.get_attr(key).design = float(data[key])

                else:
                    self.get_attr(key).design = np.nan

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        return

    def check_parameter_bounds(self):
        r"""Check parameter value limits."""
        for p in self.parameters.keys():
            data = self.get_attr(p)
            if isinstance(data, dc_cp):
                if data.val > data.max_val + ERR:
                    msg = (
                        f"Invalid value for {p}: {p} = {data.val} above "
                        f"maximum value ({data.max_val}) at component "
                        f"{self.label}."
                    )
                    logger.warning(msg)

                elif data.val < data.min_val - ERR:
                    msg = (
                        f"Invalid value for {p}: {p} = {data.val} below "
                        f"minimum value ({data.min_val}) at component "
                        f"{self.label}."
                    )
                    logger.warning(msg)

            elif isinstance(data, dc_cc) and data.is_set:
                expr = self.get_char_expr(data.param, **data.char_params)
                data.char_func.get_domain_errors(expr, self.label)

            elif isinstance(data, dc_gcc) and data.is_set:
                for char in data.elements:
                    char_data = self.get_attr(char)
                    expr = self.get_char_expr(
                        char_data.param, **char_data.char_params)
                    char_data.char_func.get_domain_errors(expr, self.label)

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

    def pressure_equality_func(self):
        r"""
        Equation for pressure equality.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = p_{in,i} - p_{out,i} \;\forall i\in\text{inlets}
        """
        residual = []
        for i in range(self.num_i):
            residual += [self.inl[i].p.val_SI - self.outl[i].p.val_SI]
        return residual

    def pressure_equality_func_doc(self, label):
        r"""
        Equation for pressure equality.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        indices = list(range(1, self.num_i + 1))
        if len(indices) > 1:
            indices = ', '.join(str(idx) for idx in indices)
        else:
            indices = str(indices[0])
        latex = (
            r'0=p_{\mathrm{in,}i}-p_{\mathrm{out,}i}'
            r'\; \forall i \in [' + indices + r']')
        return generate_latex_eq(self, latex, label)

    def pressure_equality_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the mass flow balance
            equations.
        """
        for i in range(self.num_i):
            if self.inl[i].p.is_var:
                self.jacobian[k + i, self.inl[i].p.J_col] = 1
            if self.outl[i].p.is_var:
                self.jacobian[k + i, self.outl[i].p.J_col] = -1

    def enthalpy_equality_func(self):
        r"""
        Equation for enthalpy equality.

        Returns
        -------
        residual : list
            Residual values of equations.

            .. math::

                0 = h_{in,i} - h_{out,i} \;\forall i\in\text{inlets}
        """
        residual = []
        for i in range(self.num_i):
            residual += [self.inl[i].h.val_SI - self.outl[i].h.val_SI]
        return residual

    def enthalpy_equality_func_doc(self, label):
        r"""
        Equation for enthalpy equality.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        indices = list(range(1, self.num_i + 1))
        if len(indices) > 1:
            indices = ', '.join(str(idx) for idx in indices)
        else:
            indices = str(indices[0])
        latex = (
            r'0=h_{\mathrm{in,}i}-h_{\mathrm{out,}i}'
            r'\; \forall i \in [' + indices + r']'
        )
        return generate_latex_eq(self, latex, label)

    def enthalpy_equality_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the mass flow balance
            equations.
        """
        for i in range(self.num_i):
            if self.inl[i].h.is_var:
                self.jacobian[k + i, self.inl[i].h.J_col] = 1
            if self.outl[i].h.is_var:
                self.jacobian[k + i, self.outl[i].h.J_col] = -1

    def numeric_deriv(self, func, dx, conn=None, **kwargs):
        r"""
        Calculate partial derivative of the function func to dx.

        For details see :py:func:`tespy.tools.helpers._numeric_deriv`
        """
        return _numeric_deriv(self, func, dx, conn, **kwargs)

    def pr_func(self, pr=None, inconn=0, outconn=0):
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

        Returns
        -------
        residual : float
            Residual value of function.

            .. math::

                0 = p_{in} \cdot pr - p_{out}
        """
        pr = self.get_attr(pr)
        return self.inl[inconn].p.val_SI * pr.val - self.outl[outconn].p.val_SI

    def pr_func_doc(self, label, pr='', inconn=0, outconn=0):
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

        Returns
        -------
        residual : float
            Residual value of function.
        """
        latex = (
            r'0=p_\mathrm{in,' + str(inconn + 1) + r'}\cdot ' + pr +
            r' - p_\mathrm{out,' + str(outconn + 1) + r'}'
        )
        return generate_latex_eq(self, latex, label)

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
        i = self.inl[inconn]
        o = self.outl[inconn]
        if i.p.is_var:
            self.jacobian[k, i.p.J_col] = pr.val
        if o.p.is_var:
            self.jacobian[k, o.p.J_col] = -1
        if pr.is_var:
            self.jacobian[k, self.pr.J_col] = i.p.val_SI

    def calc_zeta(self, i, o):
        if abs(i.m.val_SI) <= 1e-4:
            return 0
        else:
            return (
                (i.p.val_SI - o.p.val_SI) * math.pi ** 2
                / (4 * i.m.val_SI ** 2 * (i.vol.val_SI + o.vol.val_SI))
            )

    def zeta_func(self, zeta='', inconn=0, outconn=0):
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
            v_i = v_mix_ph(i.p.val_SI, i.h.val_SI, i.fluid_data, i.mixing_rule, T0=i.T.val_SI)
            v_o = v_mix_ph(o.p.val_SI, o.h.val_SI, o.fluid_data, o.mixing_rule, T0=o.T.val_SI)
            return (
                data.val - (i.p.val_SI - o.p.val_SI) * math.pi ** 2
                / (8 * abs(i.m.val_SI) * i.m.val_SI * (v_i + v_o) / 2)
            )

    def zeta_func_doc(self, label, zeta='', inconn=0, outconn=0):
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
        """
        inl = r'_\mathrm{in,' + str(inconn + 1) + r'}'
        outl = r'_\mathrm{out,' + str(outconn + 1) + r'}'
        latex = (
            r'0 = \begin{cases}' + '\n' +
            r'p' + inl + r'- p' + outl + r' & |\dot{m}' + inl +
            r'| < \unitfrac[0.0001]{kg}{s} \\' + '\n' +
            r'\frac{\zeta}{D^4}-\frac{(p' + inl + r'-p' + outl + r')'
            r'\cdot\pi^2}{8\cdot\dot{m}' + inl + r'\cdot|\dot{m}' + inl +
            r'|\cdot\frac{v' + inl + r' + v' + outl + r'}{2}}' +
            r'& |\dot{m}' + inl + r'| \geq \unitfrac[0.0001]{kg}{s}' + '\n'
            r'\end{cases}'
        )
        return generate_latex_eq(self, latex, label)

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
        i = self.inl[inconn]
        o = self.outl[outconn]
        kwargs = dict(zeta=zeta, inconn=inconn, outconn=outconn)
        if self.is_variable(i.m, increment_filter):
            self.jacobian[k, i.m.J_col] = self.numeric_deriv(f, 'm', i, **kwargs)
        if self.is_variable(i.p, increment_filter):
            self.jacobian[k, i.p.J_col] = self.numeric_deriv(f, 'p', i, **kwargs)
        if self.is_variable(i.h, increment_filter):
            self.jacobian[k, i.h.J_col] = self.numeric_deriv(f, 'h', i, **kwargs)
        if self.is_variable(o.p, increment_filter):
            self.jacobian[k, o.p.J_col] = self.numeric_deriv(f, 'p', o, **kwargs)
        if self.is_variable(o.h, increment_filter):
            self.jacobian[k, o.h.J_col] = self.numeric_deriv(f, 'h', o, **kwargs)
        # custom variable zeta
        if data.is_var:
            self.jacobian[k, data.J_col] = self.numeric_deriv(f, zeta, None, **kwargs)

    def dp_func(self, dp=None, inconn=None, outconn=None):
        """Calculate residual value of pressure difference function.

        Parameters
        ----------
        dp : str
            Component parameter to evaluate the dp_func on, e.g.
            :code:`dp1`.

        inconn : int
            Connection index of inlet.

        outconn : int
            Connection index of outlet.

        Returns
        -------
        residual : float
            Residual value of function.

            .. math::

                0 = p_{in} - p_{out} - dp
        """
        inlet_conn = self.inl[inconn]
        outlet_conn = self.outl[outconn]
        dp_value = self.get_attr(dp).val_SI
        return inlet_conn.p.val_SI - outlet_conn.p.val_SI - dp_value

    def dp_deriv(self, increment_filter, k, dp=None, inconn=None, outconn=None):
        r"""
        Calculate residual value of pressure difference function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.

        dp : str
            Component parameter to evaluate the dp_func on, e.g.
            :code:`dp1`.

        inconn : int
            Connection index of inlet.

        outconn : int
            Connection index of outlet.
        """
        inlet_conn = self.inl[inconn]
        outlet_conn = self.outl[outconn]
        if inlet_conn.p.is_var:
            self.jacobian[k, inlet_conn.p.J_col] = 1
        if outlet_conn.p.is_var:
            self.jacobian[k, outlet_conn.p.J_col] = -1
