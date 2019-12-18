# -*- coding: utf-8

"""Module class component.

All tespy components inherit from this class.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/components.py

SPDX-License-Identifier: MIT
"""

import numpy as np

import logging

from tespy.tools.characteristics import char_line, char_map, compressor_map
from tespy.tools.characteristics import load_default_char as ldc
from tespy.tools.data_containers import (
        data_container, dc_cc, dc_cm, dc_cp, dc_gcp, dc_simple
        )
from tespy.tools.fluid_properties import v_mix_ph

# %%


class component:
    r"""
    Class component is the base class of all TESPy components.

    Parameters
    ----------
    label : str
        The label of the component.

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
    The initialisation method (__init__), setter method (set_attr) and getter
    method (get_attr) are used for instances of class component and its
    children.

    Allowed keywords in kwargs are 'design_path', 'design' and 'offdesign'.
    Additional keywords depend on the type of component you want to create.

    Example
    -------
    Basic example for a setting up a tespy.components.components.component
    object. This example does not run a tespy calculation.

    >>> from tespy.components.components import component
    >>> comp = component('myComponent')
    >>> type(comp)
    <class 'tespy.components.components.component'>
    """

    def __init__(self, label, **kwargs):

        # check if components label is of type str and for prohibited chars
        if not isinstance(label, str):
            msg = 'Component label must be of type str!'
            logging.error(msg)
            raise ValueError(msg)

        elif len([x for x in [';', ',', '.'] if x in label]) > 0:
            msg = ('Can\'t use ' + str([';', ',', '.']) + ' in label (' +
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
        var = self.attr()

        for key in var.keys():
            self.__dict__.update({key: var[key]})

        self.set_attr(**kwargs)

    def set_attr(self, **kwargs):
        r"""
        Sets reset or unset attributes of a component for provided arguments.

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
        :func:`tespy.components.components.component.set_attr` method.
        """
        var = self.attr().keys()

        # set specified values
        for key in kwargs:
            if key in var:
                # data container specification
                if isinstance(kwargs[key], data_container):
                    if isinstance(kwargs[key], type(self.get_attr(key))):
                        self.__dict__.update({key: kwargs[key]})

                    else:
                        msg = ('The keyword ' + key + ' expects a '
                               'data_container of type ' +
                               str(type(self.get_attr(key))) +
                               ', a data_container of type ' +
                               str(type(kwargs[key])) + ' was supplied.')
                        logging.error(msg)
                        raise TypeError(msg)

                elif isinstance(self.get_attr(key), dc_cp):
                    # value specification for component properties
                    if (isinstance(kwargs[key], float) or
                            isinstance(kwargs[key], np.float64) or
                            isinstance(kwargs[key], np.int64) or
                            isinstance(kwargs[key], int)):
                        if np.isnan(kwargs[key]):
                            self.get_attr(key).set_attr(is_set=False)
                            self.get_attr(key).set_attr(is_var=False)

                        else:
                            self.get_attr(key).set_attr(val=kwargs[key])
                            self.get_attr(key).set_attr(is_set=True)
                            self.get_attr(key).set_attr(is_var=False)

                    elif kwargs[key] == 'var':
                        self.get_attr(key).set_attr(is_set=True)
                        self.get_attr(key).set_attr(is_var=True)

                    # invalid datatype for keyword
                    else:
                        msg = ('Bad datatype for keyword argument ' + key +
                               ' at ' + self.label + '.')
                        logging.error(msg)
                        raise TypeError(msg)

                elif (isinstance(self.get_attr(key), dc_cc) or
                      isinstance(self.get_attr(key), dc_cm)):
                    # value specification for characteristic lines
                    if isinstance(kwargs[key], char_line):
                        self.get_attr(key).func = kwargs[key]

                    # value specification for characteristic maps
                    elif (isinstance(kwargs[key], char_map) or
                          isinstance(kwargs[key], compressor_map)):
                        self.get_attr(key).func = kwargs[key]

                    # invalid datatype for keyword
                    else:
                        msg = ('Bad datatype for keyword argument ' + key +
                               ' at ' + self.label + '.')
                        logging.error(msg)
                        raise TypeError(msg)

                elif isinstance(self.get_attr(key), dc_gcp):
                    # value specification of grouped component parameter method
                    if isinstance(kwargs[key], str):
                        self.get_attr(key).method = kwargs[key]

                    # invalid datatype for keyword
                    else:
                        msg = ('Bad datatype for keyword argument ' + key +
                               ' at ' + self.label + '.')
                        logging.error(msg)
                        raise TypeError(msg)

                elif isinstance(self.get_attr(key), dc_simple):
                    if (isinstance(kwargs[key], float) or
                            isinstance(kwargs[key], np.float64) or
                            isinstance(kwargs[key], np.int64) or
                            isinstance(kwargs[key], int)):
                        if np.isnan(kwargs[key]):
                            self.get_attr(key).set_attr(val_set=False)

                        else:
                            self.get_attr(key).set_attr(
                                    val=kwargs[key], val_set=True)

                    else:
                        self.get_attr(key).set_attr(
                                val=kwargs[key], val_set=True)

            elif key == 'design' or key == 'offdesign':
                if not isinstance(kwargs[key], list):
                    msg = ('Please provide the ' + key + ' parameters as list '
                           'at ' + self.label + '.')
                    logging.error(msg)
                    raise TypeError(msg)
                if set(kwargs[key]).issubset(list(var)):
                    self.__dict__.update({key: kwargs[key]})

                else:
                    msg = ('Available parameters for (off-)design '
                           'specification are: ' + str(list(var)) + ' at '
                           + self.label + '.')
                    logging.error(msg)
                    raise ValueError(msg)

            elif key == 'local_design' or key == 'local_offdesign':
                if not isinstance(kwargs[key], bool):
                    msg = ('Please provide the ' + key + ' as boolean '
                           'at ' + self.label + '.')
                    logging.error(msg)
                    raise TypeError(msg)

                else:
                    self.__dict__.update({key: kwargs[key]})

            elif key == 'design_path':
                if isinstance(kwargs[key], str):
                    self.__dict__.update({key: kwargs[key]})
                    self.new_design = True

                elif np.isnan(kwargs[key]):
                    self.design_path = None
                    self.new_design = True

                else:
                    msg = ('Please provide the ' + key + ' parameter as '
                           'string or as nan.')
                    logging.error(msg)
                    raise TypeError(msg)

            elif key == 'printout':
                if not isinstance(kwargs[key], bool):
                    msg = ('Please provide the ' + key + ' as boolean.')
                    logging.error(msg)
                    raise TypeError(msg)
                else:
                    self.__dict__.update({key: kwargs[key]})

            # invalid keyword
            else:
                msg = ('Component ' + self.label + ' has no attribute ' +
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

    def comp_init(self, nw):
        r"""
        Performs component initialization in network preprocessing.

        Parameters
        ----------
        nw : tespy.networks.network
            Network this component is integrated in.
        """
        self.vars = {}
        self.num_vars = 0
        var = self.attr()
        for val in var.keys():
            if isinstance(self.attr()[val], dc_cp):
                if self.get_attr(val).is_var:
                    self.get_attr(val).var_pos = self.num_vars
                    self.num_vars += 1
                    self.vars[self.get_attr(val)] = val

        msg = ('The component ' + self.label + ' has ' + str(self.num_vars) +
               ' custom variables.')
        logging.debug(msg)

        # characteristics creation
        for key, val in var.items():
            if isinstance(val, dc_cc):
                if self.get_attr(key).func is None:
                    try:
                        self.get_attr(key).func = ldc(
                            self.component(), key, 'DEFAULT', char_line)
                    except KeyError:
                        self.get_attr(key).func = char_line(x=[0, 1], y=[1, 1])

                    if self.char_warnings is True:
                        msg = ('Created characteristic line for parameter ' +
                               key + ' at component ' + self.label + ' from '
                               'default data.\n'
                               'You can specify your own data using '
                               'component.' + key +
                               '.set_attr(func=custom_char).\n'
                               'If you want to disable these warnings use '
                               'component.char_warnings=False.')
                        logging.warning(msg)

        self.num_fl = len(nw.fluids)
        self.fluids = nw.fluids

    def attr(self):
        return {}

    def inlets(self):
        return []

    def outlets(self):
        return []

    def equations(self):
        return []

    def derivatives(self):
        return []

    def initialise_source(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's
        outlet.

        Parameters
        ----------
        c : tespy.connections.connection
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
        Returns a starting value for pressure and enthalpy at component's
        inlet.

        Parameters
        ----------
        c : tespy.connections.connection
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
        if mode == 'design' or self.local_design is True:
            self.new_design = True

        for key, dc in self.attr().items():
            if isinstance(dc, dc_cp):
                if ((mode == 'offdesign' and self.local_design is False) or
                        (mode == 'design' and self.local_offdesign is True)):
                    self.get_attr(key).design = data[key]

                else:
                    self.get_attr(key).design = np.nan

    def calc_parameters(self):
        r"""
        Postprocessing parameter calculation.
        """
        return

    def check_parameter_bounds(self):
        for p, data in self.attr().items():
            if isinstance(data, dc_cp):
                val = self.get_attr(p).val
                if val > data.max_val:
                    msg = ('Invalid value for ' + p + ': ' + p + ' = ' +
                           str(val) + ' above maximum value (' +
                           str(data.max_val) + ') at component ' +
                           self.label + '.')
                    logging.warning(msg)

                elif val < data.min_val:
                    msg = ('Invalid value for ' + p + ': ' + p + ' = ' +
                           str(val) + ' below minimum value (' +
                           str(data.min_val) + ') at component ' +
                           self.label + '.')
                    logging.warning(msg)

    def initialise_fluids(self, nw):
        return

    def convergence_check(self, nw):
        return

# %%

    def fluid_func(self):
        r"""
        Calculates the vector of residual values for component's fluid balance
        equations.

        Returns
        -------
        vec_res : list
            Vector of residual values for component's fluid balance.

            .. math::
                0 = fluid_{i,in} - fluid_{i,out} \;
                \forall i \in \mathrm{fluid}
        """
        vec_res = []

        for fluid, x in self.inl[0].fluid.val.items():
            vec_res += [x - self.outl[0].fluid.val[fluid]]
        return vec_res

    def fluid_deriv(self):
        r"""
        Calculates the partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """

        deriv = np.zeros((self.num_fl, 2 + self.num_vars, 3 + self.num_fl))
        i = 0
        for fluid in self.fluids:
            deriv[i, 0, i + 3] = 1
            deriv[i, 1, i + 3] = -1
            i += 1
        return deriv.tolist()

# %%

    def mass_flow_func(self):
        r"""
        Calculates the residual value for component's mass flow balance
        equation.

        Returns
        -------
        vec_res : list
            Vector with residual value for component's mass flow balance.

            .. math::
                0 = \sum \dot{m}_{in,i} - \sum \dot{m}_{out,j} \;
                \forall i \in inlets, \forall j \in outlets
        """

        res = 0
        for i in self.inl:
            res += i.m.val_SI
        for o in self.outl:
            res -= o.m.val_SI
        return [res]

    def mass_flow_deriv(self):
        r"""
        Calculates the partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the mass flow balance
            equations.
        """

        deriv = np.zeros((1, self.num_i + self.num_o +
                          self.num_vars, 3 + self.num_fl))
        for i in range(self.num_i):
            deriv[0, i, 0] = 1
        for j in range(self.num_o):
            deriv[0, j + i + 1, 0] = -1
        return deriv.tolist()

# %%

    def numeric_deriv(self, func, dx, pos, **kwargs):
        r"""
        Calculates partial derivative of the function func to dx at given
        connection.

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
            dp = 1
        elif dx == 'h':
            dh = 1
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

# %%

    def zeta_func(self):
        r"""
        Calculates residual value of :math:`\zeta`-function.

        Returns
        -------
        val : float
            Residual value of function.

            .. math::

                val = \begin{cases}
                p_{in} - p_{out} & |\dot{m}| < \epsilon \\
                \frac{\zeta}{D^4} - \frac{(p_{in} - p_{out}) \cdot \pi^2}{8 \cdot
                \dot{m}_{in} \cdot |\dot{m}_{in}| \cdot \frac{v_{in} +
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
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()
        if hasattr(self, 'zeta'):
            val = self.zeta.val
        else:
            val = self.zeta1.val

        if abs(i[0]) < 1e-4:
            return i[1] - o[1]

        else:
            v_i = v_mix_ph(i, T0=self.inl[0].T.val_SI)
            v_o = v_mix_ph(o, T0=self.outl[0].T.val_SI)
            return (val - (i[1] - o[1]) * np.pi ** 2 /
                    (8 * abs(i[0]) * i[0] * (v_i + v_o) / 2))

    def zeta2_func(self):
        r"""
        calculates residual value of :math:`\zeta`-function (for heat
        exchangers at lower temperature side).

        Returns
        -------
        val : float
            Residual value of function.

            .. math::

                val = \begin{cases}
                p_{in} - p_{out} & |\dot{m}| < \epsilon \\
                \frac{\zeta_2}{D^4} - \frac{(p_{2,in} - p_{2,out}) \cdot \pi^2}
                {8 \cdot \dot{m}_{2,in} \cdot |\dot{m}_{2,in}| \cdot
                \frac{v_{2,in} + v_{2,out}}{2}} &
                |\dot{m}| > \epsilon
                \end{cases}

        Note
        ----
        The zeta value is caluclated on the basis of a given pressure loss at
        a given flow rate in the design case. As the cross sectional area A
        will not change, it is possible to handle the equation in this way:

        .. math::

            \frac{\zeta_2}{D^4} =  \frac{\Delta p_2 \cdot \pi^2}
            {8 \cdot \dot{m}_2^2 \cdot v}
        """
        i = self.inl[1].to_flow()
        o = self.outl[1].to_flow()

        if abs(i[0]) < 1e-4:
            return i[1] - o[1]
        else:
            v_i = v_mix_ph(i, T0=self.inl[1].T.val_SI)
            v_o = v_mix_ph(o, T0=self.outl[1].T.val_SI)
            return (self.zeta2.val - (i[1] - o[1]) * np.pi ** 2 /
                    (8 * abs(i[0]) * i[0] * (v_i + v_o) / 2))
