# -*- coding: utf-8

"""
.. module:: components
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import numpy as np
import math

import logging

import CoolProp.CoolProp as CP

from tespy.tools.characteristics import char_map, characteristics

from tespy.tools.data_containers import (
        data_container, dc_cc, dc_cm, dc_cp, dc_gcp, dc_simple
        )
from tespy.tools.fluid_properties import (
        h_mix_pQ, h_mix_pT, dh_mix_dpQ,
        memorise,
        s_mix_ph,
        T_bp_p,
        T_mix_ph, dT_mix_dph, dT_mix_pdh, dT_mix_ph_dfluid,
        v_mix_ph,
        visc_mix_ph
        )
from tespy.tools.global_vars import molar_masses

from tespy.tools.helpers import (
        lamb, num_fluids, single_fluid, TESPyComponentError
        )

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

    >>> from tespy import cmp
    >>> comp = cmp.component('myComponent')
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
        self.interface = False

        self.new_design = True
        self.design_path = None
        self.design = []
        self.offdesign = []
        self.local_design = False
        self.local_offdesign = False

        # add container for components attributes
        var = self.attr()

        for key in var.keys():
            self.__dict__.update({key: var[key]})

        self.set_attr(**kwargs)

    def set_attr(self, **kwargs):
        r"""
        Sets, resets or unsets attributes of a component for provided keyword
        arguments.

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
                      isinstance(self.get_attr(key), dc_cm) or
                      isinstance(self.get_attr(key), dc_gcp)):
                    # value specification for component characteristics
                    if isinstance(kwargs[key], str):
                        self.get_attr(key).set_attr(method=kwargs[key])
                        if (isinstance(self.get_attr(key), dc_cc) or
                                isinstance(self.get_attr(key), dc_cm)):
                            self.get_attr(key).set_attr(func=None)

                    elif (isinstance(kwargs[key], characteristics) or
                          isinstance(kwargs[key], char_map)):
                        self.get_attr(key).func=kwargs[key]
                        self.get_attr(key).x = self.get_attr(key).func.x
                        self.get_attr(key).y = self.get_attr(key).func.y

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

            # export sources or sinks as subsystem interface
            elif key == 'interface':
                if isinstance(self, source) or isinstance(self, sink):
                    if isinstance(kwargs[key], bool):
                        self.interface = kwargs[key]
                    else:
                        msg = ('Datatype for keyword argument ' + str(key) +
                               ' must be bool at ' + self.label + '.')
                        logging.error(msg)
                        raise ValueError(msg)
                else:
                    msg = ('Only sinks and sources can be attributed with the '
                           'interface parameter at ' + self.label + ').')
                    logging.error(msg)
                    raise TESPyComponentError(msg)

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
        for var in self.attr().keys():
            if isinstance(self.attr()[var], dc_cp):
                if self.get_attr(var).is_var:
                    self.get_attr(var).var_pos = self.num_vars
                    self.num_vars += 1
                    self.vars[self.get_attr(var)] = var

        msg = ('The component ' + self.label + ' has ' + str(self.num_vars) +
               ' custom variables.')
        logging.debug(msg)

        # characteristics creation
        for key, val in self.attr().items():
            if isinstance(val, dc_cc):
                generate_char = False
                if self.get_attr(key).func is None:
                    generate_char = True
                elif (not np.array_equal(self.get_attr(key).func.x,
                                         self.get_attr(key).x) or
                      not np.array_equal(self.get_attr(key).func.y,
                                         self.get_attr(key).y)):
                    generate_char = True

                if generate_char:
                    self.get_attr(key).func = characteristics(
                            method=self.get_attr(key).method,
                            x=self.get_attr(key).x,
                            y=self.get_attr(key).y, comp=self.component())
                    self.get_attr(key).x = self.get_attr(key).func.x
                    self.get_attr(key).y = self.get_attr(key).func.y

                    msg = ('Generated characteristic line for attribute ' +
                           key + ' at component ' + self.label + '.')
                    logging.debug(msg)

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
            return (val - (i[1] - o[1]) * math.pi ** 2 /
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
            return (self.zeta2.val - (i[1] - o[1]) * math.pi ** 2 /
                    (8 * abs(i[0]) * i[0] * (v_i + v_o) / 2))

# %%


class source(component):
    r"""
    A flow originates from a source.

    Equations
        This component is unconstrained.

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Example
    -------
    >>> from tespy import cmp
    >>> so = cmp.source('some source')
    >>> so.component()
    'source'
    >>> so.label
    'some source'
    """

    def component(self):
        return 'source'

    def outlets(self):
        return ['out1']

# %%


class sink(component):
    r"""
    A flow drains in a sink.

    Equations
        This component is unconstrained.

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Example
    -------
    >>> from tespy import cmp
    >>> si = cmp.sink('some sink')
    >>> si.component()
    'sink'
    >>> si.label
    'some sink'
    """

    def component(self):
        return 'sink'

    def inlets(self):
        return ['in1']

# %%


class node(component):
    r"""
    The component node is the parent class for splitter, separator and merge.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.mass_flow_func`

        .. math::

            0 = p_{in,1} - p_{out,i}\\
            \forall i \in \mathrm{outlets}

        **additional equations**

        - :func:`tespy.components.components.splitter.additional_equations`
        - :func:`tespy.components.components.separator.additional_equations`
        - :func:`tespy.components.components.merge.additional_equations`

    Inlets/Outlets

        - specify number of outlets with :code:`num_in` (default value: 2)
        - specify number of outlets with :code:`num_out` (default value: 2)

    Image

        .. image:: _images/node.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    num_in : float/tespy.helpers.dc_simple
        Number of inlets for this component, default value: 2.

    num_out : float/tespy.helpers.dc_simple
        Number of outlets for this component, default value: 2.

    Note
    ----
    - Node: Fluid composition and enthalpy at all **outgoing** connections
      (mass flow leaves the node) is result of mixture of the properties of
      the incoming connections (mass flow enters node).
      Incoming and outgoing connections can be a result of the calculation and
      are not identical to the inlets and outlets!
    - Splitter: Fluid composition and enthalpy at all outlets is the same as
      the inlet's properties.
    - Separator: Fluid composition is variable for all outlets, temperature at
      all outlets is the same as the inlet's temperature.
    - Merge: Fluid composition and enthalpy at outlet is result of mixture of
      the inlet's properties.

    Example
    -------
    >>> from tespy import cmp, con, nwk, hlp
    >>> import numpy as np
    >>> fluid_list = ['O2', 'N2']
    >>> nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ...     h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> so1 = cmp.source('source1')
    >>> so2 = cmp.source('source2')
    >>> si1 = cmp.sink('sink1')
    >>> si2 = cmp.sink('sink2')
    >>> n = cmp.node('node', num_in=2, num_out=2)
    >>> n.component()
    'node'
    >>> inc1 = con.connection(so1, 'out1', n, 'in1')
    >>> inc2 = con.connection(so2, 'out1', n, 'in2')
    >>> outg1 = con.connection(n, 'out1', si1, 'in1')
    >>> outg2 = con.connection(n, 'out2', si2, 'in1')
    >>> nw.add_conns(inc1, inc2, outg1, outg2)
    >>> inc1.set_attr(fluid={'O2': 1, 'N2': 0}, p=1, T=20, m=2)
    >>> inc2.set_attr(fluid={'O2': 0.5, 'N2': 0.5}, T=50, m=5)
    >>> outg1.set_attr(m=3)
    >>> nw.solve('design')
    >>> (round(outg1.fluid.val['O2'], 3), round(outg1.fluid.val['N2'], 3))
    (0.643, 0.357)
    >>> inc2.set_attr(m=np.nan)
    >>> outg1.set_attr(fluid={'O2': 0.8})
    >>> nw.solve('design')
    >>> round(inc2.m.val_SI, 3)
    1.333
    """

    def component(self):
        return 'node'

    def attr(self):
        return {'num_in': dc_simple(),
                'num_out': dc_simple()}

    def inlets(self):
        if self.num_in.val_set:
            return ['in' + str(i + 1) for i in range(self.num_in.val)]
        else:
            self.set_attr(num_in=2)
            return self.inlets()

    def outlets(self):
        if self.num_out.val_set:
            return ['out' + str(i + 1) for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=2)
            return self.outlets()

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.m_deriv = self.mass_flow_deriv()
        self.p_deriv = self.pressure_deriv()

    def equations(self):
        r"""
        Calculates vector vec_res with results of equations for this component.

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # eqation for mass flow balance
        vec_res += self.mass_flow_func()

        ######################################################################
        # equations for pressure
        inl = []
        if self.num_i > 1:
            inl = self.inl[1:]
        for c in inl + self.outl:
            vec_res += [self.inl[0].p.val_SI - c.p.val_SI]

        ######################################################################
        # additional eqations
        vec_res += self.additional_equations()

        return vec_res

    def derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivative for mass flow balance equation
        mat_deriv += self.m_deriv

        ######################################################################
        # derivatives for pressure equations
        mat_deriv += self.p_deriv

        ######################################################################
        # additional derivatives
        mat_deriv += self.additional_derivatives()

        return np.asarray(mat_deriv)

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for this
        component.

        Equations

            **mandatroy equations**

            - :func:`tespy.components.components.node.fluid_func`

            .. math::

                0 = \sum_i \left(\dot{m}_{i} \cdot h_{i}\right) - h_{o} \cdot
                \sum_i \dot{m}_{i}\\
                \forall o \in \text{outgoing mass flows}\\
                \text{i: incoming mass flows}

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # check for incoming/outgoing mass flows in inlets and outlets

        loc = 0
        # total incoming enthalpy
        h = 0
        # total incoming mass flow (constant within every iteration)
        self.m_inc = 0

        self.inc = []
        self.outg = []
        for c in self.inl:
            # incoming
            if c.m.val_SI >= 0:
                self.inc += [[c, loc]]
                self.m_inc += c.m.val_SI
                h += c.m.val_SI * c.h.val_SI
            # outgoing
            else:
                self.outg += [[c, loc]]
            loc += 1

        for c in self.outl:
            # inconming
            if c.m.val_SI < 0:
                self.inc += [[c, loc]]
                self.m_inc -= c.m.val_SI
                h -= c.m.val_SI * c.h.val_SI
            # outgoing
            else:
                self.outg += [[c, loc]]
            loc += 1

        ######################################################################
        # equations for fluid composition
        vec_res += self.fluid_func()

        ######################################################################
        # equations for energy balance
        for o in self.outg:
            vec_res += [h - o[0].h.val_SI * self.m_inc]

        return vec_res

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional
        equations.

        Returns
        -------
        mat_deriv : list
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for fluid balance equations
        mat_deriv += self.fluid_deriv()

        ######################################################################
        # derivatives for energy balance equations
        deriv = np.zeros((len(self.outg), self.num_i + self.num_o,
                          self.num_fl + 3))
        k = 0
        for o in self.outg:
            deriv[k, o[1], 2] = -self.m_inc
            for i in self.inc:
                deriv[k, i[1], 0] = i[0].h.val_SI - o[0].h.val_SI
                deriv[k, i[1], 2] = abs(i[0].m.val_SI)
            k += 1
        mat_deriv += deriv.tolist()

        return mat_deriv

    def fluid_func(self):
        r"""
        Calculates the vector of residual values for component's fluid balance
        equations.

        Returns
        -------
        vec_res : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = \sum_i \left(\dot{m}_{i} \cdot x_{i,j}\right) - x_{o,j}
                \cdot  \sum_i \dot{m}_{i}\\
                \forall j \in \text{fluids}\\
                \forall o \in \text{outgoing mass flows}\\
                \text{i: incoming mass flows}
        """
        vec_res = []

        for fluid in self.fluids:
            m = 0
            for i in self.inc:
                m += abs(i[0].m.val_SI) * i[0].fluid.val[fluid]
            for o in self.outg:
                vec_res += [m - o[0].fluid.val[fluid] * self.m_inc]
        return vec_res

    def fluid_deriv(self):
        r"""
        Calculates the partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        num_o = len(self.outg)
        deriv = np.zeros((self.num_fl * num_o, self.num_i + self.num_o,
                          3 + self.num_fl))
        j = 0
        k = 0
        for fluid in self.fluids:
            for o in self.outg:
                deriv[k, o[1], j + 3] = -self.m_inc
                for i in self.inc:
                    deriv[k, i[1], 0] = -i[0].fluid.val[fluid]
                    deriv[k, i[1], j + 3] = -abs(i[0].m.val_SI)
                k += 1
            j += 1

        return deriv.tolist()

    def pressure_deriv(self):
        r"""
        Calculates the partial derivatives for all pressure equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((self.num_i + self.num_o - 1, self.num_i + self.num_o,
                          self.num_fl + 3))

        inl = []
        if self.num_i > 1:
            inl = self.inl[1:]
        for k in range(len(inl + self.outl)):
            deriv[k, 0, 1] = 1
            deriv[k, k + 1, 1] = -1
        return deriv.tolist()

    def initialise_fluids(self, nw):
        r"""
        Fluid initialisation for fluid mixture at outlet of the node.

        Parameters
        ----------
        nw : tespy.networks.network
            Network using this component object.
        """
        num_fl = {}
        for o in self.outl:
            num_fl[o] = num_fluids(o.fluid.val)

        for i in self.inl:
            num_fl[i] = num_fluids(i.fluid.val)

        ls = []
        if any(num_fl.values()) and not all(num_fl.values()):
            for conn, num in num_fl.items():
                if num == 1:
                    ls += [conn]

            for c in ls:
                for fluid in nw.fluids:
                    for o in self.outl:
                        if not o.fluid.val_set[fluid]:
                            o.fluid.val[fluid] = c.fluid.val[fluid]
                    for i in self.inl:
                        if not i.fluid.val_set[fluid]:
                            i.fluid.val[fluid] = c.fluid.val[fluid]

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
                10^5 & \text{key = 'p'}\\
                5 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            return 5e5

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
                10^5 & \text{key = 'p'}\\
                5 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            return 5e5

# %%


class splitter(node):
    r"""
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.mass_flow_func`

        .. math::

            0 = p_{in} - p_{out,i} \;
            \forall i \in \mathrm{outlets}

        **additional equations**

        - :func:`tespy.components.components.splitter.additional_equations`

    Inlets/Outlets

        - in1
        - specify number of outlets with :code:`num_out` (default value: 2)

    Image

        .. image:: _images/split.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    num_out : float/tespy.helpers.dc_simple
        Number of outlets for this component, default value: 2.

    Example
    -------
    >>> from tespy import cmp, con, nwk, hlp
    >>> import numpy as np
    >>> fluid_list = ['O2', 'N2']
    >>> nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ...     h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> so1 = cmp.source('source1')
    >>> si1 = cmp.sink('sink1')
    >>> si2 = cmp.sink('sink2')
    >>> si3 = cmp.sink('sink3')
    >>> s = cmp.splitter('splitter', num_out=3)
    >>> s.component()
    'splitter'
    >>> inc1 = con.connection(so1, 'out1', s, 'in1')
    >>> outg1 = con.connection(s, 'out1', si1, 'in1')
    >>> outg2 = con.connection(s, 'out2', si2, 'in1')
    >>> outg3 = con.connection(s, 'out3', si3, 'in1')
    >>> nw.add_conns(inc1, outg1, outg2, outg3)
    >>> inc1.set_attr(fluid={'O2': 0.23, 'N2': 0.77}, p=1, T=20, m=5)
    >>> outg1.set_attr(m=3)
    >>> outg2.set_attr(m=1)
    >>> nw.solve('design')
    >>> nw.lin_dep
    False
    >>> nw.res[-1] < 1e-3
    True
    """

    def component(self):
        return 'splitter'

    def attr(self):
        return {'num_out': dc_simple()}

    def inlets(self):
        return ['in1']

    def outlets(self):
        if self.num_out.val_set:
            return ['out' + str(i + 1) for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=2)
            return self.outlets()

    def comp_init(self, nw):

        node.comp_init(self, nw)

        self.fl_deriv = self.fluid_deriv()
        self.h_deriv = self.enthalpy_deriv()

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for
        this component.

        Equations

            **mandatroy equations**

            .. math:: 0 = fluid_{i,in} - fluid_{i,out_{j}} \;
                \forall i \in \mathrm{fluid}, \; \forall j \in outlets

            .. math::
                0 = h_{in} - h_{out,i} \;
                \forall i \in \mathrm{outlets}\\

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equations for fluid balance
        for o in self.outl:
            for fluid, x in self.inl[0].fluid.val.items():
                vec_res += [x - o.fluid.val[fluid]]

        ######################################################################
        # equations for energy balance
        for o in self.outl:
            vec_res += [self.inl[0].h.val_SI - o.h.val_SI]

        return vec_res

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional
        equations.

        Returns
        -------
        mat_deriv : list
            Matrix of partial derivatives.
        """
        ######################################################################
        # derivatives for fluid and energy balance equations are constant
        return self.fl_deriv + self.h_deriv

    def fluid_deriv(self):
        r"""
        Calculates the partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((self.num_fl * self.num_o, 1 + self.num_o,
                          3 + self.num_fl))
        k = 0
        for o in self.outl:
            i = 0
            for fluid in self.fluids:
                deriv[i + k * self.num_fl, 0, i + 3] = 1
                deriv[i + k * self.num_fl, k + 1, i + 3] = -1
                i += 1
            k += 1
        return deriv.tolist()

    def enthalpy_deriv(self):
        r"""
        Calculates matrix of partial derivatives for enthalpy balance equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((self.num_o, 1 + self.num_o, self.num_fl + 3))
        k = 0
        for o in self.outl:
            deriv[k, 0, 2] = 1
            deriv[k, k + 1, 2] = -1
            k += 1

        return deriv.tolist()

    def initialise_fluids(self, nw):
        r"""
        Fluid initialisation for fluid mixture at outlet of the node.

        Parameters
        ----------
        nw : tespy.networks.network
            Network using this component object.
        """
        return

# %%


class separator(node):
    r"""
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.mass_flow_func`

        .. math::

            0 = p_{in} - p_{out,i} \;
            \forall i \in \mathrm{outlets}

        **additional equations**

        - :func:`tespy.components.components.separator.additional_equations`

    Inlets/Outlets

        - in1
        - specify number of outlets with :code:`num_out` (default value: 2)

    Image

        .. image:: _images/split.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    TODO

        - fluid separation requires power and cooling, equations have not
          been implemented!

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    num_out : float/tespy.helpers.dc_simple
        Number of outlets for this component, default value: 2.

    Example
    -------
    >>> from tespy import cmp, con, nwk, hlp
    >>> import numpy as np
    >>> fluid_list = ['O2', 'N2']
    >>> nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ...     h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> so1 = cmp.source('source1')
    >>> si1 = cmp.sink('sink1')
    >>> si2 = cmp.sink('sink2')
    >>> s = cmp.separator('separator', num_out=2)
    >>> s.component()
    'separator'
    >>> inc1 = con.connection(so1, 'out1', s, 'in1')
    >>> outg1 = con.connection(s, 'out1', si1, 'in1')
    >>> outg2 = con.connection(s, 'out2', si2, 'in1')
    >>> nw.add_conns(inc1, outg1, outg2)
    >>> inc1.set_attr(fluid={'O2': 0.23, 'N2': 0.77}, p=1, T=20, m=5)
    >>> outg1.set_attr(fluid={'O2': 0.1, 'N2': 0.9}, m=1)
    >>> outg2.set_attr(fluid0={'O2': 0.5, 'N2': 0.5}, m0=4)
    >>> nw.solve('design')
    >>> nw.lin_dep
    False
    >>> nw.res[-1] < 1e-3
    True
    """

    def component(self):
        return 'separator'

    def attr(self):
        return {'num_out': dc_simple()}

    def inlets(self):
        return ['in1']

    def outlets(self):
        if self.num_out.val_set:
            return ['out' + str(i + 1) for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=2)
            return self.outlets()

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for
        this component.

        Equations

            **mandatroy equations**

            .. math:: 0 = fluid_{i,in} - fluid_{i,out_{j}} \;
                \forall i \in \mathrm{fluid}, \; \forall j \in outlets

            .. math::

                0 = T_{in} - T_{out,i} \;
                \forall i \in \mathrm{outlets}

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equations for fluid balance
        for fluid, x in self.inl[0].fluid.val.items():
            res = x * self.inl[0].m.val_SI
            for o in self.outl:
                res -= o.fluid.val[fluid] * o.m.val_SI
            vec_res += [res]

        ######################################################################
        # equations for energy balance
        for o in self.outl:
            vec_res += [
                    T_mix_ph(self.inl[0].to_flow(), T0=self.inl[0].T.val_SI) -
                    T_mix_ph(o.to_flow(), T0=o.T.val_SI)]

        return vec_res

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional
        equations.

        Returns
        -------
        mat_deriv : list
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for fluid balance equations
        mat_deriv += self.fluid_deriv()

        ######################################################################
        # derivatives for energy balance equations
        deriv = np.zeros((self.num_o, 1 + self.num_o, self.num_fl + 3))
        i = self.inl[0].to_flow()
        k = 0
        for o in self.outl:
            o = o.to_flow()
            deriv[k, 0, 1] = dT_mix_dph(i)
            deriv[k, 0, 2] = dT_mix_pdh(i)
            deriv[k, 0, 3:] = dT_mix_ph_dfluid(i)
            deriv[k, k + 1, 1] = -dT_mix_dph(o)
            deriv[k, k + 1, 2] = -dT_mix_pdh(o)
            deriv[k, k + 1, 3:] = -1 * dT_mix_ph_dfluid(o)
            k += 1
        mat_deriv += deriv.tolist()

        return mat_deriv

    def fluid_deriv(self):
        r"""
        Calculates the partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((self.num_fl, 1 + self.num_o, 3 + self.num_fl))
        j = 0
        for fluid in self.fluids:
            k = 0
            for o in self.outl:
                deriv[j, k + 1, 0] = -o.fluid.val[fluid]
                deriv[j, k + 1, j + 3] = -o.m.val_SI
                k += 1
            deriv[j, 0, 0] = self.inl[0].fluid.val[fluid]
            deriv[j, 0, j + 3] = self.inl[0].m.val_SI
            j += 1
        return deriv.tolist()

    def initialise_fluids(self, nw):
        r"""
        Fluid initialisation for fluid mixture at outlet of the node.

        Parameters
        ----------
        nw : tespy.networks.network
            Network using this component object.
        """
        return

# %%


class merge(node):
    r"""
    The component node is the parent class for splitter, separator and merge.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.mass_flow_func`

        .. math::

            0 = p_{in} - p_{out,i} \;
            \forall i \in \mathrm{outlets}

        **additional equations**

        - :func:`tespy.components.components.merge.additional_equations`

    Inlets/Outlets

        - specify number of outlets with :code:`num_in` (default value: 2)
        - out1

    Image

        .. image:: _images/merge.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    num_in : float/tespy.helpers.dc_simple
        Number of inlets for this component, default value: 2.

    Example
    -------
    >>> from tespy import cmp, con, nwk, hlp
    >>> fluid_list = ['O2', 'N2']
    >>> nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ...     h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> so1 = cmp.source('source1')
    >>> so2 = cmp.source('source2')
    >>> so3 = cmp.source('source3')
    >>> si1 = cmp.sink('sink1')
    >>> m = cmp.merge('merge', num_in=3)
    >>> m.component()
    'merge'
    >>> inc1 = con.connection(so1, 'out1', m, 'in1')
    >>> inc2 = con.connection(so2, 'out1', m, 'in2')
    >>> inc3 = con.connection(so3, 'out1', m, 'in3')
    >>> outg1 = con.connection(m, 'out1', si1, 'in1')
    >>> nw.add_conns(inc1, inc2, inc3, outg1)
    >>> inc1.set_attr(fluid={'O2': 0.23, 'N2': 0.77}, p=1, T=20, m=5)
    >>> inc2.set_attr(fluid={'O2': 1, 'N2':0}, T=20, m=5)
    >>> inc3.set_attr(fluid={'O2': 0, 'N2': 1}, T=20)
    >>> outg1.set_attr(fluid={'N2': 0.4})
    >>> nw.solve('design')
    >>> round(inc3.m.val_SI, 2)
    0.25
    >>> round(outg1.fluid.val['O2'], 1)
    0.6
    """

    def component(self):
        return 'merge'

    def attr(self):
        return {'num_in': dc_simple(),
                'zero_flag': dc_simple()}

    def inlets(self):
        if self.num_in.val_set:
            return ['in' + str(i + 1) for i in range(self.num_in.val)]
        else:
            self.set_attr(num_in=2)
            return self.inlets()

    def outlets(self):
        return ['out1']

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for
        this component.

        Equations

            **mandatroy equations**

            .. math::

                0 = \dot{m}_{in_{j}} \cdot fluid_{i,in_{j}} -
                    \dot {m}_{out} \cdot fluid_{i,out} \\
                \forall i \in \mathrm{fluid}, \; \forall j \in inlets

            .. math::

                0 = h_{in} - h_{out,i} \;
                \forall i \in \mathrm{outlets}\\

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equations for fluid balance
        for fluid, x in self.outl[0].fluid.val.items():
            res = -x * self.outl[0].m.val_SI
            for i in self.inl:
                res += i.fluid.val[fluid] * i.m.val_SI
            vec_res += [res]

        ######################################################################
        # equation for energy balance
        h_res = -self.outl[0].m.val_SI * self.outl[0].h.val_SI
        for i in self.inl:
            h_res += i.m.val_SI * i.h.val_SI
        vec_res += [h_res]

        return vec_res

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional
        equations.

        Returns
        -------
        mat_deriv : list
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for fluid balance equations
        mat_deriv += self.fluid_deriv()

        ######################################################################
        # derivatives for energy balance equations
        deriv = np.zeros((1, self.num_i + 1, self.num_fl + 3))
        deriv[0, self.num_i, 0] = -self.outl[0].h.val_SI
        deriv[0, self.num_i, 2] = -self.outl[0].m.val_SI
        k = 0
        for i in self.inl:
            deriv[0, k, 0] = i.h.val_SI
            deriv[0, k, 2] = i.m.val_SI
            k += 1
        mat_deriv += deriv.tolist()

        return mat_deriv

    def fluid_deriv(self):
        r"""
        Calculates the partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((self.num_fl, self.num_i + 1, 3 + self.num_fl))
        j = 0
        for fluid, x in self.outl[0].fluid.val.items():
            k = 0
            for i in self.inl:
                deriv[j, k, 0] = i.fluid.val[fluid]
                deriv[j, k, j + 3] = i.m.val_SI
                k += 1
            deriv[j, k, 0] = -x
            deriv[j, k, j + 3] = -self.outl[0].m.val_SI
            j += 1
        return deriv.tolist()

# %%


class water_electrolyzer(component):
    r"""
    Equations

        **mandatory equations**

        .. math::

            0  = x_{i,in1} - x_{i,out1} \forall i \in \text{fluids}

            \forall i \in \text{network fluids}:

            0 = \begin{cases}
                1 - x_{i,in2} & \text{i=}H_{2}O\\
                x_{i,in2} & \text{else}
            \end{cases}\\

            0 = \begin{cases}
                1 - x_{i,out2} & \text{i=}O_{2}\\
                x_{i,out2} & \text{else}
            \end{cases}\\

            0 = \begin{cases}
                1 - x_{i,out3} & \text{i=}H_{2}\\
                x_{i,out3} & \text{else}
            \end{cases}\\

            O_2 = \frac{M_{O_2}}{M_{O_2} + 2 \cdot M_{H_2}}\\

            0 = \dot{m}_{H_{2}O,in1} - \dot{m}_{H_{2}O,out1}\\
            0 = O_2 \cdot \dot{m}_{H_{2}O,in2} - \dot{m}_{O_2,out2}\\
            0 = \left(1 - O_2\right) \cdot \dot{m}_{H_{2}O,in2} -
            \dot{m}_{H_2,out3}\\

            0 = p_{H_{2}O,in2} - p_{O_2,out2}\\
            0 = p_{H_{2}O,in2} - p_{H_2,out3}

            0 = P - f_{eb}\left( \right)

        For energy balance (f_eb) calculation see
        :func:`tespy.components.components.water_electrolyzer.energy_balance`.

        .. math::

            0 = T_{O_2,out2} - T_{H_2,out3}

        **optional equations**

        .. math::

            0 = P - \dot{m}_{H_2,out3} \cdot e\\

            0 = p_{H_{2}O,in1} \cdot pr - p_{H_{2}O,out1}

        - :func:`tespy.components.components.component.zeta_func`

        .. math::

            0 = \dot{Q} - \dot{m}_{in1} \cdot \left(h_{in1} -
            h_{out1}\right) \\

            0 = P - \dot{m}_{H_2,out3} \cdot \frac{e_0}{\eta}

        - :func:`tespy.components.components.water_electrolyzer.eta_char_func`

    Inlets/Outlets

        - in1 (cooling inlet), in2 (feed water inlet)
        - out1 (cooling outlet), out2 (hydrogen outlet), out3 (oxigen outlet)

    Image

        .. image:: _images/water_electrolyzer.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    P : float/tespy.helpers.dc_cp
        Power input, :math:`P/\text{W}`.

    Q : float/tespy.helpers.dc_cp
        Heat output of cooling, :math:`Q/\text{W}`

    e : float/tespy.helpers.dc_cp
        Electrolysis specific energy consumption,
        :math:`e/(\text{J}/\text{m}^3)`.

    eta : float/tespy.helpers.dc_cp
        Electrolysis efficiency, :math:`\eta/1`.

    eta_char : str/tespy.helpers.dc_cc
        Electrolysis efficiency characteristic line.

    pr : float/tespy.helpers.dc_cp
        Cooling loop pressure ratio, :math:`pr/1`.

    zeta : float/tespy.helpers.dc_cp
        Geometry independent friction coefficient for cooling loop pressure
        drop, :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    Note
    ----

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> import shutil
    >>> fluid_list = ['O2', 'water', 'H2']
    >>> nw = nwk.network(fluids=fluid_list, T_unit='C', p_unit='bar',
    ... h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')

    >>> fw = cmp.source('feed water')
    >>> oxy = cmp.sink('oxygen sink')
    >>> hydro = cmp.sink('hydrogen sink')
    >>> cw = cmp.source('cooling water')
    >>> cw_hot = cmp.sink('cooling water out')

    >>> el = cmp.water_electrolyzer('electrolyzer 1', eta=0.8, design=['eta'],
    ... offdesign=['eta_char'])
    >>> el.component()
    'water electrolyzer'
    >>> comp = cmp.compressor('compressor', eta_s=0.9)

    >>> fw_el = con.connection(fw, 'out1', el, 'in2', m=0.1, p=10, T=15)
    >>> el_o = con.connection(el, 'out2', oxy, 'in1')
    >>> el_cmp = con.connection(el, 'out3', comp, 'in1', T=50)
    >>> cmp_h = con.connection(comp, 'out1', hydro, 'in1', p=50)
    >>> cw_el = con.connection(cw, 'out1', el, 'in1', p=5, T=15,
    ... fluid={'water': 1, 'H2': 0, 'O2': 0})
    >>> el_cw = con.connection(el, 'out1', cw_hot, 'in1', T=45, p=4.9)
    >>> nw.add_conns(fw_el, el_o, el_cmp, cmp_h, cw_el, el_cw)
    >>> nw.solve('design')
    >>> round(el.eta.val, 1)
    0.8
    >>> nw.save('tmp')
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(el.eta.val, 1)
    0.8
    >>> fw_el.set_attr(m=0.05)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(el.eta.val, 2)
    0.82
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'water electrolyzer'

    def attr(self):
        return {'P': dc_cp(min_val=0),
                'Q': dc_cp(max_val=0),
                'eta': dc_cp(min_val=0, max_val=1),
                'e': dc_cp(),
                'pr_c': dc_cp(max_val=1),
                'zeta': dc_cp(min_val=0),
                'eta_char': dc_cc(method='GENERIC'),
                'S': dc_simple()}

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1', 'out2', 'out3']

    def comp_init(self, nw):

        if not self.P.is_set:
            self.set_attr(P='var')
            msg = ('The power output of cogeneration units must be set! '
                   'We are adding the power output of component ' +
                   self.label + ' as custom variable of the system.')
            logging.info(msg)

        component.comp_init(self, nw)

        o2 = [x for x in nw.fluids if x in [a.replace(' ', '')
              for a in CP.get_aliases('O2')]]
        if len(o2) == 0:
            msg = ('Missing oxygen in network fluids, component ' +
                   self.label + ' of type ' + self.component() +
                   ' requires oxygen in network fluids.')
            logging.error(msg)
            raise ValueError(msg)
        else:
            self.o2 = o2[0]

        h2o = [x for x in nw.fluids if x in [a.replace(' ', '')
               for a in CP.get_aliases('H2O')]]
        if len(h2o) == 0:
            msg = ('Missing water in network fluids, component ' +
                   self.label + ' of type ' + self.component() +
                   ' requires water in network fluids.')
            logging.error(msg)
            raise ValueError(msg)
        else:
            self.h2o = h2o[0]

        h2 = [x for x in nw.fluids if x in [a.replace(' ', '')
              for a in CP.get_aliases('H2')]]
        if len(h2) == 0:
            msg = ('Missing hydrogen in network fluids, component ' +
                   self.label + ' of type ' + self.component() +
                   ' requires hydrogen in network fluids.')
            logging.error(msg)
            raise ValueError(msg)
        else:
            self.h2 = h2[0]

        self.e0 = self.calc_e0()

    def calc_e0(self):
        r"""
        Calculates the minimum specific energy required for electrolysis.

        Returns
        -------
        val : float
            Minimum specific energy.

            .. math::
                LHV = -\frac{\sum_i {\Delta H_f^0}_i -
                \sum_j {\Delta H_f^0}_j }
                {M_{fuel}}\\
                \forall i \in \text{reation products},\\
                \forall j \in \text{reation educts},\\
                \Delta H_f^0: \text{molar formation enthalpy}
        """

        hf = {}
        hf['H2O'] = -286
        hf['H2'] = 0
        hf['O2'] = 0
        M = molar_masses['H2']
        e0 = -(2 * hf['H2O'] - 2 * hf['H2'] + hf['O2']) / (2 * M)

        return e0 * 1000

    def equations(self):
        r"""
        Calculates vector vec_res with results of equations for this component.

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equations for fluids

        # equations for fluid composition in cooling water
        for fluid, x in self.inl[0].fluid.val.items():
            vec_res += [x - self.outl[0].fluid.val[fluid]]

        # equations to constrain fluids to inlets/outlets
        vec_res += [1 - self.inl[1].fluid.val[self.h2o]]
        vec_res += [1 - self.outl[1].fluid.val[self.o2]]
        vec_res += [1 - self.outl[2].fluid.val[self.h2]]

        # equations to ban fluids off inlets/outlets
        for fluid in self.inl[1].fluid.val.keys():
            if fluid != self.h2o:
                vec_res += [0 - self.inl[1].fluid.val[fluid]]
            if fluid != self.o2:
                vec_res += [0 - self.outl[1].fluid.val[fluid]]
            if fluid != self.h2:
                vec_res += [0 - self.outl[2].fluid.val[fluid]]

        ######################################################################
        # eqations for mass flow balance
        # equation to calculate the ratio of o2 in water
        o2 = molar_masses[self.o2] / (molar_masses[self.o2] +
                                      2 * molar_masses[self.h2])

        # equation for mass flow balance cooling water
        vec_res += [self.inl[0].m.val_SI - self.outl[0].m.val_SI]

        # equations for mass flow balance electrolyzer
        vec_res += [o2 * self.inl[1].m.val_SI - self.outl[1].m.val_SI]
        vec_res += [(1 - o2) * self.inl[1].m.val_SI - self.outl[2].m.val_SI]

        ######################################################################
        # equations for pressure to set o2 and h2 output equal
        vec_res += [self.inl[1].p.val_SI - self.outl[1].p.val_SI]
        vec_res += [self.inl[1].p.val_SI - self.outl[2].p.val_SI]

        ######################################################################
        # equation for energy balance
        vec_res += [self.P.val + self.energy_balance()]

        ######################################################################
        # temperature electrolyzer outlet
        vec_res += [T_mix_ph(self.outl[1].to_flow()) -
                    T_mix_ph(self.outl[2].to_flow())]

        ######################################################################
        # power vs hydrogen production
        if self.e.is_set:
            vec_res += [self.P.val - self.outl[2].m.val_SI * self.e.val]

        ######################################################################
        #pr_c.val = pressure ratio Druckverlust (als Faktor vorgegeben)
        if self.pr_c.is_set:
            vec_res += [self.inl[0].p.val_SI * self.pr_c.val -
                        self.outl[0].p.val_SI]

        if self.zeta.is_set:
            vec_res += [self.zeta_func()]

        # equation for heat transfer

        if self.Q.is_set:
            vec_res += [self.Q.val - self.inl[0].m.val_SI *
                        (self.inl[0].h.val_SI - self.outl[0].h.val_SI)]

        ######################################################################
        # specified efficiency (efficiency definition: e0 / e)
        if self.eta.is_set:
            vec_res += [self.P.val - self.outl[2].m.val_SI *
                        self.e0 / self.eta.val]

        ######################################################################
        # specified characteristic line for efficiency
        if self.eta_char.is_set:
            vec_res += [self.eta_char_func()]

        return vec_res

    def derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for cooling liquid composition
        deriv = np.zeros((self.num_fl, 5 + self.num_vars, self.num_fl + 3))

        j = 0
        for fluid, x in self.inl[0].fluid.val.items():
            deriv[j, 0, 3 + j] = 1
            deriv[j, 2, 3 + j] = -1
            j += 1

        mat_deriv += deriv.tolist()

        # derivatives to constrain fluids to inlets/outlets
        deriv = np.zeros((3, 5 + self.num_vars, self.num_fl + 3))

        i = 0
        for fluid in self.inl[0].fluid.val.keys():
            if fluid == self.h2o:
                deriv[0, 1, 3 + i] = -1
            elif fluid == self.o2:
                deriv[1, 3, 3 + i] = -1
            elif fluid == self.h2:
                deriv[2, 4, 3 + i] = -1
            i += 1

        mat_deriv += deriv.tolist()

        # derivatives to ban fluids off inlets/outlets
        deriv = np.zeros((3 * len(self.inl[1].fluid.val.keys()) - 3,
                          5 + self.num_vars, self.num_fl + 3))

        i = 0
        j = 0
        for fluid in self.inl[1].fluid.val.keys():
            if fluid != self.h2o:
                deriv[j, 1, 3 + i] = -1
                j += 1
            if fluid != self.o2:
                deriv[j, 3, 3 + i] = -1
                j += 1
            if fluid != self.h2:
                deriv[j, 4, 3 + i] = -1
                j += 1
            i += 1

        mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for mass balance equations

        # deritatives for mass flow balance in the heat exchanger
        deriv = np.zeros((3, 5 + self.num_vars, self.num_fl + 3))

        deriv[0, 0, 0] = 1
        deriv[0, 2, 0] = -1

        # derivatives for mass flow balance for oxygen output
        o2 = molar_masses[self.o2] / (molar_masses[self.o2] +
                                      2 * molar_masses[self.h2])
        deriv[1, 1, 0] = o2
        deriv[1, 3, 0] = -1

        # derivatives for mass flow balance for hydrogen output
        deriv[2, 1, 0] = (1 - o2)
        deriv[2, 4, 0] = -1

        mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for pressure equations

        # derivatives for pressure oxygen outlet
        deriv = np.zeros((2, 5 + self.num_vars, self.num_fl + 3))

        deriv[0, 1, 1] = 1
        deriv[0, 3, 1] = -1

        # derivatives for pressure hydrogen outlet
        deriv[1, 1, 1] = 1
        deriv[1, 4, 1] = -1

        mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for energy balance equations

        deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

        T_ref = 293.15
        p_ref = 1e5

        h_refh2o = h_mix_pT([1, p_ref, 0, self.inl[1].fluid.val], T_ref)
        h_refh2 = h_mix_pT([1, p_ref, 0, self.outl[2].fluid.val], T_ref)
        h_refo2 = h_mix_pT([1, p_ref, 0, self.outl[1].fluid.val], T_ref)

        # derivatives cooling water inlet
        deriv[0, 0, 0] = - (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        deriv[0, 0, 2] = self.inl[0].m.val_SI

        # derivatives feed water inlet
        deriv[0, 1, 0] = (self.inl[1].h.val_SI - h_refh2o)
        deriv[0, 1, 2] = self.inl[1].m.val_SI

        # derivative cooling water outlet
        deriv[0, 2, 2] = - self.inl[0].m.val_SI

        # derivatives oxygen outlet
        deriv[0, 3, 0] = - (self.outl[1].h.val_SI - h_refo2)
        deriv[0, 3, 2] = - self.outl[1].m.val_SI

        # derivatives hydrogen outlet
        deriv[0, 4, 0] = - self.e0 - (self.outl[2].h.val_SI - h_refh2)
        deriv[0, 4, 2] = - self.outl[2].m.val_SI

        # derivatives for variable P
        if self.P.is_var:
            deriv[0, 5 + self.P.var_pos, 0] = 1

        mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for temperature at gas outlets

        deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

        # derivatives for outlet 1
        deriv[0, 3, 1] = dT_mix_dph(self.outl[1].to_flow())
        deriv[0, 3, 2] = dT_mix_pdh(self.outl[1].to_flow())

        # derivatives for outlet 2
        deriv[0, 4, 1] = - dT_mix_dph(self.outl[2].to_flow())
        deriv[0, 4, 2] = - dT_mix_pdh(self.outl[2].to_flow())

        mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for power vs. hydrogen production

        if self.e.is_set:
            deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

            deriv[0, 4, 0] = - self.e.val

            # derivatives for variable P
            if self.P.is_var:
                deriv[0, 5 + self.P.var_pos, 0] = 1

            # derivatives for variable e
            if self.e.is_var:
                deriv[0, 5 + self.e.var_pos, 0] = - self.outl[2].m.val_SI

            mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for pressure ratio
        if self.pr_c.is_set:

            deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

            deriv[0, 0, 1] = self.pr_c.val
            deriv[0, 2, 1] = - 1

            mat_deriv += deriv.tolist()

        ######################################################################
        #pr_c.val = pressure ratio Druckverlust (als Faktor vorgegeben)
        # derivatives for zeta value
        if self.zeta.is_set:

            deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))
            deriv[0, 0, 0] = self.numeric_deriv(self.zeta_func, 'm', 0)
            deriv[0, 0, 1] = self.numeric_deriv(self.zeta_func, 'p', 0)
            deriv[0, 0, 2] = self.numeric_deriv(self.zeta_func, 'h', 0)
            deriv[0, 2, 1] = self.numeric_deriv(self.zeta_func, 'p', 2)
            deriv[0, 2, 2] = self.numeric_deriv(self.zeta_func, 'h', 2)

            # derivatives for variable zeta
            if self.zeta.is_var:
                deriv[0, 5 + self.zeta.var_pos, 0] = (
                        self.numeric_deriv(self.zeta_func, 'zeta', 5))

            mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for heat flow
        if self.Q.is_set:

            deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

            deriv[0, 0, 0] = - (self.inl[0].h.val_SI - self.outl[0].h.val_SI)
            deriv[0, 0, 2] = - self.inl[0].m.val_SI
            deriv[0, 2, 2] = self.inl[0].m.val_SI

            mat_deriv += deriv.tolist()

        ######################################################################
        # specified efficiency (efficiency definition: e0 / e)
        if self.eta.is_set:

            deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

            deriv[0, 4, 0] = - self.e0 / self.eta.val

            # derivatives for variable P
            if self.P.is_var:
                deriv[0, 5 + self.P.var_pos, 0] = 1

            mat_deriv += deriv.tolist()

        ######################################################################
        # specified characteristic line for efficiency
        if self.eta_char.is_set:

            mat_deriv += self.eta_char_deriv()

        ######################################################################

        return np.asarray(mat_deriv)

    def eta_char_func(self):
        r"""
        Equation for given efficiency characteristic of a water electrolyzer.
        Efficiency is linked to hydrogen production.

        Returns
        -------
        res : ndarray
            Residual value of equation.

            .. math::

                0 = P - \dot{m}_{H_2,out3} \cdot \frac{e_0}{\eta_0 \cdot
                f\left(\frac{\dot{m}_{H_2,out3}}{\dot{m}_{H_2,out3,0}} \right)}
        """
        expr = self.outl[2].m.val_SI / self.outl[2].m.design

        return (self.P.val - self.outl[2].m.val_SI * self.e0 / (
                self.eta.design * self.eta_char.func.f_x(expr)))

    def eta_char_deriv(self):
        r"""
        Calculates the matrix of partial derivatives of the efficiency
        characteristic function.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """

        deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

        deriv[0, 4, 0] = self.numeric_deriv(self.eta_char_func, 'm', 4)

        # derivatives for variable P
        if self.P.is_var:
            deriv[0, 5 + self.P.var_pos, 0] = 1

        return deriv.tolist()

    def bus_func(self, bus):
        r"""
        Calculates the residual value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        val : float
            Residual value of equation.

            .. math::

                val = \begin{cases}
                P \cdot f_{char}\left( \frac{P}{P_{ref}}\right) &
                \text{key = 'P'}\\
                \dot{Q} \cdot f_{char}\left( \frac{\dot{Q}}
                {\dot{Q}_{ref}}\right) & \text{key = 'Q'}\\
                \end{cases}\\
                \dot{Q} = - \dot{m}_{1,in} \cdot
                \left(h_{out,1} - h_{in,1} \right)\\
        """
        ######################################################################
        # equations for power on bus
        if bus.param == 'P':
            P = - self.energy_balance()
            if np.isnan(bus.P_ref):
                expr = 1
            else:
                expr = abs(P / bus.P_ref)
            return P * bus.char.f_x(expr)

        ######################################################################
        # equations for heat on bus

        elif bus.param == 'Q':
            val = - self.inl[0].m.val_SI * (self.outl[0].h.val_SI -
                                            self.inl[0].h.val_SI)
            if np.isnan(bus.P_ref):
                expr = 1
            else:
                expr = abs(val / bus.P_ref)
            return val * bus.char.f_x(expr)

        ######################################################################
        # missing/invalid bus parameter

        else:
            msg = ('The parameter ' + str(bus.param) + ' is not a valid '
                   'parameter for a component of type ' + self.component() +
                   '. Please specify a bus parameter (P/Q) for component ' +
                   self.label + '.')
            logging.error(msg)
            raise ValueError(msg)

    def bus_deriv(self, bus):
        r"""
        Calculates the matrix of partial derivatives of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

        ######################################################################
        # derivatives for power on bus
        if bus.param == 'P':
            deriv[0, 0, 0] = self.numeric_deriv(self.bus_func, 'm', 0, bus=bus)
            deriv[0, 0, 2] = self.numeric_deriv(self.bus_func, 'h', 0, bus=bus)

            deriv[0, 1, 0] = self.numeric_deriv(self.bus_func, 'm', 1, bus=bus)
            deriv[0, 1, 2] = self.numeric_deriv(self.bus_func, 'h', 1, bus=bus)

            deriv[0, 2, 2] = self.numeric_deriv(self.bus_func, 'h', 2, bus=bus)

            deriv[0, 3, 0] = self.numeric_deriv(self.bus_func, 'm', 3, bus=bus)
            deriv[0, 3, 2] = self.numeric_deriv(self.bus_func, 'h', 3, bus=bus)

            deriv[0, 4, 0] = self.numeric_deriv(self.bus_func, 'm', 4, bus=bus)
            deriv[0, 4, 2] = self.numeric_deriv(self.bus_func, 'h', 4, bus=bus)
            # variable power
            if self.P.is_var:
                deriv[0, 5 + self.P.var_pos, 0] = (
                        self.numeric_deriv(self.bus_func, 'P', 5, bus=bus))

        ######################################################################
        # derivatives for heat on bus
        elif bus.param == 'Q':

            deriv = np.zeros((1, 5 + self.num_vars, self.num_fl + 3))

            deriv[0, 0, 0] = self.numeric_deriv(self.bus_func, 'm', 0, bus=bus)
            deriv[0, 0, 2] = self.numeric_deriv(self.bus_func, 'h', 0, bus=bus)
            deriv[0, 2, 2] = self.numeric_deriv(self.bus_func, 'h', 2, bus=bus)

        ######################################################################
        # missing/invalid bus parameter

        else:
            msg = ('The parameter ' + str(bus.param) + ' is not a valid '
                   'parameter for a component of type ' + self.component() +
                   '. Please specify a bus parameter (P/Q) for component ' +
                   self.label + '.')
            logging.error(msg)
            raise ValueError(msg)

        return deriv

    def energy_balance(self):
        r"""
        Calculates the residual in energy balance of the adiabatic water
        electrolyzer. The residual is the negative to the necessary power
        input.

        Returns
        -------
        res : float
            Residual value.

            .. math::

                \begin{split}
                res = & \dot{m}_{in,2} \cdot \left( h_{in,2} - h_{in,2,ref}
                \right)\\ & - \dot{m}_{out,3} \cdot e_0\\
                & -\dot{m}_{in,1} \cdot \left( h_{out,1} - h_{in,1} \right)\\
                & - \dot{m}_{out,2} \cdot \left( h_{out,2} - h_{out,2,ref}
                \right)\\
                & - \dot{m}_{out,3} \cdot \left( h_{out,3} - h_{out,3,ref}
                \right)\\
                \end{split}

        Note
        ----
        The temperature for the reference state is set to 20 °C, thus
        the feed water must be liquid as proposed in the calculation of
        the minimum specific energy consumption for electrolysis:
        :func:`tespy.components.components.water_electrolyzer.calc_e0`.
        The part of the equation regarding the cooling water is implemented
        with negative sign as the energy for cooling is extracted from the
        reactor.

        - Reference temperature: 293.15 K.
        - Reference pressure: 1 bar.
        """
        T_ref = 293.15
        p_ref = 1e5

        # equations to set a reference point for each h2o, h2 and o2
        h_refh2o = h_mix_pT([1, p_ref, 0, self.inl[1].fluid.val], T_ref)
        h_refh2 = h_mix_pT([1, p_ref, 0, self.outl[2].fluid.val], T_ref)
        h_refo2 = h_mix_pT([1, p_ref, 0, self.outl[1].fluid.val], T_ref)

        val = (self.inl[1].m.val_SI * (self.inl[1].h.val_SI - h_refh2o) -
               self.outl[2].m.val_SI * self.e0 -
               self.inl[0].m.val_SI * (self.outl[0].h.val_SI -
                                       self.inl[0].h.val_SI) -
               self.outl[1].m.val_SI * (self.outl[1].h.val_SI - h_refo2) -
               self.outl[2].m.val_SI * (self.outl[2].h.val_SI - h_refh2))
        return val

    def initialise_fluids(self, nw):
        r"""
        Sets values to pure fluid on water inlet and gas outlets.

        Parameters
        ----------
        nw : tespy.networks.network
            Network using this component object.
        """
        self.outl[1].fluid.val[self.o2] = 1
        self.outl[2].fluid.val[self.h2] = 1
        self.inl[1].fluid.val[self.h2o] = 1

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
                5  \cdot 10^5 & \text{key = 'p'}\\
                h\left(T=323.15, p=5  \cdot 10^5\right) & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            flow = [c.m.val0, 5e5, c.h.val_SI, c.fluid.val]
            T = 50 + 273.15
            return h_mix_pT(flow, T)

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
                5  \cdot 10^5 & \text{key = 'p'}\\
                h\left(T=293.15, p=5  \cdot 10^5\right) & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            flow = [c.m.val0, 5e5, c.h.val_SI, c.fluid.val]
            T = 20 + 273.15
            return h_mix_pT(flow, T)

    def calc_parameters(self):
        r"""
        Postprocessing parameter calculation.
        """
        self.Q.val = - self.inl[0].m.val_SI * (self.outl[0].h.val_SI -
                                               self.inl[0].h.val_SI)
        self.pr_c.val = self.outl[0].p.val_SI / self.inl[0].p.val_SI
        self.e.val = self.P.val / self.outl[2].m.val_SI
        self.eta.val = self.e0 / self.e.val

        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()
        self.zeta.val = ((i[1] - o[1]) * math.pi ** 2 /
                         (8 * i[0] ** 2 * (v_mix_ph(i) + v_mix_ph(o)) / 2))

        if self.eta_char.is_set:
            # get bound errors for kA hot side characteristics
            expr = self.outl[2].m.val_SI / self.outl[2].m.design
            self.eta_char.func.get_bound_errors(expr, self.label)

        self.check_parameter_bounds()

# %%


class valve(component):
    r"""
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        .. math::

            0 = h_{in} - h_{out}

        **optional equations**

        .. math::

            0 = p_{in} \cdot pr - p_{out}

        - :func:`tespy.components.components.component.zeta_func`

    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/valve.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    pr : str/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`

    zeta : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    Example
    -------
    >>> from tespy import cmp, con, nwk, hlp
    >>> import shutil
    >>> import numpy as np
    >>> fluid_list = ['CH4']
    >>> nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C')
    >>> nw.set_printoptions(print_level='none')
    >>> so = cmp.source('source')
    >>> si = cmp.sink('sink')
    >>> v = cmp.valve('valve')
    >>> v.component()
    'valve'
    >>> so_v = con.connection(so, 'out1', v, 'in1')
    >>> v_si = con.connection(v, 'out1', si, 'in1')
    >>> nw.add_conns(so_v, v_si)
    >>> v.set_attr(pr=0.05, design=['pr'], offdesign=['zeta'])
    >>> so_v.set_attr(fluid={'CH4': 1}, m=10)
    >>> v_si.set_attr(p=2, T=10)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(v.zeta.val, 1)
    122239.1
    >>> so_v.set_attr(m=12)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(v.pr.val, 3)
    0.036
    >>> round(so_v.T.val, 1)
    33.1
    >>> so_v.set_attr(m=8)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(v.pr.val, 3)
    0.074
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'valve'

    def attr(self):
        return {'pr': dc_cp(min_val=1e-4, max_val=1),
                'zeta': dc_cp(min_val=1e-4),
                'Sirr': dc_simple()}

    def inlets(self):
        return ['in1']

    def outlets(self):
        return ['out1']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.fl_deriv = self.fluid_deriv()
        self.m_deriv = self.mass_flow_deriv()
        self.h_deriv = self.enthalpy_deriv()

    def equations(self):
        r"""
        Calculates vector vec_res with results of equations for this component.

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # eqations for fluids
        vec_res += self.fluid_func()

        ######################################################################
        # eqations for mass flow
        vec_res += self.mass_flow_func()

        ######################################################################
        # eqations for enthalpy
        vec_res += [self.inl[0].h.val_SI - self.outl[0].h.val_SI]

        ######################################################################
        # eqations for specified pressure ratio
        if self.pr.is_set:
            vec_res += [self.inl[0].p.val_SI * self.pr.val -
                        self.outl[0].p.val_SI]

        ######################################################################
        # eqations specified zeta
        if self.zeta.is_set:
            vec_res += [self.zeta_func()]

        return vec_res

    def derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives fluid composition
        mat_deriv += self.fl_deriv

        ######################################################################
        # derivatives for mass flow balance
        mat_deriv += self.m_deriv

        ######################################################################
        # derivatives for enthalpy
        mat_deriv += self.h_deriv

        ######################################################################
        # derivatives for specified pressure ratio
        if self.pr.is_set:
            deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            deriv[0, 0, 1] = self.pr.val
            deriv[0, 1, 1] = -1
            if self.pr.is_var:
                deriv[0, 2 + self.pr.var_pos, 0] = self.inl[0].p.val_SI
            mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for specified zeta
        if self.zeta.is_set:
            deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            deriv[0, 0, 0] = self.numeric_deriv(self.zeta_func, 'm', 0)
            deriv[0, 0, 1] = self.numeric_deriv(self.zeta_func, 'p', 0)
            deriv[0, 0, 2] = self.numeric_deriv(self.zeta_func, 'h', 0)
            deriv[0, 1, 1] = self.numeric_deriv(self.zeta_func, 'p', 1)
            deriv[0, 1, 2] = self.numeric_deriv(self.zeta_func, 'h', 1)
            if self.zeta.is_var:
                deriv[0, 2 + self.zeta.var_pos, 0] = (
                        self.numeric_deriv(self.zeta_func, 'zeta', 2))
            mat_deriv += deriv.tolist()

        return np.asarray(mat_deriv)

    def enthalpy_deriv(self):
        r"""
        Calculates matrix of partial derivatives for enthalpy balance equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
        deriv[0, 0, 2] = 1
        deriv[0, 1, 2] = -1
        return deriv.tolist()

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
                4 \cdot 10^5 & \text{key = 'p'}\\
                5 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 4e5
        elif key == 'h':
            return 5e5

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
                5 \cdot 10^5 & \text{key = 'p'}\\
                5 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            return 5e5

    def calc_parameters(self):
        r"""
        Postprocessing parameter calculation.
        """
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()
        self.pr.val = o[1] / i[1]
        self.zeta.val = ((i[1] - o[1]) * math.pi ** 2 /
                         (8 * i[0] ** 2 * (v_mix_ph(i) + v_mix_ph(o)) / 2))
        self.Sirr.val = i[0] * (s_mix_ph(o) - s_mix_ph(i))

        self.check_parameter_bounds()

# %%


class heat_exchanger_simple(component):
    r"""
    The component heat_exchanger_simple is the parent class for pipe and
    solar_collector.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        **optional equations**

        - :func:`tespy.components.components.heat_exchanger_simple.Q_func`

        .. math::

            0 = p_{in} \cdot pr - p_{out}

        - :func:`tespy.components.components.component.zeta_func`

        - :func:`tespy.components.components.heat_exchanger_simple.darcy_func`
          or :func:`tespy.components.components.heat_exchanger_simple.hw_func`

        **additional equations**

        - :func:`tespy.components.components.heat_exchanger_simple.additional_equations`

    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/pipe.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Q : str/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr : str/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`.

    zeta : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    D : str/float/tespy.helpers.dc_cp
        Diameter of the pipes, :math:`D/\text{m}`.

    L : str/float/tespy.helpers.dc_cp
        Length of the pipes, :math:`L/\text{m}`.

    ks : str/float/tespy.helpers.dc_cp
        Pipes roughness, :math:`ks/\text{m}` for darcy friction,
        :math:`ks/\text{1}` for hazen-williams equation.

    hydro_group : str/tespy.helpers.dc_gcp
        Parametergroup for pressure drop calculation based on pipes dimensions.
        Choose 'HW' for hazen-williams equation, else darcy friction factor is
        used.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient, provide x and y
        values or use generic values (e. g. calculated from design case).
        Standard method 'HE_COLD', Parameter 'm'.

    Tamb : float/tespy.helpers.dc_cp
        Ambient temperature, provide parameter in network's temperature
        unit.

    Tamb_ref : float/tespy.helpers.dc_cp
         Ambient temperature for reference in offdesign case, provide
         parameter in network's temperature unit.

    kA_group : tespy.helpers.dc_gcp
        Parametergroup for heat transfer calculation from ambient temperature
        and area independent heat transfer coefficient kA.

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> import shutil
    >>> fluids = ['H2O']
    >>> nw = nwk.network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> so1 = cmp.source('source 1')
    >>> si1 = cmp.sink('sink 1')
    >>> pi = cmp.pipe('pipe')
    >>> pi.component()
    'pipe'
    >>> pi.set_attr(Tamb=10, pr=0.95, design=['pr'], offdesign=['zeta', 'kA'])
    >>> inc = con.connection(so1, 'out1', pi, 'in1')
    >>> outg = con.connection(pi, 'out1', si1, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> inc.set_attr(fluid={'H2O': 1}, m=1, T=200, p=12)
    >>> outg.set_attr(T=190, design=['T'])
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(pi.Q.val, 1)
    -22252.3
    >>> inc.set_attr(m=1.2)
    >>> pi.set_attr(Tamb=-10)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(pi.kA.val, 1)
    126.5
    >>> round(pi.Q.val, 1)
    -25890.6
    >>> round(outg.T.val, 1)
    189.5
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'heat exchanger simple'

    def attr(self):
        return {'Q': dc_cp(),
                'pr': dc_cp(min_val=1e-4, max_val=1),
                'zeta': dc_cp(min_val=1e-4),
                'D': dc_cp(min_val=1e-2, max_val=2, d=1e-3),
                'L': dc_cp(min_val=1e-1, d=1e-3),
                'ks': dc_cp(val=1e-4, min_val=1e-7, max_val=1e-4, d=1e-8),
                'kA': dc_cp(min_val=0, d=1),
                'Tamb': dc_cp(),
                'kA_char': dc_cc(method='HE_HOT', param='m'),
                'SQ1': dc_simple(), 'SQ2': dc_simple(), 'Sirr': dc_simple(),
                'hydro_group': dc_gcp(), 'kA_group': dc_gcp()}

    def inlets(self):
        return ['in1']

    def outlets(self):
        return ['out1']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.fl_deriv = self.fluid_deriv()
        self.m_deriv = self.mass_flow_deriv()

        self.Tamb.val_SI = ((self.Tamb.val + nw.T[nw.T_unit][0]) *
                            nw.T[nw.T_unit][1])
        self.Tamb.design = ((self.Tamb.design + nw.T[nw.T_unit][0]) *
                            nw.T[nw.T_unit][1])

        # parameters for hydro group
        self.hydro_group.set_attr(elements=[self.L, self.ks, self.D])

        is_set = True
        for e in self.hydro_group.elements:
            if not e.is_set:
                is_set = False

        if is_set:
            self.hydro_group.set_attr(is_set=True)
            if self.hydro_group.method == 'HW':
                method = 'Hazen-Williams equation.'
            else:
                method = 'darcy friction factor.'
            msg = ('Pressure loss calculation from pipe dimensions method is '
                   'set to ' + method + '.')
            logging.debug(msg)

        elif self.hydro_group.is_set:
            msg = ('All parameters of the component group have to be '
                   'specified! This component group uses the following '
                   'parameters: L, ks, D at ' + self.label + '. '
                   'Group will be set to False.')
            logging.info(msg)
            self.hydro_group.set_attr(is_set=False)
        else:
            self.hydro_group.set_attr(is_set=False)

        # parameters for kA group
        self.kA_group.set_attr(elements=[self.kA, self.Tamb])

        is_set = True
        for e in self.kA_group.elements:
            if not e.is_set:
                is_set = False

        if is_set:
            self.kA_group.set_attr(is_set=True)
        elif self.kA_group.is_set:
            msg = ('All parameters of the component group have to be '
                   'specified! This component group uses the following '
                   'parameters: kA, Tamb at ' + self.label + '. '
                   'Group will be set to False.')
            logging.info(msg)
            self.kA_group.set_attr(is_set=False)
        else:
            self.kA_group.set_attr(is_set=False)

    def equations(self):
        r"""
        Calculates vector vec_res with results of equations for this component.

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equations for fluid balance
        vec_res += self.fluid_func()

        ######################################################################
        # equations for mass flow balance
        vec_res += self.mass_flow_func()

        ######################################################################
        # equations for specified heta transfer
        if self.Q.is_set:
            vec_res += [self.Q_func()]

        ######################################################################
        # equations for specified pressure ratio
        if self.pr.is_set:
            vec_res += [self.inl[0].p.val_SI * self.pr.val -
                        self.outl[0].p.val_SI]

        ######################################################################
        # equations for specified zeta
        if self.zeta.is_set:
            vec_res += [self.zeta_func()]

        ######################################################################
        # equation for specified hydro-group paremeters
        if self.hydro_group.is_set:
            # hazen williams equation
            if self.hydro_group.method == 'HW':
                func = self.hw_func
            # darcy friction factor
            else:
                func = self.darcy_func
            vec_res += [func()]

        ######################################################################
        # additional equations
        vec_res += self.additional_equations()

        return vec_res

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for this
        component.

        Equations

            **optional equations**

            - :func:`tespy.components.components.heat_exchanger_simple.kA_func`

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equation for specified kA-group paremeters
        if self.kA_group.is_set:
            vec_res += [self.kA_func()]

        return vec_res

    def derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for fluid balance
        mat_deriv += self.fl_deriv

        ######################################################################
        # derivatives for mass flow balance
        mat_deriv += self.m_deriv

        ######################################################################
        # derivatives for specified heat transfer
        if self.Q.is_set:
            mat_deriv += self.Q_deriv()

        ######################################################################
        # derivatives for specified pressure ratio
        if self.pr.is_set:
            pr_deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            pr_deriv[0, 0, 1] = self.pr.val
            pr_deriv[0, 1, 1] = -1
            # custom variable pr
            if self.pr.is_var:
                pr_deriv[0, 2 + self.pr.var_pos, 0] = self.inl[0].p.val_SI
            mat_deriv += pr_deriv.tolist()

        ######################################################################
        # derivatives for specified zeta
        if self.zeta.is_set:
            zeta_deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            zeta_deriv[0, 0, 0] = self.numeric_deriv(self.zeta_func, 'm', 0)
            zeta_deriv[0, 0, 1] = self.numeric_deriv(self.zeta_func, 'p', 0)
            zeta_deriv[0, 0, 2] = self.numeric_deriv(self.zeta_func, 'h', 0)
            zeta_deriv[0, 1, 1] = self.numeric_deriv(self.zeta_func, 'p', 1)
            zeta_deriv[0, 1, 2] = self.numeric_deriv(self.zeta_func, 'h', 1)
            # custom variable zeta
            if self.zeta.is_var:
                zeta_deriv[0, 2 + self.zeta.var_pos, 0] = (
                        self.numeric_deriv(self.zeta_func, 'zeta', 2))
            mat_deriv += zeta_deriv.tolist()

        ######################################################################
        # derivatives for specified hydro-group parameters
        if self.hydro_group.is_set:
            # hazen williams equation
            if self.hydro_group.method == 'HW':
                func = self.hw_func
            # darcy friction factor
            else:
                func = self.darcy_func

            deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            deriv[0, 0, 0] = self.numeric_deriv(func, 'm', 0)
            deriv[0, 0, 1] = self.numeric_deriv(func, 'p', 0)
            deriv[0, 0, 2] = self.numeric_deriv(func, 'h', 0)
            deriv[0, 1, 1] = self.numeric_deriv(func, 'p', 1)
            deriv[0, 1, 2] = self.numeric_deriv(func, 'h', 1)
            # custom variables of hydro group
            for var in self.hydro_group.elements:
                if var.is_var:
                    deriv[0, 2 + var.var_pos, 0] = (
                            self.numeric_deriv(func, self.vars[var], 2))
            mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for additional equations
        mat_deriv += self.additional_derivatives()

        return np.asarray(mat_deriv)

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional
        equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for specified kA-group paremeters
        if self.kA_group.is_set:
            deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            deriv[0, 0, 0] = self.numeric_deriv(self.kA_func, 'm', 0)
            deriv[0, 0, 1] = self.numeric_deriv(self.kA_func, 'p', 0)
            deriv[0, 0, 2] = self.numeric_deriv(self.kA_func, 'h', 0)
            deriv[0, 1, 1] = self.numeric_deriv(self.kA_func, 'p', 1)
            deriv[0, 1, 2] = self.numeric_deriv(self.kA_func, 'h', 1)
            #
            for var in self.kA_group.elements:
                if var.is_var:
                    deriv[0, 2 + var.var_pos, 0] = (
                            self.numeric_deriv(self.kA_func, self.vars[var], 2)
                            )
            mat_deriv += deriv.tolist()

        return mat_deriv

    def Q_func(self):
        r"""
        Equation for heat transfer of the simple heat exchanger.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) -
                \dot{Q}
        """
        return self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI) - self.Q.val

    def Q_deriv(self):
        r"""
        Calculates the matrix of partial derivatives for heat transfer equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
        deriv[0, 0, 0] = self.outl[0].h.val_SI - self.inl[0].h.val_SI
        deriv[0, 0, 2] = -self.inl[0].m.val_SI
        deriv[0, 1, 2] = self.inl[0].m.val_SI
        # custom variable Q
        if self.Q.is_var:
            deriv[0, 2 + self.Q.var_pos, 0] = -1

        return deriv.tolist()

    def darcy_func(self):
        r"""
        Equation for pressure drop calculation from darcy friction factor.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                Re = \frac{4 \cdot |\dot{m}_{in}|}{\pi \cdot D \cdot
                \frac{\eta_{in}+\eta_{out}}{2}}\\

                0 = p_{in} - p_{out} - \frac{8 \cdot |\dot{m}_{in}| \cdot
                \dot{m}_{in} \cdot \frac{v_{in}+v_{out}}{2} \cdot L \cdot
                \lambda\left(Re, ks, D\right)}{\pi^2 \cdot D^5}\\

                \eta: \text{dynamic viscosity}\\
                v: \text{specific volume}\\
                \lambda: \text{darcy friction factor}
        """
        i, o = self.inl[0].to_flow(), self.outl[0].to_flow()

        if abs(i[0]) < 1e-4:
            return i[1] - o[1]

        visc_i = visc_mix_ph(i, T0=self.inl[0].T.val_SI)
        visc_o = visc_mix_ph(o, T0=self.outl[0].T.val_SI)
        v_i = v_mix_ph(i, T0=self.inl[0].T.val_SI)
        v_o = v_mix_ph(o, T0=self.outl[0].T.val_SI)

        re = 4 * abs(i[0]) / (math.pi * self.D.val * (visc_i + visc_o) / 2)

        return ((i[1] - o[1]) - 8 * abs(i[0]) * i[0] * (v_i + v_o) / 2 *
                self.L.val * lamb(re, self.ks.val, self.D.val) /
                (math.pi ** 2 * self.D.val ** 5))

    def hw_func(self):
        r"""
        Equation for pressure drop calculation from Hazen-Williams equation.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = \left(p_{in} - p_{out} \right) \cdot \left(-1\right)^i -
                \frac{10.67 \cdot |\dot{m}_{in}| ^ {1.852}
                \cdot L}{ks^{1.852} \cdot D^{4.871}} \cdot g \cdot
                \left(\frac{v_{in} + v_{out}}{2}\right)^{0.852}

                i = \begin{cases}
                0 & \dot{m}_{in} \geq 0\\
                1 & \dot{m}_{in} < 0
                \end{cases}

        Note
        ----
        Gravity g is set to :math:`9.81 \frac{m}{s^2}`
        """
        i, o = self.inl[0].to_flow(), self.outl[0].to_flow()

        if abs(i[0]) < 1e-4:
            return i[1] - o[1]

        v_i = v_mix_ph(i, T0=self.inl[0].T.val_SI)
        v_o = v_mix_ph(o, T0=self.outl[0].T.val_SI)
        flow_dir = np.sign(i[0])

        return ((i[1] - o[1]) * flow_dir -
                (10.67 * abs(i[0]) ** 1.852 * self.L.val /
                 (self.ks.val ** 1.852 * self.D.val ** 4.871)) *
                (9.81 * ((v_i + v_o) / 2) ** 0.852))

    def kA_func(self):
        r"""
        Equation for heat transfer calculation from ambient conditions and heat
        transfer coefficient.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                ttd_u = \begin{cases}
                T_{in} - T_{amb} & \dot{m} \geq 0\\
                T_{out} - T_{amb} & \dot{m} < 0
                \end{cases}

                ttd_l = \begin{cases}
                T_{in} - T_{amb} & \dot{m} < 0\\
                T_{out} - T_{amb} & \dot{m} \geq 0
                \end{cases}

                0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in}\right) +
                kA \cdot f_{kA} \cdot \frac{ttd_u - ttd_l}
                {\ln{\frac{ttd_u}{ttd_l}}}

                f_{kA} = f_1\left(\frac{m_1}{m_{1,ref}}\right)

                T_{amb}: \text{ambient temperature}

            for f\ :subscript:`1` \ see class
            :func:`tespy.components.characteristics.characteristics`
        """
        i, o = self.inl[0].to_flow(), self.outl[0].to_flow()

        ttd_1 = T_mix_ph(i, T0=self.inl[0].T.val_SI) - self.Tamb.val_SI
        ttd_2 = T_mix_ph(o, T0=self.outl[0].T.val_SI) - self.Tamb.val_SI

        if ttd_1 > ttd_2:
            td_log = (ttd_1 - ttd_2) / math.log(ttd_1 / ttd_2)
        elif ttd_1 < ttd_2:
            td_log = (ttd_2 - ttd_1) / math.log(ttd_2 / ttd_1)
        else:
            td_log = 0

        fkA = 1
        if not np.isnan(self.inl[0].m.design):
            if self.kA_char.param == 'm':
                fkA = self.kA_char.func.f_x(i[0] / self.inl[0].m.design)

        return i[0] * (o[2] - i[2]) + self.kA.val * fkA * td_log

    def bus_func(self, bus):
        r"""
        Calculates the residual value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        val : float
            Residual value of equation.

            .. math::

                val = P \cdot f\left( \frac{P}{P_{ref}}\right)

                P = \dot{m}_{in} \cdot \left( h_{out} - h_{in} \right)
        """
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()

        val = i[0] * (o[2] - i[2])
        if np.isnan(bus.P_ref):
            expr = 1
        else:
            expr = abs(val / bus.P_ref)
        return val * bus.char.f_x(expr)

    def bus_deriv(self, bus):
        r"""
        Calculates the matrix of partial derivatives of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 2, len(self.inl[0].fluid.val) + 3))
        deriv[0, 0, 0] = self.numeric_deriv(self.bus_func, 'm', 0, bus=bus)
        deriv[0, 0, 2] = self.numeric_deriv(self.bus_func, 'h', 0, bus=bus)
        deriv[0, 1, 2] = self.numeric_deriv(self.bus_func, 'h', 1, bus=bus)
        return deriv

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
                \begin{cases}
                1 \cdot 10^5 \; \frac{\text{J}}{\text{kg}} & \dot{Q} < 0\\
                3 \cdot 10^5 \; \frac{\text{J}}{\text{kg}} & \dot{Q} = 0\\
                5 \cdot 10^5 \; \frac{\text{J}}{\text{kg}} & \dot{Q} > 0
                \end{cases} & \text{key = 'h'}\\
                \; \; \; \; 10^5 \text{Pa} & \text{key = 'p'}
                \end{cases}

        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            if self.Q.val < 0 and self.Q.is_set:
                return 1e5
            elif self.Q.val > 0 and self.Q.is_set:
                return 5e5
            else:
                return 3e5

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
                1 \cdot 10^5 & \text{key = 'p'}\\
                \begin{cases}
                5 \cdot 10^5 & \dot{Q} < 0\\
                3 \cdot 10^5 & \dot{Q} = 0\\
                1 \cdot 10^5 & \dot{Q} > 0
                \end{cases} & \text{key = 'h'}\\
                \end{cases}
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            if self.Q.val < 0 and self.Q.is_set:
                return 5e5
            elif self.Q.val > 0 and self.Q.is_set:
                return 1e5
            else:
                return 3e5

    def calc_parameters(self):
        r"""
        Postprocessing parameter calculation.
        """
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()
        v_i = v_mix_ph(i, T0=self.inl[0].T.val_SI)
        v_o = v_mix_ph(o, T0=self.outl[0].T.val_SI)

        self.SQ1.val = i[0] * (s_mix_ph(o) - s_mix_ph(i))
        self.Q.val = i[0] * (o[2] - i[2])
        self.pr.val = o[1] / i[1]
        self.zeta.val = ((i[1] - o[1]) * math.pi ** 2 /
                         (8 * i[0] ** 2 * (v_i + v_o) / 2))

        if self.Tamb.is_set:
            self.SQ2.val = -i[0] * (o[2] - i[2]) / self.Tamb.val_SI
            self.Sirr.val = self.SQ1.val + self.SQ2.val

            ttd_1 = T_mix_ph(i, T0=self.inl[0].T.val_SI) - self.Tamb.val_SI
            ttd_2 = T_mix_ph(o, T0=self.outl[0].T.val_SI) - self.Tamb.val_SI

            if ttd_1 > ttd_2:
                td_log = (ttd_1 - ttd_2) / math.log(ttd_1 / ttd_2)
            elif ttd_1 < ttd_2:
                td_log = (ttd_2 - ttd_1) / math.log(ttd_2 / ttd_1)
            else:
                td_log = 0

            self.kA.val = abs(i[0] * (o[2] - i[2]) / td_log)

        if self.kA.is_set:
            # get bound errors for kA characteristic line
            if self.kA_char.param == 'm':
                self.kA_char.func.get_bound_errors(i[0] / self.inl[0].m.design,
                                                   self.label)

        self.check_parameter_bounds()

# %%


class pipe(heat_exchanger_simple):
    r"""
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        **optional equations**

        - :func:`tespy.components.components.heat_exchanger_simple.Q_func`

        .. math::

            0 = p_{in} \cdot pr - p_{out}

        - :func:`tespy.components.components.component.zeta_func`

        - :func:`tespy.components.components.heat_exchanger_simple.darcy_func`
          or :func:`tespy.components.components.heat_exchanger_simple.hw_func`

        **additional equations**

        - :func:`tespy.components.components.heat_exchanger_simple.additional_equations`

    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/pipe.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Q : str/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr : str/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`.

    zeta : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    D : str/float/tespy.helpers.dc_cp
        Diameter of the pipes, :math:`D/\text{m}`.

    L : str/float/tespy.helpers.dc_cp
        Length of the pipes, :math:`L/\text{m}`.

    ks : str/float/tespy.helpers.dc_cp
        Pipes roughness, :math:`ks/\text{m}` for darcy friction,
        :math:`ks/\text{1}` for hazen-williams equation.

    hydro_group : str/tespy.helpers.dc_gcp
        Parametergroup for pressure drop calculation based on pipes dimensions.
        Choose 'HW' for hazen-williams equation, else darcy friction factor is
        used.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient, provide x and y
        values or use generic values (e. g. calculated from design case).
        Standard method 'HE_COLD', Parameter 'm'.

    Tamb : float/tespy.helpers.dc_cp
        Ambient temperature, provide parameter in network's temperature
        unit.

    Tamb_ref : float/tespy.helpers.dc_cp
         Ambient temperature for reference in offdesign case, provide
         parameter in network's temperature unit.

    kA_group : tespy.helpers.dc_gcp
        Parametergroup for heat transfer calculation from ambient temperature
        and area independent heat transfer coefficient kA.

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> import shutil
    >>> fluids = ['H2O']
    >>> nw = nwk.network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> so1 = cmp.source('source 1')
    >>> si1 = cmp.sink('sink 1')
    >>> pi = cmp.pipe('pipe')
    >>> pi.component()
    'pipe'
    >>> pi.set_attr(pr=0.95, Q=0, design=['pr'], L=100, D='var', ks=5e-5)
    >>> inc = con.connection(so1, 'out1', pi, 'in1')
    >>> outg = con.connection(pi, 'out1', si1, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> inc.set_attr(fluid={'H2O': 1}, m=1, T=100, p=12)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(pi.D.val, 3)
    0.032
    >>> round(outg.p.val, 1)
    11.4
    >>> inc.set_attr(m=1.2)
    >>> pi.set_attr(D=pi.D.val)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(outg.p.val, 2)
    11.14
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'pipe'

# %%


class solar_collector(heat_exchanger_simple):
    r"""
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        **optional equations**

        - :func:`tespy.components.components.heat_exchanger_simple.Q_func`

        .. math::

            0 = p_{in} \cdot pr - p_{out}

        - :func:`tespy.components.components.component.zeta_func`

        - :func:`tespy.components.components.heat_exchanger_simple.darcy_func`
          or :func:`tespy.components.components.heat_exchanger_simple.hw_func`

        **additional equations**

        - :func:`tespy.components.components.solar_collector.additional_equations`

    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/solar_collector.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Q : str/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr : str/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`.

    zeta : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    D : str/float/tespy.helpers.dc_cp
        Diameter of the pipes, :math:`D/\text{m}`.

    L : str/float/tespy.helpers.dc_cp
        Length of the pipes, :math:`L/\text{m}`.

    ks : str/float/tespy.helpers.dc_cp
        Pipes roughness, :math:`ks/\text{m}` for darcy friction,
        :math:`ks/\text{1}` for hazen-williams equation.

    hydro_group : str/tespy.helpers.dc_gcp
        Parametergroup for pressure drop calculation based on pipes dimensions.
        Choose 'HW' for hazen-williams equation, else darcy friction factor is
        used.

    E : str/float/tespy.helpers.dc_cp
        Absorption on the inclined surface,
        :math:`E/\frac{\text{W}}{\text{m}^2}`.

    lkf_lin : str/float/tespy.helpers.dc_cp
        Linear loss key figure,
        :math:`\alpha_1/\frac{\text{W}}{\text{K} \cdot \text{m}^2}`.

    lkf_quad : str/float/tespy.helpers.dc_cp
        Quadratic loss key figure,
        :math:`\alpha_2/\frac{\text{W}}{\text{K}^2 \cdot \text{m}^2}`.

    A : str/float/tespy.helpers.dc_cp
        Collector surface area :math:`A/\text{m}^2`.

    Tamb : float/tespy.helpers.dc_cp
        Ambient temperature, provide parameter in network's temperature unit.

    energy_group : tespy.helpers.dc_gcp
        Parametergroup for energy balance of solarthermal collector.

    Note
    ----
    The solar collector does not take optical losses into accout. The incoming
    radiation E represents the actual absorption of the solar collector.
    Optical losses should be handeled in preprocessing.

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> import shutil
    >>> fluids = ['H2O']
    >>> nw = nwk.network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> so1 = cmp.source('source 1')
    >>> si1 = cmp.sink('sink 1')
    >>> sc = cmp.solar_collector('test')
    >>> sc.component()
    'solar collector'
    >>> sc.set_attr(pr=0.95, Q=1e4, design=['pr', 'Q'], offdesign=['zeta'],
    ...     Tamb=25, A='var', lkf_lin=1, lkf_quad=0.005, E=8e2)
    >>> inc = con.connection(so1, 'out1', sc, 'in1')
    >>> outg = con.connection(sc, 'out1', si1, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> inc.set_attr(fluid={'H2O': 1}, T=40, p=3, offdesign=['m'])
    >>> outg.set_attr(T=90, design=['T'])
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(sc.A.val, 1)
    13.3
    >>> sc.set_attr(A=sc.A.val, E=5e2, Tamb=20)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(sc.Q.val, 1)
    6097.3
    >>> round(outg.T.val, 1)
    70.5
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'solar collector'

    def attr(self):
        return {'Q': dc_cp(),
                'pr': dc_cp(min_val=1e-4, max_val=1),
                'zeta': dc_cp(min_val=1e-4),
                'D': dc_cp(min_val=1e-2, max_val=2, d=1e-3),
                'L': dc_cp(min_val=1e-1, d=1e-3),
                'ks': dc_cp(val=1e-4, min_val=1e-7, max_val=1e-4, d=1e-8),
                'E': dc_cp(min_val=0),
                'lkf_lin': dc_cp(min_val=0),
                'lkf_quad': dc_cp(min_val=0),
                'A': dc_cp(min_val=0),
                'Tamb': dc_cp(),
                'SQ': dc_simple(),
                'hydro_group': dc_gcp(), 'energy_group': dc_gcp()}

    def inlets(self):
        return ['in1']

    def outlets(self):
        return ['out1']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.fl_deriv = self.fluid_deriv()
        self.m_deriv = self.mass_flow_deriv()

        self.Tamb.val_SI = ((self.Tamb.val + nw.T[nw.T_unit][0]) *
                            nw.T[nw.T_unit][1])

        # parameters for hydro group
        self.hydro_group.set_attr(elements=[self.L, self.ks, self.D])

        is_set = True
        for e in self.hydro_group.elements:
            if not e.is_set:
                is_set = False

        if is_set:
            self.hydro_group.set_attr(is_set=True)
        elif self.hydro_group.is_set:
            msg = ('All parameters of the component group have to be '
                   'specified! This component group uses the following '
                   'parameters: L, ks, D at ' + self.label + '. '
                   'Group will be set to False.')
            logging.info(msg)
            self.hydro_group.set_attr(is_set=False)
        else:
            self.hydro_group.set_attr(is_set=False)

        # parameters for kA group
        self.energy_group.set_attr(elements=[self.E, self.lkf_lin,
                                             self.lkf_quad, self.A, self.Tamb])

        is_set = True
        for e in self.energy_group.elements:
            if not e.is_set:
                is_set = False

        if is_set:
            self.energy_group.set_attr(is_set=True)
        elif self.energy_group.is_set:
            msg = ('All parameters of the component group have to be '
                   'specified! This component group uses the following '
                   'parameters: E, lkf_lin, lkf_quad, A, Tamb at ' +
                   self.label + '. Group will be set to False.')
            logging.info(msg)
            self.energy_group.set_attr(is_set=False)
        else:
            self.energy_group.set_attr(is_set=False)

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for this
        component.

        Equations

            **optional equations**

            - :func:`tespy.components.components.solar_collector.energy_func`

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equation for specified energy-group paremeters
        if self.energy_group.is_set:
            vec_res += [self.energy_func()]

        return vec_res

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional
        equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for specified energy-group paremeters
        if self.energy_group.is_set:
            deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            deriv[0, 0, 0] = self.outl[0].h.val_SI - self.inl[0].h.val_SI
            deriv[0, 0, 1] = self.numeric_deriv(self.energy_func, 'p', 0)
            deriv[0, 0, 2] = self.numeric_deriv(self.energy_func, 'h', 0)
            deriv[0, 1, 1] = self.numeric_deriv(self.energy_func, 'p', 1)
            deriv[0, 1, 2] = self.numeric_deriv(self.energy_func, 'h', 1)
            # custom variables for the energy-group
            for var in self.energy_group.elements:
                if var.is_var:
                    deriv[0, 2 + var.var_pos, 0] = (
                            self.numeric_deriv(self.energy_func,
                                               self.vars[var], 2))
            mat_deriv += deriv.tolist()

        return mat_deriv

    def energy_func(self):
        r"""
        Equation for solar collector energy balance.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                T_m = \frac{T_{out} + T_{in}}{2}\\

                \begin{split}
                0 = & \dot{m} \cdot \left( h_{out} - h_{in} \right)\\
                & - A \cdot \left[E - \alpha_1 \cdot
                \left(T_m - T_{amb} \right) - \alpha_2 \cdot
                \left(T_m - T_{amb}\right)^2 \right]
                \end{split}
        """

        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()

        T_m = (T_mix_ph(i, T0=self.inl[0].T.val_SI) +
               T_mix_ph(o, T0=self.outl[0].T.val_SI)) / 2

        return (i[0] * (o[2] - i[2]) - self.A.val * (
                self.E.val - (T_m - self.Tamb.val_SI) * self.lkf_lin.val -
                self.lkf_quad.val * (T_m - self.Tamb.val_SI) ** 2))

    def calc_parameters(self):
        r"""
        Postprocessing parameter calculation.
        """
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()

        self.SQ.val = i[0] * (s_mix_ph(o) - s_mix_ph(i))
        self.Q.val = i[0] * (o[2] - i[2])
        self.pr.val = o[1] / i[1]
        self.zeta.val = ((i[1] - o[1]) * math.pi ** 2 /
                         (8 * i[0] ** 2 * (v_mix_ph(i) + v_mix_ph(o)) / 2))

        self.check_parameter_bounds()

# %%


class heat_exchanger(component):
    r"""
    Class heat_exchanger is the parent class for condenser and desuperheater.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.heat_exchanger.fluid_func`
        - :func:`tespy.components.components.heat_exchanger.mass_flow_func`

        **heat exchanger**
        - :func:`tespy.components.components.heat_exchanger.energy_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) - \dot{Q}

        **heat exchanger**

        - :func:`tespy.components.components.heat_exchanger.kA_func`
        - :func:`tespy.components.components.heat_exchanger.ttd_u_func`
        - :func:`tespy.components.components.heat_exchanger.ttd_l_func`

        .. math::

            0 = p_{1,in} \cdot pr1 - p_{1,out}\\
            0 = p_{2,in} \cdot pr2 - p_{2,out}

        - :func:`tespy.components.components.heat_exchanger.zeta_func`
        - :func:`tespy.components.components.heat_exchanger.zeta2_func`

        **additional equations**

        - :func:`tespy.components.components.heat_exchanger.additional_equations`

    Inlets/Outlets

        - in1, in2 (index 1: hot side, index 2: cold side)
        - out1, out2 (index 1: hot side, index 2: cold side)

    Image

        .. image:: _images/heat_exchanger.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Q : str/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : str/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : str/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char1 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at hot side, provide
        x and y values or use generic values (e. g. calculated from design
        case). Standard method 'HE_HOT', Parameter 'm'.

    kA_char2 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at cold side,
        provide x and y values or use generic values (e. g. calculated from
        design case). Standard method 'HE_COLD', Parameter 'm'.

    Note
    ---
    The heat exchanger and subclasses (desuperheater, condenser) are
    countercurrent heat exchangers. Equations (kA, ttd_u, ttd_l) do not work
    for directcurrent and crosscurrent or combinations of different types.

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> import shutil
    >>> nw = nwk.network(fluids=['water'], T_unit='C', p_unit='bar',
    ...     h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> tesin = cmp.sink('TES in')
    >>> tesout = cmp.source('TES out')
    >>> hsin = cmp.sink('HS in')
    >>> hsout = cmp.source('HS out')
    >>> he = cmp.heat_exchanger('heat exchanger')
    >>> he.component()
    'heat exchanger'
    >>> tes_he = con.connection(tesout, 'out1', he, 'in2')
    >>> he_tes = con.connection(he, 'out2', tesin, 'in1')
    >>> hs_he = con.connection(hsout, 'out1', he, 'in1')
    >>> he_hs = con.connection(he, 'out1', hsin, 'in1')
    >>> nw.add_conns(tes_he, he_tes, hs_he, he_hs)
    >>> he.set_attr(pr1=0.98, pr2=0.98, ttd_u=5,
    ...     design=['pr1', 'pr2', 'ttd_u'], offdesign=['zeta1', 'zeta2', 'kA'])
    >>> hs_he.set_attr(Td_bp=-10, p=3, fluid={'water': 1})
    >>> he_hs.set_attr(T=70)
    >>> tes_he.set_attr(p=5, fluid={'water': 1})
    >>> tes_he.set_attr(T=40)
    >>> he.set_attr(Q=-80e3)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(tes_he.m.val, 2)
    0.24
    >>> round(he_tes.T.val, 1)
    118.5
    >>> he.set_attr(Q=-60e3)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(tes_he.m.val, 2)
    0.18
    >>> round(he_tes.T.val, 1)
    119.4
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'heat exchanger'

    def attr(self):
        return {'Q': dc_cp(max_val=0),
                'kA': dc_cp(min_val=0),
                'td_log': dc_cp(min_val=0),
                'ttd_u': dc_cp(min_val=0), 'ttd_l': dc_cp(min_val=0),
                'pr1': dc_cp(max_val=1), 'pr2': dc_cp(max_val=1),
                'zeta1': dc_cp(min_val=0), 'zeta2': dc_cp(min_val=0),
                'kA_char1': dc_cc(method='HE_HOT', param='m'),
                'kA_char2': dc_cc(method='HE_COLD', param='m'),
                'SQ1': dc_simple(), 'SQ2': dc_simple(), 'Sirr': dc_simple(),
                'zero_flag': dc_simple()}

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1', 'out2']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.fl_deriv = self.fluid_deriv()
        self.m_deriv = self.mass_flow_deriv()

    def equations(self):
        r"""
        Calculates vector vec_res with results of equations for this component.

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equations for fluid balance
        vec_res += self.fluid_func()

        ######################################################################
        # equations for mass flow balance
        vec_res += self.mass_flow_func()

        ######################################################################
        # equations for energy balance
        vec_res += [self.energy_func()]

        ######################################################################
        # equations for specified heat transfer
        if self.Q.is_set:
            vec_res += [self.inl[0].m.val_SI *
                        (self.outl[0].h.val_SI - self.inl[0].h.val_SI) -
                        self.Q.val]

        ######################################################################
        # equations for specified heat transfer coefficient
        if self.kA.is_set:
            vec_res += [self.kA_func()]

        ######################################################################
        # equations for specified upper terminal temperature difference
        if self.ttd_u.is_set:
            vec_res += [self.ttd_u_func()]

        ######################################################################
        # equations for specified lower terminal temperature difference
        if self.ttd_l.is_set:
            vec_res += [self.ttd_l_func()]

        ######################################################################
        # equations for specified pressure ratio at hot side
        if self.pr1.is_set:
            vec_res += [self.pr1.val * self.inl[0].p.val_SI -
                        self.outl[0].p.val_SI]

        ######################################################################
        # equations for specified pressure ratio at cold side
        if self.pr2.is_set:
            vec_res += [self.pr2.val * self.inl[1].p.val_SI -
                        self.outl[1].p.val_SI]

        ######################################################################
        # equations for specified zeta at hot side
        if self.zeta1.is_set:
            vec_res += [self.zeta_func()]

        ######################################################################
        # equations for specified zeta at cold side
        if self.zeta2.is_set:
            vec_res += [self.zeta2_func()]

        ######################################################################
        # additional equations
        vec_res += self.additional_equations()

        return vec_res

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for
        this component.

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        return []

    def derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for fluid balance equations
        mat_deriv += self.fl_deriv
        ######################################################################
        # derivatives for mass flow balance equations
        mat_deriv += self.m_deriv

        ######################################################################
        # derivatives for energy balance equation
        mat_deriv += self.energy_deriv()

        ######################################################################
        # derivatives for specified heat transfer
        if self.Q.is_set:
            deriv = np.zeros((1, 4, self.num_fl + 3))
            deriv[0, 0, 0] = self.outl[0].h.val_SI - self.inl[0].h.val_SI
            deriv[0, 0, 2] = -self.inl[0].m.val_SI
            deriv[0, 2, 2] = self.inl[0].m.val_SI
            mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for specified heat transfer coefficient
        if self.kA.is_set:
            kA_deriv = np.zeros((1, 4, self.num_fl + 3))
            kA_deriv[0, 0, 0] = self.numeric_deriv(self.kA_func, 'm', 0)
            kA_deriv[0, 1, 0] = self.numeric_deriv(self.kA_func, 'm', 1)
            for i in range(4):
                kA_deriv[0, i, 1] = self.numeric_deriv(self.kA_func, 'p', i)
                kA_deriv[0, i, 2] = self.numeric_deriv(self.kA_func, 'h', i)
            mat_deriv += kA_deriv.tolist()

        ######################################################################
        # derivatives for specified upper terminal temperature difference
        if self.ttd_u.is_set:
            mat_deriv += self.ttd_u_deriv()

        ######################################################################
        # derivatives for specified lower terminal temperature difference
        if self.ttd_l.is_set:
            mat_deriv += self.ttd_l_deriv()

        ######################################################################
        # derivatives for specified pressure ratio at hot side
        if self.pr1.is_set:
            pr1_deriv = np.zeros((1, 4, self.num_fl + 3))
            pr1_deriv[0, 0, 1] = self.pr1.val
            pr1_deriv[0, 2, 1] = -1
            mat_deriv += pr1_deriv.tolist()

        ######################################################################
        # derivatives for specified pressure ratio at cold side
        if self.pr2.is_set:
            pr2_deriv = np.zeros((1, 4, self.num_fl + 3))
            pr2_deriv[0, 1, 1] = self.pr2.val
            pr2_deriv[0, 3, 1] = -1
            mat_deriv += pr2_deriv.tolist()

        ######################################################################
        # derivatives for specified zeta at hot side
        if self.zeta1.is_set:
            zeta1_deriv = np.zeros((1, 4, self.num_fl + 3))
            zeta1_deriv[0, 0, 0] = self.numeric_deriv(self.zeta_func, 'm', 0)
            zeta1_deriv[0, 0, 1] = self.numeric_deriv(self.zeta_func, 'p', 0)
            zeta1_deriv[0, 0, 2] = self.numeric_deriv(self.zeta_func, 'h', 0)
            zeta1_deriv[0, 2, 1] = self.numeric_deriv(self.zeta_func, 'p', 2)
            zeta1_deriv[0, 2, 2] = self.numeric_deriv(self.zeta_func, 'h', 2)
            mat_deriv += zeta1_deriv.tolist()

        ######################################################################
        # derivatives for specified zeta at cold side
        if self.zeta2.is_set:
            zeta2_deriv = np.zeros((1, 4, self.num_fl + 3))
            zeta2_deriv[0, 1, 0] = self.numeric_deriv(self.zeta2_func, 'm', 1)
            zeta2_deriv[0, 1, 1] = self.numeric_deriv(self.zeta2_func, 'p', 1)
            zeta2_deriv[0, 1, 2] = self.numeric_deriv(self.zeta2_func, 'h', 1)
            zeta2_deriv[0, 3, 1] = self.numeric_deriv(self.zeta2_func, 'p', 3)
            zeta2_deriv[0, 3, 2] = self.numeric_deriv(self.zeta2_func, 'h', 3)
            mat_deriv += zeta2_deriv.tolist()

        ######################################################################
        # derivatives for additional equations
        mat_deriv += self.additional_derivatives()

        return np.asarray(mat_deriv)

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional
        equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        return []

    def fluid_func(self):
        r"""
        Calculates the vector of residual values for component's fluid balance
        equations.

        Returns
        -------
        vec_res : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = fluid_{i,in_{j}} - fluid_{i,out_{j}} \;
                \forall i \in \mathrm{fluid}, \; \forall j \in inlets/outlets
        """
        vec_res = []

        for i in range(self.num_i):
            for fluid, x in self.inl[i].fluid.val.items():
                vec_res += [x - self.outl[i].fluid.val[fluid]]
        return vec_res

    def mass_flow_func(self):
        r"""
        Calculates the residual value for component's mass flow balance
        equation.

        Returns
        -------
        vec_res : list
            Vector with residual value for component's mass flow balance.

            .. math::

                0 = \dot{m}_{in,i} - \dot{m}_{out,i} \;
                \forall i \in inlets/outlets
        """
        vec_res = []
        for i in range(self.num_i):
            vec_res += [self.inl[i].m.val_SI - self.outl[i].m.val_SI]
        return vec_res

    def fluid_deriv(self):
        r"""
        Calculates the partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((self.num_fl * 2, 4 + self.num_vars, 3 + self.num_fl))
        # hot side
        i = 0
        for fluid in self.fluids:
            deriv[i, 0, i + 3] = 1
            deriv[i, 2, i + 3] = -1
            i += 1
        # cold side
        j = 0
        for fluid in self.fluids:
            deriv[i + j, 1, j + 3] = 1
            deriv[i + j, 3, j + 3] = -1
            j += 1
        return deriv.tolist()

    def mass_flow_deriv(self):
        r"""
        Calculates the partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the mass flow balance
            equations.
        """
        deriv = np.zeros((2, 4 + self.num_vars, self.num_fl + 3))
        for i in range(self.num_i):
            deriv[i, i, 0] = 1
        for j in range(self.num_o):
            deriv[j, j + i + 1, 0] = -1
        return deriv.tolist()

    def energy_func(self):
        r"""
        Equation for heat exchanger energy balance.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{1,in} \cdot \left(h_{1,out} - h_{1,in} \right) +
                \dot{m}_{2,in} \cdot \left(h_{2,out} - h_{2,in} \right)
        """
#        if self.zero_flag.val_set:
#            c = self.zero_flag.val
#            if c[0] > 0 and c[1] < 3:
#                return self.inl[0].m.val_SI
#
#            elif ((c[0] == 0 and c[1] < 3) or
#                  (c[0] > 1 and c[1] > 2 and c[1] < 5)):
#                return self.outl[0].h.val_SI - self.inl[0].h.val_SI
#
#            elif ((c[0] < 2 and c[1] > 2 and c[1] < 5) or
#                  (c[0] == 3 and c[1] == 5)):
#                return self.inl[1].m.val_SI
#            else:
#                return self.outl[1].h.val_SI - self.inl[1].h.val_SI
#
#        else:
        return (self.inl[0].m.val_SI * (self.outl[0].h.val_SI -
                                        self.inl[0].h.val_SI) +
                self.inl[1].m.val_SI * (self.outl[1].h.val_SI -
                                        self.inl[1].h.val_SI))

    def energy_deriv(self):
        r"""
        Calculates the matrix of partial derivatives for energy balance
        equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 4, len(self.inl[0].fluid.val) + 3))

#        if self.zero_flag.val_set:
#            c = self.zero_flag.val
#            if c[0] > 0 and c[1] < 3:
#                deriv[0, 0, 0] = 1
#
#            elif ((c[0] == 0 and c[1] < 3) or
#                  (c[0] > 1 and c[1] > 2 and c[1] < 5)):
#                deriv[0, 0, 2] = -1
#                deriv[0, 2, 2] = 1
#
#            elif ((c[0] < 2 and c[1] > 2 and c[1] < 5) or
#                  (c[0] == 3 and c[1] == 5)):
#                deriv[0, 1, 0] = 1
#            else:
#                deriv[0, 1, 2] = -1
#                deriv[0, 3, 2] = 1
#
#        else:
        for k in range(2):
            deriv[0, k, 0] = self.outl[k].h.val_SI - self.inl[k].h.val_SI
            deriv[0, k, 2] = -self.inl[k].m.val_SI

        deriv[0, 2, 2] = self.inl[0].m.val_SI
        deriv[0, 3, 2] = self.inl[1].m.val_SI
        return deriv.tolist()

    def kA_func(self):
        r"""
        Equation for heat transfer from conditions on both sides of heat
        exchanger.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = \dot{m}_{1,in} \cdot \left( h_{1,out} - h_{1,in}\right) +
                kA \cdot f_{kA} \cdot \frac{T_{1,out} -
                T_{2,in} - T_{1,in} + T_{2,out}}
                {\ln{\frac{T_{1,out} - T_{2,in}}{T_{1,in} - T_{2,out}}}}

                f_{kA} = f_1\left(\frac{m_1}{m_{1,ref}}\right) \cdot
                f_2\left(\frac{m_2}{m_{2,ref}}\right)

        Note
        ----
        For standard functions f\ :subscript:`1` \ and f\ :subscript:`2` \ see
        class :func:`tespy.components.characteristics.characteristics`.

        - Calculate temperatures at inlets and outlets.
        - Perform value manipulation, if temperature levels are not physically
          feasible.
        """

#        if self.zero_flag.val_set:
#            c = self.zero_flag.val
#            if c[1] == 2 or c[1] == 4 or c[1] == 5:
#                T_i1 = T_mix_ph(self.inl[0].to_flow(), T0=self.inl[0].T.val_SI)
#                T_i2 = T_mix_ph(self.inl[1].to_flow(), T0=self.inl[1].T.val_SI)
#                T_o1 = T_mix_ph(self.outl[0].to_flow(),
#                                T0=self.outl[0].T.val_SI)
#                T_o2 = T_mix_ph(self.outl[1].to_flow(),
#                                T0=self.outl[1].T.val_SI)
#                return T_o1 - T_i2 - T_i1 + T_o2
#
#            elif c[0] < 3 and (c[1] == 1 or c[1] == 3):
#                return self.outl[1].h.val_SI - self.inl[1].h.val_SI
#
#            elif ((c[0] < 2 and c[1] == 0) or
#                  (c[0] == 3 and (c[1] == 1 or c[1] == 3))):
#                return self.inl[1].m.val_SI
#
#            else:
#                return self.outl[0].h.val_SI - self.inl[0].h.val_SI

        i1 = self.inl[0].to_flow()
        i2 = self.inl[1].to_flow()
        o1 = self.outl[0].to_flow()
        o2 = self.outl[1].to_flow()

        i1_d = self.inl[0].to_flow_design()
        i2_d = self.inl[1].to_flow_design()

        T_i1 = T_mix_ph(i1, T0=self.inl[0].T.val_SI)
        T_i2 = T_mix_ph(i2, T0=self.inl[1].T.val_SI)
        T_o1 = T_mix_ph(o1, T0=self.outl[0].T.val_SI)
        T_o2 = T_mix_ph(o2, T0=self.outl[1].T.val_SI)

#        if T_i1 <= T_o2 and self.inl[0].T.val_set is False:
        if T_i1 <= T_o2:
            T_i1 = T_o2 + 0.01
#        if T_i1 <= T_o2 and self.outl[1].T.val_set is False:
        if T_i1 <= T_o2:
            T_o2 = T_i1 - 0.01
#        if T_i1 < T_o2 and self.inl[0].T.val_set and self.outl[1].T.val_set:
#            msg = ('Infeasibility at ' + str(self.label) + ': Value for upper '
#                   'temperature difference is ' + str(round(T_i1 - T_o2)) +
#                   '.')
#            logging.error(msg)
#            raise ValueError(msg)

#        if T_i1 <= T_o2 and self.outl[1].T.val_set is False:
        if T_i1 <= T_o2:
            T_o1 = T_i2 + 0.02
#        if T_o1 <= T_i2 and self.inl[1].T.val_set is False:
        if T_o1 <= T_i2:
            T_i2 = T_o1 - 0.02
#        if T_o1 < T_i2 and self.inl[1].T.val_set and self.outl[0].T.val_set:
#            msg = ('Infeasibility at ' + str(self.label) + ': Value for lower '
#                   'temperature difference is ' + str(round(T_o1 - T_i2)) +
#                   '.')
#            logging.error(msg)
#            raise ValueError(msg)

        fkA1 = 1
        if self.kA_char1.param == 'm':
            if not np.isnan(i1_d[0]):
                if not i1[0] == 0:
                    fkA1 = self.kA_char1.func.f_x(i1[0] / i1_d[0])

        fkA2 = 1
        if self.kA_char2.param == 'm':
            if not np.isnan(i2_d[0]):
                if not i2[0] == 0:
                    fkA2 = self.kA_char2.func.f_x(i2[0] / i2_d[0])

        td_log = ((T_o1 - T_i2 - T_i1 + T_o2) /
                  math.log((T_o1 - T_i2) / (T_i1 - T_o2)))
        return i1[0] * (o1[2] - i1[2]) + self.kA.val * fkA1 * fkA2 * td_log

    def ttd_u_func(self):
        r"""
        Equation for upper terminal temperature difference.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = ttd_{u} - T_{1,in} + T_{2,out}
        """
        T_i1 = T_mix_ph(self.inl[0].to_flow(), T0=self.inl[0].T.val_SI)
        T_o2 = T_mix_ph(self.outl[1].to_flow(), T0=self.outl[1].T.val_SI)
        return self.ttd_u.val - T_i1 + T_o2

    def ttd_u_deriv(self):
        r"""
        Calculates the matrix of partial derivatives for upper temperature
        difference equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 4, len(self.inl[0].fluid.val) + 3))
        for i in range(2):
            deriv[0, i * 3, 1] = self.numeric_deriv(self.ttd_u_func,
                                                    'p', i * 3)
            deriv[0, i * 3, 2] = self.numeric_deriv(self.ttd_u_func,
                                                    'h', i * 3)
        return deriv.tolist()

    def ttd_l_func(self):
        r"""
        Equation for upper terminal temperature difference.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = ttd_{l} - T_{1,out} + T_{2,in}
        """
        i2 = self.inl[1].to_flow()
        o1 = self.outl[0].to_flow()
        return (self.ttd_l.val - T_mix_ph(o1, T0=self.outl[0].T.val_SI) +
                T_mix_ph(i2, T0=self.inl[1].T.val_SI))

    def ttd_l_deriv(self):
        r"""
        Calculates the matrix of partial derivatives for lower temperature
        difference equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 4, len(self.inl[0].fluid.val) + 3))
        for i in range(2):
            deriv[0, i + 1, 1] = self.numeric_deriv(self.ttd_l_func,
                                                    'p', i + 1)
            deriv[0, i + 1, 2] = self.numeric_deriv(self.ttd_l_func,
                                                    'h', i + 1)
        return deriv.tolist()

    def bus_func(self, bus):
        r"""
        Calculates the residual value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        val : float
            Residual value of equation.

            .. math::

                val = P \cdot f\left( \frac{P}{P_{ref}}\right)

                P = \dot{m}_{1,in} \cdot \left( h_{1,out} - h_{1,in} \right)
        """
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()

        val = i[0] * (o[2] - i[2])
        if np.isnan(bus.P_ref):
            expr = 1
        else:
            expr = abs(val / bus.P_ref)
        return val * bus.char.f_x(expr)

    def bus_deriv(self, bus):
        r"""
        Calculates the matrix of partial derivatives of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus
            TESPy bus object.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 4, len(self.inl[0].fluid.val) + 3))
        deriv[0, 0, 0] = self.numeric_deriv(self.bus_func, 'm', 0, bus=bus)
        deriv[0, 0, 2] = self.numeric_deriv(self.bus_func, 'h', 0, bus=bus)
        deriv[0, 2, 2] = self.numeric_deriv(self.bus_func, 'h', 2, bus=bus)
        return deriv

    def convergence_check(self, nw):
        r"""
        Performs a convergence check.

        Parameters
        ----------
        nw : tespy.networks.network
            The network object using this component.

        Note
        ----
        Manipulate enthalpies/pressure at inlet and outlet if not specified by
        user to match physically feasible constraints, keep fluid composition
        within feasible range and then propagates it towards the outlet.
        """
        i, o = self.inl, self.outl

        if self.ttd_l.is_set or self.ttd_u.is_set:
            fl_i1 = single_fluid(i[0].fluid.val)
            fl_i2 = single_fluid(i[1].fluid.val)
            fl_o1 = single_fluid(o[0].fluid.val)
            fl_o2 = single_fluid(o[1].fluid.val)

        if self.ttd_l.is_set:
            if isinstance(fl_o1, str):
                T_min_o1 = memorise.vrange[fl_o1][2] * 1.1
            else:
                T_min_o1 = nw.T_range_SI[0] * 1.1
            if isinstance(fl_i2, str):
                T_min_i2 = memorise.vrange[fl_i2][2] * 1.1
            else:
                T_min_i2 = nw.T_range_SI[0] * 1.1
            h_min_o1 = h_mix_pT(o[0].to_flow(), T_min_o1)
            h_min_i2 = h_mix_pT(i[1].to_flow(), T_min_i2)
            if not o[0].h.val_set and o[0].h.val_SI < h_min_o1 * 2:
                o[0].h.val_SI = h_min_o1 * 2
            if not i[1].h.val_set and i[1].h.val_SI < h_min_i2:
                i[1].h.val_SI = h_min_i2 * 1.1

        if self.ttd_u.is_set:
            if isinstance(fl_i1, str):
                T_min_i1 = memorise.vrange[fl_i1][2] * 1.1
            else:
                T_min_i1 = nw.T_range_SI[0] * 1.1
            if isinstance(fl_o2, str):
                T_min_o2 = memorise.vrange[fl_o2][2] * 1.1
            else:
                T_min_o2 = nw.T_range_SI[0] * 1.1
            h_min_i1 = h_mix_pT(i[0].to_flow(), T_min_i1)
            h_min_o2 = h_mix_pT(o[1].to_flow(), T_min_o2)
            if not i[0].h.val_set and i[0].h.val_SI < h_min_i1 * 2:
                i[0].h.val_SI = h_min_i1 * 2
            if not o[1].h.val_set and o[1].h.val_SI < h_min_o2:
                o[1].h.val_SI = h_min_o2 * 1.1

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
                4 \cdot 10^5 & \text{key = 'p'}\\
                h\left(p, 200 \text{K} \right) & \text{key = 'h' at outlet 1}\\
                h\left(p, 250 \text{K} \right) & \text{key = 'h' at outlet 2}
                \end{cases}
        """
        if key == 'p':
            return 50e5
        elif key == 'h':
            flow = [c.m.val0, c.p.val_SI, c.h.val_SI, c.fluid.val]
            if c.s_id == 'out1':
                T = 200 + 273.15
                return h_mix_pT(flow, T)
            else:
                T = 250 + 273.15
                return h_mix_pT(flow, T)

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
                4 \cdot 10^5 & \text{key = 'p'}\\
                h\left(p, 300 \text{K} \right) & \text{key = 'h' at inlet 1}\\
                h\left(p, 220 \text{K} \right) & \text{key = 'h' at outlet 2}
                \end{cases}
        """
        if key == 'p':
            return 50e5
        elif key == 'h':
            flow = [c.m.val0, c.p.val_SI, c.h.val_SI, c.fluid.val]
            if c.t_id == 'in1':
                T = 300 + 273.15
                return h_mix_pT(flow, T)
            else:
                T = 220 + 273.15
                return h_mix_pT(flow, T)

    def calc_parameters(self):
        r"""
        Postprocessing parameter calculation.
        """
        # connection information
        i1 = self.inl[0].to_flow()
        i2 = self.inl[1].to_flow()
        o1 = self.outl[0].to_flow()
        o2 = self.outl[1].to_flow()

        # temperatures
        if isinstance(self, condenser):
            T_i1 = T_bp_p(i1)
        else:
            T_i1 = T_mix_ph(i1, T0=self.inl[0].T.val_SI)
        T_i2 = T_mix_ph(i2, T0=self.inl[1].T.val_SI)
        T_o1 = T_mix_ph(o1, T0=self.outl[0].T.val_SI)
        T_o2 = T_mix_ph(o2, T0=self.outl[1].T.val_SI)

        # specific volume
        v_i1 = v_mix_ph(i1, T0=T_i1)
        v_i2 = v_mix_ph(i2, T0=T_i2)
        v_o1 = v_mix_ph(o1, T0=T_o1)
        v_o2 = v_mix_ph(o2, T0=T_o2)

        # specific entropy
        s_i1 = s_mix_ph(i1, T0=T_i1)
        s_i2 = s_mix_ph(i2, T0=T_i2)
        s_o1 = s_mix_ph(o1, T0=T_o1)
        s_o2 = s_mix_ph(o2, T0=T_o2)

        # component parameters
        self.ttd_u.val = T_i1 - T_o2
        self.ttd_l.val = T_o1 - T_i2
        self.Q.val = i1[0] * (o1[2] - i1[2])

        self.pr1.val = o1[1] / i1[1]
        self.pr2.val = o2[1] / i2[1]
        self.zeta1.val = ((i1[1] - o1[1]) * math.pi ** 2 /
                          (8 * i1[0] ** 2 * (v_i1 + v_o1) / 2))
        self.zeta2.val = ((i2[1] - o2[1]) * math.pi ** 2 /
                          (8 * i2[0] ** 2 * (v_i2 + v_o2) / 2))

        self.SQ1.val = self.inl[0].m.val_SI * (s_o1 - s_i1)
        self.SQ2.val = self.inl[1].m.val_SI * (s_o2 - s_i2)
        self.Sirr.val = self.SQ1.val + self.SQ2.val

        # kA and logarithmic temperature difference
        if T_i1 <= T_o2 or T_o1 <= T_i2:
            self.td_log.val = np.nan
            self.kA.val = np.nan
        else:
            self.td_log.val = ((T_o1 - T_i2 - T_i1 + T_o2) /
                               math.log((T_o1 - T_i2) / (T_i1 - T_o2)))
            self.kA.val = -(i1[0] * (o1[2] - i1[2]) / self.td_log.val)

        if self.kA.is_set:
            # get bound errors for kA hot side characteristics
            if self.kA_char1.param == 'm':
                i1_d = self.inl[0].to_flow_design()
                if not np.isnan(i1_d[0]):
                    if not i1[0] == 0:
                        self.kA_char1.func.get_bound_errors(i1[0] / i1_d[0],
                                                            self.label)

            # get bound errors for kA copld side characteristics
            if self.kA_char2.param == 'm':
                i2_d = self.inl[1].to_flow_design()
                if not np.isnan(i2_d[0]):
                    if not i1[0] == 0:
                        self.kA_char2.func.get_bound_errors(i2[0] / i2_d[0],
                                                            self.label)

        self.check_parameter_bounds()

# %%


class condenser(heat_exchanger):
    r"""
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.heat_exchanger.fluid_func`
        - :func:`tespy.components.components.heat_exchanger.mass_flow_func`
        - :func:`tespy.components.components.condenser.energy_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) - \dot{Q}

        - :func:`tespy.components.components.condenser.kA_func`
        - :func:`tespy.components.components.condenser.ttd_u_func`

        .. math::

            0 = p_{1,in} \cdot pr1 - p_{1,out}\\
            0 = p_{2,in} \cdot pr2 - p_{2,out}

        - :func:`tespy.components.components.heat_exchanger.zeta_func`
        - :func:`tespy.components.components.heat_exchanger.zeta2_func`

        **additional equations**

        - :func:`tespy.components.components.condenser.additional_equations`

    Inlets/Outlets

        - in1, in2 (index 1: hot side, index 2: cold side)
        - out1, out2 (index 1: hot side, index 2: cold side)

    Image

        .. image:: _images/condenser.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Q : str/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : str/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : str/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char1 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at hot side, provide
        x and y values or use generic values (e. g. calculated from design
        case). Standard method 'COND_HOT', Parameter 'm'.

    kA_char2 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at cold side,
        provide x and y values or use generic values (e. g. calculated from
        design case). Standard method 'COND_COLD', Parameter 'm'.

    subcooling : bool
        Enable/disable subcooling, default value: disabled.

    Note
    ----

    - The condenser has an additional equation for enthalpy at hot side outlet.
    - The pressure drop via zeta1 at hot side is not an offdesign parameter.
    - It has different calculation method for given heat transfer coefficient
      and upper terminal temperature difference.

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> import shutil
    >>> nw = nwk.network(fluids=['water', 'air'], T_unit='C', p_unit='bar',
    ...     h_unit='kJ / kg', m_range=[0.01, 10])
    >>> nw.set_printoptions(print_level='none')
    >>> amb_in = cmp.sink('ambient in')
    >>> amb_out = cmp.source('ambient out')
    >>> hsin = cmp.sink('HS in')
    >>> hsout = cmp.source('HS out')
    >>> he = cmp.condenser('condenser')
    >>> he.component()
    'condenser'
    >>> amb_he = con.connection(amb_out, 'out1', he, 'in2')
    >>> he_amb = con.connection(he, 'out2', amb_in, 'in1')
    >>> hs_he = con.connection(hsout, 'out1', he, 'in1')
    >>> he_hs = con.connection(he, 'out1', hsin, 'in1')
    >>> nw.add_conns(amb_he, he_amb, hs_he, he_hs)
    >>> he.set_attr(pr1=0.98, pr2=0.999, design=['pr2'],
    ...     offdesign=['zeta2', 'kA'])
    >>> hs_he.set_attr(Td_bp=20, p=1, fluid={'water': 1, 'air': 0})
    >>> amb_he.set_attr(fluid={'water': 0, 'air': 1}, T=20)
    >>> he_amb.set_attr(p=1, T=40, design=['T'])
    >>> he.set_attr(Q=-80e3)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(hs_he.m.val, 2)
    0.03
    >>> round(amb_he.m.val, 2)
    3.97
    >>> round(he_amb.T.val, 1)
    40.0
    >>> he.set_attr(Q=-60e3)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(amb_he.m.val, 2)
    2.78
    >>> round(he_amb.T.val, 1)
    41.5
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'condenser'

    def attr(self):
        return {'Q': dc_cp(max_val=0),
                'kA': dc_cp(min_val=0),
                'td_log': dc_cp(min_val=0),
                'ttd_u': dc_cp(min_val=0), 'ttd_l': dc_cp(min_val=0),
                'pr1': dc_cp(max_val=1), 'pr2': dc_cp(max_val=1),
                'zeta1': dc_cp(min_val=0), 'zeta2': dc_cp(min_val=0),
                'subcooling': dc_simple(val=False),
                'kA_char1': dc_cc(method='COND_HOT', param='m'),
                'kA_char2': dc_cc(method='COND_COLD', param='m'),
                'SQ1': dc_simple(), 'SQ2': dc_simple(), 'Sirr': dc_simple(),
                'zero_flag': dc_simple()}

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for this
        component.

        Equations

            **mandatory equations**

            .. math::

                0 = h_{1,out} - h\left(p, x=0 \right)\\
                x: \text{vapour mass fraction}

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equation for saturated liquid at hot side outlet
        if not self.subcooling.val:
            outl = self.outl
            o1 = outl[0].to_flow()
            vec_res += [o1[2] - h_mix_pQ(o1, 0)]

        return vec_res

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for saturated liquid at hot side outlet equation
        if not self.subcooling.val:
            o1 = self.outl[0].to_flow()
            x_deriv = np.zeros((1, 4, self.num_fl + 3))
            x_deriv[0, 2, 1] = -dh_mix_dpQ(o1, 0)
            x_deriv[0, 2, 2] = 1
            mat_deriv += x_deriv.tolist()

        return mat_deriv

    def energy_func(self):
        r"""
        Equation for condenser energy balance.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{1,in} \cdot \left(h_{1,out} - h_{1,in} \right) +
                \dot{m}_{2,in} \cdot \left(h_{2,out} - h_{2,in} \right)
        """
        return (self.inl[0].m.val_SI * (self.outl[0].h.val_SI -
                                        self.inl[0].h.val_SI) +
                self.inl[1].m.val_SI * (self.outl[1].h.val_SI -
                                        self.inl[1].h.val_SI))

    def energy_deriv(self):
        r"""
        Calculates the matrix of partial derivatives for energy balance
        equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 4, len(self.inl[0].fluid.val) + 3))
        for k in range(2):
            deriv[0, k, 0] = self.outl[k].h.val_SI - self.inl[k].h.val_SI
            deriv[0, k, 2] = -self.inl[k].m.val_SI

        deriv[0, 2, 2] = self.inl[0].m.val_SI
        deriv[0, 3, 2] = self.inl[1].m.val_SI
        return deriv.tolist()

    def kA_func(self):
        r"""
        Equation for heat transfer from conditions on both sides of condenser.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = \dot{m}_{1,in} \cdot \left( h_{1,out} - h_{1,in}\right) +
                kA \cdot f_{kA} \cdot \frac{T_{1,out} -
                T_{2,in} - T_s \left(p_{1,in}\right) +
                T_{2,out}}
                {\ln{\frac{T_{1,out} - T_{2,in}}
                {T_s \left(p_{1,in}\right) - T_{2,out}}}}

                f_{kA} = f_1\left(\frac{m_1}{m_{1,ref}}\right) \cdot
                f_2\left(\frac{m_2}{m_{2,ref}}\right)

        Note
        ----
        For standard functions f\ :subscript:`1` \ and f\ :subscript:`2` \ see
        class :func:`tespy.components.characteristics.characteristics`.

        - Calculate temperatures at inlets and outlets.
        - Perform value manipulation, if temperature levels are physically
          infeasible.
        """
        if self.zero_flag.val_set:
            return self.inl[0].p.val_SI - self.inl[0].p.design

        i1 = self.inl[0].to_flow()
        i2 = self.inl[1].to_flow()
        o1 = self.outl[0].to_flow()
        o2 = self.outl[1].to_flow()

        i1_d = self.inl[0].to_flow_design()
        i2_d = self.inl[1].to_flow_design()

        T_i1 = T_bp_p(i1)
        T_i2 = T_mix_ph(i2, T0=self.inl[1].T.val_SI)
        T_o1 = T_mix_ph(o1, T0=self.outl[0].T.val_SI)
        T_o2 = T_mix_ph(o2, T0=self.outl[1].T.val_SI)

        if T_i1 <= T_o2 and not self.inl[0].T.val_set:
            T_i1 = T_o2 + 0.5
        if T_i1 <= T_o2 and not self.outl[1].T.val_set:
            T_o2 = T_i1 - 0.5

        if T_o1 <= T_i2 and not self.outl[0].T.val_set:
            T_o1 = T_i2 + 1
        if T_o1 <= T_i2 and not self.inl[1].T.val_set:
            T_i2 = T_o1 - 1

        fkA1 = 1
        if self.kA_char1.param == 'm':
            if not np.isnan(i1_d[0]):
                fkA1 = self.kA_char1.func.f_x(i1[0] / i1_d[0])

        fkA2 = 1
        if self.kA_char2.param == 'm':
            if not np.isnan(i2_d[0]):
                fkA2 = self.kA_char2.func.f_x(i2[0] / i2_d[0])

        td_log = ((T_o1 - T_i2 - T_i1 + T_o2) /
                  math.log((T_o1 - T_i2) / (T_i1 - T_o2)))
        return i1[0] * (o1[2] - i1[2]) + self.kA.val * fkA1 * fkA2 * td_log

    def ttd_u_func(self):
        r"""
        Equation for upper terminal temperature difference.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                res = ttd_{u} - T_s \left(p_{1,in}\right) + T_{2,out}

        Note
        ----
        The upper terminal temperature difference ttd_u refers to boiling
        temperature at hot side inlet.
        """
        i1 = self.inl[0].to_flow()
        o2 = self.outl[1].to_flow()
        T_o2 = T_mix_ph(o2, T0=self.outl[1].T.val_SI)
        return self.ttd_u.val - T_bp_p(i1) + T_o2

# %%


class desuperheater(heat_exchanger):
    r"""
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.heat_exchanger.fluid_func`
        - :func:`tespy.components.components.heat_exchanger.mass_flow_func`
        - :func:`tespy.components.components.heat_exchanger.energy_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) - \dot{Q}

        - :func:`tespy.components.components.heat_exchanger.kA_func`
        - :func:`tespy.components.components.heat_exchanger.ttd_u_func`

        .. math::

            0 = p_{1,in} \cdot pr1 - p_{1,out}\\
            0 = p_{2,in} \cdot pr2 - p_{2,out}

        - :func:`tespy.components.components.heat_exchanger.zeta_func`
        - :func:`tespy.components.components.heat_exchanger.zeta2_func`

        **additional equations**

        - :func:`tespy.components.components.desuperheater.additional_equations`

    Inlets/Outlets

        - in1, in2 (index 1: hot side, index 2: cold side)
        - out1, out2 (index 1: hot side, index 2: cold side)

    Image

        .. image:: _images/heat_exchanger.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Q : str/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : str/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : str/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char1 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at hot side, provide
        x and y values or use generic values (e. g. calculated from design
        case). Standard method 'COND_HOT', Parameter 'm'.

    kA_char2 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at cold side,
        provide x and y values or use generic values (e. g. calculated from
        design case). Standard method 'COND_COLD', Parameter 'm'.

    Note
    ----
    The desuperheater has an additional equation for enthalpy at hot side
    outlet.

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> import shutil
    >>> nw = nwk.network(fluids=['water', 'air'], T_unit='C', p_unit='bar',
    ...     h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> amb_in = cmp.sink('ambient in')
    >>> amb_out = cmp.source('ambient out')
    >>> hsin = cmp.sink('HS in')
    >>> hsout = cmp.source('HS out')
    >>> he = cmp.desuperheater('desuperheater')
    >>> he.component()
    'desuperheater'
    >>> amb_he = con.connection(amb_out, 'out1', he, 'in2')
    >>> he_amb = con.connection(he, 'out2', amb_in, 'in1')
    >>> hs_he = con.connection(hsout, 'out1', he, 'in1')
    >>> he_hs = con.connection(he, 'out1', hsin, 'in1')
    >>> nw.add_conns(amb_he, he_amb, hs_he, he_hs)
    >>> he.set_attr(pr1=0.98, pr2=0.999, design=['pr1', 'pr2'],
    ...     offdesign=['zeta1', 'zeta2', 'kA'])
    >>> hs_he.set_attr(T=200, p=1, fluid={'water': 1, 'air': 0})
    >>> amb_he.set_attr(fluid={'water': 0, 'air': 1}, T=20)
    >>> he_amb.set_attr(p=1, T=40, design=['T'])
    >>> he.set_attr(Q=-80e3)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(hs_he.m.val, 1)
    0.4
    >>> round(amb_he.m.val, 2)
    3.97
    >>> round(he_amb.T.val, 1)
    40.0
    >>> he.set_attr(Q=-60e3)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(hs_he.m.val, 1)
    0.3
    >>> round(amb_he.m.val, 2)
    2.56
    >>> round(he_amb.T.val, 1)
    43.3
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'desuperheater'

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for this
        component.

        Equations

            **mandatory equations**

            .. math::

                0 = h_{1,out} - h\left(p, x=1 \right)\\
                x: \text{vapour mass fraction}

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equation for saturated gas at hot side outlet
        o1 = self.outl[0].to_flow()
        vec_res += [o1[2] - h_mix_pQ(o1, 1)]

        return vec_res

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional
        equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for saturated gas at hot side outlet equation
        o1 = self.outl[0].to_flow()
        deriv = np.zeros((1, 4, self.num_fl + 3))
        deriv[0, 2, 1] = -dh_mix_dpQ(o1, 1)
        deriv[0, 2, 2] = 1
        mat_deriv += deriv.tolist()

        return mat_deriv

# %%


class drum(component):
    r"""
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.drum.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        .. math::

            0 = \sum_i \left(\dot{m}_{i,in} \cdot h_{i,in} \right) -
            \sum_j \left(\dot{m}_{j,out} \cdot h_{j,out} \right)\\
            \forall i \in inlets, \; \forall j \in outlet

            0 = p_{in,1} - p_{out,i}\\
            \forall i \in \mathrm{outlets}

            0 = h_{1,out} - h\left(p, x=0 \right)

            0 = h_{2,out} - h\left(p, x=1 \right)\\
            x: \text{vapour mass fraction}

    Inlets/Outlets

        - in1, in2 (index 1: from economiser, index 2: from evaporator)
        - out1, out2 (index 1: to evaporator, index 2: to superheater)

    Image

        .. image:: _images/drum.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Note
    ----
    If you are using a drum in a network with multiple fluids, it is likely
    the fluid propagation causes trouble. If this is the case, try to
    specify the fluid composition at another connection of your network.

    This component assumes, that the fluid composition between outlet 1 and
    inlet 2 does not change, thus there is no equation for the fluid mass
    fraction at the inlet 2!

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> import shutil
    >>> nw = nwk.network(fluids=['NH3', 'air'], T_unit='C', p_unit='bar',
    ...     h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> f = cmp.source('feed')
    >>> ha = cmp.source('hot air')
    >>> ch = cmp.sink('chimney')
    >>> s = cmp.sink('steam')
    >>> dr = cmp.drum('drum')
    >>> dr.component()
    'drum'
    >>> ev = cmp.heat_exchanger('evaporator')
    >>> erp = cmp.pump('evaporator reciculation pump')
    >>> f_dr = con.connection(f, 'out1', dr, 'in1')
    >>> dr_erp = con.connection(dr, 'out1', erp, 'in1')
    >>> erp_ev = con.connection(erp, 'out1', ev, 'in2')
    >>> ev_dr = con.connection(ev, 'out2', dr, 'in2')
    >>> dr_s = con.connection(dr, 'out2', s, 'in1')
    >>> nw.add_conns(f_dr, dr_erp, erp_ev, ev_dr, dr_s)
    >>> ha_ev = con.connection(ha, 'out1', ev, 'in1')
    >>> ev_ch = con.connection(ev, 'out1', ch, 'in1')
    >>> nw.add_conns(ha_ev, ev_ch)
    >>> ev.set_attr(pr1=0.999, pr2=0.99, ttd_l=20, kA_char1='EVA_HOT',
    ...     kA_char2='EVA_COLD', design=['pr1', 'ttd_l'],
    ...     offdesign=['zeta1', 'kA'])
    >>> ev.set_attr(Q=-1e6)
    >>> erp.set_attr(eta_s=0.8)
    >>> f_dr.set_attr(p=5, T=-5)
    >>> erp_ev.set_attr(m=con.ref(f_dr, 4, 0), fluid={'air': 0, 'NH3': 1})
    >>> ha_ev.set_attr(fluid={'air': 1, 'NH3': 0}, T=100)
    >>> ev_ch.set_attr(p=1)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(ev.ttd_l.val, 1)
    20.0
    >>> round(f_dr.h.val, 1)
    320.2
    >>> round(dr_erp.h.val, 1)
    362.4
    >>> round(ev_dr.h.val, 1)
    684.7
    >>> round(f_dr.m.val, 2)
    0.78
    >>> ev.set_attr(Q=-0.75e6)
    >>> nw.solve('offdesign', init_path='tmp', design_path='tmp')
    >>> round(f_dr.m.val, 2)
    0.58
    >>> round(ev.ttd_l.val, 1)
    16.1
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'drum'

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1', 'out2']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.fl_deriv = self.fluid_deriv()
        self.m_deriv = self.mass_flow_deriv()
        self.p_deriv = self.pressure_deriv()

    def equations(self):
        r"""
        Calculates vector vec_res with results of equations for this component.

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # eqations for fluid balance
        vec_res += self.fluid_func()

        ######################################################################
        # eqations for mass flow balance
        vec_res += self.mass_flow_func()

        ######################################################################
        # eqations for pressure
        p = self.inl[0].p.val_SI
        for c in [self.inl[1]] + self.outl:
            vec_res += [p - c.p.val_SI]

        ######################################################################
        # eqations for enthalpy
        val = 0
        for i in self.inl:
            val += i.m.val_SI * i.h.val_SI
        for o in self.outl:
            val -= o.m.val_SI * o.h.val_SI
        vec_res += [val]

        ######################################################################
        # eqations for staturated fluid state at outlets
        vec_res += [h_mix_pQ(self.outl[0].to_flow(), 0) -
                    self.outl[0].h.val_SI]
        vec_res += [h_mix_pQ(self.outl[1].to_flow(), 1) -
                    self.outl[1].h.val_SI]

        return vec_res

    def derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for fluid balance equations
        mat_deriv += self.fl_deriv

        ######################################################################
        # derivatives for mass flow balance equation
        mat_deriv += self.m_deriv

        ######################################################################
        # derivatives for pressure eqauations
        mat_deriv += self.p_deriv

        ######################################################################
        # derivatives for energy balance equation
        deriv = np.zeros((1, 4, self.num_fl + 3))
        k = 0
        for i in self.inl:
            deriv[0, k, 0] = i.h.val_SI
            deriv[0, k, 2] = i.m.val_SI
            k += 1
        j = 0
        for o in self.outl:
            deriv[0, j + k, 0] = -o.h.val_SI
            deriv[0, j + k, 2] = -o.m.val_SI
            j += 1
        mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives of equations for saturated states at outlets
        x_deriv = np.zeros((2, 4, self.num_fl + 3))
        x_deriv[0, 2, 1] = dh_mix_dpQ(self.outl[0].to_flow(), 0)
        x_deriv[0, 2, 2] = -1
        x_deriv[1, 3, 1] = dh_mix_dpQ(self.outl[1].to_flow(), 1)
        x_deriv[1, 3, 2] = -1
        mat_deriv += x_deriv.tolist()

        return np.asarray(mat_deriv)

    def fluid_func(self):
        r"""
        Calculates the vector of residual values for component's fluid balance
        equations.

        Returns
        -------
        vec_res : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = fluid_{i,in_1} - fluid_{i,out_{j}}\\
                \forall i \in \mathrm{fluid}, \; \forall j \in inlets

        """
        vec_res = []

        for o in self.outl:
            for fluid, x in self.inl[0].fluid.val.items():
                vec_res += [x - o.fluid.val[fluid]]
        return vec_res

    def fluid_deriv(self):
        r"""
        Calculates the partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((2 * self.num_fl, 4, 3 + self.num_fl))
        for k in range(2):
            for i in range(self.num_fl):
                deriv[i + k * self.num_fl, 0, i + 3] = 1
                deriv[i + k * self.num_fl, k + 2, i + 3] = -1
        return deriv.tolist()

    def pressure_deriv(self):
        r"""
        Calculates the partial derivatives for pressure equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((3, 4, self.num_fl + 3))
        for k in range(3):
            deriv[k, 0, 1] = 1
            deriv[k, k + 1, 1] = -1
        return deriv.tolist()

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
                10^6 & \text{key = 'p'}\\
                h\left(p, x=0 \right) & \text{key = 'h' at outlet 1}\\
                h\left(p, x=1 \right) & \text{key = 'h' at outlet 2}
                \end{cases}
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            if c.s_id == 'out1':
                return h_mix_pQ(c.to_flow(), 0)
            else:
                return h_mix_pQ(c.to_flow(), 1)

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
                10^6 & \text{key = 'p'}\\
                h\left(p, x=0 \right) & \text{key = 'h' at inlet 1}\\
                h\left(p, x=0.7 \right) & \text{key = 'h' at inlet 2}
                \end{cases}
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            if c.t_id == 'in1':
                return h_mix_pQ(c.to_flow(), 0)
            else:
                return h_mix_pQ(c.to_flow(), 0.7)

# %%


class subsys_interface(component):
    r"""
    Equations

        **mandatory equations**

        .. math:: 0 = fluid_{i,in_{j}} - fluid_{i,out_{j}} \;
            \forall i \in \mathrm{fluid}, \; \forall j \in inlets/outlets

        .. math:: 0 = \dot{m}_{in_{j}} - \dot{m}_{out_{j}} \;
            \forall j \in inlets/outlets

        .. math:: 0 = p_{in_{j}} - p_{out_{j}} \;
            \forall j \in inlets/outlets

        .. math:: 0 = h_{in_{j}} - h_{out_{j}} \;
            \forall j \in inlets/outlets

    Inlets/Outlets

        - Specify number of inlets and outlets with :code:`num_inter`,
          predefined value: 1.

    Image

        .. image:: _images/subsys_interface.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    num_inter : float/tespy.helpers.dc_simple
        Number of interfaces for subsystem.

    Note
    ----
    This component passes all fluid properties and mass flow from its inlet to
    the outlet.

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> fluids = ['H2O', 'N2']
    >>> nw = nwk.network(fluids=fluids)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> so1 = cmp.source('source 1')
    >>> si1 = cmp.sink('sink 1')
    >>> si = cmp.subsys_interface('test', num_inter=1)
    >>> si2 = cmp.subsys_interface('test2', num_inter=np.nan)
    >>> si.component()
    'subsystem interface'
    >>> len(si.inlets()) == len(si2.inlets())
    True
    >>> inc = con.connection(so1, 'out1', si, 'in1')
    >>> outg = con.connection(si, 'out1', si1, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> inc.set_attr(fluid={'H2O': 1, 'N2': 0}, T=40, p=3, m=100)
    >>> nw.solve('design')
    >>> nw.iter
    2
    >>> nw.lin_dep
    False
    """

    def component(self):
        return 'subsystem interface'

    def attr(self):
        return {'num_inter': dc_simple()}

    def inlets(self):
        if self.num_inter.val_set:
            return ['in' + str(i + 1) for i in range(self.num_inter.val)]
        else:
            return ['in1']

    def outlets(self):
        if self.num_inter.val_set:
            return ['out' + str(i + 1) for i in range(self.num_inter.val)]
        else:
            return ['out1']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        # retrieve always constant derivatives
        self.fl_deriv = self.fluid_deriv()
        self.m_deriv = self.inout_deriv(0)
        self.p_deriv = self.inout_deriv(1)
        self.h_deriv = self.inout_deriv(2)

    def equations(self):
        r"""
        Calculates vector vec_res with results of equations for this component.

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # eqations for fluids
        for i in range(self.num_i):
            for fluid, x in self.inl[i].fluid.val.items():
                vec_res += [x - self.outl[i].fluid.val[fluid]]

        ######################################################################
        # equations for mass flow
        for i in range(self.num_i):
            vec_res += [self.inl[i].m.val_SI - self.outl[i].m.val_SI]

        ######################################################################
        # equations for pressure
        for i in range(self.num_i):
            vec_res += [self.inl[i].p.val_SI - self.outl[i].p.val_SI]

        ######################################################################
        # equations for enthalpy
        for i in range(self.num_i):
            vec_res += [self.inl[i].h.val_SI - self.outl[i].h.val_SI]

        ######################################################################

        return vec_res

    def derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        ######################################################################
        # derivatives with constant value (all for this component)
        mat_deriv = self.fl_deriv + self.m_deriv + self.p_deriv + self.h_deriv

        return np.asarray(mat_deriv)

    def fluid_deriv(self):
        r"""
        Calculates the partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((
                self.num_fl * self.num_i, 2 * self.num_i, 3 + self.num_fl))
        for i in range(self.num_i):
            for j in range(self.num_fl):
                deriv[i * self.num_fl + j, i, j + 3] = 1
                deriv[i * self.num_fl + j, self.num_i + i, j + 3] = -1
        return deriv.tolist()

    def inout_deriv(self, pos):
        r"""
        Calculates the partial derivatives for all mass flow, pressure and
        enthalpy equations.

        Parameters
        ----------
        pos : int
            Position of the variable in the matrix of derivatives.
            mass flow: 0, pressure: 1, enthalpy: 2.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((self.num_i, 2 * self.num_i, self.num_fl + 3))
        for i in range(self.num_i):
            deriv[i, i, pos] = 1
        for j in range(self.num_i):
            deriv[j, j + self.num_i, pos] = -1
        return deriv.tolist()
