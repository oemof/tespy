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

from tespy.tools.helpers import (
    num_fluids, fluid_structure, TESPyComponentError, tespy_fluid,
    v_mix_ph, h_mix_pT, h_mix_ps, s_mix_pT, s_mix_ph, T_mix_ph, visc_mix_ph,
    dT_mix_dph, dT_mix_pdh, dT_mix_ph_dfluid, h_mix_pQ, dh_mix_dpQ,
    h_ps, h_pT ,s_ph, s_pT,
    molar_mass_flow, lamb,
    molar_masses, err,
    dc_cp, dc_cc, dc_cm, dc_gcp, memorise, single_fluid, dc_simple, data_container
)
from tespy.components import characteristics as cmp_char

# %%


class component:
    r"""
    Class component is the base class of all TESPy components.

    Parameters
    ----------
    label : str
        The label of the component.

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    **kwargs :
        See the class documentation of desired component for available keywords.

    Note
    ----
    The initialisation method (__init__), setter method (set_attr) and getter method (get_attr)
    are used for instances of class component and its children.

    Allowed keywords in kwargs are 'mode', 'design' and 'offdesign'. Additional
    keywords depend on the type of component you want to create.

    Example
    -------
    Basic example for a setting up a tespy.components.components.component object.
    This example does not run a tespy calculation.

    >>> from tespy import cmp
    >>> comp = cmp.component('myComponent')
    >>> comp.set_attr(mode='man')
    >>> type(comp)
    <class 'tespy.components.components.component'>
    >>> comp.get_attr('mode')
    'man'
    """

    def __init__(self, label, **kwargs):

        # check if components label is of type str and for prohibited chars
        if not isinstance(label, str):
            msg = 'Component label must be of type str!'
            logging.error(msg)
            raise ValueError(msg)
        elif len([x for x in [';', ',', '.'] if x in label]) > 0:
            msg = ('Can\'t use ' + str([';', ',', '.']) + ' in label (' + str(self.component()) + ').')
            logging.error(msg)
            raise ValueError(msg)
        else:
            self.label = label

        self.mode = 'auto'

        # set default design and offdesign parameters
        self.design = []
        self.offdesign = []
        self.interface = False

        # add container for components attributes
        var = self.attr()

        for key in var.keys():
            self.__dict__.update({key: var[key]})

        self.set_attr(**kwargs)

    def set_attr(self, **kwargs):
        r"""
        Sets, resets or unsets attributes of a component for provided keyword arguments.

        Parameters
        ----------
        mode : str
            'auto' for automatic design to offdesign switch, 'man' for manual switch.

        design : list
            List containing design parameters (stated as String).

        offdesign : list
            List containing offdesign parameters (stated as String).

        **kwargs :
            See the class documentation of desired component for available keywords.

        Note
        ----
        Allowed keywords in kwargs are obtained from class documentation as all
        components share the :func:`tespy.components.components.component.set_attr` method.
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
                        msg = ('The keyword ' + key + ' expects a data_container of type ' + str(type(self.get_attr(key))) +
                               ', a data_container of type ' + str(type(kwargs[key])) + ' was supplied.')
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
                        msg = ('Bad datatype for keyword argument ' + key + ' at ' + self.label + '.')
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

                    # invalid datatype for keyword
                    else:
                        msg = ('Bad datatype for keyword argument ' + key + ' at ' + self.label + '.')
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
                            self.get_attr(key).set_attr(val=kwargs[key], val_set=True)
                    else:
                        self.get_attr(key).set_attr(val=kwargs[key], val_set=True)

            # export sources or sinks as subsystem interface
            elif key == 'interface':
                if isinstance(self, source) or isinstance(self, sink):
                    if isinstance(kwargs[key], bool):
                        self.interface = kwargs[key]
                    else:
                        msg = ('Datatype for keyword argument ' + str(key) + ' must be bool at ' + self.label + '.')
                        logging.error(msg)
                        raise ValueError(msg)
                else:
                    msg = ('Only sinks and sources can be attributed with the interface parameter (error at ' + self.label + ').')
                    logging.error(msg)
                    raise TESPyComponentError(msg)

            elif key == 'design' or key == 'offdesign':
                if not isinstance(kwargs[key], list):
                    msg = ('Please provide the ' + key + ' parameters as list at ' + self.label + '.')
                    logging.error(msg)
                    raise TypeError(msg)
                if set(kwargs[key]).issubset(list(var)):
                    self.__dict__.update({key: kwargs[key]})
                else:
                    msg = ('Available parameters for (off-)design specification are: ' +
                           str(list(var)) + ' at ' + self.label + '.')
                    logging.error(msg)
                    raise ValueError(msg)

            elif key == 'mode':
                if kwargs[key] in ['man', 'auto']:
                    self.__dict__.update({key: kwargs[key]})
                else:
                    msg = ('Mode must be \'man\' or \'auto\' at ' + self.label + '.')
                    logging.error(msg)
                    raise ValueError(msg)

            # invalid keyword
            else:
                msg = ('Component ' + self.label + ' has no attribute ' + str(key) + '.')
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
            msg = 'Component ' + self.label + ' has no attribute \"' + key + '\".'
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

        msg = 'Component ' + self.label + ' has ' + str(self.num_vars) + ' custom variables.'
        logging.debug(msg)


        # characteristics creation
        for key, val in self.attr().items():
            if isinstance(val, dc_cc):
                generate_char = False
                if self.get_attr(key).func is None:
                    generate_char = True
                elif (not np.array_equal(self.get_attr(key).func.x, self.get_attr(key).x) or
                      not np.array_equal(self.get_attr(key).func.y, self.get_attr(key).y)):
                    generate_char = True

                if generate_char:
                    self.get_attr(key).func = cmp_char.characteristics(
                            method=self.get_attr(key).method, x=self.get_attr(key).x,
                            y=self.get_attr(key).y, comp=self.component())
                    self.get_attr(key).x = self.get_attr(key).func.x
                    self.get_attr(key).y = self.get_attr(key).func.y

                    msg = 'Generated characteristic line for attribute ' + key + ' at component ' + self.label + '.'
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

    def bus_func(self, bus):
        return 0

    def bus_deriv(self, bus):
        return

    def initialise_source(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's outlet.

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
        Returns a starting value for pressure and enthalpy at component's inlet.

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
            Setting component design values for :code:`mode='offdesign'` and unsetting them for :code:`mode='design'`.

        df : pandas.core.series.Series
            Series containing the component parameters.
        """
        for key, dc in self.attr().items():
            if isinstance(dc, dc_cp):
                if mode == 'offdesign':
                    self.get_attr(key).design = data[key]
                else:
                    self.get_attr(key).design = np.nan

    def calc_parameters(self, mode):
        r"""
        Post and preprocessing parameter calculation/specification.

        Parameters
        ----------

        mode : str
            Pre- or postprocessing calculation.

        Note
        ----
        The method :func:`tespy.components.components.component.calc_parameters` is the base method called by every component specific method.
        This method is for preprocessing of offdesign calculations and sets all component attributes provided as offdesign parameters to their design value.

        Postprocessing is handled by the component specific methods.
        """
        if mode == 'pre':
            # set component attributes to design-value if specified as offdesign parameter
            switched = False
            msg = 'Set component attributes '
            for key, dc in self.attr().items():
                if isinstance(dc, dc_cp) and key in self.offdesign:
                    switched = True
                    self.get_attr(key).val = self.get_attr(key).design

                    msg += key + ', '

            msg = msg[:-2] + ' to design value at component ' + self.label + '.'
            if switched:
                logging.debug(msg)

    def initialise_fluids(self, nw):
        return

    def convergence_check(self, nw):
        return

# %%

    def fluid_func(self):
        r"""
        Calculates the vector of residual values for component's fluid balance equations.

        Returns
        -------
        vec_res : list
            Vector of residual values for component's fluid balance.

            .. math::
                0 = fluid_{i,in} - fluid_{i,out} \; \forall i \in \mathrm{fluid}
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
        Calculates the residual value for component's mass flow balance equation.

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
            Matrix with partial derivatives for the mass flow balance equations.
        """

        deriv = np.zeros((1, self.num_i + self.num_o + self.num_vars, 3 + self.num_fl))
        for i in range(self.num_i):
            deriv[0, i, 0] = 1
        for j in range(self.num_o):
            deriv[0, j + i + 1, 0] = -1
        return deriv.tolist()

# %%

    def numeric_deriv(self, func, dx, pos, **kwargs):
        r"""
        Calculates partial derivative of the function func to dx at given connection.

        Parameters
        ----------
        func : function
            Function :math:`f` to calculate the partial derivative for.

        dx : str
            Partial derivative.

        pos : int
            Position of connection regarding to inlets and outlet of the component,
            logic: ['in1', 'in2', ..., 'out1', ...] -> 0, 1, ..., n, n + 1, ..., n + m

        Returns
        -------
        deriv : float/list
            Partial derivative(s) of the function :math:`f` to variable(s) :math:`x`.

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

    def zeta_func(self ):
        r"""
        Calculates residual value of :math:`\zeta`-function.

        Returns
        -------
        val : float
            Residual value of function.

            .. math::

                val = \begin{cases}
                p_{in} - p_{out} & |\dot{m}| < \epsilon \\
                \zeta - \frac{(p_{in} - p_{out}) \cdot \pi^2}{8 \cdot
                \dot{m}_{in} \cdot |\dot{m}_{in}| \cdot \frac{v_{in} + v_{out}}{2}} &
                |\dot{m}| > \epsilon
                \end{cases}

        Note
        ----
        The zeta value is caluclated on the basis of a given pressure loss at a given flow rate.
        As the cross sectional area A will not change, it is possible to handle the equation in this way.

        .. math::

            \zeta = \frac{\Delta p \cdot v \cdot 2}{c^2}\\
            c = \frac{\dot{m} \cdot v}{A}
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
            return (val - (i[1] - o[1]) * math.pi ** 2 /
                    (8 * abs(i[0]) * i[0] * (v_mix_ph(i) + v_mix_ph(o)) / 2))

    def zeta2_func(self):
        r"""
        calculates residual value of :math:`\zeta`-function (for heat exchangers at lower temperature side).

        Returns
        -------
        val : float
            Residual value of function.

            .. math::

                val = \begin{cases}
                p_{in} - p_{out} & |\dot{m}| < \epsilon \\
                \zeta_2 - \frac{(p_{2,in} - p_{2,out}) \cdot \pi^2}{8 \cdot
                \dot{m}_{2,in} \cdot |\dot{m}_{2,in}| \cdot \frac{v_{2,in} + v_{2,out}}{2}} &
                |\dot{m}| > \epsilon
                \end{cases}

        Note
        ----
        The zeta value is caluclated on the basis of a given pressure loss at a given flow rate.
        As the cross sectional area A will not change, it is possible to handle the equation in this way.

        .. math::

            \zeta_2 = \frac{\Delta p_2 \cdot v_2 \cdot 2}{c_2^2}\\
            c_2 = \frac{\dot{m}_2 \cdot v_2}{A_2}
        """
        i = self.inl[1].to_flow()
        o = self.outl[1].to_flow()

        if abs(i[0]) < 1e-4:
            return i[1] - o[1]
        else:

            return (self.zeta2.val - (i[1] - o[1]) * math.pi ** 2 /
                    (8 * abs(i[0]) * i[0] * (v_mix_ph(i) + v_mix_ph(o)) / 2))

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

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).
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

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).
    """

    def component(self):
        return 'sink'

    def inlets(self):
        return ['in1']

# %%


class turbomachine(component):
    r"""
    The component turbomachine is the parent class for pump, compressor and turbine.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in} \right) - P\\
            0 = pr \cdot p_{in} - p_{out}

        isentropic efficiency equations (optional)

        - :func:`tespy.components.components.pump.eta_s_func`
        - :func:`tespy.components.components.compressor.eta_s_func`
        - :func:`tespy.components.components.turbine.eta_s_func`

        **additional equations**

        - :func:`tespy.components.components.pump.additional_equations`
        - :func:`tespy.components.components.compressor.additional_equations`
        - :func:`tespy.components.components.turbine.additional_equations`

    Inlets/Outlets

        - in1
        - out1

    Parameters
    ----------
    label : str
        The label of the component.

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    P : Sring/float/tespy.helpers.dc_cp
        Power, :math:`P/\text{W}`

    eta_s : Sring/float/tespy.helpers.dc_cp
        Isentropic efficiency, :math:`\eta_s/1`

    pr : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`

    eta_s_char : str/tespy.helpers.dc_cc
        Characteristic curve for isentropic efficiency, provide x and y values
        or use generic values (e. g. calculated from design case).
    """

    def component(self):
        return 'turbomachine'

    def attr(self):
        return {'P': dc_cp(), 'eta_s': dc_cp(), 'pr': dc_cp(),
                'eta_s_char': dc_cc(), 'Sirr': dc_cp()}

    def inlets(self):
        return ['in1']

    def outlets(self):
        return ['out1']

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
        # eqations for fluids
        vec_res += self.fluid_func()

        ######################################################################
        # eqations for mass flow balance
        vec_res += self.mass_flow_func()

        ######################################################################
        # eqations for specified power
        if self.P.is_set:
            vec_res += [self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI) - self.P.val]

        ######################################################################
        # eqations for specified pressure ratio
        if self.pr.is_set:
            vec_res += [self.pr.val * self.inl[0].p.val_SI - self.outl[0].p.val_SI]

        ######################################################################
        # eqations for specified isentropic efficiency
        if self.eta_s.is_set:
            self.eta_s_res = self.eta_s_func()
            vec_res += [self.eta_s_res]

        ######################################################################
        # additional equations
        vec_res += self.additional_equations()

        return vec_res

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for this component.

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
        # derivatives fluid composition
        mat_deriv += self.fl_deriv

        ######################################################################
        # derivatives for mass flow balance
        mat_deriv += self.m_deriv

        ######################################################################
        # derivatives for specified power
        if self.P.is_set:
            P_deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            P_deriv[0, 0, 0] = self.outl[0].h.val_SI - self.inl[0].h.val_SI
            P_deriv[0, 0, 2] = -self.inl[0].m.val_SI
            P_deriv[0, 1, 2] = self.inl[0].m.val_SI
            mat_deriv += P_deriv.tolist()

        ######################################################################
        # derivatives for specified pressure ratio
        if self.pr.is_set:
            pr_deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            pr_deriv[0, 0, 1] = self.pr.val
            pr_deriv[0, 1, 1] = -1
            mat_deriv += pr_deriv.tolist()

        ######################################################################
        # derivatives for specified isentropic efficiency
        if self.eta_s.is_set:
            mat_deriv += self.eta_s_deriv()

        ######################################################################
        # derivatives for additional equations
        mat_deriv += self.additional_derivatives()

        return np.asarray(mat_deriv)

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        return []

    def eta_s_func(self):
        r"""
        Calculates residual value of isentropic efficiency function, see subclasses.
        """
        msg = 'If you want to use eta_s as parameter, please specify which type of turbomachine you are using.'
        logging.error(msg)
        raise TESPyComponentError(msg)

    def eta_s_deriv(self):
        r"""
        Calculates partial derivatives for isentropic efficiency function, see subclasses.
        """
        msg = 'If you want to use eta_s as parameter, please specify which type of turbomachine you are using.'
        logging.error(msg)
        raise TESPyComponentError(msg)

    def h_os(self, mode):
        r"""
        Calculates the enthalpy at the outlet if compression or expansion is isentropic.

        Parameters
        ----------
        mode : str
            Determines wheather calculation is in preprocessing mode.

        Returns
        -------
        h : float
            Enthalpy after isentropic state change.

            .. math::

                h = \begin{cases}
                h\left(p_{out}, s\left(p_{in}, h_{in}\right) \right) & \text{pure fluids}\\
                h\left(p_{out}, s\left(p_{in}, T_{in}\right) \right) & \text{mixtures}\\
                \end{cases}
        """
        if mode == 'pre':
            i = self.inl[0].to_flow_design()
            o = self.outl[0].to_flow_design()
        else:
            i = self.inl[0].to_flow()
            o = self.outl[0].to_flow()

        if num_fluids(i[3]) == 1:
            for fluid, x in i[3].items():
                if x > err:
                    return h_ps(o[1], s_ph(i[1], i[2], fluid), fluid)
        else:
            T_mix = T_mix_ph(i)
            s_mix = s_mix_pT(i, T_mix)
            return h_mix_ps(o, s_mix)

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

                P = \dot{m}_{in} \cdot \left( h_{out} - h_{in} \right)\\
                val = P \cdot f_{char}\left( \frac{P}{P_{ref}}\right)
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

    def calc_parameters(self, mode):
        r"""
        Post and preprocessing parameter calculation/specification.

        Parameters
        ----------

        mode : str
            Pre- or postprocessing calculation.

        Note
        ----
        Generic preprocessing is handled by the base class. This method handles class specific pre- and postprocessing.
        """
        component.calc_parameters(self, mode)

        if mode == 'post':
            i, o = self.inl[0].to_flow(), self.outl[0].to_flow()
            self.P.val = i[0] * (o[2] - i[2])
            self.pr.val = o[1] / i[1]
            self.Sirr.val = self.inl[0].m.val_SI * (
                    s_mix_ph(self.outl[0].to_flow()) -
                    s_mix_ph(self.inl[0].to_flow()))

        else:
            self.dh_s_ref = (self.h_os(mode) - self.inl[0].h.design)

# %%


class pump(turbomachine):
    r"""
    The component turbomachine is the parent class for pump, compressor and turbine.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in} \right) - P\\
            0 = pr \cdot p_{in} - p_{out}

        - :func:`tespy.components.components.pump.eta_s_func`

        **additional equations**

        - :func:`tespy.components.components.pump.additional_equations`

    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/pump.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    P : Sring/float/tespy.helpers.dc_cp
        Power, :math:`P/\text{W}`

    eta_s : Sring/float/tespy.helpers.dc_cp
        Isentropic efficiency, :math:`\eta_s/1`

    pr : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`

    eta_s_char : str/tespy.helpers.dc_cc
        Characteristic curve for isentropic efficiency, provide x and y values
        or use generic values (e. g. calculated from design case).

    flow_char : str/tespy.helpers.dc_cc
        Characteristic curve for pressure rise vs. volumetric flow rate,
        provide data: :math:`x/\frac{\text{m}^3}{\text{s}} \,
        y/\text{Pa}`

    Example
    -------
    >>> from tespy import cmp, con, nwk, hlp
    >>> import shutil
    >>> fluid_list = ['water']
    >>> nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ...     h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> si = cmp.sink('sink')
    >>> so = cmp.source('source')
    >>> p = cmp.pump('pump')
    >>> inc = con.connection(so, 'out1', p, 'in1')
    >>> outg = con.connection(p, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> v = np.array([0, 0.4, 0.8, 1.2, 1.6, 2]) / 1000
    >>> dp = np.array([15, 14, 12, 9, 5, 0]) * 1e5
    >>> char = hlp.dc_cc(x=v, y=dp, is_set=True)
    >>> p.set_attr(pr=10, eta_s=0.8, flow_char=char, design=['eta_s'],
    ...     offdesign=['eta_s_char'])
    >>> inc.set_attr(fluid={'water': 1}, p=1, T=20)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> p.set_attr(pr=14)
    >>> round(inc.m.val_SI, 3)
    1.198
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(inc.m.val_SI, 3)
    0.599
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'pump'

    def attr(self):
        return {'P': dc_cp(), 'eta_s': dc_cp(), 'pr': dc_cp(), 'Sirr': dc_cp(),
                'eta_s_char': dc_cc(method='GENERIC'),
                'flow_char': dc_cc()}

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for pump.

        Equations

            **optional equations**

            - :func:`tespy.components.components.pump.eta_s_char_func`
            - :func:`tespy.components.components.pump.flow_char_func`

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equations for specified isentropic efficiency characteristics
        if self.eta_s_char.is_set:
            vec_res += [self.eta_s_char_func()]

        ######################################################################
        # equations for specified pressure rise vs. flowrate characteristics
        if self.flow_char.is_set:
            vec_res += [self.flow_char_func()]

        return vec_res

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional equations.

        Returns
        -------
        mat_deriv : list
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for specified isentropic efficiency characteristics
        if self.eta_s_char.is_set:
            mat_deriv += self.eta_s_char_deriv()

        ######################################################################
        # derivatives for specified pressure rise vs. flowrate characteristics
        if self.flow_char.is_set:
            mat_deriv += self.flow_char_deriv()

        return mat_deriv

    def eta_s_func(self):
        r"""
        Equation for given isentropic efficiency of a pump.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                0 = -\left( h_{out} - h_{in} \right) \cdot \eta_{s,c} +
                \left( h_{out,s} - h_{in} \right)
        """
        return (-(self.outl[0].h.val_SI - self.inl[0].h.val_SI) *
                self.eta_s.val + (self.h_os('post') - self.inl[0].h.val_SI))

    def eta_s_deriv(self):
        r"""
        Calculates the matrix of partial derivatives of the isentropic efficiency function.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))

        for i in range(2):
            deriv[0, i, 1] = self.numeric_deriv(self.eta_s_func, 'p', i)
            if i == 0:
                deriv[0, i, 2] = self.numeric_deriv(self.eta_s_func, 'h', i)
            else:
                deriv[0, i, 2] = -self.eta_s.val

        return deriv.tolist()

    def eta_s_char_func(self):
        r"""
        Equation for given isentropic efficiency characteristic of a pump.

        Returns
        -------
        res : ndarray
            Residual value of equation.

            .. math::

                0 = \left( h_{out} - h_{in} \right) \cdot \frac{\Delta h_{s,ref}}{\Delta h_{ref}}
                \cdot char\left( \frac{\dot{m}_{in} \cdot v_{in}}{\dot{m}_{in,ref} \cdot v_{in,ref}} \right) - \left( h_{out,s} - h_{in} \right)
        """
        # actual values
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()
        # design values
        i_d = self.inl[0].to_flow_design()
        o_d = self.outl[0].to_flow_design()

        expr = i[0] * v_mix_ph(i) / (i_d[0] * v_mix_ph(i_d))

        return (o[2] - i[2]) * self.dh_s_ref / (o_d[2] - i_d[2]) * self.eta_s_char.func.f_x(expr) - (self.h_os('post') - i[2])

    def eta_s_char_deriv(self):
        r"""
        Calculates the matrix of partial derivatives of the isentropic efficiency characteristic function.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))

        deriv[0, 0, 0] = self.numeric_deriv(self.eta_s_char_func, 'm', 0)
        for i in range(2):
            deriv[0, i, 1] = self.numeric_deriv(self.eta_s_char_func, 'p', i)
            deriv[0, i, 2] = self.numeric_deriv(self.eta_s_char_func, 'h', i)

        return deriv.tolist()

    def flow_char_func(self):
        r"""
        Equation for given flow characteristic of a pump.

        Returns
        -------
        res : ndarray
            Residual value of equation.

            .. math::

                0 = p_{out} - p_{in} - char\left( \dot{m}_{in} \cdot v_{in} \right)
        """
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()

        expr = i[0] * v_mix_ph(i)

        return o[1] - i[1] - self.flow_char.func.f_x(expr)

    def flow_char_deriv(self):
        r"""
        Calculates the matrix of partial derivatives of the flow characteristic of a pump.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))

        deriv[0, 0, 0] = self.numeric_deriv(self.flow_char_func, 'm', 0)
        deriv[0, 0, 2] = self.numeric_deriv(self.flow_char_func, 'h', 0)
        for i in range(2):
            deriv[0, i, 1] = self.numeric_deriv(self.flow_char_func, 'p', i)

        return deriv.tolist()

    def convergence_check(self, nw):
        r"""
        Performs a convergence check.

        Parameters
        ----------
        nw : tespy.networks.network
            The network object using this component.

        Note
        ----
        Manipulate enthalpies/pressure at inlet and outlet if not specified by user to match physically feasible constraints.
        """
        i, o = self.inl, self.outl

        if not o[0].p.val_set and o[0].p.val_SI < i[0].p.val_SI:
                o[0].p.val_SI = o[0].p.val_SI * 2
        if not i[0].p.val_set and o[0].p.val_SI < i[0].p.val_SI:
                i[0].p.val_SI = o[0].p.val_SI * 0.5

        if not o[0].h.val_set and o[0].h.val_SI < i[0].h.val_SI:
                o[0].h.val_SI = o[0].h.val_SI * 1.1
        if not i[0].h.val_set and o[0].h.val_SI < i[0].h.val_SI:
                i[0].h.val_SI = o[0].h.val_SI * 0.9

        if self.flow_char.is_set:
            expr = i[0].m.val_SI * v_mix_ph(i[0].to_flow())

            if expr > self.flow_char.func.x[-1] and not i[0].m.val_set:
                i[0].m.val_SI = self.flow_char.func.x[-1] / v_mix_ph(i[0].to_flow())
            elif expr < self.flow_char.func.x[1] and not i[0].m.val_set:
                i[0].m.val_SI = self.flow_char.func.x[0] / v_mix_ph(i[0].to_flow())
            else:
                pass

    def initialise_source(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's outlet.

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
                3 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            return 3e5

    def initialise_target(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's inlet.

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
                2.9 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            return 2.9e5

    def calc_parameters(self, mode):
        r"""
        Post and preprocessing parameter calculation/specification.

        Parameters
        ----------

        mode : str
            Pre- or postprocessing calculation.

        Note
        ----
        Generic preprocessing is handled by the base class. This method handles class specific pre- and postprocessing.
        """
        turbomachine.calc_parameters(self, mode)

        if mode == 'post':
            self.eta_s.val = ((self.h_os('post') - self.inl[0].h.val_SI) /
                              (self.outl[0].h.val_SI - self.inl[0].h.val_SI))
            if (self.eta_s.val > 1 or self.eta_s.val <= 0):
                msg = ('Invalid value for isentropic efficiency: '
                       'eta_s =' + str(self.eta_s.val) + ' at ' + self.label + '.')
                logging.error(msg)

            if self.eta_s_char.is_set:
                # get bound errors for isentropic efficiency characteristics
                i = self.inl[0].to_flow()
                i_d = self.inl[0].to_flow_design()
                expr = i[0] * v_mix_ph(i) / (i_d[0] * v_mix_ph(i_d))
                self.eta_s_char.func.get_bound_errors(expr)

            if self.flow_char.is_set:
                # get bound errors for flow characteristics
                i = self.inl[0].to_flow()
                o = self.outl[0].to_flow()
                expr = i[0] * v_mix_ph(i)
                self.flow_char.func.get_bound_errors(expr)

# %%


class compressor(turbomachine):
    r"""
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in} \right) - P\\
            0 = pr \cdot p_{in} - p_{out}

        - :func:`tespy.components.components.compressor.eta_s_func`

        **additional equations**

        - :func:`tespy.components.components.compressor.additional_equations`

    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/compressor.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    P : Sring/float/tespy.helpers.dc_cp
        Power, :math:`P/\text{W}`

    eta_s : Sring/float/tespy.helpers.dc_cp
        Isentropic efficiency, :math:`\eta_s/1`

    pr : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`

    eta_s_char : str/tespy.helpers.dc_cc
        Characteristic curve for isentropic efficiency, provide x and y values
        or use generic values (e. g. calculated from design case).

    char_map : str/tespy.helpers.dc_cm
        Characteristic map for pressure rise and isentropic efficiency vs. nondimensional mass flow,
        see tespy.components.characteristics.compressor for further information.

    igva : str/float/tespy.helpers.dc_cp
        Inlet guide vane angle, :math:`igva/^\circ`.

    Example
    -------
    >>> from tespy import cmp, con, nwk, hlp
    >>> import shutil
    >>> fluid_list = ['air']
    >>> nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ...     h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> si = cmp.sink('sink')
    >>> so = cmp.source('source')
    >>> cp = cmp.compressor('compressor')
    >>> inc = con.connection(so, 'out1', cp, 'in1')
    >>> outg = con.connection(cp, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> cp.set_attr(pr=10, eta_s=0.8, P=1e5, design=['eta_s'],
    ...     offdesign=['char_map'])
    >>> inc.set_attr(fluid={'air': 1}, p=1, T=20)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> cp.set_attr(P=9e4, igva='var')
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(cp.eta_s.val, 3)
    0.755
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'compressor'

    def attr(self):
        return {'P': dc_cp(), 'eta_s': dc_cp(), 'pr': dc_cp(),
                'igva': dc_cp(min_val=-45, max_val=45, d=1e-2, val=0),
                'Sirr': dc_cp(),
                'char_map': dc_cm(method='GENERIC'),
                'eta_s_char': dc_cc(param='m', method='GENERIC')}

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.fl_deriv = self.fluid_deriv()
        self.m_deriv = self.mass_flow_deriv()

        generate_char = False
        if self.char_map.func is None:
            generate_char = True
        elif (not np.array_equal(self.char_map.x, self.char_map.func.x) or
              not np.array_equal(self.char_map.y, self.char_map.func.y) or
              not np.array_equal(self.char_map.z1, self.char_map.func.z1) or
              not np.array_equal(self.char_map.z2, self.char_map.func.z2)):
            generate_char = True

        if generate_char:
            self.char_map.func = cmp_char.char_map(
                    x=self.char_map.x, y=self.char_map.y, z1=self.char_map.z1,
                    z2=self.char_map.z2, method=self.char_map.method, comp=self.component())
            self.char_map.x = self.char_map.func.x
            self.char_map.y = self.char_map.func.y
            self.char_map.z1 = self.char_map.func.z1
            self.char_map.z2 = self.char_map.func.z2

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for compressor.

        Equations

            **optional equations**

            - :func:`tespy.components.components.compressor.eta_s_char_func`
            - :func:`tespy.components.components.compressor.char_map_func`

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # equations for specified characteristic map
        if self.char_map.is_set:
            vec_res += self.char_map_func().tolist()

        ######################################################################
        # equation for specified isentropic efficiency characteristics
        if self.eta_s_char.is_set:
            vec_res += [self.eta_s_char_func()]

        return vec_res

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional equations.

        Returns
        -------
        mat_deriv : list
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for specified characteristic map
        if self.char_map.is_set:
            mat_deriv += self.char_map_deriv()

        ######################################################################
        # derivatives for specified isentropic efficiency characteristics
        if self.eta_s_char.is_set:
            mat_deriv += self.eta_s_char_deriv()

        return mat_deriv

    def eta_s_func(self):
        r"""
        Equation for given isentropic efficiency of a compressor.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                0 = -\left( h_{out} - h_{in} \right) \cdot \eta_{s,c} +
                \left( h_{out,s} - h_{in} \right)
        """
        return (-(self.outl[0].h.val_SI - self.inl[0].h.val_SI) *
                self.eta_s.val + (self.h_os('post') - self.inl[0].h.val_SI))

    def eta_s_deriv(self):
        r"""
        Calculates the matrix of partial derivatives of the isentropic efficiency function.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        mat_deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))

        for i in range(2):
            mat_deriv[0, i, 1] = self.numeric_deriv(self.eta_s_func, 'p', i)
            if i == 0:
                mat_deriv[0, i, 2] = self.numeric_deriv(self.eta_s_func, 'h', i)
            else:
                mat_deriv[0, i, 2] = -self.eta_s.val

        return mat_deriv.tolist()

    def eta_s_char_func(self):
        r"""
        Equation for given isentropic efficiency characteristic of a compressor.

        Returns
        -------
        res : ndarray
            Residual value of equation.

            .. math::

                0 = \left( h_{out} - h_{in} \right) \cdot \frac{\Delta h_{s,ref}}{\Delta h_{ref}}
                \cdot char\left( \dot{m}_{in} \cdot v_{in} \right) - \left( h_{out,s} - h_{in} \right)
        """
        # actual values
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()
        # design values
        i_d = self.inl[0].to_flow_design()
        o_d = self.outl[0].to_flow_design()

        expr = 1
        if self.eta_s_char.param == 'm':
            if not np.isnan(i_d[0]):
                expr = i[0] / i_d[0]
        elif self.eta_s_char.param == 'pr':
            if not np.isnan([i_d[1], o_d[1]]).any():
                expr = (o[1] * i_d[1]) / (i[1] * o_d[1])
        else:
            msg = 'Must provide a parameter for eta_s_char at component ' + self.label + '.'
            logging.error(msg)
            raise ValueError(msg)

        return (self.dh_s_ref / (o_d[2] - i_d[2]) * self.eta_s_char.func.f_x(expr) * (o[2] - i[2]) - (self.h_os('post') - i[2]))

    def eta_s_char_deriv(self):
        r"""
        Calculates the matrix of partial derivatives of the isentropic efficiency characteristic function.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        mat_deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))

        mat_deriv[0, 0, 1] = self.numeric_deriv(self.eta_s_char_func, 'p', 0)
        mat_deriv[0, 1, 1] = self.numeric_deriv(self.eta_s_char_func, 'p', 1)
        mat_deriv[0, 0, 2] = self.numeric_deriv(self.eta_s_char_func, 'h', 0)
        mat_deriv[0, 1, 2] = self.numeric_deriv(self.eta_s_char_func, 'h', 1)

        return mat_deriv.tolist()

    def char_map_func(self):
        r"""
        Equations for characteristic map of compressor.

        Parameters

            - X: speedline index (rotational speed is constant)
            - Y: nondimensional mass flow
            - Z1: pressure ratio equation
            - Z2: isentropic efficiency equation
            - igva: variable inlet guide vane angle (assumed 0° if not specified)

            .. math::

                X = \sqrt{\frac{T_\mathrm{1,ref}}{T_\mathrm{1}}}

                Y = \frac{\dot{m}_\mathrm{1} \cdot p_\mathrm{1,ref}}
                {\dot{m}_\mathrm{1,ref} \cdot p_\mathrm{1} \cdot X}

                Z1 = \frac{p_2 \cdot p_\mathrm{1,ref}}{p_1 \cdot p_\mathrm{2,ref}}-
                pr_{c}(char(m, igva))

                Z2 = \frac{\eta_\mathrm{s,c}}{\eta_\mathrm{s,c,ref}} -
                \eta_{s,c}(char(m, igva))

        Returns
        -------
        res : ndarray (Z1, Z2)
            Residual values of equations.
        """
        # actual values
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()
        # design values
        i_d = self.inl[0].to_flow_design()
        o_d = self.outl[0].to_flow_design()

        x = math.sqrt(T_mix_ph(i_d) / T_mix_ph(i))
        y = (i[0] * i_d[1]) / (i_d[0] * i[1] * x)

        pr, eta = self.char_map.func.get_pr_eta(x, y, self.igva.val)

        z1 = o[1] * i_d[1] / (i[1] * o_d[1]) - pr
        z2 = (self.h_os('post') - i[2]) / (o[2] - i[2]) / (self.dh_s_ref / (o_d[2] - i_d[2])) - eta

        return np.array([z1, z2])

    def char_map_deriv(self):
        r"""
        Calculates the matrix of partial derivatives of the compressor characteristic map function.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        m11 = self.numeric_deriv(self.char_map_func, 'm', 0)
        p11 = self.numeric_deriv(self.char_map_func, 'p', 0)
        h11 = self.numeric_deriv(self.char_map_func, 'h', 0)

        p21 = self.numeric_deriv(self.char_map_func, 'p', 1)
        h21 = self.numeric_deriv(self.char_map_func, 'h', 1)

        if self.igva.is_var:
            igva = self.numeric_deriv(self.char_map_func, 'igva', 1)

        deriv = np.zeros((2, 2 + self.num_vars, self.num_fl + 3))
        deriv[0, 0, 0] = m11[0]
        deriv[0, 0, 1] = p11[0]
        deriv[0, 0, 2] = h11[0]
        deriv[0, 1, 1] = p21[0]
        deriv[0, 1, 2] = h21[0]
        deriv[1, 0, 0] = m11[1]
        deriv[1, 0, 1] = p11[1]
        deriv[1, 0, 2] = h11[1]
        deriv[1, 1, 1] = p21[1]
        deriv[1, 1, 2] = h21[1]
        if self.igva.is_var:
            deriv[0, 2 + self.igva.var_pos, 0] = igva[0]
            deriv[1, 2 + self.igva.var_pos, 0] = igva[1]
        return deriv.tolist()

    def convergence_check(self, nw):
        r"""
        Performs a convergence check.

        Parameters
        ----------
        nw : tespy.networks.network
            The network object using this component.

        Note
        ----
        Manipulate enthalpies/pressure at inlet and outlet if not specified by user to match physically feasible constraints.
        """
        i, o = self.inl, self.outl

        if not o[0].p.val_set and o[0].p.val_SI < i[0].p.val_SI:
            o[0].p.val_SI = o[0].p.val_SI * 2
        if not i[0].p.val_set and o[0].p.val_SI < i[0].p.val_SI:
            i[0].p.val_SI = o[0].p.val_SI * 0.5

        if not o[0].h.val_set and o[0].h.val_SI < i[0].h.val_SI:
            o[0].h.val_SI = o[0].h.val_SI * 1.1
        if not i[0].h.val_set and o[0].h.val_SI < i[0].h.val_SI:
            i[0].h.val_SI = o[0].h.val_SI * 0.9

    def initialise_source(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's outlet.

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
                6 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            return 6e5

    def initialise_target(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's inlet.

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
                4 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            return 4e5

    def calc_parameters(self, mode):
        r"""
        Post and preprocessing parameter calculation/specification.

        Parameters
        ----------

        mode : str
            Pre- or postprocessing calculation.

        Note
        ----
        Generic preprocessing is handled by the base class. This method handles class specific pre- and postprocessing.
        """
        turbomachine.calc_parameters(self, mode)

        if mode == 'post':
            self.eta_s.val = ((self.h_os('post') - self.inl[0].h.val_SI) /
                              (self.outl[0].h.val_SI - self.inl[0].h.val_SI))
            if (self.eta_s.val > 1 or self.eta_s.val <= 0):
                msg = ('Invalid value for isentropic efficiency: '
                       'eta_s =' + str(self.eta_s.val) + ' at ' + self.label + '.')
                logging.error(msg)

            if self.char_map.is_set:
                # get bound errors for characteristic map
                i = self.inl[0].to_flow()
                i_d = self.inl[0].to_flow_design()
                x = math.sqrt(T_mix_ph(i_d)) / math.sqrt(T_mix_ph(i))
                y = (i[0] * i_d[1]) / (i_d[0] * i[1] * x)
                self.char_map.func.get_bound_errors(x, y, self.igva.val)

            if self.eta_s_char.is_set:
                # get bound errors for isentropic efficiency characteristics
                i = self.inl[0].to_flow()
                o = self.outl[0].to_flow()
                i_d = self.inl[0].to_flow_design()
                o_d = self.outl[0].to_flow_design()

                expr = 1
                if self.eta_s_char.param == 'm':
                    if not np.isnan(i_d[0]):
                        expr = i[0] / i_d[0]
                elif self.eta_s_char.param == 'pr':
                    if not np.isnan([i_d[1], o_d[1]]).any():
                        expr = (o[1] * i_d[1]) / (i[1] * o_d[1])

                self.eta_s_char.func.get_bound_errors(expr)

# %%


class turbine(turbomachine):
    r"""
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in} \right) - P\\
            0 = pr \cdot p_{in} - p_{out}

        - :func:`tespy.components.components.turbine.eta_s_func`

        **additional equations**

        - :func:`tespy.components.components.turbine.additional_equations`

    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/turbine.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    P : Sring/float/tespy.helpers.dc_cp
        Power, :math:`P/\text{W}`

    eta_s : Sring/float/tespy.helpers.dc_cp
        Isentropic efficiency, :math:`\eta_s/1`

    pr : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`

    eta_s_char : str/tespy.helpers.dc_cc
        Characteristic curve for isentropic efficiency, provide x and y values
        or use generic values (e. g. calculated from design case).

    cone : tespy.helpers.dc_cc
        Characteristics for stodolas cone law.

    Example
    -------
    >>> from tespy import cmp, con, nwk, hlp
    >>> import shutil
    >>> fluid_list = ['water']
    >>> nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ...     h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> si = cmp.sink('sink')
    >>> so = cmp.source('source')
    >>> t = cmp.turbine('turbine')
    >>> inc = con.connection(so, 'out1', t, 'in1')
    >>> outg = con.connection(t, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> t.set_attr(pr=0.02, eta_s=0.8, P=-1e5, design=['eta_s', 'pr'],
    ...     offdesign=['eta_s_char', 'cone'])
    >>> inc.set_attr(fluid={'water': 1}, T=600)
    >>> outg.set_attr(p=0.5)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> t.set_attr(P=-9e4)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(t.eta_s.val, 3)
    0.8
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'turbine'

    def attr(self):
        return {'P': dc_cp(), 'eta_s': dc_cp(), 'pr': dc_cp(),
                'Sirr': dc_cp(),
                'eta_s_char': dc_cc(method='GENERIC', param='m'),
                'cone': dc_cc(method='default')}

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for compressor.

        Equations

            **optional equations**

            - :func:`tespy.components.components.turbine.eta_s_char_func`
            - :func:`tespy.components.components.turbine.cone_func`

        Returns
        -------
        vec_res : list
            Vector of residual values.
        """
        vec_res = []

        ######################################################################
        # derivatives for specified isentropic efficiency characteristics
        if self.eta_s_char.is_set:
            vec_res += [self.eta_s_char_func()]

        ######################################################################
        # equation for specified cone law
        if self.cone.is_set:
            vec_res += [self.cone_func()]

        return vec_res

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional equations.

        Returns
        -------
        mat_deriv : list
            Matrix of partial derivatives.
        """
        mat_deriv = []

        ######################################################################
        # derivatives for specified isentropic efficiency characteristics
        if self.eta_s_char.is_set:
            mat_deriv += self.eta_s_char_deriv()

        ######################################################################
        # derivatives for specified cone law
        if self.cone.is_set:
            cone_deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            cone_deriv[0, 0, 0] = -1
            cone_deriv[0, 0, 1] = self.numeric_deriv(self.cone_func, 'p', 0)
            cone_deriv[0, 0, 2] = self.numeric_deriv(self.cone_func, 'h', 0)
            cone_deriv[0, 1, 2] = self.numeric_deriv(self.cone_func, 'p', 1)
            mat_deriv += cone_deriv.tolist()

        return mat_deriv

    def eta_s_func(self):
        r"""
        Equation for given isentropic efficiency of a compressor.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                0 = -\left( h_{out} - h_{in} \right) +
                \left( h_{out,s} - h_{in} \right) \cdot \eta_{s,e}
        """
        return (-(self.outl[0].h.val_SI - self.inl[0].h.val_SI) +
                (self.h_os('post') - self.inl[0].h.val_SI) *
                self.eta_s.val)

    def eta_s_deriv(self):
        r"""
        Calculates the matrix of partial derivatives of the isentropic efficiency function.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        mat_deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))

        for i in range(2):
            mat_deriv[0, i, 1] = self.numeric_deriv(self.eta_s_func, 'p', i)
            if i == 0:
                mat_deriv[0, i, 2] = self.numeric_deriv(self.eta_s_func, 'h', i)
            else:
                mat_deriv[0, i, 2] = -1

        return mat_deriv.tolist()

    def cone_func(self):
        r"""
        Equation for stodolas cone law.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                0 = \frac{\dot{m}_{in,ref} \cdot p_{in}}{p_{in,ref}} \cdot
                \sqrt{\frac{p_{in,ref} \cdot v_{in}}{p_{in} \cdot v_{in,ref}}}
                \cdot \sqrt{\frac{1 - \left(\frac{p_{out}}{p_{in}} \right)^{2}}
                {1 - \left(\frac{p_{out,ref}}{p_{in,ref}} \right)^{2}}} -
                \dot{m}_{in}
        """
        # actual values
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()
        # design values
        i_d = self.inl[0].to_flow_design()
        o_d = self.outl[0].to_flow_design()

        n = 1
        return (i_d[0] * i[1] / i_d[1] * math.sqrt(i_d[1] * v_mix_ph(i_d) / (i[1] * v_mix_ph(i))) *
                math.sqrt(abs((1 - (o[1] / i[1]) ** ((n + 1) / n)) / (1 - (o_d[1] / i_d[1]) ** ((n + 1) / n)))) - i[0])

    def eta_s_char_func(self):
        r"""
        Equation for given isentropic efficiency characteristic of a turbine.

        Returns
        -------
        res : ndarray
            Residual value of equation.

            .. math::

                0 = - \left( h_{out} - h_{in} \right) + \eta_{s,e,0} \cdot f\left(
                expr \right) \cdot \Delta h_{s}
        """
        # actual values
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()
        # design values
        i_d = self.inl[0].to_flow_design()
        o_d = self.outl[0].to_flow_design()

        if self.eta_s_char.param == 'dh_s':
            expr = math.sqrt(self.dh_s_ref / (self.h_os('post') - i[2]))
        elif self.eta_s_char.param == 'm':
            expr = i[0] / i_d[0]
        elif self.eta_s_char.param == 'v':
            expr = i[0] * v_mix_ph(i) / (i_d[0] * v_mix_ph(i_d))
        elif self.eta_s_char.param == 'pr':
            expr = (o[1] * i_d[1]) / (i[1] * o_d[1])
        else:
            msg = 'Please choose the parameter, you want to link the isentropic efficiency to.'
            logging.error(msg)
            raise ValueError(msg)

        return -(o[2] - i[2]) + (o_d[2] - i_d[2]) / self.dh_s_ref * self.eta_s_char.func.f_x(expr) * (self.h_os('post') - i[2])

    def eta_s_char_deriv(self):
        r"""
        Calculates the matrix of partial derivatives of the isentropic efficiency characteristic function.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        mat_deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))

        mat_deriv[0, 0, 0] = self.numeric_deriv(self.eta_s_char_func, 'm', 0)
        for i in range(2):
            mat_deriv[0, i, 1] = self.numeric_deriv(self.eta_s_char_func, 'p', i)
            mat_deriv[0, i, 2] = self.numeric_deriv(self.eta_s_char_func, 'h', i)

        return mat_deriv.tolist()

    def convergence_check(self, nw):
        r"""
        Performs a convergence check.

        Parameters
        ----------
        nw : tespy.networks.network
            The network object using this component.

        Note
        ----
        Manipulate enthalpies/pressure at inlet and outlet if not specified by user to match physically feasible constraints.
        """
        i, o = self.inl, self.outl

        if i[0].p.val_SI <= 1e5 and not i[0].p.val_set:
            i[0].p.val_SI = 1e5

        if i[0].p.val_SI <= o[0].p.val_SI and not o[0].p.val_set:
            o[0].p.val_SI = i[0].p.val_SI / 2

        if i[0].h.val_SI < 10e5 and not i[0].h.val_set:
            i[0].h.val_SI = 10e5

        if o[0].h.val_SI < 5e5 and not o[0].h.val_set:
            o[0].h.val_SI = 5e5

        if i[0].h.val_SI <= o[0].h.val_SI and not o[0].h.val_set:
            o[0].h.val_SI = i[0].h.val_SI * 0.75

    def initialise_source(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's outlet.

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
                5 \cdot 10^4 & \text{key = 'p'}\\
                1.5 \cdot 10^6 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 0.5e5
        elif key == 'h':
            return 1.5e6

    def initialise_target(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's inlet.

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
                2.5 \cdot 10^6 & \text{key = 'p'}\\
                2 \cdot 10^6 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 2.5e6
        elif key == 'h':
            return 2e6

    def calc_parameters(self, mode):
        r"""
        Post and preprocessing parameter calculation/specification.

        Parameters
        ----------

        mode : str
            Pre- or postprocessing calculation.

        Note
        ----
        Generic preprocessing is handled by the base class. This method handles class specific pre- and postprocessing.
        """
        turbomachine.calc_parameters(self, mode)

        if mode == 'post':
            self.eta_s.val = ((self.outl[0].h.val_SI - self.inl[0].h.val_SI) /
                              (self.h_os('post') - self.inl[0].h.val_SI))
            if (self.eta_s.val > 1 or self.eta_s.val <= 0):
                msg = ('Invalid value for isentropic efficiency: '
                       'eta_s =' + str(self.eta_s.val) + ' at ' + self.label + '.')
                logging.error(msg)

            if self.eta_s_char.is_set:
                # get bound errors for isentropic efficiency characteristics
                i = self.inl[0].to_flow()
                o = self.outl[0].to_flow()
                i_d = self.inl[0].to_flow_design()
                o_d = self.outl[0].to_flow_design()

                if self.eta_s_char.param == 'dh_s':
                    expr = math.sqrt(self.dh_s_ref / (self.h_os('post') - i[2]))
                elif self.eta_s_char.param == 'm':
                    expr = i[0] / i_d[0]
                elif self.eta_s_char.param == 'v':
                    expr = i[0] * v_mix_ph(i) / (i_d[0] * v_mix_ph(i_d))
                elif self.eta_s_char.param == 'pr':
                    expr = (o[1] * i_d[1]) / (i[1] * o_d[1])

                self.eta_s_char.func.get_bound_errors(expr)

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

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    num_in : float/tespy.helpers.dc_simple
        Number of inlets for this component.

    num_out : float/tespy.helpers.dc_simple
        Number of outlets for this component.

    Note
    ----
    - Node: Fluid composition and enthalpy at all **outgoing** connections (mass flow leaves the node) is result of mixture of the properties of the incoming connections (mass flow enters node).
      Incoming and outgoing connections can be a result of the calculation and are not identical to the inlets and outlets!
    - Splitter: Fluid composition and enthalpy at all outlets is the same as the inlet's properties.
    - Separator: Fluid composition is variable for all outlets, temperature at all outlets is the same as the inlet's temperature.
    - Merge: Fluid composition and enthalpy at outlet is result of mixture of the inlet's properties.

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
        Calculates vector vec_res with results of additional equations for this component.

        Equations

            **mandatroy equations**

            - :func:`tespy.components.components.node.fluid_func`

            .. math::

                0 = \sum_i \left(\dot{m}_{i} \cdot h_{i}\right) - h_{o} \cdot  \sum_i \dot{m}_{i}\\
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
        Calculates matrix of partial derivatives for given additional equations.

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
        deriv = np.zeros((len(self.outg), self.num_i + self.num_o, self.num_fl + 3))
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
        Calculates the vector of residual values for component's fluid balance equations.

        Returns
        -------
        vec_res : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = \sum_i \left(\dot{m}_{i} \cdot x_{i,j}\right) - x_{o,j} \cdot  \sum_i \dot{m}_{i}\\
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
        deriv = np.zeros((self.num_fl * num_o, self.num_i + self.num_o, 3 + self.num_fl))
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
        deriv = np.zeros((self.num_i + self.num_o - 1, self.num_i + self.num_o, self.num_fl + 3))

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
        Returns a starting value for pressure and enthalpy at component's outlet.

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
        Returns a starting value for pressure and enthalpy at component's inlet.

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

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    num_out : float/tespy.helpers.dc_simple
        Number of outlets for this component.

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
        Calculates vector vec_res with results of additional equations for this component.

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
        Calculates matrix of partial derivatives for given additional equations.

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
        deriv = np.zeros((self.num_fl * self.num_o, 1 + self.num_o, 3 + self.num_fl))
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

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    num_out : float/tespy.helpers.dc_simple
        Number of outlets for this component.

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
        Calculates vector vec_res with results of additional equations for this component.

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
            vec_res += [T_mix_ph(self.inl[0].to_flow()) -
                        T_mix_ph(o.to_flow())]

        return vec_res

    def additional_derivatives(self):
        r"""
        Calculates matrix of partial derivatives for given additional equations.

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

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    num_in : float/tespy.helpers.dc_simple
        Number of inlets for this component.

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
                'zero_flag': dc_cp(printout=False)}

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
        Calculates vector vec_res with results of additional equations for this component.

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
        Calculates matrix of partial derivatives for given additional equations.

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


class combustion_chamber(component):
    r"""
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.combustion_chamber.reaction_balance`
        - :func:`tespy.components.components.component.mass_flow_func`

        .. math::

            0 = p_{in,i} - p_{out} \;
            \forall i \in \mathrm{inlets}

        - :func:`tespy.components.components.combustion_chamber.energy_balance`

        **optional equations**

        - :func:`tespy.components.components.combustion_chamber.lambda_func`
        - :func:`tespy.components.components.combustion_chamber.ti_func`

    Available fuels

        - methane, ethane, propane, butane, hydrogen

    Inlets/Outlets

        - in1, in2
        - out1

    Image

        .. image:: _images/combustion_chamber.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    .. note::

        The fuel and the air components can be connected to either of the inlets.

    Parameters
    ----------
    label : str
        The label of the component.

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    fuel : str/tespy.helpers.dc_simple
        Fuel for the combustion chamber, see list of available fluids above.

    lamb : float/tespy.helpers.dc_cp
        Air to stoichiometric air ratio, :math:`\lambda/1`.

    ti : float/tespy.helpers.dc_cp
        Thermal input, (:math:`{LHV \cdot \dot{m}_f}`),
        :math:`ti/\text{W}`.

    Note
    ----

        For more information on the usage of the combustion chamber see the
        examples section on github or look for the combustion chamber tutorials
        at tespy.readthedocs.io

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> fluid_list = ['Ar', 'N2', 'O2', 'CO2', 'CH4', 'H2O']
    >>> nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ...     p_range=[0.5, 10], T_range=[10, 1200])
    >>> nw.set_printoptions(print_level='none')
    >>> amb = cmp.source('ambient')
    >>> sf = cmp.source('fuel')
    >>> fg = cmp.sink('flue gas outlet')
    >>> comb = cmp.combustion_chamber('combustion chamber')
    >>> amb_comb = con.connection(amb, 'out1', comb, 'in1')
    >>> sf_comb = con.connection(sf, 'out1', comb, 'in2')
    >>> comb_fg = con.connection(comb, 'out1', fg, 'in1')
    >>> nw.add_conns(sf_comb, amb_comb, comb_fg)
    >>> comb.set_attr(fuel='CH4', ti=50000)
    >>> amb_comb.set_attr(p=1, T=20,
    ...     fluid={'Ar': 0.0129, 'N2': 0.7553, 'H2O': 0,
    ...         'CH4': 0, 'CO2': 0.0004, 'O2': 0.2314})
    >>> sf_comb.set_attr(T=25,
    ...     fluid={'CO2': 0.04, 'Ar': 0, 'N2': 0,
    ...         'O2': 0, 'H2O': 0, 'CH4': 0.96})
    >>> comb_fg.set_attr(T=1200)
    >>> nw.solve('design')
    >>> round(comb.lamb.val, 3)
    2.009
    >>> round(comb.ti.val)
    50000.0
    """

    def component(self):
        return 'combustion chamber'

    def attr(self):
        return {'fuel': dc_simple(), 'lamb': dc_cp(), 'ti': dc_cp(),
                'S': dc_cp()}

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1']

    def fuels(self):
        return ['methane', 'ethane', 'propane', 'butane',
                'hydrogen']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.m_deriv = self.mass_flow_deriv()
        self.p_deriv = self.pressure_deriv()

        if not self.fuel.val_set:
            msg = 'Must specify fuel for component ' + self.label + '. Available fuels are: ' + str(self.fuels()) + '.'
            logging.error(msg)
            raise TESPyComponentError(msg)

        if (len([x for x in nw.fluids if x in [a.replace(' ', '') for a in
                 CP.get_aliases(self.fuel.val)]]) == 0):
            msg = 'The fuel you specified for component ' + self.label + ' does not match the fuels available within the network.'
            logging.error(msg)
            raise TESPyComponentError(msg)

        if (len([x for x in self.fuels() if x in [a.replace(' ', '') for a in
                 CP.get_aliases(self.fuel.val)]])) == 0:
            msg = 'The fuel you specified is not available for component ' + self.label + '. Available fuels are: ' + str(self.fuels()) + '.'
            logging.error(msg)
            raise TESPyComponentError(msg)

        self.fuel.val = [x for x in nw.fluids if x in [a.replace(' ', '') for a in CP.get_aliases(self.fuel.val)]][0]

        self.o2 = [x for x in nw.fluids if x in [a.replace(' ', '') for a in CP.get_aliases('O2')]][0]
        self.co2 = [x for x in nw.fluids if x in [a.replace(' ', '') for a in CP.get_aliases('CO2')]][0]
        self.h2o = [x for x in nw.fluids if x in [a.replace(' ', '') for a in CP.get_aliases('H2O')]][0]
        self.n2 = [x for x in nw.fluids if x in [a.replace(' ', '') for a in CP.get_aliases('N2')]][0]

        structure = fluid_structure(self.fuel.val)

        self.n = {}
        for el in ['C', 'H', 'O']:
            if el in structure.keys():
                self.n[el] = structure[el]
            else:
                self.n[el] = 0

        self.lhv = self.calc_lhv()
        msg = 'Combustion chamber fuel (' + self.fuel.val + ') LHV is ' + str(self.lhv) + ' for component ' + self.label + '.'
        logging.debug(msg)

    def calc_lhv(self):
        r"""
        calculates the lower heating value of the combustion chambers fuel.

        Returns
        -------
        val : float
            Lower heating value of the combustion chambers fuel.

            .. math::
                LHV = -\frac{\sum_i {\Delta H_f^0}_i -
                \sum_j {\Delta H_f^0}_j }
                {M_{fuel}}\\
                \forall i \in \text{reation products},\\
                \forall j \in \text{reation educts},\\
                \Delta H_f^0: \text{molar formation enthalpy}
        """
        hf = {}
        hf['hydrogen'] = 0
        hf['methane'] = -74.85
        hf['ethane'] = -84.68
        hf['propane'] = -103.8
        hf['butane'] = -124.51
        hf[self.o2] = 0
        hf[self.co2] = -393.5
        # water (gaseous)
        hf[self.h2o] = -241.8

        key = set(list(hf.keys())).intersection(
                set([a.replace(' ', '')
                     for a in CP.get_aliases(self.fuel.val)]))

        val = (-(self.n['H'] / 2 * hf[self.h2o] + self.n['C'] * hf[self.co2] -
                 ((self.n['C'] + self.n['H'] / 4) * hf[self.o2] +
                  hf[list(key)[0]])) /
               molar_masses[self.fuel.val] * 1000)

        return val

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
        # equations for fluids in reaction balance
        for fluid in self.inl[0].fluid.val.keys():
            vec_res += [self.reaction_balance(fluid)]

        ######################################################################
        # eqation for mass flow balance
        vec_res += self.mass_flow_func()

        ######################################################################
        # equations for pressure
        for i in self.inl:
            vec_res += [self.outl[0].p.val_SI - i.p.val_SI]

        ######################################################################
        # equation for energy balance
        vec_res += [self.energy_balance()]

        ######################################################################
        # equation for specified air to stoichiometric air ratio lamb
        if self.lamb.is_set:
            vec_res += [self.lambda_func()]

        ######################################################################
        # equation for speciified thermal input
        if self.ti.is_set:
            vec_res += [self.ti_func()]

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
        # derivatives for reaction balance
        j = 0
        deriv = np.zeros((self.num_fl, 3, self.num_fl + 3))
        for fluid in self.fluids:
            for i in range(3):
                deriv[j, i, 0] = self.rb_numeric_deriv('m', i, fluid)
                deriv[j, i, 3:] = self.rb_numeric_deriv('fluid', i, fluid)

            j += 1
        mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for mass balance equations
        mat_deriv += self.m_deriv

        ######################################################################
        # derivatives for pressure equations
        mat_deriv += self.p_deriv

        ######################################################################
        # derivatives for energy balance equations
        deriv = np.zeros((1, 3, self.num_fl + 3))
        for i in range(3):
            deriv[0, i, 0] = self.numeric_deriv(self.energy_balance, 'm', i)
            deriv[0, i, 1] = self.numeric_deriv(self.energy_balance, 'p', i)
            if i >= self.num_i:
                deriv[0, i, 2] = -(self.inl + self.outl)[i].m.val_SI
            else:
                deriv[0, i, 2] = (self.inl + self.outl)[i].m.val_SI
        mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for specified lamb
        if self.lamb.is_set:
            deriv = np.zeros((1, 3, self.num_fl + 3))
            for i in range(2):
                deriv[0, i, 0] = self.numeric_deriv(self.lambda_func, 'm', i)
                deriv[0, i, 3:] = self.numeric_deriv(self.lambda_func, 'fluid', i)
            mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for specified thermal input
        if self.ti.is_set:
            # stoichiometric combustion chamber
            if isinstance(self, combustion_chamber_stoich):
                pos = 3 + self.fluids.index('TESPy::' + self.fuel_alias.val)
                fuel = 'TESPy::' + self.fuel_alias.val
            # combustion chamber
            else:
                pos = 3 + self.fluids.index(self.fuel.val)
                fuel = self.fuel.val

            deriv = np.zeros((1, 3, self.num_fl + 3))
            for i in range(2):
                deriv[0, i, 0] = -self.inl[i].fluid.val[fuel]
                deriv[0, i, pos] = -self.inl[i].m.val_SI
            deriv[0, 2, 0] = self.outl[0].fluid.val[fuel]
            deriv[0, 2, pos] = self.outl[0].m.val_SI
            mat_deriv += (deriv * self.lhv).tolist()

        return np.asarray(mat_deriv)

    def pressure_deriv(self):
        r"""
        Calculates the partial derivatives for all pressure equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((2, 3, self.num_fl + 3))
        for k in range(2):
            deriv[k][2][1] = 1
            deriv[k][k][1] = -1
        return deriv.tolist()

    def reaction_balance(self, fluid):
        r"""
        Calculates the reaction balance for one fluid.

        - determine molar mass flows of fuel and oxygen
        - calculate excess fuel
        - calculate residual value of the fluids balance

        General equations

            .. math::

                \text{combustion chamber: } i \in [1,2], o \in [1]\\
                \text{cogeneration unit: } i \in [3,4], o \in [3]\\

                res = \sum_i \left(x_{fluid,i} \cdot \dot{m}_{i}\right) -
                \sum_j \left(x_{fluid,j} \cdot \dot{m}_{j}\right) \;
                \forall i, \; \forall j

                \dot{m}_{fluid,m} = \sum_i \frac{x_{fluid,i} \cdot \dot{m}_{i}}
                {M_{fluid}} \; \forall i

                \lambda = \frac{\dot{m}_{f,m}}{\dot{m}_{O_2,m} \cdot
                \left(n_{C,fuel} + 0.25 \cdot n_{H,fuel}\right)}

        Equation for fuel

            .. math::

                0 = res - \left(\dot{m}_{f,m} - \dot{m}_{f,exc,m}\right)
                \cdot M_{fuel}\\

                \dot{m}_{f,exc,m} = \begin{cases}
                0 & \lambda \geq 1\\
                \dot{m}_{f,m} - \frac{\dot{m}_{O_2,m}}
                {n_{C,fuel} + 0.25 \cdot n_{H,fuel}} & \lambda < 1
                \end{cases}

        Equation for oxygen

            .. math::

                0 = res - \begin{cases}
                -\frac{\dot{m}_{O_2,m} \cdot M_{O_2}}{\lambda} & \lambda \geq 1\\
                - \dot{m}_{O_2,m} \cdot M_{O_2} & \lambda < 1
                \end{cases}

        Equation for water

            .. math::

                0 = res + \left( \dot{m}_{f,m} - \dot{m}_{f,exc,m} \right)
                \cdot 0.5 \cdot n_{H,fuel} \cdot M_{H_2O}

        Equation for carbondioxide

            .. math::

                0 = res + \left( \dot{m}_{f,m} - \dot{m}_{f,exc,m} \right)
                \cdot n_{C,fuel} \cdot M_{CO_2}

        Equation for all other fluids

        .. math::

            0 = res

        Parameters
        ----------
        fluid : str
            The fluid to calculate the reation balance for.

        Returns
        -------
        res : float
            Residual value of equation.
        """

        if isinstance(self, cogeneration_unit):
            inl = self.inl[2:]
            outl = self.outl[2:]
        else:
            inl = self.inl
            outl = self.outl

        ######################################################################
        # molar mass flow for fuel and oxygen
        n_fuel = 0
        for i in inl:
            n_fuel += i.m.val_SI * i.fluid.val[self.fuel.val] / molar_masses[self.fuel.val]

        n_oxygen = 0
        for i in inl:
            n_oxygen += i.m.val_SI * i.fluid.val[self.o2] / molar_masses[self.o2]

        if n_fuel == 0:
            n_fuel = 1

        ######################################################################
        # calculate lambda if not set
        if not self.lamb.is_set:
            self.lamb.val = n_oxygen / (n_fuel * (self.n['C'] + self.n['H'] / 4))

        ######################################################################
        # calculate excess fuel if lambda is lower than 1
        n_fuel_exc = 0
        if self.lamb.val < 1:
            n_fuel_exc = n_fuel - n_oxygen / (self.n['C'] + self.n['H'] / 4)

        ######################################################################
        # equation for carbondioxide
        if fluid == self.co2:
            dm = (n_fuel - n_fuel_exc) * self.n['C'] * molar_masses[self.co2]

        ######################################################################
        # equation for water
        elif fluid == self.h2o:
            dm = (n_fuel - n_fuel_exc) * self.n['H'] / 2 * molar_masses[self.h2o]

        ######################################################################
        # equation for oxygen
        elif fluid == self.o2:
            if self.lamb.val < 1:
                dm = -n_oxygen * molar_masses[self.o2]
            else:
                dm = -n_oxygen / self.lamb.val * molar_masses[self.o2]

        ######################################################################
        # equation for fuel
        elif fluid == self.fuel.val:
            dm = -(n_fuel - n_fuel_exc) * molar_masses[self.fuel.val]

        ######################################################################
        # equation for other fluids
        else:
            dm = 0

        res = dm
        for i in inl:
            res += i.fluid.val[fluid] * i.m.val_SI
        for o in outl:
            res -= o.fluid.val[fluid] * o.m.val_SI
        return res

    def rb_numeric_deriv(self, dx, pos, fluid):
        r"""
        Calculates derivative of the reaction balance to dx at components inlet
        or outlet in position pos for the fluid fluid.

        Parameters
        ----------
        dx : str
            Partial derivative.

        pos : int
            Position of connection regarding to inlets and outlet of the component,
            logic: ['in1', 'in2', ..., 'out1', ...] -> 0, 1, ..., n, n + 1, ..., n + m

        fluid : str
            Fluid to calculate partial derivative of reaction balance for.

        Returns
        -------
        deriv : float/list
            Partial derivative(s) of the function :math:`f` to variable(s) :math:`x`.

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
                exp += self.reaction_balance(fluid)
                if (self.inl + self.outl)[pos].fluid.val[f] - 2 * df >= 0:
                    (self.inl + self.outl)[pos].fluid.val[f] -= 2 * df
                else:
                    (self.inl + self.outl)[pos].fluid.val[f] = 0
                exp -= self.reaction_balance(fluid)
                (self.inl + self.outl)[pos].fluid.val[f] = val

                deriv += [exp / (2 * (dm + dp + dh + df))]

        else:
            exp = 0
            (self.inl + self.outl)[pos].m.val_SI += dm
            (self.inl + self.outl)[pos].p.val_SI += dp
            (self.inl + self.outl)[pos].h.val_SI += dh
            exp += self.reaction_balance(fluid)

            (self.inl + self.outl)[pos].m.val_SI -= 2 * dm
            (self.inl + self.outl)[pos].p.val_SI -= 2 * dp
            (self.inl + self.outl)[pos].h.val_SI -= 2 * dh
            exp -= self.reaction_balance(fluid)
            deriv = exp / (2 * (dm + dp + dh + df))

            (self.inl + self.outl)[pos].m.val_SI += dm
            (self.inl + self.outl)[pos].p.val_SI += dp
            (self.inl + self.outl)[pos].h.val_SI += dh

        return deriv

    def energy_balance(self):
        r"""
        Calculates the energy balance of the adiabatic combustion chamber.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::
                res = \sum_i \dot{m}_{in,i} \cdot \left( h_{in,i} - h_{in,i,ref}
                \right) - \sum_j \dot{m}_{out,j} \cdot
                \left( h_{out,j} - h_{out,j,ref} \right) +\\
                H_{I,f} \cdot \left(\sum_i \dot{m}_{in,i} \cdot x_{f,i} -
                \sum_j \dot{m}_{out,j} \cdot x_{f,j} \right)\\
                \forall i \in \text{inlets}\; \forall j \in \text{outlets}

        Note
        ----
        The temperature for the reference state is set to 20 °C, thus
        the water may be liquid. In order to make sure, the state is
        referring to the lower heating value, the necessary enthalpy
        difference for evaporation is added. The stoichiometric combustion
        chamber uses a different reference, you will find it in the
        :func:`tespy.components.components.combustion_chamber_stoich.energy_balance`
        documentation.

        - Reference temperature: 293.15 K.
        - Reference pressure: 1 bar.
        """
        T_ref = 293.15
        p_ref = 1e5

        res = 0
        for i in self.inl:
            res += i.m.val_SI * (i.h.val_SI - h_mix_pT([0, p_ref, 0, i.fluid.val], T_ref))

        for o in self.outl:
            dh = 0
            n_h2o = o.fluid.val[self.h2o] / molar_masses[self.h2o]
            if n_h2o > 0:
                p = p_ref * n_h2o / molar_mass_flow(o.fluid.val)
                h = h_pT(p, T_ref, self.h2o)
                h_steam = CP.PropsSI('H', 'P', p, 'Q', 1, self.h2o)
                if h < h_steam:
                    dh = (h_steam - h) * o.fluid.val[self.h2o]

            res -= o.m.val_SI * (o.h.val_SI - h_mix_pT([0, p_ref, 0, o.fluid.val], T_ref) - dh)

        res += self.calc_ti()

        return res

    def lambda_func(self):
        r"""
        Calculates the residual for specified lambda.

        Returns
        -------
        val : float
            Residual value of function.

            .. math::

                \dot{m}_{fluid,m} = \sum_i \frac{x_{fluid,i} \cdot \dot{m}_{i}}
                {M_{fluid}}\\ \forall i \in inlets

                val = \frac{\dot{m}_{f,m}}{\dot{m}_{O_2,m} \cdot
                \left(n_{C,fuel} + 0.25 \cdot n_{H,fuel}\right)} - \lambda
        """
        if isinstance(self, cogeneration_unit):
            inl = self.inl[2:]
        else:
            inl = self.inl

        n_fuel = 0
        for i in inl:
            n_fuel += (i.m.val_SI * i.fluid.val[self.fuel.val] / molar_masses[self.fuel.val])

        n_oxygen = 0
        for i in inl:
            n_oxygen += (i.m.val_SI * i.fluid.val[self.o2] / molar_masses[self.o2])

        return (n_oxygen / (n_fuel * (self.n['C'] + self.n['H'] / 4)) - self.lamb.val)

    def ti_func(self):
        r"""
        Calculates the residual for specified thermal input.

        Returns
        -------
        val : float
            Residual value of function.

            .. math::

                val = ti - \dot{m}_f \cdot LHV
        """
        return self.ti.val - self.calc_ti()

    def calc_ti(self):
        r"""
        Calculates the thermal input of the combustion chamber.

        Returns
        -------
        ti : float
            Thermal input.

            .. math::

                ti = LHV \cdot \left[\sum_i \left(\dot{m}_{in,i} \cdot x_{f,i}
                \right) - \dot{m}_{out,1} \cdot x_{f,1} \right]
                \; \forall i \in [1,2]
        """
        m = 0
        for i in self.inl:
            m += i.m.val_SI * i.fluid.val[self.fuel.val]

        for o in self.outl:
            m -= o.m.val_SI * o.fluid.val[self.fuel.val]

        return m * self.lhv

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

                val = LHV \cdot \dot{m}_{f} \cdot f_{char}\left( \frac{\dot{m}_{f}}{\dot{m}_{f,ref}}\right)
        """
        val = self.calc_ti()
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
        deriv = np.zeros((1, 3, len(self.inl[0].fluid.val) + 3))
        for i in range(2):
            deriv[0, i, 0] = self.numeric_deriv(self.bus_func, 'm', i, bus=bus)
            deriv[0, i, 3:] = self.numeric_deriv(self.bus_func, 'fluid', i, bus=bus)

        deriv[0, 2, 0] = self.numeric_deriv(self.bus_func, 'm', 2, bus=bus)
        deriv[0, 2, 3:] = self.numeric_deriv(self.bus_func, 'fluid', 2, bus=bus)
        return deriv

    def initialise_fluids(self, nw):
        r"""
        Calculates reaction balance with given lambda of 3 for good generic starting values at the component's outlet.

        Parameters
        ----------
        nw : tespy.networks.network
            Network using this component object.
        """
        N_2 = 0.7655
        O_2 = 0.2345

        n_fuel = 1
        lamb = 3
        m_co2 = n_fuel * self.n['C'] * molar_masses[self.co2]
        m_h2o = n_fuel * self.n['H'] / 2 * molar_masses[self.h2o]

        n_o2 = (m_co2 / molar_masses[self.co2] + 0.5 * m_h2o / molar_masses[self.h2o]) * lamb

        m_air = n_o2 * molar_masses[self.o2] / O_2
        m_fuel = n_fuel * molar_masses[self.fuel.val]
        m_fg = m_air + m_fuel

        m_o2 = n_o2 * molar_masses[self.o2] * (1 - 1 / lamb)
        m_n2 = N_2 * m_air

        fg = {
            self.n2: m_n2 / m_fg,
            self.co2: m_co2 / m_fg,
            self.o2: m_o2 / m_fg,
            self.h2o: m_h2o / m_fg
        }

        for o in self.outl:
            for fluid, x in o.fluid.val.items():
                if not o.fluid.val_set[fluid] and fluid in fg.keys():
                    o.fluid.val[fluid] = fg[fluid]

    def convergence_check(self, nw):
        r"""
        Performs a convergence check.

        Parameters
        ----------
        nw : tespy.networks.network
            The network object using this component.

        Note
        ----
        Manipulate enthalpies/pressure at inlet and outlet if not specified by user to match physically feasible constraints,
        keep fluid composition within feasible range and then propagates it towards the outlet.
        """
        if isinstance(self, cogeneration_unit):
            inl = self.inl[2:]
            outl = self.outl[2:]
        else:
            inl = self.inl
            outl = self.outl

        m = 0
        for i in inl:
            if i.m.val_SI < 0 and not i.m.val_set:
                i.m.val_SI = 0.01
            m += i.m.val_SI

        ######################################################################
        # check fluid composition
        for o in outl:
            fluids = [f for f in o.fluid.val.keys() if not o.fluid.val_set[f]]
            for f in fluids:
                if f not in [self.o2, self.co2, self.h2o, self.fuel.val]:
                    m_f = 0
                    for i in inl:
                        m_f += i.fluid.val[f] * i.m.val_SI

                    if abs(o.fluid.val[f] - m_f / m) > 0.03:
                        o.fluid.val[f] = m_f / m

                elif f == self.o2:
                    if o.fluid.val[f] > 0.25:
                        o.fluid.val[f] = 0.2
                    if o.fluid.val[f] < 0.05:
                        o.fluid.val[f] = 0.05

                elif f == self.co2:
                    if o.fluid.val[f] > 0.075:
                        o.fluid.val[f] = 0.075
                    if o.fluid.val[f] < 0.02:
                        o.fluid.val[f] = 0.02

                elif f == self.h2o:
                    if o.fluid.val[f] > 0.075:
                        o.fluid.val[f] = 0.075
                    if o.fluid.val[f] < 0.02:
                        o.fluid.val[f] = 0.02

                elif f == self.fuel.val:
                    if o.fluid.val[f] > 0:
                        o.fluid.val[f] = 0

        ######################################################################
        # flue gas propagation
        for o in outl:
            if o.m.val_SI < 0 and not o.m.val_set:
                o.m.val_SI = 10
            nw.init_target(o, o.t)

            if o.h.val_SI < 7.5e5 and not o.h.val_set:
                o.h.val_SI = 1e6

        ######################################################################
        # additional checks for performance improvement
        if self.lamb.val < 2 and not self.lamb.is_set:
            for i in inl:
                fuel_set = True
                if i.fluid.val[self.fuel.val] > 0.75 and not i.m.val_set:
                    fuel_set = False
                if i.fluid.val[self.fuel.val] < 0.75:
                    air_tmp = i.m.val_SI

            if not fuel_set:
                for i in inl:
                    if i.fluid.val[self.fuel.val] > 0.75:
                        i.m.val_SI = air_tmp / 25

    def initialise_source(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's outlet.

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
                10^6 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            return 10e5

    def initialise_target(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's inlet.

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
                5 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            return 5e5

    def calc_parameters(self, mode):
        r"""
        Post and preprocessing parameter calculation/specification.

        Parameters
        ----------

        mode : str
            Pre- or postprocessing calculation.

        Note
        ----
        Generic preprocessing is handled by the base class. This method handles class specific pre- and postprocessing.
        """
        component.calc_parameters(self, mode)

        if mode == 'post':
            self.ti.val = self.calc_ti()

            n_fuel = 0
            for i in self.inl:
                n_fuel += i.m.val_SI * i.fluid.val[self.fuel.val] / molar_masses[self.fuel.val]

            n_oxygen = 0
            for i in self.inl:
                n_oxygen += i.m.val_SI * i.fluid.val[self.o2] / molar_masses[self.o2]

            self.lamb.val = n_oxygen / (n_fuel * (self.n['C'] + self.n['H'] / 4))

            val = 0
            T_ref = 293.15
            p_ref = 1e5

            for i in self.inl:
                val += i.m.val_SI * (s_mix_ph(i.to_flow()) - s_mix_pT([0, p_ref, 0, i.fluid.val], T_ref))

            for o in self.outl:
                dS = 0
                n_h2o = o.fluid.val[self.h2o] / molar_masses[self.h2o]
                if n_h2o > 0:
                    p = p_ref * n_h2o / molar_mass_flow(o.fluid.val)
                    S = s_pT(p, T_ref, self.h2o)
                    S_steam = CP.PropsSI('H', 'P', p, 'Q', 1, self.h2o)
                    if S < S_steam:
                        dS = (S_steam - S) * o.fluid.val[self.h2o]
                val -= o.m.val_SI * (s_mix_ph(o.to_flow()) - s_mix_pT([0, p_ref, 0, o.fluid.val], T_ref) - dS)

            self.S.val = val

# %%


class combustion_chamber_stoich(combustion_chamber):
    r"""
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.combustion_chamber_stoich.reaction_balance`
        - :func:`tespy.components.components.component.mass_flow_func`

        .. math::

            0 = p_{in,i} - p_{out} \;
            \forall i \in \mathrm{inlets}

        - :func:`tespy.components.components.combustion_chamber_stoich.energy_balance`

        **optional equations**

        - :func:`tespy.components.components.combustion_chamber_stoich.lambda_func`
        - :func:`tespy.components.components.combustion_chamber_stoich.ti_func`

    Available fuels

        - methane, ethane, propane, butane, hydrogen

    Inlets/Outlets

        - in1, in2
        - out1

    Image

        .. image:: _images/combustion_chamber.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    .. note::

        The fuel and the air components can be connected to either of the inlets.

    Parameters
    ----------
    label : str
        The label of the component.

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    fuel : dict
        Fuel composition, e. g. :code:`{'CH4': 0.96, 'CO2': 0.04}`.

    fuel_alias : str
        Alias for the fuel, name of fuel for usage in network will be TESPy::fuel_alias.

    air : dict
        Fresh air composition, e. g. :code:`{'N2': 0.76, 'O2': 0.23, 'Ar': 0.01}`.

    air_alias : str
        Alias for the fresh air, name of air for usage in network will be TESPy::air_alias.

    path : str
        Path to existing fluid property table.

    lamb : float/tespy.helpers.dc_cp
        Air to stoichiometric air ratio, :math:`\lambda/1`.

    ti : float/tespy.helpers.dc_cp
        Thermal input, (:math:`{LHV \cdot \dot{m}_f}`),
        :math:`ti/\text{W}`.

    Note
    ----
    This combustion chamber uses fresh air and its fuel as the only
    reactive gas components. Therefore note the following restrictions. You
    are to

    - specify the fluid composition of the fresh air,
    - fully define the fuel's fluid components,
    - provide the aliases of the fresh air and the fuel and
    - make sure, both of the aliases are part of the network fluid vector.

    If you choose 'Air' or 'air' as alias for the fresh air, TESPy will use
    the fluid properties from CoolProp's air. Else, a custom fluid
    'TESPy::yourairalias' will be created.

    The name of the flue gas will be: 'TESPy::yourfuelalias_fg'. It is also
    possible to use fluid mixtures for the fuel, e. g.
    :code:`fuel={CH4: 0.9, 'CO2': 0.1}`. If you specify a fluid mixture for
    the fuel, TESPy will automatically create a custom fluid called:
    'TESPy::yourfuelalias'. For more information see the examples section
    or look for the combustion chamber tutorials at tespy.readthedocs.io.

    Example
    -------
    >>> from tespy import con, cmp, nwk, hlp
    >>> import shutil
    >>> fluid_list = ['TESPy::myAir', 'TESPy::myFuel', 'TESPy::myFuel_fg']
    >>> nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ...     p_range=[0.001, 10], T_range=[10, 2000])
    >>> amb = cmp.source('ambient')
    >>> sf = cmp.source('fuel')
    >>> fg = cmp.sink('flue gas outlet')
    >>> comb = cmp.combustion_chamber_stoich('stoichiometric combustion chamber')
    >>> amb_comb = con.connection(amb, 'out1', comb, 'in1')
    >>> sf_comb = con.connection(sf, 'out1', comb, 'in2')
    >>> comb_fg = con.connection(comb, 'out1', fg, 'in1')
    >>> nw.add_conns(sf_comb, amb_comb, comb_fg)
    >>> comb.set_attr(fuel={'CH4': 0.96, 'CO2': 0.04},
    ...     air={'Ar': 0.0129, 'N2': 0.7553, 'H2O': 0,
    ...     'CH4': 0, 'CO2': 0.0004, 'O2': 0.2314},
    ...     fuel_alias='myFuel', air_alias='myAir',
    ...     lamb=3, ti=20000)
    >>> amb_comb.set_attr(T=20, p=1,
    ...     fluid={'TESPy::myAir': 1, 'TESPy::myFuel': 0,
    ...     'TESPy::myFuel_fg': 0})
    >>> sf_comb.set_attr(T=25,
    ...     fluid={'TESPy::myAir': 0, 'TESPy::myFuel': 1,
    ...     'TESPy::myFuel_fg': 0})
    >>> nw.set_printoptions(iterinfo=False)
    >>> nw.solve('design')
    >>> round(comb_fg.T.val, 1)
    860.2
    >>> comb.set_attr(path='./LUT')
    >>> nw.solve('design')
    >>> round(comb_fg.T.val, 1)
    860.2
    >>> shutil.rmtree('./LUT', ignore_errors=True)
    """

    def component(self):
        return 'combustion chamber stoichiometric flue gas'

    def attr(self):
        return {'fuel': dc_simple(),
                'fuel_alias': dc_simple(),
                'air': dc_simple(),
                'air_alias': dc_simple(),
                'path': dc_simple(),
                'lamb': dc_cp(), 'ti': dc_cp(), 'S': dc_cp()}

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1']

    def fuels(self):
        return ['methane', 'ethane', 'propane', 'butane',
                'hydrogen']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.m_deriv = self.mass_flow_deriv()
        self.p_deriv = self.pressure_deriv()

        if not self.fuel.val_set or not isinstance(self.fuel.val, dict):
            msg = 'Must specify fuel composition for combustion chamber.'
            logging.error(msg)
            raise TESPyComponentError(msg)

        if not self.fuel_alias.val_set:
            msg = 'Must specify fuel alias for combustion chamber.'
            logging.error(msg)
            raise TESPyComponentError(msg)
        if 'TESPy::' in self.fuel_alias.val:
            msg = 'Can not use \'TESPy::\' at this point.'
            logging.error(msg)
            raise TESPyComponentError(msg)

        if not self.air.val_set or not isinstance(self.air.val, dict):
            msg = 'Must specify air composition for combustion chamber.'
            logging.error(msg)
            raise TESPyComponentError(msg)

        if not self.air_alias.val_set:
            msg = 'Must specify air alias for combustion chamber.'
            logging.error(msg)
            raise TESPyComponentError(msg)
        if 'TESPy::' in self.air_alias.val:
            msg = 'Can not use \'TESPy::\' at this point.'
            logging.error(msg)
            raise TESPyComponentError(msg)

        # adjust the names for required fluids according to naming in the network
        # air
        for f in self.air.val.keys():
            alias = [x for x in nw.fluids if x in [a.replace(' ', '') for a in CP.get_aliases(f)]]
            if len(alias) > 0:
                self.air.val[alias[0]] = self.air.val.pop(f)

        # fuel
        for f in self.fuel.val.keys():
            alias = [x for x in self.air.val.keys() if x in [a.replace(' ', '') for a in CP.get_aliases(f)]]
            if len(alias) > 0:
                self.fuel.val[alias[0]] = self.fuel.val.pop(f)

        # list of all fluids of air and fuel
        fluids = list(self.air.val.keys()) + list(self.fuel.val.keys())

        # oxygen
        alias = [x for x in fluids if x in [a.replace(' ', '') for a in CP.get_aliases('O2')]]
        if len(alias) == 0:
            msg = 'Oxygen missing in input fluids.'
            logging.error(msg)
            raise TESPyComponentError(msg)
        else:
            self.o2 = alias[0]

        # carbondioxide
        self.co2 = [x for x in nw.fluids if x in [a.replace(' ', '') for a in CP.get_aliases('CO2')]]
        if len(self.co2) == 0:
            self.co2 = 'CO2'
        else:
            self.co2 = self.co2[0]

        # water
        self.h2o = [x for x in nw.fluids if x in [a.replace(' ', '') for a in CP.get_aliases('H2O')]]
        if len(self.h2o) == 0:
            self.h2o = 'H2O'
        else:
            self.h2o = self.h2o[0]

        for f in fluids:
            memorise.heos[f] = CP.AbstractState('HEOS', f)

        # calculate lower heating value of specified fuel
        self.lhv = self.calc_lhv()
        msg = 'Combustion chamber fuel (' + self.fuel_alias.val + ') LHV is ' + str(self.lhv) + ' for component ' + self.label + '.'
        logging.debug(msg)
        # generate fluid properties for stoichiometric flue gas
        self.stoich_flue_gas(nw)

    def calc_lhv(self):
        r"""
        calculates the lower heating value of the combustion chambers fuel.

        Returns
        -------
        val : float
            Lower heating value of the combustion chambers fuel.

            .. math::

                LHV = \sum_{fuels} \left(-\frac{\sum_i {\Delta H_f^0}_i -
                \sum_j {\Delta H_f^0}_j }
                {M_{fuel}} \cdot x_{fuel} \right)\\
                \forall i \in \text{reation products},\\
                \forall j \in \text{reation educts},\\
                \forall fuel \in \text{fuels},\\
                \Delta H_f^0: \text{molar formation enthalpy},\\
                x_{fuel}: \text{mass fraction of fuel in fuel mixture}
        """
        hf = {}
        hf['hydrogen'] = 0
        hf['methane'] = -74.85
        hf['ethane'] = -84.68
        hf['propane'] = -103.8
        hf['butane'] = -124.51
        hf['O2'] = 0
        hf['CO2'] = -393.5
        # water (gaseous)
        hf['H2O'] = -241.8

        lhv = 0

        for f, x in self.fuel.val.items():
            molar_masses[f] = CP.PropsSI('M', f)
            fl = set(list(hf.keys())).intersection(
                    set([a.replace(' ', '') for a in CP.get_aliases(f)]))
            if len(fl) == 0:
                continue

            if list(fl)[0] in self.fuels():
                structure = fluid_structure(f)

                n = {}
                for el in ['C', 'H', 'O']:
                    if el in structure.keys():
                        n[el] = structure[el]
                    else:
                        n[el] = 0

                lhv += (-(n['H'] / 2 * hf['H2O'] + n['C'] * hf['CO2'] -
                          ((n['C'] + n['H'] / 4) * hf['O2'] +
                           hf[list(fl)[0]])) / molar_masses[f] * 1000) * x

        return lhv

    def stoich_flue_gas(self, nw):
        r"""
        Calculates the fluid composition of the stoichiometric flue gas and creates a custom fluid.

        - uses one mole of fuel as reference quantity and :math:`\lambda=1`
          for stoichiometric flue gas calculation (no oxygen in flue gas)
        - calculate molar quantities of (reactive) fuel components to determine
          water and carbondioxide mass fraction in flue gas
        - calculate required molar quantity for oxygen and required fresh
          air mass
        - calculate residual mass fractions for non reactive components of
          fresh air in the flue gas
        - calculate flue gas fluid composition
        - generate custom fluid porperties



        Reactive components in fuel

            .. math::

                m_{fuel} = \frac{1}{M_{fuel}}\\
                m_{CO_2} = \sum_{i} \frac{x_{i} \cdot m_{fuel} \cdot num_{C,i}
                \cdot M_{CO_{2}}}{M_{i}}\\
                m_{H_{2}O} = \sum_{i} \frac{x_{i} \cdot m_{fuel} \cdot num_{H,i}
                \cdot M_{H_{2}O}}{2 \cdot M_{i}}\\
                \forall i \in \text{fuels in fuel vector},\\
                num = \text{number of atoms in molecule}

        Other components of fuel vector

            .. math::

                m_{fg,j} = x_{j} \cdot m_{fuel}\\
                \forall j \in \text{non fuels in fuel vecotr, e. g. } CO_2,\\
                m_{fg,j} = \text{mass of fluid component j in flue gas}

        Non-reactive components in air

            .. math::

                n_{O_2} = \left( \frac{m_{CO_2}}{M_{CO_2}} +
                \frac{m_{H_{2}O}}{0,5 \cdot M_{H_{2}O}} \right) \cdot \lambda,\\
                n_{O_2} = \text{mol of oxygen required}\\
                m_{air} = \frac{n_{O_2} \cdot M_{O_2}}{x_{O_{2}, air}},\\
                m_{air} = \text{required total air mass}\\
                m_{fg,j} = x_{j, air} \cdot m_{air}\\
                m_{fg, O_2} = 0,\\
                m_{fg,j} = \text{mass of fluid component j in flue gas}

        Flue gas composition

            .. math::

                x_{fg,j} = \frac{m_{fg, j}}{m_{air} + m_{fuel}}

        Parameters
        ----------
        nw : tespy.networks.network
            TESPy network to generate stoichiometric flue gas for.
        """
        lamb = 1
        n_fuel = 1
        m_fuel = 1 / molar_mass_flow(self.fuel.val) * n_fuel
        m_fuel_fg = m_fuel
        m_co2 = 0
        m_h2o = 0
        molar_masses[self.h2o] = CP.PropsSI('M', self.h2o)
        molar_masses[self.co2] = CP.PropsSI('M', self.co2)
        molar_masses[self.o2] = CP.PropsSI('M', self.o2)

        self.fg = {}
        self.fg[self.co2] = 0
        self.fg[self.h2o] = 0

        for f, x in self.fuel.val.items():
            fl = set(list(self.fuels())).intersection(
                    set([a.replace(' ', '') for a in CP.get_aliases(f)]))

            if len(fl) == 0:
                if f in self.fg.keys():
                    self.fg[f] += x * m_fuel
                else:
                    self.fg[f] = x * m_fuel
            else:
                n_fluid = x * m_fuel / molar_masses[f]
                m_fuel_fg -= n_fluid * molar_masses[f]
                structure = fluid_structure(f)
                n = {}
                for el in ['C', 'H', 'O']:
                    if el in structure.keys():
                        n[el] = structure[el]
                    else:
                        n[el] = 0

                m_co2 += n_fluid * n['C'] * molar_masses[self.co2]
                m_h2o += n_fluid * n['H'] / 2 * molar_masses[self.h2o]

        self.fg[self.co2] += m_co2
        self.fg[self.h2o] += m_h2o

        n_o2 = (m_co2 / molar_masses[self.co2] + 0.5 * m_h2o / molar_masses[self.h2o]) * lamb
        m_air = n_o2 * molar_masses[self.o2] / self.air.val[self.o2]

        self.air_min = m_air / m_fuel

        for f, x in self.air.val.items():
            if f != self.o2:
                if f in self.fg.keys():
                    self.fg[f] += m_air * x
                else:
                    self.fg[f] = m_air * x

        m_fg = m_fuel + m_air

        for f in self.fg.keys():
            self.fg[f] /= m_fg

        if not self.path.val_set:
            self.path.val = None
        tespy_fluid(self.fuel_alias.val, self.fuel.val, [1000, nw.p_range_SI[1]], nw.T_range_SI, path=self.path.val)
        tespy_fluid(self.fuel_alias.val + '_fg', self.fg, [1000, nw.p_range_SI[1]], nw.T_range_SI, path=self.path.val)
        msg = 'Generated lookup table for ' + self.fuel_alias.val + ' and for stoichiometric flue gas at stoichiometric combustion chamber ' + self.label + '.'
        logging.debug(msg)

        if self.air_alias.val not in ['Air', 'air']:
            tespy_fluid(self.air_alias.val, self.air.val, [1000, nw.p_range_SI[1]], nw.T_range_SI, path=self.path.val)
            msg = 'Generated lookup table for ' + self.air_alias.val + ' at stoichiometric combustion chamber ' + self.label + '.'
        else:
            msg = 'Using CoolProp air at stoichiometric combustion chamber ' + self.label + '.'
        logging.debug(msg)

    def reaction_balance(self, fluid):
        r"""
        Calculates the reaction balance for one fluid.

        - determine molar mass flows of fuel and oxygen
        - calculate excess fuel
        - calculate residual value of the fluids balance

        General equations

            .. math::

                res = \sum_i \left(x_{fluid,i} \cdot \dot{m}_{i}\right) -
                \sum_j \left(x_{fluid,j} \cdot \dot{m}_{j}\right)\\
                \forall i \in [1,2], \; \forall j \in [1]

                \dot{m}_{air,min} = \dot{m}_{fuel} \cdot air_{min}

                \lambda = \frac{\dot{m}_{air}}{\dot{m}_{air,min}}

        Equation for fuel

            .. math::

                0 = res - \left(\dot{m}_{f} - \dot{m}_{f,exc}\right)

                \dot{m}_{f,exc} = \begin{cases}
                0 & \lambda \geq 1\\
                \dot{m}_{f} - \frac{\dot{m}_{air}}
                {\lambda \cdot air_{min}} & \lambda < 1
                \end{cases}

        Equation for air

            .. math::

                0 = res - \begin{cases}
                -\dot{m}_{air,min} & \lambda \geq 1\\
                -\dot{m}_{air} & \lambda < 1
                \end{cases}

        Equation for stoichiometric flue gas

            .. math::

                0 = res + \dot{m}_{air,min} + \dot{m}_{f}

        Equation for all other fluids

        .. math::

            0 = res

        Parameters
        ----------
        fluid : str
            The fluid to calculate the reation balance for.

        Returns
        -------
        res : float
            Residual value of equation.
        """
        if self.air_alias.val in ['air', 'Air']:
            air = self.air_alias.val
        else:
            air = 'TESPy::' + self.air_alias.val
        fuel = 'TESPy::' + self.fuel_alias.val
        flue_gas = 'TESPy::' + self.fuel_alias.val + '_fg'

        ######################################################################
        # calculate fuel and air mass flow
        m_fuel = 0
        for i in self.inl:
            m_fuel += i.m.val_SI * i.fluid.val[fuel]

        m_air = 0
        for i in self.inl:
            m_air += i.m.val_SI * i.fluid.val[air]

        m_air_min = self.air_min * m_fuel

        ######################################################################
        # calculate lambda if not specified
        if not self.lamb.is_set:
            self.lamb.val = m_air / (self.air_min * m_fuel)

        ######################################################################
        # calculate excess fuel if lambda is smaller than 1
        m_fuel_exc = 0
        if self.lamb.val < 1:
            m_fuel_exc = m_fuel - m_air / (self.lamb.val * self.air_min)

        ######################################################################
        # equation for air
        if fluid == air:
            if self.lamb.val >= 1:
                dm = -m_air_min
            else:
                dm = -m_air

        ######################################################################
        # equation for fuel
        elif fluid == fuel:
            dm = -(m_fuel - m_fuel_exc)

        ######################################################################
        # equation for flue gas
        elif fluid == flue_gas:
            dm = m_air_min + m_fuel

        ######################################################################
        # equation for other components
        else:
            dm = 0

        res = dm
        for i in self.inl:
            res += i.fluid.val[fluid] * i.m.val_SI
        for o in self.outl:
            res -= o.fluid.val[fluid] * o.m.val_SI
        return res

    def energy_balance(self):
        r"""
        Calculates the energy balance of the adiabatic combustion chamber.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::
                res = \sum_i \dot{m}_{in,i} \cdot \left( h_{in,i} - h_{in,i,ref}
                \right) - \sum_j \dot{m}_{out,j} \cdot
                \left( h_{out,j} - h_{out,j,ref} \right) +
                H_{I,f} \cdot \left(\sum_i \dot{m}_{in,i} \cdot x_{f,i} -
                \sum_j \dot{m}_{out,j} \cdot x_{f,j} \right)
                \; \forall i \in \text{inlets}\; \forall j \in \text{outlets}

        Note
        ----
        The temperature for the reference state is set to 100 °C, as the
        custom fluid properties are inacurate at the dew-point of water in
        the flue gas!

        - Reference temperature: 373.15 K.
        - Reference pressure: 1 bar.
        """
        T_ref = 373.15
        p_ref = 1e5

        res = 0
        for i in self.inl:
            res += i.m.val_SI * (i.h.val_SI - h_mix_pT([0, p_ref, 0, i.fluid.val], T_ref))
        for o in self.outl:
            res -= o.m.val_SI * (o.h.val_SI - h_mix_pT([0, p_ref, 0, o.fluid.val], T_ref))

        return res + self.calc_ti()

    def lambda_func(self):
        r"""
        Calculates the residual for specified lambda.

        Returns
        -------
        val : float
            Residual value of function.

            .. math::

                val = \lambda - \frac{\dot{m}_{air}}{\dot{m}_{air,min}}
        """
        if self.air_alias.val in ['air', 'Air']:
            air = self.air_alias.val
        else:
            air = 'TESPy::' + self.air_alias.val
        fuel = 'TESPy::' + self.fuel_alias.val

        m_air = 0
        m_fuel = 0

        for i in self.inl:
            m_air += (i.m.val_SI * i.fluid.val[air])
            m_fuel += (i.m.val_SI * i.fluid.val[fuel])

        return self.lamb.val - m_air / (m_fuel * self.air_min)

    def ti_func(self):
        r"""
        Calculates the residual for specified thermal input.

        Returns
        -------
        val : float
            Residual value of function.

            .. math::

                val = ti - \dot{m}_f \cdot LHV
        """
        return self.ti.val - self.calc_ti()

    def calc_ti(self):
        r"""
        Calculates the thermal input of the combustion chamber.

        Returns
        -------
        ti : float
            Thermal input.

            .. math::

                ti = LHV \cdot \left[\sum_i \left(\dot{m}_{in,i} \cdot x_{f,i}
                \right) - \dot{m}_{out,1} \cdot x_{f,1} \right]
                \; \forall i \in [1,2]
        """
        fuel = 'TESPy::' + self.fuel_alias.val

        m = 0
        for i in self.inl:
            m += i.m.val_SI * i.fluid.val[fuel]

        for o in self.outl:
            m -= o.m.val_SI * o.fluid.val[fuel]

        return m * self.lhv

    def initialise_fluids(self, nw):
        r"""
        Calculates reaction balance with given lambda of 3 for good generic starting values at the component's outlet.

        Parameters
        ----------
        nw : tespy.networks.network
            Network using this component object.
        """
        if self.air_alias.val in ['air', 'Air']:
            air = self.air_alias.val
        else:
            air = 'TESPy::' + self.air_alias.val
        flue_gas = 'TESPy::' + self.fuel_alias.val + "_fg"

        for c in nw.comps.loc[self].o:
            if not c.fluid.val_set[air]:
                c.fluid.val[air] = 0.8
            if not c.fluid.val_set[flue_gas]:
                c.fluid.val[flue_gas] = 0.2

    def convergence_check(self, nw):
        r"""
        Performs a convergence check.

        Parameters
        ----------
        nw : tespy.networks.network
            The network object using this component.

        Note
        ----
        Manipulate enthalpies/pressure at inlet and outlet if not specified by user to match physically feasible constraints,
        keep fluid composition within feasible range and then propagates it towards the outlet.
        """
        if self.air_alias.val in ['air', 'Air']:
            air = self.air_alias.val
        else:
            air = 'TESPy::' + self.air_alias.val
        flue_gas = 'TESPy::' + self.fuel_alias.val + "_fg"
        fuel = 'TESPy::' + self.fuel_alias.val

        for c in nw.comps.loc[self].o:
            if not c.fluid.val_set[air]:
                if c.fluid.val[air] > 0.95:
                    c.fluid.val[air] = 0.95
                if c.fluid.val[air] < 0.5:
                    c.fluid.val[air] = 0.5

            if not c.fluid.val_set[flue_gas]:
                if c.fluid.val[flue_gas] > 0.5:
                    c.fluid.val[flue_gas] = 0.5
                if c.fluid.val[flue_gas] < 0.05:
                    c.fluid.val[flue_gas] = 0.05

            if not c.fluid.val_set[fuel]:
                if c.fluid.val[fuel] > 0:
                    c.fluid.val[fuel] = 0

            nw.init_target(c, c.t)

        for i in nw.comps.loc[self].i:
            if i.m.val_SI < 0 and not i.m.val_set:
                i.m.val_SI = 0.01

        for c in nw.comps.loc[self].o:
            if c.m.val_SI < 0 and not c.m.val_set:
                c.m.val_SI = 10
            nw.init_target(c, c.t)

        if self.lamb.val < 1 and not self.lamb.is_set:
            self.lamb.val = 2

    def calc_parameters(self, mode):
        r"""
        Post and preprocessing parameter calculation/specification.

        Parameters
        ----------

        mode : str
            Pre- or postprocessing calculation.

        Note
        ----
        Generic preprocessing is handled by the base class. This method handles class specific pre- and postprocessing.
        """
        component.calc_parameters(self, mode)

        if mode == 'post':

            if self.air_alias.val in ['air', 'Air']:
                air = self.air_alias.val
            else:
                air = 'TESPy::' + self.air_alias.val
            fuel = 'TESPy::' + self.fuel_alias.val

            m_fuel = 0
            for i in self.inl:
                m_fuel += i.m.val_SI * i.fluid.val[fuel]

            m_air = 0
            for i in self.inl:
                m_air += i.m.val_SI * i.fluid.val[air]

            self.lamb.val = (m_air / m_fuel) / self.air_min

            S = 0
            T_ref = 373.15
            p_ref = 1e5

            for i in self.inl:
                S += i.m.val_SI * (s_mix_ph(i.to_flow()) - s_mix_pT([0, p_ref, 0, i.fluid.val], T_ref))

            for o in self.outl:
                S -= o.m.val_SI * (s_mix_ph(o.to_flow()) - s_mix_pT([0, p_ref, 0, o.fluid.val], T_ref))

            self.S.val = S

            ti = 0
            for i in self.inl:
                ti += i.m.val_SI * i.fluid.val[fuel] * self.lhv

            self.ti.val = ti

# %%


class cogeneration_unit(combustion_chamber):
    r"""
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.cogeneration_unit.reaction_balance`
        - :func:`tespy.components.components.cogeneration_unit.fluid_func`
          (for cooling water)
        - :func:`tespy.components.components.cogeneration_unit.mass_flow_func`

        .. math::

            0 = p_{3,in} - p_{3,out}\\
            0 = p_{4,in} - p_{3,out}

        - :func:`tespy.components.components.cogeneration_unit.energy_balance`

        **optional equations**

        - :func:`tespy.components.components.cogeneration_unit.lambda_func`
        - :func:`tespy.components.components.cogeneration_unit.ti_func`
        - :func:`tespy.components.components.cogeneration_unit.Q1_func`
        - :func:`tespy.components.components.cogeneration_unit.Q2_func`

        .. math::

            0 = p_{1,in} \cdot pr1 - p_{1,out}\\
            0 = p_{2,in} \cdot pr2 - p_{2,out}

        - :func:`tespy.components.components.component.zeta_func`
        - :func:`tespy.components.components.component.zeta2_func`

    Available fuels

        - methane, ethane, propane, butane, hydrogen

    Inlets/Outlets

        - in1, in2 (cooling water), in3, in4 (air and fuel)
        - out1, out2 (cooling water), out3 (flue gas)

    Image

        .. image:: _images/cogeneration_unit.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    .. note::

        The fuel and the air components can be connected to either of the inlets.

    Parameters
    ----------
    label : str
        The label of the component.

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    fuel : str
        Fuel for the combustion chamber, see list of available fluids above.

    lamb : float/tespy.helpers.dc_cp
        Air to stoichiometric air ratio, :math:`\lambda/1`.

    ti : float/tespy.helpers.dc_cp
        Thermal input, (:math:`{LHV \cdot \dot{m}_f}`),
        :math:`ti/\text{W}`.

    P : str/float/tespy.helpers.dc_cp
        Power output, :math:`P/\text{W}`.

    Q1 : str/float/tespy.helpers.dc_cp
        Heat output 1, :math:`\dot Q/\text{W}`.

    Q2 : str/float/tespy.helpers.dc_cp
        Heat output 2, :math:`\dot Q/\text{W}`.

    Qloss : str/float/tespy.helpers.dc_cp
        Heat loss, :math:`\dot Q_{loss}/\text{W}`.

    pr1 : str/float/tespy.helpers.dc_cp
        Pressure ratio heat outlet 1, :math:`pr/1`.

    pr2 : str/float/tespy.helpers.dc_cp
        Pressure ratio heat outlet 2, :math:`pr/1`.

    zeta1 : str/float/tespy.helpers.dc_cp
        Pressure ratio heat outlet 2, :math:`\zeta/\frac{\text{Pa}}{\text{m}^4}`.

    zeta2 : str/float/tespy.helpers.dc_cp
        Pressure ratio heat outlet 2, :math:`\zeta/\frac{\text{Pa}}{\text{m}^4}`.

    tiP_char : str/tespy.helpers.dc_cc
        Characteristic line linking fuel input to power output.

    Q1_char : str/tespy.helpers.dc_cc
        Characteristic line linking heat output 1 to power output.

    Q2_char : str/tespy.helpers.dc_cc
        Characteristic line linking heat output 2 to power output.

    Qloss_char : str/tespy.helpers.dc_cc
        Characteristic line linking heat loss to power output.

    Note
    ----

        For more information on the usage of the cogeneration unit see the
        examples in the tespy_examples repository.

    Example
    -------
    >>> from tespy import cmp, con, nwk
    >>> import shutil
    >>> import numpy as np
    >>> fluid_list = ['Ar', 'N2', 'O2', 'CO2', 'CH4', 'H2O']
    >>> nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ...     p_range=[0.5, 10], T_range=[10, 1200])
    >>> nw.set_printoptions(print_level='none')
    >>> amb = cmp.source('ambient')
    >>> sf = cmp.source('fuel')
    >>> fg = cmp.sink('flue gas outlet')
    >>> cw_in1 = cmp.source('cooling water inlet1')
    >>> cw_in2 = cmp.source('cooling water inlet2')
    >>> cw_out1 = cmp.sink('cooling water outlet1')
    >>> cw_out2 = cmp.sink('cooling water outlet2')
    >>> split = cmp.splitter('splitter')
    >>> merge = cmp.merge('merge')
    >>> chp = cmp.cogeneration_unit(label='cogeneration unit')
    >>> amb_comb = con.connection(amb, 'out1', chp, 'in3')
    >>> sf_comb = con.connection(sf, 'out1', chp, 'in4')
    >>> comb_fg = con.connection(chp, 'out3', fg, 'in1')
    >>> nw.add_conns(sf_comb, amb_comb, comb_fg)
    >>> cw1_chp1 = con.connection(cw_in1, 'out1', chp, 'in1')
    >>> cw2_chp2 = con.connection(cw_in2, 'out1', chp, 'in2')
    >>> nw.add_conns(cw1_chp1, cw2_chp2)
    >>> chp1_cw = con.connection(chp, 'out1', cw_out1, 'in1')
    >>> chp2_cw = con.connection(chp, 'out2', cw_out2, 'in1')
    >>> nw.add_conns(chp1_cw, chp2_cw)
    >>> chp.set_attr(fuel='CH4', pr1=0.99, pr2=0.99, P=10e6, lamb=1.2, design=['pr1', 'pr2'], offdesign=['zeta1', 'zeta2'])
    >>> amb_comb.set_attr(p=5, T=30,
    ...     fluid={'Ar': 0.0129, 'N2': 0.7553, 'H2O': 0, 'CH4': 0,
    ...         'CO2': 0.0004, 'O2': 0.2314})
    >>> sf_comb.set_attr(T=30,
    ...     fluid={'CO2': 0, 'Ar': 0, 'N2': 0, 'O2': 0, 'H2O': 0, 'CH4': 1})
    >>> cw1_chp1.set_attr(p=3, T=60, m=50,
    ...     fluid={'CO2': 0, 'Ar': 0, 'N2': 0, 'O2': 0, 'H2O': 1, 'CH4': 0})
    >>> cw2_chp2.set_attr(p=3, T=80, m=50,
    ...     fluid={'CO2': 0, 'Ar': 0, 'N2': 0, 'O2': 0, 'H2O': 1, 'CH4': 0})
    >>> mode = 'design'
    >>> nw.solve(mode=mode)
    >>> nw.save('tmp')
    >>> round(chp.ti.val)
    22500000.0
    >>> round(chp.Q1.val)
    1743636.0
    >>> chp.set_attr(Q1=1.5e6, P=np.nan)
    >>> mode = 'offdesign'
    >>> nw.solve(mode=mode, init_path='tmp', design_path='tmp')
    >>> round(chp.ti.val)
    17427210.0
    >>> round(chp.P.val / chp.P.design, 3)
    0.747
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'cogeneration unit'

    def attr(self):
        return {'fuel': dc_simple(), 'lamb': dc_cp(), 'ti': dc_cp(),
                'P': dc_cp(val=1e6, d=1, val_min=1),
                'Q1': dc_cp(), 'Q2': dc_cp(),
                'Qloss': dc_cp(val=1e5, d=1, val_min=1),
                'pr1': dc_cp(), 'pr2': dc_cp(),
                'zeta1': dc_cp(), 'zeta2': dc_cp(),
                'tiP_char': dc_cc(method='TI'),
                'Q1_char': dc_cc(method='Q1'),
                'Q2_char': dc_cc(method='Q2'),
                'Qloss_char': dc_cc(method='QLOSS'),
                'S': dc_cp()}

    def inlets(self):
        return ['in1', 'in2', 'in3', 'in4']

    def outlets(self):
        return ['out1', 'out2', 'out3']

    def comp_init(self, nw):

        if not self.P.is_set:
            self.set_attr(P='var')
            msg = ('The power output of cogeneration units must be set! '
                   'We are adding the power output of component ' +
                   self.label + ' as custom variable of the system.')
            logging.info(msg)

        if not self.Qloss.is_set:
            self.set_attr(Qloss='var')
            msg = ('The heat loss of cogeneration units must be set! '
                   'We are adding the heat loss of component ' +
                   self.label + ' as custom variable of the system.')
            logging.info(msg)

        combustion_chamber.comp_init(self, nw)

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
        # equations for fluids in combustion chamber
        for fluid in self.inl[0].fluid.val.keys():
            vec_res += [self.reaction_balance(fluid)]

        ######################################################################
        # equations for fluids in cooling loops
        vec_res += self.fluid_func()

        ######################################################################
        # equations for mass flow
        vec_res += self.mass_flow_func()

        ######################################################################
        # equations for pressure balance in combustion
        vec_res += [self.inl[2].p.val_SI - self.outl[2].p.val_SI]
        vec_res += [self.inl[2].p.val_SI - self.inl[3].p.val_SI]

        ######################################################################
        # equation for cogeneration unit energy balance
        vec_res += [self.energy_balance()]

        ######################################################################
        # equation for power to thermal input ratio from characteristic line
        vec_res += [self.tiP_char_func()]

        ######################################################################
        # equations for heat outputs from characteristic line
        vec_res += [self.Q1_char_func()]
        vec_res += [self.Q2_char_func()]

        ######################################################################
        # equation for heat loss from characteristic line
        vec_res += [self.Qloss_char_func()]

        ######################################################################
        # equation for specified lambda
        if self.lamb.is_set:
            vec_res += [self.lambda_func()]

        ######################################################################
        # equation for specified thermal input
        if self.ti.is_set:
            vec_res += [self.ti_func()]

        ######################################################################
        # equations for specified heat ouptputs
        if self.Q1.is_set:
            vec_res += [self.Q1_func()]

        if self.Q2.is_set:
            vec_res += [self.Q2_func()]

        ######################################################################
        # equations for specified pressure ratios at cooling loops
        if self.pr1.is_set:
            vec_res += [self.pr1.val * self.inl[0].p.val_SI - self.outl[0].p.val_SI]

        if self.pr2.is_set:
            vec_res += [self.pr2.val * self.inl[1].p.val_SI - self.outl[1].p.val_SI]

        ######################################################################
        # equations for specified zeta values at cooling loops
        if self.zeta1.is_set:
            vec_res += [self.zeta_func()]

        if self.zeta2.is_set:
            vec_res += [self.zeta2_func()]

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
        # derivatives for reaction balance
        deriv = np.zeros((self.num_fl, 7 + self.num_vars, self.num_fl + 3))
        j = 0
        for fluid in self.fluids:

            # fresh air and fuel inlets
            for i in range(2):
                deriv[j, i + 2, 0] = self.rb_numeric_deriv('m', i + 2, fluid)
                deriv[j, i + 2, 3:] = self.rb_numeric_deriv('fluid', i + 2, fluid)

            # combustion outlet
            deriv[j, 6, 0] = self.rb_numeric_deriv('m', 6, fluid)
            deriv[j, 6, 3:] = self.rb_numeric_deriv('fluid', 6, fluid)
            j += 1
        mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for cooling water fluid composition and mass flow
        mat_deriv += self.fl_deriv
        mat_deriv += self.m_deriv

        ######################################################################
        # derivatives for pressure equations
        mat_deriv += self.p_deriv

        ######################################################################
        # derivatives for energy balance
        eb_deriv = np.zeros((1, 7 + self.num_vars, self.num_fl + 3))

        # mass flow cooling water
        for i in [0, 1]:
            eb_deriv[0, i, 0] = -(self.outl[i].h.val_SI - self.inl[i].h.val_SI)

        # mass flow and pressure for combustion reaction
        for i in [2, 3, 6]:
            eb_deriv[0, i, 0] = self.numeric_deriv(self.energy_balance, 'm', i)
            eb_deriv[0, i, 1] = self.numeric_deriv(self.energy_balance, 'p', i)

        # enthalpy
        for i in range(4):
            eb_deriv[0, i, 2] = self.inl[i].m.val_SI
        for i in range(3):
            eb_deriv[0, i + 4, 2] = -self.outl[i].m.val_SI

        # fluid composition
        pos = 3 + self.fluids.index(self.fuel.val)
        eb_deriv[0, 2, pos] = self.inl[2].m.val_SI * self.lhv
        eb_deriv[0, 3, pos] = self.inl[3].m.val_SI * self.lhv
        eb_deriv[0, 6, pos] = -self.outl[2].m.val_SI * self.lhv

        # power and heat loss
        if self.P.is_var:
            eb_deriv[0, 7 + self.P.var_pos, 0] = self.numeric_deriv(self.energy_balance, 'P', 7)
        if self.Qloss.is_var:
            eb_deriv[0, 7 + self.Qloss.var_pos, 0] = self.numeric_deriv(self.energy_balance, 'Qloss', 7)
        mat_deriv += eb_deriv.tolist()

        ######################################################################
        # derivatives for thermal input to power charactersitics
        tiP_deriv = np.zeros((1, 7 + self.num_vars, self.num_fl + 3))
        for i in range(2):
            tiP_deriv[0, i + 2, 0] = self.numeric_deriv(self.tiP_char_func, 'm', i + 2)
            tiP_deriv[0, i + 2, 3:] = self.numeric_deriv(self.tiP_char_func, 'fluid', i + 2)

        tiP_deriv[0, 6, 0] = self.numeric_deriv(self.tiP_char_func, 'm', 6)
        tiP_deriv[0, 6, 3:] = self.numeric_deriv(self.tiP_char_func, 'fluid', 6)

        if self.P.is_var:
            tiP_deriv[0, 7 + self.P.var_pos, 0] = self.numeric_deriv(self.tiP_char_func, 'P', 7)
        mat_deriv += tiP_deriv.tolist()

        ######################################################################
        # derivatives for heat output 1 to power charactersitics
        Q1_deriv = np.zeros((1, 7 + self.num_vars, self.num_fl + 3))
        Q1_deriv[0, 0, 0] = self.numeric_deriv(self.Q1_char_func, 'm', 0)
        Q1_deriv[0, 0, 2] = self.numeric_deriv(self.Q1_char_func, 'h', 0)
        Q1_deriv[0, 4, 2] = self.numeric_deriv(self.Q1_char_func, 'h', 4)
        for i in range(2):
            Q1_deriv[0, i + 2, 0] = self.numeric_deriv(self.Q1_char_func, 'm', i + 2)
            Q1_deriv[0, i + 2, 3:] = self.numeric_deriv(self.Q1_char_func, 'fluid', i + 2)
        Q1_deriv[0, 6, 0] = self.numeric_deriv(self.Q1_char_func, 'm', 6)
        Q1_deriv[0, 6, 3:] = self.numeric_deriv(self.Q1_char_func, 'fluid', 6)

        if self.P.is_var:
            Q1_deriv[0, 7 + self.P.var_pos, 0] = self.numeric_deriv(self.Q1_char_func, 'P', 7)
        mat_deriv += Q1_deriv.tolist()

        ######################################################################
        # derivatives for heat output 2 to power charactersitics
        Q2_deriv = np.zeros((1, 7 + self.num_vars, self.num_fl + 3))
        Q2_deriv[0, 1, 0] = self.numeric_deriv(self.Q2_char_func, 'm', 1)
        Q2_deriv[0, 1, 2] = self.numeric_deriv(self.Q2_char_func, 'h', 1)
        Q2_deriv[0, 5, 2] = self.numeric_deriv(self.Q2_char_func, 'h', 5)
        for i in range(2):
            Q2_deriv[0, i + 2, 0] = self.numeric_deriv(self.Q2_char_func, 'm', i + 2)
            Q2_deriv[0, i + 2, 3:] = self.numeric_deriv(self.Q2_char_func, 'fluid', i + 2)
        Q2_deriv[0, 6, 0] = self.numeric_deriv(self.Q2_char_func, 'm', 6)
        Q2_deriv[0, 6, 3:] = self.numeric_deriv(self.Q2_char_func, 'fluid', 6)

        if self.P.is_var:
            Q2_deriv[0, 7 + self.P.var_pos, 0] = self.numeric_deriv(self.Q2_char_func, 'P', 7)
        mat_deriv += Q2_deriv.tolist()

        ######################################################################
        # derivatives for heat loss to power charactersitics
        Ql_deriv = np.zeros((1, 7 + self.num_vars, self.num_fl + 3))
        for i in range(2):
            Ql_deriv[0, i + 2, 0] = self.numeric_deriv(self.Qloss_char_func, 'm', i + 2)
            Ql_deriv[0, i + 2, 3:] = self.numeric_deriv(self.Qloss_char_func, 'fluid', i + 2)
        Ql_deriv[0, 6, 0] = self.numeric_deriv(self.Qloss_char_func, 'm', 6)
        Ql_deriv[0, 6, 3:] = self.numeric_deriv(self.Qloss_char_func, 'fluid', 6)

        if self.P.is_var:
            Ql_deriv[0, 7 + self.P.var_pos, 0] = self.numeric_deriv(self.Qloss_char_func, 'P', 7)
        if self.Qloss.is_var:
            Ql_deriv[0, 7 + self.Qloss.var_pos, 0] = self.numeric_deriv(self.Qloss_char_func, 'Qloss', 7)
        mat_deriv += Ql_deriv.tolist()

        ######################################################################
        # derivatives for specified lambda
        if self.lamb.is_set:
            lamb_deriv = np.zeros((1, 7 + self.num_vars, self.num_fl + 3))
            for i in range(2):
                lamb_deriv[0, i + 2, 0] = self.numeric_deriv(self.lambda_func, 'm', i + 2)
                lamb_deriv[0, i + 2, 3:] = self.numeric_deriv(self.lambda_func, 'fluid', i + 2)
            mat_deriv += lamb_deriv.tolist()

        ######################################################################
        # derivatives for specified thermal input
        if self.ti.is_set:
            ti_deriv = np.zeros((1, 7 + self.num_vars, self.num_fl + 3))
            for i in range(2):
                ti_deriv[0, i + 2, 0] = self.numeric_deriv(self.ti_func, 'm', i + 2)
                ti_deriv[0, i + 2, 3:] = self.numeric_deriv(self.ti_func, 'fluid', i + 2)
            ti_deriv[0, 6, 0] = self.numeric_deriv(self.ti_func, 'm', 6)
            ti_deriv[0, 6, 3:] = self.numeric_deriv(self.ti_func, 'fluid', 6)
            mat_deriv += ti_deriv.tolist()

        ######################################################################
        # derivatives for specified heat outputs
        if self.Q1.is_set:
            Q_deriv = np.zeros((1, 7 + self.num_vars, self.num_fl + 3))
            Q_deriv[0, 0, 0] = - (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
            Q_deriv[0, 0, 2] = self.inl[0].m.val_SI
            Q_deriv[0, 4, 2] = -self.inl[0].m.val_SI
            mat_deriv += Q_deriv.tolist()

        if self.Q2.is_set:
            Q_deriv = np.zeros((1, 7 + self.num_vars, self.num_fl + 3))
            Q_deriv[0, 1, 0] = - (self.outl[1].h.val_SI - self.inl[1].h.val_SI)
            Q_deriv[0, 1, 2] = self.inl[1].m.val_SI
            Q_deriv[0, 5, 2] = -self.inl[1].m.val_SI
            mat_deriv += Q_deriv.tolist()

        ######################################################################
        # derivatives for specified pressure ratio at cooling loops
        if self.pr1.is_set:
            pr1_deriv = np.zeros((1, 7 + self.num_vars, self.num_fl + 3))
            pr1_deriv[0, 0, 1] = self.pr1.val
            pr1_deriv[0, 4, 1] = -1
            mat_deriv += pr1_deriv.tolist()

        if self.pr2.is_set:
            pr2_deriv = np.zeros((1, 7 + self.num_vars, self.num_fl + 3))
            pr2_deriv[0, 1, 1] = self.pr2.val
            pr2_deriv[0, 5, 1] = -1
            mat_deriv += pr2_deriv.tolist()

        ######################################################################
        # derivatives for specified zeta values at cooling loops
        if self.zeta1.is_set:
            zeta1_deriv = np.zeros((1, 7 + self.num_vars, self.num_fl + 3))
            zeta1_deriv[0, 0, 0] = self.numeric_deriv(self.zeta_func, 'm', 0)
            zeta1_deriv[0, 0, 1] = self.numeric_deriv(self.zeta_func, 'p', 0)
            zeta1_deriv[0, 0, 2] = self.numeric_deriv(self.zeta_func, 'h', 0)
            zeta1_deriv[0, 4, 1] = self.numeric_deriv(self.zeta_func, 'p', 4)
            zeta1_deriv[0, 4, 2] = self.numeric_deriv(self.zeta_func, 'h', 4)
            mat_deriv += zeta1_deriv.tolist()

        if self.zeta2.is_set:
            zeta2_deriv = np.zeros((1, 7 + self.num_vars, self.num_fl + 3))
            zeta2_deriv[0, 1, 0] = self.numeric_deriv(self.zeta2_func, 'm', 1)
            zeta2_deriv[0, 1, 1] = self.numeric_deriv(self.zeta2_func, 'p', 1)
            zeta2_deriv[0, 1, 2] = self.numeric_deriv(self.zeta2_func, 'h', 1)
            zeta2_deriv[0, 5, 1] = self.numeric_deriv(self.zeta2_func, 'p', 5)
            zeta2_deriv[0, 5, 2] = self.numeric_deriv(self.zeta2_func, 'h', 5)
            mat_deriv += zeta2_deriv.tolist()

        return np.asarray(mat_deriv)

    def fluid_func(self):
        r"""
        Calculates the vector of residual values for cooling loop fluid balance equations.

        Returns
        -------
        vec_res : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = fluid_{i,in_{j}} - fluid_{i,out_{j}}\\
                \forall i \in \mathrm{fluid}, \; \forall j \in [1, 2]
        """
        vec_res = []

        for i in range(2):
            for fluid, x in self.inl[i].fluid.val.items():
                vec_res += [x - self.outl[i].fluid.val[fluid]]
        return vec_res

    def mass_flow_func(self):
        r"""
        Calculates the residual value for component's mass flow balance equation.

        Returns
        -------
        vec_res : list
            Vector with residual value for component's mass flow balance.

            .. math::

                0 = \dot{m}_{in,i} - \dot{m}_{out,i}\\
                \forall i \in [1, 2]\\
                0 = \dot{m}_{in,3} + \dot{m}_{in,4} - \dot{m}_{out,3}
        """

        vec_res = []
        for i in range(2):
            vec_res += [self.inl[i].m.val_SI - self.outl[i].m.val_SI]
        vec_res += [self.inl[2].m.val_SI + self.inl[3].m.val_SI -
                    self.outl[2].m.val_SI]
        return vec_res

    def fluid_deriv(self):
        r"""
        Calculates the partial derivatives for cooling loop fluid balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((self.num_fl * 2, 7 + self.num_vars, 3 + self.num_fl))
        for i in range(self.num_fl):
            deriv[i, 0, i + 3] = 1
            deriv[i, 4, i + 3] = -1
        for j in range(self.num_fl):
            deriv[i + 1 + j, 1, j + 3] = 1
            deriv[i + 1 + j, 5, j + 3] = -1
        return deriv.tolist()

    def mass_flow_deriv(self):
        r"""
        Calculates the partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((3, 7 + self.num_vars, self.num_fl + 3))
        for i in range(2):
            deriv[i, i, 0] = 1
        for j in range(2):
            deriv[j, self.num_i + j, 0] = -1
        deriv[2, 2, 0] = 1
        deriv[2, 3, 0] = 1
        deriv[2, 6, 0] = -1
        return deriv.tolist()

    def pressure_deriv(self):
        r"""
        Calculates the partial derivatives for combustion pressure equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((2, 7 + self.num_vars, self.num_fl + 3))
        for k in range(2):
            deriv[k, 2, 1] = 1
        deriv[0, 6, 1] = -1
        deriv[1, 3, 1] = -1
        return deriv.tolist()

    def energy_balance(self):
        r"""
        Calculates the energy balance of the cogeneration unit.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \sum_i \dot{m}_{in,i} \cdot \left( h_{in,i} - h_{in,i,ref}
                \right)\\
                & - \sum_j \dot{m}_{out,3} \cdot \left( h_{out,3} - h_{out,3,ref}
                \right)\\
                & + H_{I,f} \cdot \left(\sum_i \left(\dot{m}_{in,i} \cdot x_{f,i}
                \right)- \dot{m}_{out,3} \cdot x_{f,3} \right)\\
                & - \dot{Q}_1 - \dot{Q}_2 - P - \dot{Q}_{loss}\\
                \end{split}\\
                \forall i \in [3,4]

        Note
        ----
        The temperature for the reference state is set to 20 °C, thus
        the water may be liquid. In order to make sure, the state is
        referring to the lower heating value, the necessary enthalpy
        difference for evaporation is added.

        - Reference temperature: 293.15 K.
        - Reference pressure: 1 bar.
        """
        T_ref = 293.15
        p_ref = 1e5

        res = 0
        for i in self.inl[2:]:
            res += i.m.val_SI * (i.h.val_SI - h_mix_pT([0, p_ref, 0, i.fluid.val], T_ref))

        for o in self.outl[2:]:
            dh = 0
            n_h2o = o.fluid.val[self.h2o] / molar_masses[self.h2o]
            if n_h2o > 0:
                p = p_ref * n_h2o / molar_mass_flow(o.fluid.val)
                h = h_pT(p, T_ref, self.h2o)
                h_steam = CP.PropsSI('H', 'P', p, 'Q', 1, self.h2o)
                if h < h_steam:
                    dh = (h_steam - h) * o.fluid.val[self.h2o]

            res -= o.m.val_SI * (o.h.val_SI - h_mix_pT([0, p_ref, 0, o.fluid.val], T_ref) - dh)

        res += self.calc_ti()

        # cooling water
        for i in range(2):
            res -= self.inl[i].m.val_SI * (self.outl[i].h.val_SI - self.inl[i].h.val_SI)

        # power output and heat loss
        res -= self.P.val + self.Qloss.val

        return res

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
                LHV \cdot \dot{m}_{f} \cdot f_{char}\left( \frac{LHV \cdot \dot{m}_{f}}{LHV \cdot \dot{m}_{f, ref}}\right) & \text{key = 'TI'}\\
                P \cdot f_{char}\left( \frac{P}{P_{ref}}\right) & \text{key = 'P'}\\
                \left(\dot{Q}_1 + \dot{Q}_2\right) \cdot f_{char}\left( \frac{\dot{Q}_1 + \dot{Q}_2}{\dot{Q}_{1,ref} + \dot{Q}_{2,ref}}\right)& \text{key = 'Q'}\\
                \dot{Q}_1 \cdot f_{char}\left( \frac{\dot{Q}_1}{\dot{Q}_{1,ref}}\right) & \text{key = 'Q1'}\\
                \dot{Q}_2 \cdot f_{char}\left( \frac{\dot{Q}_2}{\dot{Q}_{2,ref}}\right) & \text{key = 'Q2'}\\
                \dot{Q}_{loss} \cdot f_{char}\left( \frac{\dot{Q}_{loss}}{\dot{Q}_{loss,ref}}\right) & \text{key = 'Qloss'}
                \end{cases}

                \dot{Q}_1=\dot{m}_1 \cdot \left( h_{1,out} - h_{1,in} \right)\\
                \dot{Q}_2=\dot{m}_2 \cdot \left( h_{2,out} - h_{2,in} \right)
        """
        if bus.param == 'TI':
            ti = self.calc_ti()
            if np.isnan(bus.P_ref):
                expr = 1
            else:
                expr = abs(ti / bus.P_ref)
            return ti * bus.char.f_x(expr)

        elif bus.param == 'P':
            P = self.calc_P()
            if np.isnan(bus.P_ref):
                expr = 1
            else:
                expr = abs(P / bus.P_ref)
            return P * bus.char.f_x(expr)

        elif bus.param == 'Q':
            val = 0
            for j in range(2):
                i = self.inl[j]
                o = self.outl[j]
                val += i.m.val_SI * (o.h.val_SI - i.h.val_SI)

            if np.isnan(bus.P_ref):
                expr = 1
            else:
                expr = abs(val / bus.P_ref)
            return val * bus.char.f_x(expr)

        elif bus.param == 'Q1':
            i = self.inl[0]
            o = self.outl[0]
            val = i.m.val_SI * (o.h.val_SI - i.h.val_SI)

            if np.isnan(bus.P_ref):
                expr = 1
            else:
                expr = abs(val / bus.P_ref)
            return val * bus.char.f_x(expr)

        elif bus.param == 'Q2':
            i = self.inl[1]
            o = self.outl[1]
            val = i.m.val_SI * (o.h.val_SI - i.h.val_SI)

            if np.isnan(bus.P_ref):
                expr = 1
            else:
                expr = abs(val / bus.P_ref)
            return val * bus.char.f_x(expr)

        elif bus.param == 'Qloss':
            Q = self.calc_Qloss()
            if np.isnan(bus.P_ref):
                expr = 1
            else:
                expr = abs(Q / bus.P_ref)
            return Q * bus.char.f_x(expr)

        else:
            msg = 'The parameter ' + bus.param + 'is not a valid parameter for a ' + self.component() + '.'
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
        deriv = np.zeros((1, 7 + self.num_vars, len(self.inl[0].fluid.val) + 3))

        ######################################################################
        # derivatives for specified zeta values at cooling loops
        if bus.param == 'TI':
            for i in range(2):
                deriv[0, i + 2, 0] = self.numeric_deriv(self.bus_func, 'm', i + 2, bus=bus)
                deriv[0, i + 2, 3:] = self.numeric_deriv(self.bus_func, 'fluid', i + 2, bus=bus)
            deriv[0, 6, 0] = self.numeric_deriv(self.bus_func, 'm', 6, bus=bus)
            deriv[0, 6, 3:] = self.numeric_deriv(self.bus_func, 'fluid', 6, bus=bus)

        ######################################################################
        # derivatives for specified zeta values at cooling loops
        elif bus.param == 'P':
            for i in range(2):
                deriv[0, i + 2, 0] = self.numeric_deriv(self.bus_func, 'm', i + 2, bus=bus)
                deriv[0, i + 2, 3:] = self.numeric_deriv(self.bus_func, 'fluid', i + 2, bus=bus)

            deriv[0, 6, 0] = self.numeric_deriv(self.bus_func, 'm', 6, bus=bus)
            deriv[0, 6, 3:] = self.numeric_deriv(self.bus_func, 'fluid', 6, bus=bus)

            # variable power
            if self.P.is_var:
                deriv[0, 7 + self.P.var_pos, 0] = self.numeric_deriv(self.bus_func, 'P', 7, bus=bus)

        ######################################################################
        # derivatives for specified zeta values at cooling loops
        elif bus.param == 'Q':
            for i in range(2):
                deriv[0, i, 0] = self.numeric_deriv(self.bus_func, 'm', i, bus=bus)
                deriv[0, i, 2] = self.numeric_deriv(self.bus_func, 'h', i, bus=bus)
                deriv[0, i + 4, 2] = self.numeric_deriv(self.bus_func, 'h', i + 4, bus=bus)

        ######################################################################
        # derivatives for specified zeta values at cooling loops
        elif bus.param == 'Q1':
            deriv[0, 0, 0] = self.numeric_deriv(self.bus_func, 'm', 0, bus=bus)
            deriv[0, 0, 2] = self.numeric_deriv(self.bus_func, 'h', 0, bus=bus)
            deriv[0, 4, 2] = self.numeric_deriv(self.bus_func, 'h', 4, bus=bus)

        ######################################################################
        # derivatives for specified zeta values at cooling loops
        elif bus.param == 'Q2':
            deriv[0, 1, 0] = self.numeric_deriv(self.bus_func, 'm', 1, bus=bus)
            deriv[0, 1, 2] = self.numeric_deriv(self.bus_func, 'h', 1, bus=bus)
            deriv[0, 5, 2] = self.numeric_deriv(self.bus_func, 'h', 5, bus=bus)

        ######################################################################
        # derivatives for specified zeta values at cooling loops
        elif bus.param == 'Qloss':
            for i in range(2):
                deriv[0, i + 2, 0] = self.numeric_deriv(self.bus_func, 'm', i + 2, bus=bus)
                deriv[0, i + 2, 3:] = self.numeric_deriv(self.bus_func, 'fluid', i + 2, bus=bus)

            deriv[0, 6, 0] = self.numeric_deriv(self.bus_func, 'm', 6, bus=bus)
            deriv[0, 6, 3:] = self.numeric_deriv(self.bus_func, 'fluid', 6, bus=bus)

            # variable power
            if self.P.is_var:
                deriv[0, 7 + self.P.var_pos, 0] = self.numeric_deriv(self.bus_func, 'P', 7, bus=bus)

        else:
            msg = 'The parameter ' + bus.param + 'is not a valid parameter for a ' + self.component() + '.'
            logging.error(msg)
            raise ValueError(msg)

        return deriv

    def Q1_func(self):
        r"""
        Calculates residual value with specified Q1.

        Returns
        -------
        val : float
            Residual value of equation.

            .. math::

                val = \dot{m}_1 \cdot \left(h_{out,1} - h_{in,1} \right) - \dot{Q}_1
        """
        i = self.inl[0]
        o = self.outl[0]

        return self.Q1.val - i.m.val_SI * (o.h.val_SI - i.h.val_SI)

    def Q2_func(self):
        r"""
        Calculates residual value with specified Q2.

        Returns
        -------
        val : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_2 \cdot \left(h_{out,2} - h_{in,2} \right) - \dot{Q}_2
        """
        i = self.inl[1]
        o = self.outl[1]

        return self.Q2.val - i.m.val_SI * (o.h.val_SI - i.h.val_SI)

    def tiP_char_func(self):
        r"""
        Calculates the relation of output power and thermal input from specified characteristic line.

        Returns
        -------
        val : float
            Residual value of equation.

            .. math::

                0 = P \cdot f_{TI}\left(\frac{P}{P_{ref}}\right)- LHV \cdot
                \left[\sum_i \left(\dot{m}_{in,i} \cdot
                x_{f,i}\right) - \dot{m}_{out,3} \cdot x_{f,3} \right]
                \; \forall i \in [1,2]
        """
        if np.isnan(self.P.design):
            expr = 1
        else:
            expr = self.P.val / self.P.design

        return self.calc_ti() - self.tiP_char.func.f_x(expr) * self.P.val

    def Q1_char_func(self):
        r"""
        Calculates the relation of heat output 1 and thermal input from specified characteristic lines.

        Returns
        -------
        val : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \dot{m}_1 \cdot \left(h_{out,1} - h_{in,1} \right) \cdot
                f_{TI}\left(\frac{P}{P_{ref}}\right) \\
                & - LHV \cdot \left[\sum_i
                \left(\dot{m}_{in,i} \cdot x_{f,i}\right) -
                \dot{m}_{out,3} \cdot x_{f,3} \right] \cdot
                f_{Q1}\left(\frac{P}{P_{ref}}\right)\\
                \end{split}\\
                \forall i \in [3,4]
        """
        i = self.inl[0]
        o = self.outl[0]

        if np.isnan(self.P.design):
            expr = 1
        else:
            expr = self.P.val / self.P.design

        return (self.calc_ti() * self.Q1_char.func.f_x(expr) -
                self.tiP_char.func.f_x(expr) * i.m.val_SI * (o.h.val_SI - i.h.val_SI))

    def Q2_char_func(self):
        r"""
        Calculates the relation of heat output 2 and thermal input from specified characteristic lines.

        Returns
        -------
        val : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \dot{m}_2 \cdot \left(h_{out,2} - h_{in,2} \right) \cdot
                f_{TI}\left(\frac{P}{P_{ref}}\right) \\
                & - LHV \cdot \left[\sum_i
                \left(\dot{m}_{in,i} \cdot x_{f,i}\right) -
                \dot{m}_{out,3} \cdot x_{f,3} \right] \cdot
                f_{Q2}\left(\frac{P}{P_{ref}}\right)\\
                \end{split}\\
                \forall i \in [3,4]
        """
        i = self.inl[1]
        o = self.outl[1]

        if np.isnan(self.P.design):
            expr = 1
        else:
            expr = self.P.val / self.P.design

        return (self.calc_ti() * self.Q2_char.func.f_x(expr) -
                self.tiP_char.func.f_x(expr) * i.m.val_SI * (o.h.val_SI - i.h.val_SI))

    def Qloss_char_func(self):
        r"""
        Calculates the relation of heat loss and thermal input from specified characteristic lines.

        Returns
        -------
        val : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \dot{Q}_{loss} \cdot
                f_{TI}\left(\frac{P}{P_{ref}}\right) \\
                & - LHV \cdot \left[\sum_i
                \left(\dot{m}_{in,i} \cdot x_{f,i}\right) -
                \dot{m}_{out,3} \cdot x_{f,3} \right] \cdot
                f_{QLOSS}\left(\frac{P}{P_{ref}}\right)\\
                \end{split}\\
                \forall i \in [3,4]
        """
        if np.isnan(self.P.design):
            expr = 1
        else:
            expr = self.P.val / self.P.design

        return (self.calc_ti() * self.Qloss_char.func.f_x(expr) -
                self.tiP_char.func.f_x(expr) * self.Qloss.val)

    def calc_ti(self):
        r"""
        Calculates the thermal input of the cogeneration unit.

        Returns
        -------
        ti : float
            Thermal input.

            .. math::

                ti = LHV \cdot \left[\sum_i \left(\dot{m}_{in,i} \cdot x_{f,i}
                \right) - \dot{m}_{out,3} \cdot x_{f,3} \right]

                \forall i \in [3,4]
        """
        m = 0
        for i in self.inl[2:]:
            m += i.m.val_SI * i.fluid.val[self.fuel.val]

        for o in self.outl[2:]:
            m -= o.m.val_SI * o.fluid.val[self.fuel.val]

        return m * self.lhv

    def calc_P(self):
        r"""
        Calculates the power output of the cogeneration unit.

        Returns
        -------
        P : float
            Power output.

            .. math::

                P = \frac{LHV \cdot \dot{m}_{f}}{f_{TI}\left(\frac{P}{P_{ref}}\right)}

        """
        if np.isnan(self.P.design):
            expr = 1
        else:
            expr = self.P.val / self.P.design

        return self.calc_ti() / self.tiP_char.func.f_x(expr)

    def calc_Qloss(self):
        r"""
        Calculates the heat loss of the cogeneration unit.

        Returns
        -------
        Qloss : float
            Heat loss.

            .. math::

                \dot{Q}_{loss} = \frac{LHV \cdot \dot{m}_{f} \cdot
                f_{QLOSS}\left(\frac{P}{P_{ref}}\right)}
                {f_{TI}\left(\frac{P}{P_{ref}}\right)}
        """
        if np.isnan(self.P.design):
            expr = 1
        else:
            expr = self.P.val / self.P.design

        return (self.calc_ti() * self.Qloss_char.func.f_x(expr) / self.tiP_char.func.f_x(expr))

    def initialise_fluids(self, nw):
        r"""
        Calculates reaction balance with given lambda of 3 for good generic starting values at the combustion's outlet.

        Parameters
        ----------
        nw : tespy.networks.network
            Network using this component object.
        """
        N_2 = 0.7655
        O_2 = 0.2345

        n_fuel = 1
        lamb = 3
        m_co2 = n_fuel * self.n['C'] * molar_masses[self.co2]
        m_h2o = n_fuel * self.n['H'] / 2 * molar_masses[self.h2o]

        n_o2 = (m_co2 / molar_masses[self.co2] + 0.5 * m_h2o / molar_masses[self.h2o]) * lamb

        m_air = n_o2 * molar_masses[self.o2] / O_2
        m_fuel = n_fuel * molar_masses[self.fuel.val]
        m_fg = m_air + m_fuel

        m_o2 = n_o2 * molar_masses[self.o2] * (1 - 1 / lamb)
        m_n2 = N_2 * m_air

        fg = {
            self.n2: m_n2 / m_fg,
            self.co2: m_co2 / m_fg,
            self.o2: m_o2 / m_fg,
            self.h2o: m_h2o / m_fg
        }

        o = self.outl[2]
        for fluid, x in o.fluid.val.items():
            if not o.fluid.val_set[fluid] and fluid in fg.keys():
                o.fluid.val[fluid] = fg[fluid]

    def initialise_source(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's outlet.

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
                10^6 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            return 10e5

    def initialise_target(self, c, key):
        r"""
        Returns a starting value for pressure and enthalpy at component's inlet.

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

    def calc_parameters(self, mode):
        r"""
        Post and preprocessing parameter calculation/specification.

        Parameters
        ----------

        mode : str
            Pre- or postprocessing calculation.

        Note
        ----
        Generic preprocessing is handled by the base class. This method handles class specific pre- and postprocessing.
        """
        combustion_chamber.calc_parameters(self, mode)

        if mode == 'post':
            i1 = self.inl[0].to_flow()
            i2 = self.inl[1].to_flow()
            o1 = self.outl[0].to_flow()
            o2 = self.outl[1].to_flow()

            self.pr1.val = o1[1] / i1[1]
            self.pr2.val = o2[1] / i2[1]
            self.zeta1.val = (i1[1] - o1[1]) * math.pi ** 2 / (8 * i1[0] ** 2 * (v_mix_ph(i1) + v_mix_ph(o1)) / 2)
            self.zeta2.val = (i2[1] - o2[1]) * math.pi ** 2 / (8 * i2[0] ** 2 * (v_mix_ph(i2) + v_mix_ph(o2)) / 2)
            self.Q1.val = i1[0] * (o1[2] - i1[2])
            self.Q2.val = i2[0] * (o2[2] - i2[2])
            self.P.val = self.calc_P()
            self.Qloss.val = self.calc_Qloss()

            # get bound errors for characteristic lines
            if np.isnan(self.P.design):
                expr = 1
            else:
                expr = self.P.val / self.P.design
            self.tiP_char.func.get_bound_errors(expr)
            self.Qloss_char.func.get_bound_errors(expr)
            self.Q1_char.func.get_bound_errors(expr)
            self.Q2_char.func.get_bound_errors(expr)

# %%


class valve(component):
    r"""
    The component turbomachine is the parent class for pump, compressor and turbine.

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

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    pr : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`

    zeta : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient, :math:`\zeta/\frac{\text{Pa}}{\text{m}^4}`.

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
        return {'pr': dc_cp(min_val=1e-4),
                'zeta': dc_cp(min_val=1e-4),
                'Sirr': dc_cp()}

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
            for i in range(2):
                deriv[0, i, 1] = self.numeric_deriv(self.zeta_func, 'p', i)
                deriv[0, i, 2] = self.numeric_deriv(self.zeta_func, 'h', i)
            if self.zeta.is_var:
                deriv[0, 2 + self.zeta.var_pos, 0] = self.numeric_deriv(self.zeta_func, 'zeta', i)
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
        Returns a starting value for pressure and enthalpy at component's outlet.

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
        Returns a starting value for pressure and enthalpy at component's inlet.

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

    def calc_parameters(self, mode):
        r"""
        Post and preprocessing parameter calculation/specification.

        Parameters
        ----------

        mode : str
            Pre- or postprocessing calculation.

        Note
        ----
        Generic preprocessing is handled by the base class. This method handles class specific pre- and postprocessing.
        """
        component.calc_parameters(self, mode)

        if mode == 'post':
            i = self.inl[0].to_flow()
            o = self.outl[0].to_flow()
            self.pr.val = o[1] / i[1]
            self.zeta.val = (i[1] - o[1]) * math.pi ** 2 / (8 * i[0] ** 2 * (v_mix_ph(i) + v_mix_ph(o)) / 2)
            self.Sirr.val = i[0] * (s_mix_ph(o) - s_mix_ph(i))

# %%


class heat_exchanger_simple(component):
    r"""
    The component heat_exchanger_simple is the parent class for pipe and solar_collector.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        **optional equations**

        - :func:`tespy.components.components.heat_exchanger_simple.Q_func`

        .. math::

            0 = p_{in} \cdot pr - p_{out}

        - :func:`tespy.components.components.component.zeta_func`

        - :func:`tespy.components.components.heat_exchanger_simple.darcy_func` or
          :func:`tespy.components.components.heat_exchanger_simple.hw_func`

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

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Q : Sring/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`.

    zeta : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient, :math:`\zeta/\frac{\text{Pa}}{\text{m}^4}`.

    D : str/float/tespy.helpers.dc_cp
        Diameter of the pipes, :math:`D/\text{m}`.

    L : str/float/tespy.helpers.dc_cp
        Length of the pipes, :math:`L/\text{m}`.

    ks : str/float/tespy.helpers.dc_cp
        Pipes roughness, :math:`ks/\text{m}` for darcy friction,
        :math:`ks/\text{1}` for hazen-williams equation.

    hydro_group : Sring/tespy.helpers.dc_gcp
        Parametergroup for pressure drop calculation based on pipes dimensions.
        Choose 'HW' for hazen-williams equation, else darcy friction factor is used.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient, :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient, provide x and y values
        or use generic values (e. g. calculated from design case).
        Standard method 'HE_COLD', Parameter 'm'.

    Tamb : float/tespy.helpers.dc_cp
        Ambient temperature, provide parameter in network's temperature
        unit.

    Tamb_ref : float/tespy.helpers.dc_cp
         Ambient temperature for reference in offdesign case, provide
         parameter in network's temperature unit.

    kA_group : tespy.helpers.dc_gcp
        Parametergroup for heat transfer calculation from ambient temperature and area
        independent heat transfer coefficient kA.

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
    >>> pi = cmp.pipe('test')
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
                'pr': dc_cp(min_val=1e-4),
                'zeta': dc_cp(min_val=1e-4),
                'D': dc_cp(min_val=1e-2, max_val=2, d=1e-3),
                'L': dc_cp(min_val=1e-1, d=1e-3),
                'ks': dc_cp(min_val=1e-7, max_val=1e-4, d=1e-8),
                'kA': dc_cp(min_val=1, d=1),
                'Tamb': dc_cp(),
                'kA_char': dc_cc(method='HE_HOT', param='m'),
                'SQ1': dc_cp(), 'SQ2': dc_cp(), 'Sirr': dc_cp(),
                'hydro_group': dc_gcp(), 'kA_group': dc_gcp()}

    def inlets(self):
        return ['in1']

    def outlets(self):
        return ['out1']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.fl_deriv = self.fluid_deriv()
        self.m_deriv = self.mass_flow_deriv()

        self.Tamb.val_SI = ((self.Tamb.val + nw.T[nw.T_unit][0]) * nw.T[nw.T_unit][1])
        self.Tamb.design = ((self.Tamb.design + nw.T[nw.T_unit][0]) * nw.T[nw.T_unit][1])

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
            msg = 'Pressure loss calculation from pipe dimensions method is set to ' + method
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
            vec_res += [self.inl[0].p.val_SI * self.pr.val - self.outl[0].p.val_SI]

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
        Calculates vector vec_res with results of additional equations for this component.

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
            for i in range(2):
                zeta_deriv[0, i, 1] = self.numeric_deriv(self.zeta_func, 'p', i)
                zeta_deriv[0, i, 2] = self.numeric_deriv(self.zeta_func, 'h', i)
            # custom variable zeta
            if self.zeta.is_var:
                zeta_deriv[0, 2 + self.zeta.var_pos, 0] = self.numeric_deriv(self.zeta_func, 'zeta', i)
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
            for i in range(2):
                deriv[0, i, 1] = self.numeric_deriv(func, 'p', i)
                deriv[0, i, 2] = self.numeric_deriv(func, 'h', i)
            # custom variables of hydro group
            for var in self.hydro_group.elements:
                if var.is_var:
                    deriv[0, 2 + var.var_pos, 0] = self.numeric_deriv(func, self.vars[var], i)
            mat_deriv += deriv.tolist()

        ######################################################################
        # derivatives for additional equations
        mat_deriv += self.additional_derivatives()

        return np.asarray(mat_deriv)

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
        # derivatives for specified kA-group paremeters
        if self.kA_group.is_set:
            deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
            deriv[0, 0, 0] = self.numeric_deriv(self.kA_func, 'm', 0)
            for i in range(2):
                deriv[0, i, 1] = self.numeric_deriv(self.kA_func, 'p', i)
                deriv[0, i, 2] = self.numeric_deriv(self.kA_func, 'h', i)
            #
            for var in self.kA_group.elements:
                if var.is_var:
                    deriv[0, 2 + var.var_pos, 0] = self.numeric_deriv(self.kA_func, self.vars[var], i)
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

                res = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) - \dot{Q}
        """
        return self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI) - self.Q.val

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

                0 = p_{in} - p_{out} - \frac{8 \cdot |\dot{m}_{in}| \cdot \dot{m}_{in} \cdot
                \frac{v_{in}+v_{out}}{2} \cdot L \cdot \lambda\left(
                Re, ks, D\right)}{\pi^2 \cdot D^5}\\

                \eta: \text{dynamic viscosity}\\
                v: \text{specific volume}\\
                \lambda: \text{darcy friction factor}
        """
        i, o = self.inl[0].to_flow(), self.outl[0].to_flow()

        if abs(i[0]) < 1e-4:
            return i[1] - o[1]

        visc_i, visc_o = visc_mix_ph(i), visc_mix_ph(o)
        v_i, v_o = v_mix_ph(i), v_mix_ph(o)

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

        v_i, v_o = v_mix_ph(i), v_mix_ph(o)
        flow_dir = np.sign(i[0])

        return ((i[1] - o[1]) * flow_dir -
                (10.67 * abs(i[0]) ** 1.852 * self.L.val /
                 (self.ks.val ** 1.852 * self.D.val ** 4.871)) *
                (9.81 * ((v_i + v_o) / 2) ** 0.852))

    def kA_func(self):
        r"""
        Equation for heat transfer calculation from ambient conditions and heat transfer coefficient.

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

        ttd_1 = T_mix_ph(i) - self.Tamb.val_SI
        ttd_2 = T_mix_ph(o) - self.Tamb.val_SI

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
        Returns a starting value for pressure and enthalpy at component's outlet.

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
        Returns a starting value for pressure and enthalpy at component's inlet.

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

    def calc_parameters(self, mode):
        r"""
        Post and preprocessing parameter calculation/specification.

        Parameters
        ----------

        mode : str
            Pre- or postprocessing calculation.

        Note
        ----
        Generic preprocessing is handled by the base class. This method handles class specific pre- and postprocessing.
        """
        component.calc_parameters(self, mode)

        if mode == 'post':
            i = self.inl[0].to_flow()
            o = self.outl[0].to_flow()

            self.SQ1.val = i[0] * (s_mix_ph(o) - s_mix_ph(i))
            self.Q.val = i[0] * (o[2] - i[2])
            self.pr.val = o[1] / i[1]
            self.zeta.val = ((i[1] - o[1]) * math.pi ** 2 / (8 * i[0] ** 2 * (v_mix_ph(i) + v_mix_ph(o)) / 2))

            if self.Tamb.is_set:
                self.SQ2.val = -i[0] * (o[2] - i[2]) / self.Tamb.val_SI
                self.Sirr.val = self.SQ1.val + self.SQ2.val

                ttd_1 = T_mix_ph(i) - self.Tamb.val_SI
                ttd_2 = T_mix_ph(o) - self.Tamb.val_SI

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
                    self.kA_char.func.get_bound_errors(i[0] / self.inl[0].m.design)

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

        - :func:`tespy.components.components.heat_exchanger_simple.darcy_func` or
          :func:`tespy.components.components.heat_exchanger_simple.hw_func`

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

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Q : Sring/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`.

    zeta : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient, :math:`\zeta/\frac{\text{Pa}}{\text{m}^4}`.

    D : str/float/tespy.helpers.dc_cp
        Diameter of the pipes, :math:`D/\text{m}`.

    L : str/float/tespy.helpers.dc_cp
        Length of the pipes, :math:`L/\text{m}`.

    ks : str/float/tespy.helpers.dc_cp
        Pipes roughness, :math:`ks/\text{m}` for darcy friction,
        :math:`ks/\text{1}` for hazen-williams equation.

    hydro_group : Sring/tespy.helpers.dc_gcp
        Parametergroup for pressure drop calculation based on pipes dimensions.
        Choose 'HW' for hazen-williams equation, else darcy friction factor is used.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient, :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient, provide x and y values
        or use generic values (e. g. calculated from design case).
        Standard method 'HE_COLD', Parameter 'm'.

    Tamb : float/tespy.helpers.dc_cp
        Ambient temperature, provide parameter in network's temperature
        unit.

    Tamb_ref : float/tespy.helpers.dc_cp
         Ambient temperature for reference in offdesign case, provide
         parameter in network's temperature unit.

    kA_group : tespy.helpers.dc_gcp
        Parametergroup for heat transfer calculation from ambient temperature and area
        independent heat transfer coefficient kA.

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
    >>> pi = cmp.pipe('test')
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
    The component turbomachine is the parent class for pump, compressor and turbine.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        **optional equations**

        - :func:`tespy.components.components.heat_exchanger_simple.Q_func`

        .. math::

            0 = p_{in} \cdot pr - p_{out}

        - :func:`tespy.components.components.component.zeta_func`

        - :func:`tespy.components.components.heat_exchanger_simple.darcy_func` or
          :func:`tespy.components.components.heat_exchanger_simple.hw_func`

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

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Q : Sring/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`.

    zeta : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient, :math:`\zeta/\frac{\text{Pa}}{\text{m}^4}`.

    D : str/float/tespy.helpers.dc_cp
        Diameter of the pipes, :math:`D/\text{m}`.

    L : str/float/tespy.helpers.dc_cp
        Length of the pipes, :math:`L/\text{m}`.

    ks : str/float/tespy.helpers.dc_cp
        Pipes roughness, :math:`ks/\text{m}` for darcy friction,
        :math:`ks/\text{1}` for hazen-williams equation.

    hydro_group : Sring/tespy.helpers.dc_gcp
        Parametergroup for pressure drop calculation based on pipes dimensions.
        Choose 'HW' for hazen-williams equation, else darcy friction factor is used.

    E :
        Absorption on the inclined surface, :math:`E/\frac{\text{W}}{\text{m}^2}`.

    lkf_lin : str/float/tespy.helpers.dc_cp
        Linear loss key figure, :math:`\alpha_1/\frac{\text{W}}{\text{K} \cdot \text{m}}`.

    lkf_quad : str/float/tespy.helpers.dc_cp
        Quadratic loss key figure, :math:`\alpha_2/\frac{\text{W}}{\text{K}^2 \cdot \text{m}^2}`.

    A : str/float/tespy.helpers.dc_cp
        Collector surface area :math:`A/\text{m}^2`.

    Tamb : float/tespy.helpers.dc_cp
        Ambient temperature, provide parameter in network's temperature unit.

    energy_group : tespy.helpers.dc_gcp
        Parametergroup for energy balance of solarthermal collector.

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
    15.8
    >>> sc.set_attr(A=sc.A.val, E=5e2, Tamb=20)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(sc.Q.val, 1)
    5848.8
    >>> round(outg.T.val, 1)
    69.3
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def component(self):
        return 'solar collector'

    def attr(self):
        return {'Q': dc_cp(),
                'pr': dc_cp(min_val=1e-4),
                'zeta': dc_cp(min_val=1e-4),
                'D': dc_cp(min_val=1e-2, max_val=2, d=1e-3),
                'L': dc_cp(min_val=1e-1, d=1e-3),
                'ks': dc_cp(min_val=1e-7, max_val=1e-4, d=1e-8),
                'E': dc_cp(min_val=0), 'lkf_lin': dc_cp(), 'lkf_quad': dc_cp(),
                'A': dc_cp(min_val=0), 'Tamb': dc_cp(),
                'SQ': dc_cp(),
                'hydro_group': dc_gcp(), 'energy_group': dc_gcp()}

    def inlets(self):
        return ['in1']

    def outlets(self):
        return ['out1']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.fl_deriv = self.fluid_deriv()
        self.m_deriv = self.mass_flow_deriv()

        self.Tamb.val_SI = ((self.Tamb.val + nw.T[nw.T_unit][0]) * nw.T[nw.T_unit][1])

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
        self.energy_group.set_attr(elements=[self.E, self.lkf_lin, self.lkf_quad, self.A, self.Tamb])

        is_set = True
        for e in self.energy_group.elements:
            if not e.is_set:
                is_set = False

        if is_set:
            self.energy_group.set_attr(is_set=True)
        elif self.energy_group.is_set:
            msg = ('All parameters of the component group have to be '
                   'specified! This component group uses the following '
                   'parameters: E, lkf_lin, lkf_quad, A, Tamb at ' + self.label
                   + '. Group will be set to False.')
            logging.info(msg)
            self.energy_group.set_attr(is_set=False)
        else:
            self.energy_group.set_attr(is_set=False)

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for this component.

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
        Calculates matrix of partial derivatives for given additional equations.

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
            for i in range(2):
                deriv[0, i, 1] = self.numeric_deriv(self.energy_func, 'p', i)
                deriv[0, i, 2] = self.numeric_deriv(self.energy_func, 'h', i)
            # custom variables for the energy-group
            for var in self.energy_group.elements:
                if var.is_var:
                    deriv[0, 2 + var.var_pos, 0] = self.numeric_deriv(self.energy_func, self.vars[var], i)
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

                0 = \dot{m} \cdot \left( h_{out} - h_{in} \right) -
                A \cdot \left\{E - \left(T_m - T_{amb} \right) \cdot
                \left[ \alpha_1 + \alpha_2 \cdot A \cdot \left(\
                T_m - T_{amb}\right) \right] \right\}
        """

        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()

        T_m = (T_mix_ph(i) + T_mix_ph(o)) / 2

        return (i[0] * (o[2] - i[2]) - self.A.val * (self.E.val - (T_m - self.Tamb.val_SI) *
                (self.lkf_lin.val + self.lkf_quad.val * self.A.val * (T_m - self.Tamb.val_SI))))

    def calc_parameters(self, mode):
        r"""
        Post and preprocessing parameter calculation/specification.

        Parameters
        ----------

        mode : str
            Pre- or postprocessing calculation.

        Note
        ----
        Generic preprocessing is handled by the base class. This method handles class specific pre- and postprocessing.
        """
        component.calc_parameters(self, mode)

        if mode == 'post':
            i = self.inl[0].to_flow()
            o = self.outl[0].to_flow()

            self.SQ.val = i[0] * (s_mix_ph(o) - s_mix_ph(i))
            self.Q.val = i[0] * (o[2] - i[2])
            self.pr.val = o[1] / i[1]
            self.zeta.val = ((i[1] - o[1]) * math.pi ** 2 / (8 * i[0] ** 2 * (v_mix_ph(i) + v_mix_ph(o)) / 2))

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

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Q : Sring/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at hot side, :math:`\zeta/\frac{\text{Pa}}{\text{m}^4}`.

    zeta2 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at cold side, :math:`\zeta/\frac{\text{Pa}}{\text{m}^4}`.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient, :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char1 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at hot side, provide x and y values
        or use generic values (e. g. calculated from design case). Standard method 'HE_HOT', Parameter 'm'.

    kA_char2 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at cold side, provide x and y values
        or use generic values (e. g. calculated from design case). Standard method 'HE_COLD', Parameter 'm'.

    Note
    ---
    The heat exchanger and subclasses (desuperheater, condenser) are countercurrent heat exchangers.
    Equations (kA, ttd_u, ttd_l) do not work for directcurrent and crosscurrent or combinations of different types.

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
        # derivatives for logarithmic temperature difference not implemented
        return {'Q': dc_cp(), 'kA': dc_cp(), 'td_log': dc_cp(),
                'kA_char1': dc_cc(method='HE_HOT', param='m'),
                'kA_char2': dc_cc(method='HE_COLD', param='m'),
                'ttd_u': dc_cp(), 'ttd_l': dc_cp(),
                'pr1': dc_cp(), 'pr2': dc_cp(),
                'zeta1': dc_cp(), 'zeta2': dc_cp(),
                'SQ1': dc_cp(), 'SQ2': dc_cp(), 'Sirr': dc_cp(),
                'zero_flag': dc_cp(printout=False)}

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
            vec_res += [self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI) - self.Q.val]

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
            vec_res += [self.pr1.val * self.inl[0].p.val_SI - self.outl[0].p.val_SI]

        ######################################################################
        # equations for specified pressure ratio at cold side
        if self.pr2.is_set:
            vec_res += [self.pr2.val * self.inl[1].p.val_SI - self.outl[1].p.val_SI]

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
        Calculates vector vec_res with results of additional equations for this component.

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
        Calculates matrix of partial derivatives for given additional equations.

        Returns
        -------
        mat_deriv : ndarray
            Matrix of partial derivatives.
        """
        return []

    def fluid_func(self):
        r"""
        Calculates the vector of residual values for component's fluid balance equations.

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
        Calculates the residual value for component's mass flow balance equation.

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
            Matrix with partial derivatives for the mass flow balance equations.
        """
        deriv = np.zeros((2, 4 + self.num_vars, self.num_fl + 3))
        for i in range(self.num_i):
            deriv[i, i, 0] = 1
        for j in range(self.num_o):
            deriv[j, j + i + 1, 0] = -1
        return deriv.tolist()

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
        if self.zero_flag.is_set:
            c = self.zero_flag.val
            if c[0] > 0 and c[1] < 3:
                return self.inl[0].m.val_SI

            elif ((c[0] == 0 and c[1] < 3) or
                  (c[0] > 1 and c[1] > 2 and c[1] < 5)):
                return self.outl[0].h.val_SI - self.inl[0].h.val_SI

            elif ((c[0] < 2 and c[1] > 2 and c[1] < 5) or
                  (c[0] == 3 and c[1] == 5)):
                return self.inl[1].m.val_SI
            else:
                return self.outl[1].h.val_SI - self.inl[1].h.val_SI

        else:
            return (self.inl[0].m.val_SI * (self.outl[0].h.val_SI -
                                            self.inl[0].h.val_SI) +
                    self.inl[1].m.val_SI * (self.outl[1].h.val_SI -
                                            self.inl[1].h.val_SI))

    def energy_deriv(self):
        r"""
        Calculates the matrix of partial derivatives for energy balance equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 4, len(self.inl[0].fluid.val) + 3))

        if self.zero_flag.is_set:
            c = self.zero_flag.val
            if c[0] > 0 and c[1] < 3:
                deriv[0, 0, 0] = 1

            elif ((c[0] == 0 and c[1] < 3) or
                  (c[0] > 1 and c[1] > 2 and c[1] < 5)):
                deriv[0, 0, 2] = -1
                deriv[0, 2, 2] = 1

            elif ((c[0] < 2 and c[1] > 2 and c[1] < 5) or
                  (c[0] == 3 and c[1] == 5)):
                deriv[0, 1, 0] = 1
            else:
                deriv[0, 1, 2] = -1
                deriv[0, 3, 2] = 1

        else:
            for k in range(2):
                deriv[0, k, 0] = self.outl[k].h.val_SI - self.inl[k].h.val_SI
                deriv[0, k, 2] = -self.inl[k].m.val_SI

            deriv[0, 2, 2] = self.inl[0].m.val_SI
            deriv[0, 3, 2] = self.inl[1].m.val_SI
        return deriv.tolist()

    def kA_func(self):
        r"""
        Equation for heat transfer from conditions on both sides of heat exchanger.

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
        class :func:`tespy.component.characteristics.characteristics`.

        - Calculate temperatures at inlets and outlets.
        - Perform value manipulation, if temperature levels are not physically feasible.
        """

        if self.zero_flag.is_set:
            c = self.zero_flag.val
            if c[1] == 2 or c[1] == 4 or c[1] == 5:
                T_i1 = T_mix_ph(self.inl[0].to_flow())
                T_i2 = T_mix_ph(self.inl[1].to_flow())
                T_o1 = T_mix_ph(self.outl[0].to_flow())
                T_o2 = T_mix_ph(self.outl[1].to_flow())
                return T_o1 - T_i2 - T_i1 + T_o2

            elif c[0] < 3 and (c[1] == 1 or c[1] == 3):
                return self.outl[1].h.val_SI - self.inl[1].h.val_SI

            elif ((c[0] < 2 and c[1] == 0) or
                  (c[0] == 3 and (c[1] == 1 or c[1] == 3))):
                return self.inl[1].m.val_SI

            else:
                return self.outl[0].h.val_SI - self.inl[0].h.val_SI

        i1 = self.inl[0].to_flow()
        i2 = self.inl[1].to_flow()
        o1 = self.outl[0].to_flow()
        o2 = self.outl[1].to_flow()

        i1_d = self.inl[0].to_flow_design()
        i2_d = self.inl[1].to_flow_design()

        T_i1 = T_mix_ph(i1)
        T_i2 = T_mix_ph(i2)
        T_o1 = T_mix_ph(o1)
        T_o2 = T_mix_ph(o2)

        if T_i1 <= T_o2 and not self.inl[0].T.val_set:
            T_i1 = T_o2 + 0.01
        if T_i1 <= T_o2 and not self.outl[1].T.val_set:
            T_o2 = T_i1 - 0.01
        if T_i1 < T_o2 and self.inl[0].T.val_set and self.outl[1].T.val_set:
            msg = ('Infeasibility at ' + str(self.label) + ': Value for upper '
                   'temperature difference is ' + str(round(T_i1 - T_o2)) + '.')
            logging.error(msg)
            raise ValueError(msg)

        if T_o1 <= T_i2 and not self.outl[0].T.val_set:
            T_o1 = T_i2 + 0.02
        if T_o1 <= T_i2 and not self.inl[1].T.val_set:
            T_i2 = T_o1 - 0.02
        if T_o1 < T_i2 and self.inl[1].T.val_set and self.outl[0].T.val_set:
            msg = ('Infeasibility at ' + str(self.label) + ': Value for lower '
                   'temperature difference is ' + str(round(T_o1 - T_i2)) + '.')
            logging.error(msg)
            raise ValueError(msg)

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

        td_log = (T_o1 - T_i2 - T_i1 + T_o2) / math.log((T_o1 - T_i2) / (T_i1 - T_o2))
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
        i1 = self.inl[0].to_flow()
        o2 = self.outl[1].to_flow()
        return self.ttd_u.val - T_mix_ph(i1) + T_mix_ph(o2)

    def ttd_u_deriv(self):
        r"""
        Calculates the matrix of partial derivatives for upper temperature difference equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 4, len(self.inl[0].fluid.val) + 3))
        for i in range(2):
            deriv[0, i * 3, 1] = self.numeric_deriv(self.ttd_u_func, 'p', i * 3)
            deriv[0, i * 3, 2] = self.numeric_deriv(self.ttd_u_func, 'h', i * 3)
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
        return self.ttd_l.val - T_mix_ph(o1) + T_mix_ph(i2)

    def ttd_l_deriv(self):
        r"""
        Calculates the matrix of partial derivatives for lower temperature difference equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 4, len(self.inl[0].fluid.val) + 3))
        for i in range(2):
            deriv[0, i + 1, 1] = self.numeric_deriv(self.ttd_l_func, 'p', i + 1)
            deriv[0, i + 1, 2] = self.numeric_deriv(self.ttd_l_func, 'h', i + 1)
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
        Manipulate enthalpies/pressure at inlet and outlet if not specified by user to match physically feasible constraints,
        keep fluid composition within feasible range and then propagates it towards the outlet.
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
        Returns a starting value for pressure and enthalpy at component's outlet.

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
        Returns a starting value for pressure and enthalpy at component's inlet.

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

    def calc_parameters(self, mode):
        r"""
        Post and preprocessing parameter calculation/specification.

        Parameters
        ----------

        mode : str
            Pre- or postprocessing calculation.

        Note
        ----
        Generic preprocessing is handled by the base class. This method handles class specific pre- and postprocessing.
        """
        component.calc_parameters(self, mode)

        if mode == 'post':
            # connection information
            i1 = self.inl[0].to_flow()
            i2 = self.inl[1].to_flow()
            o1 = self.outl[0].to_flow()
            o2 = self.outl[1].to_flow()

            # temperatures
            T_i2 = T_mix_ph(i2)
            T_o1 = T_mix_ph(o1)
            T_o2 = T_mix_ph(o2)

            if isinstance(self, condenser):
                T_i1 = T_mix_ph([i1[0], i1[1], h_mix_pQ(i1, 1), i1[3]])
            else:
                T_i1 = T_mix_ph(i1)

            # component parameters
            self.ttd_u.val = T_i1 - T_o2
            self.ttd_l.val = T_o1 - T_i2
            self.Q.val = i1[0] * (o1[2] - i1[2])

            self.pr1.val = o1[1] / i1[1]
            self.pr2.val = o2[1] / i2[1]
            self.zeta1.val = (i1[1] - o1[1]) * math.pi ** 2 / (8 * i1[0] ** 2 * (v_mix_ph(i1) + v_mix_ph(o1)) / 2)
            self.zeta2.val = (i2[1] - o2[1]) * math.pi ** 2 / (8 * i2[0] ** 2 * (v_mix_ph(i2) + v_mix_ph(o2)) / 2)

            self.SQ1.val = self.inl[0].m.val_SI * (s_mix_ph(o1) - s_mix_ph(i1))
            self.SQ2.val = self.inl[1].m.val_SI * (s_mix_ph(o2) - s_mix_ph(i2))
            self.Sirr.val = self.SQ1.val + self.SQ2.val

            # kA and logarithmic temperature difference
            if T_i1 <= T_o2 or T_o1 <= T_i2:
                self.td_log.val = np.nan
                self.kA.val = np.nan
            else:
                self.td_log.val = (T_o1 - T_i2 - T_i1 + T_o2) / math.log((T_o1 - T_i2) / (T_i1 - T_o2))
                self.kA.val = -(i1[0] * (o1[2] - i1[2]) / self.td_log.val)

            if self.ttd_u.val < 0:
                msg = ('Invalid value for terminal temperature difference (upper) '
                       'at component ' + self.label + ': ttd_u = ' + str(self.ttd_u.val) + ' K.')
                logging.error(msg)

            if self.ttd_l.val < 0:
                msg = ('Invalid value for terminal temperature difference (lower) '
                       'at component ' + self.label + ': ttd_l = ' + str(self.ttd_l.val) + ' K.')
                logging.error(msg)

            if self.kA.is_set:
                # get bound errors for kA hot side characteristics
                if self.kA_char1.param == 'm':
                    i1_d = self.inl[0].to_flow_design()
                    if not np.isnan(i1_d[0]):
                        if not i1[0] == 0:
                            self.kA_char1.func.get_bound_errors(i1[0] / i1_d[0])

                # get bound errors for kA copld side characteristics
                if self.kA_char2.param == 'm':
                    i2_d = self.inl[1].to_flow_design()
                    if not np.isnan(i2_d[0]):
                        if not i1[0] == 0:
                            self.kA_char2.func.get_bound_errors(i2[0] / i2_d[0])

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

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Q : Sring/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at hot side, :math:`\zeta/\frac{\text{Pa}}{\text{m}^4}`.

    zeta2 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at cold side, :math:`\zeta/\frac{\text{Pa}}{\text{m}^4}`.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient, :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char1 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at hot side, provide x and y values
        or use generic values (e. g. calculated from design case). Standard method 'COND_HOT', Parameter 'm'.

    kA_char2 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at cold side, provide x and y values
        or use generic values (e. g. calculated from design case). Standard method 'COND_COLD', Parameter 'm'.

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
    ...     h_unit='kJ / kg')
    >>> nw.set_printoptions(print_level='none')
    >>> amb_in = cmp.sink('ambient in')
    >>> amb_out = cmp.source('ambient out')
    >>> hsin = cmp.sink('HS in')
    >>> hsout = cmp.source('HS out')
    >>> he = cmp.condenser('condenser')
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
        return {'Q': dc_cp(), 'kA': dc_cp(), 'td_log': dc_cp(),
                'kA_char1': dc_cc(method='COND_HOT', param='m'),
                'kA_char2': dc_cc(method='COND_COLD', param='m'),
                'ttd_u': dc_cp(), 'ttd_l': dc_cp(),
                'pr1': dc_cp(), 'pr2': dc_cp(),
                'zeta1': dc_cp(), 'zeta2': dc_cp(),
                'SQ1': dc_cp(), 'SQ2': dc_cp(), 'Sirr': dc_cp(),
                'zero_flag': dc_cp()}

    def additional_equations(self):
        r"""
        Calculates vector vec_res with results of additional equations for this component.

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
        Calculates the matrix of partial derivatives for energy balance equation.

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
        class :func:`tespy.component.characteristics.characteristics`.

        - Calculate temperatures at inlets and outlets.
        - Perform value manipulation, if temperature levels are not physically feasible.
        """
        if self.zero_flag.is_set:
            return self.inl[0].p.val_SI - self.inl[0].p.design

        i1 = self.inl[0].to_flow()
        i2 = self.inl[1].to_flow()
        o1 = self.outl[0].to_flow()
        o2 = self.outl[1].to_flow()

        i1_d = self.inl[0].to_flow_design()
        i2_d = self.inl[1].to_flow_design()

        T_i1 = T_mix_ph([i1[0], i1[1], h_mix_pQ(i1, 1), i1[3]])
        T_i2 = T_mix_ph(i2)
        T_o1 = T_mix_ph(o1)
        T_o2 = T_mix_ph(o2)

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

        td_log = (T_o1 - T_i2 - T_i1 + T_o2) / math.log((T_o1 - T_i2) / (T_i1 - T_o2))
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
        The upper terminal temperature difference ttd_u refers to boiling temperature at hot side inlet.
        """
        i1 = self.inl[0].to_flow()
        o2 = self.outl[1].to_flow()
        return (self.ttd_u.val - T_mix_ph([i1[0], i1[1], h_mix_pQ(i1, 1), i1[3]]) + T_mix_ph(o2))

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

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Q : Sring/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : Sring/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at hot side, :math:`\zeta/\frac{\text{Pa}}{\text{m}^4}`.

    zeta2 : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient at cold side, :math:`\zeta/\frac{\text{Pa}}{\text{m}^4}`.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient, :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char1 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at hot side, provide x and y values
        or use generic values (e. g. calculated from design case). Standard method 'COND_HOT', Parameter 'm'.

    kA_char2 : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient at cold side, provide x and y values
        or use generic values (e. g. calculated from design case). Standard method 'COND_COLD', Parameter 'm'.

    Note
    ----

    - The desuperheater has an additional equation for enthalpy at hot side outlet.

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
        Calculates vector vec_res with results of additional equations for this component.

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
        Calculates matrix of partial derivatives for given additional equations.

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

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    Note
    ----
    If you are using a drum in a network with multiple fluids, it is likely
    the fluid propagation causes trouble. If this is the case, try to
    specify the fluid composition at another connection of your network.

    This component assumes, that the fluid composition between outlet 1 and inlet 2 does not change,
    thus there is no equation for the fluid mass fraction at the inlet 2!

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
        vec_res += [h_mix_pQ(self.outl[0].to_flow(), 0) - self.outl[0].h.val_SI]
        vec_res += [h_mix_pQ(self.outl[1].to_flow(), 1) - self.outl[1].h.val_SI]

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
        Calculates the vector of residual values for component's fluid balance equations.

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
        Returns a starting value for pressure and enthalpy at component's outlet.

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
        Returns a starting value for pressure and enthalpy at component's inlet.

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

        - Specify number of inlets and outlets with :code:`num_inter`, predefined value: 1.

    Image

        .. image:: _images/subsys_interface.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the component.

    mode : str
        'auto' for automatic design to offdesign switch, 'man' for manual switch.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    num_inter : float/tespy.helpers.dc_simple
        Number of interfaces for subsystem.

    Note
    ----
    This component passes all fluid properties and mass flow from its inlet to the outlet.

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
        Calculates the partial derivatives for all mass flow, pressure and enthalpy equations.

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
