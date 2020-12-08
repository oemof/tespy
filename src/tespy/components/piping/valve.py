# -*- coding: utf-8

"""Module of class Valve.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/piping.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import Component
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.global_vars import err


class Valve(Component):
    r"""
    The Valve throttles a fluid without changing enthalpy.

    Equations

        **mandatory equations**

        - :py:meth:`tespy.components.component.Component.fluid_func`
        - :py:meth:`tespy.components.component.Component.mass_flow_func`

        .. math::

            0 = h_{in} - h_{out}

        **optional equations**

        .. math::

            0 = p_{in} \cdot pr - p_{out}

        - :py:meth:`tespy.components.component.Component.zeta_func`

        - :py:meth:`tespy.components.piping.valve.Valve.dp_char_func`


    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/Valve.svg
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

    pr : str, float, tespy.tools.data_containers.ComponentProperties
        Outlet to inlet pressure ratio, :math:`pr/1`

    zeta : str, float, tespy.tools.data_containers.ComponentProperties
        Geometry independent friction coefficient,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    dp_char : tespy.tools.characteristics.CharLine, tespy.tools.data_containers.ComponentCharacteristics
        Characteristic line for difference pressure to mass flow.

    Example
    -------
    A mass flow of 1 kg/s methane is throttled from 80 bar to 15 bar in a
    valve. The inlet temperature is at 50 °C. It is possible to determine the
    outlet temperature as the throttling does not change enthalpy.

    >>> from tespy.components import Sink, Source, Valve
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import shutil
    >>> fluid_list = ['CH4']
    >>> nw = Network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ... iterinfo=False)
    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> v = Valve('valve')
    >>> v.component()
    'valve'
    >>> so_v = Connection(so, 'out1', v, 'in1')
    >>> v_si = Connection(v, 'out1', si, 'in1')
    >>> nw.add_conns(so_v, v_si)
    >>> v.set_attr(offdesign=['zeta'])
    >>> so_v.set_attr(fluid={'CH4': 1}, m=1, T=50, p=80, design=['m'])
    >>> v_si.set_attr(p=15)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(v_si.T.val, 1)
    26.3
    >>> round(v.pr.val, 3)
    0.188

    The simulation determined the area independant zeta value
    :math:`\frac{\zeta}{D^4}`. This zeta remains constant if the cross
    sectional area of the valve opening does not change. Using the zeta value
    we can determine the pressure ratio at a different feed pressure.

    >>> so_v.set_attr(p=70)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(so_v.m.val, 1)
    0.9
    >>> round(v_si.T.val, 1)
    30.0
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'valve'

    @staticmethod
    def attr():
        return {
            'pr': dc_cp(min_val=1e-4, max_val=1),
            'zeta': dc_cp(min_val=0, max_val=1e15),
            'dp_char': dc_cc(param='m')
        }

    @staticmethod
    def inlets():
        return ['in1']

    @staticmethod
    def outlets():
        return ['out1']

    def comp_init(self, nw):

        Component.comp_init(self, nw)

        # number of mandatroy equations for
        # fluid balance: num_fl
        # mass flow: 1
        # enthalpy: 1
        self.num_eq = self.num_nw_fluids + 2
        for var in [self.pr, self.zeta, self.dp_char]:
            if var.is_set:
                self.num_eq += 1

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        pos = self.num_nw_fluids
        self.jacobian[0:pos] = self.fluid_deriv()
        self.jacobian[pos:pos + 1] = self.mass_flow_deriv()
        self.jacobian[pos + 1:pos + 2] = self.enthalpy_deriv()

    def equations(self):
        r"""Calculate residual vector with results of equations."""
        k = 0
        ######################################################################
        # eqations for fluids
        self.residual[k:k + self.num_nw_fluids] = self.fluid_func()
        k += self.num_nw_fluids

        ######################################################################
        # eqation for mass flow
        self.residual[k] = self.mass_flow_func()
        k += 1

        ######################################################################
        # eqation for enthalpy
        self.residual[k] = self.inl[0].h.val_SI - self.outl[0].h.val_SI
        k += 1

        ######################################################################
        # eqation for specified pressure ratio
        if self.pr.is_set:
            self.residual[k] = (
                self.inl[0].p.val_SI * self.pr.val - self.outl[0].p.val_SI)
            k += 1

        ######################################################################
        # eqation specified zeta
        if self.zeta.is_set:
            if (np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0 or
                    self.always_all_equations):
                self.residual[k] = self.zeta_func(zeta='zeta')
            k += 1

        ######################################################################
        # equation for specified difference pressure char
        if self.dp_char.is_set:
            if (np.absolute(self.residual[k]) > err ** 2 or self.it % 4 == 0 or
                    self.always_all_equations):
                self.residual[k] = self.dp_char_func()
            k += 1

    def derivatives(self, increment_filter):
        r"""Calculate partial derivatives for given equations."""
        ######################################################################
        # derivatives fluid, mass flow and enthalpy balance are static
        k = self.num_nw_fluids + 2

        ######################################################################
        # derivatives for specified pressure ratio
        if self.pr.is_set:
            self.jacobian[k, 0, 1] = self.pr.val
            self.jacobian[k, 1, 1] = -1
            if self.pr.is_var:
                self.jacobian[k, 2 + self.pr.var_pos, 0] = (
                    self.inl[0].p.val_SI)
            k += 1

        ######################################################################
        # derivatives for specified zeta
        if self.zeta.is_set:
            f = self.zeta_func
            if not increment_filter[0, 0]:
                self.jacobian[k, 0, 0] = self.numeric_deriv(
                    f, 'm', 0, zeta='zeta')
            if not increment_filter[0, 1]:
                self.jacobian[k, 0, 1] = self.numeric_deriv(
                    f, 'p', 0, zeta='zeta')
            if not increment_filter[0, 2]:
                self.jacobian[k, 0, 2] = self.numeric_deriv(
                    f, 'h', 0, zeta='zeta')
            if not increment_filter[1, 1]:
                self.jacobian[k, 1, 1] = self.numeric_deriv(
                    f, 'p', 1, zeta='zeta')
            if not increment_filter[1, 2]:
                self.jacobian[k, 1, 2] = self.numeric_deriv(
                    f, 'h', 1, zeta='zeta')
            if self.zeta.is_var:
                self.jacobian[k, 2 + self.zeta.var_pos, 0] = (
                    self.numeric_deriv(f, 'zeta', 2, zeta='zeta'))
            k += 1

        ######################################################################
        # derivatives for specified difference pressure
        if self.dp_char.is_set:
            if not increment_filter[0, 0]:
                self.jacobian[k, 0, 0] = self.numeric_deriv(
                    self.dp_char_func, 'm', 0)
            self.jacobian[k, 0, 1] = 1
            self.jacobian[k, 1, 1] = -1
            k += 1

    def enthalpy_deriv(self):
        r"""
        Calculate matrix of partial derivatives for enthalpy balance equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 2 + self.num_vars, self.num_nw_vars))
        deriv[0, 0, 2] = 1
        deriv[0, 1, 2] = -1
        return deriv.tolist()

    def dp_char_func(self):
        r"""
        Equation for characteristic line of difference pressure to mass flow.

        Returns
        -------
        res : ndarray
            Residual value of equation.

            .. math::

                res =  p_1 - p_2 - f \left( \dot{m} \right)
        """
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()

        return i[1] - o[1] - self.dp_char.func.evaluate(i[0])

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
                5 \cdot 10^5 & \text{key = 'p'}\\
                5 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            return 5e5

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        i = self.inl[0].to_flow()
        o = self.outl[0].to_flow()
        self.pr.val = o[1] / i[1]
        self.zeta.val = ((i[1] - o[1]) * np.pi ** 2 / (
            4 * i[0] ** 2 * (self.inl[0].vol.val_SI + self.outl[0].vol.val_SI)
            ))
        self.check_parameter_bounds()

    def entropy_balance(self):
        r"""
        Calculate entropy balance of a valve.

        Note
        ----
        The entropy balance makes the follwing parameter available:

        .. math::

            \text{S\_irr}=\dot{m} \cdot \left(s_\mathrm{out}-s_\mathrm{in}
            \right)\\
        """
        self.S_irr = self.inl[0].m.val_SI * (
            self.outl[0].s.val_SI - self.inl[0].s.val_SI)

    def exergy_balance(self, Tamb):
        r"""
        Calculate exergy balance of a valve.

        Note
        ----
         .. math::

            \dot{E}_\mathrm{P} = \text{not defined (nan)}\\
            \dot{E}_\mathrm{F} = \dot{m}_\mathrm{in} \cdot \left(
            e_\mathrm{ph,in} - e_\mathrm{ph,out}\right)
        """
        self.E_P = np.nan
        self.E_F = self.inl[0].Ex_physical - self.outl[0].Ex_physical
        self.E_D = self.E_F
        self.epsilon = np.nan

    def get_plotting_data(self):
        """Generate a dictionary containing FluProDia plotting information.

        Returns
        -------
        data : dict
            A nested dictionary containing the keywords required by the
            :code:`calc_individual_isoline` method of the
            :code:`FluidPropertyDiagram` class. First level keys are the
            connection index ('in1' -> 'out1', therefore :code:`1` etc.).
        """
        return {
            1: {
                'isoline_property': 'h',
                'isoline_value': self.inl[0].h.val,
                'isoline_value_end': self.outl[0].h.val,
                'starting_point_property': 'v',
                'starting_point_value': self.inl[0].vol.val,
                'ending_point_property': 'v',
                'ending_point_value': self.outl[0].vol.val
            }
        }
