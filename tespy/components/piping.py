# -*- coding: utf-8

"""Module for piping components.

Components in this module:

    - :func:`tespy.components.piping.pipe`
    - :func:`tespy.components.piping.valve`


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/piping.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.components import component
from tespy.components.heat_exchangers import heat_exchanger_simple

from tespy.tools.data_containers import dc_cc, dc_cp, dc_simple
from tespy.tools.fluid_properties import s_mix_ph, v_mix_ph


# %%


class pipe(heat_exchanger_simple):
    r"""
    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        **optional equations**

        - :func:`tespy.components.heat_exchangers.heat_exchanger_simple.Q_func`

        .. math::

            0 = p_{in} \cdot pr - p_{out}

        - :func:`tespy.components.components.component.zeta_func`

        - :func:`tespy.components.heat_exchangers.heat_exchanger_simple.darcy_func`
          or :func:`tespy.components.heat_exchangers.heat_exchanger_simple.hw_func`

        **additional equations**

        - :func:`tespy.components.heat_exchangers.heat_exchanger_simple.additional_equations`

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

    Q : String/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr : String/float/tespy.helpers.dc_cp
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

    hydro_group : String/tespy.helpers.dc_gcp
        Parametergroup for pressure drop calculation based on pipes dimensions.
        Choose 'HW' for hazen-williams equation, else darcy friction factor is
        used.

    kA : str/float/tespy.helpers.dc_cp
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : str/tespy.helpers.dc_cc
        Characteristic curve for heat transfer coefficient, provide x and y
        values or use generic values (e. g. calculated from design case).
        Standard parameter 'm'.

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
    A mass flow of 10 kg/s ethanol is transported in a pipeline. The pipe is
    considered adiabatic and has a length of 500 meters. We can calculate the
    diameter required at a given pressure loss of 2.5 %. After we determined
    the required diameter, we can predict pressure loss at a different mass
    flow through the pipeline.

    >>> from tespy.components import sink, source, pipe
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> import shutil
    >>> fluid_list = ['ethanol']
    >>> nw = network(fluids=fluid_list)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> so = source('source 1')
    >>> si = sink('sink 1')
    >>> pi = pipe('pipeline')
    >>> pi.component()
    'pipe'
    >>> pi.set_attr(pr=0.975, Q=0, design=['pr'], L=100, D='var', ks=5e-5)
    >>> inc = connection(so, 'out1', pi, 'in1')
    >>> outg = connection(pi, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> inc.set_attr(fluid={'ethanol': 1}, m=10, T=30, p=3)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(pi.D.val, 3)
    0.119
    >>> outg.p.val / inc.p.val == pi.pr.val
    True
    >>> inc.set_attr(m=15)
    >>> pi.set_attr(D=pi.D.val)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(pi.pr.val, 2)
    0.94
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'pipe'

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

        - :func:`tespy.components.piping.valve.dp_char_func`


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

    pr : String/float/tespy.helpers.dc_cp
        Outlet to inlet pressure ratio, :math:`pr/1`

    zeta : str/float/tespy.helpers.dc_cp
        Geometry independent friction coefficient,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    dp_char : str/tespy.helpers.dc_cc
        Characteristic line for difference pressure to mass flow.

    Example
    -------
    A mass flow of 1 kg/s methane is throttled from 80 bar to 15 bar in a
    valve. The inlet temperature is at 50 °C. It is possible to determine the
    outlet temperature as the throttling does not change enthalpy.

    >>> from tespy.components import sink, source, valve
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> import shutil
    >>> fluid_list = ['CH4']
    >>> nw = network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ... iterinfo=False)
    >>> so = source('source')
    >>> si = sink('sink')
    >>> v = valve('valve')
    >>> v.component()
    'valve'
    >>> so_v = connection(so, 'out1', v, 'in1')
    >>> v_si = connection(v, 'out1', si, 'in1')
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
        return {'pr': dc_cp(min_val=1e-4, max_val=1),
                'zeta': dc_cp(min_val=0),
                'dp_char': dc_cc(param='m'),
                'Sirr': dc_simple()}

    @staticmethod
    def inlets():
        return ['in1']

    @staticmethod
    def outlets():
        return ['out1']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        self.fl_deriv = self.fluid_deriv()
        self.m_deriv = self.mass_flow_deriv()
        self.h_deriv = self.enthalpy_deriv()

    def equations(self):
        r"""
        Calculate vector vec_res with results of equations for this component.

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
        # eqation for mass flow
        vec_res += self.mass_flow_func()

        ######################################################################
        # eqation for enthalpy
        vec_res += [self.inl[0].h.val_SI - self.outl[0].h.val_SI]

        ######################################################################
        # eqation for specified pressure ratio
        if self.pr.is_set:
            vec_res += [self.inl[0].p.val_SI * self.pr.val -
                        self.outl[0].p.val_SI]

        ######################################################################
        # eqation specified zeta
        if self.zeta.is_set:
            vec_res += [self.zeta_func()]

        ######################################################################
        # equation for specified difference pressure char
        if self.dp_char.is_set:
            vec_res += [self.dp_char_func()]

        return vec_res

    def derivatives(self):
        r"""
        Calculate matrix of partial derivatives for given equations.

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

        ######################################################################
        # derivatives for specified difference pressure
        if self.dp_char.is_set:
            mat_deriv += self.dp_char_deriv()

        return np.asarray(mat_deriv)

    def enthalpy_deriv(self):
        r"""
        Calculate matrix of partial derivatives for enthalpy balance equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))
        deriv[0, 0, 2] = 1
        deriv[0, 1, 2] = -1
        return deriv.tolist()

    def dp_char_func(self):
        r"""
        Equation for characteristic line defining difference pressure to
        mass flow.

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

    def dp_char_deriv(self):
        r"""
        Calculate the matrix of partial derivatives of the difference pressure
        to mass flow characteristic line.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((1, 2 + self.num_vars, self.num_fl + 3))

        deriv[0, 0, 0] = self.numeric_deriv(self.dp_char_func, 'm', 0)
        deriv[0, 0, 1] = 1
        deriv[0, 1, 1] = -1

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
        self.zeta.val = ((i[1] - o[1]) * np.pi ** 2 /
                         (8 * i[0] ** 2 * (v_mix_ph(i) + v_mix_ph(o)) / 2))
        self.Sirr.val = i[0] * (s_mix_ph(o) - s_mix_ph(i))

        self.check_parameter_bounds()
