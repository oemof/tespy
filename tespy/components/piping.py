# -*- coding: utf-8

"""
.. module:: turbomachine
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

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

    Q : Sring/float/tespy.helpers.dc_cp
        Heat transfer, :math:`Q/\text{W}`.

    pr : Sring/float/tespy.helpers.dc_cp
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

    hydro_group : Sring/tespy.helpers.dc_gcp
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

    pr : Sring/float/tespy.helpers.dc_cp
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
