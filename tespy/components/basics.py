# -*- coding: utf-8

"""
.. module:: turbomachine
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

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
