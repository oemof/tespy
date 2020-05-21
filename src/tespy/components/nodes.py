# -*- coding: utf-8

"""Module for components of type node.

Components in this module:

- :func:`tespy.components.nodes.droplet_separator`
- :func:`tespy.components.nodes.drum`
- :func:`tespy.components.nodes.node`
- :func:`tespy.components.nodes.merge`
- :func:`tespy.components.nodes.splitter`
- :func:`tespy.components.nodes.separator`


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/nodes.py

SPDX-License-Identifier: MIT
"""

import logging

import numpy as np

from tespy.components.components import component
from tespy.tools.data_containers import dc_simple
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import dh_mix_dpQ
from tespy.tools.fluid_properties import dT_mix_dph
from tespy.tools.fluid_properties import dT_mix_pdh
from tespy.tools.fluid_properties import dT_mix_ph_dfluid
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.helpers import num_fluids

# %%


class node(component):
    r"""
    Class node is the parent class for splitter, separator and merge.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.mass_flow_func`

        .. math::

            0 = p_{in,1} - p_{out,i}\\
            \forall i \in \mathrm{outlets}

        **additional equations**

        - :func:`tespy.components.nodes.splitter.additional_equations`
        - :func:`tespy.components.nodes.separator.additional_equations`
        - :func:`tespy.components.nodes.merge.additional_equations`

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

    design_path: str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings: boolean
        Ignore warnings on default characteristics usage for this component.

    printout: boolean
        Include this component in the network's results printout.

    num_in : float/tespy.tools.data_containers.dc_simple
        Number of inlets for this component, default value: 2.

    num_out : float/tespy.tools.data_containers.dc_simple
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
    The node can serve as merge or as splitter at the same time: The sum of
    all mass flow going into the node is identical to the mass flow leaving it.
    All incoming fluids are mixed in the node (mass flow weighted fluid mass
    fractions and enthalpy). All outgoing fluids have the composition of the
    mixture at the mixtures enthalpy/temperature.

    >>> from tespy.components import sink, source, node
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> import shutil
    >>> import numpy as np
    >>> fluid_list = ['O2', 'N2']
    >>> nw = network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ... h_unit='kJ / kg', iterinfo=False)
    >>> so1 = source('source1')
    >>> so2 = source('source2')
    >>> si1 = sink('sink1')
    >>> si2 = sink('sink2')
    >>> n = node('node', num_in=2, num_out=2)
    >>> n.component()
    'node'
    >>> inc1 = connection(so1, 'out1', n, 'in1')
    >>> inc2 = connection(so2, 'out1', n, 'in2')
    >>> outg1 = connection(n, 'out1', si1, 'in1')
    >>> outg2 = connection(n, 'out2', si2, 'in1')
    >>> nw.add_conns(inc1, inc2, outg1, outg2)

    2 kg/s of pure oxygen is mixed with 5 kg/s mixture (50 % nitrogen and
    50 % oxygen). One of the outgoing connections should be at 3 kg/s.

    >>> inc1.set_attr(fluid={'O2': 1, 'N2': 0}, p=1, T=20, m=2)
    >>> inc2.set_attr(fluid={'O2': 0.5, 'N2': 0.5}, T=50, m=5)
    >>> outg1.set_attr(m=3)
    >>> nw.solve('design')
    >>> (round(outg1.fluid.val['O2'], 3), round(outg1.fluid.val['N2'], 3))
    (0.643, 0.357)
    >>> round(outg1.T.val, 1)
    41.8

    Now calculate at what mass flow of the 50/50 mixture the total
    oxygen fraction of the nodes outlets will be at 80 %.

    >>> inc2.set_attr(m=np.nan)
    >>> outg1.set_attr(fluid={'O2': 0.8})
    >>> nw.solve('design')
    >>> round(inc2.m.val_SI, 3)
    1.333
    """

    @staticmethod
    def component():
        return 'node'

    @staticmethod
    def attr():
        return {'num_in': dc_simple(), 'num_out': dc_simple()}

    def inlets(self):
        if self.num_in.is_set:
            return ['in' + str(i + 1) for i in range(self.num_in.val)]
        else:
            self.set_attr(num_in=2)
            return self.inlets()

    def outlets(self):
        if self.num_out.is_set:
            return ['out' + str(i + 1) for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=2)
            return self.outlets()

    def comp_init(self, nw):

        component.comp_init(self, nw)

        # number of mandatroy equations for
        # mass flow: 1
        # pressure: number of inlets + number of outlets - 1
        # the number of fluid and enthalpy equations is variable for a node!

        self.num_eq = 1 + self.num_i + self.num_o - 1

        # number of fluid equations
        self.num_eq += self.num_nw_fluids * self.num_o
        # number of enthalpy equations
        self.num_eq += self.num_o

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        self.jacobian[0:1] = self.mass_flow_deriv()
        self.jacobian[1:self.num_i + self.num_o] = self.pressure_deriv()

    def equations(self):
        r"""Calculate residual vector with results of equations."""
        k = 0
        ######################################################################
        # eqation for mass flow balance
        self.residual[k] = self.mass_flow_func()
        k += 1

        ######################################################################
        # equations for pressure
        inl = []
        if self.num_i > 1:
            inl = self.inl[1:]
        for c in inl + self.outl:
            self.residual[k] = self.inl[0].p.val_SI - c.p.val_SI
            k += 1

        ######################################################################
        # additional eqations
        self.additional_equations(k)

    def derivatives(self, increment_filter):
        r"""Calculate partial derivatives for given equations."""
        ######################################################################
        # derivatives for mass and pressure are static
        k = self.num_i + self.num_o
        ######################################################################
        # additional derivatives
        self.additional_derivatives(increment_filter, k)

    def additional_equations(self, k):
        r"""
        Calculate the residual vector of the additional equations.

        Equations

            **mandatroy equations**

            .. math::

                0 = \sum_i \left(\dot{m}_{i} \cdot x_{i,j}\right) - x_{o,j}
                \cdot  \sum_i \dot{m}_{i}\\
                \forall j \in \text{fluids}\\
                \forall o \in \text{outgoing mass flows}\\
                \text{i: incoming mass flows}

            .. math::

                0 = \sum_i \dot{m}_{i} \cdot h_{i}
                - \sum_o \dot{m}_{o} \cdot h_{o}
                \forall o \in \text{outgoing mass flows}\\
                \text{i: incoming mass flows}
        """
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

        self.jacobian = self.jacobian[0:k]
        self.jacobian = np.append(self.jacobian, np.zeros((
            len(self.outg) * (self.num_nw_fluids + 1),
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars)), axis=0)

        ######################################################################
        # equations for fluid composition

        for fluid in self.nw_fluids:
            m = 0
            for i in self.inc:
                m += abs(i[0].m.val_SI) * i[0].fluid.val[fluid]
            for o in self.outg:
                self.residual[k] = m - o[0].fluid.val[fluid] * self.m_inc
                k += 1

        ######################################################################
        # equations for energy balance
        for o in self.outg:
            self.residual[k] = h - o[0].h.val_SI * self.m_inc
            k += 1

    def additional_derivatives(self, increment_filter, k):
        r"""Calculate partial derivatives for given additional equations."""
        ######################################################################
        # derivatives for fluid balance equations
        j = 0
        for fluid in self.nw_fluids:
            for o in self.outg:
                self.jacobian[k, o[1], j + 3] = -self.m_inc
                for i in self.inc:
                    self.jacobian[k, i[1], 0] = -i[0].fluid.val[fluid]
                    self.jacobian[k, i[1], j + 3] = -abs(i[0].m.val_SI)
                k += 1
            j += 1

        ######################################################################
        # derivatives for energy balance equations
        for o in self.outg:
            self.jacobian[k, o[1], 2] = -self.m_inc
            for i in self.inc:
                self.jacobian[k, i[1], 0] = i[0].h.val_SI - o[0].h.val_SI
                self.jacobian[k, i[1], 2] = abs(i[0].m.val_SI)
            k += 1

        logging.info(str(self.jacobian))

    def pressure_deriv(self):
        r"""
        Calculate partial derivatives for all pressure equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((
            self.num_i + self.num_o - 1,
            self.num_i + self.num_o,
            self.num_nw_vars))

        inl = []
        if self.num_i > 1:
            inl = self.inl[1:]
        for k in range(len(inl + self.outl)):
            deriv[k, 0, 1] = 1
            deriv[k, k + 1, 1] = -1
        return deriv

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

    @staticmethod
    def initialise_source(c, key):
        r"""
        Return a starting value for pressure and enthalpy at outlet.

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

    @staticmethod
    def initialise_target(c, key):
        r"""
        Return a starting value for pressure and enthalpy at inlet.

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


class droplet_separator(component):
    r"""
    Separate liquid phase from gas phase of a single fluid.

    Equations

        **mandatory equations**

        - :func:`tespy.components.nodes.droplet_separator.fluid_func`
        - :func:`tespy.components.components.component.mass_flow_func`

        .. math::

            0 = \dot{m}_{in} \cdot h_{in} -
            \sum_j \left(\dot{m}_{j,out} \cdot h_{j,out} \right)\\
            \forall j \in outlet

            0 = p_{in} - p_{out,i}\\
            \forall i \in \mathrm{outlets}

            0 = h_{1,out} - h\left(p, x=1 \right)

            0 = h_{2,out} - h\left(p, x=0 \right)\\
            x: \text{vapour mass fraction}

    Inlets/Outlets

        - in1
        - out1, out2 (index 1: gas phase, index 2: liquid phase)

    Image

        .. image:: _images/droplet_separator.svg
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

    design_path: str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings: boolean
        Ignore warnings on default characteristics usage for this component.

    printout: boolean
        Include this component in the network's results printout.

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
    The droplet separator separates gas from liquid phase. From a stream of
    water the liquid phase will be separated.

    >>> from tespy.components import sink, source, droplet_separator
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> from tespy.tools.fluid_properties import Q_ph, T_bp_p
    >>> import shutil
    >>> nw = network(fluids=['water'], T_unit='C', p_unit='bar',
    ... h_unit='kJ / kg', iterinfo=False)
    >>> so = source('two phase inflow')
    >>> sig = sink('gas outflow')
    >>> sil = sink('liquid outflow')
    >>> ds = droplet_separator('droplet separator')
    >>> ds.component()
    'droplet separator'
    >>> so_ds = connection(so, 'out1', ds, 'in1')
    >>> ds_sig = connection(ds, 'out1', sig, 'in1')
    >>> ds_sil = connection(ds, 'out2', sil, 'in1')
    >>> nw.add_conns(so_ds, ds_sig, ds_sil)

    We specify the fluid's state at the inlet. At the gas outflow saturated
    gas enthalpy is expected, at the liquid gas outflow saturated liquid
    enthalpy. The mass flow at the outlets is expected to split according to
    the vapor mass fraction:

    .. math::

        \dot{m}_\mathrm{out,1} = \frac{h_\mathrm{in} - h'}{h'' - h'} \cdot
        \dot{m}_\mathrm{in}

        \dot{m}_\mathrm{out,2} = \left(1 - \frac{h_\mathrm{in} - h'}{h'' - h'}
        \right) \cdot \dot{m}_\mathrm{in}

    >>> so_ds.set_attr(fluid={'water': 1}, p=1, h=1500, m=10)
    >>> nw.solve('design')
    >>> Q_in = Q_ph(so_ds.p.val_SI, so_ds.h.val_SI, 'water')
    >>> round(Q_in * so_ds.m.val_SI, 6) == round(ds_sig.m.val_SI, 6)
    True
    >>> round((1 - Q_in) * so_ds.m.val_SI, 6) == round(ds_sil.m.val_SI, 6)
    True
    >>> Q_ph(ds_sig.p.val_SI, ds_sig.h.val_SI, 'water')
    1.0
    >>> Q_ph(ds_sil.p.val_SI, ds_sil.h.val_SI, 'water')
    0.0

    In a different setup, we unset pressure and enthalpy and specify gas
    temperature and mass flow instead. The temperature specification must yield
    the corresponding boiling point pressure and the mass flow must yield the
    inlet enthalpy. The inlet vapor mass fraction must be equal to fraction of
    gas mass flow to inlet mass flow (0.95 in this example).

    >>> so_ds.set_attr(fluid={'water': 1}, p=None, h=None, T=150, m=10)
    >>> ds_sig.set_attr(m=9.5)
    >>> nw.solve('design')
    >>> round(Q_ph(so_ds.p.val_SI, so_ds.h.val_SI, 'water'), 6)
    0.95
    >>> T_boil = T_bp_p(so_ds.to_flow())
    >>> round(T_boil, 6) == round(so_ds.T.val_SI, 6)
    True
    """

    @staticmethod
    def component():
        return 'droplet separator'

    @staticmethod
    def inlets():
        return ['in1']

    @staticmethod
    def outlets():
        return ['out1', 'out2']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        # number of mandatroy equations for
        # fluid balance: num_fl * 2
        # mass flow: 1
        # pressure: 2
        # enthalpy: 1
        # saturated liquid outlet: 1
        # saturated gas outlet: 1
        self.num_eq = self.num_nw_fluids * 2 + 6

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        pos = self.num_nw_fluids * 2
        self.jacobian[0:pos] = self.fluid_deriv()
        self.jacobian[pos:pos + 1] = self.mass_flow_deriv()
        self.jacobian[pos + 1:pos + 3] = self.pressure_deriv()

    def equations(self):
        r"""Calculate residual vector with results of equations."""
        k = 0
        ######################################################################
        # eqations for fluid balance
        self.residual[k:k + self.num_nw_fluids * 2] = self.fluid_func()
        k += self.num_nw_fluids * 2

        ######################################################################
        # eqations for mass flow balance
        self.residual[k] = self.mass_flow_func()
        k += 1

        ######################################################################
        # eqations for pressure
        p = self.inl[0].p.val_SI
        for c in self.outl:
            self.residual[k] = p - c.p.val_SI
            k += 1

        ######################################################################
        # eqations for enthalpy
        val = 0
        for i in self.inl:
            val += i.m.val_SI * i.h.val_SI
        for o in self.outl:
            val -= o.m.val_SI * o.h.val_SI
        self.residual[k] = val
        k += 1

        ######################################################################
        # eqations for staturated fluid state at outlets
        self.residual[k] = h_mix_pQ(
            self.outl[0].to_flow(), 1) - self.outl[0].h.val_SI
        k += 1
        self.residual[k] = h_mix_pQ(
            self.outl[1].to_flow(), 0) - self.outl[1].h.val_SI
        k += 1

    def derivatives(self, increment_filter):
        r"""Calculate partial derivatives for given equations."""
        ######################################################################
        # derivatives fluid, mass flow and pressure balance are static
        k = self.num_nw_fluids * 2 + 3

        ######################################################################
        # derivatives for energy balance equation
        self.jacobian[k, 0, 0] = self.inl[0].h.val_SI
        self.jacobian[k, 0, 2] = self.inl[0].m.val_SI
        j = 0
        for outl in self.outl:
            self.jacobian[k, j + 1, 0] = -outl.h.val_SI
            self.jacobian[k, j + 1, 2] = -outl.m.val_SI
            j += 1
        k += 1

        ######################################################################
        # derivatives of equations for saturated states at outlets
        self.jacobian[k, 1, 1] = dh_mix_dpQ(self.outl[0].to_flow(), 0)
        self.jacobian[k, 1, 2] = -1
        k += 1
        self.jacobian[k, 2, 1] = dh_mix_dpQ(self.outl[1].to_flow(), 1)
        self.jacobian[k, 2, 2] = -1
        k += 1

    def fluid_func(self):
        r"""
        Calculate the vector of residual values for fluid balance equations.

        Returns
        -------
        residual : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = fluid_{i,in_1} - fluid_{i,out_{j}}\\
                \forall i \in \mathrm{fluid}, \; \forall j \in inlets

        """
        residual = []

        for o in self.outl:
            for fluid, x in self.inl[0].fluid.val.items():
                residual += [x - o.fluid.val[fluid]]
        return residual

    def fluid_deriv(self):
        r"""
        Calculate partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((2 * self.num_nw_fluids, 3, self.num_nw_vars))
        for k in range(2):
            for i in range(self.num_nw_fluids):
                deriv[i + k * self.num_nw_fluids, 0, i + 3] = 1
                deriv[i + k * self.num_nw_fluids, k + 1, i + 3] = -1
        return deriv

    def pressure_deriv(self):
        r"""
        Calculate partial derivatives for pressure equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the pressure equations.
        """
        deriv = np.zeros((2, 3, self.num_nw_vars))
        for k in range(2):
            deriv[k, 0, 1] = 1
            deriv[k, k + 1, 1] = -1
        return deriv

    @staticmethod
    def initialise_source(c, key):
        r"""
        Return a starting value for pressure and enthalpy at outlet.

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
                h\left(p, x=1 \right) & \text{key = 'h' at outlet 1}\\
                h\left(p, x=0 \right) & \text{key = 'h' at outlet 2}
                \end{cases}
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            if c.source_id == 'out1':
                return h_mix_pQ(c.to_flow(), 1)
            else:
                return h_mix_pQ(c.to_flow(), 0)

    @staticmethod
    def initialise_target(c, key):
        r"""
        Return a starting value for pressure and enthalpy at inlet.

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
                h\left(p, x=0.5 \right) & \text{key = 'h' at inlet 1}
                \end{cases}
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            return h_mix_pQ(c.to_flow(), 0.5)

# %%


class drum(component):
    r"""
    A drum separates saturated gas from saturated liquid.

    Equations

        **mandatory equations**

        - :func:`tespy.components.nodes.drum.fluid_func`
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

    design_path: str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings: boolean
        Ignore warnings on default characteristics usage for this component.

    printout: boolean
        Include this component in the network's results printout.

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
    The drum separates saturated gas from saturated liquid. The liquid phase is
    transported to an evaporator, the staturated gas phase is extracted from
    the drum. In this example ammonia is evaporated using ambient air. A
    characteristic function is applied for the heat transfer coefficient of the
    evaporator. It is possible to load the char_line with the function
    :code:`load_default_char` from the default lines. We want to use the
    'EVAPORATING FLUID' lines of the heat exchanger.

    >>> from tespy.components import sink, source, drum, pump, heat_exchanger
    >>> from tespy.connections import connection, ref
    >>> from tespy.networks import network
    >>> from tespy.tools.characteristics import char_line
    >>> from tespy.tools.characteristics import load_default_char as ldc
    >>> import shutil
    >>> import numpy as np
    >>> nw = network(fluids=['NH3', 'air'], T_unit='C', p_unit='bar',
    ...     h_unit='kJ / kg', iterinfo=False)
    >>> fa = source('feed ammonia')
    >>> amb_in = source('air inlet')
    >>> amb_out = sink('air outlet')
    >>> s = sink('steam')
    >>> dr = drum('drum')
    >>> dr.component()
    'drum'
    >>> ev = heat_exchanger('evaporator')
    >>> erp = pump('evaporator reciculation pump')
    >>> f_dr = connection(fa, 'out1', dr, 'in1')
    >>> dr_erp = connection(dr, 'out1', erp, 'in1')
    >>> erp_ev = connection(erp, 'out1', ev, 'in2')
    >>> ev_dr = connection(ev, 'out2', dr, 'in2')
    >>> dr_s = connection(dr, 'out2', s, 'in1')
    >>> nw.add_conns(f_dr, dr_erp, erp_ev, ev_dr, dr_s)
    >>> amb_ev = connection(amb_in, 'out1', ev, 'in1')
    >>> ev_amb = connection(ev, 'out1', amb_out, 'in1')
    >>> nw.add_conns(amb_ev, ev_amb)

    The ambient air enters the evaporator at 30 °C. The pinch point temperature
    difference (ttd_l) of the evaporator is at 5 K, and 1 MW of heat should be
    transferred. State of ammonia at the inlet is at -5 °C and 5 bar. From this
    design it is possible to calculate offdesign performance at 75 % part load.

    >>> char1 = ldc('heat exchanger', 'kA_char1', 'DEFAULT',
    ... char_line)
    >>> char2 = ldc('heat exchanger', 'kA_char2', 'EVAPORATING FLUID',
    ... char_line)
    >>> ev.set_attr(pr1=0.999, pr2=0.99, ttd_l=5, kA_char1=char1,
    ...     kA_char2=char2, design=['pr1', 'ttd_l'],
    ...     offdesign=['zeta1', 'kA_char'])
    >>> ev.set_attr(Q=-1e6)
    >>> erp.set_attr(eta_s=0.8)
    >>> f_dr.set_attr(p=5, T=-5)
    >>> erp_ev.set_attr(m=ref(f_dr, 4, 0), fluid={'air': 0, 'NH3': 1})
    >>> amb_ev.set_attr(fluid={'air': 1, 'NH3': 0}, T=30)
    >>> ev_amb.set_attr(p=1)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(ev_amb.T.val - erp_ev.T.val ,1)
    5.0
    >>> round(f_dr.h.val, 1)
    322.7
    >>> round(dr_erp.h.val, 1)
    364.9
    >>> round(ev_dr.h.val, 1)
    687.2
    >>> round(f_dr.m.val, 2)
    0.78
    >>> ev.set_attr(Q=-0.75e6)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(f_dr.m.val, 2)
    0.58
    >>> round(ev_amb.T.val - erp_ev.T.val ,1)
    3.0
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'drum'

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1', 'out2']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        # number of mandatroy equations for
        # fluid balance: num_fl * 2
        # mass flow: 1
        # pressure: 3
        # enthalpy: 1
        # saturated liquid outlet: 1
        # saturated gas outlet: 1
        self.num_eq = self.num_nw_fluids * 2 + 7

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        pos = self.num_nw_fluids * 2
        self.jacobian[0:pos] = self.fluid_deriv()
        self.jacobian[pos:pos + 1] = self.mass_flow_deriv()
        self.jacobian[pos + 1:pos + 4] = self.pressure_deriv()

    def equations(self):
        r"""Calculate residual vector with results of equations."""
        k = 0
        ######################################################################
        # eqations for fluid balance
        self.residual[k:k + self.num_nw_fluids * 2] = self.fluid_func()
        k += self.num_nw_fluids * 2

        ######################################################################
        # eqations for mass flow balance
        self.residual[k] = self.mass_flow_func()
        k += 1

        ######################################################################
        # eqations for pressure
        p = self.inl[0].p.val_SI
        for c in [self.inl[1]] + self.outl:
            self.residual[k] = p - c.p.val_SI
            k += 1

        ######################################################################
        # eqations for enthalpy
        val = 0
        for i in self.inl:
            val += i.m.val_SI * i.h.val_SI
        for o in self.outl:
            val -= o.m.val_SI * o.h.val_SI
        self.residual[k] = val
        k += 1

        ######################################################################
        # eqations for staturated fluid state at outlets
        self.residual[k] = h_mix_pQ(
            self.outl[0].to_flow(), 0) - self.outl[0].h.val_SI
        k += 1
        self.residual[k] = h_mix_pQ(
            self.outl[1].to_flow(), 1) - self.outl[1].h.val_SI
        k += 1

    def derivatives(self, increment_filter):
        r"""Calculate partial derivatives for given equations."""
        ######################################################################
        # derivatives fluid, mass flow and pressure balance are static
        k = self.num_nw_fluids * 2 + 4

        ######################################################################
        # derivatives for energy balance equation
        i = 0
        for inl in self.inl:
            self.jacobian[k, i, 0] = inl.h.val_SI
            self.jacobian[k, i, 2] = inl.m.val_SI
            i += 1
        j = 0
        for outl in self.outl:
            self.jacobian[k, j + i, 0] = -outl.h.val_SI
            self.jacobian[k, j + i, 2] = -outl.m.val_SI
            j += 1
        k += 1

        ######################################################################
        # derivatives of equations for saturated states at outlets
        self.jacobian[k, 2, 1] = dh_mix_dpQ(self.outl[0].to_flow(), 0)
        self.jacobian[k, 2, 2] = -1
        k += 1
        self.jacobian[k, 3, 1] = dh_mix_dpQ(self.outl[1].to_flow(), 1)
        self.jacobian[k, 3, 2] = -1
        k += 1

    def fluid_func(self):
        r"""
        Calculate the vector of residual values for fluid balance equations.

        Returns
        -------
        residual : list
            Vector of residual values for component's fluid balance.

            .. math::

                0 = fluid_{i,in_1} - fluid_{i,out_{j}}\\
                \forall i \in \mathrm{fluid}, \; \forall j \in inlets

        """
        residual = []

        for o in self.outl:
            for fluid, x in self.inl[0].fluid.val.items():
                residual += [x - o.fluid.val[fluid]]
        return residual

    def fluid_deriv(self):
        r"""
        Calculate partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((2 * self.num_nw_fluids, 4, self.num_nw_vars))
        for k in range(2):
            for i in range(self.num_nw_fluids):
                deriv[i + k * self.num_nw_fluids, 0, i + 3] = 1
                deriv[i + k * self.num_nw_fluids, k + 2, i + 3] = -1
        return deriv

    def pressure_deriv(self):
        r"""
        Calculate partial derivatives for pressure equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the pressure equations.
        """
        deriv = np.zeros((3, 4, self.num_nw_vars))
        for k in range(3):
            deriv[k, 0, 1] = 1
            deriv[k, k + 1, 1] = -1
        return deriv

    @staticmethod
    def initialise_source(c, key):
        r"""
        Return a starting value for pressure and enthalpy at outlet.

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
            if c.source_id == 'out1':
                return h_mix_pQ(c.to_flow(), 0)
            else:
                return h_mix_pQ(c.to_flow(), 1)

    @staticmethod
    def initialise_target(c, key):
        r"""
        Return a starting value for pressure and enthalpy at inlet.

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
            if c.target_id == 'in1':
                return h_mix_pQ(c.to_flow(), 0)
            else:
                return h_mix_pQ(c.to_flow(), 0.7)

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

    design_path: str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings: boolean
        Ignore warnings on default characteristics usage for this component.

    printout: boolean
        Include this component in the network's results printout.

    num_in : float/tespy.tools.data_containers.dc_simple
        Number of inlets for this component, default value: 2.

    Example
    -------
    The merge mixes a specified number of mass flows and has a single outlet.
    At the outlet, fluid composition and enthalpy are calculated by mass
    weighted fluid composition and enthalpy of the inlets.

    >>> from tespy.components import sink, source, merge
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> import shutil
    >>> import numpy as np
    >>> fluid_list = ['O2', 'N2']
    >>> nw = network(fluids=fluid_list, p_unit='bar', iterinfo=False)
    >>> so1 = source('source1')
    >>> so2 = source('source2')
    >>> so3 = source('source3')
    >>> si1 = sink('sink')
    >>> m = merge('merge', num_in=3)
    >>> m.component()
    'merge'
    >>> inc1 = connection(so1, 'out1', m, 'in1')
    >>> inc2 = connection(so2, 'out1', m, 'in2')
    >>> inc3 = connection(so3, 'out1', m, 'in3')
    >>> outg = connection(m, 'out1', si1, 'in1')
    >>> nw.add_conns(inc1, inc2, inc3, outg)

    A merge with three inlets mixes air (simplified) with pure nitrogen and
    pure oxygen. All gases enter the component at the same temperature. As
    mixing effects are not considered, the outlet temperature should thus be
    similar to the three inlet temperatures (difference might occur due to
    rounding in fluid property functions, let's check it for two different
    temperatures). It is e.g. possible to find the required mass flow of pure
    nitrogen given the nitrogen mass fraction in the outlet.

    >>> T = 293.15
    >>> inc1.set_attr(fluid={'O2': 0.23, 'N2': 0.77}, p=1, T=T, m=5)
    >>> inc2.set_attr(fluid={'O2': 1, 'N2':0}, T=T, m=5)
    >>> inc3.set_attr(fluid={'O2': 0, 'N2': 1}, T=T)
    >>> outg.set_attr(fluid={'N2': 0.4})
    >>> nw.solve('design')
    >>> round(inc3.m.val_SI, 2)
    0.25
    >>> abs(round((outg.T.val_SI - T) / T, 5)) < 0.01
    True
    >>> T = 173.15
    >>> inc1.set_attr(T=T)
    >>> inc2.set_attr(T=T)
    >>> inc3.set_attr(T=T)
    >>> nw.solve('design')
    >>> abs(round((outg.T.val_SI - T) / T, 5)) < 0.01
    True
    """

    @staticmethod
    def component():
        return 'merge'

    @staticmethod
    def attr():
        return {'num_in': dc_simple()}

    def inlets(self):
        if self.num_in.is_set:
            return ['in' + str(i + 1) for i in range(self.num_in.val)]
        else:
            self.set_attr(num_in=2)
            return self.inlets()

    @staticmethod
    def outlets():
        return ['out1']

    def comp_init(self, nw):

        component.comp_init(self, nw)

        # number of mandatroy equations for
        # mass flow: 1
        # pressure: number of inlets + number of outlets - 1
        # fluid: number of fluids
        # enthalpy: 1

        self.num_eq = self.num_i + self.num_o + self.num_nw_fluids + 1

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        self.jacobian[0:1] = self.mass_flow_deriv()
        end = self.num_i + self.num_o
        self.jacobian[1:end] = self.pressure_deriv()

    def additional_equations(self, k):
        r"""
        Calculate results of additional equations.

        Equations

            **mandatroy equations**

            .. math::

                0 = \dot{m}_{in_{j}} \cdot fluid_{i,in_{j}} -
                    \dot {m}_{out} \cdot fluid_{i,out} \\
                \forall i \in \mathrm{fluid}, \; \forall j \in inlets

            .. math::

                0 = h_{in} - h_{out,i} \;
                \forall i \in \mathrm{outlets}\\
        """
        ######################################################################
        # equations for fluid balance
        for fluid, x in self.outl[0].fluid.val.items():
            res = -x * self.outl[0].m.val_SI
            for i in self.inl:
                res += i.fluid.val[fluid] * i.m.val_SI
            self.residual[k] = res
            k += 1

        ######################################################################
        # equation for energy balance
        h_res = -self.outl[0].m.val_SI * self.outl[0].h.val_SI
        for i in self.inl:
            h_res += i.m.val_SI * i.h.val_SI
        self.residual[k] = h_res
        k += 1

    def additional_derivatives(self, increment_filter, k):
        r"""Calculate partial derivatives for given additional equations."""
        ######################################################################
        # derivatives for fluid balance equations
        i = 0
        for fluid, x in self.outl[0].fluid.val.items():
            j = 0
            for inl in self.inl:
                self.jacobian[k, j, 0] = inl.fluid.val[fluid]
                self.jacobian[k, j, i + 3] = inl.m.val_SI
                j += 1
            self.jacobian[k, j, 0] = -x
            self.jacobian[k, j, i + 3] = -self.outl[0].m.val_SI
            i += 1
            k += 1

        ######################################################################
        # derivatives for energy balance equations
        self.jacobian[k, self.num_i, 0] = -self.outl[0].h.val_SI
        self.jacobian[k, self.num_i, 2] = -self.outl[0].m.val_SI
        j = 0
        for i in self.inl:
            self.jacobian[k, j, 0] = i.h.val_SI
            self.jacobian[k, j, 2] = i.m.val_SI
            j += 1
        k += 1

# %%


class separator(node):
    r"""
    A separator separates fluid components from a mass flow.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.mass_flow_func`

        .. math::

            0 = p_{in} - p_{out,i} \;
            \forall i \in \mathrm{outlets}

        **additional equations**

        - :func:`tespy.components.nodes.separator.additional_equations`

    Inlets/Outlets

        - in1
        - specify number of outlets with :code:`num_out` (default value: 2)

    Image

        .. image:: _images/split.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Note
    ----
    Fluid separation requires power and cooling, equations have not been
    implemented, yet!

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

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings: boolean
        Ignore warnings on default characteristics usage for this component.

    printout: boolean
        Include this component in the network's results printout.

    num_out : float/tespy.tools.data_containers.dc_simple
        Number of outlets for this component, default value: 2.

    Example
    -------
    The separator is used to split up a single mass flow into a specified
    number of different parts at identical pressure and temperature but
    different fluid composition. Fluids can be separated from each other.

    >>> from tespy.components import sink, source, separator
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> import shutil
    >>> import numpy as np
    >>> fluid_list = ['O2', 'N2']
    >>> nw = network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ... iterinfo=False)
    >>> so = source('source')
    >>> si1 = sink('sink1')
    >>> si2 = sink('sink2')
    >>> s = separator('separator', num_out=2)
    >>> s.component()
    'separator'
    >>> inc = connection(so, 'out1', s, 'in1')
    >>> outg1 = connection(s, 'out1', si1, 'in1')
    >>> outg2 = connection(s, 'out2', si2, 'in1')
    >>> nw.add_conns(inc, outg1, outg2)

    An Air (simplified) mass flow of 5 kg/s is split up into two mass flows.
    One mass flow of 1 kg/s containing 10 % oxygen and 90 % nitrogen leaves the
    separator. It is possible to calculate the fluid composition of the second
    mass flow. Specify starting values for the second mass flow fluid
    composition for calculation stability.

    >>> inc.set_attr(fluid={'O2': 0.23, 'N2': 0.77}, p=1, T=20, m=5)
    >>> outg1.set_attr(fluid={'O2': 0.1, 'N2': 0.9}, m=1)
    >>> outg2.set_attr(fluid0={'O2': 0.5, 'N2': 0.5})
    >>> nw.solve('design')
    >>> outg2.fluid.val['O2']
    0.2625

    In the same way, it is possible to specify one of the fluid components in
    the second mass flow instead of the first mass flow. The solver will find
    the mass flows matching the desired composition. 65 % of the mass flow
    will leave the separator at the second outlet the case of 30 % oxygen
    mass fraction for this outlet.

    >>> outg1.set_attr(m=np.nan)
    >>> outg2.set_attr(fluid={'O2': 0.3})
    >>> nw.solve('design')
    >>> outg2.fluid.val['O2']
    0.3
    >>> round(outg2.m.val_SI / inc.m.val_SI, 2)
    0.65
    """

    @staticmethod
    def component():
        return 'separator'

    @staticmethod
    def attr():
        return {'num_out': dc_simple()}

    @staticmethod
    def inlets():
        return ['in1']

    def outlets(self):
        if self.num_out.is_set:
            return ['out' + str(i + 1) for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=2)
            return self.outlets()

    def comp_init(self, nw):

        component.comp_init(self, nw)

        # number of mandatroy equations for
        # mass flow: 1
        # pressure: number of inlets + number of outlets - 1
        # fluid: number of fluids
        # enthalpy: number of outlets

        self.num_eq = self.num_i + self.num_o * 2 + self.num_nw_fluids

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        self.jacobian[0:1] = self.mass_flow_deriv()
        end = self.num_i + self.num_o
        self.jacobian[1:end] = self.pressure_deriv()

    def additional_equations(self, k):
        r"""
        Calculate results of additional equations.

        Equations

            **mandatroy equations**

            .. math:: 0 = fluid_{i,in} - fluid_{i,out_{j}} \;
                \forall i \in \mathrm{fluid}, \; \forall j \in outlets

            .. math::

                0 = T_{in} - T_{out,i} \;
                \forall i \in \mathrm{outlets}
        """
        ######################################################################
        # equations for fluid balance
        for fluid, x in self.inl[0].fluid.val.items():
            res = x * self.inl[0].m.val_SI
            for o in self.outl:
                res -= o.fluid.val[fluid] * o.m.val_SI
            self.residual[k] = res
            k += 1

        ######################################################################
        # equations for energy balance
        for o in self.outl:
            self.residual[k] = (
                T_mix_ph(self.inl[0].to_flow(), T0=self.inl[0].T.val_SI) -
                T_mix_ph(o.to_flow(), T0=o.T.val_SI))
            k += 1

    def additional_derivatives(self, increment_filter, k):
        r"""Calculate partial derivatives for given additional equations."""
        ######################################################################
        # derivatives for fluid balance equations
        i = 0
        for fluid in self.nw_fluids:
            j = 0
            for o in self.outl:
                self.jacobian[k, j + 1, 0] = -o.fluid.val[fluid]
                self.jacobian[k, j + 1, i + 3] = -o.m.val_SI
                j += 1
            self.jacobian[k, 0, 0] = self.inl[0].fluid.val[fluid]
            self.jacobian[k, 0, i + 3] = self.inl[0].m.val_SI
            k += 1
            i += 1

        ######################################################################
        # derivatives for energy balance equations
        i = self.inl[0].to_flow()
        j = 0
        for o in self.outl:
            o = o.to_flow()
            self.jacobian[k, 0, 1] = dT_mix_dph(i)
            self.jacobian[k, 0, 2] = dT_mix_pdh(i)
            self.jacobian[k, 0, 3:] = dT_mix_ph_dfluid(i)
            self.jacobian[k, j + 1, 1] = -dT_mix_dph(o)
            self.jacobian[k, j + 1, 2] = -dT_mix_pdh(o)
            self.jacobian[k, j + 1, 3:] = -1 * dT_mix_ph_dfluid(o)
            j += 1
            k += 1

    @staticmethod
    def initialise_fluids(nw):
        r"""
        Fluid initialisation for fluid mixture at outlet of the node.

        Parameters
        ----------
        nw : tespy.networks.network
            Network using this component object.
        """
        return

# %%


class splitter(node):
    r"""
    Split up a mass flow in parts of identical enthalpy and fluid composition.

    Equations

        **mandatory equations**

        - :func:`tespy.components.components.component.mass_flow_func`

        .. math::

            0 = p_{in} - p_{out,i} \;
            \forall i \in \mathrm{outlets}

        **additional equations**

        - :func:`tespy.components.nodes.splitter.additional_equations`

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

    design_path: str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings: boolean
        Ignore warnings on default characteristics usage for this component.

    printout: boolean
        Include this component in the network's results printout.

    num_out : float/tespy.tools.data_containers.dc_simple
        Number of outlets for this component, default value: 2.

    Example
    -------
    A splitter is used to split up a single mass flow into a specified number
    of different parts at identical pressure, enthalpy and fluid composition.

    >>> from tespy.components import sink, source, splitter
    >>> from tespy.connections import connection
    >>> from tespy.networks import network
    >>> import shutil
    >>> import numpy as np
    >>> fluid_list = ['O2', 'N2']
    >>> nw = network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ... iterinfo=False)
    >>> so = source('source')
    >>> si1 = sink('sink1')
    >>> si2 = sink('sink2')
    >>> si3 = sink('sink3')
    >>> s = splitter('splitter', num_out=3)
    >>> s.component()
    'splitter'
    >>> inc = connection(so, 'out1', s, 'in1')
    >>> outg1 = connection(s, 'out1', si1, 'in1')
    >>> outg2 = connection(s, 'out2', si2, 'in1')
    >>> outg3 = connection(s, 'out3', si3, 'in1')
    >>> nw.add_conns(inc, outg1, outg2, outg3)

    An Air (simplified) mass flow is split up into three mass flows. The total
    incoming mass flow is 5 kg/s, 3 kg/s and 1 kg/s respectively are leaving
    the splitter into the first two outlets. The residual mass flow will
    drain in the last outlet. Temperature and fluid composition will not
    change.

    >>> inc.set_attr(fluid={'O2': 0.23, 'N2': 0.77}, p=1, T=20, m=5)
    >>> outg1.set_attr(m=3)
    >>> outg2.set_attr(m=1)
    >>> nw.solve('design')
    >>> round(outg3.m.val_SI, 1)
    1.0
    >>> round(inc.T.val, 1)
    20.0
    >>> round(outg3.T.val, 1)
    20.0
    """

    @staticmethod
    def component():
        return 'splitter'

    @staticmethod
    def attr():
        return {'num_out': dc_simple()}

    @staticmethod
    def inlets():
        return ['in1']

    def outlets(self):
        if self.num_out.is_set:
            return ['out' + str(i + 1) for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=2)
            return self.outlets()

    def comp_init(self, nw):

        component.comp_init(self, nw)

        # number of mandatroy equations for
        # mass flow: 1
        # pressure: number of inlets + number of outlets - 1
        # fluid: number of outlets * number of fluid_set
        # enthalpy: number of outlets

        self.num_eq = self.num_i + self.num_o * (2 + self.num_nw_fluids)

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        self.jacobian[0:1] = self.mass_flow_deriv()
        end = self.num_i + self.num_o
        self.jacobian[1:end] = self.pressure_deriv()
        start = end
        end = start + self.num_o * self.num_nw_fluids
        self.jacobian[start:end] = self.fluid_deriv()
        start = end
        end = start + self.num_o
        self.jacobian[start:end] = self.enthalpy_deriv()

    def additional_equations(self, k):
        r"""
        Calculate results of additional equations.

        Equations

            **mandatroy equations**

            .. math:: 0 = fluid_{i,in} - fluid_{i,out_{j}} \;
                \forall i \in \mathrm{fluid}, \; \forall j \in outlets

            .. math::
                0 = h_{in} - h_{out,i} \;
                \forall i \in \mathrm{outlets}\\
        """
        ######################################################################
        # equations for fluid balance
        for o in self.outl:
            for fluid, x in self.inl[0].fluid.val.items():
                self.residual[k] = x - o.fluid.val[fluid]
                k += 1

        ######################################################################
        # equations for energy balance
        for o in self.outl:
            self.residual[k] = self.inl[0].h.val_SI - o.h.val_SI
            k += 1

    def additional_derivatives(self, increment_filter, k):
        r"""Calculate partial derivatives for given additional equations."""
        return

    def fluid_deriv(self):
        r"""
        Calculate partial derivatives for all fluid balance equations.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        deriv = np.zeros((self.num_nw_fluids * self.num_o, 1 + self.num_o,
                          self.num_nw_vars))
        k = 0
        for o in self.outl:
            i = 0
            for fluid in self.nw_fluids:
                deriv[i + k * self.num_nw_fluids, 0, i + 3] = 1
                deriv[i + k * self.num_nw_fluids, k + 1, i + 3] = -1
                i += 1
            k += 1
        return deriv

    def enthalpy_deriv(self):
        r"""
        Calculate partial derivatives for enthalpy balance equation.

        Returns
        -------
        deriv : list
            Matrix of partial derivatives.
        """
        deriv = np.zeros((self.num_o, 1 + self.num_o, self.num_nw_vars))
        k = 0
        for o in self.outl:
            deriv[k, 0, 2] = 1
            deriv[k, k + 1, 2] = -1
            k += 1
        return deriv

    @staticmethod
    def initialise_fluids(nw):
        r"""
        Fluid initialisation for fluid mixture at outlet of the node.

        Parameters
        ----------
        nw : tespy.networks.network
            Network using this component object.
        """
        return
