# -*- coding: utf-8

"""Module of class Desuperheater.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/desuperheater.py

SPDX-License-Identifier: MIT
"""
import logging
import warnings

import numpy as np

from tespy.components.component import Component
from tespy.components.heat_exchangers.heat_exchanger import HeatExchanger
from tespy.tools.fluid_properties import dh_mix_dpQ
from tespy.tools.fluid_properties import h_mix_pQ


class Desuperheater(HeatExchanger):
    r"""
    The Desuperheater cools a fluid to the saturated gas state.

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger.HeatExchanger.energy_balance_func`
    - :py:meth:`tespy.components.heat_exchangers.desuperheater.Desuperheater.saturated_gas_func`

    **Optional Equations**

    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger.HeatExchanger.energy_balance_hot_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger.HeatExchanger.kA_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger.HeatExchanger.kA_char_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger.HeatExchanger.ttd_u_func`
    - :py:meth:`tespy.components.heat_exchangers.heat_exchanger.HeatExchanger.ttd_l_func`
    - hot side :py:meth:`tespy.components.component.Component.pr_func`
    - cold side :py:meth:`tespy.components.component.Component.pr_func`
    - hot side :py:meth:`tespy.components.component.Component.zeta_func`
    - cold side :py:meth:`tespy.components.component.Component.zeta_func`

    Inlets/Outlets

    - in1, in2 (index 1: hot side, index 2: cold side)
    - out1, out2 (index 1: hot side, index 2: cold side)

    Image

    .. image:: _images/HeatExchanger.svg
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

    Q : str, float, tespy.tools.data_containers.ComponentProperties
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : str, float, tespy.tools.data_containers.ComponentProperties
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : str, float, tespy.tools.data_containers.ComponentProperties
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : str, float, tespy.tools.data_containers.ComponentProperties
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : str, float, tespy.tools.data_containers.ComponentProperties
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    kA : str, float, tespy.tools.data_containers.ComponentProperties
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char1 : tespy.tools.characteristics.CharLine, tespy.tools.data_containers.ComponentCharacteristics
        Characteristic line for hot side heat transfer coefficient.

    kA_char2 : tespy.tools.characteristics.CharLine, tespy.tools.data_containers.ComponentCharacteristics
        Characteristic line for cold side heat transfer coefficient.

    Note
    ----
    The desuperheater has an additional equation for enthalpy at hot side
    outlet: The fluid leaves the component in saturated gas state.

    Example
    -------
    Overheated enthanol is cooled with water in a heat exchanger until it
    reaches the state of saturated gas.

    >>> from tespy.components import Sink, Source, Desuperheater
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools.fluid_properties import T_bp_p
    >>> import shutil
    >>> nw = Network(fluids=['water', 'ethanol'], T_unit='C', p_unit='bar',
    ... h_unit='kJ / kg', v_unit='l / s', m_range=[0.001, 10], iterinfo=False)
    >>> et_in = Source('ethanol inlet')
    >>> et_out = Sink('ethanol outlet')
    >>> cw_in = Source('cooling water inlet')
    >>> cw_out = Sink('cooling water outlet')
    >>> desu = Desuperheater('desuperheater')
    >>> desu.component()
    'desuperheater'
    >>> et_de = Connection(et_in, 'out1', desu, 'in1')
    >>> de_et = Connection(desu, 'out1', et_out, 'in1')
    >>> cw_de = Connection(cw_in, 'out1', desu, 'in2')
    >>> de_cw = Connection(desu, 'out2', cw_out, 'in1')
    >>> nw.add_conns(et_de, de_et, cw_de, de_cw)

    The cooling water enters the component at 15 °C. 10 l/s of ethanol is
    cooled from 100 K above boiling point. The water flow rate is at 1 l/s.
    Knowing the component's design parameters it is possible to predict
    behavior at different inlet temperatures or different volumetric flow of
    ethanol. Controlling the ethanol's state at the outlet is only possible,
    if the cooling water flow rate is adjusted accordingly.

    >>> desu.set_attr(pr1=0.99, pr2=0.98, design=['pr1', 'pr2'],
    ... offdesign=['zeta1', 'zeta2', 'kA_char'])
    >>> cw_de.set_attr(fluid={'water': 1, 'ethanol': 0}, T=15, v=1,
    ... design=['v'])
    >>> de_cw.set_attr(p=1)
    >>> et_de.set_attr(fluid={'water': 0, 'ethanol': 1}, Td_bp=100, v=10)
    >>> de_et.set_attr(p=1)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(de_cw.T.val, 1)
    15.5
    >>> round(de_et.x.val, 1)
    1.0
    >>> et_de.set_attr(v=12)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(cw_de.v.val, 2)
    1.94
    >>> et_de.set_attr(v=7)
    >>> nw.solve('offdesign', design_path='tmp', init_path='tmp')
    >>> round(cw_de.v.val, 2)
    0.41
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'desuperheater'

    def comp_init(self, nw):

        # number of mandatroy equations for
        # fluid balance: num_fl * 2
        # mass flow: 2
        # energy balance: 1
        # enthalpy hot side outlet: 1
        num_eq = len(nw.fluids) * 2 + 4
        Component.comp_init(self, nw, num_eq=num_eq)
        # constant derivatives
        pos = self.num_nw_fluids * 2
        self.jacobian[0:pos] = self.fluid_deriv()
        self.jacobian[pos:pos + 2] = self.mass_flow_deriv()

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
        k = HeatExchanger.mandatory_equations(self, doc=doc)
        ######################################################################
        # equation for saturated gas at hot side outlet
        self.residual[k] = self.saturated_gas_func()
        if doc:
            self.equation_docs[k:k + 1] = self.saturated_gas_func(doc=doc)
        k += 1
        return k

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
        k = HeatExchanger.mandatory_derivatives(self, increment_filter)
        ######################################################################
        # derivatives for saturated gas at hot side outlet equation
        o1 = self.outl[0].to_flow()
        self.jacobian[k, 2, 1] = -dh_mix_dpQ(o1, 1)
        self.jacobian[k, 2, 2] = 1
        k += 1
        return k

    def saturated_gas_func(self, doc=False):
        r"""
        Calculate hot side outlet state.

        Parameters
        ----------
        doc : boolean
            Return equation in LaTeX format instead of value.

        Returns
        -------
        residual : float
            Residual value of equation

            .. math::

                0 = h_{out,1} - h\left(p_{out,1}, x=1 \right)
        """
        if not doc:
            return self.outl[0].h.val_SI - h_mix_pQ(self.outl[0].to_flow(), 1)
        else:
            latex = r'0=h_\mathrm{out,1}-h\left(p_\mathrm{out,1}, x=1 \right)'
            return [self.generate_latex(latex, 'saturated_gas_func')]
