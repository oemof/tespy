# -*- coding: utf-8

"""Module of class Desuperheater.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/desuperheater.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.heat_exchangers.base import HeatExchanger
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.fluid_properties import dh_mix_dpQ
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import single_fluid


@component_registry
class Desuperheater(HeatExchanger):
    r"""
    The Desuperheater cools a fluid to the saturated gas state.

    **Mandatory Equations**

    - fluid: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - mass flow: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_func`
    - :py:meth:`tespy.components.heat_exchangers.desuperheater.Desuperheater.saturated_gas_func`

    **Optional Equations**

    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.kA_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.kA_char_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.ttd_u_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.ttd_l_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.ttd_min_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.eff_cold_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.eff_hot_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.eff_max_func`

    For hot and cold side individually:

    - :py:meth:`tespy.components.component.Component.pr_structure_matrix`
    - :py:meth:`tespy.components.component.Component.dp_structure_matrix`
    - :py:meth:`tespy.components.component.Component.zeta_func`

    Inlets/Outlets

    - in1, in2 (index 1: hot side, index 2: cold side)
    - out1, out2 (index 1: hot side, index 2: cold side)

    Image

    .. image:: /api/_images/HeatExchanger.svg
       :alt: flowsheet of the desuperheater
       :align: center
       :class: only-light

    .. image:: /api/_images/HeatExchanger_darkmode.svg
       :alt: flowsheet of the desuperheater
       :align: center
       :class: only-dark

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

    Q : float, dict
        Heat transfer, :math:`Q/\text{W}`.

    pr1 : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    zeta1 : float, dict, :code:`"var"`
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : float, dict, :code:`"var"`
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    ttd_l : float, dict
        Lower terminal temperature difference :math:`ttd_\mathrm{l}/\text{K}`.

    ttd_u : float, dict
        Upper terminal temperature difference :math:`ttd_\mathrm{u}/\text{K}`.

    kA : float, dict
        Area independent heat transfer coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char1 : tespy.tools.characteristics.CharLine, dict
        Characteristic line for hot side heat transfer coefficient.

    kA_char2 : tespy.tools.characteristics.CharLine, dict
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
    >>> import os
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC", "enthalpy": "kJ/kg",
    ...     "volumetric_flow": "l/s"
    ... })
    >>> et_in = Source('ethanol inlet')
    >>> et_out = Sink('ethanol outlet')
    >>> cw_in = Source('cooling water inlet')
    >>> cw_out = Sink('cooling water outlet')
    >>> desu = Desuperheater('desuperheater')
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

    >>> desu.set_attr(
    ...     pr1=0.99, pr2=0.98, design=['pr1', 'pr2'],
    ...     offdesign=['zeta1', 'zeta2', 'kA_char']
    ... )
    >>> cw_de.set_attr(fluid={'water': 1}, T=15, v=1, design=['v'])
    >>> de_cw.set_attr(p=1)
    >>> et_de.set_attr(fluid={'ethanol': 1}, td_dew=100, v=10)
    >>> de_et.set_attr(p=1)
    >>> nw.solve('design')
    >>> nw.save('tmp.json')
    >>> round(de_cw.T.val, 1)
    15.5
    >>> round(de_et.x.val, 1)
    1.0
    >>> et_de.set_attr(v=12)
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(cw_de.v.val, 2)
    1.94
    >>> et_de.set_attr(v=7)
    >>> nw.solve('offdesign', init_path='tmp.json', design_path='tmp.json')
    >>> round(cw_de.v.val, 2)
    0.41
    >>> os.remove('tmp.json')
    """

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        constraints.update({
            'saturated_gas_constraints': dc_cmc(**{
                'num_eq_sets': 1,
                'func': self.saturated_gas_func,
                'deriv': self.saturated_gas_deriv,
                'dependents': self.saturated_gas_dependents
            })
        })
        return constraints

    def saturated_gas_func(self):
        r"""
        Calculate hot side outlet state.

        Returns
        -------
        residual : float
            Residual value of equation

            .. math::

                0 = h_{out,1} - h\left(p_{out,1}, x=1 \right)
        """
        o = self.outl[0]
        return o.h.val_SI - h_mix_pQ(o.p.val_SI, 1, o.fluid_data)

    def saturated_gas_deriv(self, increment_filter, k, dependents=None):
        r"""
        Partial derivatives of saturated gas at hot side outlet function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        o = self.outl[0]
        self._partial_derivative(o.h, k, 1, increment_filter)
        if o.p.is_var:
            self._partial_derivative(
                o.p, k, -dh_mix_dpQ(o.p.val_SI, 1, o.fluid_data), increment_filter
            )

    def saturated_gas_dependents(self):
        return [
            self.outl[0].p,
            self.outl[0].h
        ]

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
        """
        if c.source_id == 'out1':
            if key == 'h':
                return h_mix_pQ(c.p.val_SI, 1, c.fluid_data, c.mixing_rule)

        return super().initialise_source(c, key)
