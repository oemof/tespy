# -*- coding: utf-8

"""Module of class CombustionChamber.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/combustion/diabatic.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components import CombustionChamber
from tespy.components.component import component_registry
from tespy.tools import logger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.helpers import convert_to_SI


@component_registry
class DiabaticCombustionChamber(CombustionChamber):
    r"""
    The class CombustionChamber is parent class of all combustion components.

    **Mandatory Equations**

    - :py:meth:`tespy.components.combustion.base.CombustionChamber.mass_flow_func`
    - :py:meth:`tespy.components.combustion.base.CombustionChamber.combustion_pressure_func`
    - :py:meth:`tespy.components.combustion.base.CombustionChamber.stoichiometry`

    **Optional Equations**

    - :py:meth:`tespy.components.combustion.base.CombustionChamber.lambda_func`
    - :py:meth:`tespy.components.combustion.base.CombustionChamber.ti_func`
    - :py:meth:`tespy.components.combustion.diabatic.DiabaticCombustionChamber.energy_balance_func`
    - :py:meth:`tespy.components.combustion.diabatic.DiabaticCombustionChamber.pr_func`

    Available fuels

    - methane, ethane, propane, butane, hydrogen

    Inlets/Outlets

    - in1, in2
    - out1

    Image

    .. image:: /api/_images/CombustionChamber.svg
       :alt: flowsheet of the combustion chamber
       :align: center
       :class: only-light

    .. image:: /api/_images/CombustionChamber_darkmode.svg
       :alt: flowsheet of the combustion chamber
       :align: center
       :class: only-dark

    .. note::

        The fuel and the air components can be connected to either of the
        inlets. The pressure of inlet 2 is disconnected from the pressure of
        inlet 1. A warning is prompted, if the pressure at inlet 2 is lower than
        the pressure at inlet 1.

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

    lamb : float, dict
        Actual oxygen to stoichiometric oxygen ratio, :math:`\lambda/1`.

    ti : float, dict
        Thermal input, (:math:`{LHV \cdot \dot{m}_f}`), :math:`ti/\text{W}`.

    eta : float, dict
        Combustion thermal efficiency, :math:`\eta`. Heat loss calculation based
        on share of thermal input.

    pr : float, dict
        Pressure ratio of outlet 1 to inlet 1, :math:`pr`.

    Note
    ----
    For more information on the usage of the combustion chamber see the
    examples section on github or look for the combustion chamber tutorials
    at tespy.readthedocs.io.

    Example
    -------
    The combustion chamber calculates energy input due to combustion as well as
    the flue gas composition based on the type of fuel and the amount of
    oxygen supplied. In this example a mixture of methane, hydrogen and
    carbondioxide is used as fuel.

    >>> from tespy.components import Sink, Source, DiabaticCombustionChamber
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools.fluid_properties import T_sat_p
    >>> import shutil
    >>> nw = Network(p_unit='bar', T_unit='C', iterinfo=False)
    >>> amb = Source('ambient air')
    >>> sf = Source('fuel')
    >>> fg = Sink('flue gas outlet')
    >>> comb = DiabaticCombustionChamber('combustion chamber')
    >>> comb.component()
    'diabatic combustion chamber'
    >>> amb_comb = Connection(amb, 'out1', comb, 'in1')
    >>> sf_comb = Connection(sf, 'out1', comb, 'in2')
    >>> comb_fg = Connection(comb, 'out1', fg, 'in1')
    >>> nw.add_conns(sf_comb, amb_comb, comb_fg)

    Specify the thermal input of the combustion chamber. At the given fluid
    compositions this determines the mass flow of the fuel. The outlet
    temperature of the flue gas determines the ratio of oxygen to fuel mass
    flow. In contrast to the simple combustion chamber, this component does
    allow for a pressure drop. Therefore the outlet pressure or the pressure
    ratio of the combustion chamber must be specified. Since the component is
    not adiabatic, an efficiency value :code:`eta` can be supplied to account
    for heat loss to the ambient. First, we specify :code:`eta=1` and expect
    identical lambda or outlet temperature as in an adiabatic combustion
    chamber.

    >>> comb.set_attr(ti=500000, pr=0.95, eta=1, lamb=1.5)
    >>> amb_comb.set_attr(p=1.2, T=20, fluid={'Ar': 0.0129, 'N2': 0.7553,
    ... 'CO2': 0.0004, 'O2': 0.2314})
    >>> sf_comb.set_attr(T=25, fluid={'CO2': 0.03, 'H2': 0.01, 'CH4': 0.96}, p=1.3)
    >>> nw.solve('design')
    >>> comb_fg.set_attr(T=1200)
    >>> comb.set_attr(lamb=None)
    >>> nw.solve('design')
    >>> round(comb.lamb.val, 3)
    2.014
    >>> round(comb_fg.p.val, 2)
    1.14

    Instead of the pressure ration, we can also specify the outlet pressure.
    The pressure ratio is the ratio or pressure at the outlet to the pressure
    at the inlet 1 (ambient air inlet in this example).

    >>> comb.set_attr(pr=None)
    >>> comb_fg.set_attr(p=1)
    >>> nw.solve('design')
    >>> round(comb.pr.val, 3)
    0.833

    We can change lambda to a specific value and unset the flue gas temperature:

    >>> comb.set_attr(lamb=2)
    >>> comb_fg.set_attr(T=None)
    >>> nw.solve('design')
    >>> round(comb_fg.T.val, 1)
    1206.5

    Now, if we change the efficiency value, e.g. to 0.9, a total of 10 % of
    heat respective to the thermal input will be transferred to the ambient.
    Note, that the heat loss :code:`Q_loss` has a negative value as it is
    extracted from the system.

    >>> eta = 0.9
    >>> comb.set_attr(eta=eta)
    >>> nw.solve('design')
    >>> round(comb.Q_loss.val, 0)
    -50000.0
    >>> round(comb.ti.val * comb.eta.val, 0)
    450000.0
    """

    @staticmethod
    def component():
        return 'diabatic combustion chamber'

    def preprocess(self, num_nw_vars):
        super().preprocess(num_nw_vars)

        if self.dp.is_set:
            self.dp.val_SI = convert_to_SI('p', self.dp.val, self.inl[0].p.unit)

    def get_parameters(self):
        return {
            'lamb': dc_cp(
                min_val=1, deriv=self.lambda_deriv, func=self.lambda_func,
                latex=self.lambda_func_doc, num_eq=1),
            'ti': dc_cp(
                min_val=0, deriv=self.ti_deriv, func=self.ti_func,
                latex=self.ti_func_doc, num_eq=1),
            'pr': dc_cp(
                min_val=0,
                deriv=self.pr_deriv,
                func=self.pr_func,
                latex=self.pr_func_doc,
                num_eq=1,
                func_params={"inconn": 0, "outconn": 0, "pr": "pr"}
                ),
            'dp': dc_cp(
                min_val=0, deriv=self.dp_deriv,
                func=self.dp_func,
                num_eq=1, func_params={"inconn": 0, "outconn": 0, "dp": "dp"}
            ),
            'eta': dc_cp(
                max_val=1, min_val=0, deriv=self.energy_balance_deriv,
                func=self.energy_balance_func,
                latex=self.energy_balance_func_doc, num_eq=1),
            'Q_loss': dc_cp(max_val=0, is_result=True)
        }

    def get_mandatory_constraints(self):
        return {
            k: v for k, v in super().get_mandatory_constraints().items()
            if k in ["mass_flow_constraints", "stoichiometry_constraints"]
        }

    def energy_balance_func(self):
        r"""
        Calculate the energy balance of the diabatic combustion chamber.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \sum_i \dot{m}_{in,i} \cdot
                \left( h_{in,i} - h_{in,i,ref} \right)\\
                & -\dot{m}_{out,2}\cdot\left( h_{out,1}-h_{out,1,ref} \right)\\
                & + LHV_{fuel} \cdot\left(\sum_i\dot{m}_{in,i}\cdot
                x_{fuel,in,i}- \dot{m}_{out,1} \cdot x_{fuel} \right)
                \cdot \eta
                \end{split}\\

                \forall i \in \text{inlets}

        Note
        ----
        The temperature for the reference state is set to 25 °C, thus
        the water may be liquid. In order to make sure, the state is
        referring to the lower heating value, the state of the water in the
        flue gas is fored to gaseous.

        - Reference temperature: 298.15 K.
        - Reference pressure: 1 bar.
        """
        T_ref = 298.15
        p_ref = 1e5

        res = 0
        for i in self.inl:
            i.build_fluid_data()
            res += i.m.val_SI * (
                i.h.val_SI
                - h_mix_pT(p_ref, T_ref, i.fluid_data, mixing_rule="forced-gas")
            )

        for o in self.outl:
            o.build_fluid_data()
            res -= o.m.val_SI * (
                o.h.val_SI
                - h_mix_pT(p_ref, T_ref, o.fluid_data, mixing_rule="forced-gas")
            )

        res += self.calc_ti() * self.eta.val
        return res

    def energy_balance_func_doc(self, label):
        r"""
        Calculate the energy balance of the diabatic combustion chamber.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = (
            r'\begin{split}' + '\n'
            r'0 = & \sum_i \dot{m}_{\mathrm{in,}i} \cdot\left( '
            r'h_{\mathrm{in,}i} - h_{\mathrm{in,}i\mathrm{,ref}} \right) -'
            r'\dot{m}_\mathrm{out,1}\cdot\left( h_\mathrm{out,1}'
            r' - h_\mathrm{out,1,ref}\right)\\' + '\n'
            r'& + LHV_{fuel} \cdot \left(\sum_i \dot{m}_{\mathrm{in,}i} '
            r'\cdot x_{fuel\mathrm{,in,}i} - \dot{m}_\mathrm{out,1} '
            r'\cdot x_{fuel\mathrm{,out,1}} \right) \cdot \eta\\' + '\n'
            r'& \forall i \in \text{inlets}\\'
            r'& T_\mathrm{ref}=\unit[298.15]{K}'
            r'\;p_\mathrm{ref}=\unit[10^5]{Pa}\\'
            '\n' + r'\end{split}'
        )
        return generate_latex_eq(self, latex, label)

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        super().calc_parameters()

        T_ref = 298.15
        p_ref = 1e5

        res = 0
        for i in self.inl:
            i.build_fluid_data()
            res += i.m.val_SI * (
                i.h.val_SI
                - h_mix_pT(p_ref, T_ref, i.fluid_data, mixing_rule="forced-gas")
            )

        for o in self.outl:
            o.build_fluid_data()
            res -= o.m.val_SI * (
                o.h.val_SI
                - h_mix_pT(p_ref, T_ref, o.fluid_data, mixing_rule="forced-gas")
            )

        self.eta.val = -res / self.ti.val
        self.Q_loss.val = -(1 - self.eta.val) * self.ti.val

        self.pr.val = self.outl[0].p.val_SI / self.inl[0].p.val_SI
        self.dp.val_SI = self.inl[0].p.val_SI - self.outl[0].p.val_SI
        self.dp.val = self.inl[0].p.val - self.outl[0].p.val
        for num, i in enumerate(self.inl):
            if i.p.val < self.outl[0].p.val:
                msg = (
                    f"The pressure at inlet {num + 1} is lower than the "
                    f"pressure at the outlet of component {self.label}."
                )
                logger.warning(msg)

    def exergy_balance(self, T0):

        self.E_P = self.outl[0].Ex_physical - (
            self.inl[0].Ex_physical + self.inl[1].Ex_physical
        )
        self.E_F = (
            self.inl[0].Ex_chemical + self.inl[1].Ex_chemical -
            self.outl[0].Ex_chemical
        )

        self.E_D = self.E_F - self.E_P
        self.epsilon = self._calc_epsilon()
        self.E_bus = {"chemical": np.nan, "physical": np.nan, "massless": np.nan}
