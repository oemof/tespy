# -*- coding: utf-8

"""Module of class DiabaticCombustionChamber.


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
from tespy.tools.fluid_properties import h_mix_pT


@component_registry
class DiabaticCombustionChamber(CombustionChamber):
    r"""
    An extension of the adiabatic combustion chamber with pressure drop and heat loss.

    .. image:: /api/_images/components/CombustionChamber.svg
       :alt: flowsheet of the diabaticcombustionchamber
       :align: center
       :class: only-light

    .. image:: /api/_images/components/CombustionChamber_darkmode.svg
       :alt: flowsheet of the diabaticcombustionchamber
       :align: center
       :class: only-dark

    Ports
    -----

    - Fluid inlets: in1, in2
    - Fluid outlets: out1

    Mandatory Equations
    -------------------

    - mass flow balance over all inflows and outflows: :py:meth:`mass_flow_func <tespy.components.combustion.base.CombustionChamber.mass_flow_func>`
    - constraints for stoichiometry of the reaction: :py:meth:`stoichiometry_func <tespy.components.combustion.base.CombustionChamber.stoichiometry_func>`

    Parameters
    ----------

    char_warnings : bool
        Ignore warnings on default characteristics usage for this component.

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the components design case.

    dp : float, dict
        Inlet 0 to outlet 0 absolute pressure change. Quantity:
        :code:`pressure_difference`.
        Equation: :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`.

    eta : float, dict
        Heat dissipation ratio relative to thermal input. Quantity:
        :code:`efficiency`.
        Equation: :py:meth:`energy_balance_func <tespy.components.combustion.diabatic.DiabaticCombustionChamber.energy_balance_func>`.

    f_nox : float, dict
        Mass-based nitric oxide (NO) generation rate in flue gas in mass of
        created NO per mass of fuel and air input. Only active if value is
        explicitly set. Quantity: :code:`ratio`.

    label : str
        The label of the component.

    lamb : float, dict
        Available oxygen to stoichiometric oxygen ratio. Quantity:
        :code:`ratio`.
        Equation: :py:meth:`lambda_func <tespy.components.combustion.base.CombustionChamber.lambda_func>`.

    local_design : bool
        Treat this component in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this component in offdesign mode in a design calculation.

    offdesign : list
        List containing offdesign parameters (stated as String).

    pr : float, dict
        Outlet 0 to inlet 0 pressure ratio. Quantity: :code:`ratio`.
        Equation: :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`.

    printout : bool
        Include this component in the network's results printout.

    Qloss : float, dict
        Heat dissipation. Quantity: :code:`heat`.

    ti : float, dict
        Thermal input of fuel: lower heating value multiplied with mass flow.
        Quantity: :code:`heat`.
        Equation: :py:meth:`ti_func <tespy.components.combustion.base.CombustionChamber.ti_func>`.

    Notes
    -----

    .. tip::

        You can add more fluids by importing :code:`COMBUSTION_FLUIDS` from
        the :code:`tespy.tools` module and passing the respective information.
        See in the example of
        :py:class:`tespy.components.combustion.base.CombustionChamber`, how to
        do that. To retrieve the fluids available by default run:

        .. code-block:: python

            from tespy.tools.global_vars import COMBUSTION_FLUIDS
            COMBUSTION_FLUIDS.fluids.keys()

    .. note::

        The fuel and the air components can be connected to either of the
        inlets. The pressure of inlet 2 is disconnected from the pressure of
        inlet 1. A warning is prompted, if the pressure at inlet 2 is lower than
        the pressure at inlet 1.

    Example
    -------
    The combustion chamber calculates energy input due to combustion as well as
    the flue gas composition based on the type of fuel and the amount of
    oxygen supplied. In this example a mixture of methane, hydrogen and
    carbon dioxide is used as fuel.

    >>> from tespy.components import Sink, Source, DiabaticCombustionChamber
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools.fluid_properties import T_sat_p
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC"
    ... })
    >>> amb = Source('ambient air')
    >>> sf = Source('fuel')
    >>> fg = Sink('flue gas outlet')
    >>> comb = DiabaticCombustionChamber('combustion chamber')
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

    >>> comb.set_attr(ti=500000, pr=0.95, eta=1)
    >>> amb_comb.set_attr(
    ...     p=1.2, T=20,
    ...     fluid={'Ar': 0.0129, 'N2': 0.7553, 'CO2': 0.0004, 'O2': 0.2314}
    ... )
    >>> sf_comb.set_attr(
    ...     p=1.3, T=25, fluid={'CO2': 0.03, 'H2': 0.01, 'CH4': 0.96}
    ... )
    >>> comb_fg.set_attr(T=1200)
    >>> nw.solve('design')
    >>> round(comb.lamb.val, 3)
    2.014
    >>> round(comb_fg.p.val, 2)
    1.14

    Instead of the pressure ratio, we can also specify the outlet pressure.
    The pressure ratio is the ratio of pressure at the outlet to the pressure
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
    Note, that the heat loss :code:`Qloss` has a negative value as it is
    extracted from the system.

    >>> eta = 0.9
    >>> comb.set_attr(eta=eta)
    >>> nw.solve('design')
    >>> round(comb.Qloss.val, 0)
    -50000.0
    >>> round(comb.ti.val * comb.eta.val, 0)
    450000.0
    """
    def get_parameters(self):
        params = super().get_parameters()
        params.update({
            'pr': dc_cp(
                min_val=0,
                num_eq_sets=1,
                structure_matrix=self.pr_structure_matrix,
                func_params={"inconn": 0, "outconn": 0, "pr": "pr"},
                quantity="ratio",
                description="outlet 0 to inlet 0 pressure ratio",
                calc=self._calc_pr
            ),
            'dp': dc_cp(
                min_val=0,
                num_eq_sets=1,
                structure_matrix=self.dp_structure_matrix,
                func_params={"inconn": 0, "outconn": 0, "dp": "dp"},
                quantity="pressure_difference",
                description="inlet 0 to outlet 0 absolute pressure change",
                calc=self._calc_dp
            ),
            'eta': dc_cp(
                max_val=1, min_val=0,
                func=self.energy_balance_func,
                dependents=self.energy_balance_dependents,
                num_eq_sets=1,
                quantity="efficiency",
                description="heat dissipation ratio relative to thermal input",
                calc=self._calc_eta, calc_deps=['ti']
            ),
            'Qloss': dc_cp(
                max_val=0, is_result=True, quantity="heat",
                description="heat dissipation",
                calc=self._calc_Qloss, calc_deps=['ti', 'eta']
            )
        })
        return params

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
        flue gas is forced to gaseous.

        - Reference temperature: 298.15 K.
        - Reference pressure: 1 bar.
        """
        T_ref = 298.15
        p_ref = 1e5

        res = 0
        for i in self.inl:
            res += i.m.val_SI * (
                i.h.val_SI
                - h_mix_pT(p_ref, T_ref, i.fluid_data, mixing_rule="forced-gas")
            )

        for o in self.outl:
            res -= o.m.val_SI * (
                o.h.val_SI
                - h_mix_pT(p_ref, T_ref, o.fluid_data, mixing_rule="forced-gas")
            )

        res += self._calc_ti() * self.eta.val_SI
        return res

    def _calc_eta(self):
        T_ref, p_ref = 298.15, 1e5
        res = sum(
            i.m.val_SI * (i.h.val_SI - h_mix_pT(p_ref, T_ref, i.fluid_data, mixing_rule="forced-gas"))
            for i in self.inl
        )
        res -= sum(
            o.m.val_SI * (o.h.val_SI - h_mix_pT(p_ref, T_ref, o.fluid_data, mixing_rule="forced-gas"))
            for o in self.outl
        )
        return -res / self.ti.val_SI

    def _calc_Qloss(self):
        return -(1 - self.eta.val_SI) * self.ti.val_SI

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        super().calc_parameters()
        for num, i in enumerate(self.inl):
            if i.p.val_SI < self.outl[0].p.val_SI:
                msg = (
                    f"The pressure at inlet {num + 1} is lower than the "
                    f"pressure at the outlet of component {self.label}."
                )
                logger.warning(msg)
