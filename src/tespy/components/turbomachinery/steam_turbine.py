# -*- coding: utf-8

"""Module of class SteamTurbine.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/turbomachinery/steam_turbine.py

SPDX-License-Identifier: MIT
"""

from scipy.optimize import brentq

from tespy.components.component import component_registry
from tespy.components.turbomachinery.turbine import Turbine
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.fluid_properties import Q_mix_ph
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import isentropic
from tespy.tools.fluid_properties.helpers import single_fluid
from tespy.tools.helpers import fluidalias_in_list
from tespy.tools.logger import logger


@component_registry
class SteamTurbine(Turbine):
    r"""
    Class for steam turbines with wet expansion.

    .. image:: /api/_images/components/Turbine.svg
       :alt: flowsheet of the steamturbine
       :align: center
       :class: only-light

    .. image:: /api/_images/components/Turbine_darkmode.svg
       :alt: flowsheet of the steamturbine
       :align: center
       :class: only-dark

    Ports
    -----

    Fluid inlets: in1

    Fluid outlets: out1

    Power outlets: power

    Mandatory Equations
    -------------------

    - mass flow equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
    - fluid composition equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`

    When a power or heat connector is attached:

    - energy_connector_balance: :py:meth:`energy_connector_balance_func <tespy.components.turbomachinery.turbine.Turbine.energy_connector_balance_func>`

    Parameters
    ----------

    alpha : float, dict
        Influence factor for wetness efficiency modifier. Quantity:
        :code:`ratio`.

    char_warnings : bool
        Ignore warnings on default characteristics usage for this component.

    cone : bool
        Cone law equation for offdesign.
        Equation: :py:meth:`cone_func <tespy.components.turbomachinery.turbine.Turbine.cone_func>`.

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the components design case.

    dp : float, dict
        Inlet to outlet absolute pressure change. Quantity:
        :code:`pressure_difference`.
        Equation: :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`.

    eta_s : float, dict
        Isentropic efficiency. Quantity: :code:`efficiency`.
        Equation: :py:meth:`eta_s_func <tespy.components.turbomachinery.turbine.Turbine.eta_s_func>`.

    eta_s_char : tespy.tools.characteristics.CharLine, dict
        Isentropic efficiency lookup table for offdesign.
        Equation: :py:meth:`eta_s_char_func <tespy.components.turbomachinery.turbine.Turbine.eta_s_char_func>`.

    eta_s_dry : float, dict
        Isentropic efficiency of dry expansion. Quantity: :code:`efficiency`.

    eta_s_dry_group : GroupedComponentProperties
        Method to apply Baumann rule. Elements: :code:`alpha`,
        :code:`eta_s_dry`.
        Equation: :py:meth:`eta_s_wet_func <tespy.components.turbomachinery.steam_turbine.SteamTurbine.eta_s_wet_func>`.

    label : str
        The label of the component.

    local_design : bool
        Treat this component in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this component in offdesign mode in a design calculation.

    offdesign : list
        List containing offdesign parameters (stated as String).

    P : float, dict
        Power input/output of the component. Quantity: :code:`power`.
        Equation: :py:meth:`energy_balance_func <tespy.components.turbomachinery.base.Turbomachine.energy_balance_func>`.

    pr : float, dict
        Outlet to inlet pressure ratio. Quantity: :code:`ratio`.
        Equation: :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`.

    printout : bool
        Include this component in the network's results printout.

    Example
    -------
    A steam turbine expands 10 kg/s of superheated steam at 550 °C and 110 bar
    to 0,5 bar at the outlet. For example, it is possible to calculate the
    power output and vapour content at the outlet for a given isentropic
    efficiency. The :code:`SteamTurbine` class follows the implementation of
    the Baumann rule :cite:`Tanuma2017`

    >>> from tespy.components import Sink, Source, SteamTurbine
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC", "enthalpy": "kJ/kg"
    ... })
    >>> si = Sink('sink')
    >>> so = Source('source')
    >>> st = SteamTurbine('steam turbine')
    >>> inc = Connection(so, 'out1', st, 'in1')
    >>> outg = Connection(st, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)

    In design conditions the isentropic efficiency is specified.
    >>> st.set_attr(eta_s=0.9)
    >>> inc.set_attr(fluid={'water': 1}, m=10, T=250, p=20)
    >>> outg.set_attr(p=0.1)
    >>> nw.solve('design')
    >>> round(st.P.val, 0)
    -7471296.0
    >>> round(outg.x.val, 3)
    0.821

    To capture the effect of liquid drop-out on the isentropic
    efficiency, the dry turbine efficiency is specified
    >>> st.set_attr(eta_s=None)
    >>> st.set_attr(eta_s_dry=0.9, alpha=1.0)
    >>> nw.solve('design')
    >>> round(st.P.val, 0)
    -7009682.0
    >>> round(outg.x.val, 3)
    0.84
    """

    def get_parameters(self):

        params = super().get_parameters()
        params["alpha"] = dc_cp(
            min_val=0.4, max_val=2.5, quantity="ratio",
            description="influence factor for wetness efficiency modifier"
        )
        params["eta_s_dry"] = dc_cp(
            min_val=0.0, max_val=1.0, quantity="efficiency",
            description="isentropic efficiency of dry expansion"
        )
        params["eta_s_dry_group"] = dc_gcp(
            num_eq_sets=1, elements=["alpha", "eta_s_dry"],
            func=self.eta_s_wet_func,
            dependents=self.eta_s_dependents,  # same dependents!
            description="method to apply Baumann rule"
        )

        return params

    def _preprocess(self, num_nw_vars):

        fluid = single_fluid(self.inl[0].fluid_data)
        if fluid is None:
            msg = "The SteamTurbine can only be used with pure fluids."
            logger.error(msg)
            raise ValueError(msg)

        if not fluidalias_in_list(fluid, ["water"]):
            msg = "The SteamTurbine is intended to be used with water only."
            logger.warning(msg)

        return super()._preprocess(num_nw_vars)

    def eta_s_wet_func(self):
        r"""
        Equation for given dry isentropic efficiency of a turbine under wet
        expansion.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = -\left( h_{out} - h_{in} \right) +
                \left( h_{out,s} - h_{in} \right) \cdot \eta_{s,e}

                \eta_{s,e} = \eta_{s,e}^{dry} \cdot \left( 1 - \alpha
                \cdot y_m \right)

                y_m = \frac{\left( 1-x_{in}\right)+ \left( 1-x_{out}
                \right)}{2}
        """
        inl = self.inl[0]
        outl = self.outl[0]

        state = inl.calc_phase()
        if state == "tp":  # two-phase or saturated vapour

            ym = 1 - (inl.calc_x() + outl.calc_x()) / 2  # average wetness
            return (
                self.calc_eta_s()
                - self.eta_s_dry.val_SI * (1 - self.alpha.val_SI * ym)
            )

        else:  # superheated vapour
            if outl.calc_phase() == "g":
                return self.calc_eta_s() - self.eta_s_dry.val_SI

            dp = inl.p.val_SI - outl.p.val_SI

            # compute the pressure and enthalpy at which the expansion
            # enters the vapour dome
            def find_sat(frac):
                psat = inl.p.val_SI - frac * dp

                # calculate enthalpy under dry expansion to psat
                hout_isen = isentropic(
                    inl.p.val_SI,
                    inl.h.val_SI,
                    psat,
                    inl.fluid_data,
                    inl.mixing_rule,
                    T0=inl.T.val_SI
                )
                hout = (
                    inl.h.val_SI - self.eta_s_dry.val_SI
                    * (inl.h.val_SI - hout_isen)
                )

                # calculate enthalpy of saturated vapour at psat
                hsat = h_mix_pQ(psat, 1, inl.fluid_data)

                return hout - hsat

            frac = 1
            if round(outl.calc_Q(), 3) != 1:
                frac = brentq(find_sat, 1, 0)

            psat = inl.p.val_SI - frac * dp
            hsat = h_mix_pQ(psat, 1, inl.fluid_data)

            # calculate the isentropic efficiency for wet expansion
            ym = 1 - (1.0 + outl.calc_x()) / 2  # average wetness
            eta_s = self.eta_s_dry.val_SI * (1 - self.alpha.val_SI * ym)

            # calculate the final outlet enthalpy
            hout_isen = isentropic(
                psat,
                hsat,
                outl.p.val_SI,
                inl.fluid_data,
                inl.mixing_rule,
                T0=inl.T.val_SI
            )
            hout = hsat - eta_s * (hsat - hout_isen)

            # return residual: outlet enthalpy = calculated outlet enthalpy
            return outl.h.val_SI - hout

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
        if key == 'p':
            return super().initialise_source(c, key)
        elif key == 'h':
            fluid = single_fluid(c.fluid_data)
            if fluid is not None:
                return h_mix_pQ(c.p.val_SI, 0.9, c.fluid_data, c.mixing_rule)
