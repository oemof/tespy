# -*- coding: utf-8

"""Module of class Turbine.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/turbomachinery/turbine.py

SPDX-License-Identifier: MIT
"""

from scipy.optimize import brentq

from tespy.components.component import component_registry
from tespy.components.turbomachinery.turbine import Turbine
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import isentropic
from tespy.tools.fluid_properties.helpers import single_fluid
from tespy.tools.helpers import fluidalias_in_list
from tespy.tools.logger import logger


@component_registry
class SteamTurbine(Turbine):
    r"""
    Class for steam turbines with wet expansion.

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_func`
    - :py:meth:`tespy.components.turbomachinery.base.Turbomachine.energy_balance_func`
    - :py:meth:`tespy.components.turbomachinery.turbine.Turbine.eta_dry_s_func`
    - :py:meth:`tespy.components.turbomachinery.turbine.Turbine.eta_s_func`
    - :py:meth:`tespy.components.turbomachinery.turbine.Turbine.eta_s_char_func`
    - :py:meth:`tespy.components.turbomachinery.turbine.Turbine.cone_func`

    Inlets/Outlets

    - in1
    - out1

    Image

    .. image:: /api/_images/Turbine.svg
       :alt: flowsheet of the turbine
       :align: center
       :class: only-light

    .. image:: /api/_images/Turbine_darkmode.svg
       :alt: flowsheet of the turbine
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

    P : float, dict
        Power, :math:`P/\text{W}`

    eta_s : float, dict
        Isentropic efficiency, :math:`\eta_s/1`

    eta_dry_s : float, dict
        Dry isentropic efficiency, :math:`\eta_s/1`

    pr : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio, :math:`pr/1`

    eta_s_char : tespy.tools.characteristics.CharLine, dict
        Characteristic curve for isentropic efficiency, provide CharLine as
        function :code:`func`.

    cone : dict
        Apply Stodola's cone law (works in offdesign only).

    Example
    -------
    A steam turbine expands 10 kg/s of superheated steam at 550 °C and 110 bar
    to 0,5 bar at the outlet. For example, it is possible to calulate the power
    output and vapour content at the outlet for a given isentropic efficiency.
    The :code:`SteamTurbine` class follows the implementation of the Baumann
    rule :cite:`Tanuma2017`

    >>> from tespy.components import Sink, Source, SteamTurbine
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> si = Sink('sink')
    >>> so = Source('source')
    >>> st = SteamTurbine('steam turbine')
    >>> st.component()
    'steam turbine'
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

    @staticmethod
    def component():
        return 'steam turbine'

    def get_parameters(self):

        params = super().get_parameters()
        params["alpha"] = dc_cp(min_val=0.4, max_val=2.5)
        params["eta_s_dry"] = dc_cp(min_val=0.0, max_val=1.0)
        params["eta_s_dry_group"] = dc_gcp(
            num_eq=1, elements=["alpha", "eta_s_dry"],
            func=self.eta_s_wet_func,
            deriv=self.eta_s_wet_deriv
        )

        return params

    def preprocess(self, num_nw_vars):

        fluid = single_fluid(self.inl[0].fluid_data)
        if fluid is None:
            msg = "The SteamTurbine can only be used with pure fluids."
            logger.error(msg)
            raise ValueError(msg)

        if not fluidalias_in_list(fluid, ["water"]):
            msg = "The SteamTurbine is intended to be used with water only."
            logger.warning(msg)

        return super().preprocess(num_nw_vars)

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
                - self.eta_s_dry.val * (1 - self.alpha.val * ym)
            )

        else:  # superheated vapour
            if outl.calc_phase() == "g":
                return self.calc_eta_s() - self.eta_s_dry.val

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
                    inl.h.val_SI - self.eta_s_dry.val
                    * (inl.h.val_SI - hout_isen)
                )

                # calculate enthalpy of saturated vapour at psat
                hsat = h_mix_pQ(psat, 1, inl.fluid_data)

                return hout - hsat

            frac = brentq(find_sat, 1, 0)
            psat = inl.p.val_SI - frac * dp
            hsat = h_mix_pQ(psat, 1, inl.fluid_data)

            # calculate the isentropic efficiency for wet expansion
            ym = 1 - (1.0 + outl.calc_x()) / 2  # average wetness
            eta_s = self.eta_s_dry.val * (1 - self.alpha.val * ym)

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

    def eta_s_wet_deriv(self, increment_filter, k):
        r"""
        Partial derivatives for dry isentropic efficiency function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        f = self.eta_s_wet_func
        i = self.inl[0]
        o = self.outl[0]
        if self.is_variable(i.p, increment_filter):
            self.jacobian[k, i.p.J_col] = self.numeric_deriv(f, "p", i)
        if self.is_variable(o.p, increment_filter):
            self.jacobian[k, o.p.J_col] = self.numeric_deriv(f, "p", o)
        if self.is_variable(i.h, increment_filter):
            self.jacobian[k, i.h.J_col] = self.numeric_deriv(f, "h", i)
        if self.is_variable(o.h, increment_filter):
            self.jacobian[k, o.h.J_col] = self.numeric_deriv(f, "h", o)
