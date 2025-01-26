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
from tespy.tools.fluid_properties import isentropic
from tespy.tools.fluid_properties import h_mix_pQ


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

    >>> from tespy.components import Sink, Source, Turbine
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools import ComponentCharacteristics as dc_cc
    >>> import shutil
    >>> nw = Network(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> si = Sink('sink')
    >>> so = Source('source')
    >>> st = SteamTurbine('steam turbine')
    >>> st.component()
    'turbine'
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

    To capture the effect of liquid drop-out on the isentropic efficiency, the dry turbine efficiency is specified
    >>> st.set_attr(eta_s=None)
    >>> st.set_attr(eta_dry_s=0.9)
    >>> nw.solve('design')
    >>> round(st.P.val, 0)
    -7009682.0
    >>> round(outg.x.val, 3)
    0.840
    """


    @staticmethod
    def component():
        return 'steam turbine'

    def get_parameters(self):

        params = super().get_parameters()

        params["alpha"] = dc_cp(
            min_val=0.4, max_val=2.5, num_eq=0)
        params["eta_dry_s"] = dc_cp(
            min_val=0, max_val=1, num_eq=1,
            func=self.eta_dry_s_func,
            deriv=self.eta_dry_s_deriv
        )

        return params

    def eta_dry_s_func(self):

        inl = self.inl[0]
        outl = self.outl[0]

        state = inl.calc_phase()
        if state == "tp":  # two-phase or saturated vapour

            ym = 1 - (inl.calc_x() + outl.calc_x()) / 2  # average wetness
            self.eta_s.val = self.eta_dry.val * (1 - self.alpha.val * ym)

            return self.eta_s_func()

        else:  # superheated vapour
            dp = inl.p.val_SI - outl.p.val_SI

            # compute the pressure and enthalpy at which the expansion enters the vapour dome
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
                hout = inl.h.val_SI - self.eta_dry_s.val * (inl.h.val_SI - hout_isen)

                # calculate enthalpy of saturated vapour at psat
                hsat = h_mix_pQ(psat, 1, inl.fluid_data)

                return hout - hsat
            frac = brentq(find_sat, 1, 0)
            psat = inl.p.val_SI - frac * dp
            hsat = h_mix_pQ(psat, 1, inl.fluid_data)

            # calculate the isentropic efficiency for wet expansion
            ym = 1 - (1.0 + outl.calc_x()) / 2  # average wetness
            eta_s = self.eta_dry_s.val * (1 - self.alpha.val * ym)

            # calculate the final outlet enthalpy
            hout_isen = isentropic(
                                psat,
                                hsat,
                                outl.p.val_SI,
                                inl.fluid_data,
                                inl.mixing_rule,
                                T0=inl.T.val_SI
                            )
            hout = hsat - eta_s * (hsat- hout_isen)

            # calculate the difference in enthalpy
            dh_isen = hout - inl.h.val_SI
            dh = outl.h.val_SI - inl.h.val_SI

            # return residual
            return -dh + dh_isen

    def eta_dry_s_deriv(self, increment_filter, k):

        f = self.eta_dry_s_func
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
