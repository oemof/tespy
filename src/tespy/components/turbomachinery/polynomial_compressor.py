
# -*- coding: utf-8

"""Module of class PolynomialCompressor.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/turbomachinery/polynomial_compressor.py

SPDX-License-Identifier: MIT
"""
from tespy.components.component import Component
from tespy.components.turbomachinery.base import Turbomachine
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import isentropic
from tespy.tools.helpers import _get_dependents


# the polynomial compressor is a fundamentally different component, therefore
# inherits from Turbomachine and not from Compressor
class PolynomialCompressor(Turbomachine):
    r"""
    Class for a compressor model following the EN12900 implementation of
    cite:`cecchinato2010`.

    See the example for the intended use of the component.

    **Mandatory Equations**

    - fluid: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - mass flow: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.dp_structure_matrix`
    - :py:meth:`tespy.components.component.Component.pr_structure_matrix`
    - :py:meth:`tespy.components.turbomachinery.polynomial_compressor.PolynomialCompressor.energy_balance_group_func`
    - :py:meth:`tespy.components.turbomachinery.polynomial_compressor.PolynomialCompressor.eta_s_group_func`
    - :py:meth:`tespy.components.turbomachinery.polynomial_compressor.PolynomialCompressor.eta_vol_group_func`

    Inlets/Outlets

    - in1
    - out1

    Optional inlets

    - power

    Image

    .. image:: /api/_images/Compressor.svg
       :alt: flowsheet of the compressor
       :align: center
       :class: only-light

    .. image:: /api/_images/Compressor_darkmode.svg
       :alt: flowsheet of the compressor
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
        Compressor power, :math:`P/\text{W}`

    Q_diss_rel : float, dict
        Relative heat loss of compressor, :math:`Q_\text{diss,rel}/1`

    eta_s : float, dict
        Isentropic efficiency, :math:`\eta_\text{s}/1`

    eta_vol : float, dict
        Volumetric efficiency, :math:`\eta_\text{vol}/1`

    reference_state: dict
        Reference state for the polynomial and displacement.

    pr : float, dict
        Outlet to inlet pressure ratio, :math:`pr/1`

    dp : float, dict
        Inlet to outlet pressure difference, :math:`dp/\text{p}_\text{unit}`
        Is specified in the Network's pressure unit

    Example
    -------
    The utilization of this component is intended to be done in two steps:

    1. Calculate the reference state isentropic and volumetric efficiency based
       on the provided polynomial formulation and the respective parameters.
    2. Set the resulting isentropic and volumetric efficiency. Under the
       assumption of constant efficiency, the outlet state will be determined
       with the volumetric flow at inlet. The volumetric flow at inlet scales
       with the rpm of the compressor.

    >>> from tespy.components import Source, Sink, PolynomialCompressor
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import numpy as np
    >>> from CoolProp.CoolProp import PropsSI
    >>> nw = Network(T_unit="C", p_unit="bar", iterinfo=False)
    >>> so = Source("from evaporator")
    >>> si = Sink("to condenser")
    >>> compressor = PolynomialCompressor("compressor")
    >>> c1 = Connection(so, "out1", compressor, "in1", label="c1")
    >>> c2 = Connection(compressor, "out1", si, "in1", label="c2")
    >>> nw.add_conns(c1, c2)

    Now, we provide a 10-coefficient polynome together with the reference
    state to the `calc_etas_from_polynome` method of the compressor. For this
    example we will take data for a compressor using "R134a" as working fluid.

    .. note::

        The 10-coefficient polynomes can be retrieved from manufacturers or
        need to be calculated based on manufacturer data. For example, Bitzer
        provides such data, which can be used to retrieve a polynomial. The
        polynomial coefficients for this example have been retrieved based on
        :cite:`bitzer2025_HSK`. This has been done under the assumption of
        5 % heat losses in the compressor power.

    >>> power_coefficients = np.array([
    ...     2.37745503e-01, 2.10760772e-01, 3.03362809e+00, 1.16180302e-02,
    ...     -1.42430044e-02, -5.56458029e-02, 3.49370564e-04, -4.33166771e-04,
    ...     2.38000000e-04, 6.49256691e-04
    ... ])
    >>> cooling_coefficients = np.array([
    ...     2.26460980e+00, 1.46392120e+01, 2.88891959e+01, 2.44045057e-01,
    ...     -3.52877642e-02, -8.11340980e-01, 1.53122397e-03, -1.41299188e-03,
    ...     -9.49473684e-04, 6.58981519e-03
    ... ])

    The reference state is available from the datasheet as well. It is defined
    through the reference evaporation and condensation temperatures, the
    respective superheat and subcooling. On top, the rpm for the polynomial and
    the displacement with its respective rpm.

    >>> reference_state = {
    ...     "T_evap": 0 + 273.15,
    ...     "T_cond": 50 + 273.15,
    ...     "T_sh": 10,
    ...     "T_sc": 0,
    ...     "rpm_poly": 50 * 60,
    ...     "displacement": 214,
    ...     "rpm_displacement": 20 * 60
    ... }

    We can now calculate the isentropic and volumetric efficiency.

    >>> polynomes = {
    ... "power": power_coefficients,
    ... "cooling": cooling_coefficients,
    ... }
    >>> eta = compressor.calc_etas_from_polynome(
    ...     "R134a", reference_state, polynomes
    ... )
    >>> round(eta["eta_s"], 2)
    np.float64(0.63)
    >>> round(eta["eta_vol"], 2)
    np.float64(0.87)

    We can take over these efficiency values and set the to the compressor
    together with the reference state and the assumption on heat dissipation.

    >>> compressor.set_attr(
    ...     eta_s=eta["eta_s"], eta_vol=eta["eta_vol"], Q_diss_rel=0.05,
    ...     reference_state=reference_state
    ... )

    First, we can impose the boundary conditions on "c1" that are equal to the
    displacement reference state. In that case, we should be able to get the
    same displacement value as inputted into the reference.

    >>> c1.set_attr(fluid={"R134a": 1}, T=0, Td_bp=10)
    >>> compressor.set_attr(rpm=1200)
    >>> p_sat = PropsSI("P", "Q", 0, "T", reference_state["T_cond"], "R134a")
    >>> c2.set_attr(p=p_sat / 1e5)
    >>> nw.solve("design")
    >>> round(c1.v.val * 3600 / eta["eta_vol"], 2)
    np.float64(214.0)

    We can also check the compressor power. It is higher than the power of an
    adiabatic compressor due to the heat dissipation. The compressor power plus
    heat dissipation will give the actual power required for isentropic
    compression. The heat dissipation is negative due to the heat leaving the
    component.

    >>> round(compressor.P.val)
    32829
    >>> round(compressor.Q_diss.val)
    -1641
    >>> round(compressor.P.val + compressor.Q_diss.val)
    31188

    It is also possible, to make the rpm a variable. This is useful, in case
    mass flow through the compressor is goverend from external. Usually, this
    could be the case, if a specific heat transfer is required to be provided
    by the condenser or from the evaporator. In this case, we just fix the
    displacement to mimic that.

    >>> compressor.set_attr(rpm="var")
    >>> c1.set_attr(v=400/3600)
    >>> nw.solve("design")
    >>> round(compressor.rpm.val)
    2569
    """

    @staticmethod
    def powerinlets():
        return ["power"]

    def _preprocess(self, num_nw_vars):
        return Component._preprocess(self, num_nw_vars)

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        if len(self.power_inl) > 0:
            constraints["energy_connector_balance"] = dc_cmc(**{
                "func": self.energy_connector_balance_func,
                "dependents": self.energy_connector_dependents,
                "num_eq_sets": 1
            })

        return constraints

    def get_parameters(self):
        return {
            "Q_diss": dc_cp(max_val=0),
            "P": dc_cp(min_val=0),
            "eta_vol": dc_cp(min_val=0, max_val=1),
            "Q_diss_rel": dc_cp(min_val=0, max_val=1),
            "rpm": dc_cp(min_val=0),
            "reference_state": dc_simple(),
            "eta_vol_group": dc_gcp(
                elements=["reference_state", "eta_vol", "rpm"],
                func=self.eta_vol_group_func,
                dependents=self.eta_vol_group_dependents,
                num_eq_sets=1
            ),
            "eta_s": dc_cp(min_val=0, max_val=1),
            "eta_s_group": dc_gcp(
                elements=["eta_s", "Q_diss_rel"],
                func=self.eta_s_group_func,
                dependents=self.eta_s_group_dependents,
                num_eq_sets=1
            ),
            # "energy_balance_group": dc_gcp(
            #     elements=["P", "Q_diss_rel"],
            #     func=self.energy_balance_group_func,
            #     dependents=self.energy_balance_group_dependents,
            #     num_eq_sets=1
            # )
        }

    # this is a bit different that in other cases, because the power cannot
    # directly be deduced from the change in enthalpy
    def energy_connector_balance_func(self):
        return (
            self.inl[0].m.val_SI
            * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
            / (1 - self.Q_diss_rel.val)
            - self.power_inl[0].E.val_SI
        )

    def energy_connector_dependents(self):
        return [
            self.power_inl[0].E, self.inl[0].m, self.inl[0].h, self.outl[0].h
        ]

    def eta_s_group_func(self):
        r"""Isentropic efficiency function featuning in a heat loss

        .. math::

            0 = \eta_\text{s} \cdot
            \frac{ h_\text{out} - h_\text{in}}{1 - \dot Q_\text{diss,rel}}
            - \left( h_\text{out,s} - h_\text{in} \right)

        .. note::

            For this, the actual enthalpy increase is further incresead by the
            relative heat loss to calculate the (virtual) compression outlet
            enthalpy. The compressor is first considered isentropic, then the
            relative heat loss (relative to compression power) is substracted
            from the outlet enthalpy afterwards.

        Returns
        -------
        float
            residual value

        """
        i = self.inl[0]
        o = self.outl[0]
        h_out_s = isentropic(
            i.p.val_SI,
            i.h.val_SI,
            o.p.val_SI,
            i.fluid_data,
            i.mixing_rule,
            T0=None
        )
        return (
            self.eta_s.val
            * (o.h.val_SI - i.h.val_SI) / (1 - self.Q_diss_rel.val)
            - (h_out_s - i.h.val_SI)
        )

    def eta_s_deriv(self, increment_filter, k, dependents=None):
        r"""
        Partial derivatives for isentropic efficiency.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        dependents = dependents["scalars"][0]
        i = self.inl[0]
        o = self.outl[0]
        f = self.eta_s_group_func

        if o.h.is_var and not i.h.is_var:
            self._partial_derivative(o.h, k, self.eta_s.val / (1 - self.Q_diss_rel.val), increment_filter)
            # remove o.h from the dependents
            dependents = dependents.difference(_get_dependents([o.h])[0])

        for dependent in dependents:
            self._partial_derivative(dependent, k, f, increment_filter)

    def eta_s_group_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
        ]

    def eta_vol_group_func(self):
        r"""Volumetric efficiency function.

        .. math::

            0 = \dot m_\text{in} -
            \eta_\text{vol} \cdot \frac{\dot V_\text{ref}}{3600} \cdot
            \frac{rpm}{rpm_\text{ref}} \cdot \frac{1}{v_\text{in}}

        Returns
        -------
        float
            residual value
        """
        displacement = (
            self.reference_state.val["displacement"] / 3600
            * self.rpm.val / self.reference_state.val["rpm_displacement"]
        )
        return (
            self.inl[0].m.val_SI
            - self.eta_vol.val * displacement / self.inl[0].calc_vol()
        )

    def eta_vol_group_dependents(self):
        return [self.inl[0].m, self.inl[0].p, self.inl[0].h, self.rpm]

    def calc_parameters(self):
        i = self.inl[0]
        o = self.outl[0]
        if self.Q_diss_rel.is_set:
            h_2 = (o.h.val_SI - i.h.val_SI * self.Q_diss_rel.val) / (1 - self.Q_diss_rel.val)
            self.P.val = i.m.val_SI * (h_2 - i.h.val_SI)
            self.Q_diss.val = i.m.val_SI * (o.h.val_SI - h_2)
            self.eta_s.val = (
                isentropic(
                    i.p.val_SI,
                    i.h.val_SI,
                    o.p.val_SI,
                    i.fluid_data,
                    i.mixing_rule,
                    T0=None
                ) - i.h.val_SI
            ) / (
                (h_2 - i.h.val_SI)
            )

    @staticmethod
    def calc_etas_from_polynome(fluid: str, reference_state: dict, polynomes: dict) -> dict:
        """Calculate the isentropic and polynomial efficiency for the reference
        state according to the procedure in :cite:`cecchinato2010`.

        Parameters
        ----------
        fluid : str
            Name of the fluid
        reference_state : dict
            Dictionary with reference state information, i.e.

            - T_evap: evaporation temperature (°C)
            - T_cond: condensation temperature (°C)
            - T_sh: superheating delta T (Kelvin)
            - T_sc: subcooling delta T (Kelvin)
            - rpm_poly: reference rpm of the polynomial
            - displancement: displacement in m3/h
            - rpm_displacement: reference rpm of the displacement

        polynomes : dict
            Dictionary containing the power and the evaporator heat polynomial

        Returns
        -------
        dict
            Dictonary with keys "eta_s" and "eta_vol" for isentropic and
            volumetric efficiency
        """
        # this method will have few calls, therefore high-level interface
        # usage is fine
        from CoolProp.CoolProp import PropsSI as PSI

        T_evap = reference_state["T_evap"]
        T_cond = reference_state["T_cond"]
        P_comp = calc_EN12900(polynomes["power"], T_evap, T_cond) * 1000
        Q_evap = calc_EN12900(polynomes["cooling"], T_evap, T_cond) * 1000

        p_evap = PSI("P", "T", T_evap, "Q", 1, fluid)
        p_cond = PSI("P", "T", T_cond, "Q", 0, fluid)
        T_sh = reference_state["T_sh"]
        T_sc = reference_state["T_sc"]

        if T_sh > 0:
            h_evap_out = PSI("H", "P", p_evap, "T", T_evap + T_sh, fluid)
        else:
            h_evap_out = PSI("H", "P", p_evap, "Q", 1, fluid)

        if T_sc > 0:
            h_cond_out = PSI("H", "P", p_cond, "T", T_cond - T_sc, fluid)
        else:
            h_cond_out = PSI("H", "P", p_cond, "Q", 0, fluid)

        s_comp_in = PSI("S", "P", p_evap, "H", h_evap_out, fluid)
        h_comp_s = PSI("H", "P", p_cond, "S", s_comp_in, fluid)

        dot_m = Q_evap / (h_evap_out - h_cond_out)
        eta_s = dot_m * (h_comp_s - h_evap_out) / P_comp

        rpm_poly = reference_state["rpm_poly"]
        rpm_displacement = reference_state["rpm_displacement"]
        displacement = (
            (reference_state["displacement"] / 3600)
            * (rpm_poly / rpm_displacement)
        )
        rho_comp_in = PSI("D", "P", p_evap, "H", h_evap_out, fluid)
        eta_vol = dot_m / (rho_comp_in * displacement)

        return {
            "eta_s": eta_s,
            "eta_vol": eta_vol
        }


def calc_EN12900(c: list, t_evap: float, t_cond: float) -> float:
    """Calculate the EN12900 polynomial

    Parameters
    ----------
    c : list
        Coefficients of the polynomial
    t_evap : float
        Reference evaporation temperature
    t_cond : float
        Reference condensation temperature

    Returns
    -------
    float
        Result of polynomial
    """
    t_evap = t_evap - 273.15
    t_cond = t_cond - 273.15
    return (
        c[0]
        + c[1] * t_evap + c[2] * t_cond
        + c[3] * t_evap**2 + c[4] * t_evap * t_cond + c[5] * t_cond ** 2
        + c[6] * t_evap**3 + c[7] * t_evap ** 2 * t_cond
        + c[8] * t_evap * t_cond ** 2 + c[9] * t_cond ** 3
    )
