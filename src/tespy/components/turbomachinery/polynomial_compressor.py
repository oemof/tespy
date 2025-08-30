
# -*- coding: utf-8

"""Module of class PolynomialCompressor.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/turbomachinery/polynomial_compressor.py

SPDX-License-Identifier: MIT
"""
import numpy as np

from tespy.components.component import Component
from tespy.components.turbomachinery.base import Turbomachine
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import T_sat_p
from tespy.tools.fluid_properties import isentropic


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

    eta_s_poly : array, dict
        Polynomial coefficients for isentropic efficiency

    eta_vol : float, dict
        Volumetric efficiency, :math:`\eta_\text{vol}/1`

    eta_vol_poly : array, dict
        Polynomial coefficients for volumetric efficiency

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

    1. Calculate the reference state isentropic and volumetric efficiency
       polynomials based on the provided manufacturer data.
    2. Set the resulting isentropic and volumetric efficiency polynomials.
       Under the assumption of isentropic and volumetric efficiency not being
       constant at variable compressor rpm, the outlet state will be determined
       with the volumetric flow at inlet. The volumetric flow at inlet scales
       linearly with the rpm of the compressor.

    >>> from tespy.components import Source, Sink, PolynomialCompressor
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import pandas as pd
    >>> from CoolProp.CoolProp import PropsSI
    >>> nw = Network(T_unit="C", p_unit="bar", iterinfo=False)
    >>> so = Source("from evaporator")
    >>> si = Sink("to condenser")
    >>> compressor = PolynomialCompressor("compressor")
    >>> c1 = Connection(so, "out1", compressor, "in1", label="c1")
    >>> c2 = Connection(compressor, "out1", si, "in1", label="c2")
    >>> nw.add_conns(c1, c2)

    Now, we can either provide

    - a 10-coefficient polynomial for power and cooling or
    - provide the respective power and cooling energy from a datasheet of a
      compressor manufacturer to generate such a polynomial

    Then we can used a precalculation method, which transforms the polynomial
    or the data into two polynomials, one for the isentropic efficiency and one
    for the volumetric efficiency, both as a function of evaporation and
    condensation temperature. Additionally information on a reference state
    have to be provided, which include

    - superheating at suction
    - subcooling after the condensation
    - the rpm belonging to the original data
    - a displacement value (kg/h) with the respective rpm for this displacement

    .. tip::

        The compressor data or the 10-coefficient polynomials can be retrieved
        from manufacturers. For example, Bitzer provides such data, which can
        be used to retrieve a polynomial. The data for this example have been
        retrieved from :cite:`bitzer2025_HSK`.

    >>> reference_state = {
    ...     "T_sh": 20,  # superheating
    ...     "T_sc": 0,  # subcooling
    ...     "rpm_poly": 50 * 60,  # rpm belonging to the polynomial data
    ...     "rpm_displacement": 20 * 60,  # rpm belonging to the displacement
    ...     "displacement": 214  # kg / h
    ... }
    >>> power = pd.DataFrame(
    ...     columns=[10,7.5,5,0,-5,-10], index=[30, 40, 50], dtype=float
    ... )
    >>> cooling = power.copy()
    >>> cooling.loc[30] = [465600,424100,385500,316700,257900,208000]
    >>> cooling.loc[40] = [418900,380400,344800,281400,227400,181600]
    >>> cooling.loc[50] = [365900,331300,299200,242100,193700,152900]
    >>> power.loc[30] = [62.0,61.8,61.8,61.8,61.7,61.3]
    >>> power.loc[40] = [78.0,78.0,78.0,78.0,77.7,76.8]
    >>> power.loc[50] = [99.2,99.2,99.2,98.9,98.1,96.5]
    >>> power = power * 1000

    .. attention::

        The data or polynomial formulations must be in SI units!

    We can now use the inbuilt method to determine the isentropic and
    volumetric efficiency polynomials. For that we need to import the
    respective method. Apart from this method, there is also the
    :py:func:`tespy.components.turbomachinery.polynomial_compressor.generate_eta_polys_from_power_and_cooling_polys`
    method, that can do the same step provided a polynomial for power and one
    for the cooling.

    >>> from tespy.components.turbomachinery.polynomial_compressor import (
    ...     generate_eta_polys_from_data
    ... )
    >>> eta_s_poly, eta_vol_poly = generate_eta_polys_from_data(
    ...     power, cooling, "R134a", reference_state
    ... )
    >>> eta_s_poly
    array([ 3.44223012e-03, -3.75139140e-02,  4.39204462e-02, -9.21644870e-04,
            1.68576190e-03, -8.97540501e-04, -7.54781107e-06,  1.61377008e-05,
           -1.53820046e-05,  5.04818089e-06])
    >>> eta_vol_poly
    array([ 5.81192914e-03, -7.18820053e-04,  7.41463587e-02,  2.84410052e-05,
            6.51372426e-05, -1.89872495e-03,  7.84206012e-07, -1.90585865e-06,
            4.52695494e-07,  1.51321175e-05])

    We can take these polynomials and set them on the compressor instance
    together with the reference state and the assumption on heat dissipation.

    >>> compressor.set_attr(
    ...     eta_s_poly=eta_s_poly, eta_vol_poly=eta_vol_poly, Q_diss_rel=0.05,
    ...     reference_state=reference_state
    ... )

    First, we can impose the boundary conditions on "c1" that are equal to the
    displacement reference state. In that case, we should be able to get the
    same displacement value as inputted into the reference.

    >>> c1.set_attr(fluid={"R134a": 1}, T=0, Td_bp=10)  # T_evap=-10°C
    >>> compressor.set_attr(rpm=1200)
    >>> p_sat = PropsSI("P", "Q", 0, "T", 50 + 273.15, "R134a")  # T_cond=50°C
    >>> c2.set_attr(p=p_sat / 1e5)
    >>> nw.solve("design")
    >>> round(c1.v.val * 3600 / compressor.eta_vol.val, 2)
    214.0
    >>> round(compressor.eta_s.val, 3)
    0.5
    >>> round(compressor.eta_vol.val, 3)
    0.814

    We can also double check our resulting isentropic and volumetric efficiency
    values with the evaluation of the polynomials.

    >>> from tespy.components.turbomachinery.polynomial_compressor import (
    ...     calc_EN12900
    ... )
    >>> round(compressor.eta_s.val, 3) == round(calc_EN12900(eta_s_poly, -10, 50), 3)
    np.True_
    >>> round(compressor.eta_vol.val, 3) == round(calc_EN12900(eta_vol_poly, -10, 50), 3)
    np.True_

    .. tip::

        You can also create polynomials for power and cooling from respective
        data. For that, import the
        :py:func:`tespy.components.turbomachinery.polynomial_compressor.fit_EN12900`
        method and pass the respective data.

    We can also check the compressor power. It is higher than the power of an
    adiabatic compressor due to the heat dissipation. The compressor power plus
    heat dissipation will give the actual power required for isentropic
    compression. The heat dissipation is negative due to the heat leaving the
    component.

    >>> round(compressor.P.val)
    38385
    >>> round(compressor.Q_diss.val)
    -1919
    >>> round(compressor.P.val + compressor.Q_diss.val)
    36466

    Now, let's see what happens, if evaporation or condensation temperature
    change:

    >>> c1.set_attr(T=20, Td_bp=10)  # T_evap=10°C
    >>> p_sat = PropsSI("P", "Q", 0, "T", 40 + 273.15, "R134a")  # T_cond=40°C
    >>> c2.set_attr(p=p_sat / 1e5)
    >>> nw.solve("design")
    >>> round(compressor.eta_s.val, 3)
    0.665
    >>> round(compressor.eta_vol.val, 3)
    0.924

    It is also possible, to make the rpm a variable. This is useful, in case
    mass flow through the compressor is governed from external. Usually, this
    could be the case, if a specific heat transfer is required to be provided
    by the condenser or from the evaporator. In this case, we just fix the
    displacement to mimic that.

    >>> compressor.set_attr(rpm="var")
    >>> c1.set_attr(v=400/3600)
    >>> nw.solve("design")
    >>> round(compressor.rpm.val)
    2427

    As final remarks: You can also set fixed isentropic and fixed volumetric
    efficiencies for these components.
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
        params = super().get_parameters()
        params.update({
            "Q_diss": dc_cp(max_val=0),
            "P": dc_cp(min_val=0),
            "eta_vol": dc_cp(min_val=0, max_val=1),
            "Q_diss_rel": dc_cp(min_val=0, max_val=1),
            "rpm": dc_cp(min_val=0),
            "reference_state": dc_simple(),
            "eta_s_poly": dc_simple(),
            "eta_vol_poly": dc_simple(),
            "eta_vol_poly_group": dc_gcp(
                elements=["reference_state", "eta_vol_poly", "rpm"],
                func=self.eta_vol_poly_group_func,
                dependents=self.eta_vol_poly_group_dependents,
                num_eq_sets=1
            ),
            "eta_vol_group": dc_gcp(
                elements=["reference_state", "eta_vol", "rpm"],
                func=self.eta_vol_group_func,
                dependents=self.eta_vol_group_dependents,
                num_eq_sets=1
            ),
            "eta_s": dc_cp(min_val=0, max_val=1),
            "eta_s_poly_group": dc_gcp(
                elements=["eta_s_poly", "Q_diss_rel"],
                func=self.eta_s_poly_group_func,
                dependents=self.eta_s_group_dependents,
                num_eq_sets=1
            ),
            "eta_s_group": dc_gcp(
                elements=["eta_s", "Q_diss_rel"],
                func=self.eta_s_group_func,
                dependents=self.eta_s_group_dependents,
                num_eq_sets=1
            ),
            "energy_balance_group": dc_gcp(
                elements=["P", "Q_diss_rel"],
                func=self.energy_balance_group_func,
                dependents=self.energy_balance_group_dependents,
                num_eq_sets=1
            )
        })
        return params

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

    def energy_balance_group_func(self):
        r"""Equation for specified power and relative heat dissipation

        .. math::

            \dot m \cdot
            \frac{ h_\text{out} - h_\text{in}}{1 - \dot Q_\text{diss,rel}}
            - P

        Returns
        -------
        float
            residual
        """
        return (
            self.inl[0].m.val_SI
            * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
            / (1 - self.Q_diss_rel.val)
            - self.P.val
        )

    def energy_balance_group_dependents(self):
        return [self.inl[0].m, self.inl[0].h, self.outl[0].h]

    def eta_s_group_func(self):
        r"""Isentropic efficiency function with a heat loss term

        .. math::

            0 = \eta_\text{s} \left(T_\text{evap}, T_\text{cond}\right)\cdot
            \frac{ h_\text{out} - h_\text{in}}{1 - \dot Q_\text{diss,rel}}
            - \left( h_\text{out,s} - h_\text{in} \right)

        .. note::

            For this, the actual enthalpy increase is further increased by the
            relative heat loss to calculate the (virtual) compression outlet
            enthalpy. The compressor is first considered isentropic, then the
            relative heat loss (relative to compression power) is subtracted
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

    def eta_s_poly_group_func(self):
        r"""Isentropic efficiency function with a heat loss term

        .. math::

            0 = \eta_\text{s} \left(T_\text{evap}, T_\text{cond}\right)\cdot
            \frac{ h_\text{out} - h_\text{in}}{1 - \dot Q_\text{diss,rel}}
            - \left( h_\text{out,s} - h_\text{in} \right)

        .. note::

            For this, the actual enthalpy increase is further increased by the
            relative heat loss to calculate the (virtual) compression outlet
            enthalpy. The compressor is first considered isentropic, then the
            relative heat loss (relative to compression power) is subtracted
            from the outlet enthalpy afterwards.

        Returns
        -------
        float
            residual value

        """
        i = self.inl[0]
        o = self.outl[0]

        t_evap = T_sat_p(i.p.val_SI, i.fluid_data)
        t_cond = T_sat_p(o.p.val_SI, o.fluid_data)
        eta_s = _calc_EN12900_SI(self.eta_s_poly.val, t_evap, t_cond)
        h_out_s = isentropic(
            i.p.val_SI,
            i.h.val_SI,
            o.p.val_SI,
            i.fluid_data,
            i.mixing_rule,
            T0=None
        )

        return (
            eta_s * (o.h.val_SI - i.h.val_SI) / (1 - self.Q_diss_rel.val)
            - (h_out_s - i.h.val_SI)
        )

    def eta_s_group_dependents(self):
        return [
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
        i = self.inl[0]

        displacement = (
            self.reference_state.val["displacement"] / 3600
            * self.rpm.val / self.reference_state.val["rpm_displacement"]
        )

        return (
            i.m.val_SI - self.eta_val.val * displacement / i.calc_vol()
        )

    def eta_vol_group_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.rpm
        ]

    def eta_vol_poly_group_func(self):
        r"""Volumetric efficiency function.

        .. math::

            0 = \dot m_\text{in} -
            \eta_\text{vol} \left(T_\text{evap}, T_\text{cond}\right) \cdot
            \frac{\dot V_\text{ref}}{3600} \cdot
            \frac{rpm}{rpm_\text{ref}} \cdot \frac{1}{v_\text{in}}

        Returns
        -------
        float
            residual value
        """
        i = self.inl[0]
        o = self.outl[0]

        t_evap = T_sat_p(i.p.val_SI, i.fluid_data)
        t_cond = T_sat_p(o.p.val_SI, o.fluid_data)
        eta_vol = _calc_EN12900_SI(self.eta_vol_poly.val, t_evap, t_cond)
        displacement = (
            self.reference_state.val["displacement"] / 3600
            * self.rpm.val / self.reference_state.val["rpm_displacement"]
        )

        return (
            i.m.val_SI - eta_vol * displacement / i.calc_vol()
        )

    def eta_vol_poly_group_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.rpm
        ]

    def calc_parameters(self):
        i = self.inl[0]
        o = self.outl[0]

        self.pr.val = o.p.val_SI / i.p.val_SI
        self.dp.val_SI = i.p.val_SI - o.p.val_SI
        self.dp.val = i.p.val - o.p.val

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

        if self.reference_state.is_set:
            displacement = (
                self.reference_state.val["displacement"] / 3600
                * self.rpm.val / self.reference_state.val["rpm_displacement"]
            )
            self.eta_vol.val = i.m.val_SI * i.calc_vol() / displacement


def _calc_etas_from_polynome(fluid: str, T_evap: float, T_cond: float, reference_state: dict, polynomes: dict) -> dict:
    """Calculate the isentropic and polynomial efficiency for the reference
    state according to the procedure in :cite:`cecchinato2010`.

    Parameters
    ----------
    fluid : str
        Name of the fluid
    T_evap : float
        Evaporation temperature in K
    T_cond : float
        Evaporation temperature in K
    reference_state : dict
        Dictionary with reference state information, i.e.
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
        Dictionary with keys "eta_s" and "eta_vol" for isentropic and
        volumetric efficiency
    """
    # this method will have few calls, therefore high-level interface
    # usage is fine
    from CoolProp.CoolProp import PropsSI as PSI

    P_comp = _calc_EN12900_SI(polynomes["power"], T_evap, T_cond)
    Q_evap = _calc_EN12900_SI(polynomes["cooling"], T_evap, T_cond)

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
    return (
        c[0]
        + c[1] * t_evap + c[2] * t_cond
        + c[3] * t_evap**2 + c[4] * t_evap * t_cond + c[5] * t_cond ** 2
        + c[6] * t_evap**3 + c[7] * t_evap ** 2 * t_cond
        + c[8] * t_evap * t_cond ** 2 + c[9] * t_cond ** 3
    )


def _calc_EN12900_SI(c: list, t_evap: float, t_cond: float) -> float:
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
    return calc_EN12900(c, t_evap - 273.15, t_cond - 273.15)


def fit_EN12900(t_evap: np.array, t_cond: np.array, data: np.array, check_diff: bool=True) -> np.array:
    """Fit the polynome coefficients of EN12900 polynome based on evaporation
    and condensation temperature and respective measurements

    Parameters
    ----------
    t_evap : np.array
        1-d array of evaporation temperatures
    t_cond : np.array
        1-d array of condensation temperatures
    data : np.array
        datasheet information

    Returns
    -------
    np.array
        1-d array of polynome coefficients
    """
    z = data.flatten()
    t_cond, t_evap = np.meshgrid(t_evap, t_cond)  # needs to be inverted here
    x = t_cond.flatten()
    y = t_evap.flatten()
    A = np.column_stack([
        np.ones_like(x),           # c0
        x,                         # c1 * t_evap
        y,                         # c2 * t_cond
        x ** 2,                    # c3 * t_evap^2
        x * y,                     # c4 * t_evap*t_cond
        y ** 2,                    # c5 * t_cond^2
        x ** 3,                    # c6 * t_evap^3
        (x ** 2) * y,              # c7 * t_evap^2*t_cond
        x * (y ** 2),              # c8 * t_evap*t_cond^2
        y ** 3                     # c9 * t_cond^3
    ])
    c, _, _, _ = np.linalg.lstsq(A, z, rcond=None)
    if check_diff:
        np.testing.assert_allclose(
            calc_EN12900(c, x, y), z,
            rtol=0.05
        )
    return c


def generate_eta_polys_from_power_and_cooling_polys(power_poly: list, cooling_poly: list, t_evap: np.array, t_cond: np.array, fluid: str, reference_state: dict) -> tuple:
    """Generate polynomials for calculation of isentropic and volumetric
    efficiency of a compressor

    Parameters
    ----------
    power_poly : list
        List of polynomial coefficients for power
    cooling_poly : list
        List of polynomial coefficients for cooling
    fluid : str
        Name of fluid
    reference_state : dict
        Dictionary with reference state information, i.e.
        - T_sh: superheating delta T (Kelvin)
        - T_sc: subcooling delta T (Kelvin)
        - rpm_poly: reference rpm of the polynomial
        - displancement: displacement in m3/h
        - rpm_displacement: reference rpm of the displacement

    Returns
    -------
    tuple
        Polynomial coefficients for isentropic and volumetric efficiency as
        function of evaporation and condensation temperature
    """
    columns = t_evap
    index = t_cond

    t_evap, t_cond = np.meshgrid(t_evap, t_cond)
    t_evap = t_evap.flatten() + 273.15
    t_cond = t_cond.flatten() + 273.15
    etas = _calc_etas_from_polynome(
        fluid,
        t_evap,
        t_cond,
        reference_state=reference_state,
        polynomes={"power": power_poly, "cooling": cooling_poly}
    )
    eta_s_poly = fit_EN12900(columns, index, etas["eta_s"].reshape(3, 6))
    eta_vol_poly = fit_EN12900(columns, index, etas["eta_vol"].reshape(3, 6))
    return eta_s_poly, eta_vol_poly


def generate_eta_polys_from_data(df_power, df_cooling, fluid: str, reference_state: dict) -> tuple:
    """Generate polynomials for calculation of isentropic and volumetric
    efficiency of a compressor

    Parameters
    ----------
    df_power : pd.DataFrame
        Power consumption data
    df_cooling : pd.DataFrame
        Cooling data
    fluid : str
        Name of fluid
    reference_state : dict
        Dictionary with reference state information, i.e.
        - T_sh: superheating delta T (Kelvin)
        - T_sc: subcooling delta T (Kelvin)
        - rpm_poly: reference rpm of the polynomial
        - displancement: displacement in m3/h
        - rpm_displacement: reference rpm of the displacement

    Returns
    -------
    tuple
        Polynomial coefficients for isentropic and volumetric efficiency as
        function of evaporation and condensation temperature
    """
    same_columns = (df_cooling.columns == df_power.columns).all()
    same_index = (df_cooling.index == df_power.index).all()
    if not same_columns or not same_index:
        msg = (
            "The dataframes provided to this method must have identical "
            "columns (evaporation temperature) and indices (condensation "
            "temperature)."
        )
        raise ValueError(msg)
    t_evap = df_power.columns
    t_cond = df_power.index
    power_poly = fit_EN12900(t_evap, t_cond, df_power.values)
    cooling_poly = fit_EN12900(t_evap, t_cond, df_cooling.values)

    return generate_eta_polys_from_power_and_cooling_polys(
        power_poly, cooling_poly, df_power.columns, df_power.index, fluid, reference_state
    )
