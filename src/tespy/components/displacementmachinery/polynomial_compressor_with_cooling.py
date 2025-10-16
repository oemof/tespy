from tespy.components.displacementmachinery.polynomial_compressor import PolynomialCompressor
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.helpers import TESPyComponentError


class PolynomialCompressorWithCooling(PolynomialCompressor):
    r"""
    Class for a compressor model following the EN12900 implementation of
    :cite:`cecchinato2010` and adding an inflow and and outflow for a cooling
    fluid.

    See the example for the intended use of the component.

    **Mandatory Equations**

    - fluid: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - mass flow: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - cooling energy balance: :py:meth:`tespy.components.displacementmachinery.polynomial_compressor_with_cooling.PolynomialCompressorWithCooling.cooling_energy_balance_func`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.dp_structure_matrix`
    - :py:meth:`tespy.components.component.Component.pr_structure_matrix`
    - :py:meth:`tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor.energy_balance_group_func`
    - :py:meth:`tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor.eta_s_group_func`
    - :py:meth:`tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor.eta_vol_group_func`
    - :py:meth:`tespy.components.component.Component.dp_structure_matrix` for cooling
    - :py:meth:`tespy.components.component.Component.pr_structure_matrix` for cooling

    Inlets/Outlets

    - in1, in2 (cooling)
    - out1, out2 (cooling)

    Optional inlets

    - power

    Image

    .. image:: /api/_images/PolynomialCompressorWithCooling.svg
       :alt: flowsheet of the compressor
       :align: center
       :class: only-light

    .. image:: /api/_images/PolynomialCompressorWithCooling_darkmode.svg
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

    dissipation_ratio : float, dict
        Relative heat loss of compressor, :math:`Q_\text{diss,rel}/1`

    eta_recovery : float, dict
        Share of heat recovered in the cooling fluid of the heat loss of
        compressor, :math:`Q_\text{diss,rel}/1`

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

    pr_cooling : float, dict
        Outlet to inlet pressure ratio for cooling, :math:`pr/1`

    dp_cooling : float, dict
        Inlet to outlet pressure difference for cooling,
        :math:`dp/\text{p}_\text{unit}`
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
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC"
    ... })
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
    :py:func:`tespy.components.displacementmachinery.polynomial_compressor.generate_eta_polys_from_power_and_cooling_polys`
    method, that can do the same step provided a polynomial for power and one
    for the cooling.

    >>> from tespy.components.displacementmachinery.polynomial_compressor import (
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
    ...     eta_s_poly=eta_s_poly, eta_vol_poly=eta_vol_poly, dissipation_ratio=0.05,
    ...     reference_state=reference_state
    ... )

    First, we can impose the boundary conditions on "c1" that are equal to the
    displacement reference state. In that case, we should be able to get the
    same displacement value as inputted into the reference.

    >>> c1.set_attr(fluid={"R134a": 1}, T=0, td_dew=10)  # T_evap=-10°C
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

    >>> from tespy.components.displacementmachinery.polynomial_compressor import (
    ...     calc_EN12900
    ... )
    >>> round(compressor.eta_s.val, 3) == round(calc_EN12900(eta_s_poly, -10, 50), 3)
    np.True_
    >>> round(compressor.eta_vol.val, 3) == round(calc_EN12900(eta_vol_poly, -10, 50), 3)
    np.True_

    .. tip::

        You can also create polynomials for power and cooling from respective
        data. For that, import the
        :py:func:`tespy.components.displacementmachinery.polynomial_compressor.fit_EN12900`
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

    >>> c1.set_attr(T=20, td_dew=10)  # T_evap=10°C
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

    def _preprocess(self, row_idx):
        if not self.eta_recovery.is_set:
            msg = (
                f"The component {self.label} of type {self.__class__.__name__}"
                "requires you to specify the share of heat recovery "
                "eta_recovery."
            )
            raise TESPyComponentError(msg)

        return super()._preprocess(row_idx)

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1', 'out2']

    def get_mandatory_constraints(self) -> dict:
        constraints = super().get_mandatory_constraints()
        # this is a dictionary
        constraints["cooling_energy_balance_constraints"] = dc_cmc(
            func=self.cooling_energy_balance_func,
            dependents=self.cooling_energy_balance_dependents,
            num_eq_sets=1
        )
        return constraints

    def get_parameters(self):
        params = super().get_parameters()
        params["eta_recovery"] = dc_cp(
            quantity="efficiency"
        )
        params["td_minimal"] = dc_cp(
            min_val=0,
            quantity="temperature_difference"
        )
        params["dp_cooling"] = dc_cp(
            min_val=0,
            structure_matrix=self.dp_structure_matrix,
            func_params={"inconn": 1, "outconn": 1, "dp": "dp_cooling"},
            quantity="pressure"
        )
        params["pr_cooling"] = dc_cp(
            min_val=0,
            structure_matrix=self.pr_structure_matrix,
            func_params={"inconn": 1, "outconn": 1, "pr": "pr_cooling"},
            quantity="ratio"
        )
        return params

    def cooling_energy_balance_func(self):
        r"""Energy balance equation for the cooling port

        Returns
        -------
        float
            residual of equation

            .. math::

                0 = \dot m_\text{in,2} \cdot \left( h_\text{out,2} - h_\text{in,2}\right)
                + \dot m_\text{in,1} \cdot \left(
                h_\text{out,1} - \frac{h_\text{out,1}}{1 - \text{diss_ratio}}
                + \frac{h_\text{in,1}\cdot\text{diss_ratio}}{1 - \text{diss_ratio}}
                \right)
        """
        residual = (
            self.inl[1].m.val_SI * (self.outl[1].h.val_SI - self.inl[1].h.val_SI)
            + self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI
                - self.outl[0].h.val_SI / (1 - self.dissipation_ratio.val_SI)
                + self.inl[0].h.val_SI * (
                    self.dissipation_ratio.val_SI / (1 - self.dissipation_ratio.val_SI)
                )
            ) * self.eta_recovery.val_SI
        )
        return residual

    def cooling_energy_balance_dependents(self):
        return [
            self.inl[0].m, self.inl[1].m,
            self.inl[0].h, self.inl[1].h,
            self.outl[0].h, self.outl[1].h
        ]

    def calc_parameters(self):
        super().calc_parameters()

        i = self.inl[0]
        o = self.outl[0]
        h_2 = (
            (o.h.val_SI - i.h.val_SI * self.dissipation_ratio.val_SI)
            / (1 - self.dissipation_ratio.val_SI)
        )
        T_max_compressor_internal = T_mix_ph(
            self.outl[0].p.val_SI,
            h_2,
            self.outl[0].fluid_data,
            self.outl[0].mixing_rule,
            T0=self.outl[0].T.val_SI
        )
        self.td_minimal.val_SI = (
            T_max_compressor_internal
            - self.outl[1].T.val_SI
        )

        self.dp_cooling.val_SI = self.inl[1].p.val_SI - self.outl[1].p.val_SI
        self.pr_cooling.val_SI = self.outl[1].p.val_SI / self.inl[1].p.val_SI
