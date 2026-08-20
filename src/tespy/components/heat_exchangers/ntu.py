# -*- coding: utf-8

"""Module of class NTUHeatExchanger.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/heat_exchangers/ntu.py

SPDX-License-Identifier: MIT
"""
import math

import numpy as np
from scipy.optimize import brentq
from scipy.optimize import minimize_scalar

from tespy.components.component import component_registry
from tespy.components.heat_exchangers.base import HeatExchanger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentCharacteristics as dc_gcc
from tespy.tools.data_containers import SimpleDataContainer as dc_simple

FLOW_ARRANGEMENTS = (
    "counterflow",
    "parallelflow",
    "crossflow_both_unmixed",
    "crossflow_hot_mixed",
    "crossflow_cold_mixed",
    "crossflow_both_mixed",
    "shell_and_tube"
)
_NTU_MAX = 1e4


def calc_epsilon(ntu, C_r, arrangement, num_shell_passes=1):
    r"""Calculate heat exchanger effectiveness from NTU and capacity ratio.

    The relations are taken from :cite:`Shah2003`, the approximation for
    crossflow with both fluids unmixed from :cite:`Bergman2011`.

    In case one side is isothermal, the capacity ratio is zero and the relation
    below is used:

    .. math::

        \varepsilon = 1 - \exp\left(-\text{NTU}\right)

    Counterflow:

    .. math::

        \varepsilon = \frac{1 - \exp\left(-\text{NTU} \cdot
        \left(1 - C_r\right)\right)}{1 - C_r \cdot \exp\left(-\text{NTU}
        \cdot \left(1 - C_r\right)\right)}

    Parallelflow:

    .. math::

        \varepsilon = \frac{1 - \exp\left(-\text{NTU} \cdot
        \left(1 + C_r\right)\right)}{1 + C_r}

    Crossflow with both fluids unmixed (approximation):

    .. math::

        \varepsilon = 1 - \exp\left(\frac{\text{NTU}^{0.22}}{C_r} \cdot
        \left(\exp\left(-C_r \cdot \text{NTU}^{0.78}\right) - 1\right)\right)

    Crossflow with the maximum capacity rate fluid mixed:

    .. math::

        \varepsilon = \frac{1 - \exp\left(-C_r \cdot \left(1 -
        \exp\left(-\text{NTU}\right)\right)\right)}{C_r}

    Crossflow with the minimum capacity rate fluid mixed:

    .. math::

        \varepsilon = 1 - \exp\left(-\frac{1 - \exp\left(-C_r \cdot
        \text{NTU}\right)}{C_r}\right)

    Crossflow with both fluids mixed:

    .. math::

        \varepsilon = \left(\frac{1}{1 - \exp\left(-\text{NTU}\right)} +
        \frac{C_r}{1 - \exp\left(-C_r \cdot \text{NTU}\right)} -
        \frac{1}{\text{NTU}}\right)^{-1}

    Shell and tube (TEMA E, one shell pass with 2, 4, ... tube passes):

    .. math::

        \varepsilon_1 = 2 \cdot \left(\left(1 + C_r\right) +
        \sqrt{1 + C_r^2} \cdot \frac{1 + \exp\left(-\text{NTU}_1 \cdot
        \sqrt{1 + C_r^2}\right)}{1 - \exp\left(-\text{NTU}_1 \cdot
        \sqrt{1 + C_r^2}\right)}\right)^{-1}

    For :math:`N` shell passes in series with :math:`\text{NTU}_1 =
    \text{NTU} / N` per shell pass:

    .. math::

        \varepsilon = \frac{\left(\frac{1 - \varepsilon_1 \cdot C_r}
        {1 - \varepsilon_1}\right)^N - 1}{\left(\frac{1 - \varepsilon_1
        \cdot C_r}{1 - \varepsilon_1}\right)^N - C_r}

    Parameters
    ----------
    ntu : float
        Number of transfer units.

    C_r : float
        Minimum to maximum capacity rate ratio.

    arrangement : str
        Flow arrangement, one of :code:`'counterflow'`,
        :code:`'parallelflow'`, :code:`'crossflow_both_unmixed'`,
        :code:`'crossflow_min_mixed'`, :code:`'crossflow_max_mixed'`,
        :code:`'crossflow_both_mixed'`, :code:`'shell_and_tube'`.

    num_shell_passes : int
        Number of shell passes for the shell and tube arrangement.

    Returns
    -------
    float
        Heat exchanger effectiveness.
    """
    if ntu <= 0:
        return 0.0
    if C_r < 1e-9:
        return 1 - math.exp(-ntu)
    if arrangement == "counterflow":
        if C_r > 1 - 1e-9:
            return ntu / (1 + ntu)
        e = math.exp(-ntu * (1 - C_r))
        return (1 - e) / (1 - C_r * e)
    elif arrangement == "parallelflow":
        return (1 - math.exp(-ntu * (1 + C_r))) / (1 + C_r)
    elif arrangement == "crossflow_both_unmixed":
        return 1 - math.exp(
            ntu ** 0.22 / C_r * (math.exp(-C_r * ntu ** 0.78) - 1)
        )
    elif arrangement == "crossflow_max_mixed":
        return (1 - math.exp(-C_r * (1 - math.exp(-ntu)))) / C_r
    elif arrangement == "crossflow_min_mixed":
        return 1 - math.exp(-(1 - math.exp(-C_r * ntu)) / C_r)
    elif arrangement == "crossflow_both_mixed":
        return 1 / (
            1 / (1 - math.exp(-ntu))
            + C_r / (1 - math.exp(-C_r * ntu))
            - 1 / ntu
        )
    elif arrangement == "shell_and_tube":
        ntu_1 = ntu / num_shell_passes
        S = (1 + C_r ** 2) ** 0.5
        e = math.exp(-ntu_1 * S)
        eps_1 = 2 / ((1 + C_r) + S * (1 + e) / (1 - e))
        if num_shell_passes == 1:
            return eps_1
        if C_r > 1 - 1e-9:
            return (
                num_shell_passes * eps_1
                / (1 + (num_shell_passes - 1) * eps_1)
            )
        R = ((1 - eps_1 * C_r) / (1 - eps_1)) ** num_shell_passes
        return (R - 1) / (R - C_r)
    else:
        msg = f"The flow arrangement '{arrangement}' is not available."
        raise ValueError(msg)


def calc_ntu(epsilon, C_r, arrangement, num_shell_passes=1):
    r"""Calculate NTU from heat exchanger effectiveness and capacity ratio.

    Inverse of :py:func:`calc_epsilon`, see there for the underlying relations
    and the available flow arrangements. For the crossflow arrangements with
    both fluids unmixed or both fluids mixed the inversion is performed
    numerically. If the effectiveness exceeds the value reachable with the
    specified arrangement, :code:`nan` is returned.

    Parameters
    ----------
    epsilon : float
        Heat exchanger effectiveness.

    C_r : float
        Minimum to maximum capacity rate ratio.

    arrangement : str
        Flow arrangement, see :py:func:`calc_epsilon`.

    num_shell_passes : int
        Number of shell passes for the shell and tube arrangement.

    Returns
    -------
    float
        Number of transfer units.
    """
    if epsilon <= 0:
        return 0.0
    try:
        if C_r < 1e-9:
            return -math.log(1 - epsilon)
        if arrangement == "counterflow":
            if C_r > 1 - 1e-9:
                return epsilon / (1 - epsilon)
            return (
                math.log((1 - epsilon * C_r) / (1 - epsilon)) / (1 - C_r)
            )
        elif arrangement == "parallelflow":
            return -math.log(1 - epsilon * (1 + C_r)) / (1 + C_r)
        elif arrangement == "crossflow_max_mixed":
            return -math.log(1 + math.log(1 - epsilon * C_r) / C_r)
        elif arrangement == "crossflow_min_mixed":
            return -math.log(1 + C_r * math.log(1 - epsilon)) / C_r
        elif arrangement in ("crossflow_both_unmixed", "crossflow_both_mixed"):
            def residual(ntu):
                return calc_epsilon(ntu, C_r, arrangement) - epsilon

            # the both mixed relation decreases again beyond its maximum,
            # invert on the ascending branch only
            ntu_upper = _NTU_MAX
            if arrangement == "crossflow_both_mixed":
                res = minimize_scalar(
                    lambda ntu: -calc_epsilon(ntu, C_r, arrangement),
                    bounds=(1e-12, _NTU_MAX), method="bounded"
                )
                ntu_upper = res.x
            if residual(ntu_upper) < 0:
                return np.nan
            return brentq(residual, 1e-12, ntu_upper)
        elif arrangement == "shell_and_tube":
            if num_shell_passes == 1:
                eps_1 = epsilon
            elif C_r > 1 - 1e-9:
                eps_1 = epsilon / (
                    num_shell_passes - (num_shell_passes - 1) * epsilon
                )
            else:
                F = (
                    (epsilon * C_r - 1) / (epsilon - 1)
                ) ** (1 / num_shell_passes)
                eps_1 = (F - 1) / (F - C_r)
            S = (1 + C_r ** 2) ** 0.5
            E = (2 / eps_1 - (1 + C_r)) / S
            if E <= 1:
                return np.nan
            return -num_shell_passes * math.log((E - 1) / (E + 1)) / S
        else:
            msg = f"The flow arrangement '{arrangement}' is not available."
            raise ValueError(msg)
    except (ValueError, ZeroDivisionError) as err:
        if "flow arrangement" in str(err):
            raise
        return np.nan


@component_registry
class NTUHeatExchanger(HeatExchanger):
    r"""
    Class for a heat exchanger with effectiveness-NTU based heat transfer.

    In contrast to the :py:class:`HeatExchanger
    <tespy.components.heat_exchangers.base.HeatExchanger>` class, which links
    the heat transfer coefficient :code:`UA` to the heat transferred through
    the logarithmic mean temperature difference, this component applies the
    effectiveness-NTU relation :cite:`Shah2003` based on the flow arrangement
    specified by the user:

    .. math::

        \dot{Q} = \varepsilon \cdot C_\text{min} \cdot
        \left(T_{in,1} - T_{in,2}\right) \text{ with }
        \varepsilon = f\left(\text{NTU}, C_r\right) \text{ and }
        \text{NTU} = \frac{UA}{C_\text{min}}

    The capacity rates :math:`C` are calculated from the enthalpy and
    temperature changes of the respective side, e.g.
    :math:`C_\text{hot}=\dot{m}_{in,1} \cdot \frac{h_{in,1} - h_{out,1}}
    {T_{in,1} - T_{out,1}}`. In case of equal inlet and outlet temperature the
    capacity rate is set to infinite (e.g. on a phase change of a pure fluid).
    The available relations are listed in :py:func:`calc_epsilon
    <tespy.components.heat_exchangers.ntu.calc_epsilon>`. For the crossflow
    arrangements with one fluid mixed the mixed side is specified physically
    (hot or cold), the assignment to the minimum or maximum capacity rate
    relation happens automatically.

    .. note::

        The effectiveness-NTU relations assume constant capacity rates along
        the heat exchanger. If multiple phases are present on one side of the
        heat exchanger use the
        :py:class:`MovingBoundaryHeatExchanger
        <tespy.components.heat_exchangers.movingboundary.MovingBoundaryHeatExchanger>`
        instead.

    .. note::

        The implemented equations assume isothermal phase changes. A pressure
        drop on a condensing or evaporating side lowers its saturation
        temperature. Especially in evaporation processes this leads to a large
        error. Therefore, either specify no pressure drop on a phase changing
        side of this component, or use the
        :py:class:`MovingBoundaryHeatExchanger
        <tespy.components.heat_exchangers.movingboundary.MovingBoundaryHeatExchanger>`
        or :py:class:`SectionedHeatExchanger
        <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger>`.

        Corrections to correct the NTU approach for mixtures with temperature
        glide and pressure drops exist, e.g. in :cite:`Qiao2025`, but are not
        implemented in tespy yet.

    .. attention::

        The terminal temperature differences :code:`ttd_u` and :code:`ttd_l`
        will output the results in the standard counterflow convention
        independent of flow arrangement specified here.

    .. image:: /api/_images/components/HeatExchanger.svg
       :alt: flowsheet of the ntuheatexchanger
       :align: center
       :class: only-light

    .. image:: /api/_images/components/HeatExchanger_darkmode.svg
       :alt: flowsheet of the ntuheatexchanger
       :align: center
       :class: only-dark

    Ports
    -----

    - Fluid inlets: in1, in2
    - Fluid outlets: out1, out2

    Mandatory Equations
    -------------------

    - mass flow equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
    - fluid composition equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
    - hot side to cold side heat transfer equation: :py:meth:`energy_balance_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_func>`

    Parameters
    ----------

    C_r : float, dict
        Minimum to maximum capacity rate ratio. Quantity: :code:`ratio`.

    char_warnings : bool
        Ignore warnings on default characteristics usage for this component.

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the components design case.

    dp1 : float, dict
        Hot side inlet to outlet absolute pressure change. Quantity:
        :code:`pressure_difference`.
        Equation: :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`.

    dp2 : float, dict
        Cold side inlet to outlet absolute pressure change. Quantity:
        :code:`pressure_difference`.
        Equation: :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`.

    eff_cold : float, dict
        Heat exchanger effectiveness for cold side. Quantity:
        :code:`efficiency`.
        Equation: :py:meth:`eff_cold_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_cold_func>`.

    eff_hot : float, dict
        Heat exchanger effectiveness for hot side. Quantity: :code:`efficiency`.
        Equation: :py:meth:`eff_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_hot_func>`.

    eff_max : float, dict
        Maximum heat exchanger effectiveness. Quantity: :code:`efficiency`.
        Equation: :py:meth:`eff_max_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_max_func>`.

    eps : float, dict
        Heat exchanger effectiveness, heat transfer over maximum transferable
        heat flow at constant capacity rates. Quantity: :code:`efficiency`.
        Equation: :py:meth:`eps_func <tespy.components.heat_exchangers.ntu.NTUHeatExchanger.eps_func>`.

    flow_arrangement : str
        Flow arrangement for the effectiveness-NTU relation, mandatory, one of
        :code:`'counterflow'`, :code:`'parallelflow'`,
        :code:`'crossflow_both_unmixed'`, :code:`'crossflow_hot_mixed'`,
        :code:`'crossflow_cold_mixed'`, :code:`'crossflow_both_mixed'`,
        :code:`'shell_and_tube'`.

    kA : float, dict
        Deprecated, use :code:`UA` instead. Quantity:
        :code:`heat_transfer_coefficient`.

    kA_char : GroupedComponentCharacteristics
        Deprecated, use :code:`UA_char` instead. Elements: :code:`kA_char1`,
        :code:`kA_char2`.

    kA_char1 : tespy.tools.characteristics.CharLine, dict
        Deprecated, use :code:`UA_char1` instead.

    kA_char2 : tespy.tools.characteristics.CharLine, dict
        Deprecated, use :code:`UA_char2` instead.

    label : str
        The label of the component.

    lmtd : float, dict
        Effective logarithmic mean temperature difference :code:`Q/UA`.
        Quantity: :code:`temperature_difference`.

    lmtd_per_section : numpy.ndarray
        Logarithmic mean temperature difference in each section. Quantity:
        :code:`temperature_difference`. Result only - populated by the network
        after each solve.

    local_design : bool
        Treat this component in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this component in offdesign mode in a design calculation.

    NTU : float, dict
        Number of transfer units :code:`UA` over minimum capacity rate.
        Equation: :py:meth:`NTU_func <tespy.components.heat_exchangers.ntu.NTUHeatExchanger.NTU_func>`.

    num_shell_passes : int
        Number of shell passes for the shell and tube flow arrangement (with 2,
        4, 6, ... tube passes per shell pass).

    offdesign : list
        List containing offdesign parameters (stated as String).

    pr1 : float, dict
        Hot side outlet to inlet pressure ratio. Quantity: :code:`ratio`.
        Equation: :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`.

    pr2 : float, dict
        Cold side outlet to inlet pressure ratio. Quantity: :code:`ratio`.
        Equation: :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`.

    printout : bool
        Include this component in the network's results printout.

    Q : float, dict
        Heat transfer from hot side. Quantity: :code:`heat`.
        Equation: :py:meth:`energy_balance_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func>`.

    Q_per_section : numpy.ndarray
        Heat transferred from hot to cold side in each section. Quantity:
        :code:`heat`. Result only - populated by the network after each solve.

    Q_sections : numpy.ndarray
        Cumulative heat transferred from hot to cold side up to each section
        boundary. Quantity: :code:`heat`. Result only - populated by the network
        after each solve.

    T_cold_sections : numpy.ndarray
        Cold side temperature at each section boundary. Quantity:
        :code:`temperature`. Result only - populated by the network after each
        solve.

    T_hot_sections : numpy.ndarray
        Hot side temperature at each section boundary. Quantity:
        :code:`temperature`. Result only - populated by the network after each
        solve.

    td_log : float, dict
        Deprecated, use :code:`lmtd` instead. Quantity:
        :code:`temperature_difference`.

    ttd_l : float, dict
        Terminal temperature difference at hot side outlet to cold side inlet.
        Quantity: :code:`temperature_difference`.
        Equation: :py:meth:`ttd_l_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_l_func>`.

    ttd_min : float, dict
        Minimum terminal temperature difference. Quantity:
        :code:`temperature_difference`.
        Equation: :py:meth:`ttd_min_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_min_func>`.

    ttd_u : float, dict
        Terminal temperature difference at hot side inlet to cold side outlet.
        Quantity: :code:`temperature_difference`.
        Equation: :py:meth:`ttd_u_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_u_func>`.

    UA : float, dict
        Heat transfer coefficient linking heat transfer and effectiveness
        through the NTU relation. Quantity: :code:`heat_transfer_coefficient`.
        Equation: :py:meth:`UA_func <tespy.components.heat_exchangers.ntu.NTUHeatExchanger.UA_func>`.

    UA_char : GroupedComponentCharacteristics
        Equation for heat transfer based on UA and modification factor.
        Elements: :code:`UA_char1`, :code:`UA_char2`.
        Equation: :py:meth:`UA_char_func <tespy.components.heat_exchangers.ntu.NTUHeatExchanger.UA_char_func>`.

    UA_char1 : tespy.tools.characteristics.CharLine, dict
        Hot side UA modification lookup table for offdesign.

    UA_char2 : tespy.tools.characteristics.CharLine, dict
        Cold side UA modification lookup table for offdesign.

    zeta1 : float, dict
        Deprecated, use :code:`zeta1_d4` instead.

    zeta1_d4 : float, dict
        Hot side geometry-independent friction coefficient zeta/D^4 for pressure
        loss calculation.
        Equation: :py:meth:`zeta_d4_func <tespy.components.component.Component.zeta_d4_func>`.

    zeta2 : float, dict
        Deprecated, use :code:`zeta2_d4` instead.

    zeta2_d4 : float, dict
        Cold side geometry-independent friction coefficient zeta/D^4 for
        pressure loss calculation.
        Equation: :py:meth:`zeta_d4_func <tespy.components.component.Component.zeta_d4_func>`.

    Example
    -------
    Hot air is cooled with water in a crossflow heat exchanger with both fluids
    unmixed. In the design case we specify both mass flows as well as both
    inlet temperatures and the air outlet temperature. The water outflow
    temperature, :code:`NTU`, effectiveness :code:`eps` and :code:`UA` are
    results.

    >>> from tespy.components import Sink, Source, NTUHeatExchanger
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC", "heat_transfer_coefficient": "kW/K"
    ... })
    >>> air_in = Source("air inlet")
    >>> air_out = Sink("air outlet")
    >>> water_in = Source("water inlet")
    >>> water_out = Sink("water outlet")
    >>> he = NTUHeatExchanger("crossflow cooler")
    >>> he.set_attr(flow_arrangement="crossflow_both_unmixed")
    >>> c1 = Connection(air_in, "out1", he, "in1", label="c1")
    >>> c2 = Connection(he, "out1", air_out, "in1", label="c2")
    >>> c3 = Connection(water_in, "out1", he, "in2", label="c3")
    >>> c4 = Connection(he, "out2", water_out, "in1", label="c4")
    >>> nw.add_conns(c1, c2, c3, c4)
    >>> c1.set_attr(fluid={"air": 1}, m=1, T=85, p=1)
    >>> c2.set_attr(T=45)
    >>> c3.set_attr(fluid={"water": 1}, m=3, T=25, p=3)
    >>> he.set_attr(pr1=0.98, pr2=0.98)
    >>> nw.solve("design")

    We can check the NTU specific results:

    >>> round(he.eps.val, 3)
    0.667
    >>> round(he.NTU.val, 4)
    1.1486
    >>> round(he.UA.val, 4)
    1.1581

    Now, we can calculate the heat transfer by specifying any of these three
    results instead of specifying the air outlet temperature. All three are
    alternative formulations of the same equation, therefore we achieve
    identical results:

    >>> c2.set_attr(T=None)
    >>> he.set_attr(NTU=1.1486)
    >>> nw.solve("design")
    >>> round(c2.T.val, 1)
    45.0

    >>> he.set_attr(NTU=None, eps=0.6667)
    >>> nw.solve("design")
    >>> round(c2.T.val, 1)
    45.0

    >>> he.set_attr(eps=None, UA=1.1581)
    >>> nw.solve("design")
    >>> round(c2.T.val, 1)
    45.0

    For an offdesign calculation :code:`UA` remains specified as input.
    Consequently, :code:`NTU` and :code:`eps` are results. When reducing air
    mass flow the air outlet temperature will decrease and the effectiveness
    will increase.

    >>> design_state = nw.save(as_dict=True)
    >>> c1.set_attr(m=0.75)
    >>> nw.solve("offdesign", design_path=design_state)
    >>> round(he.eps.val, 3) > 0.667
    True
    >>> round(c2.T.val, 1)
    38.8
    """

    def get_parameters(self):
        params = super().get_parameters()
        params["flow_arrangement"] = dc_simple(
            dtype="str",
            description=(
                "flow arrangement for the effectiveness-NTU relation, "
                "mandatory, one of "
                + ", ".join(f":code:`'{a}'`" for a in FLOW_ARRANGEMENTS)
            )
        )
        params["num_shell_passes"] = dc_simple(
            val=1, dtype="int",
            description=(
                "number of shell passes for the shell and tube flow "
                "arrangement (with 2, 4, 6, ... tube passes per shell pass)"
            )
        )
        params["UA"] = dc_cp(
            min_val=0, num_eq_sets=1,
            func=self.UA_func,
            dependents=self.UA_dependents,
            quantity="heat_transfer_coefficient",
            description=(
                "heat transfer coefficient linking heat transfer and "
                "effectiveness through the NTU relation"
            ),
            calc=self._calc_UA, calc_deps=["Q"]
        )
        params["UA_char"] = dc_gcc(
            elements=["UA_char1", "UA_char2"],
            num_eq_sets=1,
            func=self.UA_char_func,
            dependents=self.UA_char_dependents,
            description=(
                "equation for heat transfer based on UA and modification "
                "factor"
            )
        )
        params["NTU"] = dc_cp(
            min_val=0, num_eq_sets=1,
            func=self.NTU_func,
            dependents=self.UA_dependents,
            description=(
                "number of transfer units :code:`UA` over minimum capacity "
                "rate"
            ),
            calc=self._calc_NTU, calc_deps=["UA"]
        )
        params["C_r"] = dc_cp(
            min_val=0, max_val=1, is_result=True, quantity="ratio",
            description="minimum to maximum capacity rate ratio",
            calc=self._calc_C_r
        )
        params["eps"] = dc_cp(
            min_val=0, max_val=1, num_eq_sets=1, quantity="efficiency",
            func=self.eps_func,
            dependents=self.UA_dependents,
            description=(
                "heat exchanger effectiveness, heat transfer over maximum "
                "transferable heat flow at constant capacity rates"
            ),
            calc=self._calc_eps, calc_deps=["Q"]
        )
        return params

    def _preprocess(self, row_idx):
        if not self.flow_arrangement.is_set:
            msg = (
                f"The flow arrangement of {self.label} is not specified. The "
                "effectiveness-NTU relation depends on it, you need to set it "
                "explicitly. Select from '"
                + "', '".join(FLOW_ARRANGEMENTS) + "'."
            )
            raise ValueError(msg)
        if self.flow_arrangement.val not in FLOW_ARRANGEMENTS:
            msg = (
                f"The flow arrangement '{self.flow_arrangement.val}' of "
                f"component {self.label} is not available. Available flow "
                "arrangements are: '" + "', '".join(FLOW_ARRANGEMENTS) + "'."
            )
            raise ValueError(msg)
        if (
                self.flow_arrangement.val == "shell_and_tube"
                and self.num_shell_passes.val < 1
        ):
            msg = (
                f"The number of shell passes of component {self.label} must "
                "be at least 1."
            )
            raise ValueError(msg)
        super()._preprocess(row_idx)

    def _calc_capacity_rates(self):
        r"""Calculate the capacity rates of both sides.

        The capacity rates are effective values based on the enthalpy and
        temperature change of the respective side. A side with (nearly)
        equal inlet and outlet temperature is considered isothermal, its
        capacity rate is infinite.

        Returns
        -------
        tuple
            Hot side inlet temperature, cold side inlet temperature, hot
            side capacity rate, cold side capacity rate.
        """
        i1, o1 = self.inl[0], self.outl[0]
        i2, o2 = self.inl[1], self.outl[1]
        T_i1 = i1.calc_T()
        T_o1 = o1.calc_T()
        T_i2 = i2.calc_T()
        T_o2 = o2.calc_T()
        if abs(T_i1 - T_o1) > 1e-6:
            C_hot = abs(
                i1.m.val_SI * (i1.h.val_SI - o1.h.val_SI) / (T_i1 - T_o1)
            )
        else:
            C_hot = math.inf
        if abs(T_o2 - T_i2) > 1e-6:
            C_cold = abs(
                i2.m.val_SI * (o2.h.val_SI - i2.h.val_SI) / (T_o2 - T_i2)
            )
        else:
            C_cold = math.inf
        return T_i1, T_i2, C_hot, C_cold

    def _resolve_arrangement(self, C_hot, C_cold):
        arrangement = self.flow_arrangement.val
        if arrangement == "crossflow_hot_mixed":
            if C_hot <= C_cold:
                return "crossflow_min_mixed"
            return "crossflow_max_mixed"
        elif arrangement == "crossflow_cold_mixed":
            if C_cold <= C_hot:
                return "crossflow_min_mixed"
            return "crossflow_max_mixed"
        return arrangement

    def _ntu_residual(self, UA):
        T_i1, T_i2, C_hot, C_cold = self._calc_capacity_rates()
        Q = (
            self.inl[0].m.val_SI
            * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        )
        if math.isinf(C_hot) and math.isinf(C_cold):
            # both sides isothermal: constant temperature difference
            return Q + UA * (T_i1 - T_i2)
        C_min = min(C_hot, C_cold)
        if C_min == 0:
            return Q
        C_max = max(C_hot, C_cold)
        C_r = 0.0 if math.isinf(C_max) else C_min / C_max
        eps = calc_epsilon(
            UA / C_min, C_r, self._resolve_arrangement(C_hot, C_cold),
            self.num_shell_passes.val
        )
        return Q + eps * C_min * (T_i1 - T_i2)

    def UA_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in,1} \cdot \left(h_{out,1} - h_{in,1}\right)
                + \varepsilon \cdot C_\text{min} \cdot
                \left(T_{in,1} - T_{in,2}\right)

                \varepsilon = f\left(\text{NTU}, C_r\right),\;
                \text{NTU} = \frac{UA}{C_\text{min}},\;
                C_r = \frac{C_\text{min}}{C_\text{max}}
        """
        return self._ntu_residual(self.UA.val_SI)

    def NTU_func(self):
        r"""
        Calculate heat transfer from the number of transfer units

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in,1} \cdot \left(h_{out,1} - h_{in,1}\right)
                + \varepsilon \cdot C_\text{min} \cdot
                \left(T_{in,1} - T_{in,2}\right)

                \varepsilon = f\left(\text{NTU}, C_r\right),\;
                UA = \text{NTU} \cdot C_\text{min}
        """
        _, _, C_hot, C_cold = self._calc_capacity_rates()
        C_min = min(C_hot, C_cold)
        if math.isinf(C_min):
            # both sides isothermal, the number of transfer units is not
            # defined: drive the duty towards zero to leave that state
            return (
                self.inl[0].m.val_SI
                * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
            )
        return self._ntu_residual(self.NTU.val_SI * C_min)

    def eps_func(self):
        r"""
        Calculate heat transfer from the heat exchanger effectiveness.

        :code:`eff_hot` and :code:`eff_cold` relate the heat transferred to the
        maximum possible enthalpy change, this equation is the effectiveness of
        the effectiveness-NTU method.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in,1} \cdot \left(h_{out,1} - h_{in,1}\right)
                + \varepsilon \cdot C_\text{min} \cdot
                \left(T_{in,1} - T_{in,2}\right)
        """
        T_i1, T_i2, C_hot, C_cold = self._calc_capacity_rates()
        Q = (
            self.inl[0].m.val_SI
            * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)
        )
        if math.isinf(C_hot) and math.isinf(C_cold):
            return Q
        C_min = min(C_hot, C_cold)
        return Q + self.eps.val_SI * C_min * (T_i1 - T_i2)

    def UA_dependents(self):
        return [
            self.inl[0].m,
            self.inl[1].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
            self.inl[1].p,
            self.inl[1].h,
            self.outl[1].p,
            self.outl[1].h,
        ]

    def UA_char_func(self):
        r"""
        Calculate heat transfer from heat transfer coefficient characteristic.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_{in,1} \cdot \left(h_{out,1} - h_{in,1}\right)
                + \varepsilon \cdot C_\text{min} \cdot
                \left(T_{in,1} - T_{in,2}\right)

                \varepsilon = f\left(\frac{UA_\text{design} \cdot f_{UA}}
                {C_\text{min}}, C_r\right),\;
                f_{UA} = \frac{2}{\frac{1}{f_1\left(expr_1\right)} +
                \frac{1}{f_2\left(expr_2\right)}}

        Note
        ----
        For standard functions f\ :subscript:`1` \ and f\ :subscript:`2` \
        see module :ref:`tespy.data <data_label>`.
        """
        p1 = self.UA_char1.param
        p2 = self.UA_char2.param
        if self.local_offdesign:
            design_value = self._connection_offdesign[self.inl[0].label][p1]
            actual_value = getattr(self.inl[0], p1).val_SI
            f1 = actual_value / design_value

            design_value = self._connection_offdesign[self.inl[1].label][p2]
            actual_value = getattr(self.inl[1], p2).val_SI
            f2 = actual_value / design_value
        else:
            f1 = self.get_char_expr(p1, **self.UA_char1.char_params)
            f2 = self.get_char_expr(p2, **self.UA_char2.char_params)

        fUA1 = self.UA_char1.char_func.evaluate(f1)
        fUA2 = self.UA_char2.char_func.evaluate(f2)
        fUA = 2 / (1 / fUA1 + 1 / fUA2)

        return self._ntu_residual(self.UA.design * fUA)

    def UA_char_dependents(self):
        return self.UA_dependents()

    def _calc_eps(self):
        T_i1, T_i2, C_hot, C_cold = self._calc_capacity_rates()
        if math.isinf(C_hot) and math.isinf(C_cold):
            return np.nan
        C_min = min(C_hot, C_cold)
        return -self.Q.val_SI / (C_min * (T_i1 - T_i2))

    def _calc_C_r(self):
        _, _, C_hot, C_cold = self._calc_capacity_rates()
        if math.isinf(C_hot) and math.isinf(C_cold):
            return np.nan
        C_max = max(C_hot, C_cold)
        return 0.0 if math.isinf(C_max) else min(C_hot, C_cold) / C_max

    def _calc_NTU(self):
        _, _, C_hot, C_cold = self._calc_capacity_rates()
        C_min = min(C_hot, C_cold)
        if math.isinf(C_min):
            return np.nan
        return self.UA.val_SI / C_min

    def _calc_UA(self):
        T_i1, T_i2, C_hot, C_cold = self._calc_capacity_rates()
        if math.isinf(C_hot) and math.isinf(C_cold):
            return -self.Q.val_SI / (T_i1 - T_i2)
        C_min = min(C_hot, C_cold)
        C_max = max(C_hot, C_cold)
        C_r = 0.0 if math.isinf(C_max) else C_min / C_max
        eps = -self.Q.val_SI / (C_min * (T_i1 - T_i2))
        ntu = calc_ntu(
            eps, C_r, self._resolve_arrangement(C_hot, C_cold),
            self.num_shell_passes.val
        )
        return ntu * C_min
