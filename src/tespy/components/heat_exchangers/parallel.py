import math

import numpy as np

from tespy.components.component import component_registry

from .base import HeatExchanger


@component_registry
class ParallelFlowHeatExchanger(HeatExchanger):
    r"""
    Class for parallel flow heat exchanger.

    **Mandatory Equations**

    - fluid: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - mass flow: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_func`

    **Optional Equations**

    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.kA_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.kA_char_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.ttd_u_func`
    - :py:meth:`tespy.components.heat_exchangers.base.HeatExchanger.ttd_l_func`

    For hot and cold side individually:

    - :py:meth:`tespy.components.component.Component.pr_structure_matrix`
    - :py:meth:`tespy.components.component.Component.dp_structure_matrix`
    - :py:meth:`tespy.components.component.Component.zeta_func`

    Inlets/Outlets

    - in1, in2 (index 1: hot side, index 2: cold side)
    - out1, out2 (index 1: hot side, index 2: cold side)

    Image

    .. image:: /api/_images/HeatExchanger.svg
       :alt: flowsheet of the heat exchanger
       :align: center
       :class: only-light

    .. image:: /api/_images/HeatExchanger_darkmode.svg
       :alt: flowsheet of the heat exchanger
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

    pr1 : float, dict
        Outlet to inlet pressure ratio at hot side, :math:`pr/1`.

    pr2 : float, dict
        Outlet to inlet pressure ratio at cold side, :math:`pr/1`.

    dp1 : float, dict,
        Inlet to outlet pressure delta at hot side, :math:`dp/\text{Pa}`

    dp2 : float, dict
        Inlet to outlet pressure delta at cold side, :math:`dp\text{Pa}`.

    zeta1 : float, dict
        Geometry independent friction coefficient at hot side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    zeta2 : float, dict
        Geometry independent friction coefficient at cold side,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    ttd_l : float, dict
        Initial terminal temperature difference, referring to the temperature
        difference between the two inlets of the heat exchanger,
        :math:`ttd_\mathrm{l}/\text{K}`.

    ttd_u : float, dict
        Final terminal temperature difference, referring to the temperature
        difference between the two outlets of the heat exchanger,
        :math:`ttd_\mathrm{u}/\text{K}`.

    kA : float, dict
        Area independent heat transfer coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : dict
        Area independent heat transfer coefficient characteristic.

    kA_char1 : tespy.tools.characteristics.CharLine, dict
        Characteristic line for hot side heat transfer coefficient.

    kA_char2 : tespy.tools.characteristics.CharLine, dict
        Characteristic line for cold side heat transfer coefficient.

    Note
    ----
    The :code:`ParallelFlowHeatExchanger` implements parallel flow of both
    streams, meaning the streams enter with a high temperature difference and
    then gradually reduce their temperature difference to each other. The
    initial temperature difference is the maximum temperature difference, the
    final temperature difference is the minimum temperature difference.

    Example
    -------
    Water at 75 °C is used to heat up an air stream 2500 l/s from 10 °C to
    35 °C.

    >>> from tespy.components import Sink, Source, ParallelFlowHeatExchanger
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import os
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC", "enthalpy": "kJ/kg",
    ...     "volumetric_flow": "l/s", "heat_transfer_coefficient": "kW/K"
    ... })
    >>> feed_water = Source("Feed water inlet")
    >>> return_water = Sink("Water outlet")
    >>> air_inlet = Source("Fresh air inlet")
    >>> air_warm = Sink("Air outlet")
    >>> he = ParallelFlowHeatExchanger("heat exchanger")
    >>> c1 = Connection(feed_water, 'out1', he, 'in1')
    >>> c2 = Connection(he, 'out1', return_water, 'in1')
    >>> c3 = Connection(air_inlet, 'out1', he, 'in2')
    >>> c4 = Connection(he, 'out2', air_warm, 'in1')
    >>> nw.add_conns(c1, c2, c3, c4)

    We assume pressure losses of 10 mbar on the air side, and 100 mbar on the
    water side. Depending on our specifications we can calculate:

    - What is the temperature of the water leaving the heat exchanger, or
    - What is the required mass flow of water to heat up the air

    Let's first assume, that a final pinch :code:`ttd_u` of 7.5 K is desired,
    meaning, the water should leave the heat exchanger with a temperature
    higher than the air by 7.5 K. With that, we will get the water flow, and
    also the heat transfer coefficient :code:`kA` as a result.

    .. note::

        Note, that for specification of initial or final pinch temperature
        differences, it is often important to start the calculation with a good
        initial guess. In this case, we know that the water will be in liquid
        state, therefore we can

        - either specify the outlet enthalpy initial guess to be liquid,
        - the state at the water outlet to be :code:`"l"`,
        - or use :code:`"INCOMP::Water"` as fluid, since this uses liquid phase
          only.

    >>> he.set_attr(dp1=0.1, dp2=0.01, ttd_u=7.5)
    >>> c1.set_attr(fluid={"INCOMP::Water": 1}, T=70, p=1.3)
    >>> c3.set_attr(fluid={"air": 1}, T=10, p=1.02, v=2500)
    >>> c4.set_attr(T=35)
    >>> nw.solve("design")
    >>> round(c1.v.val, 2)
    0.7
    >>> round(he.kA.val, 2)
    3.13

    Now, it might be interesting to see what happens under different operation
    conditions after we have designed the system. For that, we can assume that
    the heat transfer coefficient is constant. First we just fix the :code:`kA`
    value instead of the final pinch and then resolve again.

    >>> he.set_attr(design=["ttd_u"], offdesign=["kA"])
    >>> nw.save("design.json")
    >>> nw.solve("offdesign", design_path="design.json")
    >>> round(he.kA.val_SI / he.kA.design, 1)
    1.0

    Now, let's see what happens under different operating conditions. First
    we change the air volumetric flow, then we change the air temperature to
    check what happens to the outflow temperature of the water.

    >>> c3.set_attr(v=2000)
    >>> nw.solve("offdesign", design_path="design.json")
    >>> round(c2.T.val, 2)
    38.69
    >>> c3.set_attr(v=2500, T=8)
    >>> nw.solve("offdesign", design_path="design.json")
    >>> round(c2.T.val, 2)
    44.0
    >>> os.remove("design.json")
    """
    def get_parameters(self):
        params = super().get_parameters()
        del params["ttd_min"]
        del params["eff_hot"]
        del params["eff_cold"]
        del params["eff_max"]
        return params

    def ttd_l_func(self):
        T_i1 = self.inl[0].calc_T()
        T_i2 = self.inl[1].calc_T()
        return self.ttd_l.val_SI - T_i1 + T_i2

    def ttd_l_dependents(self):
        return [var for c in self.inl for var in [c.p, c.h]]

    def ttd_u_func(self):
        T_o1 = self.outl[0].calc_T()
        T_o2 = self.outl[1].calc_T()
        return self.ttd_u.val_SI - T_o1 + T_o2

    def ttd_u_dependents(self):
        return [var for c in self.outl for var in [c.p, c.h]]

    def calculate_td_log(self):
        i1 = self.inl[0]
        i2 = self.inl[1]
        o1 = self.outl[0]
        o2 = self.outl[1]

        # temperature value manipulation for convergence stability
        T_i1 = i1.calc_T()
        T_i2 = i2.calc_T()
        T_o1 = o1.calc_T()
        T_o2 = o2.calc_T()

        if T_i1 <= T_i2:
            T_i1 = T_i2 + 0.01
        if T_o1 <= T_o2:
            T_o2 = T_o1 - 0.01

        ttd_u = T_o1 - T_o2
        ttd_l = T_i1 - T_i2

        if round(ttd_u, 6) == round(ttd_l, 6):
            td_log = ttd_l
        else:
            td_log = (ttd_l - ttd_u) / math.log((ttd_l) / (ttd_u))

        return td_log

    def calc_parameters(self):

        self.Q.val_SI = self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI
        )
        self.ttd_u.val_SI = self.outl[0].T.val_SI - self.outl[1].T.val_SI
        self.ttd_l.val_SI = self.inl[0].T.val_SI - self.inl[1].T.val_SI

        # pr and zeta
        for i in range(2):
            self.get_attr(f'pr{i + 1}').val_SI = (
                self.outl[i].p.val_SI / self.inl[i].p.val_SI
            )
            self.get_attr(f'zeta{i + 1}').val_SI = self.calc_zeta(
                self.inl[i], self.outl[i]
            )
            self.get_attr(f'dp{i + 1}').val_SI = (
                self.inl[i].p.val_SI - self.outl[i].p.val_SI
            )

        # kA and logarithmic temperature difference
        if self.ttd_u.val_SI < 0 or self.ttd_l.val_SI < 0:
            self.td_log.val_SI = np.nan
        elif round(self.ttd_l.val_SI, 6) == round(self.ttd_u.val_SI, 6):
            self.td_log.val_SI = self.ttd_l.val_SI
        elif round(self.ttd_l.val_SI, 6) == 0 or round(self.ttd_u.val_SI, 6) == 0:
            self.td_log.val_SI = np.nan
        else:
            self.td_log.val_SI = (
                (self.ttd_l.val_SI - self.ttd_u.val_SI)
                / math.log(self.ttd_l.val_SI / self.ttd_u.val_SI)
            )

        self.kA.val_SI = -self.Q.val_SI / self.td_log.val_SI
