# -*- coding: utf-8

"""Module of class SolarCollector.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/solar_collector.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.heat_exchangers.simple import SimpleHeatExchanger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp


@component_registry
class SolarCollector(SimpleHeatExchanger):
    r"""
    The solar collector calculates heat output from irradiance.

    **Mandatory Equations**

    - fluid: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - mass flow: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_structure_matrix`
    - :py:meth:`tespy.components.component.Component.dp_structure_matrix`
    - :py:meth:`tespy.components.component.Component.zeta_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.darcy_func`
    - :py:meth:`tespy.components.heat_exchangers.simple.SimpleHeatExchanger.hazen_williams_func`
    - :py:meth:`tespy.components.heat_exchangers.solar_collector.SolarCollector.energy_group_func`

    Inlets/Outlets

    - in1
    - out1

    Image

    .. image:: /api/_images/SolarCollector.svg
       :alt: flowsheet of the solar collector
       :align: center
       :class: only-light

    .. image:: /api/_images/SolarCollector_darkmode.svg
       :alt: flowsheet of the solar collector
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

    Q : float, dict, :code:`"var"`
        Heat transfer, :math:`Q/\text{W}`.

    pr : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio, :math:`pr/1`.

    zeta : float, dict, :code:`"var"`
        Geometry independent friction coefficient,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    D : float, dict, :code:`"var"`
        Diameter of the pipes, :math:`D/\text{m}`.

    L : float, dict, :code:`"var"`
        Length of the pipes, :math:`L/\text{m}`.

    ks : float, dict, :code:`"var"`
        Pipe's roughness, :math:`ks/\text{m}`.

    darcy_group : str, dict
        Parametergroup for pressure drop calculation based on pipes dimensions
        using darcy weissbach equation.

    ks_HW : float, dict, :code:`"var"`
        Pipe's roughness, :math:`ks/\text{1}`.

    hw_group : str, dict
        Parametergroup for pressure drop calculation based on pipes dimensions
        using hazen williams equation.

    E : float, dict, :code:`"var"`
        irradiance at tilted collector surface area,
        :math:`E/\frac{\text{W}}{\text{m}^2}`.

    eta_opt : float, dict, :code:`"var"`
        optical loss at surface cover,
        :math:`\eta_{opt}`.

    lkf_lin : float, dict, :code:`"var"`
        Linear thermal loss key figure,
        :math:`\alpha_1/\frac{\text{W}}{\text{K} \cdot \text{m}^2}`.

    lkf_quad : float, dict, :code:`"var"`
        Quadratic thermal loss key figure,
        :math:`\alpha_2/\frac{\text{W}}{\text{K}^2 \cdot \text{m}^2}`.

    A : float, dict, :code:`"var"`
        Collector surface area :math:`A/\text{m}^2`.

    Tamb : float, dict
        Ambient temperature, provide parameter in network's temperature unit.

    energy_group : str, dict
        Parametergroup for energy balance of solarthermal collector.

    Example
    -------
    The solar collector is used to calculate heat transferred to the heating
    system from irradiance on a tilted plane. For instance, it is possible to
    calculate the collector surface area required to transfer a specific amount
    of heat at a given irradiance. The collector parameters are the linear and
    the quadratic loss keyfigure as well as the optical effifiency.

    >>> from tespy.components import Sink, Source, SolarCollector
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import os
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC", "enthalpy": "kJ/kg"
    ... })
    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> sc = SolarCollector('solar collector')
    >>> sc.set_attr(pr=0.95, Q=1e4, design=['pr', 'Q'], offdesign=['zeta'],
    ...     Tamb=25, A='var', eta_opt=0.92, lkf_lin=1, lkf_quad=0.005, E=8e2)
    >>> inc = Connection(so, 'out1', sc, 'in1')
    >>> outg = Connection(sc, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)

    The outlet temperature should be at 90 °C at a constant mass flow, which
    is determined in the design calculation. In offdesign operation (at a
    different irradiance) using the calculated surface area and mass flow, it
    is possible to predict the outlet temperature. It would instead be
    possible to calulate the change in mass flow required to hold the
    specified outlet temperature, too.

    >>> inc.set_attr(fluid={'H2O': 1}, T=40, p=3, offdesign=['m'])
    >>> outg.set_attr(T=90, design=['T'])
    >>> nw.solve('design')
    >>> nw.save('tmp.json')
    >>> round(sc.A.val, 1)
    14.5
    >>> sc.set_attr(A=sc.A.val, E=5e2, Tamb=20)
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(sc.Q.val, 1)
    6083.8
    >>> round(outg.T.val, 1)
    70.5
    >>> os.remove('tmp.json')
    """

    def get_parameters(self):
        data = super().get_parameters()
        for k in ["kA_group", "kA_char_group", "kA", "kA_char"]:
            del data[k]

        data.update({
            'E': dc_cp(min_val=0, quantity="heat"),
            'A': dc_cp(min_val=0, quantity="area"),
            'eta_opt': dc_cp(min_val=0, max_val=1, quantity="efficiency"),
            'lkf_lin': dc_cp(min_val=0),
            'lkf_quad': dc_cp(min_val=0),
            'Tamb': dc_cp(quantity="temperature"),
            'Q_loss': dc_cp(max_val=0, _val=0, quantity="temperature"),
            'energy_group': dc_gcp(
                elements=['E', 'eta_opt', 'lkf_lin', 'lkf_quad', 'A', 'Tamb'],
                num_eq_sets=1,
                func=self.energy_group_func,
                dependents=self.energy_group_dependents
            )
        })
        return data

    def energy_group_func(self):
        r"""
        Equation for solar collector energy balance.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \dot{m} \cdot \left( h_{out} - h_{in} \right)\\
                & - A \cdot \left[E \cdot \eta_{opt} - \alpha_1 \cdot
                \left(T_m - T_{amb} \right) - \alpha_2 \cdot
                \left(T_m - T_{amb}\right)^2 \right]\\
                T_m = & \frac{T_{out} + T_{in}}{2}\\
                \end{split}

            Reference: :cite:`Quaschning2013`.
        """
        i = self.inl[0]
        o = self.outl[0]

        T_m = 0.5 * (i.calc_T() + o.calc_T())

        return (
            i.m.val_SI * (o.h.val_SI - i.h.val_SI)
            - self.A.val_SI * (
                self.E.val_SI * self.eta_opt.val_SI
                - (T_m - self.Tamb.val_SI) * self.lkf_lin.val_SI
                - self.lkf_quad.val_SI * (T_m - self.Tamb.val_SI) ** 2
            )
        )

    def energy_group_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
        ] + [self.get_attr(element) for element in self.energy_group.elements]

    def convergence_check(self):
        pass

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        i = self.inl[0]
        o = self.outl[0]

        self.Q.val_SI = i.m.val_SI * (o.h.val_SI - i.h.val_SI)
        self.pr.val_SI = o.p.val_SI / i.p.val_SI
        self.dp.val_SI = i.p.val_SI - o.p.val_SI
        self.zeta.val_SI = self.calc_zeta(i, o)

        if self.energy_group.is_set:
            self.Q_loss.val_SI = -(self.E.val_SI * self.A.val_SI - self.Q.val_SI)
            self.Q_loss.is_result = True
        else:
            self.Q_loss.is_result = False
