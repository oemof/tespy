# -*- coding: utf-8

"""Module of class Turbine.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/turbomachinery/turbine.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import component_registry
from tespy.components.turbomachinery.base import Turbomachine
from tespy.tools import logger
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import isentropic
from tespy.tools.fluid_properties import single_fluid
from tespy.tools.helpers import _get_dependents


@component_registry
class Turbine(Turbomachine):
    r"""
    Class for gas or steam turbines.

    **Mandatory Equations**

    - fluid: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - mass flow: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.dp_structure_matrix`
    - :py:meth:`tespy.components.component.Component.pr_structure_matrix`
    - :py:meth:`tespy.components.turbomachinery.base.Turbomachine.energy_balance_func`
    - :py:meth:`tespy.components.turbomachinery.turbine.Turbine.eta_s_func`
    - :py:meth:`tespy.components.turbomachinery.turbine.Turbine.eta_s_char_func`
    - :py:meth:`tespy.components.turbomachinery.turbine.Turbine.cone_func`

    Inlets/Outlets

    - in1
    - out1

    Optional outlets

    - power

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

    pr : float, dict
        Outlet to inlet pressure ratio, :math:`pr/1`

    dp : float, dict
        Inlet to outlet pressure difference, :math:`dp/\text{p}_\text{unit}`
        Is specified in the Network's pressure unit

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
    >>> import os
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC", "enthalpy": "kJ/kg",
    ...     "mass_flow": "t/h"
    ... })
    >>> si = Sink('sink')
    >>> so = Source('source')
    >>> t = Turbine('turbine')
    >>> inc = Connection(so, 'out1', t, 'in1')
    >>> outg = Connection(t, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)

    In design conditions the isentropic efficiency is specified. For offdesign
    a characteristic function will be applied, together with Stodola's cone
    law coupling the turbine mass flow to inlet pressure.

    >>> t.set_attr(eta_s=0.9, design=['eta_s'],
    ... offdesign=['eta_s_char', 'cone'])
    >>> inc.set_attr(fluid={'water': 1}, m=36, T=550, p=110, design=['p'])
    >>> outg.set_attr(p=0.5)
    >>> nw.solve('design')
    >>> nw.save('tmp.json')
    >>> round(t.P.val, 0)
    -10452574.0
    >>> round(outg.x.val, 3)
    0.914
    >>> inc.set_attr(m=28.8)
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(t.eta_s.val, 3)
    0.898
    >>> round(inc.p.val, 1)
    88.6
    >>> os.remove('tmp.json')
    """

    @staticmethod
    def poweroutlets():
        return ["power"]

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        if len(self.power_outl) > 0:
            constraints["energy_connector_balance"] = dc_cmc(**{
                "func": self.energy_connector_balance_func,
                "dependents": self.energy_connector_dependents,
                "num_eq_sets": 1
            })

        return constraints

    def energy_connector_balance_func(self):
        r"""
        (optional) energy balance equation connecting the power connector to
        the component's power

        Returns
        -------
        residual : float
            Residual value of equation

            .. math::

                0=\dot E + \dot{m}_{in}\cdot\left(h_{out}-h_{in}\right)
        """
        return self.power_outl[0].E.val_SI + self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI
        )

    def energy_connector_dependents(self):
        return [self.power_outl[0].E, self.inl[0].m, self.outl[0].h, self.inl[0].h]

    def get_parameters(self):
        parameters = super().get_parameters()
        parameters["P"].max_val = 0
        parameters["pr"].max_val = 1
        parameters["pr"].min_val = 0
        parameters["dp"].min_val = 0
        parameters.update({
            'eta_s': dc_cp(
                min_val=0, max_val=1, num_eq_sets=1,
                func=self.eta_s_func,
                dependents=self.eta_s_dependents,
                deriv=self.eta_s_deriv,
                quantity="efficiency"
            ),
            'eta_s_char': dc_cc(
                param='m', num_eq_sets=1,
                func=self.eta_s_char_func,
                dependents=self.eta_s_char_dependents
            ),
            'cone': dc_simple(
                num_eq_sets=1,
                func=self.cone_func,
                dependents=self.cone_dependents
            )
        })
        return parameters

    def eta_s_func(self):
        r"""
        Equation for given isentropic efficiency of a turbine.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = -\left( h_{out} - h_{in} \right) +
                \left( h_{out,s} - h_{in} \right) \cdot \eta_{s,e}
        """
        inl = self.inl[0]
        outl = self.outl[0]
        return (
            -(outl.h.val_SI - inl.h.val_SI)
            + (
                isentropic(
                    inl.p.val_SI,
                    inl.h.val_SI,
                    outl.p.val_SI,
                    inl.fluid_data,
                    inl.mixing_rule,
                    T0=inl.T.val_SI
                )
                - inl.h.val_SI
            ) * self.eta_s.val_SI
        )

    def eta_s_deriv(self, increment_filter, k, dependents=None):
        r"""
        Partial derivatives for isentropic efficiency function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        dependents = dependents["scalars"][0]
        f = self.eta_s_func
        i = self.inl[0]
        o = self.outl[0]

        if o.h.is_var and not i.h.is_var:
            self._partial_derivative(o.h, k, -1, increment_filter)
            # remove o.h from the dependents
            dependents = dependents.difference(_get_dependents([o.h])[0])

        for dependent in dependents:
            self._partial_derivative(dependent, k, f, increment_filter)

    def eta_s_dependents(self):
        return [
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
        ]

    def cone_func(self):
        r"""
        Equation for stodolas cone law.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \frac{\dot{m}_{in,ref} \cdot p_{in}}{p_{in,ref}} \cdot
                \sqrt{\frac{p_{in,ref} \cdot v_{in}}{p_{in} \cdot v_{in,ref}}}
                \cdot \sqrt{\frac{1 - \left(\frac{p_{out}}{p_{in}} \right)^{2}}
                {1 - \left(\frac{p_{out,ref}}{p_{in,ref}} \right)^{2}}} -
                \dot{m}_{in}
        """
        n = 1
        i = self.inl[0]
        o = self.outl[0]
        vol = i.calc_vol(T0=i.T.val_SI)
        residual = (
            - i.m.val_SI + i.m.design * i.p.val_SI / i.p.design
            * (i.p.design * i.vol.design / (i.p.val_SI * vol)) ** 0.5
            * abs(
                    (1 - (o.p.val_SI / i.p.val_SI) ** ((n + 1) / n))
                    / (1 - (self.pr.design) ** ((n + 1) / n))
            ) ** 0.5
        )
        return residual

    def cone_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
        ]

    def eta_s_char_func(self):
        r"""
        Equation for given isentropic efficiency characteristic.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = - \left( h_\mathrm{out} - h_\mathrm{in} \right) +
                \eta_\mathrm{s,design} \cdot f\left( expr \right) \cdot
                \left(h_\mathrm{out,s}-h_\mathrm{in}\right)
        """
        p = self.eta_s_char.param
        expr = self.get_char_expr(p)
        if not expr:
            msg = (
                "Please choose a valid parameter, you want to link the "
                f"isentropic efficiency to at component {self.label}."
            )
            logger.error(msg)
            raise ValueError(msg)

        inl = self.inl[0]
        outl = self.outl[0]
        return (
            -(outl.h.val_SI - inl.h.val_SI)
            + self.eta_s.design * self.eta_s_char.char_func.evaluate(expr)
            * (
                isentropic(
                    inl.p.val_SI,
                    inl.h.val_SI,
                    outl.p.val_SI,
                    inl.fluid_data,
                    inl.mixing_rule,
                    T0=inl.T.val_SI
                ) - inl.h.val_SI
            )
        )

    def eta_s_char_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
        ]

    def calc_eta_s(self):
        inl = self.inl[0]
        outl = self.outl[0]
        return (
            (outl.h.val_SI - inl.h.val_SI)
            / (isentropic(
                    inl.p.val_SI,
                    inl.h.val_SI,
                    outl.p.val_SI,
                    inl.fluid_data,
                    inl.mixing_rule,
                    T0=inl.T.val_SI
                ) - inl.h.val_SI
            )
        )

    def convergence_check(self):
        r"""
        Perform a convergence check.

        Note
        ----
        Manipulate enthalpies/pressure at inlet and outlet if not specified by
        user to match physically feasible constraints.
        """
        i, o = self.inl[0], self.outl[0]

        if i.h.val_SI <= o.h.val_SI and o.h.is_var:
            o.h.set_reference_val_SI(i.h.val_SI - 100e3)

        if i.p.val_SI <= o.p.val_SI and o.p.is_var:
            o.p.set_reference_val_SI(i.p.val_SI * 2 /3)

    @staticmethod
    def initialise_source(c, key):
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
            fluid = single_fluid(c.fluid_data)
            if fluid is not None:
                return c.fluid.wrapper[fluid]._p_crit / 2
            else:
                return 1e5
        elif key == 'h':
            fluid = single_fluid(c.fluid_data)
            if fluid is not None:
                if c.p.val_SI >= c.fluid.wrapper[fluid]._p_crit:
                    temp = c.fluid.wrapper[fluid]._T_crit * 1.2
                    return h_mix_pT(c.p.val_SI, temp, c.fluid_data)
                else:
                    return h_mix_pQ(c.p.val_SI, 1, c.fluid_data, c.mixing_rule)
            else:
                temp = 1000
                return h_mix_pT(c.p.val_SI, temp, c.fluid_data, c.mixing_rule)

    @staticmethod
    def initialise_target(c, key):
        r"""
        Return a starting value for pressure and enthalpy at inlet.

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
            fluid = single_fluid(c.fluid_data)
            if fluid is not None:
                return c.fluid.wrapper[fluid]._p_crit / 4 * 3
            else:
                return 10e5
        elif key == 'h':
            fluid = single_fluid(c.fluid_data)
            if fluid is not None:
                if c.p.val_SI >= c.fluid.wrapper[fluid]._p_crit:
                    temp = c.fluid.wrapper[fluid]._T_crit * 1.4
                    return h_mix_pT(c.p.val_SI, temp, c.fluid_data)
                else:
                    return h_mix_pQ(c.p.val_SI, 1, c.fluid_data, c.mixing_rule) + 1e5
            else:
                temp = 500
                return h_mix_pT(c.p.val_SI, temp, c.fluid_data, c.mixing_rule)

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        super().calc_parameters()
        self.eta_s.val_SI = self.calc_eta_s()

    def exergy_balance(self, T0):
        r"""
        Calculate exergy balance of a turbine.

        Parameters
        ----------
        T0 : float
            Ambient temperature T0 / K.

        Note
        ----
        .. math::

            \dot{E}_\mathrm{P} =
            \begin{cases}
            -P & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
            -P + \dot{E}_\mathrm{out}^\mathrm{T}
            & T_\mathrm{in} > T_0 \geq T_\mathrm{out}\\
            -P +\dot{E}_\mathrm{out}^\mathrm{T}- \dot{E}_\mathrm{in}^\mathrm{T}
            & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
            \end{cases}

           \dot{E}_\mathrm{F} =
           \begin{cases}
           \dot{E}_\mathrm{in}^\mathrm{PH} - \dot{E}_\mathrm{out}^\mathrm{PH}
           & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
           \dot{E}_\mathrm{in}^\mathrm{T} + \dot{E}_\mathrm{in}^\mathrm{M} -
           \dot{E}_\mathrm{out}^\mathrm{M}
           & T_\mathrm{in} > T_0 \geq T_\mathrm{out}\\
           \dot{E}_\mathrm{in}^\mathrm{M} - \dot{E}_\mathrm{out}^\mathrm{M}
           & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
           \end{cases}

           \dot{E}_\mathrm{bus} = -P
        """
        if self.inl[0].T.val_SI >= T0 and self.outl[0].T.val_SI >= T0:
            self.E_P = -self.P.val
            self.E_F = self.inl[0].Ex_physical - self.outl[0].Ex_physical
        elif self.inl[0].T.val_SI > T0 and self.outl[0].T.val_SI <= T0:
            self.E_P = -self.P.val + self.outl[0].Ex_therm
            self.E_F = self.inl[0].Ex_therm + (
                self.inl[0].Ex_mech - self.outl[0].Ex_mech)
        elif self.inl[0].T.val_SI <= T0 and self.outl[0].T.val_SI <= T0:
            self.E_P = -self.P.val + (
                self.outl[0].Ex_therm - self.inl[0].Ex_therm)
            self.E_F = self.inl[0].Ex_mech - self.outl[0].Ex_mech
        else:
            msg = ('Exergy balance of a turbine, where outlet temperature is '
                   'larger than inlet temperature is not implmented.')
            logger.warning(msg)
            self.E_P = np.nan
            self.E_F = np.nan

        self.E_bus = {"chemical": 0, "physical": 0, "massless": -self.P.val}
        self.E_D = self.E_F - self.E_P
        self.epsilon = self._calc_epsilon()
