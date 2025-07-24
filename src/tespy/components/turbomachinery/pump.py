# -*- coding: utf-8

"""Module of class Pump.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/turbomachinery/pump.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import component_registry
from tespy.components.turbomachinery.base import Turbomachine
from tespy.tools import logger
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import isentropic
from tespy.tools.fluid_properties import single_fluid
from tespy.tools.helpers import _get_dependents


@component_registry
class Pump(Turbomachine):
    r"""
    Class for axial or radial pumps.

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_func`
    - :py:meth:`tespy.components.turbomachinery.base.Turbomachine.energy_balance_func`
    - :py:meth:`tespy.components.turbomachinery.pump.Pump.eta_s_func`
    - :py:meth:`tespy.components.turbomachinery.pump.Pump.eta_s_char_func`
    - :py:meth:`tespy.components.turbomachinery.pump.Pump.flow_char_func`

    Inlets/Outlets

    - in1
    - out1

    Image

    .. image:: /api/_images/Pump.svg
       :alt: flowsheet of the pump
       :align: center
       :class: only-light

    .. image:: /api/_images/Pump_darkmode.svg
       :alt: flowsheet of the pump
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

    pr : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio, :math:`pr/1`

    eta_s_char : tespy.tools.characteristics.CharLine, dict
        Characteristic curve for isentropic efficiency, provide CharLine as
        function :code:`func`.

    flow_char : tespy.tools.characteristics.CharLine, dict
        Characteristic curve for pressure rise as function of volumetric flow
        :math:`x/\frac{\text{m}^3}{\text{s}} \, y/\text{Pa}`.

    Example
    -------
    A pump with a known pump curve (difference pressure as function of
    volumetric flow) pumps 1,5 l/s of water in design conditions. E.g. for a
    given isentropic efficiency it is possible to calculate power consumption
    and pressure at the pump.

    >>> from tespy.components import Sink, Source, Pump
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools.characteristics import CharLine
    >>> import os
    >>> nw = Network(p_unit='bar', T_unit='C', h_unit='kJ / kg', v_unit='l / s',
    ... iterinfo=False)
    >>> si = Sink('sink')
    >>> so = Source('source')
    >>> pu = Pump('pump')
    >>> inc = Connection(so, 'out1', pu, 'in1')
    >>> outg = Connection(pu, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)

    After that we calculate offdesign performance using the pump curve and a
    characteristic function for the pump efficiency. We can calulate the
    offdesign efficiency and the volumetric flow, if the difference pressure
    changed. The default characteristc lines are to be found in the
    :py:mod:`tespy.data` module. Of course you are able to specify your own
    characteristcs, like done for the :code:`flow_char`. More information on
    how to specify characteristic functions are given in the corresponding part
    of the online documentation.

    >>> v = np.array([0, 0.4, 0.8, 1.2, 1.6, 2]) / 1000
    >>> dp = np.array([15, 14, 12, 9, 5, 0]) * 1e5
    >>> char = CharLine(x=v, y=dp)
    >>> pu.set_attr(eta_s=0.8, flow_char={'char_func': char, 'is_set': True},
    ... design=['eta_s'], offdesign=['eta_s_char'])
    >>> inc.set_attr(fluid={'water': 1}, p=1, T=20, v=1.5, design=['v'])
    >>> nw.solve('design')
    >>> nw.save('tmp.json')
    >>> round(pu.pr.val, 0)
    7.0
    >>> round(outg.p.val - inc.p.val, 0)
    6.0
    >>> round(pu.P.val, 0)
    1125.0
    >>> outg.set_attr(p=12)
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(pu.eta_s.val, 2)
    0.71
    >>> round(inc.v.val, 1)
    0.9
    >>> os.remove('tmp.json')
    """

    @staticmethod
    def powerinlets():
        return ["power"]

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
        parameters = super().get_parameters()
        parameters["P"].min_val = 0
        parameters["pr"].min_val = 1
        parameters["dp"].max_val = 0
        parameters.update({
            'eta_s': dc_cp(
                min_val=0, max_val=1, num_eq_sets=1,
                func=self.eta_s_func,
                dependents=self.eta_s_dependents,
                deriv=self.eta_s_deriv
            ),
            'eta_s_char': dc_cc(
                param='v', num_eq_sets=1,
                func=self.eta_s_char_func,
                dependents=self.eta_s_char_dependents
            ),
            'flow_char': dc_cc(
                param='v', num_eq_sets=1,
                func=self.flow_char_func,
                dependents=self.flow_char_dependents,
                char_params={'type': 'abs', 'inconn': 0, 'outconn': 0}
            )
        })
        return parameters

    def energy_connector_balance_func(self):
        r"""
        (optional) energy balance equation connecting the power connector to
        the component's power

        Returns
        -------
        residual : float
            Residual value of equation

            .. math::

                0=\dot E - \dot{m}_{in}\cdot\left(h_{out}-h_{in}\right)
        """
        return self.power_inl[0].E.val_SI - self.inl[0].m.val_SI * (
            self.outl[0].h.val_SI - self.inl[0].h.val_SI
        )

    def energy_connector_dependents(self):
        return [self.power_inl[0].E, self.inl[0].m, self.outl[0].h, self.inl[0].h]

    def eta_s_func(self):
        r"""
        Equation for given isentropic efficiency.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = -\left( h_{out} - h_{in} \right) \cdot \eta_{s} +
                \left( h_{out,s} - h_{in} \right)
        """
        i = self.inl[0]
        o = self.outl[0]
        return (
            (o.h.val_SI - i.h.val_SI) * self.eta_s.val - (
                isentropic(
                    i.p.val_SI,
                    i.h.val_SI,
                    o.p.val_SI,
                    i.fluid_data,
                    i.mixing_rule,
                    T0=None
                ) - self.inl[0].h.val_SI
            )
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
        f = self.eta_s_func

        if o.h.is_var and not i.h.is_var:
            self._partial_derivative(o.h, k, self.eta_s.val, increment_filter)
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

    def eta_s_char_func(self):
        r"""
        Equation for given isentropic efficiency characteristic.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \left(h_{out}-h_{in}\right) \cdot \eta_{s,design}
                \cdot f\left( expr \right) -\left( h_{out,s} - h_{in} \right)
        """
        p = self.eta_s_char.param
        expr = self.get_char_expr(p, **self.eta_s_char.char_params)
        if not expr:
            msg = ('Please choose a valid parameter, you want to link the '
                   'isentropic efficiency to at component ' + self.label + '.')
            logger.error(msg)
            raise ValueError(msg)

        i = self.inl[0]
        o = self.outl[0]
        return (
            (o.h.val_SI - i.h.val_SI)
            * self.eta_s.design * self.eta_s_char.char_func.evaluate(expr)
            - (
                isentropic(
                    i.p.val_SI,
                    i.h.val_SI,
                    o.p.val_SI,
                    i.fluid_data,
                    i.mixing_rule,
                    T0=None
                ) - i.h.val_SI
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

    def flow_char_func(self):
        r"""
        Equation for given flow characteristic of a pump.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = p_{out} - p_{in} - f\left( expr \right)
        """
        p = self.flow_char.param
        expr = self.get_char_expr(p, **self.flow_char.char_params)
        return (
            self.outl[0].p.val_SI - self.inl[0].p.val_SI
            - self.flow_char.char_func.evaluate(expr)
        )

    def flow_char_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
        ]

    def convergence_check(self):
        r"""
        Perform a convergence check.

        Note
        ----
        Manipulate enthalpies/pressure at inlet and outlet if not specified by
        user to match physically feasible constraints.
        """
        i = self.inl[0]
        o = self.outl[0]

        if o.h.is_var:
            fluid = single_fluid(i.fluid_data)
            if fluid is not None:
                if i.fluid_data[fluid]["wrapper"].back_end != "INCOMP":
                    h_max = h_mix_pQ(i.p.val_SI, 0, i.fluid_data)

                    if o.h.val_SI > h_max:
                        o.h.set_reference_val_SI(h_max - 20e3)

        if i.h.is_var and o.h.val_SI < i.h.val_SI:
            i.h.set_reference_val_SI(o.h.val_SI - 10e3)

        if o.p.is_var and o.p.val_SI < i.p.val_SI:
            o.p.set_reference_val_SI(i.p.val_SI * 1.5)

        if o.h.is_var and o.h.val_SI < i.h.val_SI:
            o.h.set_reference_val_SI(i.h.val_SI + 10e3)

        if i.p.is_var and o.p.val_SI < i.p.val_SI:
            i.p.set_reference_val_SI(o.p.val_SI * 2 / 3)
            i.p.val_SI = o.p.val_SI * 0.9

        if i.m.is_var and self.flow_char.is_set:
            vol = i.calc_vol(T0=i.T.val_SI)
            expr = i.m.val_SI * vol
            if expr > self.flow_char.char_func.x[-1]:
                i.m.set_reference_val_SI(self.flow_char.char_func.x[-1] / vol)
            elif expr < self.flow_char.char_func.x[1]:
                i.m.set_reference_val_SI(self.flow_char.char_func.x[0] / vol)

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
        fluid = single_fluid(c.fluid_data)
        if key == 'p':
            fluid = single_fluid(c.fluid_data)
            if fluid is not None:
                return c.fluid.wrapper[fluid]._p_crit / 2
            else:
                return 10e5
        elif key == 'h':
            fluid = single_fluid(c.fluid_data)
            if fluid is not None:
                temp = c.fluid.wrapper[fluid]._T_crit
                if temp is None:
                    temp = c.fluid.wrapper[fluid]._T_max

                dT = temp - c.fluid.wrapper[fluid]._T_min

                temp = temp - dT * 0.89
            else:
                # a pump with a mixture does not really make a lot of sense
                temp = 300
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
        fluid = single_fluid(c.fluid_data)
        if key == 'p':
            fluid = single_fluid(c.fluid_data)
            if fluid is not None:
                return c.fluid.wrapper[fluid]._p_crit / 2
            else:
                return 1e5
        elif key == 'h':
            fluid = single_fluid(c.fluid_data)
            if fluid is not None:
                temp = c.fluid.wrapper[fluid]._T_crit
                if temp is None:
                    temp = c.fluid.wrapper[fluid]._T_max

                dT = temp - c.fluid.wrapper[fluid]._T_min

                temp = temp - dT * 0.9
            else:
                # a pump with a mixture does not really make a lot of sense
                temp = 300
            return h_mix_pT(c.p.val_SI, temp, c.fluid_data, c.mixing_rule)

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        super().calc_parameters()

        i = self.inl[0]
        o = self.outl[0]
        self.eta_s.val =  (
            isentropic(
                i.p.val_SI,
                i.h.val_SI,
                o.p.val_SI,
                i.fluid_data,
                i.mixing_rule,
                T0=None
            ) - self.inl[0].h.val_SI
        ) / (o.h.val_SI - i.h.val_SI)

    def exergy_balance(self, T0):
        r"""
        Calculate exergy balance of a pump.

        Parameters
        ----------
        T0 : float
            Ambient temperature T0 / K.

        Note
        ----
        .. math::

            \dot{E}_\mathrm{P} =
            \begin{cases}
            \dot{E}_\mathrm{out}^\mathrm{PH} - \dot{E}_\mathrm{in}^\mathrm{PH}
            & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
            \dot{E}_\mathrm{out}^\mathrm{T} + \dot{E}_\mathrm{out}^\mathrm{M} -
            \dot{E}_\mathrm{in}^\mathrm{M}
            & T_\mathrm{out} > T_0 \leq T_\mathrm{in}\\
            \dot{E}_\mathrm{out}^\mathrm{M} - \dot{E}_\mathrm{in}^\mathrm{M}
            & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
            \end{cases}

            \dot{E}_\mathrm{F} =
            \begin{cases}
            P & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
            P + \dot{E}_\mathrm{in}^\mathrm{T}
            & T_\mathrm{out} > T_0 \leq T_\mathrm{in}\\
            P + \dot{E}_\mathrm{in}^\mathrm{T} -\dot{E}_\mathrm{out}^\mathrm{T}
            & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
            \end{cases}

            \dot{E}_\mathrm{bus} = P
        """
        if self.inl[0].T.val_SI >= T0 and self.outl[0].T.val_SI >= T0:
            self.E_P = self.outl[0].Ex_physical - self.inl[0].Ex_physical
            self.E_F = self.P.val
        elif self.inl[0].T.val_SI <= T0 and self.outl[0].T.val_SI > T0:
            self.E_P = self.outl[0].Ex_therm + (
                self.outl[0].Ex_mech - self.inl[0].Ex_mech)
            self.E_F = self.P.val + self.inl[0].Ex_therm
        elif self.inl[0].T.val_SI <= T0 and self.outl[0].T.val_SI <= T0:
            self.E_P = self.outl[0].Ex_mech - self.inl[0].Ex_mech
            self.E_F = self.P.val + (
                self.inl[0].Ex_therm - self.outl[0].Ex_therm)
        else:
            msg = ('Exergy balance of a pump, where outlet temperature is '
                   'smaller than inlet temperature is not implmented.')
            logger.warning(msg)
            self.E_P = np.nan
            self.E_F = np.nan

        self.E_bus = {
            "chemical": 0, "physical": 0, "massless": self.P.val
        }
        self.E_D = self.E_F - self.E_P
        self.epsilon = self._calc_epsilon()
