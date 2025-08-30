# -*- coding: utf-8

"""Module of class Compressor.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/turbomachinery/compressor.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import component_registry
from tespy.components.turbomachinery.base import Turbomachine
from tespy.tools import logger
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import isentropic
from tespy.tools.fluid_properties import single_fluid
from tespy.tools.helpers import _get_dependents


@component_registry
class Compressor(Turbomachine):
    r"""
    Class for axial or radial compressor.

    **Mandatory Equations**

    - fluid: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - mass flow: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.dp_structure_matrix`
    - :py:meth:`tespy.components.component.Component.pr_structure_matrix`
    - :py:meth:`tespy.components.turbomachinery.base.Turbomachine.energy_balance_func`
    - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.eta_s_func`
    - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.eta_s_char_func`
    - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.char_map_eta_s_func`
    - :py:meth:`tespy.components.turbomachinery.compressor.Compressor.char_map_pr_func`

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

    char_map : tespy.tools.characteristics.CharMap, dict
        Characteristic map for pressure rise and isentropic efficiency vs.
        nondimensional mass flow, see
        :py:class:`tespy.tools.characteristics.CharMap` for further
        information. Provide a CompressorMap as function :code:`func`.

    igva : float, dict, :code:`"var"`
        Inlet guide vane angle, :math:`igva/^\circ`.

    Example
    -------
    Create an air compressor model and calculate the power required for
    compression of 50 l/s of ambient air to 5 bars. Using a generic compressor
    map how does the efficiency change in different operation mode (e.g. 90 %
    of nominal volumetric flow)?

    >>> from tespy.components import Sink, Source, Compressor
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import os
    >>> nw = Network(p_unit='bar', T_unit='C', h_unit='kJ / kg', v_unit='l / s',
    ... iterinfo=False)
    >>> si = Sink('sink')
    >>> so = Source('source')
    >>> comp = Compressor('compressor')
    >>> inc = Connection(so, 'out1', comp, 'in1')
    >>> outg = Connection(comp, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)

    Specify the compressor parameters: nominal efficiency and pressure ratio.
    For offdesign mode the characteristic map is selected instead of the
    isentropic efficiency. For offdesign, the inlet guide vane angle should be
    variable in order to maintain the same pressure ratio at a different
    volumetric flow.

    >>> comp.set_attr(pr=5, eta_s=0.8, design=['eta_s'],
    ... offdesign=['char_map_pr', 'char_map_eta_s'])
    >>> inc.set_attr(fluid={'air': 1}, p=1, T=20, v=50)
    >>> nw.solve('design')
    >>> nw.save('tmp.json')
    >>> round(comp.P.val, 0)
    12772.0
    >>> round(comp.eta_s.val, 2)
    0.8
    >>> inc.set_attr(v=45)
    >>> comp.set_attr(igva='var')
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(comp.eta_s.val, 2)
    0.77
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
                deriv=self.eta_s_deriv,
                dependents=self.eta_s_dependents,
                quantity="efficiency"
            ),
            'eta_s_char': dc_cc(
                param='m', num_eq_sets=1,
                func=self.eta_s_char_func,
                dependents=self.eta_s_char_dependents,
            ),
            'igva': dc_cp(min_val=-90, max_val=90, d=1e-4, val=0),
            'char_map_eta_s': dc_cm(),
            'char_map_eta_s_group': dc_gcp(
                elements=['char_map_eta_s', 'igva'], num_eq_sets=1,
                func=self.char_map_eta_s_func,
                dependents=self.char_map_dependents
            ),
            'char_map_pr': dc_cm(),
            'char_map_pr_group': dc_gcp(
                elements=['char_map_pr', 'igva'],
                num_eq_sets=1,
                func=self.char_map_pr_func,
                dependents=self.char_map_dependents
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
        Equation for given isentropic efficiency of a compressor.

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
            (o.h.val_SI - i.h.val_SI) * self.eta_s.val_SI - (
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
            self._partial_derivative(o.h, k, self.eta_s.val_SI, increment_filter)
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

    def char_map_pr_func(self):

        r"""
        Calculate pressure ratio from characteristic map.

        Returns
        -------
        residual : float
            Residual value of equations.

        Note
        ----
        - X: speedline index (rotational speed is constant)
        - Y: nondimensional mass flow
        - igva: variable inlet guide vane angle for value manipulation
          according to :cite:`GasTurb2018`.

        .. math::

            X = \sqrt{\frac{T_\mathrm{in,design}}{T_\mathrm{in}}}\\
            Y = \frac{\dot{m}_\mathrm{in} \cdot p_\mathrm{in,design}}
            {\dot{m}_\mathrm{in,design} \cdot p_\mathrm{in} \cdot X}\\
            \vec{Y} = f\left(X,Y\right)\cdot\left(1-\frac{igva}{100}\right)\\
            \vec{Z} = f\left(X,Y\right)\cdot\left(1-\frac{igva}{100}\right)\\
            0 = \frac{p_{out} \cdot p_{in,design}}
            {p_\mathrm{in} \cdot p_\mathrm{out,design}}-
            f\left(Y,\vec{Y},\vec{Z}\right)
        """
        i = self.inl[0]
        o = self.outl[0]

        beta = np.sqrt(i.T.design / i.calc_T())
        y = (i.m.val_SI * i.p.design) / (i.m.design * i.p.val_SI * beta)

        yarr, zarr = self.char_map_pr.char_func.evaluate_x(beta)
        # value manipulation with igva
        yarr *= (1 - self.igva.val_SI / 100)
        zarr *= (1 - self.igva.val_SI / 100)
        pr = self.char_map_pr.char_func.evaluate_y(y, yarr, zarr)

        return (o.p.val_SI / i.p.val_SI) - pr * self.pr.design

    def char_map_eta_s_func(self):
        r"""
        Calculate isentropic efficiency from characteristic map.

        Returns
        -------
        residual : float
            Residual value of equation.

        Note
        ----
        - X: speedline index (rotational speed is constant)
        - Y: nondimensional mass flow
        - igva: variable inlet guide vane angle for value manipulation
          according to :cite:`GasTurb2018`.

        .. math::

            X = \sqrt{\frac{T_\mathrm{in,design}}{T_\mathrm{in}}}\\
            Y = \frac{\dot{m}_\mathrm{in} \cdot p_\mathrm{in,design}}
            {\dot{m}_\mathrm{in,design} \cdot p_\mathrm{in} \cdot X}\\
            \vec{Y} = f\left(X,Y\right)\cdot\left(1-\frac{igva}{100}\right)\\
            \vec{Z}=f\left(X,Y\right)\cdot\left(1-\frac{igva^2}{10000}\right)\\
            0 = \frac{\eta_\mathrm{s}}{\eta_\mathrm{s,design}} -
            f\left(Y,\vec{Y},\vec{Z}\right)
        """
        i = self.inl[0]
        o = self.outl[0]

        x = np.sqrt(i.T.design / i.calc_T())
        y = (i.m.val_SI * i.p.design) / (i.m.design * i.p.val_SI * x)

        yarr, zarr = self.char_map_eta_s.char_func.evaluate_x(x)
        # value manipulation with igva
        yarr *= (1 - self.igva.val_SI / 100)
        zarr *= (1 - self.igva.val_SI ** 2 / 10000)
        eta = self.char_map_eta_s.char_func.evaluate_y(y, yarr, zarr)

        return (
            (
            isentropic(
                i.p.val_SI,
                i.h.val_SI,
                o.p.val_SI,
                i.fluid_data,
                i.mixing_rule,
                T0=i.T.val_SI
            ) - i.h.val_SI)
            / (o.h.val_SI - i.h.val_SI) - eta * self.eta_s.design
        )

    def char_map_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].p,
            self.outl[0].h,
            self.igva
        ]

    def convergence_check(self):
        r"""
        Perform a convergence check.

        Note
        ----
        Manipulate enthalpies/pressure at inlet and outlet if not specified by
        user to match physically feasible constraints.
        """
        i, o = self.inl[0], self.outl[0]

        if o.p.is_var and o.p.val_SI < i.p.val_SI:
            o.p.set_reference_val_SI(i.p.val_SI * 1.5)

        if o.h.is_var and o.h.val_SI < i.h.val_SI:
            o.h.set_reference_val_SI(i.h.val_SI + 100e3)

        if i.p.is_var and o.p.val_SI < i.p.val_SI:
            i.p.set_reference_val_SI(o.p.val_SI * 2 / 3)
            i.p.val_SI = o.p.val_SI * 0.9

        if i.h.is_var and o.h.val_SI < i.h.val_SI:
            i.h.set_reference_val_SI(o.h.val_SI - 100e3)

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
                return 10e5
        elif key == 'h':
            fluid = single_fluid(c.fluid_data)
            if fluid is not None:
                if c.p.val_SI < c.fluid.wrapper[fluid]._p_crit:
                    return h_mix_pQ(c.p.val_SI, 1, c.fluid_data, c.mixing_rule) + 2e5
                else:
                    temp = c.fluid.wrapper[fluid]._T_crit
                    return h_mix_pT(
                        c.p.val_SI, temp * 1.2, c.fluid_data, c.mixing_rule
                    ) + 2e5
            else:
                temp = 450
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
                return c.fluid.wrapper[fluid]._p_crit / 3
            else:
                return 1e5
        elif key == 'h':
            fluid = single_fluid(c.fluid_data)
            if fluid is not None:
                return h_mix_pQ(c.p.val_SI, 1, c.fluid_data, c.mixing_rule)
            else:
                temp = 350
                return h_mix_pT(c.p.val_SI, temp, c.fluid_data, c.mixing_rule)

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        super().calc_parameters()

        i = self.inl[0]
        o = self.outl[0]
        self.eta_s.val_SI =  (
            isentropic(
                i.p.val_SI,
                i.h.val_SI,
                o.p.val_SI,
                i.fluid_data,
                i.mixing_rule,
                T0=None
            ) - self.inl[0].h.val_SI
        ) / (o.h.val_SI - i.h.val_SI)

    def check_parameter_bounds(self):
        r"""Check parameter value limits."""
        _no_limit_violations = super().check_parameter_bounds()

        for data in [self.char_map_pr, self.char_map_eta_s]:
            if data.is_set:
                x = np.sqrt(self.inl[0].T.design / self.inl[0].T.val_SI)
                y = (self.inl[0].m.val_SI * self.inl[0].p.design) / (
                    self.inl[0].m.design * self.inl[0].p.val_SI * x)
                yarr = data.char_func.get_domain_errors_x(x, self.label)
                yarr *= (1 - self.igva.val_SI / 100)
                data.char_func.get_domain_errors_y(y, yarr, self.label)

        return _no_limit_violations

    def exergy_balance(self, T0):
        r"""
        Calculate exergy balance of a compressor.

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
            msg = ('Exergy balance of a compressor, where outlet temperature '
                   'is smaller than inlet temperature is not implmented.')
            logger.warning(msg)
            self.E_P = np.nan
            self.E_F = np.nan

        self.E_bus = {
            "chemical": 0, "physical": 0, "massless": self.P.val
        }
        self.E_D = self.E_F - self.E_P
        self.epsilon = self._calc_epsilon()
