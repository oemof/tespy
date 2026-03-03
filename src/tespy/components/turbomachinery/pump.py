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
from tespy.tools.data_containers import ComponentCharacteristicMaps as dc_cm
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import isentropic
from tespy.tools.fluid_properties import single_fluid
from tespy.tools.global_vars import GRAVITY
from tespy.tools.helpers import _get_dependents


@component_registry
class Pump(Turbomachine):
    r"""
    Class for axial or radial pumps.

    **Mandatory Equations**

    - fluid: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - mass flow: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.pr_structure_matrix`
    - :py:meth:`tespy.components.component.Component.dp_structure_matrix`
    - :py:meth:`tespy.components.turbomachinery.base.Turbomachine.energy_balance_func`
    - :py:meth:`tespy.components.turbomachinery.pump.Pump.eta_s_func`
    - :py:meth:`tespy.components.turbomachinery.pump.Pump.eta_s_char_func`
    - :py:meth:`tespy.components.turbomachinery.pump.Pump.flow_char_func`

    Inlets/Outlets

    - in1
    - out1

    Optional inlets

    - power

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

    pr : float, dict
        Outlet to inlet pressure ratio, :math:`pr/1`

    dp : float, dict
        Inlet to outlet pressure difference, :math:`dp/\text{p}_\text{unit}`
        Is specified in the Network's pressure unit

    eta_s_char : tespy.tools.characteristics.CharLine, dict
        Characteristic curve for isentropic efficiency, provide CharLine as
        function :code:`func`.

    flow_char : tespy.tools.characteristics.CharLine, dict
        Characteristic curve for pressure rise as function of volumetric flow
        :math:`x/\frac{\text{m}^3}{\text{s}} \, y/\text{Pa}`.

    Example
    -------
    Below you will find two separate examples. First, a simple pump using a
    pressure difference to volumetric flow characteristic in both design and
    offdesign mode. Second, a pump implementing a pump characteristic map for
    hydraulic head over flow at variable frequencies and a characteristic map
    for the efficiency over flow at variable frequencies.

    A pump with a known pump curve (difference pressure as function of
    volumetric flow) pumps 1,5 l/s of water in design conditions. E.g. for a
    given isentropic efficiency it is possible to calculate power consumption
    and pressure at the pump.

    >>> from tespy.components import Sink, Source, Pump
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> from tespy.tools.characteristics import CharLine
    >>> import os
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC", "volumetric_flow": "l/s", "enthalpy": "kJ/kg"
    ... })
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
    :ref:`tespy.data <data_label>` module. Of course you are able to
    specify your own characteristcs, like done for the :code:`flow_char`. More
    information on how to specify characteristic functions are given in the
    corresponding part of the online documentation.

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

    In the second example we model a pump with characteristic maps. These can
    be retrieved from manufacturers, e.g., this one: :cite:`grundfos`. To keep
    the example simple, the pump maps here are more basic. To make use of the
    maps, we need the following information:

    - A numpy array of frequency values
    - Per frequency value a numpy array for volumetric flow
    - Per frequency value a corresponding numpy array for hydraulic head and
      for efficiency

    We can use the following data in the example:

    >>> import numpy as np
    >>> from tespy.tools import CharMap
    >>> frequencies = np.array([1500, 1750, 2000, 2250]) / 60  # SI units!
    >>> flows = np.array([
    ...     [0.   , 0.715, 1.431],  # corresponds to 1500
    ...     [0.   , 1.141, 2.284],  # corresponds to 1750...
    ...     [0.   , 1.568, 3.137],
    ...     [0.   , 1.995, 3.772]
    ... ])

    .. attention::

        The values need to be in SI units!

    Then we need to have the respective hydraulic head and efficiency for each
    of those values.

    >>> head = np.array([
    ...     [ 3.229,  3.023,  1.786],
    ...     [ 8.239,  7.712,  4.559],
    ...     [15.554, 14.558,  8.602],
    ...     [20.   , 20.   , 12.431]
    ... ])

    >>> efficiency = np.array([
    ...     [0. , 0.45, 0.4 ],
    ...     [0. , 0.5 , 0.45],
    ...     [0. , 0.52, 0.5 ],
    ...     [0. , 0.5 , 0.45]
    ... ])

    We can create the maps next:

    >>> pump_H_map = CharMap(
    ...     x=frequencies,
    ...     y=flows,
    ...     z=head
    ... )
    >>> pump_eta_map = CharMap(
    ...     x=frequencies,
    ...     y=flows,
    ...     z=efficiency
    ... )

    Now we create the model and impose the maps

    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar",
    ...     "temperature": "degC",
    ...     "volumetric_flow": "m3/s",
    ...     "enthalpy": "kJ/kg",
    ...     "frequency": "1/min"
    ... })
    >>> si = Sink('sink')
    >>> so = Source('source')
    >>> pu = Pump('pump')
    >>> inc = Connection(so, 'out1', pu, 'in1')
    >>> outg = Connection(pu, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> pu.set_attr(
    ...     head_flow_map={'char_func': pump_H_map, 'is_set': True},
    ...     eta_flow_map={'char_func': pump_eta_map, 'is_set': True},
    ...     frequency=2000
    ... )

    We can impose the inlet volumetric flow, for example, an actual set point
    on the map.

    >>> inc.set_attr(fluid={'water': 1}, p=1, T=20, v=1.568)
    >>> nw.solve('design')
    >>> round(pu.head.val, 3)
    14.558

    The process can of course be inverted, e.g., find the volumetric flow
    corresponding to a specific pressure at the outlet.

    >>> inc.set_attr(v=None)
    >>> outg.set_attr(p=2.2)
    >>> nw.solve("design")
    >>> round(inc.v.val, 3)
    2.174

    .. attention::

        Be aware that if the pump runs outside of the characteristic map it can
        easily happen that derivatives become zero because the interpolation
        returns the same values independent of the variables. A good initial
        guess can (but does not necessarily) help.

    >>> outg.set_attr(p=1.2)
    >>> nw.solve("design")
    >>> nw.status == 3  # check if status is linear dependency
    True

    Finally, it is possible to make the frequency variable to match a certain
    volumetric flow and pressure increase:

    >>> inc.set_attr(v=1.5)
    >>> outg.set_attr(p=2)
    >>> pu.set_attr(frequency="var")
    >>> nw.solve("design")
    >>> nw.status == 0  # check if converged
    True
    >>> round(pu.frequency.val)
    1862

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

    def _preprocess(self, row_idx):
        if self.frequency.is_var and self.head_flow_map.is_set:
            min_frequency = min(self.head_flow_map.char_func.x)
            max_frequency = max(self.head_flow_map.char_func.x)
            self.frequency.min_val = min_frequency
            self.frequency.max_val = max_frequency

        return super()._preprocess(row_idx)

    def get_parameters(self):
        parameters = super().get_parameters()
        parameters["P"].min_val = 0
        parameters["pr"].min_val = 1
        parameters["dp"].max_val = 0
        parameters.update({
            'eta': dc_cp(
                min_val=0, max_val=1, is_result=True,
                quantity="efficiency",
                description=(
                    "efficiency defined as specific incompressible flow work "
                    "over increase of enthalpy"
                )
            ),
            'eta_s': dc_cp(
                min_val=0, max_val=1, num_eq_sets=1,
                func=self.eta_s_func,
                dependents=self.eta_s_dependents,
                deriv=self.eta_s_deriv,
                quantity="efficiency",
                description="isentropic efficiency"
            ),
            'eta_s_char': dc_cc(
                param='v', num_eq_sets=1,
                func=self.eta_s_char_func,
                dependents=self.eta_s_char_dependents,
                description="isentropic efficiency lookup table for offdesign"
            ),
            'flow_char': dc_cc(
                param='v', num_eq_sets=1,
                func=self.flow_char_func,
                dependents=self.flow_char_dependents,
                char_params={'type': 'abs', 'inconn': 0, 'outconn': 0},
                description="pressure rise over volumetric flow lookup table"
            ),
            "frequency": dc_cp(
                min_val=0, max_val=10000, _potential_var=True,
                quantity="frequency"
            ),
            "eta_flow_map": dc_cm(
                description=(
                    "2D lookup table for pump efficiency over volumetric "
                    "flow and frequency"
                )
            ),
            "eta_flow_group": dc_gcp(
                elements=["eta_flow_map", "frequency"],
                num_eq_sets=1,
                func=self.eta_flow_frequency_group_func,
                dependents=self.eta_flow_frequency_group_dependents,
                description=(
                    "map function for efficiency over volumetric flow and "
                    "frequency"
                )
            ),
            "head_flow_map": dc_cm(
                description=(
                    "2D lookup table for hydraulic head over volumetric "
                    "flow and frequency"
                )
            ),
            "head_flow_map_group": dc_gcp(
                elements=["head_flow_map", "frequency"],
                num_eq_sets=1,
                func=self.head_flow_frequency_group_func,
                dependents=self.head_flow_frequency_group_dependents,
                description=(
                    "map function for efficiency over volumetric flow and "
                    "frequency"
                )
            ),
            "head": dc_cp(
                min_val=0, max_val=10000, num_eq_sets=1,
                func=self.hydraulic_head_func,
                dependents=self.hydraulic_head_dependents,
                quantity="length"
            )
        })
        return parameters

    def energy_connector_balance_func(self):
        r"""
        (optional) energy balance equation connecting the power connector to
        the component's power

        Returns
        -------
        float
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
        float
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
                    T0=i.T.val_SI,
                    T0_out=o.T.val_SI
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
        float
            Residual value of equation.

            .. math::

                0 = \left(h_{out}-h_{in}\right) \cdot \eta_{s,design}
                \cdot f\left( expr \right) -\left( h_{out,s} - h_{in} \right)
        """
        p = self.eta_s_char.param
        expr = self.get_char_expr(p, **self.eta_s_char.char_params)
        if not expr:
            msg = (
                "Please choose a valid parameter, you want to link the "
                f"isentropic efficiency to at component {self.label}."
            )
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
                    T0=i.T.val_SI,
                    T0_out=o.T.val_SI
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
        float
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

    def eta_flow_frequency_group_func(self):
        r"""
        Equation for efficiency :math:`\eta` over flow and frequency
        :math:`\omega` map

        Returns
        -------
        float
            Residual value of equation.

            .. math::

                0 = \eta \left(\omega,\dot{m}_\text{in}\cdot v_\text{in} \right)
                \cdot (h_\text{out} - h_\text{in})
                - v_\text{in}\cdot\left(p_\text{out} - p_\text{in}\right)
        """
        i = self.inl[0]
        o = self.outl[0]
        vol = i.calc_vol()
        flow = i.m.val_SI * vol
        frequency = self.frequency.val_SI
        eta = self.eta_flow_map.char_func.evaluate(frequency, flow)
        return (
            eta * (o.h.val_SI - i.h.val_SI) - vol * (o.p.val_SI - i.p.val_SI)
        )

    def eta_flow_frequency_group_dependents(self):
        i = self.inl[0]
        o = self.outl[0]
        return [i.m, i.p, i.h, o.p, o.h, self.frequency]

    def head_flow_frequency_group_func(self):
        r"""
        Equation for hydraulic head :math:`H` over flow and frequency
        :math:`\omega` map

        Returns
        -------
        float
            Residual value of equation.

            .. math::

                0 = H\left(\omega,\dot{m}_\text{in}\right)
                - \frac{\left( p_\text{out} - p_\text{in} \right) \cdot v}{g}
        """
        i = self.inl[0]
        flow = i.m.val_SI * i.calc_vol()
        frequency = self.frequency.val_SI
        head = self.head_flow_map.char_func.evaluate(frequency, flow)
        return head - self.calc_hydraulic_head()

    def head_flow_frequency_group_dependents(self):
        i = self.inl[0]
        o = self.outl[0]
        return [i.m, i.p, i.h, o.p, self.frequency]

    def hydraulic_head_func(self):
        r"""
        Equation for given hydraulic head :math:`H`

        Returns
        -------
        float
            Residual value of equation.

            .. math::

                0 = \frac{H \cdot g}{v}
                - \left( p_\text{out} - p_\text{in} \right)
        """
        return self.head.val_SI - self.calc_hydraulic_head()

    def hydraulic_head_dependents(self):
        i = self.inl[0]
        o = self.outl[0]
        return [i.p, i.h, o.p]

    def calc_hydraulic_head(self):
        r"""
        Calculate the hydraulic head :math:`H`

        Returns
        -------
        float
            Residual value of equation.

            .. math::

                H = \frac{\left( p_\text{out} - p_\text{in} \right) \cdot v}{g}
        """
        i = self.inl[0]
        o = self.outl[0]
        return (
            (o.p.val_SI - i.p.val_SI) * i.calc_vol() / GRAVITY
        )

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
                if i.p.val_SI < i.fluid_data[fluid]["wrapper"]._p_crit:
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
        self.eta_s.val_SI =  (
            isentropic(
                i.p.val_SI,
                i.h.val_SI,
                o.p.val_SI,
                i.fluid_data,
                i.mixing_rule,
                T0=i.T.val_SI,
                T0_out=o.T.val_SI
            ) - self.inl[0].h.val_SI
        ) / (o.h.val_SI - i.h.val_SI)

        self.head.val_SI = self.calc_hydraulic_head()
        self.eta.val_SI = (
            i.vol.val_SI * (o.p.val_SI - i.p.val_SI)
            / (o.h.val_SI - i.h.val_SI)
        )

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
