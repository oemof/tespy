# -*- coding: utf-8

"""Module of class Drum.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/nodes/drum.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import component_registry
from tespy.components.nodes.droplet_separator import DropletSeparator
from tespy.tools.fluid_properties import h_mix_pQ


@component_registry
class Drum(DropletSeparator):
    r"""
    A drum separates saturated gas from saturated liquid.

    **Mandatory Equations**

    - :py:meth:`tespy.components.nodes.base.NodeBase.mass_flow_func`
    - :py:meth:`tespy.components.nodes.base.NodeBase.pressure_equality_func`
    - :py:meth:`tespy.components.nodes.droplet_separator.DropletSeparator.fluid_func`
    - :py:meth:`tespy.components.nodes.droplet_separator.DropletSeparator.energy_balance_func`
    - :py:meth:`tespy.components.nodes.droplet_separator.DropletSeparator.outlet_states_func`

    Inlets/Outlets

    - in1, in2 (index 1: from economiser, index 2: from evaporator)
    - out1, out2 (index 1: saturated liquid, index 2: saturated gas)

    Image

    .. image:: /api/_images/Drum.svg
       :alt: flowsheet of the drum
       :align: center
       :class: only-light

    .. image:: /api/_images/Drum_darkmode.svg
       :alt: flowsheet of the drum
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

    Note
    ----
    If you are using a drum in a network with multiple fluids, it is likely
    the fluid propagation causes trouble. If this is the case, try to
    specify the fluid composition at another connection of your network.

    This component assumes, that the fluid composition between outlet 1 and
    inlet 2 does not change, thus there is no equation for the fluid mass
    fraction at the inlet 2!

    Example
    -------
    The drum separates saturated gas from saturated liquid. The liquid phase is
    transported to an evaporator, the staturated gas phase is extracted from
    the drum. In this example ammonia is evaporated using ambient air. A
    characteristic function is applied for the heat transfer coefficient of the
    evaporator. It is possible to load the CharLine with the function
    :code:`load_default_char` from the default lines. We want to use the
    'EVAPORATING FLUID' lines of the heat exchanger.

    >>> from tespy.components import Sink, Source, Drum, Pump, HeatExchanger
    >>> from tespy.connections import Connection, Ref
    >>> from tespy.networks import Network
    >>> from tespy.tools.characteristics import CharLine
    >>> from tespy.tools.characteristics import load_default_char as ldc
    >>> import shutil
    >>> nw = Network(T_unit='C', p_unit='bar', h_unit='kJ / kg', iterinfo=False)
    >>> fa = Source('feed ammonia')
    >>> amb_in = Source('air inlet')
    >>> amb_out = Sink('air outlet')
    >>> s = Sink('steam')
    >>> dr = Drum('drum')
    >>> dr.component()
    'drum'
    >>> ev = HeatExchanger('evaporator')
    >>> erp = Pump('evaporator reciculation pump')
    >>> f_dr = Connection(fa, 'out1', dr, 'in1')
    >>> dr_erp = Connection(dr, 'out1', erp, 'in1')
    >>> erp_ev = Connection(erp, 'out1', ev, 'in2')
    >>> ev_dr = Connection(ev, 'out2', dr, 'in2')
    >>> dr_s = Connection(dr, 'out2', s, 'in1')
    >>> nw.add_conns(f_dr, dr_erp, erp_ev, ev_dr, dr_s)
    >>> amb_ev = Connection(amb_in, 'out1', ev, 'in1')
    >>> ev_amb = Connection(ev, 'out1', amb_out, 'in1')
    >>> nw.add_conns(amb_ev, ev_amb)

    The ambient air enters the evaporator at 30 °C. The pinch point temperature
    difference (ttd_l) of the evaporator is at 5 K, and 1 MW of heat should be
    transferred. State of ammonia at the inlet is at -5 °C and 5 bar. From this
    design it is possible to calculate offdesign performance at 75 % part load.

    >>> char1 = ldc('heat exchanger', 'kA_char1', 'DEFAULT',
    ... CharLine)
    >>> char2 = ldc('heat exchanger', 'kA_char2', 'EVAPORATING FLUID',
    ... CharLine)
    >>> ev.set_attr(pr1=0.999, pr2=0.99, ttd_l=5, kA_char1=char1,
    ...     kA_char2=char2, design=['pr1', 'ttd_l'],
    ...     offdesign=['zeta1', 'kA_char'])
    >>> ev.set_attr(Q=-1e6)
    >>> erp.set_attr(eta_s=0.8)
    >>> f_dr.set_attr(p=5, T=-5)
    >>> erp_ev.set_attr(m=Ref(f_dr, 4, 0), fluid={'NH3': 1})
    >>> amb_ev.set_attr(fluid={'air': 1}, T=30)
    >>> ev_amb.set_attr(p=1)
    >>> nw.solve('design')
    >>> nw.save('tmp.json')
    >>> round(ev_amb.T.val - erp_ev.T.val ,1)
    5.0
    >>> round(f_dr.h.val, 1)
    322.7
    >>> round(dr_erp.h.val, 1)
    364.9
    >>> round(ev_dr.h.val, 1)
    687.2
    >>> round(f_dr.m.val, 2)
    0.78
    >>> ev.set_attr(Q=-0.75e6)
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(f_dr.m.val, 2)
    0.58
    >>> round(ev_amb.T.val - erp_ev.T.val ,1)
    3.0
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'drum'

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1', 'out2']

    def get_mandatory_constraints(self):
        num_mass_eq = 1
        if self.inl[1].m == self.outl[0].m:
            num_mass_eq = 0
        return {
            'mass_flow_constraints': {
                'func': self.mass_flow_func, 'deriv': self.mass_flow_deriv,
                'constant_deriv': True, 'latex': self.mass_flow_func_doc,
                'num_eq': num_mass_eq},
            'energy_balance_constraints': {
                'func': self.energy_balance_func,
                'deriv': self.energy_balance_deriv,
                'constant_deriv': False, 'latex': self.energy_balance_func_doc,
                'num_eq': 1},
            'pressure_constraints': {
                'func': self.pressure_equality_func,
                'deriv': self.pressure_equality_deriv,
                'constant_deriv': True,
                'latex': self.pressure_equality_func_doc,
                'num_eq': self.num_i + self.num_o - 1},
            'outlet_constraints': {
                'func': self.outlet_states_func,
                'deriv': self.outlet_states_deriv,
                'constant_deriv': False,
                'latex': self.outlet_states_func_doc,
                'num_eq': 2}
        }

    def preprocess(self, num_nw_vars):
        super().preprocess(num_nw_vars)
        self._propagation_start = False

    def mass_flow_func(self):
        r"""
        Calculate the residual value for mass flow balance equation.

        Returns
        -------
        res : float
            Residual value of equation.

            .. math::

                0 = \sum \dot{m}_{in,i} - \sum \dot{m}_{out,j} \;
                \forall i \in inlets, \forall j \in outlets
        """
        if self.inl[1].m == self.outl[0].m:
            return self.inl[0].m.val_SI - self.outl[1].m.val_SI
        else:
            res = 0
            for i in self.inl:
                res += i.m.val_SI
            for o in self.outl:
                res -= o.m.val_SI

            return res

    def mass_flow_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives for mass flow equation.

        Returns
        -------
        deriv : list
            Matrix with partial derivatives for the fluid equations.
        """
        if self.inl[1].m == self.outl[0].m:
            if self.inl[0].m.is_var:
                self.jacobian[k, self.inl[0].m.J_col] = 1
            if self.outl[1].m.is_var:
                self.jacobian[k, self.outl[1].m.J_col] = -1

        else:
            for i in self.inl:
                if i.m.is_var:
                    self.jacobian[k, i.m.J_col] = 1
            for o in self.outl:
                if o.m.is_var:
                    self.jacobian[k, o.m.J_col] = -1

    def energy_balance_func(self):
        r"""
        Calculate energy balance.

        Returns
        -------
        residual : float
            Residual value of energy balance.

            .. math::

                0 = \sum_i \left(\dot{m}_{in,i} \cdot h_{in,i} \right) -
                \sum_j \left(\dot{m}_{out,j} \cdot h_{out,j} \right)\\
                \forall i \in \text{inlets} \; \forall j \in \text{outlets}
        """
        if self.inl[1].m == self.outl[0].m:
            res = (
                (self.inl[1].h.val_SI - self.outl[0].h.val_SI)
                * self.outl[0].m.val_SI
                + (self.inl[0].h.val_SI - self.outl[1].h.val_SI)
                * self.inl[0].m.val_SI
            )
        else:
            res = 0
            for i in self.inl:
                res += i.m.val_SI * i.h.val_SI
            for o in self.outl:
                res -= o.m.val_SI * o.h.val_SI

        return res

    def energy_balance_deriv(self, increment_filter, k):
        r"""
        Calculate partial derivatives of energy balance.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        # due to topology reduction this is the case quite often
        if self.inl[1].m == self.outl[0].m:
            if self.outl[0].m.is_var:
                self.jacobian[k, self.outl[0].m.J_col] = (self.inl[1].h.val_SI - self.outl[0].h.val_SI)
            if self.inl[1].h.is_var:
                self.jacobian[k, self.inl[1].h.J_col] = self.outl[0].m.val_SI
            if self.outl[0].h.is_var:
                self.jacobian[k, self.outl[0].h.J_col] = -self.outl[0].m.val_SI

            if self.inl[0].m.is_var:
                self.jacobian[k, self.inl[0].m.J_col] = self.inl[0].h.val_SI - self.outl[1].h.val_SI
            if self.inl[0].h.is_var:
                self.jacobian[k, self.inl[0].h.J_col] = self.inl[0].m.val_SI
            if self.outl[1].h.is_var:
                self.jacobian[k, self.outl[1].h.J_col] = -self.outl[1].m.val_SI
        else:
            super().energy_balance_deriv(increment_filter, k)

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

            .. math::

                val = \begin{cases}
                10^6 & \text{key = 'p'}\\
                h\left(p, x=0 \right) & \text{key = 'h' at inlet 1}\\
                h\left(p, x=0.7 \right) & \text{key = 'h' at inlet 2}
                \end{cases}
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            if c.target_id == 'in1':
                return h_mix_pQ(c.p.val_SI, 0, c.fluid_data)
            else:
                return h_mix_pQ(c.p.val_SI, 0.7, c.fluid_data)

    def propagate_wrapper_to_target(self, branch):
        return super().propagate_wrapper_to_target(branch)

    def propagate_to_target(self, branch):

        if branch["connections"][-1].target_id == "in2":
            return

        outconn = self.outl[0]
        subbranch = {
            "connections": [outconn],
            "components": [self, outconn.target],
            "subbranches": {}
        }
        outconn.target.propagate_to_target(subbranch)
        branch["subbranches"][outconn.label] = subbranch

        outconn = self.outl[1]
        if subbranch["components"][-1] == self:

            branch["connections"] += [outconn]
            branch["components"] += [outconn.target]

            outconn.target.propagate_to_target(branch)

        else:
            subbranch = {
                "connections": [outconn],
                "components": [self, outconn.target],
                "subbranches": {}
            }
            outconn.target.propagate_to_target(subbranch)
            branch["subbranches"][outconn.label] = subbranch

    def exergy_balance(self, T0):
        r"""
        Calculate exergy balance of a merge.

        Parameters
        ----------
        T0 : float
            Ambient temperature T0 / K.

        Note
        ----
        Please note, that the exergy balance accounts for physical exergy only.

        .. math::

            \dot{E}_\mathrm{P} = \sum \dot{E}_{\mathrm{out,}j}^\mathrm{PH}\\
            \dot{E}_\mathrm{F} = \sum \dot{E}_{\mathrm{in,}i}^\mathrm{PH}
        """
        self.E_P = self.outl[0].Ex_physical + self.outl[1].Ex_physical
        self.E_F = self.inl[0].Ex_physical + self.inl[1].Ex_physical

        self.E_bus = {
            "chemical": np.nan, "physical": np.nan, "massless": np.nan
        }
        self.E_D = self.E_F - self.E_P
        self.epsilon = self._calc_epsilon()

    def get_plotting_data(self):
        """
        Generate a dictionary containing FluProDia plotting information.

        Returns
        -------
        data : dict
            A nested dictionary containing the keywords required by the
            :code:`calc_individual_isoline` method of the
            :code:`FluidPropertyDiagram` class. The keys :code:`1` and
            :code:`2` connect the saturated liquid-vapor mixture of 'in1' with
            the saturated liquid ('out1') and saturated vapor ('out2'), while
            the keys :code:`3` and :code:`4` connect the (superheated) gas of
            'in2' with the same.
            The key :code:`5` connects both saturated states.
        """
        return {
            1: {
                'isoline_property': 'p',
                'isoline_value': self.inl[0].p.val,
                'isoline_value_end': self.outl[0].p.val,
                'starting_point_property': 'v',
                'starting_point_value': self.inl[0].vol.val,
                'ending_point_property': 'v',
                'ending_point_value': self.outl[0].vol.val
            },
            2: {
                'isoline_property': 'p',
                'isoline_value': self.inl[0].p.val,
                'isoline_value_end': self.outl[1].p.val,
                'starting_point_property': 'v',
                'starting_point_value': self.inl[0].vol.val,
                'ending_point_property': 'v',
                'ending_point_value': self.outl[1].vol.val
            },
            3: {
                'isoline_property': 'p',
                'isoline_value': self.inl[1].p.val,
                'isoline_value_end': self.outl[0].p.val,
                'starting_point_property': 'v',
                'starting_point_value': self.inl[1].vol.val,
                'ending_point_property': 'v',
                'ending_point_value': self.outl[0].vol.val
            },
            4: {
                'isoline_property': 'p',
                'isoline_value': self.inl[1].p.val,
                'isoline_value_end': self.outl[1].p.val,
                'starting_point_property': 'v',
                'starting_point_value': self.inl[1].vol.val,
                'ending_point_property': 'v',
                'ending_point_value': self.outl[1].vol.val
            },
            5: {
                'isoline_property': 'p',
                'isoline_value': self.outl[0].p.val,
                'isoline_value_end': self.outl[1].p.val,
                'starting_point_property': 'v',
                'starting_point_value': self.outl[0].vol.val,
                'ending_point_property': 'v',
                'ending_point_value': self.outl[1].vol.val
            }
        }
