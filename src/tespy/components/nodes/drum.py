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
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.fluid_properties import h_mix_pQ


@component_registry
class Drum(DropletSeparator):
    r"""
    A drum separates saturated gas from saturated liquid.

    **Mandatory Equations**

    - :py:meth:`tespy.components.nodes.base.NodeBase.mass_flow_func`
    - :py:meth:`tespy.components.nodes.base.NodeBase.pressure_structure_matrix`
    - :py:meth:`tespy.components.nodes.droplet_separator.DropletSeparator.fluid_structure_matrix`
    - :py:meth:`tespy.components.nodes.droplet_separator.DropletSeparator.energy_balance_func`
    - saturated liquid: :py:meth:`tespy.components.nodes.droplet_separator.DropletSeparator.saturated_outlet_func`
    - saturated gas: :py:meth:`tespy.components.nodes.droplet_separator.DropletSeparator.saturated_outlet_func`

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
    >>> import os
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC", "enthalpy": "kJ/kg"
    ... })
    >>> fa = Source('feed ammonia')
    >>> amb_in = Source('air inlet')
    >>> amb_out = Sink('air outlet')
    >>> s = Sink('steam')
    >>> dr = Drum('drum')
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

    >>> char1 = ldc('HeatExchanger', 'kA_char1', 'DEFAULT', CharLine)
    >>> char2 = ldc('HeatExchanger', 'kA_char2', 'EVAPORATING FLUID', CharLine)
    >>> ev.set_attr(pr1=0.999, pr2=0.99, ttd_l=5, kA_char1=char1,
    ...     kA_char2=char2, design=['pr1', 'ttd_l'],
    ...     offdesign=['zeta1', 'kA_char']
    ... )
    >>> ev.set_attr(Q=-1e6)
    >>> erp.set_attr(eta_s=0.8)
    >>> f_dr.set_attr(p=5, T=-5)
    >>> erp_ev.set_attr(m=Ref(f_dr, 4, 0), fluid={'NH3': 1})
    >>> amb_ev.set_attr(fluid={'air': 1}, T=30)
    >>> ev_amb.set_attr(p=1)
    >>> nw.solve('design')
    >>> nw.assert_convergence()
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
    >>> os.remove('tmp.json')
    """

    @staticmethod
    def inlets():
        return ['in1', 'in2']

    @staticmethod
    def outlets():
        return ['out1', 'out2']

    def initialise_target(self, c, key):
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
           return super().initialise_target(c, key)
        elif key == 'h':
            if c.target_id == 'in1':
                return h_mix_pQ(c.p.val_SI, 0, c.fluid_data)
            else:
                return h_mix_pQ(c.p.val_SI, 0.7, c.fluid_data)

    def propagate_wrapper_to_target(self, branch):
        return super().propagate_wrapper_to_target(branch)

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
