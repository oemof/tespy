# -*- coding: utf-8

"""Module of class Motor.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/power/motor.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.energy._converter import _EnergyConverter


@component_registry
class Motor(_EnergyConverter):
    r"""
    A motor converts electrical energy into mechanical energy.

    **Mandatory Equations**

    - None

    **Optional Equations**

    - :py:meth:`tespy.components.power.motor.Motor.eta_func`
    - :py:meth:`tespy.components.power.motor.Motor.delta_power_func`
    - :py:meth:`tespy.components.power.motor.Motor.eta_char_func`

    Inlets/Outlets

    - None

    Optional inlets/outlets

    - power_in
    - power_out

    Image

    .. image:: /api/_images/Motor.svg
       :alt: flowsheet of the motor
       :align: center
       :class: only-light

    .. image:: /api/_images/Motor_darkmode.svg
       :alt: flowsheet of the motor
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

    eta : float, dict
        Outlet to inlet efficiency, :math:`\eta/1`

    delta_power : float, dict
        Fixed power offset, :math:`\text{delta_power}/\text{W}`

    eta_char : tespy.tools.characteristics.CharLine, dict
        Characteristic line for efficiency to power as function of design
        efficiency.

    Example
    -------
    A compressor provides compressed air which is used in a compressed air
    distribution system. The energy is provided by an electrical motor.

    >>> from tespy.components import Sink, Source, Compressor, Motor, PowerSource
    >>> from tespy.connections import Connection, PowerConnection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC"
    ... })
    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> compressor = Compressor('compressor')

    Ambient air flows into the compressor and is ejected at 4 bar. We can set
    the system up without the use of any of the power components.

    >>> c1 = Connection(so, 'out1', compressor, 'in1')
    >>> c2 = Connection(compressor, 'out1', si, 'in1')
    >>> nw.add_conns(c1, c2)
    >>> c1.set_attr(fluid={'air': 1}, T=25, p=1, m=1)
    >>> c2.set_attr(p=4)
    >>> compressor.set_attr(eta_s=0.8)
    >>> nw.solve('design')

    We can add the Motor and a PowerSource and then connect these parts to
    the compressor.

    >>> motor = Motor('motor')
    >>> power_source = PowerSource('power source')
    >>> e1 = PowerConnection(power_source, 'power', motor, 'power_in')
    >>> e2 = PowerConnection(motor, 'power_out', compressor, 'power')
    >>> nw.add_conns(e1, e2)

    Now we have added two variables to our problem (the power flows of e1 and
    e2), but only one equation (the power balance for the compressor). The
    connection between the two power flows can be made through specifying the
    efficiency of the motor:

    >>> motor.set_attr(eta=.98)
    >>> nw.solve('design')
    >>> nw.assert_convergence()
    >>> round(e2.E.val_SI) == round(compressor.P.val)
    True
    >>> round(e1.E.val_SI) == round(compressor.P.val / 0.98)
    True

    We could also specify the electrical energy instead of fixing the air
    mass flow to calculate the resulting air mass flow:

    >>> e1.set_attr(E=1e5)
    >>> c1.set_attr(m=None)
    >>> nw.solve('design')
    >>> round(c1.m.val, 3)
    0.539

    Or, fix both (electrical and mechanical power flows) and leave open the
    motor efficiency:

    >>> e2.set_attr(E=0.9e5)
    >>> motor.set_attr(eta=None)
    >>> nw.solve('design')
    >>> round(motor.eta.val, 2)
    0.9

    >>> e2.set_attr(E=None)
    >>> motor.set_attr(delta_power=5e3)
    >>> nw.solve('design')
    >>> round(motor.eta.val, 3)
    0.95
    """
    @classmethod
    def port_schema(cls):
        return {
            "inlets": {"type": "fixed", "ports": []},
            "outlets": {"type": "fixed", "ports": []},
            "powerinlets": {"type": "fixed", "ports": ["power_in"]},
            "poweroutlets": {"type": "fixed", "ports": ["power_out"]},
            "heatinlets": {"type": "fixed", "ports": []},
            "heatoutlets": {"type": "fixed", "ports": []},
        }
