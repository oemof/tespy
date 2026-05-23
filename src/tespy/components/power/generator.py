# -*- coding: utf-8

"""Module of class Generator.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/power/generator.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.energy._converter import _EnergyConverter


@component_registry
class Generator(_EnergyConverter):
    r"""
    A generator converts mechanical energy into electrical energy.

    .. image:: /api/_images/components/Generator.svg
       :alt: flowsheet of the generator
       :align: center
       :class: only-light

    .. image:: /api/_images/components/Generator_darkmode.svg
       :alt: flowsheet of the generator
       :align: center
       :class: only-dark

    Ports
    -----

    - Power inlets: power_in
    - Power outlets: power_out

    Mandatory Equations
    -------------------

    None

    Parameters
    ----------

    char_warnings : bool
        Ignore warnings on default characteristics usage for this component.

    delta_power : float, dict
        Inlet to outlet power difference. Quantity: :code:`power`.
        Equation: :py:meth:`delta_power_func <tespy.components.energy._converter._EnergyConverter.delta_power_func>`.

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the components design case.

    eta : float, dict
        Efficiency. Quantity: :code:`efficiency`.
        Equation: :py:meth:`eta_func <tespy.components.energy._converter._EnergyConverter.eta_func>`.

    eta_char : tespy.tools.characteristics.CharLine, dict
        Efficiency lookup table for offdesign.
        Equation: :py:meth:`eta_char_func <tespy.components.energy._converter._EnergyConverter.eta_char_func>`.

    label : str
        The label of the component.

    local_design : bool
        Treat this component in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this component in offdesign mode in a design calculation.

    offdesign : list
        List containing offdesign parameters (stated as String).

    printout : bool
        Include this component in the network's results printout.

    Example
    -------
    A turbine generates mechanical power which is used to generate electrical
    power by the generator.

    >>> from tespy.components import Sink, Source, Turbine, Generator, PowerSink
    >>> from tespy.connections import Connection, PowerConnection
    >>> from tespy.networks import Network
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "pressure_difference": "bar",
    ...     "temperature": "degC"
    ... })
    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> turbine = Turbine('turbine')

    Steam flows through the turbine and we can set it up as we are used to for
    systems without power components.

    >>> c1 = Connection(so, 'out1', turbine, 'in1')
    >>> c2 = Connection(turbine, 'out1', si, 'in1')
    >>> nw.add_conns(c1, c2)
    >>> c1.set_attr(fluid={'water': 1}, T=500, p=50, m=1)
    >>> c2.set_attr(p=5)
    >>> turbine.set_attr(eta_s=0.9)
    >>> nw.solve('design')

    We can add the Generator and a PowerSink and then connect these parts to
    the turbine.

    >>> generator = Generator('generator')
    >>> power_sink = PowerSink('power sink')
    >>> e1 = PowerConnection(turbine, 'power', generator, 'power_in')
    >>> e2 = PowerConnection(generator, 'power_out', power_sink, 'power')
    >>> nw.add_conns(e1, e2)

    Now we have added two variables to our problem (the power flows of e1 and
    e2), but only one equation (the power balance for the turbine). The
    connection between the two power flows can be made through specifying the
    efficiency of the generator:

    >>> generator.set_attr(eta=.98)
    >>> nw.solve('design')
    >>> nw.assert_convergence()
    >>> round(e1.E.val_SI) == -round(turbine.P.val)
    True
    >>> round(e2.E.val_SI) == -round(turbine.P.val * 0.98)
    True

    We could also specify the electrical energy instead of fixing the steam
    mass flow to calculate the resulting steam mass flow:

    >>> e2.set_attr(E=1e6)
    >>> c1.set_attr(m=None)
    >>> nw.solve('design')
    >>> round(c1.m.val, 3)
    1.837

    Or, fix both (electrical and mechanical power flows) and leave open the
    generator efficiency:

    >>> e1.set_attr(E=1.1e6)
    >>> generator.set_attr(eta=None)
    >>> nw.solve('design')
    >>> round(generator.eta.val, 2)
    0.91

    >>> e1.set_attr(E=None)
    >>> generator.set_attr(delta_power=50e3)
    >>> nw.solve('design')
    >>> round(generator.eta.val, 3)
    0.952
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
