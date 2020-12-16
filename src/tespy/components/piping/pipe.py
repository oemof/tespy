# -*- coding: utf-8

"""Module of class Pipe.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/piping/pipe.py

SPDX-License-Identifier: MIT
"""

from tespy.components.heat_exchangers.heat_exchanger_simple import HeatExchangerSimple


class Pipe(HeatExchangerSimple):
    r"""
    The Pipe is a subclass of a HeatExchangerSimple.

    Equations

        **mandatory equations**

        - :py:meth:`tespy.components.component.Component.fluid_func`
        - :py:meth:`tespy.components.component.Component.mass_flow_func`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) -
            \dot{Q}

            0 = p_{in} \cdot pr - p_{out}

        - :py:meth:`tespy.components.component.Component.zeta_func`

        - :py:meth:`tespy.components.heat_exchangers.heat_exchanger_simple.HeatExchangerSimple.darcy_func`
          or :py:meth:`tespy.components.heat_exchangers.heat_exchanger_simple.HeatExchangerSimple.hw_func`

        **additional equations**

        - :py:meth:`tespy.components.heat_exchangers.heat_exchanger_simple.HeatExchangerSimple.additional_equations`

    Inlets/Outlets

        - in1
        - out1

    Image

        .. image:: _images/Pipe.svg
           :alt: alternative text
           :align: center

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

    Q : str, float, tespy.tools.data_containers.ComponentProperties
        Heat transfer, :math:`Q/\text{W}`.

    pr : str, float, tespy.tools.data_containers.ComponentProperties
        Outlet to inlet pressure ratio, :math:`pr/1`.

    zeta : str, float, tespy.tools.data_containers.ComponentProperties
        Geometry independent friction coefficient,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    D : str, float, tespy.tools.data_containers.ComponentProperties
        Diameter of the pipes, :math:`D/\text{m}`.

    L : str, float, tespy.tools.data_containers.ComponentProperties
        Length of the pipes, :math:`L/\text{m}`.

    ks : str, float, tespy.tools.data_containers.ComponentProperties
        Pipe's roughness, :math:`ks/\text{m}` for darcy friction,
        :math:`ks/\text{1}` for hazen-williams equation.

    hydro_group : str, tespy.tools.data_containers.GroupedComponentProperties
        Parametergroup for pressure drop calculation based on pipes dimensions.
        Choose 'HW' for hazen-williams equation, else darcy friction factor is
        used.

    kA : str, float, tespy.tools.data_containers.ComponentProperties
        Area independent heat transition coefficient,
        :math:`kA/\frac{\text{W}}{\text{K}}`.

    kA_char : tespy.tools.characteristics.CharLine, tespy.tools.data_containers.ComponentCharacteristics
        Characteristic line for heat transfer coefficient.

    Tamb : float, tespy.tools.data_containers.ComponentProperties
        Ambient temperature, provide parameter in network's temperature
        unit.

    kA_group : tespy.tools.data_containers.GroupedComponentProperties
        Parametergroup for heat transfer calculation from ambient temperature
        and area independent heat transfer coefficient kA.

    Example
    -------
    A mass flow of 10 kg/s ethanol is transported in a pipeline. The pipe is
    considered adiabatic and has a length of 500 meters. We can calculate the
    diameter required at a given pressure loss of 2.5 %. After we determined
    the required diameter, we can predict pressure loss at a different mass
    flow through the pipeline.

    >>> from tespy.components import Sink, Source, Pipe
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import shutil
    >>> fluid_list = ['ethanol']
    >>> nw = Network(fluids=fluid_list)
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> so = Source('source 1')
    >>> si = Sink('sink 1')
    >>> pi = Pipe('pipeline')
    >>> pi.component()
    'pipe'
    >>> pi.set_attr(pr=0.975, Q=0, design=['pr'], L=100, D='var', ks=5e-5)
    >>> inc = Connection(so, 'out1', pi, 'in1')
    >>> outg = Connection(pi, 'out1', si, 'in1')
    >>> nw.add_conns(inc, outg)
    >>> inc.set_attr(fluid={'ethanol': 1}, m=10, T=30, p=3)
    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> round(pi.D.val, 3)
    0.119
    >>> outg.p.val / inc.p.val == pi.pr.val
    True
    >>> inc.set_attr(m=15)
    >>> pi.set_attr(D=pi.D.val)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(pi.pr.val, 2)
    0.94
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    @staticmethod
    def component():
        return 'pipe'
