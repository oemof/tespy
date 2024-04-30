# -*- coding: utf-8

"""Module for class SubsystemInterface.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/basics/subsystem_interface.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import Component
from tespy.components.component import component_registry
from tespy.tools.data_containers import SimpleDataContainer as dc_simple


@component_registry
class SubsystemInterface(Component):
    r"""
    The subsystem interface does not change fluid properties.

    **Mandatory Equations**

    - :py:meth:`tespy.components.component.Component.fluid_func`
    - :py:meth:`tespy.components.component.Component.mass_flow_func`
    - :py:meth:`tespy.components.component.Component.pressure_equality_func`
    - :py:meth:`tespy.components.component.Component.enthalpy_equality_func`

    Inlets/Outlets

    - Specify number of inlets and outlets with :code:`num_inter`,
      predefined value: 1.

    Image

    .. image:: /api/_images/SubsystemInterface.svg
       :alt: flowsheet of the subsystem interface
       :align: center
       :class: only-light

    .. image:: /api/_images/SubsystemInterface_darkmode.svg
       :alt: flowsheet of the subsystem interface
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

    num_inter : float, dict
        Number of interfaces for subsystem.

    Note
    ----
    This component passes all fluid properties and mass flow from its inlet to
    the outlet.

    Example
    -------
    As connections can only connect a component with a different
    component, the subsystem interface is used to connect subsystems with the
    rest of your network. It is necessary to specify the number of interfaces
    of the subsystem interface, if you want any number other than 1. We will
    not go in depth of subsystem usage in this example. Please refer to
    :ref:`this section <tespy_subsystems_label>` for more information on
    building your own subsystems.

    >>> from tespy.components import Sink, Source, SubsystemInterface
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> nw = Network()
    >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg', iterinfo=False)
    >>> so1 = Source('source 1')
    >>> si1 = Sink('sink 1')
    >>> so2 = Source('source 2')
    >>> si2 = Sink('sink 2')
    >>> IF = SubsystemInterface('subsystem interface', num_inter=2)
    >>> IF.component()
    'subsystem interface'
    >>> len(IF.inlets())
    2

    The interface does not change the fluid properties in any way.

    >>> inc1 = Connection(so1, 'out1', IF, 'in1')
    >>> outg1 = Connection(IF, 'out1', si1, 'in1')
    >>> inc2 = Connection(so2, 'out1', IF, 'in2')
    >>> outg2 = Connection(IF, 'out2', si2, 'in1')
    >>> nw.add_conns(inc1, outg1, inc2, outg2)
    >>> inc1.set_attr(fluid={'H2O': 1}, T=40, p=3, m=100)
    >>> inc2.set_attr(fluid={'N2': 1}, T=60, p=1, v=10)
    >>> nw.solve('design')
    >>> inc1.m.val_SI == outg1.m.val_SI
    True
    >>> inc2.m.val_SI == outg2.m.val_SI
    True
    >>> inc1.h.val_SI == outg1.h.val_SI
    True
    >>> inc2.h.val_SI == outg2.h.val_SI
    True
    """

    @staticmethod
    def component():
        return 'subsystem interface'

    def get_mandatory_constraints(self):
        return {
            'pressure_equality_constraints': {
                'func': self.pressure_equality_func,
                'deriv': self.pressure_equality_deriv,
                'constant_deriv': True,
                'latex': self.pressure_equality_func_doc,
                'num_eq': self.num_i},
            'enthalpy_equality_constraints': {
                'func': self.enthalpy_equality_func,
                'deriv': self.enthalpy_equality_deriv,
                'constant_deriv': True,
                'latex': self.enthalpy_equality_func_doc,
                'num_eq': self.num_i}
        }

    @staticmethod
    def get_parameters():
        return {'num_inter': dc_simple()}

    def inlets(self):
        if self.num_inter.is_set:
            return ['in' + str(i + 1) for i in range(self.num_inter.val)]
        else:
            return ['in1']

    def outlets(self):
        if self.num_inter.is_set:
            return ['out' + str(i + 1) for i in range(self.num_inter.val)]
        else:
            return ['out1']
