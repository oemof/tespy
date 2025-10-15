# -*- coding: utf-8

"""Module of class Generator.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/power/generator.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import Component
from tespy.components.component import component_registry
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp


@component_registry
class Generator(Component):
    r"""
    A generator converts mechanical energy into electrical energy.

    **Mandatory Equations**

    - None

    **Optional Equations**

    - :py:meth:`tespy.components.power.generator.Generator.eta_func`
    - :py:meth:`tespy.components.power.generator.Generator.delta_power_func`
    - :py:meth:`tespy.components.power.generator.Generator.eta_char_func`

    Inlets/Outlets

    - None

    Optional inlets/outlets

    - power_in
    - power_out

    Image

    .. image:: /api/_images/Generator.svg
       :alt: flowsheet of the generator
       :align: center
       :class: only-light

    .. image:: /api/_images/Generator_darkmode.svg
       :alt: flowsheet of the generator
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
    A turbine generates mechanical power which is used to generate electrical
    power by the generator.

    >>> from tespy.components import Sink, Source, Turbine, Generator, PowerSink
    >>> from tespy.connections import Connection, PowerConnection
    >>> from tespy.networks import Network
    >>> import os
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC"
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

    def powerinlets(self):
        return ["power_in"]

    def poweroutlets(self):
        return ["power_out"]

    def get_parameters(self):
        return {
            "eta": dc_cp(**{
                "structure_matrix": self.eta_structure_matrix,
                "func": self.eta_func,
                "dependents": self.eta_dependents,
                "num_eq_sets": 1,
                "max_val": 1,
                "min_val": 0,
                "quantity": "efficiency"
            }),
            "delta_power": dc_cp(**{
                "structure_matrix": self.delta_power_structure_matrix,
                "func": self.delta_power_func,
                "dependents": self.delta_power_dependents,
                "num_eq_sets": 1,
                "min_val": 0,
                "quantity": "power"
            }),
            "eta_char": dc_cc(**{
                "func": self.eta_char_func,
                "dependents": self.eta_char_dependents,
                "num_eq_sets": 1
            })
        }

    def eta_func(self):
        r"""
        Equation for efficiency of the component

        Returns
        -------
        residual : float
            Residual value of equation

            .. math::

                0=\dot E_\text{in} \cdot \eta - \dot E_\text{out}
        """
        return (
            self.power_inl[0].E.val_SI * self.eta.val_SI
            - self.power_outl[0].E.val_SI
        )

    def eta_structure_matrix(self, k):
        self._structure_matrix[k, self.power_inl[0].E.sm_col] = self.eta.val_SI
        self._structure_matrix[k, self.power_outl[0].E.sm_col] = -1

    def eta_dependents(self):
        return [self.power_inl[0].E, self.power_outl[0].E]

    def delta_power_func(self):
        r"""
        Equation for power delta of the component

        Returns
        -------
        residual : float
            Residual value of equation

            .. math::

                0=\dot E_\text{in} - \dot E_\text{out} - \Delta \dot E
        """
        return (
            self.power_inl[0].E.val_SI - self.power_outl[0].E.val_SI
            - self.delta_power.val_SI
        )

    def delta_power_structure_matrix(self, k):
        self._structure_matrix[k, self.power_inl[0].E.sm_col] = 1
        self._structure_matrix[k, self.power_outl[0].E.sm_col] = -1
        self._rhs[k] = self.delta_power.val_SI

    def delta_power_dependents(self):
        return [self.power_inl[0].E, self.power_outl[0].E]

    def eta_char_func(self):
        r"""
        Equation for efficiency characteristics of the component

        Returns
        -------
        residual : float
            Residual value of equation

            .. math::

                0=\dot E_\text{in} \cdot \eta_\text{design} \cdot
                f\left(\frac{\dot E_\text{out}}{\dot E_\text{out,design}}\right)
                - \dot E_\text{out}
        """
        expr = self.power_outl[0].E.val_SI / self.power_outl[0].E.design
        f = self.eta_char.char_func.evaluate(expr)
        return (
            self.power_inl[0].E.val_SI * self.eta.design * f
            - self.power_outl[0].E.val_SI
        )

    def eta_char_dependents(self):
        return [self.power_inl[0].E, self.power_outl[0].E]

    def calc_parameters(self):
        self.eta.val_SI = self.power_outl[0].E.val_SI / self.power_inl[0].E.val_SI
        self.delta_power.val_SI = (
            self.power_inl[0].E.val_SI - self.power_outl[0].E.val_SI
        )
