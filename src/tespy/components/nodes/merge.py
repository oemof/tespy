# -*- coding: utf-8

"""Module of class Merge.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/nodes/merge.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import Component
from tespy.components.nodes.node import Node
from tespy.tools.data_containers import dc_simple


class Merge(Node):
    r"""
    The component Node is the parent class for Splitter, Separator and Merge.

    Equations

        **mandatory equations**

        - :py:meth:`tespy.components.component.Component.mass_flow_func`

        .. math::

            0 = p_{in} - p_{out,i} \;
            \forall i \in \mathrm{outlets}

        **additional equations**

        - :py:meth:`tespy.components.basics.merge.additional_equations`

    Inlets/Outlets

        - specify number of outlets with :code:`num_in` (default value: 2)
        - out1

    Image

        .. image:: _images/merge.svg
           :scale: 100 %
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

    design_path: str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings: boolean
        Ignore warnings on default characteristics usage for this component.

    printout: boolean
        Include this component in the network's results printout.

    num_in : float/tespy.tools.data_containers.dc_simple
        Number of inlets for this component, default value: 2.

    Example
    -------
    The merge mixes a specified number of mass flows and has a single outlet.
    At the outlet, fluid composition and enthalpy are calculated by mass
    weighted fluid composition and enthalpy of the inlets.

    >>> from tespy.components import Sink, Source, Merge
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import shutil
    >>> import numpy as np
    >>> fluid_list = ['O2', 'N2']
    >>> nw = Network(fluids=fluid_list, p_unit='bar', iterinfo=False)
    >>> so1 = Source('source1')
    >>> so2 = Source('source2')
    >>> so3 = Source('source3')
    >>> si1 = Sink('sink')
    >>> m = Merge('merge', num_in=3)
    >>> m.component()
    'merge'
    >>> inc1 = Connection(so1, 'out1', m, 'in1')
    >>> inc2 = Connection(so2, 'out1', m, 'in2')
    >>> inc3 = Connection(so3, 'out1', m, 'in3')
    >>> outg = Connection(m, 'out1', si1, 'in1')
    >>> nw.add_conns(inc1, inc2, inc3, outg)

    A merge with three inlets mixes air (simplified) with pure nitrogen and
    pure oxygen. All gases enter the component at the same temperature. As
    mixing effects are not considered, the outlet temperature should thus be
    similar to the three inlet temperatures (difference might occur due to
    rounding in fluid property functions, let's check it for two different
    temperatures). It is e.g. possible to find the required mass flow of pure
    nitrogen given the nitrogen mass fraction in the outlet.

    >>> T = 293.15
    >>> inc1.set_attr(fluid={'O2': 0.23, 'N2': 0.77}, p=1, T=T, m=5)
    >>> inc2.set_attr(fluid={'O2': 1, 'N2':0}, T=T, m=5)
    >>> inc3.set_attr(fluid={'O2': 0, 'N2': 1}, T=T)
    >>> outg.set_attr(fluid={'N2': 0.4})
    >>> nw.solve('design')
    >>> round(inc3.m.val_SI, 2)
    0.25
    >>> abs(round((outg.T.val_SI - T) / T, 5)) < 0.01
    True
    >>> T = 173.15
    >>> inc1.set_attr(T=T)
    >>> inc2.set_attr(T=T)
    >>> inc3.set_attr(T=T)
    >>> nw.solve('design')
    >>> abs(round((outg.T.val_SI - T) / T, 5)) < 0.01
    True
    """

    @staticmethod
    def component():
        return 'merge'

    @staticmethod
    def attr():
        return {'num_in': dc_simple()}

    def inlets(self):
        if self.num_in.is_set:
            return ['in' + str(i + 1) for i in range(self.num_in.val)]
        else:
            self.set_attr(num_in=2)
            return self.inlets()

    @staticmethod
    def outlets():
        return ['out1']

    def comp_init(self, nw):

        Component.comp_init(self, nw)

        # number of mandatroy equations for
        # mass flow: 1
        # pressure: number of inlets + number of outlets - 1
        # fluid: number of fluids
        # enthalpy: 1

        self.num_eq = self.num_i + self.num_o + self.num_nw_fluids + 1

        self.jacobian = np.zeros((
            self.num_eq,
            self.num_i + self.num_o + self.num_vars,
            self.num_nw_vars))

        self.residual = np.zeros(self.num_eq)
        self.jacobian[0:1] = self.mass_flow_deriv()
        end = self.num_i + self.num_o
        self.jacobian[1:end] = self.pressure_deriv()

    def additional_equations(self, k):
        r"""
        Calculate results of additional equations.

        Equations

            **mandatroy equations**

            .. math::

                0 = \dot{m}_{in_{j}} \cdot fluid_{i,in_{j}} -
                    \dot {m}_{out} \cdot fluid_{i,out} \\
                \forall i \in \mathrm{fluid}, \; \forall j \in inlets

            .. math::

                0 = h_{in} - h_{out,i} \;
                \forall i \in \mathrm{outlets}\\
        """
        ######################################################################
        # equations for fluid balance
        for fluid, x in self.outl[0].fluid.val.items():
            res = -x * self.outl[0].m.val_SI
            for i in self.inl:
                res += i.fluid.val[fluid] * i.m.val_SI
            self.residual[k] = res
            k += 1

        ######################################################################
        # equation for energy balance
        h_res = -self.outl[0].m.val_SI * self.outl[0].h.val_SI
        for i in self.inl:
            h_res += i.m.val_SI * i.h.val_SI
        self.residual[k] = h_res
        k += 1

    def additional_derivatives(self, increment_filter, k):
        r"""Calculate partial derivatives for given additional equations."""
        ######################################################################
        # derivatives for fluid balance equations
        i = 0
        for fluid, x in self.outl[0].fluid.val.items():
            j = 0
            for inl in self.inl:
                self.jacobian[k, j, 0] = inl.fluid.val[fluid]
                self.jacobian[k, j, i + 3] = inl.m.val_SI
                j += 1
            self.jacobian[k, j, 0] = -x
            self.jacobian[k, j, i + 3] = -self.outl[0].m.val_SI
            i += 1
            k += 1

        ######################################################################
        # derivatives for energy balance equations
        self.jacobian[k, self.num_i, 0] = -self.outl[0].h.val_SI
        self.jacobian[k, self.num_i, 2] = -self.outl[0].m.val_SI
        j = 0
        for i in self.inl:
            self.jacobian[k, j, 0] = i.h.val_SI
            self.jacobian[k, j, 2] = i.m.val_SI
            j += 1
        k += 1

    def get_plotting_data(self):
        """Generate a dictionary containing FluProDia plotting information.

        Returns
        -------
        data : dict
            A nested dictionary containing the keywords required by the
            :code:`calc_individual_isoline` method of the
            :code:`FluidPropertyDiagram` class. First level keys are the
            connection index ('in1' -> 'out1', therefore :code:`1` etc.).
        """
        return {
            i + 1: {
                'isoline_property': 'p',
                'isoline_value': self.inl[i].p.val,
                'isoline_value_end': self.outl[0].p.val,
                'starting_point_property': 'v',
                'starting_point_value': self.inl[i].vol.val,
                'ending_point_property': 'v',
                'ending_point_value': self.outl[0].vol.val
            } for i in range(self.num_i)}
