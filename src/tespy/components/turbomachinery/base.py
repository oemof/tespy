# -*- coding: utf-8

"""Module of class Turbomachine.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/turbomachinery/base.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import Component
from tespy.components.component import component_registry
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.helpers import _numeric_deriv


@component_registry
class Turbomachine(Component):
    r"""
    Parent class for compressor, pump and turbine.

    Ports
    -----

    - Fluid inlets: in1
    - Fluid outlets: out1

    Mandatory Equations
    -------------------

    - mass flow equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
    - fluid composition equality constraint(s): :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`

    Parameters
    ----------

    char_warnings : bool
        Ignore warnings on default characteristics usage for this component.

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the components design case.

    dp : float, dict
        Inlet to outlet absolute pressure change. Quantity:
        :code:`pressure_difference`.
        Equation: :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`.

    label : str
        The label of the component.

    local_design : bool
        Treat this component in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this component in offdesign mode in a design calculation.

    offdesign : list
        List containing offdesign parameters (stated as String).

    P : float, dict
        Power input/output of the component. Quantity: :code:`power`.
        Equation: :py:meth:`energy_balance_func <tespy.components.turbomachinery.base.Turbomachine.energy_balance_func>`.

    pr : float, dict
        Outlet to inlet pressure ratio. Quantity: :code:`ratio`.
        Equation: :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`.

    printout : bool
        Include this component in the network's results printout.

    Example
    -------
    For an example please refer to:

    - :class:`tespy.components.turbomachinery.compressor.Compressor`
    - :class:`tespy.components.turbomachinery.pump.Pump`
    - :class:`tespy.components.turbomachinery.turbine.Turbine`
    - :class:`tespy.components.turbomachinery.steam_turbine.SteamTurbine`
    """
    def _calc_P(self):
        return self.inl[0].m.val_SI * (self.outl[0].h.val_SI - self.inl[0].h.val_SI)

    def get_parameters(self):
        return {
            'P': dc_cp(
                num_eq_sets=1,
                func=self.energy_balance_func,
                dependents=self.energy_balance_dependents,
                quantity="power",
                description="power input/output of the component",
                calc=self._calc_P
            ),
            'pr': dc_cp(
                num_eq_sets=1,
                func_params={'pr': 'pr'},
                structure_matrix=self.pr_structure_matrix,
                quantity="ratio",
                description="outlet to inlet pressure ratio",
                calc=self._calc_pr
            ),
            'dp': dc_cp(
                num_eq_sets=1,
                structure_matrix=self.dp_structure_matrix,
                func_params={'dp': 'dp'},
                quantity="pressure_difference",
                description="inlet to outlet absolute pressure change",
                calc=self._calc_dp
            )
        }

    def get_bypass_constraints(self):
        return {
            'mass_flow_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': self.num_i,
                'func_params': {'variable': 'm'}
            }),
            'pressure_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': self.num_i,
                'func_params': {'variable': 'p'}
            }),
            'enthalpy_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': self.num_i,
                'func_params': {'variable': 'h'}
            }),
            'fluid_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': self.num_i,
                'func_params': {'variable': 'fluid'}
            })
        }

    @staticmethod
    def inlets():
        return ['in1']

    @staticmethod
    def outlets():
        return ['out1']

    def energy_balance_func(self):
        r"""
        Calculate energy balance of a turbomachine.

        Returns
        -------
        residual : float
            Residual value of turbomachine energy balance

            .. math::

                0=\dot{m}_{in}\cdot\left(h_{out}-h_{in}\right)-P
        """
        return self._calc_P() - self.P.val_SI

    def energy_balance_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].h,
            self.outl[0].h,
        ]

    def entropy_balance(self):
        r"""
        Calculate entropy balance of turbomachine.

        Note
        ----
        The entropy balance makes the following parameter available:

        .. math::

            \text{S\_irr}=\dot{m} \cdot \left(s_\text{out}-s_\text{in}
            \right)\\
        """
        self.S_irr = self.inl[0].m.val_SI * (
            self.outl[0].s.val_SI - self.inl[0].s.val_SI
        )

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
            1: {
                'isoline_property': 's',
                'isoline_value': self.inl[0].s.val,
                'isoline_value_end': self.outl[0].s.val,
                'starting_point_property': 'vol',
                'starting_point_value': self.inl[0].vol.val,
                'ending_point_property': 'vol',
                'ending_point_value': self.outl[0].vol.val
            }
        }
