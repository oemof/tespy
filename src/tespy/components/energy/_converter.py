# -*- coding: utf-8

"""Private base class for power-domain energy converter components.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/energy/_converter.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import Component
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp


class _EnergyConverter(Component):
    """Base for components that convert power from one flow to another.

    Both Motor and Generator have one ``power_in`` and one ``power_out``
    connection, an efficiency equation and an optional characteristic curve.
    The characteristic is always normalised by the outlet flow.
    """

    @staticmethod
    def powerinlets():
        return ["power_in"]

    @staticmethod
    def poweroutlets():
        return ["power_out"]

    @staticmethod
    def get_mandatory_constraints():
        return {}

    def get_parameters(self):
        return {
            "eta": dc_cp(**{
                "structure_matrix": self.eta_structure_matrix,
                "func": self.eta_func,
                "dependents": self.eta_dependents,
                "num_eq_sets": 1,
                "max_val": 1,
                "min_val": 0,
                "quantity": "efficiency",
                "description": "efficiency"
            }),
            "delta_power": dc_cp(**{
                "structure_matrix": self.delta_power_structure_matrix,
                "func": self.delta_power_func,
                "dependents": self.delta_power_dependents,
                "num_eq_sets": 1,
                "min_val": 0,
                "quantity": "power",
                "description": "inlet to outlet power difference"
            }),
            "eta_char": dc_cc(**{
                "func": self.eta_char_func,
                "dependents": self.eta_char_dependents,
                "num_eq_sets": 1,
                "description": "efficiency lookup table for offdesign"
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
        Equation for efficiency characteristics of the component.

        The characteristic is normalised by the outlet flow.

        Returns
        -------
        residual : float
            Residual value of equation

            .. math::

                0=\dot E_\text{in} \cdot \eta_\text{design} \cdot
                f\left(\frac{\dot E_\text{out}}{\dot E_\text{out,design}}\right)
                - \dot E_\text{out}
        """
        expr = (
            self.power_outl[0].E.val_SI
            / self._conn_design(self.power_outl[0], "E")
        )
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
