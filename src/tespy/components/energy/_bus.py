# -*- coding: utf-8

"""Private base class for energy bus (distribution) components.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/energy/_bus.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import Component
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import SimpleDataContainer as dc_simple


class _EnergyBus(Component):
    """Base for components that distribute energy across many flows.

    Subclasses set :attr:`_energy_port` to :code:`"power"` or :code:`"heat"`.
    Port names follow the pattern :code:`<domain>_in<n>` / :code:`<domain>_out<n>`,
    e.g. :code:`power_in1` for a :class:`PowerBus` and :code:`heat_in1` for a
    :class:`HeatBus`.
    """

    _energy_port: str = None

    @classmethod
    def port_schema(cls):
        is_power = cls._energy_port == "power"
        return {
            "inlets": {"type": "fixed", "ports": []},
            "outlets": {"type": "fixed", "ports": []},
            "powerinlets": {"type": "variable", "parameter": "num_in", "pattern": "power_in{n}", "min": 1} if is_power else {"type": "fixed", "ports": []},
            "poweroutlets": {"type": "variable", "parameter": "num_out", "pattern": "power_out{n}", "min": 1} if is_power else {"type": "fixed", "ports": []},
            "heatinlets": {"type": "fixed", "ports": []} if is_power else {"type": "variable", "parameter": "num_in", "pattern": "heat_in{n}", "min": 1},
            "heatoutlets": {"type": "fixed", "ports": []} if is_power else {"type": "variable", "parameter": "num_out", "pattern": "heat_out{n}", "min": 1},
        }

    def powerinlets(self):
        if self._energy_port == "power":
            return [f"power_in{i + 1}" for i in range(self.num_in.val)]
        return []

    def poweroutlets(self):
        if self._energy_port == "power":
            return [f"power_out{i + 1}" for i in range(self.num_out.val)]
        return []

    def heatinlets(self):
        if self._energy_port == "heat":
            return [f"heat_in{i + 1}" for i in range(self.num_in.val)]
        return []

    def heatoutlets(self):
        if self._energy_port == "heat":
            return [f"heat_out{i + 1}" for i in range(self.num_out.val)]
        return []

    @property
    def _energy_inl(self):
        return self.power_inl if self._energy_port == "power" else self.heat_inl

    @property
    def _energy_outl(self):
        return self.power_outl if self._energy_port == "power" else self.heat_outl

    def get_parameters(self):
        return {
            "num_in": dc_simple(val=0, dtype="int", description="number of inlets"),
            "num_out": dc_simple(val=0, dtype="int", description="number of outlets"),
        }

    def get_mandatory_constraints(self):
        return {
            "energy_balance_constraint": dc_cmc(**{
                "func": self.energy_balance_func,
                "dependents": self.energy_balance_dependents,
                "num_eq_sets": 1,
                "description": "energy balance over all inflows and outflows",
            })
        }

    def energy_balance_func(self):
        r"""
        Equation for energy balance of the component

        Returns
        -------
        residual : float
            Residual value of equation

            .. math::

                0=\sum_{i} \dot E_\text{i} - \sum_{o} \dot E_\text{o}\\
                \forall i \in \text{inlets}, o \in \text{outlets}
        """
        return (
            sum(i.E.val_SI for i in self._energy_inl)
            - sum(o.E.val_SI for o in self._energy_outl)
        )

    def energy_balance_dependents(self):
        return [c.E for c in self._energy_inl + self._energy_outl]
