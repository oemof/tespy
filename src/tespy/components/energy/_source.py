# -*- coding: utf-8

"""Private base class for energy source boundary components.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/energy/_source.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import Component


class _EnergySource(Component):
    """Base for boundary components that emit a single energy flow.

    Subclasses set :attr:`_energy_port` to :code:`"power"` or :code:`"heat"` to
    declare which connection domain they participate in.
    """

    _energy_port: str = None

    @classmethod
    def port_schema(cls):
        is_power = cls._energy_port == "power"
        return {
            "inlets": {"type": "fixed", "ports": []},
            "outlets": {"type": "fixed", "ports": []},
            "powerinlets": {"type": "fixed", "ports": []},
            "poweroutlets": {"type": "fixed", "ports": ["power"] if is_power else []},
            "heatinlets": {"type": "fixed", "ports": []},
            "heatoutlets": {"type": "fixed", "ports": [] if is_power else ["heat"]},
        }

    def poweroutlets(self):
        return ["power"] if self._energy_port == "power" else []

    def heatoutlets(self):
        return ["heat"] if self._energy_port == "heat" else []

    @staticmethod
    def get_mandatory_constraints():
        return {}
