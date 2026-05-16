# -*- coding: utf-8

"""Private base class for energy sink boundary components.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/energy/_sink.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import Component


class _EnergySink(Component):
    """Base for boundary components that absorb a single energy flow.

    Subclasses set :attr:`_energy_port` to :code:`"power"` or :code:`"heat"` to
    declare which connection domain they participate in.
    """

    _energy_port: str = None

    def powerinlets(self):
        return ["power"] if self._energy_port == "power" else []

    def heatinlets(self):
        return ["heat"] if self._energy_port == "heat" else []

    @staticmethod
    def get_mandatory_constraints():
        return {}
