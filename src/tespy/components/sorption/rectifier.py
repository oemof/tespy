# -*- coding: utf-8

"""Module of the rectifier component.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/sorption/rectifier.py

SPDX-License-Identifier: MIT
"""

from tespy.components.component import component_registry
from tespy.components.sorption.base import SorptionBase
from tespy.tools.data_containers import ComponentParameter
from tespy.tools.data_containers import MandatoryConstraint


@component_registry
class Rectifier(SorptionBase):
    r"""
    Ideal rectifier (partial condenser) of absorption cycles with a volatile
    absorbent, e.g. NH3-water.

    The rectifier purifies the vapor leaving the desorber (:code:`in1`) by
    partial condensation. The purified vapor (:code:`out1`) leaves saturated
    at its dew point, the condensate is returned as reflux (:code:`out2`).
    The reflux is the saturated liquid in equilibrium with the entering
    vapor, i.e. it has the composition and temperature of the liquid at the
    dew point of the entering vapor. The heat of condensation is removed,
    :code:`Q` is negative. The composition of the purified vapor is
    typically specified on :code:`out1` (:code:`xi`).

    The reflux returns to the desorber, which closes a pressure loop. The
    rectifier therefore does not impose its pressure on the reflux, the
    reflux pressure follows from the component it enters, e.g. a
    :class:`.Merge` in front of the desorber, or from a specification.

    **Mandatory Equations**

    - mass balance :py:meth:`mass_flow_func`
    - key component balance :py:meth:`key_balance_func`
    - pressure equality of the vapor ports
      :py:meth:`pressure_structure_matrix`
    - purified vapor at its dew point :py:meth:`saturation_func`
    - reflux composition in equilibrium with the entering vapor
      :py:meth:`reflux_composition_func`
    - reflux at its bubble point :py:meth:`reflux_state_func`
    - fluid closure of both outlets :py:meth:`fluid_balance_func`,
      :py:meth:`reflux_fluid_balance_func`

    **Optional Equations**

    - heat flow :py:meth:`heat_func`

    Inlets/Outlets

    - in1: vapor from the desorber
    - out1: purified vapor
    - out2: reflux

    Parameters
    ----------
    label : str
        The label of the component.

    Q : float, dict
        Heat transferred into the rectifier :math:`\dot Q`/W, negative.
    """

    def get_parameters(self):
        return {
            "Q": ComponentParameter(
                max_val=0,
                func=self.heat_func,
                dependents=self.heat_dependents,
                num_eq_sets=1,
                quantity="heat",
                description="heat transferred into the component",
            ),
        }

    @staticmethod
    def inlets():
        return ["in1"]

    @staticmethod
    def outlets():
        return ["out1", "out2"]

    def _solution_inlet(self):
        return self.inl[0]

    def _solution_outlet(self):
        return self.outl[0]

    def _reflux(self):
        return self.outl[1]

    def get_mandatory_constraints(self):
        return {
            "mass_flow_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.mass_flow_func,
                dependents=self.mass_flow_dependents,
                description="mass balance",
            ),
            "pressure_constraints": MandatoryConstraint(
                num_eq_sets=1,
                structure_matrix=self.pressure_structure_matrix,
                description="pressure equality of the vapor ports",
            ),
            "key_balance_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.key_balance_func,
                dependents=self.key_balance_dependents,
                description="key component mass balance",
            ),
            "saturation_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.saturation_func,
                dependents=self.saturation_dependents,
                description="purified vapor at its dew point",
            ),
            "reflux_composition_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.reflux_composition_func,
                dependents=self.reflux_dependents,
                description="reflux in equilibrium with the entering vapor",
            ),
            "reflux_state_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.reflux_state_func,
                dependents=self.reflux_dependents,
                description="reflux at its bubble point",
            ),
            "fluid_balance_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.fluid_balance_func,
                dependents=self.fluid_balance_dependents,
                description="fluid closure of the purified vapor",
            ),
            "reflux_fluid_balance_constraints": MandatoryConstraint(
                num_eq_sets=1,
                func=self.reflux_fluid_balance_func,
                dependents=self.reflux_fluid_balance_dependents,
                description="fluid closure of the reflux",
            ),
        }

    def pressure_structure_matrix(self, k):
        r"""
        Equation for pressure equality of the vapor ports.

        .. math::

            0 = p_\text{in,1} - p_\text{out,1}
        """
        self._structure_matrix[k, self.inl[0].p.sm_col] = 1
        self._structure_matrix[k, self.outl[0].p.sm_col] = -1

    def saturation_func(self):
        r"""
        Equation for the purified vapor at its dew point.

        .. math::

            0 = h_\text{out,1} - h\left(p, Q=1, y_\text{out,1}\right)
        """
        c = self.outl[0]
        return c.h.val_SI - self.pair.h_pQ(c.p.val_SI, 1, self._key_fraction(c))

    def reflux_composition_func(self):
        r"""
        Equation for the reflux composition in equilibrium with the entering
        vapor.

        .. math::

            0 = w_\text{out,2} - w_\text{eq}\left(p, y_\text{in,1}\right)
        """
        vapor = self.inl[0]
        return (
            self._key_fraction(self._reflux())
            - self.pair.w_eq(vapor.p.val_SI, self._key_fraction(vapor))
        )

    def reflux_state_func(self):
        r"""
        Equation for the reflux at its bubble point.

        .. math::

            0 = h_\text{out,2} - h\left(p, Q=0, w_\text{out,2}\right)
        """
        c = self._reflux()
        return c.h.val_SI - self.pair.h_pQ(c.p.val_SI, 0, self._key_fraction(c))

    def reflux_dependents(self):
        vapor = self.inl[0]
        reflux = self._reflux()
        vectors = {}
        vectors.update(self._key_dependents(vapor))
        vectors.update(self._key_dependents(reflux))
        return {
            "scalars": [[vapor.p, reflux.p, reflux.h]],
            "vectors": [vectors],
        }

    def reflux_fluid_balance_func(self):
        r"""
        Equation for the fluid closure of the reflux.

        .. math::

            0 = 1 - \sum_i x_{i,\text{out,2}}
        """
        return sum(self._reflux().fluid.val.values()) - 1

    def reflux_fluid_balance_dependents(self):
        c = self._reflux()
        return {"scalars": [[]], "vectors": [{c.fluid: c.fluid.is_var}]}

    def _update_num_eq(self):
        self.constraints["fluid_balance_constraints"].num_eq = int(
            bool(self.outl[0].fluid.is_var)
        )
        self.constraints["reflux_fluid_balance_constraints"].num_eq = int(
            bool(self._reflux().fluid.is_var)
        )

    def _initial_affine_edges(self):
        first = self.inl[0]
        return [(first.p, c.p, 1.0, 0.0) for c in self.outl]

    def _sorption_inlets(self):
        return self.inl

    def _sorption_outlets(self):
        return self.outl

    def _initial_temperature_edges(self):
        return [
            (self.inl[0], self.outl[0], -50.0, 0.3),
            (self.inl[0], self._reflux(), 0.0, 0.3),
        ]

    def calc_parameters(self):
        self.Q.val_SI = self.calc_Q()
