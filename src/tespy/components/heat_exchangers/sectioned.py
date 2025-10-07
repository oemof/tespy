# -*- coding: utf-8

"""Module of class MovingBoundaryCondenser.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/heat_exchangers/sectioned.py

SPDX-License-Identifier: MIT
"""
import math

import numpy as np

from tespy.components.heat_exchangers.base import HeatExchanger
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import T_mix_ph


class SectionedHeatExchanger(HeatExchanger):

    def get_parameters(self):
        params = super().get_parameters()
        params.update({
            'UA': dc_cp(
                min_val=0, num_eq_sets=1,
                func=self.UA_func,
                dependents=self.UA_dependents,
                quantity="heat_transfer_coefficient"
            ),
            'refrigerant_index': dc_simple(val=0),
            're_exp_r': dc_cp(),
            're_exp_sf': dc_cp(),
            'alpha_ratio': dc_cp(),
            'area_ratio': dc_cp(),
            'UA_cecchinato': dc_gcp(
                elements=['re_exp_r', 're_exp_sf', 'alpha_ratio', 'area_ratio'],
                num_eq_sets=1,
                func=self.UA_cecchinato_func,
                dependents=self.UA_dependents,
            ),
            'td_pinch': dc_cp(
                min_val=0, num_eq_sets=1,
                func=self.td_pinch_func,
                dependents=self.td_pinch_dependents,
                quantity="temperature_difference"
            )
        })
        return params

    @staticmethod
    def _get_h_steps(c1, c2, num_sections=200):
        """Get the steps for enthalpy at the boundaries of phases during the
        change of enthalpy from one state to another

        Parameters
        ----------
        c1 : tespy.connections.connection.Connection
            Inlet connection.

        c2 : tespy.connections.connection.Connection
            Outlet connection.

        Returns
        -------
        list
            Steps of enthalpy of the specified connections
        """
        if c1.fluid.val != c2.fluid.val:
            msg = (
                "Both connections need to utilize the same fluid data: "
                f"{c1.fluid.val}, {c2.fluid.val}"
            )
            raise ValueError(msg)

        # change the order of connections to have c1 as the lower enthalpy
        # connection (enthalpy will be rising in the list)
        if c1.h.val_SI > c2.h.val_SI:
            c1, c2 = c2, c1

        h_at_steps = np.linspace(c1.h.val_SI, c2.h.val_SI, num_sections)

        return h_at_steps

    @staticmethod
    def _get_Q_sections(h_at_steps, mass_flow):
        """Calculate the heat exchange of every section given steps of
        enthalpy and mass flow.

        Parameters
        ----------
        h_at_steps : list
            Enthalpy values at sections (inlet, phase change points, outlet)
        mass_flow : float
            Mass flow value

        Returns
        -------
        float
            Heat exchanged between defined steps of enthalpy.
        """
        return np.diff(h_at_steps) * mass_flow

    def _assign_sections(self):
        """Assign the sections of the heat exchanger

        Returns
        -------
        list
            List of cumulative sum of heat exchanged defining the heat exchanger
            sections.
        """
        h_steps_hot = self._get_h_steps(self.inl[0], self.outl[0])
        Q_sections_hot = self._get_Q_sections(h_steps_hot, self.inl[0].m.val_SI)
        Q_sections_hot = np.insert(np.cumsum(Q_sections_hot)[:-1], 0, 0.0)

        return Q_sections_hot

    def _get_T_at_steps(self, Q_sections):
        """Calculate the temperature values for the provided sections.

        Parameters
        ----------
        Q_sections : list
            Cumulative heat exchanged from the hot side to the cold side
            defining the sections of the heat exchanger.

        Returns
        -------
        tuple
            Lists of cold side and hot side temperature
        """
        # now put the Q_sections back on the h_steps on both sides
        # Since Q_sections is defined increasing we have to start back from the
        # outlet of the hot side
        h_steps_hot = self.outl[0].h.val_SI + Q_sections / self.inl[0].m.val_SI
        h_steps_cold = self.inl[1].h.val_SI + Q_sections / self.inl[1].m.val_SI
        self.num_sections = 200
        p_steps_hot = np.linspace(self.outl[0].p.val_SI, self.inl[0].p.val_SI, self.num_sections)
        p_steps_cold = np.linspace(self.inl[1].p.val_SI, self.outl[1].p.val_SI, self.num_sections)
        T_steps_hot = np.array([
            T_mix_ph(p, h, self.inl[0].fluid_data, self.inl[0].mixing_rule)
            for p, h in zip(p_steps_hot, h_steps_hot)
        ])
        T_steps_cold = np.array([
            T_mix_ph(p, h, self.inl[1].fluid_data, self.inl[1].mixing_rule)
            for p, h in zip(p_steps_cold, h_steps_cold)
        ])
        return T_steps_hot, T_steps_cold

    @staticmethod
    def _calc_td_log_per_section(T_steps_hot, T_steps_cold, postprocess=False):
        """Calculate the logarithmic temperature difference values per section
        of heat exchanged.

        Parameters
        ----------
        T_steps_hot : list
            Temperature hot side at beginning and end of sections.

        T_steps_cold : list
            Temperature cold side at beginning and end of sections.

        Returns
        -------
        list
            Lists of temperature differences per section of heat exchanged.
        """
        if postprocess:
            td_at_steps = T_steps_hot - T_steps_cold
            if (td_at_steps <= 0).any():
                return np.ones(len(td_at_steps) - 1) * np.nan
        # the temperature ranges both come with increasing values
        td_at_steps = np.abs(T_steps_hot - T_steps_cold)

        return np.array([
            (td_at_steps[i + 1] - td_at_steps[i])
            / math.log(td_at_steps[i + 1] / td_at_steps[i])
            # round is required because tiny differences may cause
            # inconsistencies due to rounding errors
            if round(td_at_steps[i + 1], 6) != round(td_at_steps[i], 6)
            else td_at_steps[i + 1]
            for i in range(len(td_at_steps) - 1)
        ])

    def calc_sections(self, postprocess=True):
        """Calculate the sections of the heat exchanger.

        Returns
        -------
        tuple
            Cumulated heat transfer over sections, temperature at steps hot
            side, temperature at steps cold side, heat transfer per section
        """
        Q_sections = self._assign_sections()
        T_steps_hot, T_steps_cold = self._get_T_at_steps(Q_sections)
        Q_per_section = np.diff(Q_sections)
        td_log_per_section = self._calc_td_log_per_section(T_steps_hot, T_steps_cold, postprocess)
        return Q_sections, T_steps_hot, T_steps_cold, Q_per_section, td_log_per_section

    def calc_UA(self, sections):
        """Calculate the sum of UA for all sections in the heat exchanger

        Returns
        -------
        float
            Sum of UA values of all heat exchanger sections.
        """
        _, _, _, Q_per_section, td_log_per_section = sections
        UA_sections = Q_per_section / td_log_per_section
        return sum(UA_sections)

    def UA_func(self, **kwargs):
        r"""
        Residual method for fixed heat transfer coefficient UA.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = UA - \sum UA_\text{i}
        """
        sections = self.calc_sections(False)
        return self.UA.val_SI - self.calc_UA(sections)

    def UA_cecchinato_func(self):
        r"""
        Method to calculate heat transfer via UA design with modification
        for part load according to :cite:`cecchinato2010`. UA is determined
        over the UA values of the sections of the heat exchanger.

        .. note::

            You need to specify a couple of parameters to use this method. The
            values depend on the context, as they define relations between the
            refrigerant and the secondary fluid. For an evaporator the
            refrigerant is the cold side, for a condenser it is the hot side.
            You can check the linked publication for reference values.

        - alpha_ratio: Ratio of secondary fluid to refrigerant heat transfer
          coefficient
        - area_ratio: Ratio of secondary fluid to refrigerant area
        - re_exp_sf: Reynolds exponent for the secondary fluid mass flow
        - re_exp_r: Reynolds exponent for the refrigerant mass flow

        The modification factor for UA is calculated as follows

        .. math::

            f_\text{UA}=\frac{
                1 + \frac{\alpha_\text{sf}}{\alpha_\text{r}}
                \cdot\frac{A_\text{sf}}{A_\text{r}}
            }{
                \frac{\dot m_\text{sf}}{\dot m_\text{sf,ref}}^{-Re_\text{sf}} +
                \frac{\alpha_\text{sf}}{\alpha_\text{r}}
                \cdot\frac{A_\text{sf}}{A_\text{r}}
                \cdot\frac{\dot m_\text{r}}{\dot m_\text{r,ref}}^{-Re_\text{r}}
            }

        Returns
        -------
        float
            residual value of equation

            .. math::

                0 = UA_\text{ref} \cdot f_\text{UA} - \sum UA_\text{i}
        """
        alpha_ratio = self.alpha_ratio.val_SI
        area_ratio = self.area_ratio.val_SI
        re_exp_sf = self.re_exp_sf.val_SI
        re_exp_r = self.re_exp_r.val_SI
        refrigerant_index = self.refrigerant_index.val

        if refrigerant_index == 0:
            secondary_index = 1
        else:
            secondary_index = 0

        m_r = self.inl[refrigerant_index].m
        m_ratio_r = m_r.val_SI / m_r.design
        m_sf = self.inl[secondary_index].m
        m_ratio_sf = m_sf.val_SI / m_sf.design

        fUA = (
            (1 + alpha_ratio * area_ratio)
            / (
                m_ratio_sf ** -re_exp_sf
                + alpha_ratio * area_ratio * m_ratio_r ** -re_exp_r
            )
        )
        sections = self.calc_sections(False)
        return self.UA.design * fUA - self.calc_UA(sections)

    def UA_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].h,
            self.inl[1].m,
            self.inl[1].p,
            self.inl[1].h,
            self.outl[1].h
        ]

    def calc_td_pinch(self, sections):
        """Calculate the pinch point temperature difference

        Returns
        -------
        float
            Value of the pinch point temperature difference
        """
        _, T_steps_hot, T_steps_cold, _, _ = sections

        return min(T_steps_hot - T_steps_cold)

    def td_pinch_func(self):
        r"""
        Equation for minimal pinch temperature difference of sections.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = td_\text{pinch} - min(td_\text{i})
        """
        sections = self.calc_sections(False)
        return self.td_pinch.val_SI - self.calc_td_pinch(sections)

    def td_pinch_dependents(self):
        return [
            self.inl[0].m,
            self.inl[0].p,
            self.inl[0].h,
            self.outl[0].h,
            self.inl[1].m,
            self.inl[1].p,
            self.inl[1].h,
            self.outl[1].h
        ]

    def calc_parameters(self):
        super().calc_parameters()

        sections = self.calc_sections()
        self.UA.val_SI = self.calc_UA(sections)
        self.td_pinch.val_SI = self.calc_td_pinch(sections)
