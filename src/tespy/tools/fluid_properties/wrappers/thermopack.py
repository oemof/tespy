# -*- coding: utf-8

"""thermopack wrapper for fluid property calls.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/wrappers/thermopack.py

SPDX-License-Identifier: MIT
"""
from math import isfinite


from .base import FluidPropertyWrapper
from .base import wrapper_registry


@wrapper_registry
class ThermopackWrapper(FluidPropertyWrapper):

    def __init__(self, fluid, back_end=None, **kwargs) -> None:
        """Wrapper for the thermopack fluid property library

        Thermopack works on a molar basis, all inputs and outputs of this
        wrapper are mass based SI values in line with the other wrappers.
        Mixtures follow the same fluid string convention as the
        :code:`CoolPropWrapper`, e.g. :code:`"CO2[0.9]&N2[0.1]|molar"`, with
        component names from the thermopack component database. Intermediate
        vapor qualities in :code:`h_pQ`, :code:`h_QT`, :code:`s_QT`,
        :code:`d_QT` and :code:`p_sat_TQ` are linear interpolations between
        the bubble and dew state, which is exact for pure fluids only.
        Transport properties are not available in thermopack.

        Parameters
        ----------
        fluid : str
            Name of the fluid or mixture.
        back_end : str, optional
            Equation of state, by default "PR". Available are the cubic
            equations of state "SRK", "PR", "VdW", "RK", "SW" and "PT",
            the SAFT variants "PC-SAFT", "sPC-SAFT", "PCP-SAFT" and
            "SAFT-VR-MIE", the CPA variants "CPA-SRK" and "CPA-PR",
            Lee-Kesler "LK" as well as the multiparameter back ends
            "NIST_MEOS", "MBWR32" and "MBWR19".
        T_min, T_max, p_min, p_max : float, optional
            Override the temperature and pressure limits of the wrapper in
            SI units. Thermopack's default limits reach into regions where
            its flash routines are unstable, restricting them to the range
            of the application improves solver robustness. Pass them
            through the :code:`fluid_wrapper_kwargs` keyword on the
            connection.
        """
        # avoid unnecessary loading time if not used
        # a bare 'import thermopack' can succeed as an empty namespace
        # package, e.g. with a thermopack source checkout on the path
        try:
            import thermopack.thermo  # noqa: F401
        except ModuleNotFoundError:
            msg = (
                "To use the thermopack fluid properties you need to install "
                "thermopack, e.g. with 'uv add thermopack' or "
                "'pip install thermopack'."
            )
            raise ModuleNotFoundError(msg)

        if back_end is None:
            back_end = "PR"
        super().__init__(fluid, back_end)

        self._identify_mixture()
        self._create_AS()
        self._set_composition()
        self._set_constants(**kwargs)
        self._last_ph = None

    def _create_AS(self):
        comps = self.fluid.replace("&", ",")
        if self.back_end in ["SRK", "PR", "VdW", "RK", "SW", "PT"]:
            from thermopack.cubic import cubic
            self.AS = cubic(comps, self.back_end)
        elif self.back_end in ["PC-SAFT", "sPC-SAFT", "PCP-SAFT"]:
            from thermopack.pcsaft import pcsaft
            self.AS = pcsaft(
                comps,
                simplified=self.back_end == "sPC-SAFT",
                polar=self.back_end == "PCP-SAFT"
            )
        elif self.back_end == "SAFT-VR-MIE":
            from thermopack.saftvrmie import saftvrmie
            self.AS = saftvrmie(comps)
        elif self.back_end in ["CPA-SRK", "CPA-PR"]:
            from thermopack.cpa import cpa
            self.AS = cpa(comps, self.back_end.split("-")[1])
        elif self.back_end == "LK":
            from thermopack.lee_kesler import lee_kesler
            self.AS = lee_kesler(comps)
        elif self.back_end in ["NIST_MEOS", "MBWR32", "MBWR19"]:
            from thermopack.multiparameter import multiparam
            self.AS = multiparam(comps, self.back_end)
        else:
            msg = f"The specified back_end {self.back_end} is not available."
            raise NotImplementedError(msg)

    def _set_composition(self):
        if self.mixture_type == "volume":
            msg = (
                "Volume based composition is not available for the "
                f"{self.__class__.__name__}."
            )
            raise ValueError(msg)

        # compmoleweight returns g/mol in all thermopack versions
        self._component_molar_masses = [
            self.AS.compmoleweight(i + 1) * 1e-3
            for i in range(len(self.fluid.split("&")))
        ]

        if len(self.fractions) == 0:
            self._z = [1.0]
        elif self.mixture_type == "mass":
            molar = [
                w / M
                for w, M in zip(self.fractions, self._component_molar_masses)
            ]
            self._z = [n / sum(molar) for n in molar]
        else:
            self._z = [x / sum(self.fractions) for x in self.fractions]

    def _set_constants(self, **kwargs):
        # thermopack's default temperature minimum of 80 K lies inside a
        # region where the flash routines of many fluids stop the whole
        # interpreter instead of raising an error, restricting the
        # temperature and pressure ranges to the application of the model
        # protects the solver from stepping into that region
        if "T_min" in kwargs:
            self.AS.set_tmin(kwargs["T_min"])
        if "T_max" in kwargs:
            self.AS.set_tmax(kwargs["T_max"])
        if "p_min" in kwargs:
            self.AS.set_pmin(kwargs["p_min"])
        if "p_max" in kwargs:
            self.AS.set_pmax(kwargs["p_max"])
        self._T_min = self.AS.get_tmin()
        self._T_max = self.AS.get_tmax()
        self._p_min = self.AS.get_pmin()
        self._p_max = self.AS.get_pmax()
        self._molar_mass = sum(
            x * M for x, M in zip(self._z, self._component_molar_masses)
        )

        if len(self._z) == 1:
            self._T_crit = self.AS.critical_temperature(1)
            self._p_crit = self.AS.critical_pressure(1)
        else:
            try:
                self._T_crit, _, self._p_crit = self.AS.critical(self._z)
            except Exception:
                # the mixture critical point solver does not always converge,
                # the molar average keeps the phase region checks functional
                self._T_crit = sum(
                    x * self.AS.critical_temperature(i + 1)
                    for i, x in enumerate(self._z)
                )
                self._p_crit = sum(
                    x * self.AS.critical_pressure(i + 1)
                    for i, x in enumerate(self._z)
                )

    @staticmethod
    def _scalar(value):
        # thermopack v2 returns single element tuples for Tp-properties
        return value[0] if isinstance(value, tuple) else value

    def _single_phase(self, flash):
        # Tp-property methods require a concrete phase root, which the
        # flash does not deliver for e.g. supercritical states
        if flash.phase in [self.AS.LIQPH, self.AS.VAPPH]:
            return flash.phase
        return self.AS.guess_phase(flash.T, flash.p, flash.z)

    @staticmethod
    def _check_finite(*values):
        # thermopack's Fortran routines stop the whole interpreter on
        # non-finite inputs instead of raising an error
        if not all(map(isfinite, values)):
            msg = (
                "Thermopack property calls require finite input values, "
                f"got {values}."
            )
            raise ValueError(msg)

    def _check_flash_inputs(self, p, T=None):
        # out of bounds inputs, e.g. zero pressure during solver
        # initialization, stop the whole interpreter as well
        if not self._p_min <= p <= self._p_max:
            msg = (
                f"Pressure {p} Pa is out of the wrapper bounds "
                f"[{self._p_min}, {self._p_max}] Pa."
            )
            raise ValueError(msg)
        if T is not None and not self._T_min <= T <= self._T_max:
            msg = (
                f"Temperature {T} K is out of the wrapper bounds "
                f"[{self._T_min}, {self._T_max}] K."
            )
            raise ValueError(msg)

    def _check_sat_p(self, p):
        self._check_finite(p)
        if not self._p_min <= p <= self._p_crit:
            msg = (
                f"Saturation temperature requires a pressure between "
                f"{self._p_min} Pa and the critical pressure "
                f"{self._p_crit} Pa, got {p} Pa."
            )
            raise ValueError(msg)

    def _check_sat_T(self, T):
        self._check_finite(T)
        if not self._T_min <= T <= self._T_crit:
            msg = (
                f"Saturation pressure requires a temperature between "
                f"{self._T_min} K and the critical temperature "
                f"{self._T_crit} K, got {T} K."
            )
            raise ValueError(msg)

    def _update_ph(self, p, h):
        if self._last_ph != (p, h):
            self._check_finite(p, h)
            self._check_flash_inputs(p)
            self._flash = self.AS.two_phase_phflash(
                p, self._z, h * self._molar_mass
            )
            self._last_ph = (p, h)
        return self._flash

    def _tp_flash(self, p, T):
        self._check_finite(p, T)
        self._check_flash_inputs(p, T)
        return self.AS.two_phase_tpflash(T, p, self._z)

    def _ps_flash(self, p, s):
        self._check_finite(p, s)
        self._check_flash_inputs(p)
        return self.AS.two_phase_psflash(p, self._z, s * self._molar_mass)

    def _h_flash(self, flash):
        if flash.phase == self.AS.TWOPH:
            h = (
                flash.betaV * self._scalar(
                    self.AS.enthalpy(flash.T, flash.p, flash.y, self.AS.VAPPH)
                )
                + flash.betaL * self._scalar(
                    self.AS.enthalpy(flash.T, flash.p, flash.x, self.AS.LIQPH)
                )
            )
        else:
            h = self._scalar(
                self.AS.enthalpy(
                    flash.T, flash.p, flash.z, self._single_phase(flash)
                )
            )
        return h / self._molar_mass

    def _s_flash(self, flash):
        if flash.phase == self.AS.TWOPH:
            s = (
                flash.betaV * self._scalar(
                    self.AS.entropy(flash.T, flash.p, flash.y, self.AS.VAPPH)
                )
                + flash.betaL * self._scalar(
                    self.AS.entropy(flash.T, flash.p, flash.x, self.AS.LIQPH)
                )
            )
        else:
            s = self._scalar(
                self.AS.entropy(
                    flash.T, flash.p, flash.z, self._single_phase(flash)
                )
            )
        return s / self._molar_mass

    def _d_flash(self, flash):
        if flash.phase == self.AS.TWOPH:
            v = (
                flash.betaV * self._scalar(
                    self.AS.specific_volume(
                        flash.T, flash.p, flash.y, self.AS.VAPPH
                    )
                )
                + flash.betaL * self._scalar(
                    self.AS.specific_volume(
                        flash.T, flash.p, flash.x, self.AS.LIQPH
                    )
                )
            )
        else:
            v = self._scalar(
                self.AS.specific_volume(
                    flash.T, flash.p, flash.z, self._single_phase(flash)
                )
            )
        return self._molar_mass / v

    def _is_below_T_critical(self, T):
        return T < self._T_crit

    def isentropic(self, p_1, h_1, p_2):
        return self.h_ps(p_2, self.s_ph(p_1, h_1))

    def T_ph(self, p, h):
        return self._update_ph(p, h).T

    def T_ps(self, p, s):
        return self._ps_flash(p, s).T

    def h_pT(self, p, T):
        return self._h_flash(self._tp_flash(p, T))

    def h_ps(self, p, s):
        return self._h_flash(self._ps_flash(p, s))

    def h_pQ(self, p, Q):
        if Q == 0:
            h = self.AS.enthalpy(self.T_bubble(p), p, self._z, self.AS.LIQPH)
        elif Q == 1:
            h = self.AS.enthalpy(self.T_dew(p), p, self._z, self.AS.VAPPH)
        else:
            return (1 - Q) * self.h_pQ(p, 0) + Q * self.h_pQ(p, 1)
        return self._scalar(h) / self._molar_mass

    def h_QT(self, Q, T):
        if Q == 0:
            h = self.AS.enthalpy(T, self.p_bubble(T), self._z, self.AS.LIQPH)
        elif Q == 1:
            h = self.AS.enthalpy(T, self.p_dew(T), self._z, self.AS.VAPPH)
        else:
            return (1 - Q) * self.h_QT(0, T) + Q * self.h_QT(1, T)
        return self._scalar(h) / self._molar_mass

    def s_QT(self, Q, T):
        if Q == 0:
            s = self.AS.entropy(T, self.p_bubble(T), self._z, self.AS.LIQPH)
        elif Q == 1:
            s = self.AS.entropy(T, self.p_dew(T), self._z, self.AS.VAPPH)
        else:
            return (1 - Q) * self.s_QT(0, T) + Q * self.s_QT(1, T)
        return self._scalar(s) / self._molar_mass

    def d_QT(self, Q, T):
        if Q == 0:
            v = self.AS.specific_volume(
                T, self.p_bubble(T), self._z, self.AS.LIQPH
            )
        elif Q == 1:
            v = self.AS.specific_volume(
                T, self.p_dew(T), self._z, self.AS.VAPPH
            )
        else:
            return 1 / ((1 - Q) / self.d_QT(0, T) + Q / self.d_QT(1, T))
        return self._molar_mass / self._scalar(v)

    def T_sat(self, p):
        self._check_sat_p(p)
        return self.AS.bubble_temperature(p, self._z)[0]

    def T_dew(self, p):
        self._check_sat_p(p)
        return self.AS.dew_temperature(p, self._z)[0]

    def p_sat(self, T):
        self._check_sat_T(T)
        return self.AS.bubble_pressure(T, self._z)[0]

    def p_dew(self, T):
        self._check_sat_T(T)
        return self.AS.dew_pressure(T, self._z)[0]

    def p_sat_TQ(self, T, Q):
        if Q == 0:
            return self.p_bubble(T)
        elif Q == 1:
            return self.p_dew(T)
        return (1 - Q) * self.p_bubble(T) + Q * self.p_dew(T)

    def Q_ph(self, p, h):
        flash = self._update_ph(p, h)
        if flash.phase == self.AS.TWOPH:
            # betaV is the molar vapor fraction
            M_vap = sum(
                y * M for y, M in zip(flash.y, self._component_molar_masses)
            )
            M_liq = sum(
                x * M for x, M in zip(flash.x, self._component_molar_masses)
            )
            return (
                flash.betaV * M_vap
                / (flash.betaV * M_vap + flash.betaL * M_liq)
            )
        elif flash.phase == self.AS.LIQPH:
            return 0
        elif flash.phase == self.AS.VAPPH:
            return 1
        else:
            return -1

    def phase_ph(self, p, h):
        flash = self._update_ph(p, h)
        if flash.phase == self.AS.TWOPH:
            return "tp"
        elif flash.phase == self.AS.LIQPH:
            return "l"
        elif flash.phase == self.AS.VAPPH:
            return "g"
        else:
            return "sc"

    def d_ph(self, p, h):
        return self._d_flash(self._update_ph(p, h))

    def d_pT(self, p, T):
        return self._d_flash(self._tp_flash(p, T))

    def s_ph(self, p, h):
        return self._s_flash(self._update_ph(p, h))

    def s_pT(self, p, T):
        return self._s_flash(self._tp_flash(p, T))
