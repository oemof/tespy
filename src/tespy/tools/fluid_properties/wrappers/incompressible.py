# -*- coding: utf-8

"""Incompressible fluid wrapper based on tabular data.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/wrappers/incompressible.py

SPDX-License-Identifier: MIT
"""
import numpy as np
from scipy.optimize import brentq

from tespy.tools.fluid_properties.helpers import fit_incompressible_linear
from tespy.tools.fluid_properties.helpers import fit_incompressible_viscosity

from .base import FluidPropertyWrapper
from .base import wrapper_registry


@wrapper_registry
class IncompressibleFluidWrapper(FluidPropertyWrapper):
    """Class to represent a fluid in TESPy using tabular data

    Parameters
    ----------
    fluid : str
        Name of fluid
    back_end : str, optional
        Name of the back end in context of CoolProp, by default None
    temperature_data : np.ndarray
        Array of temperature measurements in SI units (Kelvin)
    density_data : np.ndarray
        Array of corresponding density values in SI units (kg/m3)
    heat_capacity_data : np.ndarray
        Array of corresponding heat capacity values in SI units (J/kg)
    viscosity_data : np.ndarray
        Array of corresponding **dynamic** viscosity values in SI units (Pas)
    conductivity_data : np.ndarray
        Array of corresponding thermal conductivity values in SI units
        (W/mK)
    """

    def __init__(self, fluid, back_end=None, **kwargs):
        """Class to represent a fluid in TESPy using tabular data

        Parameters
        ----------
        fluid : str
            Name of fluid
        back_end : str, optional
            Name of the back end in context of CoolProp, by default None
        temperature_data : np.ndarray
            Array of temperature measurements in SI units (Kelvin)
        density_data : np.ndarray
            Array of corresponding density values in SI units (kg/m3)
        heat_capacity_data : np.ndarray
            Array of corresponding heat capacity values in SI units (J/kg)
        viscosity_data : np.ndarray
            Array of corresponding **dynamic** viscosity values in SI units
            (Pas)
        conductivity_data : np.ndarray
            Array of corresponding thermal conductivity values in SI units
            (W/mK)
        """
        super().__init__(fluid, back_end, **kwargs)

        self.temperature_data = None
        self.heat_capacity_data = None
        self.density_data = None
        self.viscosity_data = None

        for key in ["temperature", "heat_capacity", "density", "viscosity"]:
            value = kwargs.get(f"{key}_data")
            if value is None:
                msg = (
                    f"The {self.__class__.__name__} requires specification of "
                    f"the '{key}_data' keyword in the form of a numpy array."
                )
                raise KeyError(msg)
            else:
                setattr(self, f"{key}_data", value)

        self.conductivity_data = kwargs.get("conductivity_data")

        self._T_ref = kwargs.get("T_ref", min(self.temperature_data))
        self._p_ref = kwargs.get("p_ref", 1e5)

        self._fit_data()
        self._set_constants()

    def _fit_data(self):
        A, B = fit_incompressible_linear(
            self.temperature_data, self.heat_capacity_data
        )
        self._heat_capacity = {
            "A": A,
            "B": B
        }

        A, B = fit_incompressible_linear(
            self.temperature_data, self.density_data
        )
        self._density = {
            "A": A,
            "B": B
        }

        if self.conductivity_data is not None:
            A, B = fit_incompressible_linear(
                self.temperature_data, self.conductivity_data
            )
        else:
            A, B = np.nan, np.nan

        self._conductivity = {
            "A": A,
            "B": B
        }

        A, B, C, D = fit_incompressible_viscosity(
            self.temperature_data, self.viscosity_data
        )
        self._viscosity = {
            "A": A,
            "B": B,
            "C": C,
            "D": D
        }

    def _set_constants(self):
        # evaluate h at T=T_ref
        self._h_ref = self._h_pT(None, self._T_ref)

        self._T_min = self._T_ref
        self._T_max = max(self.temperature_data)

        self._molar_mass = 1
        self._p_min = 100
        self._p_max = 10000000
        self._p_crit = self._p_max

        self._T_crit = None

    def get_fitting_report(self):
        import matplotlib.pyplot as plt

        def plot_property(ax, temperature, measurements, evaluation):

            _fit, = ax.plot(temperature, evaluation, "-", color="red")
            _data = ax.scatter(temperature, measurements, marker="x", c="blue")

            ax_err = ax.twinx()

            _err = ax_err.scatter(
                temperature, (evaluation - measurements) / measurements * 100,
                c="#0000ff66"
            )

            ax_err.set_ylabel("Deviation between fit and data in %")

            return [_data, _fit, _err]

        fig, ax = plt.subplots(2, 2, figsize=(10, 10), sharex=True)


        ax[0, 0].set_title("Heat capacity")
        ax[0, 1].set_title("Density")
        ax[1, 0].set_title("Viscosity")
        ax[1, 1].set_title("Thermal conductivity")

        temperature_data = self.temperature_data
        heat_capacity_data = self.heat_capacity_data
        density_data = self.density_data
        viscosity_data = self.viscosity_data
        conductivity_data = self.conductivity_data

        d = 0.001
        heat_capacity_eval = (
            self.h_pT(None, temperature_data + d)
            - self.h_pT(None, temperature_data - d)
        ) / (2 * d)

        density_eval = self.d_pT(None, temperature_data)
        viscosity_eval = self.viscosity_pT(None, temperature_data)
        conductivity_eval = self.conductivity_pT(None, temperature_data)

        lines = plot_property(ax[0, 0], temperature_data, heat_capacity_data, heat_capacity_eval)
        labels = ["datapoints", "fitted function", "deviation"]

        plot_property(ax[0, 1], temperature_data, density_data, density_eval)

        plot_property(ax[1, 0], temperature_data, viscosity_data, viscosity_eval)
        ax[1, 0].set_yscale("log")

        plot_property(ax[1, 1], temperature_data, conductivity_data, conductivity_eval)


        ax[0, 0].set_ylabel("Heat capacity in J/kgK")
        ax[0, 1].set_ylabel("Density in kg/m3")
        ax[1, 0].set_ylabel("Viscosity in Pas")
        ax[1, 1].set_ylabel("Thermal conductivity in W/mK")

        ax[1, 0].set_xlabel("Temperature in K")
        ax[1, 1].set_xlabel("Temperature in K")

        fig.legend(
            lines, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.05)
        )

        plt.tight_layout()

        return fig, ax

    def T_ph(self, p, h):
        # Inverse function of h_pT, using the quadratic formula in its
        # numerically stable conjugate form: the standard form divides by
        # the polynomial slope A, which suffers catastrophic cancellation
        # for (near-)constant heat capacity data with A close to zero
        return (
            2 * (h + self._h_ref)
            / (
                self._heat_capacity["B"]
                + (
                    self._heat_capacity["B"] ** 2
                    + 2 * self._heat_capacity["A"] * (h + self._h_ref)
                ) ** 0.5
            )
        )

    def h_pT(self, p, T):
        return self._h_pT(p, T) - self._h_ref

    def _h_pT(self, p, T):
        # h = integral cp(T) dT
        return (
            0.5 * self._heat_capacity["A"] * T ** 2
            + self._heat_capacity["B"] * T
        )

    def h_ps(self, p, s):
        return self.h_pT(p, self.T_ps(p, s))

    def s_ph(self, p, h):
        return self.s_pT(p, self.T_ph(p, h))

    def s_pT(self, p, T):
        # s0 = 0
        return (
            self._heat_capacity["B"] * np.log(T / self._T_ref)
            + self._heat_capacity["A"] * (T - self._T_ref)
            - self.d_pT(p, T) * (p - self._p_ref)
        )

    def isentropic(self, p_1, h_1, p_2):
        # assumption that temperature barely changes
        T = self.T_ph(p_1, h_1)
        return h_1 + (p_2 - p_1) / self.d_pT(p_1, T)

    def _inverse_s_pT(self, T, p, s):
        return s - self.s_pT(p, T)

    def T_ps(self, p, s):
        return brentq(
            self._inverse_s_pT,
            self._T_min,
            self._T_max,
            args=(p, s)
        )

    def conductivity_ph(self, p, h):
        return self.conductivity_pT(p, self.T_ph(p, h))

    def conductivity_pT(self, p, T):
        return self._conductivity["A"] * T + self._conductivity["B"]

    def d_ph(self, p, h):
        return self.d_pT(p, self.T_ph(p, h))

    def d_pT(self, p, T):
        return self._density["A"] * T + self._density["B"]

    def phase_ph(self, p, h):
        return "l"

    def viscosity_ph(self, p, h):
        return self.viscosity_pT(p, self.T_ph(p, h))

    def viscosity_pT(self, p, T):
        return np.exp(
            self._viscosity["A"] /  T ** 3
            + self._viscosity["B"] /  T ** 2
            + self._viscosity["C"] /  T
            + self._viscosity["D"]
        )
