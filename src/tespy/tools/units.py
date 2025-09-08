# -*- coding: utf-8

"""Module for Units class.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/units.py

SPDX-License-Identifier: MIT
"""
import warnings

import pint


class Units:

    @classmethod
    def from_json(cls, default_units):
        instance = cls()
        instance.set_defaults(**default_units)
        return instance

    def __init__(self):
        self.default = {
            "temperature": "kelvin",
            "temperature_difference": "delta_degC",
            "enthalpy": "J/kg",
            "specific_energy": "J/kg",
            "entropy": "J/kg/K",
            "pressure": "Pa",
            "mass_flow": "kg/s",
            "volumetric_flow": "m3/s",
            "specific_volume": "m3/kg",
            "power": "W",
            "heat": "W",
            "quality": "1",
            "vapor_mass_fraction": "1",  # backwards compatibility network import
            "efficiency": "1",
            "ratio": "1",
            "length": "m",
            "speed": "m/s",
            "area": "m2",
            "thermal_conductivity": "W/m/K",
            "heat_transfer_coefficient": "W/K",
            # None is the default if not quantity is supplied
            None: "1"
        }
        # cannot use the setter here because we have to define m3 first!
        self._ureg = pint.UnitRegistry(cache_folder=":auto:")
        # m3 is not in standard ureg
        self.ureg.define("m3 = m ** 3")
        self.ureg.define("m2 = m ** 2")
        self.ureg.define("kgK = kg * K")
        self._quantities = {
            k: self.ureg.Quantity(1, v) for k, v in self.default.items()
        }

    def set_defaults(self, **kwargs):
        """Set the default units

        Parameters
        ----------
        temperature : str
            Default unit: "kelvin"
        temperature_difference : str
            Default unit: "delta_degC"
        enthalpy : str
            Default unit: "J/kg"
        specific_energy : str
            Default unit: "J/kg"
        entropy : str
            Default unit: "J/kg/K"
        pressure : str
            Default unit: "Pa"
        mass_flow : str
            Default unit: "kg/s"
        volumetric_flow : str
            Default unit: "m3/s"
        specific_volume : str
            Default unit: "m3/kg"
        power : str
            Default unit: "W"
        heat : str
            Default unit: "W"
        quality : str
            Default unit: "1"
        efficiency : str
            Default unit: "1"
        ratio : str
            Default unit: "1"
        length : str
            Default unit: "m"
        speed : str
            Default unit: "m/s"
        area : str
            Default unit: "m2"
        thermal_conductivity : str
            Default unit: "W/m/K"
        heat_transfer_coefficient : str
            Default unit: "W/K"
        """
        for key, value in kwargs.items():
            self._check_quantity_exists(key)
            if value == "-":
                value = "1"
            elif value == "C":
                value = "degC"
                msg = (
                    "The unit 'C' is used for 'Coulomb' in pint. For "
                    "backwards compatibility it will be parsed as degC for "
                    "now. Please use '°C' (or correct pint aliases) instead  "
                    "as it will stop working with the next major release"
                )
                warnings.warn(msg, FutureWarning)
            if self._is_compatible(key, value):
                self.default[key] = value
            else:
                msg = f"Unit {value} is not compatible with quantity {key}"
                raise ValueError(msg)

    def _is_compatible(self, quantity, unit):
        if quantity == "temperature_difference":
            if unit.startswith("delta_"):
                return self._quantities[quantity].is_compatible_with(unit)
            else:
                _units = self.ureg._units
                kelvin = list(_units["K"].aliases) + [_units["K"].name]
                rankine = list(_units["rankine"].aliases) + [_units["rankine"].name]
                return unit in kelvin or unit in rankine
        else:
            return self._quantities[quantity].is_compatible_with(unit)

    def get_default(self, quantity):
        self._check_quantity_exists(quantity)
        return self.default[quantity]

    def _check_quantity_exists(self, quantity):
        if quantity not in self.default:
            msg = (
                f"The quantity {quantity} is unknown. Please specify any of "
                f"the following: "
                f"{', '.join([key for key in self.default if key is not None])}."
            )
            raise KeyError(msg)


    def set_ureg(self, ureg):
        """Replace default ureg with a custom one

        Parameters
        ----------
        ureg : pint.UnitRegistry
            Change the pint.UnitRegistry to a custom one
        """
        self._ureg = ureg
        self.ureg.define("m3 = m ** 3")
        self.ureg.define("m2 = m ** 2")
        self._quantities = {
            k: self.ureg.Quantity(1, v) for k, v in self.default.items()
        }

    def get_ureg(self):
        return self._ureg

    def _serialize(self):
        return {k: v for k, v in self.default.items() if k is not None}

    ureg = property(get_ureg, set_ureg)


_UNITS = Units()
SI_UNITS = _UNITS.default.copy()
