import pint


class Units:

    def __init__(self):
        self.default = {
            "temperature": "kelvin",
            "temperature_difference": "delta_degC",
            "enthalpy": "J/kg",
            "entropy": "J/kg/K",
            "pressure": "Pa",
            "mass_flow": "kg/s",
            "volumetric_flow": "m3/s",
            "specific_volume": "m3/kg",
            "power": "W",
            "heat": "W",
            "quality": "1",
            "vapor_mass_fraction": "1",  # for backwards compatibility
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
        self._ureg = pint.UnitRegistry()
        # m3 is not in standard ureg
        self.ureg.define("m3 = m ** 3")
        self.ureg.define("m2 = m ** 2")
        self._quantities = {
            k: self.ureg.Quantity(1, v) for k, v in self.default.items()
        }

    def set_default_units(self, **kwargs):
        for key, value in kwargs.items():
            if self._quantities[key].is_compatible_with(value):
                self.default[key] = value
            else:
                msg = f"Unit {value} is not compatible with quantity {key}"
                raise ValueError(msg)

    def set_ureg(self, ureg):
        self._ureg = ureg
        self.ureg.define("m3 = m ** 3")
        self.ureg.define("m2 = m ** 2")
        self._quantities = {
            k: self.ureg.Quantity(1, v) for k, v in self.default.items()
        }

    def get_ureg(self):
        return self._ureg

    ureg = property(get_ureg, set_ureg)


_UNITS = Units()
