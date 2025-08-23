import pint


class DefaultUnits:

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
            "efficiency": "1",
            "ratio": "1",
            # None is the default if not quantity is supplied
            None: "1"
        }
        # cannot use the setter here because we have to define m3 first!
        self._ureg = pint.UnitRegistry()
        # m3 is not in standard ureg
        self.ureg.define("m3 = meter ** 3")
        self._quantities = {
            k: self.ureg.Quantity(1, v) for k, v in self.default.items()
        }

    def set_default_units(self, **kwargs):
        for key, value in kwargs.items():
            if self._quantities[key].is_compatible_with(value):
                self.default[key] = value

    def set_ureg(self, ureg):
        self._ureg = ureg
        self._quantities = {
            k: self.ureg.Quantity(1, v) for k, v in self.default.items()
        }


    def get_ureg(self):
        return self._ureg

    ureg = property(get_ureg, set_ureg)


UNITS = DefaultUnits()
