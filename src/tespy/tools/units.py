import pint


class DefaultUnits:

    def __init__(self):
        self.default = {
            "temperature": "degC",
            "enthalpy": "kJ/kg",
            "entropy": "kJ/kg/K",
            "pressure": "bar",
            "mass_flow": "kg/s",
            # "volumetric_flow": "m3/s",  # not in ureg
            "power": "MW",
            "heat": "MW",
            "quality": "1",
            "efficiency": "1",
            "ratio": "1"
        }
        self.ureg = pint.UnitRegistry()
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
