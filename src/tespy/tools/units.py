import pint


UREG = pint.UnitRegistry()


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
        self._quantities = {
            k: UREG.Quantity(1, v) for k, v in self.default.items()
        }

    def set_default_units(self, **kwargs):
        for key, value in kwargs.items():
            if self._quantities[key].is_compatible_with(value):
                self.default[key] = value


UNITS = DefaultUnits()
