from CoolProp.CoolProp import PropsSI
from CoolProp.CoolProp import AbstractState
import CoolProp as CP
from tespy.tools.helpers import newton

PRECISION = 1e-6


def _is_larger_than_precision(value):
    return value > PRECISION


def T_mix_ph(p, h, fluid_data, mixing_rule=None, T0=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["property_object"].T_ph(p, h)
    else:
        if mixing_rule == "ideal":
            return T_mix_ph_ideal(p, h, fluid_data, T0=T0)
        elif mixing_rule == "ideal-cond":
            return T_mix_ph_ideal_cond(p, h, fluid_data, T0=T0)
        elif mixing_rule == "incompressible":
            return T_mix_ph_incompressible(p, h, fluid_data, T0=T0)
        else:
            raise ValueError()


def h_mix_pT(p, T, fluid_data, mixing_rule=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["property_object"].h_pT(p, T)
    else:
        if mixing_rule == "ideal":
            return h_mix_pT_ideal(p, T, fluid_data)
        elif mixing_rule == "ideal-cond":
            return h_mix_pT_ideal_cond(p, T, fluid_data)
        elif mixing_rule == "incompressible":
            return h_mix_pT_incompressible(p, T, fluid_data)
        else:
            raise ValueError()


def T_mix_ps(p, s, fluid_data, mixing_rule=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["property_object"].T_ps(p, s)
    else:
        if mixing_rule == "ideal":
            return T_mix_ps_ideal(p, s, fluid_data)
        elif mixing_rule == "ideal-cond":
            return T_mix_ps_ideal_cond(p, s, fluid_data)
        elif mixing_rule == "incompressible":
            return T_mix_ps_incompressible(p, s, fluid_data)
        else:
            raise ValueError()


def v_mix_ph(p, h, fluid_data, mixing_rule=None, T0=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return 1 / pure_fluid["property_object"].d_ph(p, h)
    else:
        if mixing_rule == "ideal":
            return v_mix_pT_ideal(p, T_mix_ph(p, h, fluid_data, mixing_rule), fluid_data, T0=T0)
        elif mixing_rule == "ideal-cond":
            return v_mix_pT_ideal_cond(p, T_mix_ph(p, h, fluid_data, mixing_rule), fluid_data, T0=T0)
        elif mixing_rule == "incompressible":
            raise NotImplementedError()
        else:
            raise ValueError()


def viscosity_mix_ph(p, h, fluid_data, mixing_rule=None, T0=None):
    if get_number_of_fluids(fluid_data) == 1:
        pure_fluid = get_pure_fluid(fluid_data)
        return pure_fluid["property_object"].viscosity_ph(p, h)
    else:
        if mixing_rule == "ideal" or mixing_rule == "ideal-cond":
            return viscosity_mix_pT_ideal(p, T_mix_ph(p, h, fluid_data, mixing_rule), fluid_data, T0=T0)
        elif mixing_rule == "incompressible":
            raise NotImplementedError()
        else:
            raise ValueError()


def get_number_of_fluids(fluid_data):
    return sum([1 for f in fluid_data.values() if _is_larger_than_precision(f["mass_fraction"])])


def get_pure_fluid(fluid_data):
    for f in fluid_data.values():
        if _is_larger_than_precision(f["mass_fraction"]):
            return f


def get_molar_fractions(fluid_data):
    molarflow = {
        key: value["mass_fraction"] / value["property_object"]._molar_mass
        for key, value in fluid_data.items()
    }
    molarflow_sum = sum(molarflow.values())
    return {key: value / molarflow_sum for key, value in molarflow.items()}


def h_mix_pT_ideal(p=None, T=None, fluid_data=None, **kwargs):
    molar_fractions = get_molar_fractions(fluid_data)

    h = 0
    for fluid, data in fluid_data.items():

        if _is_larger_than_precision(data["mass_fraction"]):
            pp = p * molar_fractions[fluid]
            h += data["property_object"].h_pT(pp, T) * data["mass_fraction"]

    return h


def h_mix_pT_ideal_cond(p=None, T=None, fluid_data=None, **kwargs):

    water_alias = _water_in_mixture(fluid_data)
    if water_alias:
        water_alias = next(iter(water_alias))
        mass_fractions_gas, molar_fraction_gas, mass_liquid, molar_liquid = cond_check(p, T, fluid_data, water_alias)
        if mass_liquid == 0:
            return h_mix_pT_ideal(p, T, fluid_data, **kwargs)
        h = 0
        for fluid, data in fluid_data.items():
            if _is_larger_than_precision(data["mass_fraction"]):
                if fluid == water_alias:
                    h += fluid_data[water_alias]["property_object"].h_QT(0, T) * mass_liquid
                    h += fluid_data[water_alias]["property_object"].h_QT(1, T) * mass_fractions_gas[fluid] * (1 - mass_liquid)
                else:
                    pp = p * molar_fraction_gas[fluid]
                    h += data["property_object"].h_pT(pp, T) * mass_fractions_gas[fluid] * (1 - mass_liquid)
        return h
    else:
        return h_mix_pT_ideal(p, T, fluid_data, **kwargs)


def h_mix_pT_incompressible(p, T, fluid_data, **kwargs):

    h = 0
    for data in fluid_data.values():
        if _is_larger_than_precision(data["mass_fraction"]):
            h += data["property_object"].h_pT(p, T) * data["mass_fraction"]

    return h



def s_mix_pT_ideal(p=None, T=None, fluid_data=None, **kwargs):
    molar_fractions = get_molar_fractions(fluid_data)

    s = 0
    for fluid, data in fluid_data.items():

        if _is_larger_than_precision(data["mass_fraction"]):
            pp = p * molar_fractions[fluid]
            s += data["property_object"].s_pT(pp, T) * data["mass_fraction"]

    return s


def s_mix_pT_ideal_cond(p=None, T=None, fluid_data=None, **kwargs):

    water_alias = _water_in_mixture(fluid_data)
    if water_alias:
        water_alias = next(iter(water_alias))
        mass_fractions_gas, molar_fraction_gas, mass_liquid, molar_liquid = cond_check(p, T, fluid_data, water_alias)
        if mass_liquid == 0:
            return s_mix_pT_ideal(p, T, fluid_data, **kwargs)
        s = 0
        for fluid, data in fluid_data.items():
            if _is_larger_than_precision(data["mass_fraction"]):
                if fluid == water_alias:
                    s += fluid_data[water_alias]["property_object"].s_QT(0, T) * mass_liquid
                    s += fluid_data[water_alias]["property_object"].s_QT(1, T) * mass_fractions_gas[fluid] * (1 - mass_liquid)
                else:
                    pp = p * molar_fraction_gas[fluid]
                    s += data["property_object"].s_pT(pp, T) * mass_fractions_gas[fluid] * (1 - mass_liquid)
        return s
    else:
        return s_mix_pT_ideal(p, T, fluid_data, **kwargs)


def s_mix_pT_incompressible(p=None, T=None, fluid_data=None, **kwargs):

    s = 0
    for data in fluid_data.values():

        if _is_larger_than_precision(data["mass_fraction"]):
            s += data["property_object"].s_pT(p, T) * data["mass_fraction"]

    return s


def v_mix_pT_ideal(p=None, T=None, fluid_data=None, **kwargs):
    molar_fractions = get_molar_fractions(fluid_data)

    d = 0
    for fluid, data in fluid_data.items():

        if _is_larger_than_precision(data["mass_fraction"]):
            pp = p * molar_fractions[fluid]
            d += data["property_object"].d_pT(pp, T)

    return 1 / d


def v_mix_pT_ideal_cond(p=None, T=None, fluid_data=None, **kwargs):

    water_alias = _water_in_mixture(fluid_data)
    if water_alias:
        water_alias = next(iter(water_alias))
        mass_fractions_gas, molar_fraction_gas, mass_liquid, molar_liquid = cond_check(p, T, fluid_data, water_alias)
        if mass_liquid == 0:
            return v_mix_pT_ideal(p, T, fluid_data, **kwargs)
        d = 0
        for fluid, data in fluid_data.items():
            if _is_larger_than_precision(data["mass_fraction"]):
                if fluid == water_alias:
                    d += fluid_data[water_alias]["property_object"].d_QT(0, T) * mass_liquid
                    d += fluid_data[water_alias]["property_object"].d_QT(1, T) * (1 - mass_liquid)
                else:
                    pp = p * molar_fraction_gas[fluid]
                    d += data["property_object"].d_pT(pp, T) * (1 - mass_liquid)
        return 1 / d
    else:
        return v_mix_pT_ideal(p, T, fluid_data, **kwargs)


def v_mix_pT_incompressible(p=None, T=None, fluid_data=None, **kwargs):

    v = 0
    for data in fluid_data.values():
        if _is_larger_than_precision(data["mass_fraction"]):
            v += 1 / data["property_object"].d_pT(p, T) * data["mass_fraction"]

    return v


def viscosity_mix_pT_ideal(p=None, T=None, fluid_data=None, **kwargs):
    r"""
    Calculate dynamic viscosity from pressure and temperature.

    Parameters
    ----------
    flow : list
        Fluid property vector containing mass flow, pressure, enthalpy and
        fluid composition.

    T : float
        Temperature T / K.

    Returns
    -------
    visc : float
        Dynamic viscosity visc / Pa s.

    Note
    ----
    Calculation for fluid mixtures.

    .. math::

        \eta_{mix}(p,T)=\frac{\sum_{i} \left( \eta(p,T,fluid_{i}) \cdot y_{i}
        \cdot \sqrt{M_{i}} \right)}
        {\sum_{i} \left(y_{i} \cdot \sqrt{M_{i}} \right)}\;
        \forall i \in \text{fluid components}\\
        y: \text{volume fraction}\\
        M: \text{molar mass}

    Reference: :cite:`Herning1936`.
    """
    molar_fractions = get_molar_fractions(fluid_data)

    a = 0
    b = 0
    for fluid, data in fluid_data.items():
        if _is_larger_than_precision(data["mass_fraction"]):
            bi = molar_fractions[fluid] * data["property_object"]._molar_mass ** 0.5
            b += bi
            a += bi * data["property_object"].viscosity_pT(p, T)

    return a / b


def viscosity_mix_pT_incompressible(p=None, T=None, fluid_data=None, **kwargs):

    viscosity = 0
    for data in fluid_data.values():
        if _is_larger_than_precision(data["mass_fraction"]):
            viscosity += data["property_object"].viscosity_pT(p, T) * data["mass_fraction"]

    return viscosity


def _water_in_mixture(fluid_data):
    water_aliases = set(CP.CoolProp.get_aliases("H2O"))
    return water_aliases & set([f for f in fluid_data if _is_larger_than_precision(fluid_data[f]["mass_fraction"])])


def T_mix_ph_ideal(p=None, h=None, fluid_data=None, T0=None):
    # calculate the fluid properties for fluid mixtures
    valmin, valmax = get_mixture_temperature_range(fluid_data)
    if T0 is None:
        T0 = (valmin + valmax) / 2

    function_kwargs = {
        "p": p, "fluid_data": fluid_data, "T": T0,
        "function": h_mix_pT_ideal, "parameter": "T" , "delta": 0.01
    }
    return newton_with_kwargs(
        central_difference,
        h,
        val0=T0,
        valmin=valmin,
        valmax=valmax,
        **function_kwargs
    )


def T_mix_ph_ideal_cond(p=None, h=None, fluid_data=None, T0=None):
    # calculate the fluid properties for fluid mixtures
    valmin, valmax = get_mixture_temperature_range(fluid_data)
    if T0 is None:
        T0 = (valmin + valmax) / 2

    function_kwargs = {
        "p": p, "fluid_data": fluid_data, "T": T0,
        "function": h_mix_pT_ideal_cond, "parameter": "T" , "delta": 0.01
    }
    return newton_with_kwargs(
        central_difference,
        h,
        val0=T0,
        valmin=valmin,
        valmax=valmax,
        tol_rel=1e-4,
        **function_kwargs
    )


def T_mix_ph_incompressible(p=None, h=None, fluid_data=None, T0=None):
    # calculate the fluid properties for fluid mixtures
    valmin, valmax = get_mixture_temperature_range(fluid_data)
    if T0 is None:
        T0 = (valmin + valmax) / 2

    function_kwargs = {
        "p": p, "fluid_data": fluid_data, "T": T0,
        "function": h_mix_pT_incompressible, "parameter": "T" , "delta": 0.01
    }
    return newton_with_kwargs(
        central_difference,
        s,
        val0=T0,
        valmin=valmin,
        valmax=valmax,
        **function_kwargs
    )


def T_mix_ps_ideal(p=None, s=None, fluid_data=None, T0=None):
    # calculate the fluid properties for fluid mixtures
    valmin, valmax = get_mixture_temperature_range(fluid_data)
    if T0 is None:
        T0 = (valmin + valmax) / 2

    function_kwargs = {
        "p": p, "fluid_data": fluid_data, "T": T0,
        "function": s_mix_pT_ideal, "parameter": "T" , "delta": 0.01
    }
    return newton_with_kwargs(
        central_difference,
        s,
        val0=T0,
        valmin=valmin,
        valmax=valmax,
        **function_kwargs
    )


def T_mix_ps_ideal_cond(p=None, s=None, fluid_data=None, T0=None):
    # calculate the fluid properties for fluid mixtures
    valmin, valmax = get_mixture_temperature_range(fluid_data)
    if T0 is None:
        T0 = (valmin + valmax) / 2

    function_kwargs = {
        "p": p, "fluid_data": fluid_data, "T": T0,
        "function": s_mix_pT_ideal_cond, "parameter": "T" , "delta": 0.01
    }
    return newton_with_kwargs(
        central_difference,
        s,
        val0=T0,
        valmin=valmin,
        valmax=valmax,
        tol_rel=1e-4,
        **function_kwargs
    )


def T_mix_ps_incompressible(p=None, s=None, fluid_data=None, T0=None):
    # calculate the fluid properties for fluid mixtures
    valmin, valmax = get_mixture_temperature_range(fluid_data)
    if T0 is None:
        T0 = (valmin + valmax) / 2

    function_kwargs = {
        "p": p, "fluid_data": fluid_data, "T": T0,
        "function": s_mix_pT_incompressible, "parameter": "T" , "delta": 0.01
    }
    return newton_with_kwargs(
        central_difference,
        s,
        val0=T0,
        valmin=valmin,
        valmax=valmax,
        **function_kwargs
    )


def get_mixture_temperature_range(fluid_data):
    valmin = max(
        [v["property_object"]._T_min for v in fluid_data.values() if _is_larger_than_precision(v["mass_fraction"])]
    ) + 0.1
    valmax = min(
        [v["property_object"]._T_max for v in fluid_data.values() if _is_larger_than_precision(v["mass_fraction"])]
    ) - 0.1
    return valmin, valmax


def newton_with_kwargs(derivative, target_value, val0=300, valmin=70, valmax=3000, max_iter=10, tol_rel=PRECISION, tol_abs=PRECISION, tol_mode="rel", **function_kwargs):

    # start newton loop
    iteration = 0
    expr = True
    x = val0
    parameter = function_kwargs["parameter"]
    function = function_kwargs["function"]

    while expr:
        # calculate function residual and new value
        function_kwargs[parameter] = x
        residual = target_value - function(**function_kwargs)
        x += residual / derivative(**function_kwargs)

        # check for value ranges
        if x < valmin:
            x = valmin
        if x > valmax:
            x = valmax
        iteration += 1

        if iteration > max_iter:
            msg = (
                'The Newton algorithm was not able to find a feasible value '
                f'for function {function}. Current value with x={x} is '
                f'{function(**function_kwargs)}, target value is '
                f'{target_value}, residual is {residual} after {iteration} '
                'iterations.'
            )
            print(msg)
            # logging.debug(msg)

            break
        if tol_mode == 'abs':
            expr = abs(residual) >= tol_abs
        elif tol_mode == 'rel':
            expr = abs(residual / target_value) >= tol_rel
        else:
            expr = abs(residual / target_value) >= tol_rel or abs(residual) >= tol_abs

    return x


def central_difference(function=None, parameter=None, delta=None, **kwargs):
    upper = kwargs.copy()
    upper[parameter] += delta
    lower = kwargs
    lower[parameter] -= delta
    return (function(**upper) - function(**lower)) / (2 * delta)


def cond_check(p, T, fluid_data, water_alias):
    """Check if water is partially condensing in gaseous mixture.

    Parameters
    ----------
    y_i : dict
        Mass specific fluid composition.
    x_i : dict
        Mole specific fluid composition.
    p : float
        Pressure of mass flow.
    n : float
        Molar mass flow.
    T : float
        Temperature of mass flow.

    Returns
    -------
    tuple
        Tuple containing gas phase mass specific and molar specific
        compositions and overall liquid water mass fraction.
    """
    molar_fractions = get_molar_fractions(fluid_data)
    molar_fractions_gas = molar_fractions
    mass_fractions_gas = {f: v["mass_fraction"] for f, v in fluid_data.items()}
    water_mass_liquid = 0
    water_molar_liquid = 0

    if fluid_data[water_alias]["property_object"]._is_below_T_critical(T):
        p_sat = fluid_data[water_alias]["property_object"].p_boiling(T)
        pp_water = p * molar_fractions[water_alias]

        if p_sat < pp_water:
            water_molar_gas = (1 - molar_fractions[water_alias]) / (p / p_sat - 1)
            water_molar_liquid = molar_fractions[water_alias] - water_molar_gas
            x_gas_sum = 1 - water_molar_liquid

            molar_fractions_gas = {f: x / x_gas_sum for f, x in molar_fractions.items()}
            molar_fractions_gas[water_alias] = water_molar_gas / x_gas_sum

            water_mass_liquid = (
                water_molar_liquid
                * fluid_data[water_alias]["property_object"]._molar_mass
                / calc_molar_mass_mixture(fluid_data, molar_fractions)
            )

            molar_mass_mixture = calc_molar_mass_mixture(fluid_data, molar_fractions_gas)
            mass_fractions_gas = {
                fluid: (
                    x / molar_mass_mixture
                    * fluid_data[fluid]["property_object"]._molar_mass
                )
                for fluid, x in molar_fractions_gas.items()
            }

    return mass_fractions_gas, molar_fractions_gas, water_mass_liquid, water_molar_liquid


def calc_molar_mass_mixture(fluid_data, molar_fractions):
    return sum([x * fluid_data[fluid]["property_object"]._molar_mass for fluid, x in molar_fractions.items()])


class CoolPropWrapper:

    def __init__(self, fluid, backend=None) -> None:
        self.fluid = fluid
        if backend is None:
            backend = "HEOS"

        self.AS = AbstractState(backend, fluid)
        self._set_constants()

    def _set_constants(self):
        self._p_crit = self.AS.trivial_keyed_output(CP.iP_critical)
        self._T_crit = self.AS.trivial_keyed_output(CP.iT_critical)
        self._p_min = self.AS.trivial_keyed_output(CP.iP_min)
        self._p_max = self.AS.trivial_keyed_output(CP.iP_max)
        self._T_min = self.AS.trivial_keyed_output(CP.iT_min)
        self._T_max = self.AS.trivial_keyed_output(CP.iT_max)
        self._molar_mass = self.AS.trivial_keyed_output(CP.imolar_mass)

    def _is_below_T_critical(self, T):
        return T < self._T_crit

    def _make_p_subcritical(self, p):
        if p > self._p_crit:
            p = self._p_crit * 0.99
        return p

    def T_ph(self, p, h):
        self.AS.update(CP.HmassP_INPUTS, h, p)
        return self.AS.T()

    def T_ps(self, p, s):
        self.AS.update(CP.PSmass_INPUTS, p, s)
        return self.AS.T()

    def h_pT(self, p, T):
        self.AS.update(CP.PT_INPUTS, p, T)
        return self.AS.hmass()

    def h_QT(self, Q, T):
        self.AS.update(CP.QT_INPUTS, Q, T)
        return self.AS.hmass()

    def s_QT(self, Q, T):
        self.AS.update(CP.QT_INPUTS, Q, T)
        return self.AS.smass()

    def T_boiling(self, p):
        p = self._make_p_subcritical(p)
        self.AS.update(CP.PQ_INPUTS, p, 1)
        return self.AS.T()

    def p_boiling(self, T):
        if T > self._T_crit:
            T = self._T_crit * 0.99

        self.AS.update(CP.QT_INPUTS, 1, T)
        return self.AS.p()

    def Q_ph(self, p, h):
        p = self._make_p_subcritical(p)
        self.AS.update(CP.HmassP_INPUTS, h, p)
        return self.AS.Q()

    def d_ph(self, p, h):
        self.AS.update(CP.HmassP_INPUTS, h, p)
        return self.AS.rhomass()

    def d_pT(self, p, T):
        self.AS.update(CP.PT_INPUTS, p, T)
        return self.AS.rhomass()

    def d_QT(self, Q, T):
        self.AS.update(CP.QT_INPUTS, Q, T)
        return self.AS.rhomass()

    def viscosity_ph(self, p, h):
        self.AS.update(CP.HmassP_INPUTS, h, p)
        return self.AS.viscosity()

    def viscosity_pT(self, p, T):
        self.AS.update(CP.PT_INPUTS, p, T)
        return self.AS.viscosity()

    def s_ph(self, p, h):
        self.AS.update(CP.HmassP_INPUTS, h, p)
        return self.AS.smass()

    def s_pT(self, p, T):
        self.AS.update(CP.PT_INPUTS, p, T)
        return self.AS.smass()



h2o = CoolPropWrapper("H2O", "HEOS")
n2 = CoolPropWrapper("N2", "HEOS")
o2 = CoolPropWrapper("O2", "HEOS")


fluid_data = {
    "H2O": {
        "mass_fraction": 0.1,
        "property_object": h2o
    },
    "N2": {
        "mass_fraction": 0.7,
        "property_object": n2
    },
    "O2": {
        "mass_fraction": 0.2,
        "property_object": o2
    }
}


print(h_mix_pT(1e5, 327.2409, fluid_data, "ideal"))
print(h_mix_pT(1e5, 327.2409, fluid_data, "ideal-cond"))
print(T_mix_ph(1e5, 2.5e5, fluid_data, "ideal"))
print(viscosity_mix_pT_ideal(1e5, 350, fluid_data))

h = h_mix_pT_ideal_cond(1e5, 320, fluid_data)
T = T_mix_ph(1e5, h, fluid_data, "ideal-cond")
print(T)

s = s_mix_pT_ideal_cond(1e5, 300, fluid_data)
print(s)
T = T_mix_ps(1e5, s, fluid_data, "ideal-cond")

print(T)
# print(cond_check(1e5, T, fluid_data, "H2O"))