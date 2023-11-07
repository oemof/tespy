from tespy.tools.global_vars import ERR
from tespy.tools.logger import logger


def _is_larger_than_precision(value):
    return value > ERR


def _check_mixing_rule(mixing_rule, mixing_functions, propertyfunction):
    if mixing_rule not in mixing_functions:
        msg = (
            f"The mixing rule '{mixing_rule}' is not available for "
            f"the fluid property functions for {propertyfunction}. Available "
            f"rules are '" + "', '".join(mixing_functions.keys()) + "'."
        )
        raise KeyError(msg)


def get_number_of_fluids(fluid_data):
    return sum([1 for f in fluid_data.values() if _is_larger_than_precision(f["mass_fraction"])])


def get_pure_fluid(fluid_data):
    for f in fluid_data.values():
        if _is_larger_than_precision(f["mass_fraction"]):
            return f


def get_molar_fractions(fluid_data):
    molarflow = {
        key: value["mass_fraction"] / value["wrapper"]._molar_mass
        for key, value in fluid_data.items()
    }
    molarflow_sum = sum(molarflow.values())
    return {key: value / molarflow_sum for key, value in molarflow.items()}


def newton_with_kwargs(derivative, target_value, val0=300, valmin=70, valmax=3000, max_iter=10, tol_rel=ERR, tol_abs=ERR, tol_mode="rel", **function_kwargs):

    # start newton loop
    iteration = 0
    expr = True
    x = val0
    parameter = function_kwargs["parameter"]
    function = function_kwargs["function"]
    relax = 1

    while expr:
        # calculate function residual and new value
        function_kwargs[parameter] = x
        residual = target_value - function(**function_kwargs)
        x += residual / derivative(**function_kwargs) * relax

        # check for value ranges
        if x < valmin:
            x = valmin
        if x > valmax:
            x = valmax

        iteration += 1
        # relaxation to help convergence in case of jumping
        if iteration == 5:
            relax = 0.75
            max_iter = 12

        if iteration > max_iter:
            msg = (
                'The Newton algorithm was not able to find a feasible value '
                f'for function {function}. Current value with x={x} is '
                f'{function(**function_kwargs)}, target value is '
                f'{target_value}, residual is {residual} after {iteration} '
                'iterations.'
            )
            logger.debug(msg)

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


def inverse_temperature_mixture(p=None, target_value=None, fluid_data=None, T0=None, f=None, **function_kwargs):
    # calculate the fluid properties for fluid mixtures
    valmin, valmax = get_mixture_temperature_range(fluid_data)
    if T0 is None:
        T0 = (valmin + valmax) / 2.0
    T0 = max(valmin,min(valmax,T0)) 

    function_kwargs.update({
        "p": p, "fluid_data": fluid_data, "T": T0,
        "function": f, "parameter": "T" , "delta": 0.01
    })
    return newton_with_kwargs(
        central_difference,
        target_value,
        val0=T0,
        valmin=valmin,
        valmax=valmax,
        **function_kwargs
    )


def get_mixture_temperature_range(fluid_data):
    valmin = max(
        [v["wrapper"]._T_min for v in fluid_data.values() if _is_larger_than_precision(v["mass_fraction"])]
    ) + 0.1
    valmax = min(
        [v["wrapper"]._T_max for v in fluid_data.values() if _is_larger_than_precision(v["mass_fraction"])]
    ) - 0.1
    return valmin, valmax


def calc_molar_mass_mixture(fluid_data, molar_fractions):
    return sum([x * fluid_data[fluid]["wrapper"]._molar_mass for fluid, x in molar_fractions.items()])
