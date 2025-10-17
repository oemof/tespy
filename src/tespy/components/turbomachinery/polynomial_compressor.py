import warnings

from tespy.components.displacementmachinery.polynomial_compressor import calc_EN12900 as cEN12900
from tespy.components.displacementmachinery.polynomial_compressor import fit_EN12900 as fEN12900
from tespy.components.displacementmachinery.polynomial_compressor import \
    generate_eta_polys_from_data as gen_polys_from_data
from tespy.components.displacementmachinery.polynomial_compressor import \
    generate_eta_polys_from_power_and_cooling_polys as gen_polys_from_polys
from tespy.tools.logger import logger


def calc_EN12900(c: list, t_evap: float, t_cond: float) -> float:
    msg = (
        "This method has moved. Import it from the "
        "tespy.components.displacementmachinery.polynomial_compressor module "
        "instead!"
    )
    logger.warning(msg)
    warnings.warn(msg, FutureWarning)
    return cEN12900(c, t_evap, t_cond)


def fit_EN12900(t_evap: np.array, t_cond: np.array, data: np.array) -> np.array:
    msg = (
        "This method has moved. Import it from the "
        "tespy.components.displacementmachinery.polynomial_compressor module "
        "instead!"
    )
    logger.warning(msg)
    warnings.warn(msg, FutureWarning)
    return fEN12900(t_evap, t_cond, data)


def generate_eta_polys_from_power_and_cooling_polys(power_poly: list, cooling_poly: list, t_evap: np.array, t_cond: np.array, fluid: str, reference_state: dict) -> tuple:
    msg = (
        "This method has moved. Import it from the "
        "tespy.components.displacementmachinery.polynomial_compressor module "
        "instead!"
    )
    logger.warning(msg)
    warnings.warn(msg, FutureWarning)
    return gen_polys_from_polys(power_poly, cooling_poly, t_evap, t_cond, fluid, reference_state)


def generate_eta_polys_from_data(df_power, df_cooling, fluid: str, reference_state: dict) -> tuple:
    msg = (
        "This method has moved. Import it from the "
        "tespy.components.displacementmachinery.polynomial_compressor module "
        "instead!"
    )
    logger.warning(msg)
    warnings.warn(msg, FutureWarning)
    return gen_polys_from_data(df_power, df_cooling, fluid, fluid, reference_state)
