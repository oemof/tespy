# -*- coding: utf-8
import importlib.resources
import os
import sys
import warnings

if sys.version_info < (3, 11):
    warnings.warn(
        f"Python {sys.version_info.major}.{sys.version_info.minor} is no "
        f"longer supported as of the next major release of tespy. "
        f"Please upgrade to Python 3.11 or newer.",
        FutureWarning
    )

__datapath__ = os.path.join(importlib.resources.files("tespy"), "data")
__version__ = '0.9.16.post2 - Kelvin\'s Kingdom'

# tespy data and connections import
from . import connections  # noqa: F401
from . import data  # noqa: F401
# tespy components imports
from .components import basics  # noqa: F401
from .components import combustion  # noqa: F401
from .components import component  # noqa: F401
from .components import heat_exchangers  # noqa: F401
from .components import nodes  # noqa: F401
from .components import piping  # noqa: F401
from .components import reactors  # noqa: F401
from .components import subsystem  # noqa: F401
from .components import turbomachinery  # noqa: F401
# tespy networks imports
from .networks import network  # noqa: F401
# tespy tools imports
from .tools import characteristics  # noqa: F401
from .tools import data_containers  # noqa: F401
from .tools import fluid_properties  # noqa: F401
from .tools import global_vars  # noqa: F401
from .tools import helpers  # noqa: F401
from .tools import logger  # noqa: F401
