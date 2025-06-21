# -*- coding: utf-8
import importlib.resources
import os
import sys
import warnings

if sys.version_info[1] < 10:
    msg = (
        "Supprt for python versions below 3.10 will be dropped with the "
        "next major release"
    )
    warnings.warn(FutureWarning(msg))

__datapath__ = os.path.join(importlib.resources.files("tespy"), "data")
__version__ = '0.8.3 - Newton\'s Nature'

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
from .networks import network_reader  # noqa: F401
# tespy tools imports
from .tools import characteristics  # noqa: F401
from .tools import data_containers  # noqa: F401
from .tools import fluid_properties  # noqa: F401
from .tools import global_vars  # noqa: F401
from .tools import helpers  # noqa: F401
from .tools import logger  # noqa: F401
