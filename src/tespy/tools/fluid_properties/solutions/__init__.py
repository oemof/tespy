# -*- coding: utf-8

"""Package for the working pairs of absorption cycles.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/tools/fluid_properties/solutions/__init__.py

SPDX-License-Identifier: MIT
"""
from .base import SolutionComponentWrapper  # noqa: F401
from .base import SolutionPair  # noqa: F401
from .base import get_solution_pair  # noqa: F401
from .base import register_solution_pair  # noqa: F401
from .base import solution_pair_names  # noqa: F401
from .base import solution_pair_of_component  # noqa: F401
from .libr_water import LiBrWater
from .nh3_water import AmmoniaWater

register_solution_pair(LiBrWater())
register_solution_pair(AmmoniaWater())
