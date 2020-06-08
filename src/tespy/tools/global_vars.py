# -*- coding: utf-8

"""Module for global variables used by other modules of the tespy package.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/global_vars.py

SPDX-License-Identifier: MIT
"""

err = 1e-6
molar_masses = {}
gas_constants = {}
gas_constants['uni'] = 8.3144598
# Define colors for highlighting values in result table
coloring = {
    'end': '\033[0m',
    'set': '\033[94m',
    'err': '\033[31m',
    'var': '\033[32m'
}
