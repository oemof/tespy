# -*- coding: utf-8

"""
Extend pyvalem.Reaction class with custom methods.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/reactors/water_electrolyzer.py

SPDX-License-Identifier: MIT
"""

import pyromat as pm
from pyvalem.reaction import Reaction


def get_reactant_coefficients(self):
    return {k: -v for k, v in self.reactants_text_count_map.items()}


Reaction.get_reactant_coefficients = get_reactant_coefficients


def get_stoich_coefficients(self):
    d = self.get_reactant_coefficients().copy()
    d.update(self.products_text_count_map)
    return d


Reaction.get_stoich_coefficients = get_stoich_coefficients


def get_limiting_component(self, n0):
    stoich = self.reactants_text_count_map
    stoich_norm = {k: v / max(stoich.values()) for k, v in stoich.items()}
    X_max = {c: n0[c] / stoich_norm[c] for c in stoich.keys()}
    return min(X_max, key=X_max.get)


Reaction.get_limiting_component = get_limiting_component


def calculate_composition(self, n0, X=1):
    stoich = self.get_stoich_coefficients()
    lim_comp = self.get_limiting_component(n0)
    stoich_norm = {k: -v / stoich[lim_comp] for k, v in stoich.items()}
    return {
        c: n0[c] + X * stoich_norm[c] * n0[lim_comp]
        for c in stoich_norm.keys()
    }


Reaction.calculate_composition = calculate_composition


def get_formation_enthalpy(name, per_mole=True, T0=298.15):
    r"""
    Get enthalpy of formation (Δ_{f}H^{0}) of a molecule *in gaseous phase*.
    :param name: name of the molecule
    :param per_mole: True: return enthalpy in J/mol; False: return J/g
    :param T0: reference temperature in K (default: 298.15 K)
    :return: enthalpy of formation
    """
    ig = pm.get("ig." + name)
    h0 = ig.h(T0)[0]  # J/g
    if per_mole is True:
        h0 = h0 * ig.mw()  # J/mol
    return h0


def calculate_reaction_enthalpy(self, n0=None, X=1, T0=298.15):
    # assumes that all substances are in gas phase!
    stoich = self.get_stoich_coefficients()
    form_enthalpies = {
        k: get_formation_enthalpy(k, per_mole=True, T0=T0)
        for k in stoich.keys()
    }
    delta_h0 = sum(
        [stoich[k] * form_enthalpies[k] for k in stoich.keys()]
    )  # J/mol
    # calculate absolute enthalpy change if initial composition and
    # yield are given
    if n0 is not None and X is not None:
        lim_comp = self.get_limiting_component(n0)
        delta_h0 = X * n0[lim_comp] / -stoich[lim_comp] * delta_h0  # J
    return delta_h0


Reaction.calculate_reaction_enthalpy = calculate_reaction_enthalpy
