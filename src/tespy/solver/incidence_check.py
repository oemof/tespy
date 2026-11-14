# -*- coding: utf-8

"""Module for checking the incidence against the observed dependence.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/solver/incidence_check.py

SPDX-License-Identifier: MIT
"""
import numpy as np

from tespy.tools import helpers as hlp

FLUID_STEP = 1e-4
SCALAR_STEP = 1e-4
CHANGE_TOLERANCE = 1e-12


def _residuals(problem):
    residuals = {}
    for obj in problem._all_objects():
        for data in obj.equations.values():
            values = np.atleast_1d(
                np.asarray(data.func(**data.func_params), dtype=float)
            )
            for i, value in enumerate(values[:data.num_eq]):
                residuals[data._first_eq_index + i] = value
    return residuals


def _step(data, value):
    if data["variable"] == "fluid":
        return -FLUID_STEP if value > 0.5 else FLUID_STEP
    return max(abs(value), 1) * SCALAR_STEP


def check_dependents(problem):
    r"""Compare the incidence of a problem with the observed dependence.

    Every solver variable is perturbed once and the equations whose
    residual changes are recorded. The observation only holds for the
    current state: a dependence vanishing there, e.g. on the composition
    of saturated humid air, is not seen.

    Parameters
    ----------
    problem : tespy.solver.problem.Problem
        Problem with its variable space loaded.

    Returns
    -------
    dict
        Per entry a dictionary of equation labels (object label and equation
        name) to the short labels of the variables concerned:

        - :code:`"missing"`: dependence observed but not part of the
          incidence, neither declared nor inferred. The block order can be
          wrong.
        - :code:`"unused_inferred"`: dependence on the composition inferred
          but not observed.
        - :code:`"unused_declared"`: dependence declared but not observed.

        :code:`"not_probed"` lists the variables whose perturbation raised
        an error in an equation evaluation.
    """
    incidence = problem._incidence_matrix
    inferred = problem._inferred_incidence
    base = _residuals(problem)

    observed = {eq_num: set() for eq_num in base}
    not_probed = []
    for col, data in problem.variables_dict.items():
        value = hlp.get_variable_value(data)
        hlp.set_variable_value(data, value + _step(data, value))
        try:
            perturbed = _residuals(problem)
        except Exception:
            not_probed.append(problem.format_var_label(col))
            continue
        finally:
            hlp.set_variable_value(data, value)
        for eq_num, residual in perturbed.items():
            scale = max(1, abs(residual), abs(base[eq_num]))
            if abs(residual - base[eq_num]) > CHANGE_TOLERANCE * scale:
                observed[eq_num].add(col)

    report = {"missing": {}, "unused_inferred": {}, "unused_declared": {}}
    for eq_num, columns in observed.items():
        label = problem._equation_lookup[eq_num]
        assumed = set(incidence.get(eq_num, []))
        assumed_by_inference = set(inferred.get(eq_num, []))
        groups = {
            "missing": columns - assumed,
            "unused_inferred": assumed_by_inference - columns,
            "unused_declared": assumed - assumed_by_inference - columns,
        }
        for key, cols in groups.items():
            if cols:
                report[key][label] = [
                    problem.format_var_label(col) for col in sorted(cols)
                ]
    report["not_probed"] = not_probed
    return report
