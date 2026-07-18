# -*- coding: utf-8

"""Module for the block solve strategies.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/solver/strategies.py

SPDX-License-Identifier: MIT
"""
import numpy as np
from numpy.linalg import norm

from tespy.components.component import Component
from tespy.connections.connection import ConnectionBase
from tespy.tools import helpers as hlp
from tespy.tools import logger
from tespy.tools.global_vars import ERR


class NewtonStrategy:
    r"""
    Solve a single block with the newton algorithm restricted to the block.

    The residuals and derivatives of all objects owning the block's
    equations are evaluated, the linear system is restricted to the block's
    equation rows and variable columns and only the block's variables are
    updated. Variables of preceding blocks act as constants.

    If the block does not converge or an equation evaluation raises an
    error, the block's variables are restored to their values at entry, so
    a failed attempt does not degrade the starting point of any subsequent
    solution attempt.
    """

    def solve(self, problem, block):
        """Solve the block and return its status (0, 2 or 3)."""
        snapshot = {}
        for col in block.variables:
            data = problem.variables_dict[col]
            if data["variable"] == "fluid":
                snapshot[col] = data["obj"].val[data["fluid"]]
            else:
                snapshot[col] = data["obj"]._val_SI

        try:
            status = self._iterate(problem, block)
        except hlp.TESPyNetworkError:
            # informative errors, e.g. on NaN residuals pointing at the
            # equation that produced them, must reach the user
            raise
        except Exception as e:
            msg = f"Solving the block raised an error: {e}"
            logger.warning(msg)
            status = 2

        if status != 0:
            for col, value in snapshot.items():
                data = problem.variables_dict[col]
                if data["variable"] == "fluid":
                    data["obj"].val[data["fluid"]] = value
                else:
                    data["obj"]._val_SI = value

        return status

    def _iterate(self, problem, block):
        equations = block.equations
        variables = block.variables
        variable_set = set(variables)

        objects = {problem._equation_obj_lookup[eq] for eq in equations}
        for col in variables:
            for sm_col in problem.variables_dict[col]["_represents"]:
                objects.add(problem._variable_lookup[sm_col]["object"])

        connections = [o for o in objects if isinstance(o, ConnectionBase)]
        components = [o for o in objects if isinstance(o, Component)]
        equation_objects = [
            problem._equation_obj_lookup[eq] for eq in equations
        ]
        equation_objects = list(dict.fromkeys(equation_objects))

        block.residual_history = []
        # seed the increment so the increment filter does not suppress the
        # first derivative evaluation of the block's variables
        problem.increment[variables] = 1
        # stabilization measures of the equations are staged on the objects'
        # iteration counters, every block starts a fresh count
        for obj in equation_objects:
            obj.it = 0
        problem._prev_residual = None
        for block_iter in range(problem.max_iter):
            problem.iter = block_iter
            problem.increment_filter = (
                np.absolute(problem.increment) < ERR ** 2
            )
            problem._solve_equations(equation_objects)
            problem._invert_jacobian(equations, variables)
            if problem.lin_dep:
                return 3

            if problem.oscillation_damping:
                problem._dampen_oscillating_increments()
            self._modify_increment(problem, block)
            problem._update_variables(variable_set)
            modified = problem._adapt_to_variable_bounds(
                connections, components
            )
            problem._prev_residual = problem.residual.copy()

            block.residual_history.append(
                float(norm(problem.residual[equations]))
            )
            # accept once the residual and the increment are small in an
            # iteration without interference of the value heuristics and at
            # least one newton update has been applied before the residual
            # evaluation. A deeply converged residual is accepted right
            # away (exactly solved blocks), otherwise two consecutive small
            # residuals are required, forcing another newton iteration
            # while convergence is still in progress
            deep = block.residual_history[-1] < ERR
            confirmed = (
                len(block.residual_history) >= 2
                and (
                    np.array(block.residual_history[-2:]) < ERR ** 0.5
                ).all()
            )
            if (
                    block_iter >= 1
                    and not modified
                    and (deep or confirmed)
                    and (
                        np.abs(problem.increment[variables]) < ERR ** 0.25
                    ).all()
                ):
                return 0

        return 2

    def _modify_increment(self, problem, block):
        pass


class ScalarBracketingStrategy(NewtonStrategy):
    r"""
    Solve a scalar block with newton steps and bracketing as fallback.

    When the block's residual changes sign between two iterations, the
    newton step overshot the root - but the two iterates then enclose it.
    The increment is replaced by a refinement of the root within that known
    bracket (Brent's method), which restores monotone convergence for
    non-smooth residuals. Only if the refinement fails, e.g. because
    heuristics modified the state between the iterations, a bracketing
    search on the equation is performed instead.
    """

    def __init__(self):
        self._prev_residual = None
        self._prev_value = None

    def _modify_increment(self, problem, block):
        row = block.equations[0]
        col = block.variables[0]
        residual = problem.residual[row]
        data = problem.variables_dict[col]
        if data["variable"] == "fluid":
            value = data["obj"].val[data["fluid"]]
        else:
            value = data["obj"]._val_SI

        if (
                self._prev_residual is not None
                and self._prev_residual * residual < 0
            ):
            step = problem._bracketed_step(row, col, self._prev_value)
            if step is None:
                step = problem._search_reducing_step(row, col)
            if step is not None:
                problem.increment[col] = step

        self._prev_residual = residual
        self._prev_value = value
