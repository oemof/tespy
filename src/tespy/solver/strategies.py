# -*- coding: utf-8

"""Module for the block solve strategies.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/solver/strategies.py

SPDX-License-Identifier: MIT
"""
import numpy as np
from numpy.linalg import norm
from scipy.optimize import brentq

from tespy.components.component import Component
from tespy.connections.connection import ConnectionBase
from tespy.solver.decomposition import Block
from tespy.tools import helpers as hlp
from tespy.tools import logger
from tespy.tools.global_vars import ERR
from tespy.tools.global_vars import INCREMENT_TOLERANCE
from tespy.tools.global_vars import RESIDUAL_TOLERANCE


def scalar_residual_closure(problem, row, col):
    """Build an evaluator of equation row over variable col.

    The evaluator computes the equation's residual at a given value of
    the variable and restores the current value afterwards. Works for
    both scalar variables (m, h, p, E) and vector variables (fluid mass
    fractions), where the largest other variable fluid is adjusted in
    the opposite direction to maintain the sum of 1.

    Returns the evaluator and the current variable value, or None in
    case the equation cannot be evaluated standalone.
    """
    obj = problem._equation_obj_lookup.get(row)
    if obj is None:
        return None
    _, (param_name, sub_idx) = problem._equation_lookup[row]
    if param_name not in obj.equations:
        return None
    data = obj.equations[param_name]

    var_data = problem.variables_dict[col]
    container = var_data["obj"]

    if var_data["variable"] == "fluid":
        fluid_key = var_data["fluid"]
        x0 = container.val[fluid_key]
        # Maintain sum=1 by adjusting the largest other variable fluid by
        # the same delta in the opposite direction.
        other_var_fluids = [f for f in container.is_var if f != fluid_key]
        if other_var_fluids:
            companion = max(other_var_fluids, key=lambda f: container.val.get(f, 0))
            companion_x0 = container.val[companion]
        else:
            companion = None
            companion_x0 = None

        def set_x(v):
            container.val[fluid_key] = v
            if companion is not None:
                container.val[companion] = companion_x0 - (v - x0)
    else:
        x0 = container._val_SI

        def set_x(v):
            container._val_SI = v

    def eval_r(x):
        set_x(x)
        try:
            result = data.func(**data.func_params)
        except Exception:
            return None
        finally:
            set_x(x0)
        if hasattr(result, '__iter__'):
            result = list(result)
            return result[sub_idx] if sub_idx < len(result) else result[0]
        return result

    return eval_r, x0


def bracketed_step(problem, row, col, bracket_value):
    """Refine the root of a scalar equation within a known bracket.

    The bracket value is the variable value of a preceding newton
    iteration whose residual had the opposite sign of the current one,
    so the root is enclosed between the two values and can be refined
    directly without searching for a bracket first.

    Returns the step to add to the variable, or None in case the
    refinement fails, e.g. because the residual function changed
    between the iterations through other variables or heuristics.
    """
    closure = scalar_residual_closure(problem, row, col)
    if closure is None:
        return None
    eval_r, x0 = closure

    # each iteration costs a single equation evaluation, so the root is
    # refined to close to machine precision right away
    try:
        tol = max(abs(x0) * 1e-12, 1e-12)
        x_root = brentq(
            eval_r, min(x0, bracket_value), max(x0, bracket_value),
            xtol=tol, maxiter=60
        )
        return x_root - x0
    except Exception:
        return None


def search_reducing_step(problem, row, col):
    """Find the increment for variable col that reduces equation row's
    residual.

    Searches both +/- directions with geometrically growing step sizes
    (x2 per iteration, up to 20 iterations each). Prefers the side that
    produces a sign change in the residual, which guarantees a root in
    the bracket [x0, x0±d] by the IVT, and refines its location with
    Brent's method. If both sides bracket a root, the tighter one
    (smaller |r| at the probe point) is used. Falls back to a secant
    step if brentq raises, and to the lower-magnitude heuristic when
    neither side yields a sign change.

    Returns the step to add to the variable, or None if neither direction
    improves the residual.
    """
    closure = scalar_residual_closure(problem, row, col)
    if closure is None:
        return None
    eval_r, x0 = closure

    r0 = problem.residual[row]
    abs_r0 = abs(r0)

    # Guard against x0 == 0 producing a zero initial step
    d = max(abs(x0) * 0.1, 1e-3)
    found_plus = None
    found_minus = None
    # outermost probes whose residual still equals the starting residual:
    # they tighten the bracket, e.g. cutting off most of a residual plateau
    # the search walked across
    same_plus = 0.0
    same_minus = 0.0

    for _ in range(20):
        if found_plus is None:
            r = eval_r(x0 + d)
            if r is not None and r != r0:
                found_plus = (d, r)
            elif r is not None:
                same_plus = d

        if found_minus is None:
            r = eval_r(x0 - d)
            if r is not None and r != r0:
                found_minus = (d, r)
            elif r is not None:
                same_minus = d

        if found_plus is not None and found_minus is not None:
            break
        d *= 2

    plus_sign_change = found_plus is not None and r0 * found_plus[1] < 0
    minus_sign_change = found_minus is not None and r0 * found_minus[1] < 0

    if plus_sign_change or minus_sign_change:
        # Both sides bracket a root: prefer the tighter probe (smaller |r|)
        if plus_sign_change and minus_sign_change:
            plus_d, plus_r = found_plus
            minus_d, minus_r = found_minus
            if abs(plus_r) <= abs(minus_r):
                sign, step_d, r_val = +1, plus_d, plus_r
            else:
                sign, step_d, r_val = -1, minus_d, minus_r
        elif plus_sign_change:
            sign, step_d, r_val = +1, found_plus[0], found_plus[1]
        else:
            sign, step_d, r_val = -1, found_minus[0], found_minus[1]

        a = x0 + sign * (same_plus if sign > 0 else same_minus)
        b = x0 + sign * step_d
        try:
            tol = max(abs(x0) * 1e-6, 1e-10)
            # residual plateaus (e.g. temperatures inside the two phase
            # region) degrade brentq to bisection, which needs the iterations
            x_root = brentq(
                eval_r, min(a, b), max(a, b), xtol=tol, maxiter=60
            )
            return x_root - x0
        except Exception:
            pass

        # Secant fallback: linear interpolation between x0 and the probe
        return sign * step_d * (-r0) / (r_val - r0)

    # No sign change found - fall back to lower-magnitude direction
    if found_plus is None and found_minus is None:
        return None
    if found_plus is None:
        step_d, r_val = found_minus
        return -step_d if abs(r_val) < abs_r0 else None
    if found_minus is None:
        step_d, r_val = found_plus
        return +step_d if abs(r_val) < abs_r0 else None

    plus_d, plus_r = found_plus
    minus_d, minus_r = found_minus
    if abs(plus_r) <= abs(minus_r):
        return +plus_d if abs(plus_r) < abs_r0 else None
    else:
        return -minus_d if abs(minus_r) < abs_r0 else None


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
        snapshot = {
            col: hlp.get_variable_value(problem.variables_dict[col])
            for col in block.variables
        }

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
                hlp.set_variable_value(problem.variables_dict[col], value)

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
                    np.array(block.residual_history[-2:])
                    < RESIDUAL_TOLERANCE
                ).all()
            )
            if (
                    block_iter >= 1
                    and not modified
                    and (deep or confirmed)
                    and (
                        np.abs(problem.increment[variables])
                        < INCREMENT_TOLERANCE
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
        value = hlp.get_variable_value(problem.variables_dict[col])

        if (
                self._prev_residual is not None
                and self._prev_residual * residual < 0
            ):
            step = bracketed_step(problem, row, col, self._prev_value)
            if step is None:
                step = search_reducing_step(problem, row, col)
            if step is not None:
                problem.increment[col] = step

        self._prev_residual = residual
        self._prev_value = value


class BlockDriver:
    r"""
    Orchestrate the block wise solution of a problem.

    The blocks of the decomposition are solved in precedence order and the
    variables of every solved block are frozen for the subsequent ones.
    Failures escalate in stages: a failed block restores its variables and
    the remaining blocks are solved as one coupled system; if that fails as
    well, the full system is solved simultaneously from its initial state,
    reproducing the plain simultaneous solution exactly. After a successful
    block phase the full residual vector is verified and a simultaneous
    polish runs in case couplings outside of the declared incidence left a
    remaining residual.
    """

    def __init__(self, problem):
        self.problem = problem
        self._frozen = []

    def solve(self, iterinfo, print_results=True):
        problem = self.problem
        if self._variable_composition():
            msg = (
                "The fluid composition is part of the variable space. The "
                "dependency declarations of the fluid property equations do "
                "not cover the composition, so a block ordering is not "
                "reliable in this case: solving the full system "
                "simultaneously instead."
            )
            logger.info(msg)
            problem._solve_simultaneous(iterinfo, print_results)
            return

        decomposition = problem.decompose()
        if decomposition.defective_blocks:
            problem.singularity_msg = (
                "The problem is structurally singular, block-wise solving "
                f"is not possible.{problem._structural_report()}"
            )
            problem.status = 3
            return

        self._log_decomposition(decomposition)

        # snapshot to restore the starting state in the final escalation
        # stage, so it behaves exactly like the default simultaneous solver
        snapshot = {
            key: hlp.get_variable_value(data)
            for key, data in problem.variables_dict.items()
        }

        num_blocks = len(decomposition.blocks)
        for position, block in enumerate(decomposition.blocks):
            self._solve_block(block, num_blocks, iterinfo, print_results)

            if block.status != 0:
                # the preceding blocks stay solved: their equations do not
                # depend on any variable of the failed or later blocks. The
                # failed block's variables were restored by the strategy, so
                # the coupled solution of everything remaining starts from
                # the solved blocks plus the original starting values.
                remainder = self._remainder_of(decomposition, position)
                eq_str = ", ".join(
                    f"{lbl}.{problem.format_eq_name(name)}"
                    for lbl, name in block.equation_labels
                )
                var_str = ", ".join(block.variable_labels)
                msg = (
                    f"Block {block.id} did not converge, solving the "
                    f"remaining {num_blocks - position} blocks "
                    "simultaneously.\n"
                    f"  Equations: {eq_str}\n"
                    f"  Variables: {var_str}"
                )
                logger.warning(msg)

                self._solve_block(
                    remainder, num_blocks, iterinfo, print_results
                )
                if remainder.status != 0:
                    msg = (
                        "The remaining system did not converge either, "
                        "restarting with the simultaneous solution of the "
                        "full system from its initial state."
                    )
                    logger.warning(msg)
                    self._restart(snapshot, iterinfo, print_results)
                    return
                break

            self._freeze(block)

        if not self._verified():
            self._unfreeze()
            problem.increment = np.ones([problem.variable_counter])
            problem._solve_simultaneous(iterinfo, print_results)
            return

        self._unfreeze()
        problem.status = 0

    def _variable_composition(self):
        return any(
            data["variable"] == "fluid"
            for data in self.problem.variables_dict.values()
        )

    def _log_decomposition(self, decomposition):
        kinds = {}
        for block in decomposition.blocks:
            kinds[block.kind] = kinds.get(block.kind, 0) + 1
        kind_str = ", ".join(f"{num} {kind}" for kind, num in kinds.items())
        msg = (
            f"Block decomposition: {len(decomposition.blocks)} blocks "
            f"({kind_str})."
        )
        logger.debug(msg)
        for block in decomposition.blocks:
            eq_str = ", ".join(
                f"{lbl}.{self.problem.format_eq_name(name)}"
                for lbl, name in block.equation_labels
            )
            var_str = ", ".join(block.variable_labels)
            msg = (
                f"Block {block.id} ({block.kind}): equations [{eq_str}], "
                f"variables [{var_str}]"
            )
            logger.debug(msg)

    def _solve_block(self, block, num_blocks, iterinfo, print_results):
        if block.kind == "scalar":
            strategy = ScalarBracketingStrategy()
        else:
            strategy = NewtonStrategy()
        block.status = strategy.solve(self.problem, block)

        residual = (
            block.residual_history[-1] if block.residual_history else 0.0
        )
        self.problem.residual_history = np.append(
            self.problem.residual_history, residual
        )

        if iterinfo:
            if block.kind == "remainder":
                progress = 100
                detail = f"{len(block.equations)} equations"
            else:
                progress = 100 * (block.id + 1) // num_blocks
                detail = ", ".join(
                    f"{lbl}.{self.problem.format_eq_name(name)}"
                    for lbl, name in block.equation_labels
                )
            msg = (
                f" block {block.id:3d} | {block.kind:15s} | "
                f"iterations: {len(block.residual_history):3d} | "
                f"residual: {residual:.2e} | {detail}"
            )
            logger.progress(progress, msg)
            if print_results and not logger.console_logging_enabled():
                print(msg)

    def _remainder_of(self, decomposition, position):
        remaining = decomposition.blocks[position:]
        remainder = Block(
            id=decomposition.blocks[position].id,
            kind="remainder",
            equations=sorted(eq for b in remaining for eq in b.equations),
            variables=sorted(col for b in remaining for col in b.variables),
        )
        equations = self.problem.get_equations()
        remainder.equation_labels = [
            equations[eq] for eq in remainder.equations
        ]
        remainder.variable_labels = [
            self.problem.format_var_label(col) for col in remainder.variables
        ]
        return remainder

    def _freeze(self, block):
        # the is_var guards of the convergence heuristics and of the
        # derivative computation skip frozen variables automatically in all
        # subsequent blocks
        for col in block.variables:
            container = self.problem.variables_dict[col]["obj"]
            container.is_var = False
            self._frozen.append(container)

    def _unfreeze(self):
        for container in self._frozen:
            container.is_var = True
        self._frozen = []

    def _restart(self, snapshot, iterinfo, print_results):
        problem = self.problem
        self._unfreeze()
        for key, value in snapshot.items():
            hlp.set_variable_value(problem.variables_dict[key], value)
        # iteration counters stage stabilization measures of the equations,
        # the restart must begin from counter zero like a fresh solve
        network = problem.network
        to_reset = (
            network.comps["object"].tolist()
            + network.conns["object"].tolist()
            + list(network.user_defined_eq.values())
        )
        for obj in to_reset:
            obj.it = 0
        problem.increment = np.ones([problem.variable_counter])
        problem._solve_simultaneous(iterinfo, print_results)

    def _verified(self):
        """Evaluate the full residual vector as convergence certificate.

        Couplings acting outside of the declared incidence (e.g. equations
        internally manipulating their residual during iteration) can leave
        a global residual the per-block convergence cannot see.
        """
        problem = self.problem
        if problem.variable_counter == 0:
            return True

        problem._solve_equations(residual_only=True)
        residual_norm = norm(problem.residual)
        if residual_norm > RESIDUAL_TOLERANCE:
            msg = (
                "Block-wise solving finished with a remaining global "
                f"residual of {residual_norm:.2e}, continuing with the "
                "simultaneous solution of the full system."
            )
            logger.info(msg)
            return False

        msg = (
            "Block solution verified, global residual "
            f"{residual_norm:.2e}."
        )
        logger.debug(msg)
        return True
