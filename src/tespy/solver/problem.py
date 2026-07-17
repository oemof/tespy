# -*- coding: utf-8

"""Module for the solver problem representation.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/solver/problem.py

SPDX-License-Identifier: MIT
"""
import math
from time import time

import numpy as np
from numpy.linalg import norm
from scipy.optimize import brentq

from tespy.tools import helpers as hlp
from tespy.tools import logger
from tespy.tools.data_containers import VectorVariable as dc_vecvar
from tespy.tools.global_vars import ERR

# Only require cupy if Cuda shall be used
try:
    import cupy as cu
except ModuleNotFoundError:
    cu = None


class Problem:
    r"""
    Solver representation of a prepared network.

    The problem holds the variable space, the equation lookups and the state
    of the newton algorithm (residual vector, jacobian and increment). It is
    created by the network in the preprocessing stage of every :code:`solve`
    call and remains accessible through :code:`Network.problem` afterwards
    for inspection.
    """

    def __init__(self, network):
        self.network = network

        self.variable_counter = 0
        self.variables_dict = {}

        self.num_comp_eq = 0
        self.num_conn_eq = 0
        self.num_ude_eq = 0

        self._equation_lookup = {}
        self._equation_obj_lookup = {}
        self._incidence_matrix = {}
        self._incidence_matrix_dense = None

        self.status = 99
        self.iter = 0
        self.residual = None
        self.jacobian = None
        self.increment = None
        self.lin_dep = False
        self.singularity_msg = ""

    def solve_loop(self, max_iter, min_iter, use_cuda, robust_relax,
                   oscillation_damping, iterinfo, print_results=True):
        r"""Loop of the newton algorithm."""
        self.max_iter = max_iter
        self.min_iter = min_iter
        self.use_cuda = use_cuda
        self.robust_relax = robust_relax
        self.oscillation_damping = oscillation_damping

        # parameter definitions
        self.residual_history = np.array([])
        self.residual = np.zeros([self.variable_counter])
        self.increment = np.ones([self.variable_counter])
        self.jacobian = np.zeros((self.variable_counter, self.variable_counter))
        self._prev_residual = None

        self.start_time = time()

        if iterinfo:
            self._print_iterinfo_head(print_results)

        for self.iter in range(self.max_iter):
            self.increment_filter = np.absolute(self.increment) < ERR ** 2
            self._solve_iteration()
            self.residual_history = np.append(
                self.residual_history, norm(self.residual)
            )
            if iterinfo:
                self._print_iterinfo_body(print_results)

            if self.lin_dep:
                self.status = 3
                break

            elif self.iter > 40:
                if (
                    all(
                        self.residual_history[(self.iter - 3):] >= self.residual_history[-3] * 0.95
                    ) and self.residual_history[-1] >= self.residual_history[-2] * 0.95
                ):
                    self.status = 2
                    break

            if (
                    self.iter >= self.min_iter - 1
                    and (self.residual_history[-2:] < ERR ** 0.5).all()
                    # the increment should also be small, but it does not need
                    # to be that small
                    and (abs(self.increment) < ERR ** 0.25).all()
                ):
                self.status = 0
                break

        self.end_time = time()

        if iterinfo:
            self._print_iterinfo_tail(print_results)

        if self.iter == self.max_iter - 1:
            msg = (
                f"Reached maximum iteration count ({self.max_iter}), "
                "calculation stopped. Residual value is "
                "{:.2e}. ".format(norm(self.residual)) +
                "\nPossible reasons include:\n"
                " - fluid properties moving outside the valid range of the "
                "property database (consider adjusting p_range or h_range),\n"
                " - an impossible constraint that can never be satisfied \n"
                " - bad starting values causing the Newton solver to diverge.\n"
                "Use nw.print_residuals() to identify which equations have "
                "the largest residuals."
            )
            logger.warning(msg)
            self.status = 2

    def _print_iterinfo_head(self, print_results=True):
        """Print head of convergence progress."""
        # Start with defining the format here
        self.iterinfo_fmt = ' {iter:5s} | {residual:10s} | {progress:10s} '
        self.iterinfo_fmt += '| {massflow:10s} | {pressure:10s} | {enthalpy:10s} '
        self.iterinfo_fmt += '| {fluid:10s} | {energy:10s} | {component:10s} '
        # Use the format to create the first logging entry
        msg = self.iterinfo_fmt.format(
            iter='iter',
            residual='residual',
            progress='progress',
            massflow='massflow',
            pressure='pressure',
            enthalpy='enthalpy',
            fluid='fluid',
            energy='energy',
            component='component'
        )
        logger.progress(0, msg)
        msg2 = '-' * 7 + '+------------' * 8

        logger.progress(0, msg2)
        if print_results:
            print('\n' + msg + '\n' + msg2)

    def _print_iterinfo_body(self, print_results=True):
        """Print convergence progress."""
        m = [k for k, v in self.variables_dict.items() if v["variable"] == "m"]
        p = [k for k, v in self.variables_dict.items() if v["variable"] == "p"]
        h = [k for k, v in self.variables_dict.items() if v["variable"] == "h"]
        fl = [k for k, v in self.variables_dict.items() if v["variable"] == "fluid"]
        e = [k for k, v in self.variables_dict.items() if v["variable"] == "E"]
        cp = [k for k in self.variables_dict if k not in m + p + h + fl + e]

        iter_str = str(self.iter + 1)
        residual_norm = norm(self.residual)
        residual = 'NaN'
        progress = 'NaN'
        massflow = 'NaN'
        pressure = 'NaN'
        enthalpy = 'NaN'
        fluid = 'NaN'
        energy = 'NaN'
        component = 'NaN'

        progress_val = -1

        if not np.isnan(residual_norm):
            residual = '{:.2e}'.format(residual_norm)

            if not self.lin_dep:
                massflow = '{:.2e}'.format(norm(self.increment[m]))
                pressure = '{:.2e}'.format(norm(self.increment[p]))
                enthalpy = '{:.2e}'.format(norm(self.increment[h]))
                fluid = '{:.2e}'.format(norm(self.increment[fl]))
                energy  = '{:.2e}'.format(norm(self.increment[e]))
                component  = '{:.2e}'.format(norm(self.increment[cp]))

            # This should not be hardcoded here.
            if residual_norm > np.finfo(float).eps * 100:
                progress_min = math.log(ERR)
                progress_max = math.log(ERR ** 0.5) * -1
                progress_val = math.log(max(residual_norm, ERR)) * -1
                # Scale to 0-1
                progress_scaled = (
                    (progress_val - progress_min)
                    / (progress_max - progress_min)
                )
                progress_val = max(0, min(1, progress_scaled))
                # Scale to 100%
                progress_val = int(progress_val * 100)
            else:
                progress_val = 100

            progress = '{:d} %'.format(progress_val)

        msg = self.iterinfo_fmt.format(
            iter=iter_str,
            residual=residual,
            progress=progress,
            massflow=massflow,
            pressure=pressure,
            enthalpy=enthalpy,
            fluid=fluid,
            energy=energy,
            component=component
        )
        logger.progress(progress_val, msg)
        if print_results:
            print(msg)

    def _print_iterinfo_tail(self, print_results=True):
        """Print tail of convergence progress."""
        num_iter = self.iter + 1
        clc_time = self.end_time - self.start_time
        num_ips = num_iter / clc_time if clc_time > 1e-10 else np.inf
        msg = '-' * 7 + '+------------' * 7
        logger.progress(100, msg)
        msg = (
            "Total iterations: {0:d}, Calculation time: {1:.2f} s, "
            "Iterations per second: {2:.2f}"
        ).format(num_iter, clc_time, num_ips)
        logger.debug(msg)
        if print_results:
            print(msg)
        return

    def _search_reducing_step(self, row, col):
        """Find the increment for variable col that reduces equation row's
        residual.

        Searches both +/- directions with geometrically growing step sizes
        (x2 per iteration, up to 20 iterations each). Works for both scalar
        variables (m, h, p, E) and vector variables (fluid mass fractions).
        Prefers the side that produces a sign change in the residual, which
        guarantees a root in the bracket [x0, x0±d] by the IVT, and refines
        its location with Brent's method. If both sides bracket a root, the
        tighter one (smaller |r| at the probe point) is used. Falls back to a
        secant step if brentq raises, and to the lower-magnitude heuristic
        when neither side yields a sign change.

        Returns the step to add to the variable, or None if neither direction
        improves the residual.
        """
        obj = self._equation_obj_lookup.get(row)
        if obj is None:
            return None
        _, (param_name, sub_idx) = self._equation_lookup[row]
        if param_name not in obj.equations:
            return None
        data = obj.equations[param_name]

        var_data = self.variables_dict[col]
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

        r0 = self.residual[row]
        abs_r0 = abs(r0)

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

        # Guard against x0 == 0 producing a zero initial step
        d = max(abs(x0) * 0.1, 1e-3)
        found_plus = None
        found_minus = None

        for _ in range(20):
            if found_plus is None:
                r = eval_r(x0 + d)
                if r is not None and r != r0:
                    found_plus = (d, r)

            if found_minus is None:
                r = eval_r(x0 - d)
                if r is not None and r != r0:
                    found_minus = (d, r)

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

            a = x0
            b = x0 + sign * step_d
            try:
                tol = max(abs(x0) * 1e-6, 1e-10)
                x_root = brentq(
                    eval_r, min(a, b), max(a, b), xtol=tol, maxiter=5
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

    def _fill_jacobian_surrogates(self):
        """Restore invertibility for all-zero rows and find better steps.

        For each row that is entirely zero but expected to have non-zero
        entries (per the incidence matrix), inserts 1 in the expected positions
        so the Jacobian can be inverted for all other variables.
        Subsequently searches value of associated variable(s) to find the
        increment for the affected variable(s) that reduces that equation's
        residual.

        Returns a dict {col: step} of increment overrides to apply after the
        inversion.
        """
        overrides = {}
        for row in self._check_all_zero_rows(self.jacobian):
            for col in self._incidence_matrix.get(row, []):
                if self.jacobian[row, col] == 0.0:
                    self.jacobian[row, col] = 1.0
                    step = self._search_reducing_step(row, col)
                    if step is not None:
                        overrides[col] = step
        return overrides

    def _invert_jacobian(self):
        """Compute Newton step, storing result in self.increment. Sets self.lin_dep."""
        self.lin_dep = False
        self.increment = self.residual * 0
        if len(self.variables_dict) == 0:
            return

        self._check_residual_and_jacobian_for_nan()

        overrides = self._fill_jacobian_surrogates()

        try:
            if self.use_cuda:
                self.increment = cu.asnumpy(cu.dot(
                    cu.linalg.inv(cu.asarray(self.jacobian)),
                    -cu.asarray(self.residual)
                ))
            else:
                row_scales = np.abs(self.jacobian).max(axis=1)
                row_scales[row_scales == 0] = 1.0
                J_eq = self.jacobian / row_scales[:, None]

                col_scales = np.maximum(np.abs(J_eq).max(axis=0), 1e-10)
                J_sc = J_eq / col_scales[None, :]
                r_eq = -self.residual / row_scales

                self.increment = np.linalg.solve(J_sc, r_eq) / col_scales
        except np.linalg.LinAlgError:
            self.lin_dep = True
            return

        for col, step in overrides.items():
            self.increment[col] = step

    def _check_residual_and_jacobian_for_nan(self):
        """Raise an informative error if a NaN or inf entry is found.

        A single NaN in the residual or Jacobian typically contaminates the
        whole solution vector once the linear system is solved, so the bad
        value has to be located here, before the solve, to be able to point
        at the equation that produced it. Otherwise it propagates silently
        into the next iteration's variable values and only surfaces several
        iterations later as an unrelated and confusing error.
        """
        nan_rows = set(np.where(~np.isfinite(self.residual))[0])
        nan_rows |= set(np.where(~np.isfinite(self.jacobian).any(axis=1))[0])
        if len(nan_rows) == 0:
            return

        eq_str = ", ".join(
            f"{lbl}.{self._format_eq_name(eq_name)}"
            for lbl, eq_name in self._get_equations_by_number(nan_rows).values()
        )
        msg = (
            "The residual or one of the partial derivatives of the "
            f"following equation(s) is NaN or inf: {eq_str}. Please check "
            "the inputs to and the implementation of these equations, e.g. "
            "a UserDefinedEquation."
        )
        logger.error(msg)
        raise hlp.TESPyNetworkError(msg)

    def _diagnose_singularity(self):
        """Build singularity_msg after a failed matrix solve."""
        if self.iter == 0 and np.linalg.matrix_rank(self._incidence_matrix_dense) < self._incidence_matrix_dense.shape[0]:
            self.singularity_msg = (
                "Detected singularity in Jacobian matrix. This singularity "
                "is most likely caused by the parametrization of your "
                "problem and NOT a numerical issue. Double check your "
                "setup.\n"
            )
            self._find_linear_dependencies(self.jacobian)
            return

        expected_entries = self._incidence_matrix_dense.astype(bool)
        actual_entries = self.jacobian.astype(bool)
        rows, cols = np.where(expected_entries != actual_entries)

        missing_entries = []
        for row, col in zip(rows, cols):
            lbl, eq_name = self._equation_lookup[row]
            eq_str = f"{lbl}.{self._format_eq_name(eq_name)}"
            var_str = self._format_var_label(col)
            missing_entries += [f"{eq_str}: {var_str}"]

        entries_str = ", ".join(missing_entries)
        self.singularity_msg = (
            "Found singularity in Jacobian matrix, calculation aborted! "
            "The setup of your problem seems to be solvable. It failed "
            "due to partial derivatives in the Jacobian being zero where "
            "a non-zero was expected, or vice versa. This usually lies in "
            "starting value selection or bad convergence.\n"
            "  The following equation/variable pairs may have an "
            f"unexpected zero/non-zero partial derivative:  {entries_str}\n"
        )
        self._find_linear_dependencies(self.jacobian)

    def _find_linear_dependencies(self, matrix):
        all_zero_cols = self._check_all_zero_columns(matrix)
        all_zero_rows = self._check_all_zero_rows(matrix)
        if len(all_zero_cols) + len(all_zero_rows) == 0:
            eq_indices = self._cauchy_schwarz_inequality(matrix)
            eq_str = ", ".join(
                f"{lbl}.{self._format_eq_name(eq_name)}"
                for lbl, eq_name in self._get_equations_by_number(eq_indices).values()
            )
            self.singularity_msg += (
                "The following equations form a linear dependency:\n"
                f"  {eq_str}\n"
            )
        else:
            if len(all_zero_cols) > 0:
                var_str = ", ".join(self._format_var_label(i) for i in all_zero_cols)
                self.singularity_msg += (
                    "The following variables are not associated with any equation:\n"
                    f"  {var_str}\n"
                )
            if len(all_zero_rows) > 0:
                eq_str = ", ".join(
                    f"{lbl}.{self._format_eq_name(eq_name)}"
                    for lbl, eq_name in self._get_equations_by_number(all_zero_rows).values()
                )
                self.singularity_msg += (
                    "The following equations do not depend on any variable:\n"
                    f"  {eq_str}\n"
                )

    def _check_all_zero_columns(self, matrix):
        return np.where((matrix == 0).all(axis=0))[0]

    def _check_all_zero_rows(self, matrix):
        return np.where((matrix == 0).all(axis=1))[0]

    def _cauchy_schwarz_inequality(self, matrix):
        n = matrix.shape[0]
        dependent_equations = []
        for i in range(n):
            for j in range(n):
                if i != j:
                    inner_product = np.inner(
                        matrix[i,:],
                        matrix[j,:]
                    )
                    norm_i = np.linalg.norm(matrix[i,:])
                    norm_j = np.linalg.norm(matrix[j,:])

                    if np.abs(inner_product - norm_j * norm_i) < 1e-5:
                        dependent_equations += [i]
        return list(set(dependent_equations))

    def _dampen_oscillating_increments(self):
        """Halve increments for variables whose residuals changed sign since the last iteration.

        A sign change means the Newton step overshot the zero of that equation.
        Halving the dependent columns' increments keeps the bracket from growing
        and makes subsequent iterations behave like bisection for those variables.
        """
        if self._prev_residual is None:
            return
        sign_flips = np.where(self._prev_residual * self.residual < 0)[0]
        for row in sign_flips:
            cols = np.nonzero(self.jacobian[row, :])[0]
            self.increment[cols] *= 0.5

    def _update_variables(self):
        # cast dtype to float from numpy float64
        # this is necessary to keep the doctests running and note make them
        # look ugly all over the place
        # I have yet to come up with a better idea, or vectorize all operations
        # which requires major changes in tespy
        increment = [float(val) for val in self.increment]
        # the J_cols here point to actual variables, no need to call to
        # get_J_col yet
        relax = 1
        if self.robust_relax:
            relax = 0.05 + 0.95 * min(1, self.iter / (0.25 * self.max_iter))

        for _, data in self.variables_dict.items():
            if data["variable"] in ["m", "h", "E"]:
                container = data["obj"]
                container._val_SI += increment[container.J_col] * relax
            elif data["variable"] in ["p"]:
                container = data["obj"]
                p_relax = max(
                    1, -2 * increment[container.J_col] / container.val_SI
                )
                container._val_SI += increment[container.J_col] / p_relax
            elif data["variable"] == "fluid":
                container = data["obj"]
                inc = increment[container.J_col[data["fluid"]]]
                value = container.val[data["fluid"]]
                if value > ERR:
                    f_relax = max(1, -2 * inc / value)
                else:
                    f_relax = 1

                container.val[data["fluid"]] += inc / f_relax

                if container.val[data["fluid"]] < ERR:
                    container.val[data["fluid"]] = 0
                elif container.val[data["fluid"]] > 1 - ERR:
                    container.val[data["fluid"]] = 1
            else:
                # add increment
                data["obj"]._val_SI += increment[data["obj"].J_col] * relax

                # keep value within specified value range
                if data["obj"].val_SI < data["obj"].min_val:
                    data["obj"].val_SI = data["obj"].min_val
                elif data["obj"].val_SI > data["obj"].max_val:
                    data["obj"].val_SI = data["obj"].max_val

    def _adapt_to_variable_bounds(self):

        # this could be in a different place, its kind of in between
        # network and connection
        if self.iter < 10:
            for data in self.variables_dict.values():
                if type(data["obj"]) == dc_vecvar:
                    total_mass_fractions = sum(data["obj"].val.values())
                    for fluid in data["obj"].is_var:
                        data["obj"]._val[fluid] /= total_mass_fractions

        if norm(self.increment) > 1e-1:
            for c in self.network.conns['object']:
                # check the fluid properties for physical ranges
                c._adjust_to_property_limits(self.network)

            for cp in self.network.comps['object']:
                cp._adjust_to_property_limits()

        # second check based on component heuristics
        # - for first three iterations
        # - only if the increment is sufficiently large
        # - only in design case
        if (
                self.iter < 3
                and norm(self.increment) > 1e-1
                and self.network.mode == "design"
            ):
            for cp in self.network.comps['object']:
                cp.convergence_check()

            for c in self.network.conns['object']:
                c._adjust_to_property_limits(self.network)

    def _solve_iteration(self):
        r"""
        Control iteration step of the newton algorithm.

        - Calculate the residual value for each equation
        - Calculate the jacobian matrix
        - Calculate new values for variables
        - Restrict fluid properties to value ranges
        - Check component parameters for consistency
        """
        self._solve_equations()
        self._invert_jacobian()

        if self.lin_dep:
            self._diagnose_singularity()
            return

        if self.oscillation_damping:
            self._dampen_oscillating_increments()

        self._update_variables()
        self._adapt_to_variable_bounds()
        self._prev_residual = self.residual.copy()

    def _solve_equations(self):
        r"""
        Calculate the residual and derivatives of all equations.
        """
        to_solve = (
            self.network.comps["object"].tolist()
            + self.network.conns["object"].tolist()
            + list(self.network.user_defined_eq.values())
        )
        for obj in to_solve:
            hlp.solve(obj, self.increment_filter)
            if len(obj.jacobian) > 0:
                rows = list(obj.residual.keys())
                data = list(obj.residual.values())
                self.residual[rows] = data

                rows = [k[0] for k in obj.jacobian]
                columns = [k[1] for k in obj.jacobian]
                data = list(obj.jacobian.values())
                self.jacobian[rows, columns] = data

            obj.it += 1

    def _format_var_label(self, v_idx):
        v_data = self.variables_dict[v_idx]
        v_type = v_data["variable"]
        if v_type == "fluid" and v_data["fluid"] is not None:
            return f"{v_data['fluid']}{v_idx}"
        return f"{v_type}{v_idx}"

    @staticmethod
    def _format_eq_name(eq_name):
        if isinstance(eq_name, tuple):
            name, sub_idx = eq_name
            return f"{name}{{{sub_idx}}}" if sub_idx > 0 else name
        return eq_name

    def _get_equations_by_number(self, number_list) -> dict:
        """Get the actual equations after presolving the problem by equation
        number

        Returns
        -------
        dict
            Lookup with equation number as index and tuple of label and
            parameter defining the equation. In case one parameter defines
            multiple equations, the same equation is repeated.
        """
        return {
            k: v for k, v in self._equation_lookup.items() if k in number_list
        }

    def get_sorted_residual_index(self) -> list[int]:
        """Get the sorted array of residual indices.

        Returns
        -------
        list[int]
            List of variable numbers, the index values.
        """
        return list(np.argsort(np.abs(self.residual))[::-1])
