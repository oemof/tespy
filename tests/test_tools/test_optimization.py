# -*- coding: utf-8

"""Tests of the OptimizationProblem with the flat parameter interface.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tests/test_tools/test_optimization.py

SPDX-License-Identifier: MIT
"""
import pandas as pd
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.soo.nonconvex.de import DE
from pymoo.optimize import minimize
from pytest import approx

from tespy.tools import OptimizationProblem


class FakeOptimizationModel:
    """Model with the flat interface: parameters are addressed by name."""

    def __init__(self):
        self.vars = {}
        self.obj1 = None
        self.obj2 = None

    def solve_model(self, **kwargs):
        self.vars = {name: float(value) for name, value in kwargs.items()}
        x = self.vars.get("x", 0.0)
        y = self.vars.get("y", 0.0)
        z = self.vars.get("z", 0.0)
        self.obj1 = (x - 2.0) ** 2 + (y - 1.0) ** 2 + z ** 2
        self.obj2 = x + y + z

    def get_parameter(self, name):
        if name in ("obj1", "obj2"):
            return getattr(self, name)
        return self.vars.get(name)

    def get_objectives(self, objective_list):
        return [self.get_parameter(obj) for obj in objective_list]

    def penalize(self, fitness, constraint_evaluations):
        penalty = sum(c ** 2 for c in constraint_evaluations if c > 0)
        return [f + penalty * 1000 for f in fitness]


def _optimize(problem, algorithm):
    return minimize(
        problem=problem, algorithm=algorithm, termination=("n_gen", 10),
        seed=42
    )


def test_single_objective_single_variable_unconstrained():
    problem = OptimizationProblem(
        FakeOptimizationModel(),
        variables={"x": {"min": 0, "max": 10}},
        objective=["obj1"],
        minimize=[True],
        kpi=["x", "y", "z"],
    )
    res = _optimize(problem, DE(pop_size=20))

    df_log = pd.DataFrame(problem.log)
    best = df_log.sort_values(by=["obj1"], ascending=True).iloc[0]
    assert approx(best["obj1"]) == res.F[0]


def test_single_objective_single_variable_constrained():
    problem = OptimizationProblem(
        FakeOptimizationModel(),
        variables={"x": {"min": 0, "max": 10}},
        constraints={"x": {"max": 5}},
        objective=["obj1"],
        minimize=[True],
    )
    res = _optimize(problem, DE(pop_size=20))

    df_log = pd.DataFrame(problem.log)
    mask = df_log["x<=5"] <= 0
    best = df_log.loc[mask].sort_values(by=["obj1"], ascending=True).iloc[0]
    assert approx(best["obj1"]) == res.F[0]


def test_single_objective_multi_variable_reference_unconstrained():
    problem = OptimizationProblem(
        FakeOptimizationModel(),
        variables={
            "x": {"min": 0, "max": 10},
            "y": {"min": -5, "max": 100},
            "z": {"min": 0, "max": 50},
        },
        objective=["obj1"],
        minimize=[True],
        kpi=["y"],
    )
    res = _optimize(problem, DE(pop_size=20))

    df_log = pd.DataFrame(problem.log)
    best = df_log.sort_values(by=["obj1"], ascending=True).iloc[0]
    assert approx(best["obj1"]) == res.F[0]


def test_single_objective_multi_variable_reference_constrained():
    problem = OptimizationProblem(
        FakeOptimizationModel(),
        variables={
            "x": {"min": 0, "max": 10},
            "y": {"min": -5, "max": 100},
            "z": {"min": 0, "max": 50},
        },
        constraints={"x": {"max": "y"}},
        objective=["obj1"],
        minimize=[True],
        kpi=["y"],
    )
    res = _optimize(problem, DE(pop_size=20))

    df_log = pd.DataFrame(problem.log)
    mask = df_log["x<=y"] <= 0
    best = df_log.loc[mask].sort_values(by=["obj1"], ascending=True).iloc[0]
    assert approx(best["obj1"]) == res.F[0]


def test_multi_objective_multi_variable_reference_unconstrained():
    problem = OptimizationProblem(
        FakeOptimizationModel(),
        variables={
            "x": {"min": 0, "max": 10},
            "y": {"min": -5, "max": 100},
            "z": {"min": 0, "max": 50},
        },
        objective=["obj1", "obj2"],
        minimize=[True, False],
        kpi=["y"],
    )
    res = _optimize(problem, NSGA2(pop_size=20))

    df_log = pd.DataFrame(problem.log)
    assert approx(res.F[:, 0].min()) == df_log["obj1"].min()
    assert approx(-res.F[:, 1].min()) == df_log["obj2"].max()


def test_multi_objective_multi_variable_reference_constrained():
    problem = OptimizationProblem(
        FakeOptimizationModel(),
        variables={
            "x": {"min": 0, "max": 10},
            "y": {"min": -5, "max": 100},
            "z": {"min": 0, "max": 50},
        },
        constraints={"x": {"max": "y"}, "z": {"min": 1}},
        objective=["obj1", "obj2"],
        minimize=[True, False],
        kpi=["y"],
    )
    res = _optimize(problem, NSGA2(pop_size=20))

    df_log = pd.DataFrame(problem.log)
    mask = (df_log["x<=y"] <= 0) & (df_log["z>=1"] <= 0)
    assert approx(res.F[:, 0].min()) == df_log.loc[mask, "obj1"].min()
    assert approx(-res.F[:, 1].min()) == df_log.loc[mask, "obj2"].max()


def test_multi_objective_multi_variable_reference_constrained_with_penalty():
    problem = OptimizationProblem(
        FakeOptimizationModel(),
        variables={
            "x": {"min": 0, "max": 10},
            "y": {"min": -5, "max": 100},
            "z": {"min": 0, "max": 50},
        },
        constraints={"x": {"max": "y"}, "z": {"min": 1}},
        objective=["obj1", "obj2"],
        minimize=[True, False],
        kpi=["y"],
        penalty_instead_of_constraints=True,
    )
    res = _optimize(problem, NSGA2(pop_size=20))

    df_log = pd.DataFrame(problem.log)
    assert problem.n_ieq_constr == 0
    assert len(res.F) > 0
    assert set(df_log.columns) >= {"x", "y", "z", "obj1", "obj2", "x<=y", "z>=1"}
