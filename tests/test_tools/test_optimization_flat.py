import pandas as pd
import pytest
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.soo.nonconvex.de import DE
from pymoo.optimize import minimize
from pytest import approx

from tespy.tools import OptimizationProblem


class FakeModel:
    """Minimal model compatible with the OptimizationProblem flat API."""

    def __init__(self):
        self.x = 0.0
        self.y = 0.0
        self.obj1 = None
        self.obj2 = None

    def solve_model(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, float(v))
        self.obj1 = (self.x - 2.0) ** 2 + (self.y - 1.0) ** 2
        self.obj2 = self.x + self.y

    def get_parameter(self, name):
        return getattr(self, name)

    def get_objectives(self, objective_list):
        return [self.get_parameter(obj) for obj in objective_list]

    def penalize(self, fitness, constraint_evaluations):
        penalty = sum(c ** 2 for c in constraint_evaluations if c > 0)
        return [f + penalty * 1000 for f in fitness]


# --- single objective ---

def test_single_objective_unconstrained():
    model = FakeModel()
    problem = OptimizationProblem(
        model,
        variables={"x": {"min": 0, "max": 10}, "y": {"min": 0, "max": 10}},
        objective=["obj1"],
        minimize=[True],
    )
    algorithm = DE(pop_size=20)
    res = minimize(problem, algorithm, termination=("n_gen", 10), seed=42)

    df_log = pd.DataFrame(problem.log)
    best = df_log.sort_values("obj1").iloc[0]
    assert approx(best["obj1"]) == res.F[0]


def test_single_objective_numeric_constraint():
    model = FakeModel()
    problem = OptimizationProblem(
        model,
        variables={"x": {"min": 0, "max": 10}, "y": {"min": 0, "max": 10}},
        constraints={"x": {"max": 3.0}},
        objective=["obj1"],
        minimize=[True],
    )
    algorithm = DE(pop_size=20)
    res = minimize(problem, algorithm, termination=("n_gen", 10), seed=42)

    df_log = pd.DataFrame(problem.log)
    mask = df_log["x<=3.0"] <= 0
    best = df_log.loc[mask].sort_values("obj1").iloc[0]
    assert approx(best["obj1"]) == res.F[0]


def test_single_objective_cross_parameter_constraint():
    model = FakeModel()
    problem = OptimizationProblem(
        model,
        variables={"x": {"min": 0, "max": 10}, "y": {"min": 0, "max": 10}},
        constraints={"x": {"min": "y"}},
        objective=["obj1"],
        minimize=[True],
    )
    algorithm = DE(pop_size=20)
    res = minimize(problem, algorithm, termination=("n_gen", 10), seed=42)

    df_log = pd.DataFrame(problem.log)
    mask = df_log["x>=y"] <= 0
    best = df_log.loc[mask].sort_values("obj1").iloc[0]
    assert approx(best["obj1"]) == res.F[0]


# --- multi objective ---

def test_multi_objective_unconstrained():
    model = FakeModel()
    problem = OptimizationProblem(
        model,
        variables={"x": {"min": 0, "max": 10}, "y": {"min": 0, "max": 10}},
        objective=["obj1", "obj2"],
        minimize=[True, False],
    )
    algorithm = NSGA2(pop_size=20)
    res = minimize(problem, algorithm, termination=("n_gen", 10), seed=42)

    df_log = pd.DataFrame(problem.log)
    assert approx(res.F[:, 0].min()) == df_log["obj1"].min()
    assert approx(-res.F[:, 1].min()) == df_log["obj2"].max()


def test_multi_objective_constrained():
    model = FakeModel()
    problem = OptimizationProblem(
        model,
        variables={"x": {"min": 0, "max": 10}, "y": {"min": 0, "max": 10}},
        constraints={"x": {"max": "y"}},
        objective=["obj1", "obj2"],
        minimize=[True, False],
    )
    algorithm = NSGA2(pop_size=20)
    res = minimize(problem, algorithm, termination=("n_gen", 10), seed=42)

    df_log = pd.DataFrame(problem.log)
    mask = df_log["x<=y"] <= 0
    df_log = df_log.loc[mask]
    assert approx(res.F[:, 0].min()) == df_log["obj1"].min()
    assert approx(-res.F[:, 1].min()) == df_log["obj2"].max()


def test_multi_objective_constrained_with_penalty():
    model = FakeModel()
    problem = OptimizationProblem(
        model,
        variables={"x": {"min": 0, "max": 10}, "y": {"min": 0, "max": 10}},
        constraints={"x": {"max": "y"}},
        objective=["obj1", "obj2"],
        minimize=[True, False],
        penalty_instead_of_constraints=True,
    )
    assert problem.n_ieq_constr == 0
    algorithm = NSGA2(pop_size=20)
    minimize(problem, algorithm, termination=("n_gen", 10), seed=42)
    assert len(problem.log) > 0


# --- kpi and logging ---

def test_kpi_logged():
    model = FakeModel()
    problem = OptimizationProblem(
        model,
        variables={"x": {"min": 0, "max": 10}},
        objective=["obj1"],
        minimize=[True],
        kpi=["y"],
    )
    algorithm = DE(pop_size=10)
    minimize(problem, algorithm, termination=("n_gen", 5), seed=42)

    df_log = pd.DataFrame(problem.log)
    assert "y" in df_log.columns


# --- default argument isolation ---

def test_none_defaults_do_not_mutate():
    """Passing no optional args must not share state across instances."""
    model = FakeModel()
    p1 = OptimizationProblem(model, variables={"x": {"min": 0, "max": 1}}, objective=["obj1"])
    p2 = OptimizationProblem(model, variables={"x": {"min": 0, "max": 1}}, objective=["obj1"])
    p1.log.append({"x": 1})
    assert p2.log == []


# --- argument validation ---

def test_objective_not_list_raises():
    model = FakeModel()
    with pytest.raises(TypeError, match="list"):
        OptimizationProblem(model, variables={"x": {"min": 0, "max": 1}}, objective="obj1")


def test_minimize_length_mismatch_raises():
    model = FakeModel()
    with pytest.raises(ValueError):
        OptimizationProblem(
            model,
            variables={"x": {"min": 0, "max": 1}},
            objective=["obj1", "obj2"],
            minimize=[True],
        )
