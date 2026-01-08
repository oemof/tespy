from pytest import approx
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.soo.nonconvex.de import DE
from pymoo.optimize import minimize
from tespy.tools import OptimizationProblem
import pandas as pd


class FakeOptimizationModel:

    def __init__(self):
        self.vars = {}
        self.obj1 = None
        self.obj2 = None

    def get_param(self, obj, label, param):
        return self.vars.get(label)

    def set_parameters(self, **kwargs):
        """
        Flattens the nested variable structure into self.vars
        """
        self.vars.clear()

        for group, content in kwargs.items():
            if group in ("Connections", "Components"):
                for label, data in content.items():
                    self.vars[label] = float(list(data.values())[0])
            elif group == "Customs":
                for label, value in content.items():
                    self.vars[label] = float(value)

    def solve_model(self, **kwargs):
        self.set_parameters(**kwargs)

        # Pull variables (default to 0 if missing)
        x = self.vars.get("x", 0.0)
        y = self.vars.get("y", 0.0)
        z = self.vars.get("z", 0.0)

        # Define objectives
        self.obj1 = (x - 2.0)**2 + (y - 1.0)**2 + z**2
        self.obj2 = x + y + z

    def get_objectives(self, objective_list):
        return [self.get_objective(obj) for obj in objective_list]

    def get_objective(self, obj):
        return getattr(self, obj)

    def penalize(self, fitness, constraint_evaluations):
        """
        constraint_evaluations = list of violation magnitudes
        """
        penalty = sum(
            (c ** 2) for c in constraint_evaluations if c > 0
        )

        return [f + penalty * 1000 for f in fitness]


def test_single_objective_single_variable_unconstrained():

    model = FakeOptimizationModel()

    variables = {
        "Connections": {"x": {"param": {"min": 0, "max": 10}}},
    }
    objective = ["obj1"]
    constraints = {}
    sense = [True]
    kpi = {
        "Connections": {"x": {"param"}},
        "Components": {"y": {"param"}},
        "Customs": {"z"}
    }

    problem = OptimizationProblem(
        model,
        variables=variables,
        constraints=constraints,
        objective=objective,
        minimize=sense,
        kpi=kpi,
        penalty_instead_of_constraints=False
    )

    algorithm = DE(pop_size=20)
    res = minimize(
        problem=problem,
        algorithm=algorithm,
        termination=("n_gen", 10),
        seed=42
    )

    df_log = pd.DataFrame(problem.log)
    best = df_log.sort_values(by=["obj1"], ascending=True).iloc[0]

    assert approx(best["obj1"]) == res.F[0]


def test_single_objective_single_variable_constrained():

    model = FakeOptimizationModel()

    variables = {
        "Connections": {"x": {"param": {"min": 0, "max": 10}}},
    }
    objective = ["obj1"]
    constraints = {
        "upper limits": {
            "Connections": {"x": {"param": 5}}
        }
    }
    sense = [True]
    kpi = {}

    problem = OptimizationProblem(
        model,
        variables=variables,
        constraints=constraints,
        objective=objective,
        minimize=sense,
        kpi=kpi,
        penalty_instead_of_constraints=False
    )

    algorithm = DE(pop_size=20)
    res = minimize(
        problem=problem,
        algorithm=algorithm,
        termination=("n_gen", 10),
        seed=42
    )

    df_log = pd.DataFrame(problem.log)
    mask = df_log["Connections-x-param<=5"] <= 0
    best = df_log.loc[mask].sort_values(by=["obj1"], ascending=True).iloc[0]

    assert approx(best["obj1"]) == res.F[0]


def test_single_objective_multi_variable_reference_unconstrained():

    model = FakeOptimizationModel()

    variables = {
        "Connections": {"x": {"param": {"min": 0, "max": 10}}},
        "Components": {"y": {"param": {"min": - 5, "max": 100}}},
        "Customs": {"z": {"min": 0, "max": 50}}
    }
    objective = ["obj1"]
    constraints = {}
    sense = [True]
    kpi = {"Components": {"y": {"param"}}}

    problem = OptimizationProblem(
        model,
        variables=variables,
        constraints=constraints,
        objective=objective,
        minimize=sense,
        kpi=kpi,
        penalty_instead_of_constraints=False
    )

    algorithm = DE(pop_size=20)
    res = minimize(
        problem=problem,
        algorithm=algorithm,
        termination=("n_gen", 10),
        seed=42
    )

    df_log = pd.DataFrame(problem.log)
    best = df_log.sort_values(by=["obj1"], ascending=True).iloc[0]

    assert approx(best["obj1"]) == res.F[0]


def test_single_objective_multi_variable_reference_constrained():

    model = FakeOptimizationModel()

    variables = {
        "Connections": {"x": {"param": {"min": 0, "max": 10}}},
        "Components": {"y": {"param": {"min": - 5, "max": 100}}},
        "Customs": {"z": {"min": 0, "max": 50}}
    }
    objective = ["obj1"]
    constraints = {
        "upper limits": {
            "Connections": {"x": {"param": "ref1"}}
        },
        "ref1": ["Components", "y", "param"]
    }
    sense = [True]
    kpi = {"Components": {"y": {"param"}}}

    problem = OptimizationProblem(
        model,
        variables=variables,
        constraints=constraints,
        objective=objective,
        minimize=sense,
        kpi=kpi,
        penalty_instead_of_constraints=False
    )

    algorithm = DE(pop_size=20)
    res = minimize(
        problem=problem,
        algorithm=algorithm,
        termination=("n_gen", 10),
        seed=42
    )

    df_log = pd.DataFrame(problem.log)
    mask = df_log["Connections-x-param<=Components-y-param"] <= 0
    best = df_log.loc[mask].sort_values(by=["obj1"], ascending=True).iloc[0]

    assert approx(best["obj1"]) == res.F[0]


def test_multi_objective_multi_variable_reference_unconstrained():

    model = FakeOptimizationModel()

    variables = {
        "Connections": {"x": {"param": {"min": 0, "max": 10}}},
        "Components": {"y": {"param": {"min": - 5, "max": 100}}},
        "Customs": {"z": {"min": 0, "max": 50}}
    }
    objective = ["obj1", "obj2"]
    constraints = {}
    sense = [True, False]
    kpi = {"Components": {"y": {"param"}}}

    problem = OptimizationProblem(
        model,
        variables=variables,
        constraints=constraints,
        objective=objective,
        minimize=sense,
        kpi=kpi,
        penalty_instead_of_constraints=False
    )

    algorithm = NSGA2(pop_size=20)
    res = minimize(
        problem=problem,
        algorithm=algorithm,
        termination=("n_gen", 10),
        seed=42
    )

    df_log = pd.DataFrame(problem.log)

    assert approx(res.F[:, 0].min()) == df_log["obj1"].min()
    assert approx(-res.F[:, 1].min()) == df_log["obj2"].max()


def test_multi_objective_multi_variable_reference_constrained():

    model = FakeOptimizationModel()

    variables = {
        "Connections": {"x": {"param": {"min": 0, "max": 10}}},
        "Components": {"y": {"param": {"min": - 5, "max": 100}}},
        "Customs": {"z": {"min": 0, "max": 50}}
    }
    objective = ["obj1", "obj2"]
    constraints = {
        "upper limits": {
            "Connections": {"x": {"param": "ref1"}}
        },
        "ref1": ["Components", "y", "param"]
    }
    sense = [True, False]
    kpi = {"Components": {"y": {"param"}}}

    problem = OptimizationProblem(
        model,
        variables=variables,
        constraints=constraints,
        objective=objective,
        minimize=sense,
        kpi=kpi,
        penalty_instead_of_constraints=False
    )

    algorithm = NSGA2(pop_size=20)
    res = minimize(
        problem=problem,
        algorithm=algorithm,
        termination=("n_gen", 10),
        seed=42
    )

    df_log = pd.DataFrame(problem.log)
    mask = df_log["Connections-x-param<=Components-y-param"] <= 0
    df_log = df_log.loc[mask]

    assert approx(res.F[:, 0].min()) == df_log["obj1"].min()
    assert approx(-res.F[:, 1].min()) == df_log["obj2"].max()


def test_multi_objective_multi_variable_reference_constrained_with_penalty():

    model = FakeOptimizationModel()

    variables = {
        "Connections": {"x": {"param": {"min": 0, "max": 10}}},
        "Components": {"y": {"param": {"min": - 5, "max": 100}}},
        "Customs": {"z": {"min": 0, "max": 50}}
    }
    objective = ["obj1", "obj2"]
    constraints = {
        "upper limits": {
            "Connections": {"x": {"param": "ref1"}}
        },
        "ref1": ["Components", "y", "param"]
    }
    sense = [True, False]
    kpi = {"Components": {"y": {"param"}}}

    problem = OptimizationProblem(
        model,
        variables=variables,
        constraints=constraints,
        objective=objective,
        minimize=sense,
        kpi=kpi,
        penalty_instead_of_constraints=True
    )

    algorithm = NSGA2(pop_size=20)
    res = minimize(
        problem=problem,
        algorithm=algorithm,
        termination=("n_gen", 10),
        seed=42
    )

    df_log = pd.DataFrame(problem.log)
    mask = df_log["Connections-x-param<=Components-y-param"] <= 0
    df_log = df_log.loc[mask]

    assert approx(res.F[:, 0].min()) == df_log["obj1"].min()
    assert approx(-res.F[:, 1].min()) == df_log["obj2"].max()
