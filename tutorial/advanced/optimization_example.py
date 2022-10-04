# %%[sec_1]
import numpy as np
import pygmo as pg

from tespy.components import CycleCloser
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Condenser
from tespy.components import Desuperheater
from tespy.components import HeatExchangerSimple
from tespy.components import Merge
from tespy.components import Splitter
from tespy.components import Valve
from tespy.components import Pump
from tespy.components import Turbine
from tespy.connections import Bus
from tespy.connections import Connection
from tespy.networks import Network

from tespy.tools.optimization import OptimizationProblem


class SamplePlant:
    """Class template for TESPy model usage in optimization module."""
    def __init__(self):

        self.nw = Network(fluids=["water"])
        self.nw.set_attr(
            p_unit="bar", T_unit="C", h_unit="kJ / kg", iterinfo=False
        )
        # components
        # main cycle
        sg = HeatExchangerSimple("steam generator")
        cc = CycleCloser("cycle closer")
        hpt = Turbine("high pressure turbine")
        sp1 = Splitter("splitter 1", num_out=2)
        mpt = Turbine("mid pressure turbine")
        sp2 = Splitter("splitter 2", num_out=2)
        lpt = Turbine("low pressure turbine")
        con = Condenser("condenser")
        pu1 = Pump("feed water pump")
        fwh1 = Condenser("feed water preheater 1")
        fwh2 = Condenser("feed water preheater 2")
        dsh = Desuperheater("desuperheater")
        me2 = Merge("merge2", num_in=2)
        pu2 = Pump("feed water pump 2")
        pu3 = Pump("feed water pump 3")
        me = Merge("merge", num_in=2)

        # cooling water
        cwi = Source("cooling water source")
        cwo = Sink("cooling water sink")

        # connections
        # main cycle
        c0 = Connection(sg, "out1", cc, "in1", label="0")
        c1 = Connection(cc, "out1", hpt, "in1", label="1")
        c2 = Connection(hpt, "out1", sp1, "in1", label="2")
        c3 = Connection(sp1, "out1", mpt, "in1", label="3", state="g")
        c4 = Connection(mpt, "out1", sp2, "in1", label="4")
        c5 = Connection(sp2, "out1", lpt, "in1", label="5")
        c6 = Connection(lpt, "out1", con, "in1", label="6")
        c7 = Connection(con, "out1", pu1, "in1", label="7", state="l")
        c8 = Connection(pu1, "out1", fwh1, "in2", label="8", state="l")
        c9 = Connection(fwh1, "out2", me, "in1", label="9", state="l")
        c10 = Connection(me, "out1", fwh2, "in2", label="10", state="l")
        c11 = Connection(fwh2, "out2", dsh, "in2", label="11", state="l")
        c12 = Connection(dsh, "out2", me2, "in1", label="12", state="l")
        c13 = Connection(me2, "out1", sg, "in1", label="13", state="l")

        self.nw.add_conns(
            c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13
        )

        # preheating
        c21 = Connection(sp1, "out2", dsh, "in1", label="21")
        c22 = Connection(dsh, "out1", fwh2, "in1", label="22")
        c23 = Connection(fwh2, "out1", pu2, "in1", label="23")
        c24 = Connection(pu2, "out1", me2, "in2", label="24")

        c31 = Connection(sp2, "out2", fwh1, "in1", label="31")
        c32 = Connection(fwh1, "out1", pu3, "in1", label="32")
        c33 = Connection(pu3, "out1", me, "in2", label="33")

        self.nw.add_conns(
            c21, c22, c23, c24,
            c31, c32, c33
        )

        # cooling water
        c41 = Connection(cwi, "out1", con, "in2", label="41")
        c42 = Connection(con, "out2", cwo, "in1", label="42")

        self.nw.add_conns(c41, c42)

        # busses
        # power bus
        self.power = Bus("power")
        self.power.add_comps(
            {"comp": hpt, "char": -1}, {"comp": mpt, "char": -1},
            {"comp": lpt, "char": -1}, {"comp": pu1, "char": -1},
            {"comp": pu2, "char": -1}, {"comp": pu3, "char": -1}
        )

        # heating bus
        self.heat = Bus("heat")
        self.heat.add_comps({"comp": sg, "char": 1})

        self.nw.add_busses(self.power, self.heat)

        # parametrization
        # components
        hpt.set_attr(eta_s=0.9)
        mpt.set_attr(eta_s=0.9)
        lpt.set_attr(eta_s=0.9)

        pu1.set_attr(eta_s=0.8)
        pu2.set_attr(eta_s=0.8)
        pu3.set_attr(eta_s=0.8)

        sg.set_attr(pr=0.92)

        con.set_attr(pr1=1, pr2=0.99, ttd_u=5)
        fwh1.set_attr(pr1=1, pr2=0.99, ttd_u=5)
        fwh2.set_attr(pr1=1, pr2=0.99, ttd_u=5)
        dsh.set_attr(pr1=0.99, pr2=0.99)

        c1.set_attr(m=200, T=650, p=100, fluid={"water": 1})
        c2.set_attr(p=20)
        c4.set_attr(p=3)

        c41.set_attr(T=20, p=3, fluid={"water": 1})
        c42.set_attr(T=28)

        self.nw.solve("design")
        self.stable = "_stable"
        self.nw.save(self.stable)
        self.solved = True

    # %%[sec_2]

    def get_param(self, obj, label, parameter):
        """Get the value of a parameter in the network"s unit system.

        Parameters
        ----------
        obj : str
            Object to get parameter for (Components/Connections).

        label : str
            Label of the object in the TESPy model.

        parameter : str
            Name of the parameter of the object.

        Returns
        -------
        value : float
            Value of the parameter.
        """
        if obj == "Components":
            return self.nw.get_comp(label).get_attr(parameter).val
        elif obj == "Connections":
            return self.nw.get_conn(label).get_attr(parameter).val

    def set_params(self, **kwargs):

        if "Connections" in kwargs:
            for c, params in kwargs["Connections"].items():
                self.nw.get_conn(c).set_attr(**params)

        if "Components" in kwargs:
            for c, params in kwargs["Components"].items():
                self.nw.get_comp(c).set_attr(**params)

    def solve_model(self, **kwargs):
        """
        Solve the TESPy model given the the input parameters
        """
        self.set_params(**kwargs)

        self.solved = False
        try:
            self.nw.solve("design")
            if self.nw.res[-1] >= 1e-3 or self.nw.lin_dep:
                self.nw.solve("design", init_only=True, init_path=self.stable)
            else:
                # might need more checks here!
                if (
                        any(self.nw.results["Condenser"]["Q"] > 0)
                        or any(self.nw.results["Desuperheater"]["Q"] > 0)
                        or any(self.nw.results["Turbine"]["P"] > 0)
                        or any(self.nw.results["Pump"]["P"] < 0)
                    ):
                    self.solved = False
                else:
                    self.solved = True
        except ValueError as e:
            self.nw.lin_dep = True
            self.nw.solve("design", init_only=True, init_path=self.stable)

    def get_objective(self, objective=None):
        """
        Get the current objective function evaluation.

        Parameters
        ----------
        objective : str
            Name of the objective function.

        Returns
        -------
        objective_value : float
            Evaluation of the objective function.
        """
        if self.solved:
            if objective == "efficiency":
                return 1 / (
                    self.nw.busses["power"].P.val /
                    self.nw.busses["heat"].P.val
                )
            else:
                msg = f"Objective {objective} not implemented."
                raise NotImplementedError(msg)
        else:
            return np.nan

    # %%[sec_3]

plant = SamplePlant()
plant.get_objective("efficiency")
variables = {
    "Connections": {
        "2": {"p": {"min": 1, "max": 40}},
        "4": {"p": {"min": 1, "max": 40}}
    }
}
constraints = {
    "lower limits": {
        "Connections": {
            "2": {"p": "ref1"}
        },
    },
    "ref1": ["Connections", "4", "p"]
}

optimize = OptimizationProblem(
    plant, variables, constraints, objective="efficiency"
)
# %%[sec_4]
num_ind = 10
num_gen = 100

# for algorithm selection and parametrization please consider the pygmo
# documentation! The number of generations indicated in the algorithm is
# the number of evolutions we undertake within each generation defined in
# num_gen
algo = pg.algorithm(pg.ihs(gen=3, seed=42))
# create starting population
pop = pg.population(pg.problem(optimize), size=num_ind, seed=42)

optimize.run(algo, pop, num_ind, num_gen)
# %%[sec_5]
# To access the results
print(optimize.individuals)
# check pygmo documentation to see, what you can get from the population
pop
# plot the results
import matplotlib.pyplot as plt


# make text reasonably sized
plt.rc("font", **{"size": 18})

fig, ax = plt.subplots(1, figsize=(16, 8))

filter_valid_constraint = optimize.individuals["valid"].values
filter_valid_result = ~np.isnan(optimize.individuals["efficiency"].values)
data = optimize.individuals.loc[filter_valid_constraint & filter_valid_result]

sc = ax.scatter(
    data["Connections-2-p"],
    data["Connections-4-p"],
    c=1 / data["efficiency"] * 100,
    s=100
)
cbar = plt.colorbar(sc)
cbar.set_label("Thermal efficiency in %")

ax.set_axisbelow(True)
ax.set_xlabel("Pressure at connection 2 in bar")
ax.set_ylabel("Pressure at connection 4 in bar")
plt.tight_layout()

fig.savefig("pygmo_optimization.svg")
print(data.loc[data["efficiency"].values == data["efficiency"].min()])
# %%[sec_6]
