# %%[sec_1]
import numpy as np
import pygmo as pg

from tespy.components import CycleCloser
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Condenser
from tespy.components import Desuperheater
from tespy.components import SimpleHeatExchanger
from tespy.components import Merge
from tespy.components import Splitter
from tespy.components import PowerBus
from tespy.components import PowerSink
from tespy.components import PowerSource
from tespy.components import Pump
from tespy.components import Turbine
from tespy.connections import Connection
from tespy.connections import PowerConnection
from tespy.networks import Network

from tespy.tools.optimization import OptimizationProblem


class SamplePlant:
    """Class template for TESPy model usage in optimization module."""
    def __init__(self):
        self._create_network()

    def _create_network(self):

        self.nw = Network()
        self.nw.set_attr(
            p_unit="bar", T_unit="C", h_unit="kJ / kg", iterinfo=False
        )
        # components
        # main cycle
        sg = SimpleHeatExchanger("steam generator")
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

        self.nw.add_conns(c21, c22, c23, c24, c31, c32, c33)

        # cooling water
        c41 = Connection(cwi, "out1", con, "in2", label="41")
        c42 = Connection(con, "out2", cwo, "in1", label="42")

        self.nw.add_conns(c41, c42)

        electricity = PowerBus("electricity bus", num_in=3, num_out=4)
        grid = PowerSink("grid")

        e1 = PowerConnection(hpt, "power", electricity, "power_in1", label="e1")
        e2 = PowerConnection(mpt, "power", electricity, "power_in2", label="e2")
        e3 = PowerConnection(lpt, "power", electricity, "power_in3", label="e3")
        e4 = PowerConnection(electricity, "power_out1", pu1, "power", label="e4")
        e5 = PowerConnection(electricity, "power_out2", pu2, "power", label="e5")
        e6 = PowerConnection(electricity, "power_out3", pu3, "power", label="e6")
        e7 = PowerConnection(electricity, "power_out4", grid, "power", label="e7")

        # heating bus
        sg.set_attr(power_connector_location="inlet")
        heat_source = PowerSource("heat source")

        h1 = PowerConnection(heat_source, "power", sg, "heat", label="h1")

        self.nw.add_conns(e1, e2, e3, e4, e5, e6, e7, h1)

        hpt.set_attr(eta_s=0.9)
        mpt.set_attr(eta_s=0.9)
        lpt.set_attr(eta_s=0.9)

        pu1.set_attr(eta_s=0.8)
        pu2.set_attr(eta_s=0.8)
        pu3.set_attr(eta_s=0.8)

        sg.set_attr(pr=0.92)

        con.set_attr(pr1=1, pr2=0.99)
        fwh1.set_attr(pr1=1, pr2=0.99, ttd_u=5)
        fwh2.set_attr(pr1=1, pr2=0.99, ttd_u=5)
        dsh.set_attr(pr1=0.99, pr2=0.99)

        c1.set_attr(m=200, T=650, p=100, fluid={"water": 1})
        c2.set_attr(p=20)
        c4.set_attr(p=3)
        c6.set_attr(p=0.05)

        c41.set_attr(T=20, p=3, fluid={"INCOMP::Water": 1})
        c42.set_attr(T=28)

        # parametrization
        # components
        self.nw.solve("design")
        con.set_attr(ttd_u=5)
        c6.set_attr(p=None)

        self.nw.solve("design")
        self.stable = "_stable.json"
        self.nw.save(self.stable)
        self.solved = True
        self.nw.print_results()

    def _reset_component_boundary_conditions(self):
        self.nw.get_comp("feed water preheater 1").set_attr(ttd_u=5)
        self.nw.get_comp("feed water preheater 2").set_attr(ttd_u=5)
        self.nw.get_comp("condenser").set_attr(ttd_u=5)
        self.nw.get_comp("feed water pump").set_attr(eta_s=0.8)
        self.nw.get_comp("feed water pump").set_attr(eta_s=0.8)
        self.nw.get_comp("feed water pump").set_attr(eta_s=0.8)
        self.nw.get_comp("high pressure turbine").set_attr(eta_s=0.9)
        self.nw.get_comp("mid pressure turbine").set_attr(eta_s=0.9)
        self.nw.get_comp("low pressure turbine").set_attr(eta_s=0.9)

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
            if self.nw.status != 0:
                self._reset_component_boundary_conditions()
                self.nw.solve("design", init_only=True, init_path=self.stable)
            else:
                self.solved = True
        except ValueError as e:
            self._reset_component_boundary_conditions()
            self.nw.solve("design", init_only=True, init_path=self.stable)

    def get_objectives(self, objective_list):
        """Get the objective values

        Parameters
        ----------
        objective_list : list
            Names of the objectives

        Returns
        -------
        list
            Values of the objectives
        """
        return [self.get_objective(obj) for obj in objective_list]

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
                return (
                    self.nw.get_conn("e7").E.val
                    / self.nw.get_conn("h1").E.val
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
    plant, variables, constraints, objective=["efficiency"], minimize=[False]
)
# %%[sec_4]
num_ind = 40
num_gen = 15

# for algorithm selection and parametrization please consider the pygmo
# documentation! The number of generations indicated in the algorithm is
# the number of evolutions we undertake within each generation defined in
# num_gen
algo = pg.algorithm(pg.ihs(gen=num_gen, seed=42))
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
plt.style.use("dark_background")


# make text reasonably sized
plt.rc("font", **{"size": 18})

fig, ax = plt.subplots(1, figsize=(16, 8))

mask_constraint = optimize.individuals["Connections-2-p>=Connections-4-p"] < 0
mask_objective = ~np.isnan(optimize.individuals["efficiency"].values)
data = optimize.individuals.loc[mask_constraint & mask_objective]

sc = ax.scatter(
    data["Connections-2-p"],
    data["Connections-4-p"],
    c=data["efficiency"] * 100,
    s=100
)
cbar = plt.colorbar(sc)
cbar.set_label("Thermal efficiency in %")

ax.set_axisbelow(True)
ax.set_xlabel("Pressure at connection 2 in bar")
ax.set_ylabel("Pressure at connection 4 in bar")
plt.tight_layout()

fig.savefig("pygmo_optimization_darkmode.svg")
print(data.loc[data["efficiency"].values == data["efficiency"].max()])
# %%[sec_6]
