# %%[sec_1]
import numpy as np
from pymoo.algorithms.soo.nonconvex.de import DE  # https://pymoo.org/algorithms/index.html

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
from tespy.components import HeatSource
from tespy.components import Pump
from tespy.components import Turbine
from tespy.connections import Connection
from tespy.connections import HeatConnection
from tespy.connections import PowerConnection
from tespy.models import ModelTemplate


class SamplePlant(ModelTemplate):

    def _create_network(self):
        super()._create_network()

        self.nw.iterinfo = False
        self.nw.units.set_defaults(**{
            "pressure": "bar", "temperature": "degC", "enthalpy": "kJ/kg",
            "pressure_difference": "bar"
        })
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

        # heating
        heat_source = HeatSource("heat source")

        h1 = HeatConnection(heat_source, "heat", sg, "heat", label="h1")

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

        self.nw.solve("design")
        con.set_attr(ttd_u=5)
        c6.set_attr(p=None)

        self.nw.solve("design")
        self._stable_solution = self.nw.save(as_dict=True)
        self._solved = True
        self.nw.print_results()
# %%[sec_2]
    def _parameter_lookup(self):
        return {
            "extraction pressure 1": ["Connections", "2", "p"],
            "extraction pressure 2": ["Connections", "4", "p"],
            "hpt power": ["Components", "high pressure turbine", "P"],
            "hpt pressure ratio": ["Components", "high pressure turbine", "pr"],
            "efficiency": {"get": self.calc_efficiency},
        }

    def calc_efficiency(self):
        if self._solved:
            return (
                self.nw.get_conn("e7").E.val
                / self.nw.get_conn("h1").E.val
            )
        return np.nan

    def solve_model(self, **kwargs):
        self.solve_model_design(**kwargs)
# %%[sec_3]
plant = SamplePlant()
plant.get_parameter("efficiency")

num_evo = 20
# %%[sec_4]
import os

if os.getenv("GITHUB_ACTIONS") == "true" or "PYTEST_CURRENT_TEST" in os.environ:
    num_evo = 2
# %%[sec_5]
algorithm = DE(pop_size=20)

result = plant.optimize(
    algorithm=algorithm,
    termination=("n_gen", num_evo),
    variables={
        "extraction pressure 1": {"min": 1, "max": 40},
        "extraction pressure 2": {"min": 1, "max": 40},
    },
    constraints={
        "extraction pressure 1": {"min": "extraction pressure 2"},
    },
    objective=["efficiency"],
    minimize_flags=[False],
    kpi=["hpt power", "hpt pressure ratio"],
)
# %%[sec_6]
print(result)

# plot the results
import matplotlib.pyplot as plt

plt.rc("font", **{"size": 18})

fig, ax = plt.subplots(1, figsize=(16, 8))

mask_constraint = result["extraction pressure 1>=extraction pressure 2"] < 0
mask_objective = ~np.isnan(result["efficiency"].values)
data = result.loc[mask_constraint & mask_objective]

sc = ax.scatter(
    data["extraction pressure 1"],
    data["extraction pressure 2"],
    c=data["efficiency"] * 100,
    s=100,
)
best = data.loc[data["efficiency"].values == data["efficiency"].max()]
ax.scatter(
    x=best["extraction pressure 1"],
    y=best["extraction pressure 2"],
    c="red",
    marker="x",
)
cbar = plt.colorbar(sc)
cbar.set_label("Thermal efficiency in %")

ax.set_axisbelow(True)
ax.set_xlabel("Extraction pressure 1 in bar")
ax.set_ylabel("Extraction pressure 2 in bar")
plt.tight_layout()

fig.savefig("optimization_result.svg")
print(best)
# %%[sec_7]
