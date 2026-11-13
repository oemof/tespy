# %%[sec_1]
from CoolProp.CoolProp import PropsSI

from tespy.components import Absorber
from tespy.components import CycleCloser
from tespy.components import Desorber
from tespy.components import HeatExchanger
from tespy.components import Pump
from tespy.components import SimpleHeatExchanger
from tespy.components import Valve
from tespy.connections import Connection
from tespy.connections import SolutionConnection
from tespy.networks import Network

nw = Network(iterinfo=False)
nw.units.set_defaults(temperature="degC", pressure="bar", heat="kW")

generator = Desorber("generator")
condenser = SimpleHeatExchanger("condenser")
refrigerant_valve = Valve("refrigerant valve")
evaporator = SimpleHeatExchanger("evaporator")
absorber = Absorber("absorber")
closer = CycleCloser("solution cycle closer")
pump = Pump("solution pump")
shx = HeatExchanger("solution heat exchanger")
solution_valve = Valve("solution valve")
# %%[sec_2]
s0 = SolutionConnection(
    absorber, "out1", closer, "in1", label="s0", mixing_rule="libr-water"
)
s1 = SolutionConnection(closer, "out1", pump, "in1", label="s1")
s2 = SolutionConnection(pump, "out1", shx, "in2", label="s2")
s3 = SolutionConnection(shx, "out2", generator, "in1", label="s3")
s4 = SolutionConnection(generator, "out1", shx, "in1", label="s4")
s5 = SolutionConnection(shx, "out1", solution_valve, "in1", label="s5")
s6 = SolutionConnection(solution_valve, "out1", absorber, "in1", label="s6")
r1 = Connection(generator, "out2", condenser, "in1", label="r1")
r2 = Connection(condenser, "out1", refrigerant_valve, "in1", label="r2")
r3 = Connection(refrigerant_valve, "out1", evaporator, "in1", label="r3")
r4 = Connection(evaporator, "out1", absorber, "in2", label="r4")
nw.add_conns(s0, s1, s2, s3, s4, s5, s6, r1, r2, r3, r4)
# %%[sec_3]
condenser.set_attr(dp=0)
evaporator.set_attr(dp=0, Q=100)
pump.set_attr(eta_s=0.7)
shx.set_attr(dp1=0, dp2=0, ttd_l=10)

r1.set_attr(fluid={"H2O": 1})
r2.set_attr(x=0, T=45)
r4.set_attr(x=1, T=5)
s0.set_attr(T=40)
s1.set_attr(fluid_balance=True)
s4.set_attr(T=130)

nw.solve("design")
nw.print_results()
COP = evaporator.Q.val / generator.Q.val
print(f"COP: {COP:.3f}")
# %%[sec_4]
generator.set_attr(vapor_state="outlet-equilibrium")
nw.solve("design")
print(f"COP with the vapor in equilibrium with the leaving solution: "
      f"{evaporator.Q.val / generator.Q.val:.3f}")
# %%[sec_5]
from tespy.tools.fluid_properties.solutions import register_solution_pair
from tespy.tools.fluid_properties.solutions.libr_water import LiBrWater


class ShiftedLiBrWater(LiBrWater):
    """LiBr-water with a boiling point elevation 1 K higher than Patek and
    Klomfar, as an example of a user defined working pair."""

    name = "libr-water-shifted"

    def T_bubble(self, p, w):
        return super().T_bubble(p, w) + 1

    def p_bubble(self, T, w):
        return super().p_bubble(T - 1, w)

    def w_bubble(self, p, T):
        return super().w_bubble(p, T - 1)


register_solution_pair(ShiftedLiBrWater())

pair = ShiftedLiBrWater()
p = PropsSI("P", "T", 278.15, "Q", 0, "water")
print(f"bubble point of 55 % LiBr at 5 degC water saturation pressure: "
      f"{pair.T_bubble(p, 0.55) - 273.15:.2f} degC")
# %%[sec_6]
