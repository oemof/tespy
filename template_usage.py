from tespy.components import Compressor
from tespy.components import CycleCloser
from tespy.components import HeatExchanger
from tespy.components import Motor
from tespy.components import PowerBus
from tespy.components import PowerSource
from tespy.components import SectionedHeatExchanger
from tespy.components import SimpleHeatExchanger
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Valve
from tespy.connections import Connection
from tespy.connections import PowerConnection
from tespy.models import ModelTemplate


class MyConcreteModel(ModelTemplate):

    def _parameter_lookup(self):
        return {
            "ihx pinch": ["Components", "internal heat exchanger", "td_pinch"],
            "heat": ["Components", "condenser high", "Q"],
            "high stage compressor efficiency": ["Components", "compressor high", "eta_s"],
        }

    def _subcycle_mapping(self):
        return {
            "lower": "a1",
            "upper": "b1"
        }

    def _create_network(self):
        super()._create_network()

        self.nw.units.set_defaults(
            temperature="degC",
            pressure="bar",
            power="kW"
        )

        cc_low = CycleCloser("cc low")
        compressor_low = Compressor("compressor low")
        condenser_low = SectionedHeatExchanger("internal heat exchanger")
        valve_low = Valve("valve low")
        # evaporator_low = SimpleHeatExchanger("evaporator low")
        evaporator_low = SectionedHeatExchanger("evaporator low")
        he_in = Source("source")
        he_out = Sink("sink")

        cc_high = CycleCloser("cc high")
        compressor_high = Compressor("compressor high")
        condenser_high = SimpleHeatExchanger("condenser high")
        valve_high = Valve("valve high")

        a1 = Connection(cc_low, "out1", compressor_low, "in1", label="a1")
        a2 = Connection(compressor_low, "out1", condenser_low, "in1", label="a2")
        a3 = Connection(condenser_low, "out1", valve_low, "in1", label="a3")
        a4 = Connection(valve_low, "out1", evaporator_low, "in2", label="a4")
        a5 = Connection(evaporator_low, "out2", cc_low, "in1", label="a5")

        b1 = Connection(cc_high, "out1", compressor_high, "in1", label="b1")
        b2 = Connection(compressor_high, "out1", condenser_high, "in1", label="b2")
        b3 = Connection(condenser_high, "out1", valve_high, "in1", label="b3")
        b4 = Connection(valve_high, "out1", condenser_low, "in2", label="b4")
        b5 = Connection(condenser_low, "out2", cc_high, "in1", label="b5")

        c1 = Connection(he_in, "out1", evaporator_low, "in1", label="c1")
        c2 = Connection(evaporator_low, "out1", he_out, "in1", label="c2")

        self.nw.add_conns(a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, c1, c2)

        a1.set_attr(fluid={"R600a": 1}, T_dew=0, td_dew=10, m=1)
        a3.set_attr(T_bubble=55, td_bubble=5)
        b1.set_attr(fluid={"R600a": 1}, td_dew=10)
        b3.set_attr(T_bubble=100, td_bubble=2)

        condenser_high.set_attr(dp=0)
        condenser_low.set_attr(dp1=0, dp2=0, td_pinch=5)
        # evaporator_low.set_attr(dp=0)
        evaporator_low.set_attr(dp1=0, dp2=0)
        c1.set_attr(p=1, T=20, fluid={"air": 1})
        c2.set_attr(T=10)

        compressor_low.set_attr(eta_s=0.8)
        compressor_high.set_attr(eta_s=0.8)

        self.nw.solve("design")


model = MyConcreteModel()
print(model.nw)

model.plot_logph_diagram_matplotlib("upper", ".")
model.plot_Ts_diagram_matplotlib("upper", ".")
model.plot_QT_diagram_matplotlib("internal heat exchanger", ".")
model.plot_QT_diagram_matplotlib("evaporator low", ".")
# Sensitivity analysis
param_dict = {
    "ihx pinch": [3, 3, 5, 5, 10, 10],
    "heat": [-100, -200, -100, -200, -100, -200],
    "high stage compressor efficiency": [0.7, 0.7, 0.7, 0.7, 0.7]
}
model.sensitivity_analysis(param_dict=param_dict)
# model.sensitivity_analysis()
