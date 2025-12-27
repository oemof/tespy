
from tespy.tools.pinch_analysis import TesypPinchAnalysis


# integrating a tespy model of a heat pump into a given process

# example Pinch Analysis of a manual workflow without tespy components
Example_Analysis = TesypPinchAnalysis("Example_Process_1")

# setting the minimum temperature difference for the analysis
Example_Analysis.set_minimum_temperature_difference(10)

# add all the streams manually
# Example: Kemp 2007 p. 20, reduced temperature by 50 degC to fit the heat pump example
Example_Analysis.add_cold_stream_manually(-230, 20-50, 135-50)
Example_Analysis.add_hot_stream_manually(330, 170-50, 60-50)
Example_Analysis.add_cold_stream_manually(-240, 80-50, 140-50)
Example_Analysis.add_hot_stream_manually(180, 150-50, 30-50)
# additional latent streams as shown by Arpagaus 2019 p. 99
Example_Analysis.add_hot_stream_manually(60,35,35)
Example_Analysis.add_cold_stream_manually(-40,80,80)

# with all streams added, get the results as in pina
cold_cc_data = Example_Analysis.get_cold_cc()
hot_cc_data = Example_Analysis.get_hot_cc()
shifted_cold_cc_data = Example_Analysis.get_shifted_cold_cc()
shifted_hot_cc_data = Example_Analysis.get_shifted_hot_cc()
gcc_data = Example_Analysis.get_gcc()

# this data can be used for other aspects like combining with tespy heat exchangers
# to plot the results use the following functions
# the composite curves
Example_Analysis.plot_cc_diagram()
# the shifted composite curves touching in the pinch point
Example_Analysis.plot_shifted_cc_diagram()
# the grand composite curve
Example_Analysis.plot_gcc_diagram()


# setting up the heat pump system
from tespy.networks import Network
from tespy.connections import Connection
from tespy.components import SimpleHeatExchanger, Valve, Compressor, CycleCloser

# network
nw = Network()

# taken from example heat pump
nw.units.set_defaults(
    temperature="degC", pressure="bar", enthalpy="kJ/kg", heat="kW", power="kW"
)

# components
condenser = SimpleHeatExchanger("Condenser")
evaporator = SimpleHeatExchanger("Evaporator")
desuperheater = SimpleHeatExchanger("Desuperheater")
expansion_valve = Valve("Expansion Valve")
compressor = Compressor("Compressor")
cycle_closer = CycleCloser("Cycle Closer")

# connections
c1 = Connection(evaporator, "out1", compressor,"in1", label = "connection 1")
c2 = Connection(compressor, "out1", desuperheater, "in1", "connection 2") 
c3 = Connection(desuperheater, "out1", condenser, "in1", label = "connection 3")
c4 = Connection(condenser, "out1", expansion_valve, "in1", label = "connection 4")
c5 = Connection(expansion_valve, "out1",cycle_closer, "in1", label = "connection 5")
c6 = Connection(cycle_closer, "out1", evaporator, "in1", label = "connection 6")
nw.add_conns(c1,c2,c3,c4,c5,c6)

# set up general parameters of heat pump
compressor.set_attr(eta_s = 0.7)
condenser.set_attr(dp=0)
desuperheater.set_attr(dp=0)
evaporator.set_attr(dp=0)
c1.set_attr(fluid={"R290": 1})

# set up parameters to show specific case, the desuperheater reduces the temperature to the dewline
# in that case only the condensation is part of the condenser forming a horizontal line in the GCC as
# the most simple example case
c1.set_attr(m=0.1, p=5, x=1)
c3.set_attr(p=20, x=1)
c4.set_attr(x=0)

# solve design
nw.solve("design")

# reference heat pump components for plotting in the GCC
Example_Analysis.show_heat_pump_in_gcc(condenser=condenser,evaporator=evaporator)



# second example using a moving boundary heat exchanger as the evaporator and a sectioned heat exchanger as the condenser
# as the heat exchangers now include the internal H,T-Data, the desuperheater is not included anymore.
# However, the connection numbers are kept to clarifiy the positions.

# setting up the heat pump system
from tespy.components import  MovingBoundaryHeatExchanger, SectionedHeatExchanger, Sink, Source

# network
nw_2 = Network()

# taken from example heat pump
nw_2.units.set_defaults(
    temperature="degC", pressure="bar", enthalpy="kJ/kg", heat="kW", power="kW"
)

# components
condenser_2 = SectionedHeatExchanger("Condenser_2")
evaporator_2 = MovingBoundaryHeatExchanger("Evaporator_2")
expansion_valve_2 = Valve("Expansion Valve_2")
compressor_2 = Compressor("Compressor_2")
cycle_closer_2 = CycleCloser("Cycle Closer_2")
# Sinks and Sources for the secondary media
HeatSinkIn = Source("HeatSinkIn")
HeatSinkOut = Sink("HeatSinkOut")
HeatSourceIn = Source("HeatSourceIn")
HeatSourceOut = Sink("HeatSourceOut")

# connections
c1_2 = Connection(evaporator_2, "out2", compressor_2,"in1", label = "connection 1_2")
c2_2 = Connection(compressor_2, "out1", condenser_2, "in1", "connection 2_2") 
c4_2 = Connection(condenser_2, "out1", expansion_valve_2, "in1", label = "connection 4_2")
c5_2 = Connection(expansion_valve_2, "out1",cycle_closer_2, "in1", label = "connection 5_2")
c6_2 = Connection(cycle_closer_2, "out1", evaporator_2, "in2", label = "connection 6_2")
nw_2.add_conns(c1_2,c2_2,c4_2,c5_2,c6_2)

# add the connections for secondary media
c7_2 = Connection(HeatSinkIn, "out1", condenser_2, "in2", label = "connection 7_2")
c8_2 = Connection(condenser_2, "out2", HeatSinkOut, "in1", label = "connection 8_2")
c9_2 = Connection(HeatSourceIn, "out1", evaporator_2, "in1", label = "connection 9_2")
c10_2 = Connection(evaporator_2, "out1", HeatSourceOut, "in1", label = "connection 10_2")
nw_2.add_conns(c7_2,c8_2,c9_2,c10_2)

# set up general parameters of heat pump
compressor_2.set_attr(eta_s = 0.7)
condenser_2.set_attr(dp1=0, dp2=0, td_pinch=2)
evaporator_2.set_attr(dp1=0, dp2=0, td_pinch=2)
c1_2.set_attr(fluid={"R290": 1})
# media, temperatures and pressure of heat source and sink
c7_2.set_attr(fluid={"Water": 1}, T=50, p=1)
c8_2.set_attr(T=55)
c9_2.set_attr(fluid={"Water": 1}, T=10, p=1)
c10_2.set_attr(T=5)

# set up parameters to show specific case
c1_2.set_attr(m=0.1, x=1)
c4_2.set_attr(x=0)

# solve design
nw_2.solve("design")

# reference heat pump components for plotting in the GCC of the same pinch analysis
Example_Analysis.show_heat_pump_in_gcc(condenser=condenser_2,evaporator=evaporator_2)