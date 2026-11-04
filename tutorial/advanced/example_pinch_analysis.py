
from tespy.tools.pinch_analysis import TesypPinchAnalysis
# setting up the heat pump system
from tespy.networks import Network
from tespy.connections import Connection
from tespy.components import SimpleHeatExchanger, Valve, Compressor, CycleCloser


# MARK: Example 1

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


# network
nw = Network()

# taken from example heat pump
nw.units.set_defaults(
    temperature="degC", pressure="bar", enthalpy="kJ/kg", heat="kW", power="kW", pressure_difference = "bar"
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
c3.set_attr(p=21, x=1)
c4.set_attr(x=0)

# solve design
nw.solve("design")

# reference heat pump components for plotting in the GCC
Example_Analysis.show_heat_pump_in_gcc(condenser=condenser,evaporator=evaporator)


# MARK: Example 2
# second example using a moving boundary heat exchanger as the evaporator and a sectioned heat exchanger as the condenser
# as the heat exchangers now include the internal H,T-Data, the desuperheater is not included anymore.
# However, the connection numbers are kept to clarifiy the positions.

# setting up the heat pump system
from tespy.components import  MovingBoundaryHeatExchanger, SectionedHeatExchanger, Sink, Source

# network
nw_2 = Network()

# taken from example heat pump
nw_2.units.set_defaults(
    temperature="degC", pressure="bar", enthalpy="kJ/kg", heat="kW", power="kW", pressure_difference = "bar"
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
c7_2.set_attr(fluid={"Water": 1}, T=55, p=1)
c8_2.set_attr(T=60)
c9_2.set_attr(fluid={"Water": 1}, T=10, p=1)
c10_2.set_attr(T=5)

# set up parameters to show specific case
c1_2.set_attr(m=0.1, x=1)
c4_2.set_attr(x=0)

# solve design
nw_2.solve("design")

# integrating a tespy model of a heat pump into a given process

# example Pinch Analysis of a manual workflow without tespy components
Example_Analysis2 = TesypPinchAnalysis("Example_Process_2")

# setting the minimum temperature difference for the analysis
Example_Analysis2.set_minimum_temperature_difference(10)

# add all the streams manually
# Example: Kemp 2007 p. 20, reduced temperature by 50 degC to fit the heat pump example
Example_Analysis2.add_cold_stream_manually(-230, 20-50, 135-50)
Example_Analysis2.add_hot_stream_manually(330, 170-50, 60-50)
Example_Analysis2.add_cold_stream_manually(-240, 80-50, 140-50)
Example_Analysis2.add_hot_stream_manually(180, 150-50, 30-50)
# additional latent streams as shown by Arpagaus 2019 p. 99
Example_Analysis2.add_hot_stream_manually(60,35,35)
Example_Analysis2.add_cold_stream_manually(-40,80,80)

# reference heat pump components for plotting in the GCC of the same pinch analysis
Example_Analysis2.show_heat_pump_in_gcc(condenser=condenser_2,evaporator=evaporator_2)



# MARK: Example 3
# additional example to test HeatExchanger in pinch analysis
from tespy.components import  HeatExchanger

# network
nw_3 = Network()

# taken from example heat pump
nw_3.units.set_defaults(
    temperature="degC", pressure="bar", enthalpy="kJ/kg", heat="kW", power="kW", pressure_difference = "bar"
)

# components
condenser_3 = HeatExchanger("Condenser_3")
evaporator_3 = HeatExchanger("Evaporator_3")
expansion_valve_3 = Valve("Expansion Valve_3")
compressor_3 = Compressor("Compressor_3")
cycle_closer_3 = CycleCloser("Cycle Closer_3")
# Sinks and Sources for the secondary media
HeatSinkIn = Source("HeatSinkIn")
HeatSinkOut = Sink("HeatSinkOut")
HeatSourceIn = Source("HeatSourceIn")
HeatSourceOut = Sink("HeatSourceOut")

# connections
c1_3 = Connection(evaporator_3, "out2", compressor_3,"in1", label = "connection 1_3")
c2_3 = Connection(compressor_3, "out1", condenser_3, "in1", "connection 2_3") 
c4_3 = Connection(condenser_3, "out1", expansion_valve_3, "in1", label = "connection 4_3")
c5_3 = Connection(expansion_valve_3, "out1",cycle_closer_3, "in1", label = "connection 5_3")
c6_3 = Connection(cycle_closer_3, "out1", evaporator_3, "in2", label = "connection 6_3")
nw_3.add_conns(c1_3,c2_3,c4_3,c5_3,c6_3)

# add the connections for secondary media
c7_3 = Connection(HeatSinkIn, "out1", condenser_3, "in2", label = "connection 7_3")
c8_3 = Connection(condenser_3, "out2", HeatSinkOut, "in1", label = "connection 8_3")
c9_3 = Connection(HeatSourceIn, "out1", evaporator_3, "in1", label = "connection 9_3")
c10_3 = Connection(evaporator_3, "out1", HeatSourceOut, "in1", label = "connection 10_3")
nw_3.add_conns(c7_3,c8_3,c9_3,c10_3)

# set up general parameters of heat pump
compressor_3.set_attr(eta_s = 0.7)
condenser_3.set_attr(dp1=0, dp2=0, ttd_l=2)
evaporator_3.set_attr(dp1=0, dp2=0, ttd_l=2)
c1_3.set_attr(fluid={"R290": 1})
# media, temperatures and pressure of heat source and sink
c7_3.set_attr(fluid={"Water": 1}, T=55, p=1)
c8_3.set_attr(T=60)
c9_3.set_attr(fluid={"Water": 1}, T=10, p=1)
c10_3.set_attr(T=5)

# set up parameters to show specific case
c1_3.set_attr(m=0.1, x=1)
c4_3.set_attr(x=0)

# solve design
nw_3.solve("design")

# integrating a tespy model of a heat pump into a given process

# example Pinch Analysis of a manual workflow without tespy components
Example_Analysis3 = TesypPinchAnalysis("Example_Process_3")

# setting the minimum temperature difference for the analysis
Example_Analysis3.set_minimum_temperature_difference(10)

# add all the streams manually
# Example: Kemp 2007 p. 20, reduced temperature by 50 degC to fit the heat pump example
Example_Analysis3.add_cold_stream_manually(-230, 20-50, 135-50)
Example_Analysis3.add_hot_stream_manually(330, 170-50, 60-50)
Example_Analysis3.add_cold_stream_manually(-240, 80-50, 140-50)
Example_Analysis3.add_hot_stream_manually(180, 150-50, 30-50)
# additional latent streams as shown by Arpagaus 2019 p. 99
Example_Analysis3.add_hot_stream_manually(60,35,35)
Example_Analysis3.add_cold_stream_manually(-40,80,80)

# reference heat pump components for plotting in the GCC of the same pinch analysis
Example_Analysis3.show_heat_pump_in_gcc(condenser=condenser_3,evaporator=evaporator_3)



# MARK: Example 4
# testing: ParallelFlowHeatExchanger, Desuperheater
# additional example to test HeatExchanger in pinch analysis
from tespy.components import  ParallelFlowHeatExchanger, Desuperheater

# network
nw_4 = Network()

# taken from example heat pump
nw_4.units.set_defaults(
    temperature="degC", pressure="bar", enthalpy="kJ/kg", heat="kW", power="kW", pressure_difference = "bar"
)

# components
condenser_4 = Desuperheater("Condenser_4")
evaporator_4 = ParallelFlowHeatExchanger("Evaporator_4")
expansion_valve_4 = Valve("Expansion Valve_4")
compressor_4 = Compressor("Compressor_4")
cycle_closer_4 = CycleCloser("Cycle Closer_4")
# Sinks and Sources for the secondary media
HeatSinkIn = Source("HeatSinkIn")
HeatSinkOut = Sink("HeatSinkOut")
HeatSourceIn = Source("HeatSourceIn")
HeatSourceOut = Sink("HeatSourceOut")

# connections
c1_4 = Connection(evaporator_4, "out2", compressor_4,"in1", label = "connection 1_4")
c2_4 = Connection(compressor_4, "out1", condenser_4, "in1", "connection 2_4") 
c4_4 = Connection(condenser_4, "out1", expansion_valve_4, "in1", label = "connection 4_4")
c5_4 = Connection(expansion_valve_4, "out1",cycle_closer_4, "in1", label = "connection 5_4")
c6_4 = Connection(cycle_closer_4, "out1", evaporator_4, "in2", label = "connection 6_4")
nw_4.add_conns(c1_4,c2_4,c4_4,c5_4,c6_4)

# add the connections for secondary media
c7_4 = Connection(HeatSinkIn, "out1", condenser_4, "in2", label = "connection 7_4")
c8_4 = Connection(condenser_4, "out2", HeatSinkOut, "in1", label = "connection 8_4")
c9_4 = Connection(HeatSourceIn, "out1", evaporator_4, "in1", label = "connection 9_4")
c10_4 = Connection(evaporator_4, "out1", HeatSourceOut, "in1", label = "connection 10_4")
nw_4.add_conns(c7_4,c8_4,c9_4,c10_4)

# set up general parameters of heat pump
compressor_4.set_attr(eta_s = 0.7)
condenser_4.set_attr(dp1=0, dp2=0, ttd_l=2)
evaporator_4.set_attr(dp1=0, dp2=0, ttd_l=2)
c1_4.set_attr(fluid={"R290": 1})
# media, temperatures and pressure of heat source and sink
c7_4.set_attr(fluid={"Water": 1}, T=55, p=1)
c8_4.set_attr(T=60)
c9_4.set_attr(fluid={"Water": 1}, T=10, p=1)
c10_4.set_attr(T=5)

# set up parameters to show specific case
c1_4.set_attr(m=0.1) # no x specified, as the desuperheater is used
c4_4.set_attr(x=0)

# solve design
nw_4.solve("design")

# integrating a tespy model of a heat pump into a given process

# example Pinch Analysis of a manual workflow without tespy components
Example_Analysis4 = TesypPinchAnalysis("Example_Process_4")

# setting the minimum temperature difference for the analysis
Example_Analysis4.set_minimum_temperature_difference(10)

# add all the streams manually
# Example: Kemp 2007 p. 20, reduced temperature by 50 degC to fit the heat pump example
Example_Analysis4.add_cold_stream_manually(-230, 20-50, 135-50)
Example_Analysis4.add_hot_stream_manually(330, 170-50, 60-50)
Example_Analysis4.add_cold_stream_manually(-240, 80-50, 140-50)
Example_Analysis4.add_hot_stream_manually(180, 150-50, 30-50)
# additional latent streams as shown by Arpagaus 2019 p. 99
Example_Analysis4.add_hot_stream_manually(60,35,35)
Example_Analysis4.add_cold_stream_manually(-40,80,80)

# reference heat pump components for plotting in the GCC of the same pinch analysis
Example_Analysis4.show_heat_pump_in_gcc(condenser=condenser_3,evaporator=evaporator_3)



# MARK: Example 5
# testing: Condenser
# additional example to test HeatExchanger in pinch analysis
from tespy.components import Condenser

# network
nw_5 = Network()

# taken from example heat pump
nw_5.units.set_defaults(
    temperature="degC", pressure="bar", enthalpy="kJ/kg", heat="kW", power="kW", pressure_difference = "bar"
)

# components
condenser_5 = Condenser("Condenser_5")
evaporator_5 = HeatExchanger("Evaporator_5")
expansion_valve_5 = Valve("Expansion Valve_5")
compressor_5 = Compressor("Compressor_5")
cycle_closer_5 = CycleCloser("Cycle Closer_5")
# Sinks and Sources for the secondary media
HeatSinkIn = Source("HeatSinkIn")
HeatSinkOut = Sink("HeatSinkOut")
HeatSourceIn = Source("HeatSourceIn")
HeatSourceOut = Sink("HeatSourceOut")

# connections
c1_5 = Connection(evaporator_5, "out2", compressor_5,"in1", label = "connection 1_5")
c2_5 = Connection(compressor_5, "out1", condenser_5, "in1", "connection 2_5") 
c4_5 = Connection(condenser_5, "out1", expansion_valve_5, "in1", label = "connection 4_5")
c5_5 = Connection(expansion_valve_5, "out1",cycle_closer_5, "in1", label = "connection 5_5")
c6_5 = Connection(cycle_closer_5, "out1", evaporator_5, "in2", label = "connection 6_5")
nw_5.add_conns(c1_5,c2_5,c4_5,c5_5,c6_5)

# add the connections for secondary media
c7_5 = Connection(HeatSinkIn, "out1", condenser_5, "in2", label = "connection 7_5")
c8_5 = Connection(condenser_5, "out2", HeatSinkOut, "in1", label = "connection 8_5")
c9_5 = Connection(HeatSourceIn, "out1", evaporator_5, "in1", label = "connection 9_5")
c10_5 = Connection(evaporator_5, "out1", HeatSourceOut, "in1", label = "connection 10_5")
nw_5.add_conns(c7_5,c8_5,c9_5,c10_5)

# set up general parameters of heat pump
compressor_5.set_attr(eta_s = 0.7)
condenser_5.set_attr(dp1=0, dp2=0, ttd_l=6)
evaporator_5.set_attr(dp1=0, dp2=0, ttd_l=2)
c1_5.set_attr(fluid={"R290": 1})
# media, temperatures and pressure of heat source and sink
c7_5.set_attr(fluid={"Water": 1}, T=55, p=1)
c8_5.set_attr(T=60)
c9_5.set_attr(fluid={"Water": 1}, T=10, p=1)
c10_5.set_attr(T=5)

# set up parameters to show specific case
c1_5.set_attr(m=0.1, x=1)


# solve design
nw_5.solve("design")

# integrating a tespy model of a heat pump into a given process

# example Pinch Analysis of a manual workflow without tespy components
Example_Analysis5 = TesypPinchAnalysis("Example_Process_5")

# setting the minimum temperature difference for the analysis
Example_Analysis5.set_minimum_temperature_difference(10)

# add all the streams manually
# Example: Kemp 2007 p. 20, reduced temperature by 50 degC to fit the heat pump example
Example_Analysis5.add_cold_stream_manually(-230, 20-50, 135-50)
Example_Analysis5.add_hot_stream_manually(330, 170-50, 60-50)
Example_Analysis5.add_cold_stream_manually(-240, 80-50, 140-50)
Example_Analysis5.add_hot_stream_manually(180, 150-50, 30-50)
# additional latent streams as shown by Arpagaus 2019 p. 99
Example_Analysis5.add_hot_stream_manually(60,35,35)
Example_Analysis5.add_cold_stream_manually(-40,80,80)

# reference heat pump components for plotting in the GCC of the same pinch analysis
Example_Analysis5.show_heat_pump_in_gcc(condenser=condenser_3,evaporator=evaporator_3)