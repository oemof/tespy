
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
# in that case only the condensation is part of the conndenser forming a horizontal line in the GCC as
# the most simple example case
c1.set_attr(m=0.1, p=5, x=1)
c3.set_attr(p=20, x=1)
c4.set_attr(x=0)

# solve design
nw.solve("design")

# reference heat pump components for plotting in the GCC
Example_Analysis.show_heat_pump_in_gcc(condenser=condenser,evaporator=evaporator)