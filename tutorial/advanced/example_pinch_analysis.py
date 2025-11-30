
from tespy.tools.pinch_analysis import TesypPinchAnalysis


# example Pinch Analysis of a manual workflow without tespy components
Example_Analysis = TesypPinchAnalysis("Example_Process_1")


# setting the minimum temperature difference for the analysis
Example_Analysis.set_minimum_temperature_difference(10)


# add all the streams manually
# Example: Kemp 2007 p. 20
Example_Analysis.add_cold_stream_manually(-230, 20, 135)
Example_Analysis.add_hot_stream_manually(330, 170, 60)
Example_Analysis.add_cold_stream_manually(-240, 80, 140)
Example_Analysis.add_hot_stream_manually(180, 150, 30)
# additional latent streams as shown by Arpagaus 2019 p. 99
Example_Analysis.add_hot_stream_manually(60,85,85)
Example_Analysis.add_cold_stream_manually(-40,130,130)


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