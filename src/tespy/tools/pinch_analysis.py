# using pina 0.1.1 for pinch related tasks
from pina import PinchAnalyzer, make_stream

# Useful literature: 
# Arpagaus 2019 Hochtemperatur-Wärmepumpen (p.99)
# Kemp 2006 (p.20)


class PinchTools():
    
    def __init__(self):
        self.min_dT = None
        self.temp_shift = None
        self.streams = []
        self.analyzer = None
    
    # including the general functions of pina to expand them later

    # adding cold streams to the used streams, colsd streams have to be heated and form the cold composite curvve (cold CC)
    def add_cold_streams_manually(self, enthalpy_difference:float, T_inlet:float, T_outlet:float):
        # check if the stream has a positive enthalpy flow difference, pina needs positive sign for hot streams, way to check user
        if enthalpy_difference < 0:
            print("Got positive enthaply difference, expected negative enthalpy for cold streams. stream not added")
        else:
            self.streams.append(make_stream(enthalpy_difference,T_inlet,T_outlet))

    # adding cold streams to the used streams, colsd streams have to be heated and form the cold composite curvve (hot CC)
    def add_hot_streams_manually(self,  enthalpy_difference:float, T_inlet:float, T_outlet:float):
        # check if the stream has a negative enthalpy flow difference, pina needs negative sign for cold streams, way to check user
        if enthalpy_difference > 0:
            print("Got positive enthaply difference, expected negative enthalpy for cold streams. stream not added")
        else:
            self.streams.append(make_stream(enthalpy_difference,T_inlet,T_outlet))


    # setting the value for shifting the composit curces (CCs)
    def set_minimum_temperature_difference(self, minimum_temperature_difference:float):
        # minimum temperature difference of streams 
        self.min_dT = minimum_temperature_difference
        # the hot and cold composit curve (CC) are shifted by half the minimum temperature difference to meet at the pinch point
        self.temp_shift = self.min_dT / 2


    # create the pinch analyzer from pina
    def _create_analyzer(self):
        self.analyzer = PinchAnalyzer(self.temp_shift)
        self.analyzer.add_streams(*self.streams) # add a way to add streams later


    # get datapoints of hot composite curve
    def get_hot_cc(self):
        if self.analyzer is None:
            self._create_analyzer()
        [self.hot_cc_data_enthalpy, self.hot_cc_data_temperature] = self.analyzer.hot_composite_curve
    

    # get datapoints of cold composite curve
    def get_cold_cc(self):
        [self.cold_cc_data_enthalpy, self.cold_cc_data_temperature] = self.analyzer.cold_composite_curve


    # get datapoints of shifted hot composite curve
    def get_shifted_hot_cc(self):
        [self.shifted_hot_cc_data_enthalpy, self.shifted_hot_cc_data_temperature] = self.analyzer.shifted_hot_composite_curve


    # get datapoints of shifted cold composite curve
    def get_shifted_cold_cc(self):
        [self.shifted_cold_cc_data_enthalpy, self.shifted_cold_cc_data_temperature] = self.analyzer.shifted_cold_composite_curve


    # get datapoints of grand composite curve
    def get_gcc(self):
        [self.gcc_data_enthalpy, self.gcc_data_shifted_temperature] = self.analyzer.grand_composite_curve


    # add: def get_pinch_point(self):
        # needed to check e.g. pinch rules for heat pump integration automatically

    # add: def plot_cc_diagram(self):
    
    # add: def plot_shifted_cc_diagram(self):

    # add: def plot_gcc_diagram(self):

    # include heat cascades later (for more than one point in time)


    # adding components from tespy models

    # add: def show_heat_pump_in_gcc(self):
    
        # get plot data of heat exchangers at heat sink
        # as used in e.g. the heat exchanger example plots

        # get plot data of heat exchangers at heat source

        # show (only display) by adding the plot data to the gcc

    # add: def add_heat_exchanger_to_streams(self):
        # function to read streams from given heat exchangers, later use to iterate all heat exchangers of network with possibility to exclude streams / heat exchangers
