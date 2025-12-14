# using pina 0.1.1 for pinch related tasks
from pina import PinchAnalyzer, make_stream
import matplotlib.pyplot as plt

# Useful literature: 
# Arpagaus 2019 Hochtemperatur-Wärmepumpen (ISBN 978-3-8007-4550-0)
# Kemp 2006 (p.20)
# Walden et al. 2023 (https://doi.org/10.1016/j.apenergy.2023.121933)


class TesypPinchAnalysis():
    
    def __init__(self, label:str):
        self.min_dT = None
        self.temp_shift = None
        self.streams = []
        self.analyzer = None
        self.label = label

    
    # including the general functions of pina to expand them later


    # setting the value for shifting the composit curces (CCs)
    def set_minimum_temperature_difference(self, minimum_temperature_difference:float):
        # minimum temperature difference of streams 
        self.min_dT = minimum_temperature_difference
        # the hot and cold composit curve (CC) are shifted by half the minimum temperature difference to meet at the pinch point
        self.temp_shift = self.min_dT / 2


    # adding cold streams to the used streams, colsd streams have to be heated and form the cold composite curvve (cold CC)
    def add_cold_stream_manually(self, enthalpy_difference:float, T_inlet:float, T_outlet:float):
        # check if the stream has a positive enthalpy flow difference, pina needs negative sign for hot streams, way to check user
        if enthalpy_difference > 0:
            print("Got positive enthaply difference, expected negative enthalpy for cold streams. stream not added")
        else:
            self.streams.append(make_stream(enthalpy_difference,T_inlet,T_outlet))


    # adding cold streams to the used streams, colsd streams have to be heated and form the cold composite curvve (hot CC)
    def add_hot_stream_manually(self,  enthalpy_difference:float, T_inlet:float, T_outlet:float):
        # check if the stream has a negative enthalpy flow difference, pina needs positive sign for cold streams, way to check user
        if enthalpy_difference < 0:
            print("Got positive enthaply difference, expected negative enthalpy for cold streams. stream not added")
        else:
            self.streams.append(make_stream(enthalpy_difference,T_inlet,T_outlet))


    # add: functions to read lists for hot and cold later


    # add: functions to remove streams later


    # create the pinch analyzer from pina
    def _create_analyzer(self):
        # check if necessary input is set
        if self.temp_shift is not None:
            self.analyzer = PinchAnalyzer(self.temp_shift)
        else:
            print("set minimum temperature difference first")
            return
        self.analyzer.add_streams(*self.streams) # add a way to add streams later
        self._get_analyzer_data()


    # get datapoints of hot composite curve
    def get_hot_cc(self):
        if self.analyzer is None:
            self._create_analyzer()
        [self.hot_cc_data_enthalpy, self.hot_cc_data_temperature] = self.analyzer.hot_composite_curve
    

    # get datapoints of cold composite curve
    def get_cold_cc(self):
        if self.analyzer is None:
            self._create_analyzer()
        [self.cold_cc_data_enthalpy, self.cold_cc_data_temperature] = self.analyzer.cold_composite_curve


    # get datapoints of shifted hot composite curve
    def get_shifted_hot_cc(self):
        if self.analyzer is None:
            self._create_analyzer()
        [self.shifted_hot_cc_data_enthalpy, self.shifted_hot_cc_data_temperature] = self.analyzer.shifted_hot_composite_curve


    # get datapoints of shifted cold composite curve
    def get_shifted_cold_cc(self):
        if self.analyzer is None:
            self._create_analyzer()
        [self.shifted_cold_cc_data_enthalpy, self.shifted_cold_cc_data_temperature] = self.analyzer.shifted_cold_composite_curve


    # get datapoints of grand composite curve
    def get_gcc(self):
        if self.analyzer is None:
            self._create_analyzer()
        [self.gcc_data_enthalpy, self.gcc_data_shifted_temperature] = self.analyzer.grand_composite_curve


    def _get_analyzer_data(self):  
        # get additional analysis data from pina
        self.T_pinch = self.analyzer.pinch_temps[0]
        self.cold_utility = self.analyzer.cold_utility_target
        self.hot_utility = self.analyzer.hot_utility_target
        self.heat_recovery = self.analyzer.heat_recovery_target

        # needed to check e.g. pinch rules for heat pump integration automatically
           

    def plot_cc_diagram(self, save_fig:bool = True, show_fig:bool = False, return_fig:bool = False):
        fig, ax = plt.subplots()
        # activate minor ticks
        ax.minorticks_on()
        # plot subplots with same axes limits
        ax.plot(self.hot_cc_data_enthalpy, self.hot_cc_data_temperature, color = "red")
        ax.plot(self.cold_cc_data_enthalpy, self.cold_cc_data_temperature, color = "blue")
        # add visualization of sections
        # get minimum temperature
        T_min_CC = min(min(self.hot_cc_data_temperature), min(self.cold_cc_data_temperature))
        ax.plot([0,self.cold_utility,self.heat_recovery+self.cold_utility,
                 self.heat_recovery+self.cold_utility+self.hot_utility],[T_min_CC-5]*4, "o-", color="black")
        # set x scale on lowest subplot
        ax.tick_params(axis='x', which='major', labelsize=10, rotation = 0)
        # set x label
        ax.set_xlabel("$\dot{H}$ [kW]", loc='center', fontsize="10")
        ax.xaxis.label.set_color("black")
        # set y label
        ax.set_ylabel("T [°C]", color = "black")
        # set grid
        ax.grid(visible=True, which='major', color='lightgrey', linewidth = 0.5)
        ax.grid(visible=True, which='minor', color='lightgrey', linestyle='dotted', linewidth = 0.5)
        # set aspect ratio of subplot
        ax.set_box_aspect(1)
        ax.set_title(f"Composite Curves of \"{self.label}\"",color = "black", fontsize = 10)
        
        if save_fig:
            fig.savefig(f"Composite_Curves_{self.label}.svg")
        if show_fig:
            fig.show()
        if return_fig:
            return fig
    

    def plot_shifted_cc_diagram(self, save_fig:bool = True, show_fig:bool = False, return_fig:bool = False):
        fig, ax = plt.subplots()
        # activate minor ticks
        ax.minorticks_on()
        # plot subplots with same axes limits
        ax.plot(self.shifted_hot_cc_data_enthalpy, self.shifted_hot_cc_data_temperature, color = "red")
        ax.plot(self.shifted_cold_cc_data_enthalpy, self.shifted_cold_cc_data_temperature, color = "blue")
        # add visualization of sections
        # get minimum temperature
        T_min_shifted_CC = min(min(self.shifted_hot_cc_data_temperature), min(self.shifted_cold_cc_data_temperature))
        ax.plot([0,self.cold_utility,self.heat_recovery+self.cold_utility,
                 self.heat_recovery+self.cold_utility+self.hot_utility],[T_min_shifted_CC-5]*4, "o-", color="black")
        # set x scale on lowest subplot
        ax.tick_params(axis='x', which='major', labelsize=10, rotation = 0)
        # set x label
        ax.set_xlabel("$\dot{H}$ [kW]", loc='center', fontsize="10")
        ax.xaxis.label.set_color("black")
        # set y label
        ax.set_ylabel("shifted T* [°C]", color = "black")
        # set grid
        ax.grid(visible=True, which='major', color='lightgrey', linewidth = 0.5)
        ax.grid(visible=True, which='minor', color='lightgrey', linestyle='dotted', linewidth = 0.5)
        # set aspect ratio of subplot
        ax.set_box_aspect(1)
        ax.set_title(f"Shifted Composite Curves of \"{self.label}\"",color = "black", fontsize = 10)
        
        if save_fig:
            fig.savefig(f"Shifted_Composite_Curves_{self.label}.svg")
        if show_fig:
            fig.show()
        if return_fig:
            return fig


    def plot_gcc_diagram(self, save_fig:bool = True, show_fig:bool = False, return_fig:bool = False):
        # plot
        fig, ax = plt.subplots()
        # activate minor ticks
        ax.minorticks_on()
        # plot subplots with same axes limits
        ax.plot(self.gcc_data_enthalpy, self.gcc_data_shifted_temperature, color = "black")
        # add pinch point
        ax.plot(0, self.T_pinch, "o", color = "black")
        # set x scale on lowest subplot
        ax.tick_params(axis='x', which='major', labelsize=10, rotation = 0)
        # set x label
        ax.set_xlabel("$\Delta\dot{H}$ [kW]", loc='center', fontsize="10")
        ax.xaxis.label.set_color("black")
        # set x limit to 0
        ax.set_xlim(xmin=0)
        # set y label
        ax.set_ylabel("shifted T* [°C]", color = "black")
        # set grid
        ax.grid(visible=True, which='major', color='lightgrey', linewidth = 0.5)
        ax.grid(visible=True, which='minor', color='lightgrey', linestyle='dotted', linewidth = 0.5)
        # set aspect ratio of subplot
        ax.set_box_aspect(1)
        ax.set_title(f"Grand Composite Curve of \"{self.label}\"",color = "black", fontsize = 10)
        
        self.gcc_fig, self.gcc_ax = fig, ax

        if save_fig:
            fig.savefig(f"Grand_Composite_Curve_{self.label}.svg")
        if show_fig:
            fig.show()
        if return_fig:
            return fig


    # add: include heat cascades later (for more than one point in time / investigated interval)


    # adding components from tespy models

    # adding a heat pump in the GCC (by referencing evaporator and condenser)
    
    def show_heat_pump_in_gcc(self, evaporator, condenser):
        from tespy.components import SimpleHeatExchanger
        
        # get the GCC
        fig = self.gcc_fig
        ax = self.gcc_ax

        # get the plotting data of the heat exchangers
        if isinstance(condenser,SimpleHeatExchanger):
            # get plot data of heat exchangers at heat sink (taken from user meeting example of heat exchangers)
            condenser_Q_vals = [0, abs(condenser.Q.val)]
            condenser_T_vals = [condenser.outl[0].T.val,condenser.inl[0].T.val]
        if isinstance(condenser,SimpleHeatExchanger):
            # get plot data of heat exchangers at heat source (taken from user meeting example of heat exchangers)
            evaporator_Q_vals = [0, abs(evaporator.Q.val)]         
            evaporator_T_vals = [evaporator.inl[0].T.val,evaporator.outl[0].T.val]

        # add: expand these conditions in future for other types


        # show (only display) by adding the plot data of the heat exchangers to the GCC
        ax.plot(condenser_Q_vals, condenser_T_vals, "-",color="red") # as in heat exchanger example
        ax.plot(evaporator_Q_vals, evaporator_T_vals, "-",color="blue")

        # save figure
        fig.savefig("GCC_with_heat_pump.svg")


    # add: check heat pump integration by using the integration rules
    # see e.g. Arpagaus 2019, Walden et al. 2023

    # add: def add_heat_exchanger_to_streams(self):
        # function to read streams from given heat exchangers, later use to iterate all heat exchangers of network with possibility to exclude streams / heat exchangers