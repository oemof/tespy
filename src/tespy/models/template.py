import numpy as np
from fluprodia import FluidPropertyDiagram
from matplotlib import pyplot as plt

from tespy.networks import Network
from tespy.tools import get_plotting_data
from tespy.tools.fluid_properties import single_fluid
from tespy.tools.helpers import merge_dicts


class ModelTemplate():

    _DIAGRAM_CACHE = {}

    def __init__(self) -> None:
        self.parameter_lookup = self._parameter_lookup()
        self._create_network()

    def _create_network(self) -> None:
        self.nw = Network()

    def _parameter_lookup(self) -> dict:
        """
        Example


        Returns
        -------
        dict
            return a mapping between single labels and their wiring to the
            internal model, e.g.

            {
                "evaporator pinch": ["Component", "evaporator", "td_pinch"],
                "evaporation temperature": ["Connection", "b1", "T"]
            }
        """
        return {}

    def _subcycle_mapping(self) -> dict:
        """Method to extract subcycles based on a label, which maps to
        internal connection labels for plotting cycle diagrams

        Returns
        -------
        dict
            mapping of labels to connection labels in the model
        """
        return {}

    def _map_parameter(self, parameter: str) -> tuple:
        return self.parameter_lookup[parameter]

    def _map_to_input_dict(self, **kwargs) -> dict:
        input_dict = {}
        for param, value in kwargs.items():
            if param not in self.parameter_lookup:
                msg = (
                    f"The parameter {param} is not mapped to any input of the "
                    "model. The following parameters are available:\n"
                    f"{', '.join(self.parameter_lookup)}."
                )
                raise KeyError(msg)
            key = self._map_parameter(param)
            input_dict = merge_dicts(
                input_dict,
                {key[0]: {key[1]: {key[2]: value}}}
            )
        return input_dict

    def get_parameter(self, parameter: str) -> float:
        mapped = self._map_parameter(parameter)
        if mapped[0] == "Connections":
            return self.nw.get_conn(mapped[1]).get_attr(mapped[2]).val

        elif mapped[0] == "Components":
            return self.nw.get_comp(mapped[1]).get_attr(mapped[2]).val

    def set_parameters(self, **kwargs) -> None:
        input_dict = self._map_to_input_dict(**kwargs)
        if "Connections" in input_dict:
            for c, params in input_dict["Connections"].items():
                self.nw.get_conn(c).set_attr(**params)

        if "Components" in input_dict:
            for c, params in input_dict["Components"].items():
                self.nw.get_comp(c).set_attr(**params)

    def solve_model_design(self, **kwargs) -> None:
        self.set_parameters(**kwargs)

        self._solved = False
        self.nw.solve("design")

        if self.nw.status == 0:
            self._solved = True
        # is not required in this example, but could lead to handling some
        # stuff
        elif self.nw.status == 1:
            self._solved = False
        elif self.nw.status in [2, 3, 99]:
            # in this case model is very likely corrupted!!
            # fix it by running a presolve using the stable solution
            self._solved = False
            self.nw.solve("design", init_only=True, init_path=self._stable_solution)

    def _get_diagram(self, fluid_name):
        if fluid_name in self._DIAGRAM_CACHE:
            return self._DIAGRAM_CACHE[fluid_name]

        else:
            diagram = FluidPropertyDiagram(fluid_name)
            diagram.set_unit_system(self.nw.units)
            diagram.set_isolines_subcritical(-20, 200)
            diagram.calc_isolines()
            self._DIAGRAM_CACHE[fluid_name] = diagram

            return diagram


    def _prepare_diagram_and_proccess_data(self, subcycle, fig, ax, figsize):
        connection_label = self._subcycle_mapping().get(subcycle)

        if connection_label is None:
            raise ValueError("subcycle is unknown")

        fluid_name = single_fluid(self.nw.get_conn(connection_label).fluid_data)

        if fluid_name is None:
            raise ValueError("Fluid is mixture")

        diagram = self._get_diagram(fluid_name)

        processes, points = get_plotting_data(self.nw, connection_label)
        processes = {
            key: diagram.calc_individual_isoline(**value)
            for key, value in processes.items()
            if value is not None
        }

        if fig is None or ax is None:
            if figsize is None:
                figsize = (10, 6)
            fig, ax = plt.subplots(figsize=figsize)

        return fig, ax, processes, points, diagram

    def _plot_processes_and_states(self, ax, processes, points, x_property, y_property):

        for label, values in processes.items():
            _ = ax.plot(values[x_property], values[y_property], label=label, color="tab:red", zorder=10000)

        for label, point in points.items():
            _ = ax.scatter(point[x_property], point[y_property], label=label, color="tab:red", zorder=10000)


    def plot_Ts_diagram_matplotlib(self, subcycle=None, save_path=None, fig=None, ax=None, x_min=None, x_max=None, y_min=None, y_max=None, figsize=None):

        fig, ax, processes, points, diagram = self._prepare_diagram_and_proccess_data(
            subcycle, fig, ax, figsize
        )

        if x_min is None or x_max is None:
            x_min, x_max = self._make_cycle_plot_limits(points, "s", "lin")
        if y_min is None or y_max is None:
            y_min, y_max = self._make_cycle_plot_limits(points, "T", "lin")

        diagram.draw_isolines(
            fig, ax, "Ts", x_min, x_max, y_min, y_max,
            isoline_data={
                "Q": {"values": np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])},
            }
        )
        self._plot_processes_and_states(ax, processes, points, "s", "T")

        if save_path:
            fig.savefig(f"{save_path}/ts_diagram.svg", bbox_inches="tight")

        return fig, ax

    def plot_logph_diagram_matplotlib(self, subcycle=None, save_path=None, fig=None, ax=None, x_min=None, x_max=None, y_min=None, y_max=None, figsize=None):

        fig, ax, processes, points, diagram = self._prepare_diagram_and_proccess_data(
            subcycle, fig, ax, figsize
        )

        if x_min is None or x_max is None:
            x_min, x_max = self._make_cycle_plot_limits(points, "h", "lin")
        if y_min is None or y_max is None:
            y_min, y_max = self._make_cycle_plot_limits(points, "p", "log")

        diagram.draw_isolines(
            fig, ax, "logph", x_min, x_max, y_min, y_max,
            isoline_data={
                "s": {"values": np.array([])},
                "Q": {"values": np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])},
            }
        )
        self._plot_processes_and_states(ax, processes, points, "h", "p")

        if save_path:
            fig.savefig(f"{save_path}/log_ph_diagram.svg", bbox_inches="tight")
        return fig, ax


    def plot_Ts_diagram_plotly(self, subcycle=None):
        pass

    def plot_logph_diagram_plotly(self, subcycle=None):
        pass

    def plot_QT_diagram_matplotlib(self, heatexchanger_label=None, save_path=None):
        heatex = self.nw.get_comp(heatexchanger_label)
        heat, T_hot, T_cold, heat_per_section, td_log_per_section = heatex.calc_sections()
        fig, ax = plt.subplots(1)

        ax.plot(heat, T_hot, "o-", color="red")
        ax.plot(heat, T_cold, "o-", color="blue")
        if save_path:
            fig.savefig(f"{save_path}/qt_diagram.svg", bbox_inches="tight")
        return fig, ax

    def plot_QT_diagram_plotly(self, heatexchanger_label=None):
        pass

    def _make_cycle_plot_limits(self, states: list, quantity: str, scale: str, padding_rel=0.1) -> tuple:
        """Automatically retrieve the limits for an axes based on the process
        point limits in one axis

        Parameters
        ----------
        states : list
            List of process states
        quantity : str
            Name of the quantity, e.g. :code:`T`, :code:`h`
        scale : str
            Scale of the axis to plot on
        padding_rel : float, optional
            relative difference to overall distance between min and max value,
            by default 0.1

        Returns
        -------
        tuple
            minimum and maximum value for axis
        """
        all_values = [point[quantity] for point in states.values()]
        min_val = min(all_values)
        max_val = max(all_values)

        if scale == 'lin':
            delta_val = max_val - min_val
            ax_min_val = min_val - padding_rel * delta_val
            ax_max_val = max_val + padding_rel * delta_val
        elif scale == 'log':
            delta_val = np.log10(max_val) - np.log10(min_val)
            ax_min_val = 10 ** (np.log10(min_val) - padding_rel * delta_val)
            ax_max_val = 10 ** (np.log10(max_val) + padding_rel * delta_val)

        return ax_min_val, ax_max_val
