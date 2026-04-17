import numpy as np
import pandas as pd
from fluprodia import FluidPropertyDiagram
from matplotlib import pyplot as plt
from scipy.spatial.distance import cdist

from tespy.networks import Network
from tespy.tools import OptimizationProblem
from tespy.tools import get_plotting_data
from tespy.tools.fluid_properties import single_fluid


class ModelTemplate():

    def __init__(self) -> None:
        self._diagram_cache = {}
        self.parameter_lookup = self._parameter_lookup()
        self._create_network()

    def _create_network(self) -> None:
        self.nw = Network()

    def _parameter_lookup(self) -> dict:
        """
        Define the mapping between flat parameter names and their location in
        the model. Three entry forms are supported:

        - :code:`["Connections", "label", "attr"]` or
          :code:`["Components", "label", "attr"]`: maps to a network object
          attribute. The parameter is both settable and readable.

        - :code:`{"get": callable}`: read-only. The callable takes no arguments
          and returns the current value.

        - :code:`{"set": callable}`: write-only. The callable accepts a single
          value argument. Use this for parameters that cannot be associated
          with a single Connection or Component attribute, e.g. a
          UserDefinedEquation target value.

        - :code:`{"get": callable, "set": callable}`: custom readable and
          settable, e.g. a UDE target that can also be read back.

        Returns
        -------
        dict
            e.g.

            {
                "evaporator pinch": ["Components", "evaporator", "td_pinch"],
                "cop": {"get": self.calc_cop},
                "ude target": {"get": self.get_ude_target, "set": self.set_ude_target},
            }
        """
        return {}

    def get_results(self, labels):
        return {label: self.get_parameter(label) for label in labels}

    def _subcycle_mapping(self) -> dict:
        """Method to extract subcycles based on a label, which maps to
        internal connection labels for plotting cycle diagrams

        Returns
        -------
        dict
            mapping of labels to connection labels in the model
        """
        return {}

    def _map_parameter(self, parameter: str) -> list:
        mapped = self.parameter_lookup[parameter]
        if isinstance(mapped, dict):
            raise TypeError(
                f"'{parameter}' uses a custom get/set spec and cannot be "
                "translated to a Connections/Components path."
            )
        return mapped

    def get_parameter(self, parameter: str) -> float:
        if parameter not in self.parameter_lookup:
            raise KeyError(f"'{parameter}' is not in parameter_lookup.")
        mapped = self.parameter_lookup[parameter]
        if isinstance(mapped, dict):
            if "get" not in mapped:
                raise AttributeError(f"'{parameter}' is write-only.")
            return mapped["get"]()
        if mapped[0] == "Connections":
            return self.nw.get_conn(mapped[1]).get_attr(mapped[2]).val
        elif mapped[0] == "Components":
            return self.nw.get_comp(mapped[1]).get_attr(mapped[2]).val

    def set_parameters(self, **kwargs) -> None:
        conn_params = {}
        comp_params = {}
        for param, value in kwargs.items():
            if param not in self.parameter_lookup:
                raise KeyError(
                    f"The parameter '{param}' is not in parameter_lookup. "
                    f"Available: {', '.join(self.parameter_lookup)}."
                )
            mapped = self.parameter_lookup[param]
            if isinstance(mapped, dict):
                if "set" not in mapped:
                    raise AttributeError(f"'{param}' is read-only.")
                mapped["set"](value)
            else:
                obj_type, label, attr = mapped
                if obj_type == "Connections":
                    conn_params.setdefault(label, {})[attr] = value
                elif obj_type == "Components":
                    comp_params.setdefault(label, {})[attr] = value

        for label, params in conn_params.items():
            self.nw.get_conn(label).set_attr(**params)
        for label, params in comp_params.items():
            self.nw.get_comp(label).set_attr(**params)

    def solve_model_design(self, **kwargs) -> None:
        self.set_parameters(**kwargs)

        self._solved = False
        self.nw.solve("design")

        if self.nw.status == 0:
            self._solved = True
            # path =
            # self.nw.export("stable_solution")
            # self._stable_solution = path
            # save the stable solution for later use in case of model corruption
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
        if fluid_name in self._diagram_cache:
            return self._diagram_cache[fluid_name]

        else:
            diagram = FluidPropertyDiagram(fluid_name)
            diagram.set_unit_system(self.nw.units)
            diagram.set_isolines_subcritical(-20, 200)
            diagram.calc_isolines()
            self._diagram_cache[fluid_name] = diagram

            return diagram


    def _prepare_diagram_and_process_data(self, subcycle, fig, ax, figsize):
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

        fig, ax, processes, points, diagram = self._prepare_diagram_and_process_data(
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

        fig, ax, processes, points, diagram = self._prepare_diagram_and_process_data(
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

    def solve_model_offdesign(self, **kwargs) -> None:
        self.set_parameters(**kwargs)

        self._solved = False
        # Check whether the design path is available
        self.nw.solve("offdesign", design_path=self._design_path)

        if self.nw.status == 0:
            self._solved = True
        elif self.nw.status in [1, 2, 3, 99]:
            # in this case model is very likely corrupted!!
            # fix it by running a presolve using the stable solution
            self._solved = False
            # check whether the design path and stable solution path are available
            self.nw.solve("design", init_only=True, design_path=self._design_path, init_path=self._stable_solution)

    def plot_Ts_diagram_plotly(self, subcycle=None):
        pass

    def plot_logph_diagram_plotly(self, subcycle=None):
        pass

    def plot_QT_diagram_matplotlib(self, heatexchanger_label=None, save_path=None):
        heatex = self.nw.get_comp(heatexchanger_label)
        heat, T_hot, T_cold, _, _ = heatex.calc_sections()

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

    def sensitivity_analysis(self, param_dict=None, result_param_list=None) -> pd.DataFrame:
        """
        1. Check the parameter lengths
        2. Use the order_min_change method
        3. Solve design or offdesign
        4. Deal with the large step changes
        5. Save the results - What results are needed?
            - Function needs to be passed - check for this and raise exception
            - Use the method for the objective function

        Parameters
        ----------
        param_dict : dict
            A dictionary of parameter names and lists of values to be used in the
            sensitivity analysis. All lists must have the same length, which
            determines the number of simulations to be run.

        result_func : function -> dict
            This function will be called after each simulation step and should
            return a dictionary. Its contents will be appended to a pandas
            DataFrame, which is returned after the sensitivity analyses
            finishes. Therefore, the function should return a dictionary of
            string column name and (numeric) result value pairs.

        Returns
        -------
        pandas.core.frame.DataFrame
            DataFrame with input and specified output values of the model
        """
        if param_dict is None:
            raise ValueError(
                "Parameters need to be provided for the sensitivity analysis."
            )
        # never enter kwarg=[] in function signature!
        if result_param_list is None:
            result_param_list = []

        self._check_parameter_lengths(param_dict)

        result_rows = []
        keys = list(param_dict.keys())
        input_values = np.array(list(param_dict.values())).T

        # check current state of model and get nearest to current state
        current_values = np.array([self.get_parameter(k) for k in keys])
        start_idx = int(np.argmin(cdist([current_values], input_values)[0]))

        # sort following simulations and force starting with the nearest to
        # current state
        order = self._order_min_change(input_values, start_idx=start_idx)
        sorted_input = input_values[order]

        # Sensitivity analysis loop:
        for row_num, row in enumerate(sorted_input):
            input_dict = {keys[i]: row[i] for i in range(len(keys))}
            self.solve_model_design(**input_dict)

            if self.nw.status > 1:
                # TODO: get an example, which actually would trigger this
                # and then implement a respective logic
                # previous_input = {
                #     keys[i]: sorted_input[row_num - 1][i]
                #     for i in range(len(keys))
                # }
                # self._intermediate_simulations(previous_input, input_dict)

                # if self.nw.status > 1:
                result_rows.append({
                    **input_dict,
                    **{param: None for param in result_param_list}
                })
                continue

            result_rows.append({
                **input_dict,
                **self.get_results(result_param_list)
            })

        results = pd.DataFrame(result_rows)
        results["_idx"] = order
        return results.sort_values(by="_idx").drop(columns="_idx").reset_index(drop=True)

    # Method for checking the parameter lenghts
    def _check_parameter_lengths(self, param_dict=None):
        lengths = [len(v) for v in param_dict.values()]
        if len(set(lengths)) != 1:
            raise ValueError(
                "All parameters in the sensitivity dictionary must have "
                "the same number of values."
            )

    # Method for ordering function -
    @staticmethod
    def _order_min_change(points: np.ndarray, start_idx: int = 0) -> np.ndarray:
        """Greedy heuristic: always go to the nearest unvisited point."""
        n = len(points)
        dist = cdist(points, points)
        order = [start_idx]
        visited = set(order)
        while len(order) < n:
            last = order[-1]
            # pick nearest unvisited
            candidates = [(i, dist[last, i]) for i in range(n) if i not in visited]
            next_idx = min(candidates, key=lambda x: x[1])[0]
            order.append(next_idx)
            visited.add(next_idx)
        return order

    # Method for handling large step changes
    def _intermediate_simulations(self, start: dict, end: dict) -> None:

        """Run simulation at intermediate design points for more stability."""
        n=0
        num_steps= 2
        max_n=10
        while n< max_n:
            intermediate_points = {
                param: np.linspace(start[param], end[param], num_steps).tolist()
                for param in range(len(start))
                }
            for step in range(num_steps-1):
                try:
                    self.solve_model_design(
                        {key: value[step] for key, value in intermediate_points.items()}
                        )

                except Exception as e:
                    print(f"Error solving model at intermediate point {step}: {e}")
                    for param, value in enumerate(end):
                        start[param]= intermediate_points[param][step]
                        num_steps*=2
                        n+=1
                    continue
                else:
                    continue

        else:
            raise Exception("Simulation failed. Max iteration reached for intermediate simulation step increase.")
