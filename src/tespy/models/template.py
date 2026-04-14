import pandas as pd

from tespy.networks import Network
from tespy.tools.helpers import merge_dicts
from scipy.spatial.distance import cdist
import numpy as np


class ModelTemplate():

    def __init__(self) -> None:
        self.parameter_lookup = self._parameter_lookup()
        self._create_network()

    def _create_network(self) -> None:
        self.nw = Network()

    def _parameter_lookup(self) -> dict:
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

    def solve_model_offdesign(self, **kwargs) -> None:
        pass

    def sensitivity_analysis(self, param_dict=None, result_func=None) -> None:
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
        """
        if param_dict is None:
            raise ValueError(
                "Parameters need to be provided for the sesitivity analysis."
            )

        if result_func is None:
            raise ValueError(
                "No 'result_func' keyword argument was passed. It is necessary"
                + " to extract results for the sensitivity analysis."
            )

        self._check_parameter_lengths(param_dict)

        result_rows = []

        # Sensitivity analysis loop:
        for i in range(len(list(param_dict.values())[0])):
            try:
                self.solve_model_design(
                        **{key: value[i] for key, value in param_dict.items()}
                        )
            except:

                self._intermediate_simulations( 
                    {key: value[i-1] for key, value in param_dict.items()}, 
                    {key: value[i] for key, value in param_dict.items()})
                try:
                    self.solve_model_design(
                            **{key: value[i] for key, value in param_dict.items()}
                            )
                except Exception as e:
                    raise e
                else:
                    result_rows.append(result_func())
                    continue
            else:
                result_rows.append(result_func())
                continue
        results = pd.DataFrame(result_rows)

        return results

    # Method for checking the parameter lenghts
    def _check_parameter_lengths(self, param_dict=None):
        lengths = [len(v) for v in param_dict.values()]
        if len(set(lengths)) != 1:
            raise ValueError(
                "All parameters in the sensitivity dictionary must have "
                "the same number of values."
            )

    # Method for ordering function - 
    def _order_min_change(points: np.ndarray) -> np.ndarray:
        """Greedy heuristic: always go to the nearest unvisited point."""
        n = len(points)
        dist = cdist(points, points)
        order = [0]
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



    def plot_Ts_diagram_matplotlib(self, subcycle=None):
        if subcycle is not None:
            relevant_connection = self._subcycle_mapping[subcycle]
        else:
            pass
        # if no subcycle is provided all of them are plotted in individual figures

    def plot_logph_diagram_matplotlib(self, subcycle=None):
        pass

    def plot_Ts_diagram_plotly(self, subcycle=None):
        pass

    def plot_logph_diagram_plotly(self, subcycle=None):
        pass

    def plot_QT_diagram_matplotlib(self, heatexchanger_label=None):
        pass

    def plot_QT_diagram_plotly(self, heatexchanger_label=None):
        pass
