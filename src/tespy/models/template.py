from tespy.networks import Network
from tespy.tools.helpers import merge_dicts


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

    def sensitivity_analysis(self, **kwargs) -> None:
        """
        1. Check the parameter lenghts
        2. Use the order_min_change method
        3. Deal with the large step changes
        4. Save the results - What results are needed?
            - Function needs to be passed - check for this and raise exception
            - Use the method for the objective function
        """
        pass
    # Method for objective function - 

    def order_min_change(points: np.ndarray) -> np.ndarray:
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
