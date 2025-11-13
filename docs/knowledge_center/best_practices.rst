.. _best_practices_label:

Best Practices
~~~~~~~~~~~~~~

High-level API for tespy models
###############################

Requirements for models in workflows

- models will be executed more than once
- specific results/processings/... will be retrieved frequently
- model crashes should be avoided or at least be handled in some way without
  leaving side effects (e.g. corrupting future models)
- there should be a way to exchange data (input and output) with the model
  in a structured way
- there might be the need for various tespy models doing similar things
  - various topologies
  - various working fluids
  - ...
- model inputs might change a lot (keep in mind: EVERY simulation has side effects!)

(A possible) solution

disclaimer: This is not universal truth, this is what is working well for me.
Feedback/suggestions for improvements greatly appreciated!

- Structure
  - each model is wrapped in a class (make use of inheritance!)
  - the class has methods
    - to set up the Network and create a initial stable solution
    - to input component or connection parameters
    - to solve the model in design mode
    - to retrieve component or connection parameters
  - optionally, there are methods
    - to automatically generate cycle diagrams, QT diagrams
    - to solve the model in offdesign mode

- Parameter input and output
  - via dictionaries/json/yaml
  - keys are mapped to specific parameters in the model, e.g.
    "evaporator_pinch" mapped to the parameter **td_pinch** of the
    **component** **evaporator**
  - there is a second structure in between: nested dictionaries to access
    components and connections
    - user specifies {"evaporator_pinch": 10}
    - parameter lookup -> "evaporator_pinch": ["Components", "evaporator", "td_pinch"]
    - nested dict -> {"Components": {"evaporator": {"td_pinch": 10}}}
    - internally this dict is used to set values for model or retrieve results

- Solving methods
  - solve the model
  - check the status variable of the model
  - handle issues if status is not 0 or 1
  - handle what happens with converged simulation violating physical limits (status 1)
  - e.g. set a flat if the model has been solved successfully

Template

.. code-block:: python

        >>> from tespy.tools.helpers import merge_dicts
        >>> from tespy.networks import Network


        >>> class ModelTemplate():
        ...
        ...    def __init__(self, config) -> None:
        ...        self.config = config
        ...        self.parameter_lookup = self._parameter_lookup()
        ...        self._create_network()
        ...
        ...    def _create_network(self) -> None:
        ...        self.nw = Network()
        ...        self.nw.units.set_defaults(
        ...            **self.config["units"]
        ...        )
        ...
        ...    def _parameter_lookup(self) -> dict:
        ...        return {
        ...            "evaporator_pinch": ["Components", "evaporator", "td_pinch"],
        ...            "T_geo": ["Connections", "a1", "T"],
        ...            "m_geo": ["Connections", "a1", "m"]
        ...        }
        ...
        ...    def _map_parameter(self, parameter: str) -> tuple:
        ...        return self.parameter_lookup[parameter]
        ...
        ...    def _map_to_input_dict(self, **kwargs) -> dict:
        ...        input_dict = {}
        ...        for param, value in kwargs.items():
        ...            if param not in self.parameter_lookup:
        ...                msg = (
        ...                    f"The parameter {param} is not mapped to any input of the "
        ...                    "model. The following parameters are available:\n"
        ...                    f"{', '.join(self.parameter_lookup)}."
        ...                )
        ...                raise KeyError(msg)
        ...            key = self._map_parameter(param)
        ...            input_dict = merge_dicts(
        ...                input_dict,
        ...                {key[0]: {key[1]: {key[2]: value}}}
        ...            )
        ...        return input_dict
        ...
        ...    def get_parameter(self, parameter: str) -> float:
        ...        mapped = self._map_parameter(parameter)
        ...        if mapped[0] == "Connections":
        ...            return self.nw.get_conn(mapped[1]).get_attr(mapped[2]).val
        ...
        ...        elif mapped[0] == "Components":
        ...            return self.nw.get_comp(mapped[1]).get_attr(mapped[2]).val
        ...
        ...    def set_parameters(self, **kwargs) -> None:
        ...        input_dict = self._map_to_input_dict(**kwargs)
        ...        if "Connections" in input_dict:
        ...            for c, params in input_dict["Connections"].items():
        ...                self.nw.get_conn(c).set_attr(**params)
        ...
        ...        if "Components" in input_dict:
        ...            for c, params in input_dict["Components"].items():
        ...                self.nw.get_comp(c).set_attr(**params)
        ...
        ...    def solve_model(self, **kwargs) -> None:
        ...        self.set_parameters(**kwargs)
        ...
        ...        self._solved = False
        ...        self.nw.solve("design")
        ...
        ...        if self.nw.status == 0:
        ...            self._solved = True
        ...        # is not required in this example, but could lead to handling some
        ...        # stuff
        ...        elif self.nw.status == 1:
        ...            self._solved = False
        ...        elif self.nw.status in [2, 3, 99]:
        ...            # in this case model is very likely corrupted!!
        ...            # fix it by running a presolve using the stable solution
        ...            self._solved = False
        ...            self.nw.solve("design", init_only=True, init_path=self._stable_solution)
        ...
        ...    def solve_model_offdesign(self, design, **kwargs) -> None:
        ...        self.set_parameters(**kwargs)
        ...
        ...        self._solved = False
        ...        self.nw.solve("offdesign", design_path=design)
        ...
        ...        if self.nw.status == 0:
        ...            self._solved = True
        ...        # is not required in this example, but could lead to handling some
        ...        # stuff
        ...        elif self.nw.status == 1:
        ...            self._solved = False
        ...        elif self.nw.status in [2, 3, 99]:
        ...            # in this case model is very likely corrupted!!
        ...            # fix it by running a presolve using the stable solution
        ...            self._solved = False
        ...            self.nw.solve("offdesign", init_only=True, init_path=self._stable_solution, design_path=design)
