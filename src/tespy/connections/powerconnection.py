import numpy as np

from tespy.connections.connection import ConnectionBase
from tespy.connections.connection import connection_registry
from tespy.tools.data_containers import DataContainer as dc
from tespy.tools.data_containers import FluidProperties as dc_prop
from tespy.tools.helpers import TESPyConnectionError
from tespy.tools.logger import logger


@connection_registry
class PowerConnection(ConnectionBase):

    def __init__(self, source, outlet_id, target, inlet_id, label=None, **kwargs):
        self._check_types(source, target)
        self._check_self_connect(source, target)
        self._check_connector_id(source, outlet_id, source.poweroutlets())
        self._check_connector_id(target, inlet_id, target.powerinlets())

        self.label = f"{source.label}:{outlet_id}_{target.label}:{inlet_id}"
        if label is not None:
            self.label = label
            if not isinstance(label, str):
                msg = "Please provide the label as string."
                logger.error(msg)
                raise TypeError(msg)

        # set specified values
        self.source = source
        self.source_id = outlet_id
        self.target = target
        self.target_id = inlet_id

        # defaults
        self.new_design = True
        self.design_path = None
        self.design = []
        self.offdesign = []
        self.local_design = False
        self.local_offdesign = False
        self.printout = True

        # set default values for kwargs
        self.property_data = self.get_parameters()
        self.property_data0 = [x + '0' for x in self.property_data.keys()]
        self.parameters = {
            k: v for k, v in self.get_parameters().items()
            if hasattr(v, "func") and v.func is not None
        }
        self.__dict__.update(self.property_data)
        msg = (
            f"Created connection from {self.source.label} ({self.source_id}) "
            f"to {self.target.label} ({self.target_id})."
        )
        logger.debug(msg)
        self.set_attr(**kwargs)

    def set_attr(self, **kwargs):
        r"""
        Set, reset or unset attributes of a connection.

        Parameters
        ----------
        e : float
            Energy flow specification.

        design : list
            List containing design parameters (stated as string).

        offdesign : list
            List containing offdesign parameters (stated as string).

        design_path : str
            Path to individual design case for this connection.

        local_offdesign : boolean
            Treat this connection in offdesign mode in a design calculation.

        local_design : boolean
            Treat this connection in design mode in an offdesign calculation.

        printout : boolean
            Include this connection in the network's results printout.
        """
        # set specified values
        for key in kwargs:
            if key == 'label':
                msg = 'Label can only be specified on instance creation.'
                logger.error(msg)
                raise TESPyConnectionError(msg)

            elif key in self.property_data or key in self.property_data0:
                self._parameter_specification(key, kwargs[key])

            # design/offdesign parameter list
            elif key in ['design', 'offdesign']:
                if not isinstance(kwargs[key], list):
                    msg = f"Please provide the {key} parameters as list!"
                    logger.error(msg)
                    raise TypeError(msg)
                elif set(kwargs[key]).issubset(self.property_data.keys()):
                    self.__dict__.update({key: kwargs[key]})
                else:
                    params = ', '.join(self.property_data.keys())
                    msg = (
                        "Available parameters for (off-)design specification "
                        f"are: {params}."
                    )
                    logger.error(msg)
                    raise ValueError(msg)

            # design path
            elif key == 'design_path':
                self.__dict__.update({key: kwargs[key]})
                self.new_design = True

            # other boolean keywords
            elif key in ['printout', 'local_design', 'local_offdesign']:
                if not isinstance(kwargs[key], bool):
                    msg = ('Please provide the ' + key + ' as boolean.')
                    logger.error(msg)
                    raise TypeError(msg)
                else:
                    self.__dict__.update({key: kwargs[key]})

            # invalid keyword
            else:
                msg = f"Connection has no attribute {key}."
                logger.error(msg)
                raise KeyError(msg)

    def _precalc_guess_values(self):
        pass

    def _presolve(self):
        return []

    def _reset_design(self, redesign):
        for value in self.get_variables().values():
            value.design = np.nan

        self.new_design = True

        # switch connections to design mode
        if redesign:
            for var in self.design:
                self.get_attr(var).is_set = True

            for var in self.offdesign:
                self.get_attr(var).is_set = False

    def get_variables(self):
        return {"E": self.E}

    def get_parameters(self):
        return {"E": dc_prop(d=1e-4, quantity="power")}

    def calc_results(self):
        self.E.set_val_from_SI()
        self.E.set_val0_from_SI()
        return True

    def _set_design_params(self, data):
        for var in self._result_attributes():
            self.get_attr(var).design = float(data[var])

    def _set_starting_values(self, data):
        for prop in self.get_variables():
            var = self.get_attr(prop)
            var.val0 = float(data[prop])
            var.unit = data[prop + '_unit']

    @classmethod
    def _print_attributes(cls):
        return ["E"]

    @classmethod
    def _result_attributes(cls):
        return ["E"]

    @classmethod
    def _get_result_cols(cls, all_fluids):
        return ["E", "E_unit"]

    def collect_results(self, all_fluids):
        return [self.E.val, self.E.unit]

    def _deserialize(self, data, all_connections):
        arglist = [
            _ for _ in data
            if _ not in ["source", "source_id", "target", "target_id", "label", "fluid"]
            and "ref" not in _
        ]

        for arg in arglist:
            container = self.get_attr(arg)
            if isinstance(container, dc):
                container.set_attr(**data[arg])
            else:
                self.set_attr(**{arg: data[arg]})

    def _to_exerpy(self, pamb, Tamb):
        connection_json = {}

        if self.source.__class__.__name__ in ["Motor", "Generator"]:
            source_connector = 0
        elif self.source.__class__.__name__ in ["Turbine"]:
            source_connector = 1
        elif self.source.__class__.__name__ in ["SimpleHeatExchanger"]:
            source_connector = 1
        elif self.source.__class__.__name__ in ["PowerBus"]:
            if self.source_id.startswith("power_out"):
                s_id = self.source_id.removeprefix("power_out")
                source_connector = 0 if s_id == "" else int(s_id) - 1
            elif self.source_id.startswith("power_in"):
                s_id = self.source_id.removeprefix("power_in")
                source_connector = 0 if s_id == "" else int(s_id) - 1
        else:
            source_connector = 999

        if self.target.__class__.__name__ in ["Motor", "Generator"]:
            target_connector = 0
        elif self.target.__class__.__name__ in ["Compressor", "Pump"]:
            target_connector = 1
        elif self.target.__class__.__name__ in ["SimpleHeatExchanger"]:
            target_connector = 1
        elif self.target.__class__.__name__ in ["PowerBus"]:
            if self.target_id.startswith("power_in"):
                t_id = self.target_id.removeprefix("power_in")
                target_connector = 0 if t_id == "" else int(t_id) - 1
            elif self.target_id.startswith("power_out"):
                t_id = self.target_id.removeprefix("power_out")
                target_connector = 0 if t_id == "" else int(t_id) - 1
        else:
            target_connector = 999

        connection_json[self.label] = {
            "source_component": self.source.label,
            "source_connector": source_connector,
            "target_component": self.target.label,
            "target_connector": target_connector
        }
        if self.source_id == "heat" or self.target_id == "heat":
            kind = "heat"
        else:
            kind = "power"

        connection_json[self.label].update({
            "kind": kind,
            "energy_flow": self.E.val_SI
        })

        return connection_json
