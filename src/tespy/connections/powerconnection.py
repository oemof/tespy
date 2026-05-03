from tespy.connections.connection import ConnectionBase
from tespy.connections.connection import connection_registry
from tespy.tools.data_containers import DataContainer as dc
from tespy.tools.data_containers import FluidProperties as dc_prop
from tespy.tools.helpers import TESPyConnectionError
from tespy.tools.logger import logger
from tespy.tools.units import SI_UNITS


@connection_registry
class PowerConnection(ConnectionBase):

    def _source_outlets(self, source):
        return source.poweroutlets()

    def _target_inlets(self, target):
        return target.powerinlets()

    def __init__(self, source, outlet_id, target, inlet_id, label=None, **kwargs):
        self._check_types(source, target)
        self._check_self_connect(source, target)
        self._check_connector_id(source, outlet_id, self._source_outlets(source))
        self._check_connector_id(target, inlet_id, self._target_inlets(target))
        self._init_common(source, outlet_id, target, inlet_id, label, **kwargs)

    def set_attr(self, **kwargs):
        r"""
        Set, reset or unset attributes of a connection.

        Parameters
        ----------
        E : float
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
        for key, value in kwargs.items():
            if key == 'label':
                msg = 'Label can only be specified on instance creation.'
                logger.error(msg)
                raise TESPyConnectionError(msg)
            elif key in self.property_data or key in self.property_data0:
                self._parameter_specification(key, value)
            elif key in ('design', 'offdesign'):
                self._set_design_list(key, value)
            elif key == 'design_path':
                self._set_path_attr(value)
            elif key in ('printout', 'local_design', 'local_offdesign'):
                self._set_bool_attr(key, value)
            else:
                msg = f"Connection has no attribute {key}."
                logger.error(msg)
                raise KeyError(msg)

    def _guess_starting_values(self, units):
        if self.E.is_var and not self.good_starting_values:
            self.E.set_reference_val_SI(0.0)

    def get_variables(self):
        return {"E": self.E}

    def get_parameters(self):
        return {"E": dc_prop(d=1e-4, quantity="power")}

    def calc_results(self, units, skip_postprocess):
        self.E.set_val_from_SI(units)
        self.E.set_val0_from_SI(units)
        return True

    def _get_design_state_SI(self, data, units):
        state = {}
        for var in self._result_attributes():
            param = self.get_attr(var)
            state[var] = units.ureg.Quantity(
                float(data[var]), data[f"{var}_unit"]
            ).to(SI_UNITS[param.quantity]).magnitude
        return state

    def _set_design_params(self, data, units):
        for var, val in self._get_design_state_SI(data, units).items():
            self.get_attr(var).design = val

    def _set_starting_values(self, data, units):
        for prop in self.get_variables():
            var = self.get_attr(prop)
            var.val0 = units.ureg.Quantity(
                float(data[prop]),
                data[f"{prop}_unit"]
            )

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
