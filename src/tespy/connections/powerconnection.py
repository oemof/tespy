from tespy.connections.connection import ConnectionBase
from tespy.connections.connection import connection_registry
from tespy.tools.data_containers import FluidProperties as dc_prop
from tespy.tools.helpers import TESPyConnectionError
from tespy.tools.logger import logger


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
        if self.E.is_var:
            if not self.good_starting_values:
                self.E.set_reference_val_SI(0.0)
            else:
                # Push the value preserved across the previous solve's
                # `detach()` into the freshly (re-)built reference
                # container, mirroring what `Connection._guess_starting_values`
                # does for m/p/h. Without this, presolve's fresh
                # `_reference_container` for `E` starts uninitialized on
                # every solve call after the first, discarding the warm
                # start even when good_starting_values is True.
                self.E.set_reference_val_SI(self.E._val_SI)

    def get_variables(self):
        return {"E": self.E}

    def get_parameters(self):
        return {"E": dc_prop(d=1e-4, quantity="power")}

    def calc_results(self, units, skip_postprocess):
        self.E.set_val_from_SI(units)
        self.E.set_val0_from_SI(units)
        return True

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
