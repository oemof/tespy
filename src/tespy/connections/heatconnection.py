from tespy.connections.connection import connection_registry
from tespy.connections.powerconnection import PowerConnection
from tespy.tools.data_containers import ConnectionParameter
from tespy.tools.data_containers import ReferenceEquation


@connection_registry
class HeatConnection(PowerConnection):
    r"""
    Connection for heat flows, e.g. between a heat exchanger and a heat
    source or sink component.

    The connection holds the heat flow :code:`E` as its only variable and
    links the heat outlet of the source to the heat inlet of the target
    component.

    Parameters
    ----------

    design : list
        List containing design parameters (stated as String).

    design_path : str
        Path to the individual design case for this connection.

    E : float, Ref
        Heat flow (system variable).

    label : str
        The label of the connection.

    local_design : bool
        Treat this connection in design mode in an offdesign calculation.

    local_offdesign : bool
        Treat this connection in offdesign mode in a design calculation.

    offdesign : list
        List containing offdesign parameters (stated as String).

    printout : bool
        Include this connection in the network's results printout.
    """

    def _source_outlets(self, source):
        return source.heatoutlets()

    def _target_inlets(self, target):
        return target.heatinlets()

    def get_parameters(self):
        return {
            "E": ConnectionParameter(
                d=1e-4, quantity="heat",
                description="heat flow (system variable)",
                reference=ReferenceEquation(
                    num_eq=1,
                    func_params={"variable": "E"},
                    structure_matrix=self.primary_ref_structure_matrix,
                    quantity="heat",
                    description="equation for linear relationship between two heat flows"
                )
            )
        }
