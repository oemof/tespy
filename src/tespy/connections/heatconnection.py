from tespy.connections.connection import connection_registry
from tespy.connections.powerconnection import PowerConnection
from tespy.tools.data_containers import ConnectionParameter
from tespy.tools.data_containers import ReferenceEquation


@connection_registry
class HeatConnection(PowerConnection):

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
