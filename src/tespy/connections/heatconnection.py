from tespy.connections.powerconnection import PowerConnection
from tespy.connections.connection import connection_registry
from tespy.tools.data_containers import FluidProperties as dc_prop


@connection_registry
class HeatConnection(PowerConnection):

    def _source_outlets(self, source):
        return source.heatoutlets()

    def _target_inlets(self, target):
        return target.heatinlets()

    def get_parameters(self):
        return {"E": dc_prop(d=1e-4, quantity="heat")}
