from tespy.connections.powerconnection import PowerConnection
from tespy.connections.connection import connection_registry


@connection_registry
class HeatConnection(PowerConnection):

    def _source_outlets(self, source):
        return source.heatoutlets()

    def _target_inlets(self, target):
        return target.heatinlets()
