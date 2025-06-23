from tespy.components.component import Component
from tespy.components.component import component_registry
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import SimpleDataContainer as dc_simple


@component_registry
class PowerBus(Component):

    def powerinlets(self):
        return [f"power_in{i + 1}" for i in range(self.num_inlets.val)]

    def poweroutlets(self):
        return [f"power_out{i + 1}" for i in range(self.num_outlets.val)]

    def get_parameters(self):
        return {
            "num_inlets": dc_simple(),
            "num_outlets": dc_simple()
        }

    def get_mandatory_constraints(self):
        return {
            "energy_balance_constraint": dc_cmc(**{
                "func": self.energy_balance_func,
                "dependents": self.energy_balance_dependents,
                "num_eq_sets": 1
            })
        }

    def energy_balance_func(self):
        residual = 0
        for i in self.power_inl:
            residual += i.e.val_SI
        for o in self.power_outl:
            residual -= o.e.val_SI
        return residual

    def energy_balance_dependents(self):
        return [c.e for c in self.power_inl + self.power_outl]