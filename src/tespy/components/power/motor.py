from tespy.components.component import Component
from tespy.components.component import component_registry
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp


@component_registry
class Motor(Component):

    def powerinlets(self):
        return ["power_in"]

    def poweroutlets(self):
        return ["power_out"]

    def get_parameters(self):
        return {
            "eta": dc_cp(**{
                "structure_matrix": self.eta_structure_matrix,
                "func": self.eta_func,
                "dependents": self.eta_dependents,
                "num_eq_sets": 1
            }),
            "eta_char": dc_cc(**{
                "func": self.eta_char_func,
                "dependents": self.eta_char_dependents,
                "num_eq_sets": 1
            })
        }

    def eta_func(self):
        return (
            self.power_inl[0].e.val_SI * self.eta.val
            - self.power_outl[0].e.val_SI
        )

    def eta_structure_matrix(self, k):
        self._structure_matrix[k, self.power_inl[0].e.sm_col] = self.eta.val
        self._structure_matrix[k, self.power_outl[0].e.sm_col] = -1

    def eta_dependents(self):
        return [self.power_inl[0].e, self.power_outl[0].e]

    def eta_char_func(self):
        expr = self.power_inl[0].e.val_SI / self.power_inl[0].e.design
        f = self.eta_char.char_func.evaluate(expr)
        return (
            self.power_inl[0].e.val_SI * self.eta.design * f
            - self.power_outl[0].e.val_SI
        )

    def eta_char_dependents(self):
        return [self.power_inl[0].e, self.power_outl[0].e]

    def calc_parameters(self):
        self.eta.val = self.power_outl[0].e.val_SI / self.power_inl[0].e.val_SI
