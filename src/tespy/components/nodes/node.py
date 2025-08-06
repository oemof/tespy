from tespy.components import Splitter, Merge
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc


class Node(Splitter, Merge):

    @staticmethod
    def get_parameters():
        return {
            'num_out': dc_simple(),
            'num_in': dc_simple()
        }

    def _update_num_eq(self):
        self.variable_fluids = set(
            [fluid for c in self.inl + self.outl for fluid in c.fluid.is_var]
        )
        set_fluids = set(
            [fluid for c in self.inl + self.outl for fluid in c.fluid.is_set]
        )
        self.all_fluids = self.variable_fluids | set_fluids
        if len(self.variable_fluids) == 0 and len(set_fluids) == 0:
            fluid_eq = 0
            self.constraints["mass_flow_constraints"].num_eq = 1
        elif len(self.variable_fluids) == 0:
            fluid_eq = len(self.all_fluids)
            self.constraints["mass_flow_constraints"].num_eq = 0
        else:
            fluid_eq = len(self.variable_fluids)
        self.constraints["fluid_constraints"].num_eq = fluid_eq

    def get_mandatory_constraints(self):
        return {
            'mass_flow_constraints': dc_cmc(**{
                'func': self.mass_flow_func,
                'dependents': self.mass_flow_dependents,
                'num_eq_sets': 1
            }),
            'pressure_constraints': dc_cmc(**{
                'structure_matrix': self.pressure_structure_matrix,
                'num_eq_sets': self.num_i + self.num_o - 1
            }),
            'outlet_enthalpy_constraints': dc_cmc(**{
                'structure_matrix': self.enthalpy_structure_matrix,
                'num_eq_sets': self.num_o - 1
            }),
            'outlet_fluid_constraints': dc_cmc(**{
                'structure_matrix': self.fluid_structure_matrix,
                'num_eq_sets': self.num_o - 1
            }),
            'fluid_constraints': dc_cmc(**{
                'num_eq_sets': 1,
                'func': self.fluid_func,
                'deriv': self.fluid_deriv,
                'dependents': self.fluid_dependents
            }),
            'energy_balance_constraints': dc_cmc(**{
                'num_eq_sets': 1,
                'func': self.energy_balance_func,
                'dependents': self.energy_balance_dependents,
            })
        }

    def inlets(self):
        if self.num_in.is_set:
            return [f'in{i + 1}' for i in range(self.num_in.val)]
        else:
            self.set_attr(num_in=2)
            return self.inlets()

    def outlets(self):
        if self.num_out.is_set:
            return [f'out{i + 1}' for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=2)
            return self.outlets()

    def enthalpy_structure_matrix(self, k):
        r"""
        Set up the structure matrix for the outlet enthalpy constraints:

        All outlet enthalpy values must be identical.
        """
        for eq, conn in enumerate(self.outl[1:]):
            self._structure_matrix[k + eq, self.outl[0].h.sm_col] = 1
            self._structure_matrix[k + eq, conn.h.sm_col] = -1

    def fluid_structure_matrix(self, k):
        r"""
        Set up the structure matrix for the outlet fluid constraints:

        All outlet fluid compositions must be identical.
        """
        for eq, conn in enumerate(self.outl[1:]):
            self._structure_matrix[k + eq, self.outl[0].fluid.sm_col] = 1
            self._structure_matrix[k + eq, conn.fluid.sm_col] = -1
