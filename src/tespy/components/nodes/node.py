from tespy.components.nodes.merge import Merge
from tespy.components.nodes.splitter import Splitter
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import SimpleDataContainer as dc_simple


class Node(Splitter, Merge):
    r"""
    Class for combined merge and splitting points with multiple inflows and
    outflows.

    **Mandatory Equations**

    - :py:meth:`tespy.components.nodes.base.NodeBase.mass_flow_func`
    - :py:meth:`tespy.components.nodes.base.NodeBase.pressure_structure_matrix`
    - :py:meth:`tespy.components.nodes.node.Node.enthalpy_structure_matrix`
    - :py:meth:`tespy.components.nodes.node.Node.fluid_structure_matrix`
    - :py:meth:`tespy.components.nodes.merge.Merge.fluid_func`
    - :py:meth:`tespy.components.nodes.merge.Merge.energy_balance_func`
    - :py:meth:`tespy.components.nodes.splitter.Splitter.fluid_func`
    - :py:meth:`tespy.components.nodes.splitter.Splitter.energy_balance_func`

    Inlets/Outlets

    - specify number of inlets with :code:`num_in` (default value: 2)
    - specify number of outlets with :code:`num_in` (default value: 2)

    Image

    .. image:: /api/_images/Node.svg
       :alt: flowsheet of the node
       :align: center
       :class: only-light

    .. image:: /api/_images/Node_darkmode.svg
       :alt: flowsheet of the node
       :align: center
       :class: only-dark

    Parameters
    ----------
    label : str
        The label of the component.

    design : list
        List containing design parameters (stated as String).

    offdesign : list
        List containing offdesign parameters (stated as String).

    design_path : str
        Path to the components design case.

    local_offdesign : boolean
        Treat this component in offdesign mode in a design calculation.

    local_design : boolean
        Treat this component in design mode in an offdesign calculation.

    char_warnings : boolean
        Ignore warnings on default characteristics usage for this component.

    printout : boolean
        Include this component in the network's results printout.

    num_in : float
        Number of inlets for this component, default value: 2.

    num_out : float
        Number of outlets for this component, default value: 2.

    Example
    -------
    """

    @staticmethod
    def get_parameters():
        return {
            'num_out': dc_simple(),
            'num_in': dc_simple()
        }

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
