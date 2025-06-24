from tespy.components.component import Component
from tespy.components.component import component_registry
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp


@component_registry
class Motor(Component):
    """
    A motor converts electrical energy into mechanical energy.

    **Mandatory Equations**

    - None

    **Optional Equations**

    - :py:meth:`tespy.components.power.motor.Motor.eta_func`
    - :py:meth:`tespy.components.power.motor.Motor.delta_power_func`
    - :py:meth:`tespy.components.power.motor.Motor.eta_char_func`

    Inlets/Outlets

    - None

    Optional inlets/outlets

    - power_in
    - power_out

    Image

    .. image:: /api/_images/Motor.svg
       :alt: flowsheet of the motor
       :align: center
       :class: only-light

    .. image:: /api/_images/Motor_darkmode.svg
       :alt: flowsheet of the motor
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

    eta : float, dict
        Outlet to inlet efficiency, :math:`\eta/1`

    delta_power : float, dict
        Fixed power offset, :math:`\text{delta_power}/\text{W}`

    eta_char : tespy.tools.characteristics.CharLine, dict
        Characteristic line for efficiency to power as function of design
        efficiency.

    Example
    -------
    A compressor provides compressed air which is used in a compressed air
    distribution system. The energy is provided by an electrical motor.

    >>> from tespy.components import Sink, Source, Compressor, Motor, PowerSource
    >>> from tespy.connections import Connection, PowerConnection
    >>> from tespy.networks import Network
    >>> import os
    >>> nw = Network(p_unit='bar', T_unit='C', iterinfo=False)
    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> compressor = Compressor('compressor')

    Ambient air flows into the compressor and is ejected at 4 bar. We can set
    the system up without the use of any of the power components.

    >>> c1 = Connection(so, 'out1', compressor, 'in1')
    >>> c2 = Connection(compressor, 'out1', si, 'in1')
    >>> nw.add_conns(c1, c2)
    >>> c1.set_attr(fluid={'air': 1}, T=25, p=1, m=1)
    >>> c2.set_attr(p=4)
    >>> compressor.set_attr(eta_s=0.8)
    >>> nw.solve('design')

    We can add the Motor and a PowerSource and then connect these parts to
    the compressor.

    >>> motor = Motor('motor')
    >>> power_source = PowerSource('power source')
    >>> e1 = PowerConnection(power_source, 'power', motor, 'power_in')
    >>> e2 = PowerConnection(motor, 'power_out', compressor, 'power')
    >>> nw.add_conns(e1, e2)

    Now we have added two variables to our problem (the power flows of e1 and
    e2), but only one equation (the power balance for the compressor). The
    connection between the two power flows can be made through specifying the
    efficiency of the motor:

    >>> motor.set_attr(eta=.98)
    >>> nw.solve('design')
    >>> nw.assert_convergence()
    >>> round(e2.E.val_SI) == round(compressor.P.val)
    True
    >>> round(e1.E.val_SI) == round(compressor.P.val / 0.98)
    True

    We could also specify the electrical energy instead of fixing the air
    mass flow to calculate the resulting air mass flow:

    >>> e1.set_attr(E=1e5)
    >>> c1.set_attr(m=None)
    >>> nw.solve('design')
    >>> round(c1.m.val, 3)
    0.539

    Or, fix both (electrical and mechanical power flows) and leave open the
    motor efficiency:

    >>> e2.set_attr(E=0.9e5)
    >>> motor.set_attr(eta=None)
    >>> nw.solve('design')
    >>> round(motor.eta.val, 2)
    0.9

    >>> e2.set_attr(E=None)
    >>> motor.set_attr(delta_power=5e3)
    >>> nw.solve('design')
    >>> round(motor.eta.val, 3)
    0.95
    """

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
                "num_eq_sets": 1,
                "max_val": 1,
                "min_val": 0
            }),
            "delta_power": dc_cp(**{
                "structure_matrix": self.delta_power_structure_matrix,
                "func": self.delta_power_func,
                "dependents": self.delta_power_dependents,
                "num_eq_sets": 1,
                "min_val": 0
            }),
            "eta_char": dc_cc(**{
                "func": self.eta_char_func,
                "dependents": self.eta_char_dependents,
                "num_eq_sets": 1
            })
        }

    def eta_func(self):
        return (
            self.power_inl[0].E.val_SI * self.eta.val
            - self.power_outl[0].E.val_SI
        )

    def eta_structure_matrix(self, k):
        self._structure_matrix[k, self.power_inl[0].E.sm_col] = self.eta.val
        self._structure_matrix[k, self.power_outl[0].E.sm_col] = -1

    def eta_dependents(self):
        return [self.power_inl[0].E, self.power_outl[0].E]

    def delta_power_func(self):
        return (
            self.power_inl[0].E.val_SI - self.power_outl[0].E.val_SI
            - self.delta_power.val
        )

    def delta_power_structure_matrix(self, k):
        self._structure_matrix[k, self.power_inl[0].E.sm_col] = 1
        self._structure_matrix[k, self.power_outl[0].E.sm_col] = -1
        self._rhs[k] = self.delta_power.val

    def delta_power_dependents(self):
        return [self.power_inl[0].E, self.power_outl[0].E]

    def eta_char_func(self):
        expr = self.power_inl[0].E.val_SI / self.power_inl[0].E.design
        f = self.eta_char.char_func.evaluate(expr)
        return (
            self.power_inl[0].E.val_SI * self.eta.design * f
            - self.power_outl[0].E.val_SI
        )

    def eta_char_dependents(self):
        return [self.power_inl[0].E, self.power_outl[0].E]

    def calc_parameters(self):
        self.eta.val = self.power_outl[0].E.val_SI / self.power_inl[0].E.val_SI
        self.delta_power.val = (
            self.power_inl[0].E.val_SI - self.power_outl[0].E.val_SI
        )
