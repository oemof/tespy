# -*- coding: utf-8

"""Module of class CombustionEngine.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/combustion/engine.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.combustion.base import CombustionChamber
from tespy.components.component import component_registry
from tespy.tools import logger
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import s_mix_ph
from tespy.tools.fluid_properties import s_mix_pT
from tespy.tools.helpers import _get_dependents
from tespy.tools.helpers import _numeric_deriv


@component_registry
class CombustionEngine(CombustionChamber):
    r"""
    An internal combustion engine supplies power and heat cogeneration.

    The combustion engine produces power and heat in cogeneration from fuel
    combustion. The combustion properties are identical to the combustion
    chamber. Thermal input and power output, heat output and heat losses are
    linked with an individual characteristic line for each property.

    **Mandatory Equations**

    - :py:meth:`tespy.components.combustion.base.CombustionChamber.mass_flow_func`
    - :py:meth:`tespy.components.combustion.base.CombustionChamber.combustion_pressure_structure_matrix`
    - :py:meth:`tespy.components.combustion.base.CombustionChamber.stoichiometry`
    - :py:meth:`tespy.components.combustion.engine.CombustionEngine.energy_balance_func`
    - :py:meth:`tespy.components.combustion.engine.CombustionEngine.tiP_char_func`
    - :py:meth:`tespy.components.combustion.engine.CombustionEngine.Q1_char_func`
    - :py:meth:`tespy.components.combustion.engine.CombustionEngine.Q2_char_func`
    - :py:meth:`tespy.components.combustion.engine.CombustionEngine.Qloss_char_func`

    - for each cooling loop:

      - mass flow: :py:meth:`tespy.components.combustion.engine.CombustionEngine.variable_equality_structure_matrix`
      - fluid: :py:meth:`tespy.components.combustion.engine.CombustionEngine.variable_equality_structure_matrix`

    **Optional Equations**

    - :py:meth:`tespy.components.combustion.base.CombustionChamber.lambda_func`
    - :py:meth:`tespy.components.combustion.base.CombustionChamber.ti_func`
    - :py:meth:`tespy.components.combustion.engine.CombustionEngine.Q1_func`
    - :py:meth:`tespy.components.combustion.engine.CombustionEngine.Q2_func`

    - for each cooling loop:

      - :py:meth:`tespy.components.component.Component.dp_structure_matrix`
      - :py:meth:`tespy.components.component.Component.pr_structure_matrix`
      - :py:meth:`tespy.components.component.Component.zeta_func`

    Available fuels

    - methane, ethane, propane, butane, hydrogen, carbon monoxide, nDodecane

    Inlets/Outlets

    - in1, in2 (cooling water), in3, in4 (air and fuel)
    - out1, out2 (cooling water), out3 (flue gas)

    Power outlets

    - power

    Image

    .. image:: /api/_images/CombustionEngine.svg
       :alt: flowsheet of the combustion engine
       :align: center
       :class: only-light

    .. image:: /api/_images/CombustionEngine_darkmode.svg
       :alt: flowsheet of the combustion engine
       :align: center
       :class: only-dark

    .. note::

        The fuel and the air components can be connected to either of the
        inlets.

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

    lamb : float, dict
        Air to stoichiometric air ratio, :math:`\lambda/1`.

    ti : float, dict
        Thermal input, (:math:`{LHV \cdot \dot{m}_f}`), :math:`ti/\text{W}`.

    P : float, dict, :code:`"var"`
        Power output, :math:`P/\text{W}`.

    Q1 : float, dict
        Heat output 1, :math:`\dot Q/\text{W}`.

    Q2 : float, dict
        Heat output 2, :math:`\dot Q/\text{W}`.

    Qloss : float, dict, :code:`"var"`
        Heat loss, :math:`\dot Q_{loss}/\text{W}`.

    pr1 : float, dict, :code:`"var"`
        Pressure ratio heat outlet 1, :math:`pr/1`.

    pr2 : float, dict, :code:`"var"`
        Pressure ratio heat outlet 2, :math:`pr/1`.

    zeta1 : float, dict, :code:`"var"`
        Geometry independent friction coefficient heating loop 1,
        :math:`\zeta_1/\frac{1}{\text{m}^4}`.

    zeta2 : float, dict, :code:`"var"`
        Geometry independent friction coefficient heating loop 2,
        :math:`\zeta_2/\frac{1}{\text{m}^4}`.

    tiP_char : tespy.tools.characteristics.CharLine, dict
        Characteristic line linking fuel input to power output.

    Q1_char : tespy.tools.characteristics.CharLine, dict
        Characteristic line linking heat output 1 to power output.

    Q2_char : tespy.tools.characteristics.CharLine, dict
        Characteristic line linking heat output 2 to power output.

    Qloss_char : tespy.tools.characteristics.CharLine, dict
        Characteristic line linking heat loss to power output.

    eta_mech : float
        Value of internal efficiency of the combustion engine. This value is
        required to determine the (virtual) thermodynamic temperature of heat
        inside the combustion engine for the entropy balance calculation.
        Default value is 0.85.

    Note
    ----
    Parameters available through entropy balance are listed in the respective
    method:

    - :py:meth:`tespy.components.combustion.engine.CombustionEngine.entropy_balance`

    Example
    -------
    The combustion chamber calculates energy input due to combustion as well as
    the flue gas composition based on the type of fuel and the amount of
    oxygen supplied. In this example a mixture of methane, hydrogen and
    carbondioxide is used as fuel. There are two cooling ports, the cooling
    water will flow through them in parallel.

    >>> from tespy.components import (Sink, Source, CombustionEngine, Merge,
    ... Splitter)
    >>> from tespy.connections import Connection, Ref
    >>> from tespy.networks import Network
    >>> import os
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC"
    ... })
    >>> amb = Source('ambient')
    >>> sf = Source('fuel')
    >>> fg = Sink('flue gas outlet')
    >>> cw_in = Source('cooling water inlet')
    >>> sp = Splitter('cooling water splitter', num_out=2)
    >>> me = Merge('cooling water merge', num_in=2)
    >>> cw_out = Sink('cooling water outlet')
    >>> chp = CombustionEngine(label='internal combustion engine')
    >>> amb_comb = Connection(amb, 'out1', chp, 'in3')
    >>> sf_comb = Connection(sf, 'out1', chp, 'in4')
    >>> comb_fg = Connection(chp, 'out3', fg, 'in1')
    >>> nw.add_conns(sf_comb, amb_comb, comb_fg)
    >>> cw_sp = Connection(cw_in, 'out1', sp, 'in1')
    >>> sp_chp1 = Connection(sp, 'out1', chp, 'in1')
    >>> sp_chp2 = Connection(sp, 'out2', chp, 'in2')
    >>> chp1_me = Connection(chp, 'out1', me, 'in1')
    >>> chp2_me = Connection(chp, 'out2', me, 'in2')
    >>> me_cw = Connection(me, 'out1', cw_out, 'in1')
    >>> nw.add_conns(cw_sp, sp_chp1, sp_chp2, chp1_me, chp2_me, me_cw)

    The combustion engine produces a power output of 10 MW the oxygen to
    stoichiometric oxygen ratio is set to 1. Only pressure ratio 1 is set as
    we reconnect both cooling water streams. At the merge all pressure values
    will be identical automatically. Reference the mass flow at the splitter
    to be split in half.

    >>> chp.set_attr(pr1=0.99, P=-10e6, lamb=1.0,
    ... design=['pr1'], offdesign=['zeta1'])
    >>> amb_comb.set_attr(p=5, T=30, fluid={'Ar': 0.0129, 'N2': 0.7553,
    ... 'CO2': 0.0004, 'O2': 0.2314})
    >>> sf_comb.set_attr(T=30, fluid={'CH4': 1})
    >>> cw_sp.set_attr(p=3, T=60, m=50, fluid={'H2O': 1})
    >>> sp_chp2.set_attr(m=Ref(sp_chp1, 1, 0))
    >>> mode = 'design'
    >>> nw.solve(mode=mode)
    >>> nw.save('tmp.json')
    >>> round(chp.ti.val, 0)
    25300000.0
    >>> round(chp.Q1.val, 0)
    -4980000.0
    >>> chp.set_attr(Q1=-4e6, P=None)
    >>> mode = 'offdesign'
    >>> nw.solve(mode=mode, init_path='tmp.json', design_path='tmp.json')
    >>> round(chp.ti.val, 0)
    17794554.0
    >>> round(chp.P.val / chp.P.design, 3)
    0.617
    >>> chp.set_attr(P=chp.P.design * 0.75, Q1=None)
    >>> mode = 'offdesign'
    >>> nw.solve(mode=mode, init_path='tmp.json', design_path='tmp.json')
    >>> round(chp.ti.val, 0)
    20550000.0
    >>> round(chp.P.val / chp.P.design, 3)
    0.75
    >>> os.remove('tmp.json')
    """

    def get_parameters(self):
        params = super().get_parameters()
        params.update({
            'P': dc_cp(_val=-1e6, d=1, max_val=-1, quantity="power"),
            'Q1': dc_cp(
                max_val=-1,
                num_eq_sets=1,
                dependents=self.Q1_dependents,
                func=self.Q1_func,
                quantity="heat"
            ),
            'Q2': dc_cp(
                max_val=-1,
                num_eq_sets=1,
                dependents=self.Q2_dependents,
                func=self.Q2_func,
                quantity="heat"
            ),
            'Qloss': dc_cp(_val=-1e5, d=1, max_val=-1, quantity="heat"),
            'pr1': dc_cp(
                min_val=1e-4, max_val=1, num_eq_sets=1,
                structure_matrix=self.pr_structure_matrix,
                func_params={'pr': 'pr1'},
                quantity="ratio"
            ),
            'pr2': dc_cp(
                min_val=1e-4, max_val=1, num_eq_sets=1,
                structure_matrix=self.pr_structure_matrix,
                func_params={'pr': 'pr2', 'inconn': 1, 'outconn': 1},
                quantity="ratio"
            ),
            'dp1': dc_cp(
                min_val=0, num_eq_sets=1,
                structure_matrix=self.dp_structure_matrix,
                func_params={"inconn": 0, "outconn": 0, "dp": "dp1"},
                quantity="pressure"
            ),
            'dp2': dc_cp(
                min_val=0, num_eq_sets=1,
                structure_matrix=self.dp_structure_matrix,
                func_params={"inconn": 1, "outconn": 1, "dp": "dp2"},
                quantity="pressure"
            ),
            'zeta1': dc_cp(
                min_val=0, max_val=1e15, num_eq_sets=1,
                dependents=self.zeta_dependents,
                func=self.zeta_func,
                func_params={'zeta': 'zeta1'}
            ),
            'zeta2': dc_cp(
                min_val=0, max_val=1e15, num_eq_sets=1,
                dependents=self.zeta_dependents,
                func=self.zeta_func,
                func_params={'zeta': 'zeta2', 'inconn': 1, 'outconn': 1}
            ),
            'tiP_char': dc_cc(),
            'Q1_char': dc_cc(),
            'Q2_char': dc_cc(),
            'Qloss_char': dc_cc(),
            'eta_mech': dc_simple(_val=0.85),
            'T_v_inner': dc_simple()
        })
        return params

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        constraints["energy_balance_constraints"].deriv = self.energy_balance_deriv
        constraints.update({
            'power_constraints': dc_cmc(**{
                'func': self.tiP_char_func,
                'dependents': self.tiP_char_dependents,
                'constant_deriv': False,
                'num_eq_sets': 1,
            }),
            'heat1_constraints': dc_cmc(**{
                'func': self.Q1_char_func,
                'dependents': self.Q1_char_dependents,
                'constant_deriv': False,
                'num_eq_sets': 1,
            }),
            'heat2_constraints': dc_cmc(**{
                'func': self.Q2_char_func,
                'dependents': self.Q2_char_dependents,
                'constant_deriv': False,
                'num_eq_sets': 1,
            }),
            'heatloss_constraints': dc_cmc(**{
                'func': self.Qloss_char_func,
                'dependents': self.Qloss_char_dependents,
                'constant_deriv': False,
                'num_eq_sets': 1,
            }),
            'mass_flow_cooling_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': 2,
                'func_params': {'variable': 'm'}
            }),
            'fluid_cooling_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': 2,
                'func_params': {'variable': 'fluid'}
            })
        })
        if len(self.power_outl) > 0:
            constraints["energy_connector_balance"] = dc_cmc(**{
                "func": self.energy_connector_balance_func,
                "dependents": self.energy_connector_dependents,
                "num_eq_sets": 1
            })
        return constraints

    @staticmethod
    def inlets():
        return ['in1', 'in2', 'in3', 'in4']

    @staticmethod
    def outlets():
        return ['out1', 'out2', 'out3']

    @staticmethod
    def poweroutlets():
        return ["power"]

    def propagate_wrapper_to_target(self, branch):
        inl, _ = self._get_combustion_connections()
        inconn = branch["connections"][-1]
        if inconn in inl:
            if self in branch["components"]:
                return

            outconn = self.outl[2]
        else:
            conn_idx = self.inl.index(inconn)
            outconn = self.outl[conn_idx]

        branch["connections"] += [outconn]
        branch["components"] += [self]

        outconn.target.propagate_wrapper_to_target(branch)

    def get_variables(self):
        if not self.P.is_set:
            self.set_attr(P='var')
            msg = (
                f'The power output of combustion engines ({self.label}) must '
                'either be a set value or part of the system variables. Since '
                'it has not been set to a fixed value it will be added to the '
                'system\'s variables.'
            )
            logger.info(msg)

        if not self.Qloss.is_set:
            self.set_attr(Qloss='var')
            msg = (
                f'The heat loss of combustion engines ({self.label}) must '
                'either be a set value or part of the system variables. Since '
                'it has not been set to a fixed value it will be added to the '
                'system\'s variables.'
            )
            logger.info(msg)

        return super().get_variables()

    def _get_combustion_connections(self):
        return (self.inl[2:], [self.outl[2]])

    def variable_equality_structure_matrix(self, k, **kwargs):
        variable = kwargs.get("variable")
        for count, (i, o) in enumerate(zip(self.inl[:2], self.outl[:2])):
            self._structure_matrix[k + count, i.get_attr(variable).sm_col] = 1
            self._structure_matrix[k + count, o.get_attr(variable).sm_col] = -1

    def energy_connector_balance_func(self):
        return self.power_outl[0].E.val_SI + self.P.val_SI

    def energy_connector_dependents(self):
        return [self.power_outl[0].E, self.P]


    def energy_balance_func(self):
        r"""
        Calculate the energy balance of the combustion engine.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \sum_i \dot{m}_{in,i} \cdot
                \left( h_{in,i} - h_{in,i,ref} \right)\\
                & - \sum_j \dot{m}_{out,3} \cdot
                \left( h_{out,3} - h_{out,3,ref} \right)\\
                & + LHV_{fuel} \cdot
                \left(\sum_i \left(\dot{m}_{in,i} \cdot x_{fuel,i} \right)-
                \dot{m}_{out,3} \cdot x_{fuel,3} \right)\\
                & - \dot{Q}_1 - \dot{Q}_2 - P - \dot{Q}_{loss}\\
                \end{split}\\
                \forall i \in [3,4]

        Note
        ----
        The temperature for the reference state is set to 25 °C, thus
        the water may be liquid. In order to make sure, the state is
        referring to the lower heating value, the necessary enthalpy
        difference for evaporation is added.

        - Reference temperature: 298.15 K.
        - Reference pressure: 1 bar.
        """
        res = super().energy_balance_func()

        # cooling
        for i in range(2):
            res -= self.inl[i].m.val_SI * (
                self.outl[i].h.val_SI - self.inl[i].h.val_SI
            )

        # power output and heat loss
        res += self.P.val_SI + self.Qloss.val_SI

        return res

    def energy_balance_deriv(self, increment_filter, k, dependents=None):
        """
        Calculate partial derivatives of energy balance function.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of equation in Jacobian matrix.
        """
        f = self.energy_balance_func
        # mass flow cooling water
        i0 = self.inl[0]
        i1 = self.inl[1]
        o0 = self.outl[0]
        o1 = self.outl[1]
        if i0.m.is_var:
            if i0.m._reference_container != i1.m._reference_container:
                self.jacobian[k, i0.m.J_col] = -(o0.h.val_SI - i0.h.val_SI)
                self.jacobian[k, i1.m.J_col] = -(o1.h.val_SI - i1.h.val_SI)
            else:
                self.jacobian[k, i0.m.J_col] = (
                    -(o0.h.val_SI - i0.h.val_SI)
                    -(o1.h.val_SI - i1.h.val_SI)
                )

        # mass flow and pressure for combustion reaction
        inl, outl = self._get_combustion_connections()
        variables = _get_dependents([var for c in inl + outl for var in (c.m, c.p)])[0]
        for variable in variables:
            self._partial_derivative(variable, k, f, increment_filter)

        # enthalpy all connections
        for i in self.inl:
            if i.h.is_var:
                self.jacobian[k, i.h.J_col] = i.m.val_SI

        for o in self.outl:
            if o.h.is_var:
                self.jacobian[k, o.h.J_col] = -o.m.val_SI

        # fluid composition
        for c in inl:
            for fl in (self.fuel_list & c.fluid.is_var):
                self.jacobian[k, c.fluid.J_col[fl]] = c.m.val_SI * self.fuels[fl]['LHV']

        c = outl[0]
        for fl in (self.fuel_list & c.fluid.is_var):
            self.jacobian[k, c.fluid.J_col[fl]] = -c.m.val_SI * self.fuels[fl]['LHV']

        # power and heat loss
        if self.P.is_var:
            self.jacobian[k, self.P.J_col] = 1
        if self.Qloss.is_var:
            self.jacobian[k, self.Qloss.J_col] = 1

    def energy_balance_dependents(self):
        inl, outl = self._get_combustion_connections()
        return {
            "scalars": [
                var for c in inl + outl for var in [c.m, c.p, c.h]
            ] + [
                var for c in self.inl[:2] + self.outl[:2] for var in [c.m, c.h]
            ] + [self.P, self.Qloss],
            "vectors": [{
                c.fluid: self.fuel_list & c.fluid.is_var for c in inl + outl
            }]
        }

    def Q1_func(self):
        r"""
        Calculate residual value of primary heat loop function.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_1 \cdot \left(h_{out,1} +
                h_{in,1} \right) + \dot{Q}_1
        """
        i = self.inl[0]
        o = self.outl[0]
        return i.m.val_SI * (o.h.val_SI - i.h.val_SI) + self.Q1.val_SI

    def Q1_dependents(self):
        return [self.inl[0].m, self.inl[0].h, self.outl[0].h]

    def Q2_func(self):
        r"""
        Calculate residual value of secondary heat loop function.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = \dot{m}_2 \cdot \left(h_{out,2} - h_{in,2} \right) +
                \dot{Q}_2
        """
        i = self.inl[1]
        o = self.outl[1]
        return i.m.val_SI * (o.h.val_SI - i.h.val_SI) + self.Q2.val_SI

    def Q2_dependents(self):
        return [self.inl[1].m, self.inl[1].h, self.outl[1].h]

    def tiP_char_func(self):
        r"""
        Calculate the relation of output power and thermal input.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                0 = P \cdot f_{TI}\left(\frac{P}{P_{design}}\right)+ LHV \cdot
                \left[\sum_i \left(\dot{m}_{in,i} \cdot
                x_{f,i}\right) - \dot{m}_{out,3} \cdot x_{f,3} \right]
                \; \forall i \in [3,4]
        """
        expr = self._calc_char_expr()
        return (
            self.calc_ti()
            + self.tiP_char.char_func.evaluate(expr) * self.P.val_SI
        )

    def tiP_char_dependents(self):
        inl, outl = self._get_combustion_connections()
        return {
            "scalars": [c.m for c in inl + outl] + [self.P],
            "vectors": [{
                c.fluid: self.fuel_list & c.fluid.is_var for c in inl + outl
            }]
        }

    def Q1_char_func(self):
        r"""
        Calculate the relation of heat output 1 and thermal input.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \dot{m}_1 \cdot \left(h_{out,1} - h_{in,1} \right) \cdot
                f_{TI}\left(\frac{P}{P_{design}}\right) \\
                & - LHV \cdot \left[\sum_i
                \left(\dot{m}_{in,i} \cdot x_{f,i}\right) -
                \dot{m}_{out,3} \cdot x_{f,3} \right] \cdot
                f_{Q1}\left(\frac{P}{P_{ref}}\right)\\
                \end{split}\\
                \forall i \in [3,4]
        """
        i = self.inl[0]
        o = self.outl[0]

        expr = self._calc_char_expr()
        return (
            self.calc_ti() * self.Q1_char.char_func.evaluate(expr)
            - self.tiP_char.char_func.evaluate(expr) * i.m.val_SI
            * (o.h.val_SI - i.h.val_SI)
        )

    def Q1_char_dependents(self):
        inl, outl = self._get_combustion_connections()
        return {
            "scalars": [
                c.m for c in inl + outl + [self.inl[0]]
            ] + [c.h for c in [self.inl[0], self.outl[0]]] + [self.P],
            "vectors": [{
                c.fluid: self.fuel_list & c.fluid.is_var for c in inl + outl
            }]
        }

    def Q2_char_func(self):
        r"""
        Calculate the relation of heat output 2 and thermal input.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \dot{m}_2 \cdot \left(h_{out,2} - h_{in,2} \right) \cdot
                f_{TI}\left(\frac{P}{P_{design}}\right) \\
                & - LHV \cdot \left[\sum_i
                \left(\dot{m}_{in,i} \cdot x_{f,i}\right) -
                \dot{m}_{out,3} \cdot x_{f,3} \right] \cdot
                f_{Q2}\left(\frac{P}{P_{ref}}\right)\\
                \end{split}\\
                \forall i \in [3,4]
        """
        i = self.inl[1]
        o = self.outl[1]

        expr = self._calc_char_expr()
        return (
            self.calc_ti() * self.Q2_char.char_func.evaluate(expr)
            - self.tiP_char.char_func.evaluate(expr) * i.m.val_SI
            * (o.h.val_SI - i.h.val_SI)
        )

    def Q2_char_dependents(self):
        inl, outl = self._get_combustion_connections()
        return {
            "scalars": [
                c.m for c in inl + outl + [self.inl[1]]
            ] + [c.h for c in [self.inl[1], self.outl[1]]] + [self.P],
            "vectors": [{
                c.fluid: self.fuel_list & c.fluid.is_var for c in inl + outl
            }]
        }

    def Qloss_char_func(self):
        r"""
        Calculate the relation of heat loss and thermal input.

        Returns
        -------
        residual : float
            Residual value of equation.

            .. math::

                \begin{split}
                0 = & \dot{Q}_{loss} \cdot
                f_{TI}\left(\frac{P}{P_{design}}\right) \\
                & + LHV \cdot \left[\sum_i
                \left(\dot{m}_{in,i} \cdot x_{f,i}\right) -
                \dot{m}_{out,3} \cdot x_{f,3} \right] \cdot
                f_{QLOSS}\left(\frac{P}{P_{ref}}\right)\\
                \end{split}\\
                \forall i \in [3,4]
        """
        expr = self._calc_char_expr()
        return (
            self.calc_ti() * self.Qloss_char.char_func.evaluate(expr)
            + self.tiP_char.char_func.evaluate(expr) * self.Qloss.val_SI
        )

    def Qloss_char_dependents(self):
        inl, outl = self._get_combustion_connections()
        return {
            "scalars": [c.m for c in inl + outl] + [self.P, self.Qloss],
            "vectors": [{
                c.fluid: self.fuel_list & c.fluid.is_var for c in inl + outl
            }]
        }

    def calc_P(self):
        r"""
        Calculate the power output of the combustion engine.

        Returns
        -------
        P : float
            Power output.

            .. math::

                P = -\frac{LHV \cdot \dot{m}_{f}}
                {f_{TI}\left(\frac{P}{P_{ref}}\right)}

        """
        expr = self._calc_char_expr()
        return -self.calc_ti() / self.tiP_char.char_func.evaluate(expr)

    def calc_Qloss(self):
        r"""
        Calculate the heat loss of the combustion engine.

        Returns
        -------
        Qloss : float
            Heat loss.

            .. math::

                \dot{Q}_{loss} = -\frac{LHV \cdot \dot{m}_{f} \cdot
                f_{QLOSS}\left(\frac{P}{P_{ref}}\right)}
                {f_{TI}\left(\frac{P}{P_{ref}}\right)}
        """
        expr = self._calc_char_expr()
        return (
            -self.calc_ti() * self.Qloss_char.char_func.evaluate(expr)
            / self.tiP_char.char_func.evaluate(expr)
        )

    def _calc_char_expr(self):
        if np.isnan(self.P.design):
            return 1
        else:
            return self.P.val_SI / self.P.design

    def bus_func(self, bus):
        r"""
        Calculate the value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        residual : float
            Value of energy transfer :math:`\dot{E}`. This value is passed to
            :py:meth:`tespy.components.component.Component.calc_bus_value`
            for value manipulation according to the specified characteristic
            line of the bus.

            .. math::

                \dot{E} = \begin{cases}
                LHV \cdot \dot{m}_{f} & \text{key = 'TI'}\\
                P & \text{key = 'P'}\\
                -\dot{m}_1 \cdot \left( h_{1,out} - h_{1,in} \right)
                -\dot{m}_2 \cdot \left( h_{2,out} - h_{2,in} \right) &
                \text{key = 'Q'}\\
                -\dot{m}_1 \cdot \left( h_{1,out} - h_{1,in} \right) &
                \text{key = 'Q1'}\\
                -\dot{m}_2 \cdot \left( h_{2,out} - h_{2,in} \right) &
                \text{key = 'Q2'}\\
                \dot{Q}_{loss} & \text{key = 'Qloss'}
                \end{cases}
        """
        ######################################################################
        # value for bus parameter of thermal input (TI)
        if bus['param'] == 'TI':
            val = self.calc_ti()

        ######################################################################
        # value for bus parameter of power output (P)
        elif bus['param'] == 'P':
            val = self.calc_P()

        ######################################################################
        # value for bus parameter of total heat production (Q)
        elif bus['param'] == 'Q':
            val = 0
            for j in range(2):
                i = self.inl[j]
                o = self.outl[j]
                val -= i.m.val_SI * (o.h.val_SI - i.h.val_SI)

        ######################################################################
        # value for bus parameter of heat production 1 (Q1)
        elif bus['param'] == 'Q1':
            i = self.inl[0]
            o = self.outl[0]
            val = -i.m.val_SI * (o.h.val_SI - i.h.val_SI)

        ######################################################################
        # value for bus parameter of heat production 2 (Q2)
        elif bus['param'] == 'Q2':
            i = self.inl[1]
            o = self.outl[1]
            val = -i.m.val_SI * (o.h.val_SI - i.h.val_SI)

        ######################################################################
        # value for bus parameter of heat loss (Qloss)
        elif bus['param'] == 'Qloss':
            val = self.calc_Qloss()

        ######################################################################
        # missing/invalid bus parameter
        else:
            msg = (
                f'The parameter {bus["param"]} is not a valid parameter for a '
                f'component of type {self.__class__.__name__}.'
            )
            logger.error(msg)
            raise ValueError(msg)

        return val

    def bus_deriv(self, bus):
        r"""
        Calculate the matrix of partial derivatives of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        deriv : ndarray
            Matrix of partial derivatives.
        """
        f = self.calc_bus_value
        b = bus.comps.loc[self]

        ######################################################################
        # derivatives for bus parameter of thermal input (TI)
        if b['param'] == 'TI':
            # ti deriv is identical to super ti deriv
            super().bus_deriv(bus)

        ######################################################################
        # derivatives for bus parameter of power production (P) or
        # heat loss (Qloss)
        elif b['param'] == 'P' or b['param'] == 'Qloss':
            super().bus_deriv(bus)

            # variable power
            if self.P.is_var:
                if self.P.J_col not in bus.jacobian:
                    bus.jacobian[self.P.J_col] = 0
                bus.jacobian[self.P.J_col] -= _numeric_deriv(self.P._reference_container, f, bus=bus)

        ######################################################################
        # derivatives for bus parameter of total heat production (Q)
        elif b['param'] == 'Q':
            for i, o in  zip(self.inl[:2], self.outl[:2]):
                if i.m.is_var:
                    if i.m.J_col not in bus.jacobian:
                        bus.jacobian[i.m.J_col] = 0
                    bus.jacobian[i.m.J_col] -= _numeric_deriv(i.m._reference_container, f, bus=bus)
                if i.h.is_var:
                    if i.h.J_col not in bus.jacobian:
                        bus.jacobian[i.h.J_col] = 0
                    bus.jacobian[i.h.J_col] -= _numeric_deriv(i.h._reference_container, f, bus=bus)

                if o.h.is_var:
                    if o.h.J_col not in bus.jacobian:
                        bus.jacobian[o.h.J_col] = 0
                    bus.jacobian[o.h.J_col] -= _numeric_deriv(o.h._reference_container, f, bus=bus)

        ######################################################################
        # derivatives for bus parameter of heat production 1 and 2 (Q1, Q2)
        elif b['param'] in ['Q1', 'Q2']:
            i = self.inl[int(b["param"][-1]) - 1]
            o = self.outl[int(b["param"][-1]) - 1]

            if i.m.is_var:
                if i.m.J_col not in bus.jacobian:
                    bus.jacobian[i.m.J_col] = 0
                bus.jacobian[i.m.J_col] -= _numeric_deriv(i.m._reference_container, f, bus=bus)
            if i.h.is_var:
                if i.h.J_col not in bus.jacobian:
                    bus.jacobian[i.h.J_col] = 0
                bus.jacobian[i.h.J_col] -= _numeric_deriv(i.h._reference_container, f, bus=bus)

            if o.h.is_var:
                if o.h.J_col not in bus.jacobian:
                    bus.jacobian[o.h.J_col] = 0
                bus.jacobian[o.h.J_col] -= _numeric_deriv(o.h._reference_container, f, bus=bus)

    @staticmethod
    def initialise_source(c, key):
        r"""
        Return a starting value for pressure and enthalpy at outlet.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                5 \cdot 10^5 & \text{key = 'p'}\\
                10^6 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            return 10e5

    @staticmethod
    def initialise_target(c, key):
        r"""
        Return a starting value for pressure and enthalpy at inlet.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to perform initialisation on.

        key : str
            Fluid property to retrieve.

        Returns
        -------
        val : float
            Starting value for pressure/enthalpy in SI units.

            .. math::

                val = \begin{cases}
                5 \cdot 10^5 & \text{key = 'p'}\\
                5 \cdot 10^5 & \text{key = 'h'}
                \end{cases}
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            return 5e5

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        # Q, pr and zeta
        for i in range(2):
            self.get_attr(f'Q{i + 1}').val_SI = -self.inl[i].m.val_SI * (
                self.outl[i].h.val_SI - self.inl[i].h.val_SI
            )
            self.get_attr(f'dp{i + 1}').val_SI = (
                self.inl[i].p.val_SI - self.outl[i].p.val_SI
            )
            self.get_attr(f'zeta{i + 1}').val_SI = self.calc_zeta(
                self.inl[i], self.outl[i]
            )

        self.P.val_SI = self.calc_P()
        self.Qloss.val_SI = self.calc_Qloss()

        super().calc_parameters()

    def check_parameter_bounds(self):
        r"""Check parameter value limits."""
        _no_limit_violations = super().check_parameter_bounds()
        # get bound errors for characteristic lines
        expr = self._calc_char_expr()
        self.tiP_char.char_func.get_domain_errors(expr, self.label)
        self.Qloss_char.char_func.get_domain_errors(expr, self.label)
        self.Q1_char.char_func.get_domain_errors(expr, self.label)
        self.Q2_char.char_func.get_domain_errors(expr, self.label)

        return _no_limit_violations

    def entropy_balance(self):
        r"""
        Calculate entropy balance of combustion engine.

        For the entropy balance of a combustion engine two additional
        parameters need to be specified:

        - virtual inner temperature :code:`T_v_inner` that is used to determine
          the entropy of heat transferred from the hot side.
        - mechanical efficiency :code:`eta_mech` describing the ratio of power
          output :code:`P` to reversible power of the motor
          :cite:`Zahoransky2019`. It is used to determine the irreversibilty
          inside the motor.

          .. math::

              P_\mathrm{irr,inner}=\left(1 - \frac{1}{\eta_\mathrm{mech}}
              \right) \cdot P

        The default values are:

        - :code:`T_v_inner`: flue gas temperature (result of calculation)
        - :code:`eta_mech`: 0.85

        Note
        ----
        The entropy balance makes the following parameter available:

        - :code:`T_mcomb`: Thermodynamic temperature of heat of combustion
        - :code:`S_comb`: Entropy production due to combustion
        - :code:`T_mQ1`: Thermodynamic temperature of heat at cold side of
          heater 1
        - :code:`S_Q11`: Entropy transport at hot side of heater 1
        - :code:`S_Q12`: Entropy transport at cold side of heater 1
        - :code:`S_Q1irr`: Entropy production due to heat transfer at heater 1
        - :code:`S_irr1`: Entropy production due to pressure losses at heater 1
        - :code:`T_mQ2`: Thermodynamic temperature of heat at cold side of
          heater 2
        - :code:`S_Q21`: Entropy transport at hot side of heater 2
        - :code:`S_Q22`: Entropy transport at cold side of heater 2
        - :code:`S_Q2irr`: Entropy production due to heat transfer at heater 2
        - :code:`S_irr2`: Entropy production due to pressure losses at heater 2
        - :code:`S_irr_i`: Entropy production due to internal irreversibilty
        - :code:`S_Qloss`: Entropy transport with heat loss to ambient
        - :code:`S_Qcomb`: Virtual entropy transport of heat to revert
          combustion gases to reference state
        - :code:`S_irr`: Total entropy production due to irreversibilty

        The methodology for entropy analysis of combustion processes is derived
        from :cite:`Tuschy2001`. Similar to the energy balance of a combustion
        reaction, we need to define the same reference state for the entropy
        balance of the combustion. The temperature for the reference state is
        set to 25 °C and reference pressure is 1 bar. As the water in the flue
        gas may be liquid but the thermodynmic temperature of heat of
        combustion refers to the lower heating value, the water is forced to
        gas at the reference point by considering evaporation.

        - Reference temperature: 298.15 K.
        - Reference pressure: 1 bar.

        .. math::

            \begin{split}
            T_\mathrm{m,comb}= & \frac{\dot{m}_\mathrm{fuel} \cdot LHV}
            {\dot{S}_\mathrm{comb}}\\
            \dot{S}_\mathrm{comb} =&\dot{S}_\mathrm{Q,comb}-\left(
            \dot{S}_\mathrm{Q,11} + \dot{S}_\mathrm{Q,21} +
            \dot{S}_\mathrm{Q,loss} +\dot{S}_\mathrm{irr,i}\right)\\
            \dot{S}_\mathrm{Q,comb}= & \dot{m}_\mathrm{fluegas} \cdot
            \left(s_\mathrm{fluegas}-s_\mathrm{fluegas,ref}\right)\\
            & - \sum_{i=3}^4 \dot{m}_{\mathrm{in,}i} \cdot
            \left( s_{\mathrm{in,}i} - s_{\mathrm{in,ref,}i} \right)\\
            \dot{S}_\mathrm{Q,11}= & \frac{\dot{Q}_1}{T_\mathrm{v,inner}}\\
            \dot{S}_\mathrm{Q,21}= & \frac{\dot{Q}_2}{T_\mathrm{v,inner}}\\
            \dot{S}_\mathrm{Q,loss}= & \frac{\dot{Q}_\mathrm{loss}}
            {T_\mathrm{v,inner}}\\
            \dot{S}_\mathrm{irr,i}= & \frac{\left(1 -
            \frac{1}{\eta_\mathrm{mech}}\right) \cdot P}{T_\mathrm{v,inner}}\\
            T_\mathrm{Q,12} = &\frac{-\dot{Q}_1}{\dot{m}_1 \cdot \left(
            s_\mathrm{out,1} - s_\mathrm{in,1}\right)}\\
            T_\mathrm{Q,22} = &\frac{-\dot{Q}_2}{\dot{m}_2 \cdot \left(
            s_\mathrm{out,2} - s_\mathrm{in,2}\right)}\\
            \dot{S}_\mathrm{irr} = &\sum \dot{S}_\mathrm{irr}\\
            \end{split}\\
        """
        T_ref = 298.15
        p_ref = 1e5
        o = self.outl[2]
        self.S_Qcomb = o.m.val_SI * (
            o.s.val_SI - s_mix_pT(p_ref, T_ref, o.fluid_data, "forced-gas")
        )

        for i in self.inl[2:]:
            self.S_Qcomb -= i.m.val_SI * (
                i.s.val_SI - s_mix_pT(p_ref, T_ref, i.fluid_data, "forced-gas")
            )

        # (virtual) thermodynamic temperature of combustion, use default value
        # if not specified
        if not self.T_v_inner.is_set:
            self.T_v_inner.val = o.T.val_SI

        for i in range(2):
            inl = self.inl[i]
            out = self.outl[i]
            p_star = inl.p.val_SI * (
                self.get_attr('pr' + str(i + 1)).val) ** 0.5
            s_i_star = s_mix_ph(
                p_star, inl.h.val_SI, inl.fluid_data, inl.mixing_rule,
                T0=inl.T.val_SI
            )
            s_o_star = s_mix_ph(
                p_star, out.h.val_SI, out.fluid_data, out.mixing_rule,
                T0=out.T.val_SI
            )

            setattr(
                self, 'S_Q' + str(i + 1) + '2',
                inl.m.val_SI * (s_o_star - s_i_star)
            )
            S_Q = self.get_attr('S_Q' + str(i + 1) + '2')
            setattr(
                self, 'S_irr' + str(i + 1),
                inl.m.val_SI * (out.s.val_SI - inl.s.val_SI) - S_Q
            )
            setattr(
            self, 'T_mQ' + str(i + 1),
            inl.m.val_SI * (out.h.val_SI - inl.h.val_SI) / S_Q
        )

        # internal irreversibilty
        self.P_irr_i = (1 / self.eta_mech.val_SI - 1) * self.P.val_SI

        # internal entropy flow and production
        self.S_Q11 = self.Q1.val_SI / self.T_v_inner.val_SI
        self.S_Q21 = self.Q2.val_SI / self.T_v_inner.val_SI
        self.S_Qloss = self.Qloss.val_SI / self.T_v_inner.val_SI
        self.S_irr_i = self.P_irr_i / self.T_v_inner.val_SI

        # entropy production of heaters due to heat transfer
        self.S_Q1irr = self.S_Q12 - self.S_Q11
        self.S_Q2irr = self.S_Q22 - self.S_Q21

        # calculate entropy production of combustion
        self.S_comb = (
            self.S_Qcomb - self.S_Q11 - self.S_Q21 - self.S_Qloss
            - self.S_irr_i
        )

        # thermodynamic temperature of heat input
        self.T_mcomb = self.calc_ti() / self.S_comb
        # total irreversibilty production
        self.S_irr = (
            self.S_irr_i + self.S_irr2 + self.S_irr1
            + self.S_Q1irr + self.S_Q2irr
        )

    def exergy_balance(self, T0):

        self.E_P = (
            self.outl[2].Ex_physical - (self.inl[3].Ex_physical + self.inl[2].Ex_physical)
            - self.P.val_SI + (self.outl[1] - self.inl[1]) + (self.outl[0] - self.inl[0])
        )
        self.E_F = (
            self.inl[3].Ex_chemical + self.inl[2].Ex_chemical
            - self.outl[2].Ex_chemical
        )
        self.E_D = self.E_F - self.E_P
        self.epsilon = self._calc_epsilon()
        self.E_bus = {
            "chemical": np.nan, "physical": np.nan, "massless": -self.P.val_SI
        }
