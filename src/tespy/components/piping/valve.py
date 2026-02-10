# -*- coding: utf-8

"""Module of class Valve.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/piping.py

SPDX-License-Identifier: MIT
"""

import numpy as np

from tespy.components.component import Component
from tespy.components.component import component_registry
from tespy.tools import logger
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentMandatoryConstraints as dc_cmc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import GroupedComponentProperties as dc_gcp
from tespy.tools.data_containers import SimpleDataContainer as dc_simple


@component_registry
class Valve(Component):
    r"""
    The Valve throttles a fluid without changing enthalpy.

    **Mandatory Equations**

    - fluid: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`
    - mass flow: :py:meth:`tespy.components.component.Component.variable_equality_structure_matrix`

    **Optional Equations**

    - :py:meth:`tespy.components.component.Component.dp_structure_matrix`
    - :py:meth:`tespy.components.component.Component.pr_structure_matrix`
    - :py:meth:`tespy.components.component.Component.zeta_func`
    - :py:meth:`tespy.components.piping.valve.Valve.dp_char_func`

    Inlets/Outlets

    - in1
    - out1

    Image

    .. image:: /api/_images/Valve.svg
       :alt: flowsheet of the valve
       :align: center
       :class: only-light

    .. image:: /api/_images/Valve_darkmode.svg
       :alt: flowsheet of the valve
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

    pr : float, dict, :code:`"var"`
        Outlet to inlet pressure ratio, :math:`pr/1`

    zeta : float, dict, :code:`"var"`
        Geometry independent friction coefficient,
        :math:`\frac{\zeta}{D^4}/\frac{1}{\text{m}^4}`.

    dp_char : tespy.tools.characteristics.CharLine, dict
        Characteristic line for difference pressure to mass flow.

    Example
    -------
    A mass flow of 1 kg/s methane is throttled from 80 bar to 15 bar in a
    valve. The inlet temperature is at 50 °C. It is possible to determine the
    outlet temperature as the throttling does not change enthalpy.

    >>> from tespy.components import Sink, Source, Valve
    >>> from tespy.connections import Connection
    >>> from tespy.networks import Network
    >>> import os
    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC"
    ... })
    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> v = Valve('valve')
    >>> so_v = Connection(so, 'out1', v, 'in1')
    >>> v_si = Connection(v, 'out1', si, 'in1')
    >>> nw.add_conns(so_v, v_si)
    >>> v.set_attr(offdesign=['zeta'])
    >>> so_v.set_attr(fluid={'CH4': 1}, m=1, T=50, p=80, design=['m'])
    >>> v_si.set_attr(p=15)
    >>> nw.solve('design')
    >>> nw.save('tmp.json')
    >>> round(v_si.T.val, 1)
    26.3
    >>> round(v.pr.val, 3)
    0.188

    The simulation determined the area independent zeta value
    :math:`\frac{\zeta}{D^4}`. This zeta remains constant if the cross
    sectional area of the valve opening does not change. Using the zeta value
    we can determine the pressure ratio at a different feed pressure.

    >>> so_v.set_attr(p=70)
    >>> nw.solve('offdesign', design_path='tmp.json')
    >>> round(so_v.m.val, 1)
    0.9
    >>> round(v_si.T.val, 1)
    30.0
    >>> os.remove('tmp.json')

    You can also specify the flow coefficient of the valve :code:`Kv` which is
    used in context of liquids. For this there are several methods available:

    - direct specification with :code:`Kv`
    - lookup table specification :math:`Kv=f\left(opening\right)` with
      :code:`Kv_char` and :code:`opening`
    - arbitrary function specification :math:`Kv=f\left(opening, params\right)`
      with :code:`Kv_analytical` and :code:`opening`

    >>> nw = Network(iterinfo=False)
    >>> nw.units.set_defaults(**{
    ...     "pressure": "bar", "temperature": "degC"
    ... })
    >>> so = Source('source')
    >>> si = Sink('sink')
    >>> v = Valve('valve')
    >>> so_v = Connection(so, 'out1', v, 'in1')
    >>> v_si = Connection(v, 'out1', si, 'in1')
    >>> nw.add_conns(so_v, v_si)
    >>> so_v.set_attr(fluid={'water': 1}, T=50, p=5)
    >>> v_si.set_attr(p=4)

    First we specify Kv:

    >>> v.set_attr(Kv=10)
    >>> nw.solve('design')
    >>> round(so_v.v.val, 4)
    0.0028

    Then, for example an analytical function:

    >>> def analytical(opening, *params):
    ...     return params[0] * opening
    >>> v.set_attr(
    ...     Kv=None,
    ...     opening=0.5,
    ...     Kv_analytical={"method": analytical, "params": [10]}
    ... )
    >>> nw.solve("design")
    >>> round(v.Kv.val_SI, 1)
    5.0

    Or, use the :code:`Kv_char`, which is a :code:`CharLine` instance:

    >>> from tespy.tools import CharLine
    >>> kv_data = np.array([
    ...     0.09,0.63,1.1,2.1,3.1,4.2,5.2,6.2,7.2,8.2,9.2,10.3,11.3
    ... ])
    >>> opening_data = np.array([0,5,10,20,30,40,50,60,70,80,90,100,110]) / 100
    >>> Kv_char = {
    ...     "char_func": CharLine(x=opening_data, y=kv_data) , "is_set": True
    ... }
    >>> v.set_attr(Kv_char=Kv_char, Kv_analytical=None)
    >>> nw.solve("design")
    >>> round(v.Kv.val, 1)
    5.2
    """
    def get_parameters(self):
        return {
            'pr': dc_cp(
                min_val=1e-4, max_val=1, num_eq_sets=1,
                structure_matrix=self.pr_structure_matrix,
                func_params={'pr': 'pr'},
                quantity="ratio",
                description="outlet ot inlet pressure ratio"
            ),
            'dp': dc_cp(
                min_val=0,
                num_eq_sets=1,
                structure_matrix=self.dp_structure_matrix,
                func_params={"inconn": 0, "outconn": 0, "dp": "dp"},
                quantity="pressure",
                description="inlet to outlet absolute pressure change"
            ),
            'zeta': dc_cp(
                min_val=0, max_val=1e15, num_eq_sets=1,
                func=self.zeta_func,
                dependents=self.zeta_dependents,
                func_params={'zeta': 'zeta'},
                description="non-dimensional friction coefficient for pressure loss calculation"
            ),
            'dp_char': dc_cc(
                param='m', num_eq_sets=1,
                dependents=self.dp_char_dependents,
                func=self.dp_char_func,
                char_params={'type': 'abs'},
                description="inlet to outlet absolute pressure change as function of mass flow lookup table"
            ),
            'Kv': dc_cp(
                min_val=0, max_val=1e15, num_eq_sets=1,
                func=self.Kv_func,
                dependents=self.Kv_dependents,
                description="flow coefficient in m3/h",
            ),
            'Kv_char': dc_cc(
                description="lookup-table data for flow coefficient as function of opening"
            ),
            'opening': dc_cp(
                # opening can be more than 100 % sometimes
                min_val=0, max_val=1.1,
                _potential_var=True,
                description="opening ratio of the valve",
                quantity="ratio"
            ),
            'Kv_char_group': dc_gcp(
                num_eq_sets=1,
                elements=["Kv_char", "opening"],
                func=self.Kv_char_func,
                dependents=self.Kv_char_dependents,
                description="equation for flow coefficient over opening"
            ),
            'Kv_analytical': dc_simple(
                description=(
                    "fitting parameters and method for the analytical Kv "
                    "evaluation provided in a dictionary with keys 'method' "
                    "(callable) and 'params' (list)"
                )
            ),
            'Kv_char_analytical_group': dc_gcp(
                num_eq_sets=1,
                elements=["Kv_analytical", "opening"],
                func=self.Kv_char_analytical_func,
                dependents=self.Kv_char_analytical_dependents
            ),
        }

    def get_mandatory_constraints(self):
        constraints = super().get_mandatory_constraints()
        constraints.update({
            'enthalpy_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': 1,
                'func_params': {'variable': 'h'},
                "description": "equation for enthalpy equality"
            })
        })
        return constraints

    def get_bypass_constraints(self):
        return {
            'mass_flow_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': self.num_i,
                'func_params': {'variable': 'm'}
            }),
            'pressure_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': self.num_i,
                'func_params': {'variable': 'p'}
            }),
            'enthalpy_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': self.num_i,
                'func_params': {'variable': 'h'}
            }),
            'fluid_constraints': dc_cmc(**{
                'structure_matrix': self.variable_equality_structure_matrix,
                'num_eq_sets': self.num_i,
                'func_params': {'variable': 'fluid'}
            })
        }

    @staticmethod
    def inlets():
        return ['in1']

    @staticmethod
    def outlets():
        return ['out1']

    def dp_char_func(self):
        r"""
        Equation for characteristic line of difference pressure to mass flow.

        Returns
        -------
        float
            Residual value of equation.

            .. math::

                0=p_\mathrm{in}-p_\mathrm{out}-f\left( expr \right)
        """
        p = self.dp_char.param
        expr = self.get_char_expr(p, **self.dp_char.char_params)
        if not expr:
            msg = (
                "Please choose a valid parameter for the usage of the "
                f"'dp_char_func' of the component {self.label}."
            )
            logger.error(msg)
            raise ValueError(msg)

        return (
            self.inl[0].p.val_SI - self.outl[0].p.val_SI
            - self.dp_char.char_func.evaluate(expr)
        )

    def dp_char_dependents(self):
        dependents = [
            self.inl[0].m,
            self.inl[0].p,
            self.outl[0].p,
        ]
        if self.dp_char.param == 'v':
            dependents += [self.inl[0].h]
        return dependents

    def _Kv_eq(self, Kv):
        # 1000 * delta p (bar) is 1000 * delta p / 100000, simplified to
        # delta p / 100
        # (vol * m * 3600) ** 2 / vol simplified to
        # (m * 3600) ** 2 * vol
        return (
            Kv ** 2 * (self.inl[0].p.val_SI - self.outl[0].p.val_SI) / 1e2
            - self.inl[0].calc_vol() * (self.inl[0].m.val_SI * 3600) ** 2
        )

    def Kv_func(self):
        r"""
        Equation for Kv value of a Valve

        The equation is as follows:

        .. math::

            K_v=\dot V \cdot \sqrt{\frac{\rho}{1000\cdot \Delta p}}

        The residual is reformulated as below:

        Returns
        -------
        float
            Residual value of equation.

            .. math::

                0=K_v ^ 2 \cdot \frac{\Delta p}{100}
                -\frac{\left(3600 \cdot \dot m \right) ^ 2}{\rho}
        """
        Kv = self.Kv.val_SI
        return self._Kv_eq(Kv)

    def Kv_dependents(self):
        return [self.inl[0].m, self.inl[0].p, self.inl[0].h, self.outl[0].p]

    def Kv_char_func(self):
        r"""
        Equation for Kv characteristic of a Valve opening
        :math:`K_v=f\left(opening\right)`

        Kv is determined from the degree of opening with a lookup table, the
        Kv equation is then applied:

        .. math::

            K_v=\dot V \cdot \sqrt{\frac{\rho}{1000\cdot \Delta p}}

        The residual is reformulated as below:

        Returns
        -------
        float
            Residual value of equation.

            .. math::

                0=K_v ^ 2 \cdot \frac{\Delta p}{100}
                -\frac{\left(3600 \cdot \dot m \right) ^ 2}{\rho}
        """
        Kv = self.Kv_char.char_func.evaluate(self.opening.val_SI)
        return self._Kv_eq(Kv)

    def Kv_char_dependents(self):
        return [
            self.inl[0].m, self.inl[0].p, self.inl[0].h, self.outl[0].p,
            self.opening
        ]

    def Kv_char_analytical_func(self):
        r"""
        Equation for Kv characteristic of a Valve opening
        :math:`K_v=f\left(opening\right)`

        Kv is determined from the method supplied by the user in the
        :code:`Kv_analytical` specification and the additional parameters next
        to the opening. The method must accept parameters in the following way:
        :code:`method(opening, *params)`

        .. math::

            K_v=\dot V \cdot \sqrt{\frac{\rho}{1000\cdot \Delta p}}

        The residual is reformulated as below:

        Returns
        -------
        float
            Residual value of equation.

            .. math::

                0=K_v ^ 2 \cdot \frac{\Delta p}{100}
                -\frac{\left(3600 \cdot \dot m \right) ^ 2}{\rho}
        """
        params = self.Kv_analytical.val["params"]
        method = self.Kv_analytical.val["method"]
        opening = self.opening.val_SI
        Kv = method(opening, *params)
        return self._Kv_eq(Kv)

    def Kv_char_analytical_dependents(self):
        return [
            self.inl[0].m, self.inl[0].p, self.inl[0].h, self.outl[0].p,
            self.opening
        ]

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        i = self.inl[0]
        o = self.outl[0]
        self.pr.val_SI = o.p.val_SI / i.p.val_SI
        self.dp.val_SI = i.p.val_SI - o.p.val_SI
        self.zeta.val_SI = self.calc_zeta(i, o)
        if self.dp.val_SI > 0 and i.calc_phase() == "l":
            self.Kv.val_SI = (
                i.v.val_SI * 3600 * (
                    100 / (i.vol.val_SI * self.dp.val_SI)
                ) ** 0.5
            )

    def entropy_balance(self):
        r"""
        Calculate entropy balance of a valve.

        Note
        ----
        The entropy balance makes the following parameter available:

        .. math::

            \text{S\_irr}=\dot{m} \cdot \left(s_\mathrm{out}-s_\mathrm{in}
            \right)\\
        """
        self.S_irr = self.inl[0].m.val_SI * (
            self.outl[0].s.val_SI - self.inl[0].s.val_SI
        )

    def exergy_balance(self, T0):
        r"""
        Calculate exergy balance of a valve.

        Parameters
        ----------
        T0 : float
            Ambient temperature T0 / K.

        Note
        ----
        .. math::

            \dot{E}_\mathrm{P} =
            \begin{cases}
            \text{not defined (nan)} & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
            \dot{E}_\mathrm{out}^\mathrm{T}
            & T_\mathrm{in} > T_0 \geq T_\mathrm{out}\\
            \dot{E}_\mathrm{out}^\mathrm{T} - \dot{E}_\mathrm{in}^\mathrm{T}
            & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
            \end{cases}

            \dot{E}_\mathrm{F} =
            \begin{cases}
            \dot{E}_\mathrm{in}^\mathrm{PH} - \dot{E}_\mathrm{out}^\mathrm{PH}
            & T_\mathrm{in}, T_\mathrm{out} \geq T_0\\
            \dot{E}_\mathrm{in}^\mathrm{T} + \dot{E}_\mathrm{in}^\mathrm{M}-
            \dot{E}_\mathrm{out}^\mathrm{M}
            & T_\mathrm{in} > T_0 \geq T_\mathrm{out}\\
            \dot{E}_\mathrm{in}^\mathrm{M} - \dot{E}_\mathrm{out}^\mathrm{M}
            & T_0 \geq T_\mathrm{in}, T_\mathrm{out}\\
            \end{cases}
        """
        if self.inl[0].T.val_SI > T0 and self.outl[0].T.val_SI > T0:
            self.E_P = np.nan
            self.E_F = self.inl[0].Ex_physical - self.outl[0].Ex_physical
        elif self.outl[0].T.val_SI <= T0 and self.inl[0].T.val_SI > T0:
            self.E_P = self.outl[0].Ex_therm
            self.E_F = self.inl[0].Ex_therm + (
                self.inl[0].Ex_mech - self.outl[0].Ex_mech)
        elif self.inl[0].T.val_SI <= T0 and self.outl[0].T.val_SI <= T0:
            self.E_P = self.outl[0].Ex_therm - self.inl[0].Ex_therm
            self.E_F = self.inl[0].Ex_mech - self.outl[0].Ex_mech
        else:
            msg = ('Exergy balance of a valve, where outlet temperature is '
                   'larger than inlet temperature is not implmented.')
            logger.warning(msg)
            self.E_P = np.nan
            self.E_F = np.nan

        self.E_bus = {
            "chemical": np.nan, "physical": np.nan, "massless": np.nan
        }
        if np.isnan(self.E_P):
            self.E_D = self.E_F
        else:
            self.E_D = self.E_F - self.E_P
        self.epsilon = self._calc_epsilon()

    def get_plotting_data(self):
        """Generate a dictionary containing FluProDia plotting information.

        Returns
        -------
        data : dict
            A nested dictionary containing the keywords required by the
            :code:`calc_individual_isoline` method of the
            :code:`FluidPropertyDiagram` class. First level keys are the
            connection index ('in1' -> 'out1', therefore :code:`1` etc.).
        """
        return {
            1: {
                'isoline_property': 'h',
                'isoline_value': self.inl[0].h.val,
                'isoline_value_end': self.outl[0].h.val,
                'starting_point_property': 'v',
                'starting_point_value': self.inl[0].vol.val,
                'ending_point_property': 'v',
                'ending_point_value': self.outl[0].vol.val
            }
        }
