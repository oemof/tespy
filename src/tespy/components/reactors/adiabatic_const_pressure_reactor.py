# -*- coding: utf-8

"""Module of class AdiabaticConstPressureReactor.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tespy/components/reactors/water_electrolyzer.py

SPDX-License-Identifier: MIT
"""

import logging

import CoolProp.CoolProp as CP
import numpy as np

from tespy.components.component import Component
from tespy.tools.data_containers import ComponentCharacteristics as dc_cc
from tespy.tools.data_containers import ComponentProperties as dc_cp
from tespy.tools.data_containers import DataContainerSimple as dc_simple
from tespy.tools.document_models import generate_latex_eq
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import dT_mix_dph
from tespy.tools.fluid_properties import dT_mix_pdh
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.global_vars import molar_masses
from tespy.tools.global_vars import fluid_property_data as fpd
from tespy.tools.helpers import TESPyComponentError
from tespy.tools.reaction import Reaction

# %%


class AdiabaticConstPressureReactor(Component):
    r"""
    The adiabatic constant pressure reactor class models arbitrary chemical reactions
    based on a stoichiometric formula.

    **Mandatory Equations**

    - :py:meth:`tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.fluid_func`
    - :py:meth:`tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.mass_flow_func`
    - :py:meth:`tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.reactor_pressure_func`
    - :py:meth:`tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.energy_balance_func`

    Inlets/Outlets

    - in1 (reactants inlet)
    - out1 (products outlet)

    Parameters
    ----------
    label : str
        The label of the component.

    printout : boolean
        Include this component in the network's results printout.

    X : float
        Conversion yield relative to limited reactant :math:`X/1`.

    Example
    -------
    (todo)

    """

    @staticmethod
    def component():
        return "AdiabaticConstPressureReactor"

    def get_variables(self):
        return {
            "X": dc_cp(min_val=0, max_val=1),
            "formula": dc_simple(),
        }

    def get_mandatory_constraints(self):
        return {
            "mass_flow_constraints": {
                "func": self.mass_flow_func,
                "deriv": self.mass_flow_deriv,
                "constant_deriv": True,
                "latex": self.mass_flow_func_doc,
                "num_eq": 1,
            },
            "fluid_constraints": {
                "func": self.fluid_func,
                "deriv": self.fluid_deriv,
                "constant_deriv": True,
                "latex": self.fluid_func_doc,
                "num_eq": len(self.outl[0].fluid.val),
            },
            "energy_balance_constraints": {
                "func": self.energy_balance_func,
                "deriv": self.energy_balance_deriv,
                "constant_deriv": False,
                "latex": self.energy_balance_func_doc,
                "num_eq": 1,
            },
            "reactor_pressure_constraints": {
                "func": self.reactor_pressure_func,
                "deriv": self.reactor_pressure_deriv,
                "constant_deriv": True,
                "latex": self.reactor_pressure_func_doc,
                "num_eq": 1,
            },
        }

    @staticmethod
    def inlets():
        return ["in1"]

    @staticmethod
    def outlets():
        return ["out1"]

    def comp_init(self, nw):

        if not self.X.is_set:
            msg = "The conversion yield (X) of must be set!"
            logging.error(msg)

        if not self.formula.is_set:
            msg = "The reaction formula of must be set!"
            logging.error(msg)

        # parse reaction equation
        self.reac = Reaction(self.formula.val)

        # calculate outlet composition
        m = self.inl[0].m.val_SI  # kg/s
        molar_masses = {
            k: CP.PropsSI("M", k)
            for k in self.reac.get_stoich_coefficients().keys()
        }  # kg/mol
        self.n_in = {
            k: m * x / molar_masses[k]
            for k, x in self.inl[0].fluid.val.items()
        }  # mol/s
        n_out = self.reac.calculate_composition(self.n_in, self.X.val)  # mol/s
        m_out = {k: n * molar_masses[k] for k, n in n_out.items()}  # kg/s
        self.x_out = {
            k: m_i / sum(m_out.values()) for k, m_i in m_out.items()
        }  # mass fractions

        # calculate reaction enthalpy
        T_in = self.inl[0].T.val_SI
        self.reaction_enthalpy_SI = self.reac.calculate_reaction_enthalpy(
            n0=self.n_in, X=self.X.val, T0=T_in
        )  # J/s

        Component.comp_init(self, nw)

    def energy_balance_func(self):
        r"""
        Calculate the residual in energy balance.

        Returns
        -------
        residual : float
            Residual value of energy balance equation.

            .. math::

                \begin{split}
                0=&P + \dot{m}_\mathrm{in,2}\cdot\left(h_\mathrm{in,2}-
                h_\mathrm{in,2,ref}\right)\\
                &-\dot{m}_\mathrm{in,1}\cdot\left( h_\mathrm{out,1} -
                h_\mathrm{in,1} \right)\\
                & -\dot{m}_\mathrm{out,2} \cdot \left( h_\mathrm{out,2} -
                h_\mathrm{out,2,ref} \right)\\
                & +\dot{m}_\mathrm{out,3} \cdot \left( h_\mathrm{out,3} -
                h_\mathrm{out,3,ref} + e_0\right)\\
                \end{split}

            - Reference temperature: 298.15 K.
            - Reference pressure: 1 bar.
        """
        return (
            self.outl[0].m.val_SI * self.outl[0].h.val_SI
            - self.inl[0].m.val_SI * self.inl[0].h.val_SI
        )

    def energy_balance_func_doc(self, label):
        r"""
        Calculate the residual in energy balance.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = (
            r"\begin{split}" + "\n"
            r"0=&P + \dot{m}_\mathrm{in,2}\cdot\left(h_\mathrm{in,2}-"
            r"h_\mathrm{in,2,ref}\right)\\" + "\n"
            r"&-\dot{m}_\mathrm{in,1}\cdot\left( h_\mathrm{out,1} -"
            r"h_\mathrm{in,1} \right)\\" + "\n"
            r"& - \dot{m}_\mathrm{out,2} \cdot \left( h_\mathrm{out,2} -"
            r"h_\mathrm{out,2,ref} \right)\\" + "\n"
            r"& + \dot{m}_\mathrm{out,3} \cdot \left( h_\mathrm{out,3} -"
            r"h_\mathrm{out,3,ref} + e_0\right)\\" + "\n"
            r"&p_\mathrm{ref}=\unit[1]{bar},"
            r"\;T_\mathrm{ref}=\unit[25]{^\circ C}\\" + "\n"
            r"\end{split}"
        )
        return generate_latex_eq(self, latex, label)

    def energy_balance_deriv(self, increment_filter, k):
        r"""
        Partial derivatives for reactor energy balance.

        Parameters
        ----------
        increment_filter : ndarray
            Matrix for filtering non-changing variables.

        k : int
            Position of derivatives in Jacobian matrix (k-th equation).
        """
        # inlet
        self.jacobian[k, 0, 0] = -self.inl[0].h.val_SI
        self.jacobian[k, 0, 2] = -self.inl[0].m.val_SI

        # outlet
        self.jacobian[k, 1, 0] = self.outl[0].h.val_SI
        self.jacobian[k, 1, 2] = self.outl[0].m.val_SI

    def fluid_func(self):
        r"""
        Equations for fluid composition.

        Returns
        -------
        residual : list
            Residual values of equation.

            .. math::

                0  = x_\mathrm{i,in,1} - x_\mathrm{i,out,1}
                \forall i \in \text{network fluids}

                0 = \begin{cases}
                    1 - x_\mathrm{i,in2} & \text{i=}H_{2}O\\
                    x_\mathrm{i,in2} & \text{else}
                \end{cases} \forall i \in \text{network fluids}

                0 = \begin{cases}
                    1 - x_\mathrm{i,out,2} & \text{i=}O_{2}\\
                    x_\mathrm{i,out,2} & \text{else}
                \end{cases} \forall i \in \text{network fluids}

                0 = \begin{cases}
                    1 - x_\mathrm{i,out,3} & \text{i=}H_{2}\\
                    x_\mathrm{i,out,3} & \text{else}
                \end{cases} \forall i \in \text{network fluids}
        """
        expected_composition = self.x_out
        residual = [
            x - expected_composition.get(fluid, 0)
            for fluid, x in self.outl[0].fluid.val.items()
        ]
        return residual

    def fluid_func_doc(self, label):
        r"""
        Equations for fluid composition.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = (
            r"\begin{split}" + "\n"
            r"0 = &x_\mathrm{i,in,1} - x_\mathrm{i,out,1}\\" + "\n"
            r"0 = &\begin{cases}" + "\n"
            r"1 - x_\mathrm{i,in2} & \text{i=}H_{2}O\\" + "\n"
            r"x_\mathrm{i,in2} & \text{else}\\" + "\n"
            r"\end{cases}\\" + "\n"
            r"0 =&\begin{cases}" + "\n"
            r"1 - x_\mathrm{i,out,2} & \text{i=}O_{2}\\" + "\n"
            r"x_\mathrm{i,out,2} & \text{else}\\" + "\n"
            r"\end{cases}\\" + "\n"
            r"0 =&\begin{cases}" + "\n"
            r"1 - x_\mathrm{i,out,3} & \text{i=}H_{2}\\" + "\n"
            r"x_\mathrm{i,out,3} & \text{else}\\" + "\n"
            r"\end{cases}\\" + "\n"
            r"&\forall i \in \text{network fluids}" + "\n"
            r"\end{split}"
        )
        return generate_latex_eq(self, latex, label)

    def fluid_deriv(self):
        r"""
        Calculate the partial derivatives for cooling loop fluid balance.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the fluid equations.
        """

        deriv = np.zeros(
            (
                len(self.inl[0].fluid.val),
                len(self.conn_loc) + self.num_vars,
                self.num_nw_vars,
            )
        )

        out_id = self.outl[0].conn_loc
        for i in range(len(self.inl[0].fluid.val())):
            deriv[i, out_id, 3 + i] = 1

        return deriv

    def mass_flow_func(self):
        r"""
        Equations for mass conservation.

        Returns
        -------
        residual : list
            Residual values of equation.

            .. math::

                0 =\dot{m}_\mathrm{in,1}-\dot{m}_\mathrm{out,1}\\
        """
        return self.inl[0].m.val_SI - self.outl[0].m.val_SI

    def mass_flow_func_doc(self, label):
        r"""
        Equations for mass conservation.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = r"0 =\dot{m}_\mathrm{in,1}-\dot{m}_\mathrm{out,1}"
        return generate_latex_eq(self, latex, label)

    def mass_flow_deriv(self):
        r"""
        Calculate the partial derivatives for all mass flow balance equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the mass flow equations.
        """
        # deritatives for mass flow balance in the heat exchanger
        deriv = np.zeros(
            (1, len(self.conn_loc) + self.num_vars, self.num_nw_vars)
        )
        deriv[0, 0, 0] = 1
        deriv[0, 1, 0] = -1

        return deriv

    def reactor_pressure_func(self):
        r"""
        Equations for reactor pressure balance.

        Returns
        -------
        residual : list
            Residual values of equation.

            .. math::

                0 = p_\mathrm{in,1} - p_\mathrm{out,1}\\
        """
        return self.inl[0].p.val_SI - self.outl[0].p.val_SI

    def reactor_pressure_func_doc(self, label):
        r"""
        Equations for reactor pressure balance.

        Parameters
        ----------
        label : str
            Label for equation.

        Returns
        -------
        latex : str
            LaTeX code of equations applied.
        """
        latex = r"0 = p_\mathrm{in,1} - p_\mathrm{out,1}"
        return generate_latex_eq(self, latex, label)

    def reactor_pressure_deriv(self):
        r"""
        Calculate the partial derivatives for combustion pressure equations.

        Returns
        -------
        deriv : ndarray
            Matrix with partial derivatives for the pressure equations.
        """
        deriv = np.zeros(
            (1, len(self.conn_loc) + self.num_vars, self.num_nw_vars)
        )
        deriv[0, 0, 1] = 1
        deriv[0, 1, 1] = -1

        return deriv

    def bus_func(self, bus):
        r"""
        Calculate the value of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        val : float
            Value of energy transfer :math:`\dot{E}`. This value is passed to
            :py:meth:`tespy.components.component.Component.calc_bus_value`
            for value manipulation according to the specified characteristic
            line of the bus.

            .. math::

                \dot{E} = \begin{cases}
                P & \text{key = 'P'}\\
                - \dot{m}_{in,1} \cdot \left(h_{out,1} - h_{in,1} \right) &
                \text{key = 'Q'}\\
                \end{cases}
        """
        ######################################################################
        # equations for power on bus
        if bus["param"] == "P":
            val = 0

        ######################################################################
        # equations for heat on bus

        elif bus["param"] == "Q":
            val = -self.inl[0].m.val_SI * (
                self.outl[0].h.val_SI - self.inl[0].h.val_SI
            )

        ######################################################################
        # missing/invalid bus parameter

        else:
            msg = (
                "The parameter " + str(bus["param"]) + " is not a valid "
                "parameter for a component of type "
                + self.component()
                + ". Please specify a bus parameter (P/Q) for component "
                + self.label
                + "."
            )
            logging.error(msg)
            raise ValueError(msg)

        return val

    def bus_func_doc(self, bus):
        r"""
        Return LaTeX string of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        latex : str
            LaTeX string of bus function.
        """
        if bus["param"] == "P":
            return r"P_\mathrm{el}"
        elif bus["param"] == "Q":
            return (
                r"-\dot{m}_\mathrm{in,1} \cdot \left(h_\mathrm{out,1} - "
                r"h_\mathrm{in,1} \right)"
            )

    def bus_deriv(self, bus):
        r"""
        Calculate partial derivatives of the bus function.

        Parameters
        ----------
        bus : tespy.connections.bus.Bus
            TESPy bus object.

        Returns
        -------
        deriv : ndarray
            Matrix of partial derivatives.
        """
        deriv = np.zeros(
            (1, len(self.conn_loc) + self.num_vars, self.num_nw_vars)
        )
        f = self.calc_bus_value
        b = bus.comps.loc[self]

        ######################################################################
        # derivatives for power on bus
        if b["param"] == "P":
            deriv[0, 0, 0] = self.numeric_deriv(f, "m", 0, bus=bus)
            deriv[0, 0, 2] = self.numeric_deriv(f, "h", 0, bus=bus)

            deriv[0, 1, 0] = self.numeric_deriv(f, "m", 1, bus=bus)
            deriv[0, 1, 2] = self.numeric_deriv(f, "h", 1, bus=bus)

        ######################################################################
        # derivatives for heat on bus
        elif b["param"] == "Q":

            deriv[0, 0, 0] = self.numeric_deriv(f, "m", 0, bus=bus)
            deriv[0, 0, 2] = self.numeric_deriv(f, "h", 0, bus=bus)
            deriv[0, 1, 2] = self.numeric_deriv(f, "h", 1, bus=bus)

        ######################################################################
        # missing/invalid bus parameter

        else:
            msg = (
                "The parameter " + str(b["param"]) + " is not a valid "
                "parameter for a component of type "
                + self.component()
                + ". Please specify a bus parameter (P/Q) for component "
                + self.label
                + "."
            )
            logging.error(msg)
            raise ValueError(msg)

        return deriv

    def initialise_fluids(self):
        for c in self.outl:
            c.target.propagate_fluid_to_target(c, c.target)
        for c in self.inl:
            c.source.propagate_fluid_to_source(c, c.source)

    def initialise_source(self, c, key):
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
                5  \cdot 10^5 & \text{key = 'p'}\\
                h\left(T=323.15, p=5  \cdot 10^5\right) & \text{key = 'h'}
                \end{cases}
        """
        if key == "p":
            return 5e5
        elif key == "h":
            flow = c.get_flow()
            T = 50 + 273.15
            return h_mix_pT(flow, T)

    def initialise_target(self, c, key):
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
                5  \cdot 10^5 & \text{key = 'p'}\\
                h\left(T=293.15, p=5  \cdot 10^5\right) & \text{key = 'h'}
                \end{cases}
        """
        if key == "p":
            return 5e5
        elif key == "h":
            flow = c.get_flow()
            T = 20 + 273.15
            return h_mix_pT(flow, T)

    def propagate_fluid_to_target(self, inconn, start):
        r"""
        Propagate the fluids towards connection's target in recursion.

        Parameters
        ----------
        inconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        if inconn == self.inl[0]:
            outconn = self.outl[0]

            for fluid, x in inconn.fluid.val.items():
                if (
                    outconn.fluid.val_set[fluid] is False
                    and outconn.good_starting_values is False
                ):
                    outconn.fluid.val[fluid] = x

            outconn.target.propagate_fluid_to_target(outconn, start)

    def propagate_fluid_to_source(self, outconn, start):
        r"""
        Propagate the fluids towards connection's source in recursion.

        Parameters
        ----------
        outconn : tespy.connections.connection.Connection
            Connection to initialise.

        start : tespy.components.component.Component
            This component is the fluid propagation starting point.
            The starting component is saved to prevent infinite looping.
        """
        if outconn == self.outl[0]:
            inconn = self.inl[0]

            for fluid, x in outconn.fluid.val.items():
                if (
                    inconn.fluid.val_set[fluid] is False
                    and inconn.good_starting_values is False
                ):
                    inconn.fluid.val[fluid] = x

            inconn.source.propagate_fluid_to_source(inconn, start)

    def calc_parameters(self):
        r"""Postprocessing parameter calculation."""
        pass
