# -*- coding: utf-8

"""Module for thermodynamic analyses.

The analyses module provides thermodynamic analysis tools for your simulation.
Different analysis classes are available:

- :py:class:`tespy.tools.analyses.ExergyAnalysis`


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/tools/analyses.py

SPDX-License-Identifier: MIT
"""
import numpy as np
import pandas as pd
from tabulate import tabulate

from tespy.tools import helpers as hlp
from tespy.tools import logger
from tespy.tools.fluid_properties import single_fluid
from tespy.tools.global_vars import ERR
from tespy.tools.global_vars import combustion_gases

idx = pd.IndexSlice


def categorize_fluids(conn):
    fluid = single_fluid(conn.fluid_data)
    if fluid is None:
        cat = "non-combustion-gas"
        for f, x in conn.fluid.val.items():
            if x > ERR :
                try:
                    if hlp.fluidalias_in_list(f, combustion_gases):
                        cat = "combustion-gas"
                        break
                except RuntimeError:
                    # CoolProp cannot call aliases on incompressibles
                    pass
    else:
        is_incompressible = False
        try:
            is_combustion_gas = hlp.fluidalias_in_list(fluid, combustion_gases)
        except RuntimeError:
            # CoolProp cannot call aliases on incompressibles
            is_incompressible = True
            is_combustion_gas = False
        if is_combustion_gas:
            cat = "combustion-gas"
        elif is_incompressible:
            cat = "incompressible"
        else:
            cat = "two-phase-fluid"
    return cat


class ExergyAnalysis:
    r"""Class for exergy analysis of TESPy models."""
    exergy_cats = ["chemical", "physical", "massless"]

    def __init__(self, network, E_F, E_P, E_L=[], internal_busses=[]):
        r"""
        Parameters
        ----------
        E_F : list
            List containing busses which represent fuel exergy input of the
            network, e.g. heat exchangers of the steam generator.

        E_P : list
            List containing busses which represent exergy production of the
            network, e.g. the motors and generators of a power plant.

        E_L : list
            Optional: List containing busses which represent exergy loss
            streams of the network to the ambient, e.g. flue gases of a gas
            turbine.

        internal_busses : list
            Optional: List containing internal busses that represent exergy
            transfer within your network but neither exergy production or
            exergy fuel, e.g. a steam turbine driven feed water pump. The
            conversion factors of the bus are applied to calculate exergy
            destruction which is allocated to the respective components.

        Note
        ----
        The nomenclature of the variables used in the exergy analysis is
        according to :cite:`Tsatsaronis2007`.

        The analysis is carried out by the
        :py:meth:`tespy.tools.analyses.ExergyAnalysis.analyse` method. Given
        the ambient state (pressure and temperature), it will

        - Calculate the values of physical exergy on all connections
        - Calculate exergy balance for all components. The individual exergy
          balance methods are documented in the API-documentation of the
          respective components.

          - Components for which no exergy balance has yet been implemented,
            :code:`nan` (not defined) is assigned for fuel and product
            exergy as well as exergy destruction and exergetic efficiency.
          - Dissipative components do not have product exergy (:code:`nan`) per
            definition.

        - Calculate exergy balances for busses passed to ExergyAnalysis class
          instance.
        - Calculate network fuel exergy, product exergy as well as exergy loss
          from data provided by the busses passed to the instance.
        - Calculate network exergetic efficiency.
        - Calculate exergy destruction ratios for components.

          - :math:`y_\mathrm{D}` compare the rate of exergy destruction in a
            component to the exergy rate of the fuel provided to the overall
            system.
          - :math:`y^*_\mathrm{D}` compare the component exergy destruction
            rate to the total exergy destruction rate within the system.

        .. math::

            \begin{split}
            \dot{E}_{\mathrm{D},comp} = \dot{E}_{\mathrm{F},comp} - \dot{E}_{\mathrm{P},comp}
            \;& \\
            \varepsilon_{\mathrm{comp}} =
            \frac{\dot{E}_{\mathrm{P},comp}}{\dot{E}_{\mathrm{F},comp}} \;& \\
            \dot{E}_{\mathrm{D}} = \sum_{comp} \dot{E}_{\mathrm{D},comp} \;&
            \forall comp \in \text{ network components}\\
            \dot{E}_{\mathrm{P}} = \sum_{comp} \dot{E}_{\mathrm{P},comp} \;&
            \forall comp \in
            \text{ components of busses in E\_P if 'base': 'component'}\\
            \dot{E}_{\mathrm{P}} = \dot{E}_{\mathrm{P}} - \sum_{comp} \dot{E}_{\mathrm{F},comp}
            \;& \forall comp \in
            \text{ components of busses in E\_P if 'base': 'bus'}\\
            \dot{E}_{\mathrm{F}} = \sum_{comp} \dot{E}_{\mathrm{F},comp} \;&
            \forall comp \in
            \text{ components of busses in E\_F if 'base': 'bus'}\\
            \dot{E}_{\mathrm{F}} = \dot{E}_{\mathrm{F}} - \sum_{comp} \dot{E}_{\mathrm{P},comp}
            \;& \forall comp \in
            \text{ components of busses in E\_F if 'base': 'component'}\\
            \dot{E}_{\mathrm{L}} = \sum_{comp} \dot{E}_{\mathrm{D},comp} \;&
            \forall comp \in
            \text{ sinks of network components if parameter exergy='loss'}
            \end{split}

        The exergy balance of the network must be closed, meaning fuel exergy
        minus product exergy, exergy destruction and exergy losses must be
        zero (:math:`\Delta \dot{E}_\text{max}=0.001`). If the balance is violated a
        warning message is prompted.

        .. math::

            |\dot{E}_{\text{F}} - \dot{E}_{\text{P}} - \dot{E}_{\text{L}} - \dot{E}_{\text{D}}| \leq
            \Delta \dot{E}_\text{max}\\

            \varepsilon = \frac{\dot{E}_{\text{P}}}{\dot{E}_{\text{F}}}

            y_{\text{D},comp} =
            \frac{\dot{E}_{\text{D},comp}}{\dot{E}_{\text{F}}}\\
            y^*_{\text{D},comp} =
            \frac{\dot{E}_{\text{D},comp}}{\dot{E}_{\text{D}}}

        Example
        -------
        In this example a simple clausius rankine cycle is set up and an
        exergy analysis is performed after simulation of the power plant.
        Start by defining ambient state and genereral network setup.

        >>> from tespy.components import (CycleCloser, SimpleHeatExchanger,
        ... Merge, Splitter, Valve, Compressor, Pump, Turbine)
        >>> from tespy.connections import Bus
        >>> from tespy.connections import Connection
        >>> from tespy.networks import Network
        >>> from tespy.tools import ExergyAnalysis

        >>> Tamb = 20
        >>> pamb = 1
        >>> nw = Network()
        >>> nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg',
        ... iterinfo=False)

        In order to show all functionalities available we use a feed water pump
        that is not driven electrically by a motor but instead internally by
        an own steam turbine. Therefore we split up the live steam from the
        steam generator and merge the streams after both steam turbines. For
        simplicity the steam generator and the condenser are modeled as simple
        heat exchangers.

        >>> cycle_close = CycleCloser('cycle closer')
        >>> splitter1 = Splitter('splitter 1')
        >>> merge1 = Merge('merge 1')
        >>> turb = Turbine('turbine')
        >>> fwp_turb = Turbine('feed water pump turbine')
        >>> condenser = SimpleHeatExchanger('condenser')
        >>> fwp = Pump('pump')
        >>> steam_generator = SimpleHeatExchanger('steam generator')

        >>> fs_in = Connection(cycle_close, 'out1', splitter1, 'in1')
        >>> fs_fwpt = Connection(splitter1, 'out1', fwp_turb, 'in1')
        >>> fs_t = Connection(splitter1, 'out2', turb, 'in1')
        >>> fwpt_ws = Connection(fwp_turb, 'out1', merge1, 'in1')
        >>> t_ws = Connection(turb, 'out1', merge1, 'in2')
        >>> ws = Connection(merge1, 'out1', condenser, 'in1')
        >>> cond = Connection(condenser, 'out1', fwp, 'in1')
        >>> fw = Connection(fwp, 'out1', steam_generator, 'in1')
        >>> fs_out = Connection(steam_generator, 'out1', cycle_close, 'in1')
        >>> nw.add_conns(fs_in, fs_fwpt, fs_t, fwpt_ws, t_ws, ws, cond,
        ... fw, fs_out)

        Next step is to set up the busses to later pass them according to the
        convetions in the list below:

        - E_F for fuel exergy
        - E_P for product exergy
        - internal_busses for internal energy transport
        - E_L for exergy loss streams to the ambient (sources and sinks go
          here, in case you use e.g. flue gases or air input)

        The first bus is for output power, which is only represented by the
        main steam turbine. The efficiency is set to 0.97. This bus will
        represent the product exergy.

        >>> power = Bus('power_output')
        >>> power.add_comps({'comp': turb, 'char': 0.97})

        The second bus is for driving the feed water pump. The total power of
        this bus is specified to be 0 in order to make sure, the power genrated
        by the secondary steam turbine is transferred to the feed water pump.
        For mechanical efficiency we choose 0.985 for both components, but
        we need to make sure, the :code:`'base'` of the feed water pump is
        :code:`'bus'` as the energy from the turbine drives the feed water
        pump.

        >>> fwp_power = Bus('feed water pump power', P=0)
        >>> fwp_power.add_comps(
        ... {'comp': fwp_turb, 'char': 0.985},
        ... {'comp': fwp, 'char': 0.985, 'base': 'bus'})

        The fuel exergy is the exergy input into the network which is
        represented by the heat input bus. Here again, as we have an energy
        input from outside of the network, the :code:`'base'` keyword must be
        specified to :code:`'bus'`.

        >>> heat = Bus('heat_input')
        >>> heat.add_comps({'comp': steam_generator, 'base': 'bus'})
        >>> nw.add_busses(power, fwp_power, heat)

        After setting up the busses, we specify the parameters for components
        and connections and start the simulation.

        >>> turb.set_attr(eta_s=0.9)
        >>> fwp_turb.set_attr(eta_s=0.87)
        >>> condenser.set_attr(pr=0.98)
        >>> fwp.set_attr(eta_s=0.75)
        >>> steam_generator.set_attr(pr=0.89)
        >>> fs_in.set_attr(m=10, p=120, T=600, fluid={'water': 1})
        >>> cond.set_attr(T=Tamb + 3, x=0)
        >>> nw.solve('design')

        To evaluate the exergy balance of the network, we create an instance of
        class :py:class:`tespy.tools.analyses.ExergyAnalysis`
        passing the network to analyse as well as the respective busses. To run
        the analysis, the ambient state is passed subsequently. The results of
        the analysis can be printed using the
        :py:meth:`tespy.tools.analyses.ExergyAnalysis.print_results` method.
        The exergy balance should be closed, if you set up your network
        analysis correctly. If not, an error is prompted.

        >>> ean = ExergyAnalysis(network=nw, E_F=[heat], E_P=[power],
        ... internal_busses=[fwp_power])
        >>> ean.analyse(pamb=pamb, Tamb=Tamb)
        >>> abs(round(ean.network_data['E_F'] - ean.network_data['E_P'] -
        ... ean.network_data['E_L'] - ean.network_data['E_D'], 3))
        0.0
        >>> ();ean.print_results();() # doctest: +ELLIPSIS
        (...)

        The exergy data of the passed busses, the network's components and
        connections as well as the network itself are stored as dataframes and
        therefore accessible for further investigation.

        >>> busses = ean.bus_data
        >>> components = ean.component_data
        >>> connections = ean.connection_data
        >>> network = ean.network_data

        Additionally, if you defined component groups for your components, the
        exergy data of these groups are accumulated and saved in an own
        DataFrame. Please note, that the exergy destruction values of the
        busses are allocated to the component groups in this case.

        >>> groups = ean.group_data

        Lastly, the tool offers an automatically generated data input for the
        sankey plotting functionalities of the plotly library to create a
        Grassmann diagram of your network. For more information on the usage
        of the functionality look up the corresponding method documentation:
        :py:meth:`tespy.tools.analyses.ExergyAnalysis.generate_plotly_sankey_input`.
        The method returns a dictionary containting the links for the sankey
        as well as a list of the nodes.

        >>> links, nodes = ean.generate_plotly_sankey_input()
        """
        if len(E_F) == 0:
            msg = ('Missing fuel exergy E_F of network.')
            logger.error(msg)
            raise hlp.TESPyNetworkError(msg)

        if len(E_P) == 0:
            msg = ('Missing product exergy E_P of network.')
            logger.error(msg)
            raise hlp.TESPyNetworkError(msg)

        self.nw = network
        self.E_F = E_F
        self.E_P = E_P
        self.E_L = E_L
        self.internal_busses = internal_busses

        bus_labels = [b.label for b in internal_busses + E_F + E_P + E_L]
        key_exergy_labels = ['E_P', 'E_F', 'E_D', 'E_L']
        self.reserved_fkt_groups = key_exergy_labels + bus_labels
        if len(set(bus_labels).intersection(key_exergy_labels)) > 0:
            msg = (
                "None of your busses may have the label '"
                + "', '".join(key_exergy_labels) + "' when performing the "
                + "exergy analysis."
            )
            raise ValueError(msg)

    def analyse(self, pamb, Tamb, Chem_Ex=None):
        """Run the exergy analysis.

        Parameters
        ----------
        pamb : float
            Ambient pressure value for analysis, provide value in network's
            pressure unit.

        Tamb : float
            Ambient temperature value for analysis, provide value in network's
            temperature unit.
        """
        pamb_SI = hlp.convert_to_SI('p', pamb, self.nw.p_unit)
        Tamb_SI = hlp.convert_to_SI('T', Tamb, self.nw.T_unit)

        # reset data
        dtypes = {
            "E_F": float,
            "E_P": float,
            "E_D": float,
            "epsilon": float,
            "group": str
        }
        self.component_data = pd.DataFrame(
            columns=list(dtypes.keys())
        ).astype(dtypes)

        self.bus_data = self.component_data.copy()
        self.bus_data["base"] = np.nan
        self.bus_data["base"] = self.bus_data["base"].astype(str)
        conn_exergy_data_cols = ['e_PH', 'e_T', 'e_M', 'E_PH', 'E_T', 'E_M']

        if Chem_Ex is not None:
            conn_exergy_data_cols += ['e_CH', 'E_CH']

        self.connection_data = pd.DataFrame(
            columns=conn_exergy_data_cols,
            dtype='float64'
        )

        self.network_data = pd.Series(
            index=['E_F', 'E_P', 'E_D', 'E_L'], dtype='float64'
        )
        self.network_data[:] = 0

        # physical exergy of connections
        for conn in self.nw.conns['object']:
            conn.get_physical_exergy(pamb_SI, Tamb_SI)
            conn.get_chemical_exergy(pamb_SI, Tamb_SI, Chem_Ex)
            conn_exergy_data = [
                conn.ex_physical, conn.ex_therm, conn.ex_mech,
                conn.Ex_physical, conn.Ex_therm, conn.Ex_mech
            ]
            if Chem_Ex is not None:
                conn_exergy_data += [conn.ex_chemical, conn.Ex_chemical]

            self.connection_data.loc[conn.label] = conn_exergy_data

        # todo: überprüfen der sankey data + massless exergy
        self.sankey_data = {}
        sankey_columns_dtypes = {
            'chemical': float,
            'physical': float,
            'massless': float
        }
        for label in self.reserved_fkt_groups:
            self.sankey_data[label] = pd.DataFrame(
                columns=sankey_columns_dtypes.keys(),
                index=pd.MultiIndex(
                    levels=[[], []],
                    names=["target_group", "category"],
                    codes=[[], []]
                )
            ).astype(sankey_columns_dtypes)

        # exergy balance of components
        for cp in self.nw.comps['object']:
            # save component information
            cp.exergy_balance(Tamb_SI)
            self.component_data.loc[cp.label] = [
                cp.E_F, cp.E_P, cp.E_D, cp.epsilon, cp.fkt_group
            ]

            if cp.fkt_group in self.reserved_fkt_groups:
                msg = (
                    'The labels ' + ', '.join(self.reserved_fkt_groups) + ' '
                    'cannot be used by components (if no group was assigned) '
                    'or component groups in the exergy analysis. Found '
                    'component/group with name ' + cp.fkt_group + '.'
                )
                raise ValueError(msg)
            elif cp.fkt_group not in self.sankey_data:
                self.sankey_data[cp.fkt_group] = pd.DataFrame(
                    columns=sankey_columns_dtypes.keys(),
                    index=pd.MultiIndex(
                        levels=[[], []],
                        names=["target_group", "category"],
                        codes=[[], []]
                    )
                ).astype(sankey_columns_dtypes)

            self.evaluate_busses(cp)

        # create a table that includes exergy destruction attributed to the
        # components
        bus_based = self.bus_data[self.bus_data['base'] == 'bus'].index
        component_based = self.bus_data[
            self.bus_data['base'] == 'component'
        ].index

        # add aggregated components with respective buses data
        self.aggregation_data = self.component_data.copy()
        # E_D is sum of both E_D
        self.aggregation_data.loc[self.bus_data.index, 'E_D'] = (
            self.component_data.loc[self.bus_data.index, 'E_D'] +
            self.bus_data['E_D']
        )
        # E_F for bus based components is higher by E_D of bus
        self.aggregation_data.loc[bus_based, 'E_F'] += (
            self.bus_data.loc[bus_based, 'E_D']
        )
        # E_P of component based components is lower by E_D of bus
        self.aggregation_data.loc[component_based, 'E_P'] -= (
            self.bus_data.loc[component_based, 'E_D']
        )

        # calculate network results
        self.network_data.loc['E_D'] = (
            self.component_data['E_D'].sum() + self.bus_data['E_D'].sum())
        self.network_data.loc['E_F'] = abs(self.network_data.loc['E_F'])
        self.network_data.loc['E_P'] = abs(self.network_data.loc['E_P'])
        self.network_data.loc['epsilon'] = (
            self.network_data.loc['E_P'] / self.network_data.loc['E_F']
        )

        # calculate exergy destruction ratios for components/busses
        E_F = self.network_data.loc['E_F']
        E_D = self.network_data.loc['E_D']

        for d in [self.component_data, self.bus_data, self.aggregation_data]:
            d['y_Dk'] = d['E_D'] / E_F
            d['y*_Dk'] = d['E_D'] / E_D
            d['epsilon'] = d['E_P'] / d['E_F']

        residual = abs(
            self.network_data.loc['E_F'] - self.network_data.loc['E_P'] -
            self.network_data.loc['E_D'] - self.network_data.loc['E_L']
        )

        if residual >= ERR ** 0.5:
            msg = (
                'The exergy balance of your network is not closed (residual '
                'value is ' + str(round(residual, 6)) + ', but should be '
                'smaller than ' + str(ERR ** 0.5) + '), you should check the '
                'component and network exergy data and check, if network is '
                'properly setup for the exergy analysis.')
            logger.error(msg)

        self.create_group_data()

    """+F+F+F+F++++START++++F+F+F+F+"""
    def evaluate_exergoeconomics(self, Exe_Eco_Costs, Tamb):
        Tamb_SI = hlp.convert_to_SI('T', Tamb, self.nw.T_unit)
        # TO DO: need to add checks and error messages

        # DETERMINE NUMBER OF VARIABLES FOR MATRIX AND THEIR POSITION
        col_number = 0
        for comp in self.nw.comps["object"]:
            comp.Ex_C_col = {}
            comp.C_bus =  np.nan

        for bus in self.E_F + self.E_P:
            bus.Ex_C_col = {}
            bus.E_external = {}
            bus.C_external = {}

            if not any([comp.component() in ["source", "sink"] for comp in bus.comps.index]):
                bus.Ex_C_col["external_bus"] = col_number
                col_number += 1
                print("bus: ", bus.label, bus.Ex_C_col)
                # calculate exergy transportet to/from external
                bus.E_external["bus"] = 0
                for comp, comp_data in bus.comps.iterrows():
                    if comp_data["base"] == "bus":
                        bus.E_external["bus"] += comp.E_bus["massless"] #/ comp.calc_bus_efficiency(bus)
                        #print(f"{bus.label} to {comp.label}")
                        #print(comp.label, "E_bus[massless] component based: ", comp.E_bus["massless"])
                        #print(comp.label, "comp.E_bus[massless] bus based: ", comp.E_bus["massless"] / comp.calc_bus_efficiency(bus))
                        # / efficiency -> is that correct?? is comp.E_bus the component value of the bus
                        # E_external > 0 means bus gets energy from external source
                    else:
                        bus.E_external["bus"] -= comp.E_bus["massless"] #* comp.calc_bus_efficiency(bus)
                        #print(f"{comp.label} to {bus.label}")
                        #print(comp.label, "E_bus[massless] component based: ", comp.E_bus["massless"])
                        #print(comp.label, "E_bus[massless] bus based: ", comp.E_bus["massless"] * comp.calc_bus_efficiency(bus))
                        # is comp.E_bus the component value of the bus
                        # E_external < 0 means bus gives energy to external sink
                # print("E_bus_external: ", bus.E_external["bus"], bus.label)
                # print("bus.P.val: ", bus.P.val, bus.label) # if those were the same, no need for E_external

                # WRITE GIVEN POWER AND HEAT COSTS INTO OBJECTS
                if f"{bus.label}_c" in Exe_Eco_Costs:
                    bus.C_external["bus"] = Exe_Eco_Costs[f"{bus.label}_c"] * abs(bus.E_external["bus"])
                    #print("bus external known costs: ", bus.C_external["bus"], bus.label)
                elif bus.E_external["bus"] > 0: # standard values
                    bus.C_external["bus"] = 0.02 * abs(bus.E_external["bus"])
                    print("no bus costs provided")
            for comp in bus.comps.index:                # 2 variables per component (bus based and component based)
                if comp.component() not in ["source", "sink"]:
                    comp.Ex_C_col[f"{bus.label}_bus"] = col_number
                    col_number += 1
                    comp.Ex_C_col[f"{bus.label}_component"] = col_number
                    col_number += 1
                    print("comp bus: ", comp.label, bus.label, comp.Ex_C_col, comp.E_bus["massless"])

        for comp in self.nw.comps['object']:
            if comp.dissipative.val:
                comp.Z_col = col_number      # ändern zu comp.Ex_C_col["_dissipative"], dann auch in valve, heat exchanger
                print("dissipative: ", comp.label, col_number)
                col_number += 1
        num_variables = (
            len(self.nw.conns) * 3
            + col_number
        )

        for conn_num, conn in enumerate(self.nw.conns["object"]):
            conn.Ex_C_col = {
                "therm": conn_num * 3 + col_number,
                "mech": conn_num * 3 + 1 + col_number,
                "chemical": conn_num * 3 + 2 + col_number
            }
            print(conn.label, conn.Ex_C_col)

        # WRITE GIVEN COMPONENT AND SOURCE COSTS INTO OBJECTS
        for bus in self.E_F + self.E_P + self.E_L + self.internal_busses:
            for i, comp in enumerate(bus.comps.index):
                if comp.component() == "source":
                    c_id = f"{comp.label}_c"                                    # material costs
                    if c_id in Exe_Eco_Costs:
                        comp.set_source_costs(Exe_Eco_Costs[f"{c_id}"])
                    else:
                        comp.set_source_costs_standard()
                        print("no source costs entered")

        for comp in self.nw.comps['object']:
            Z_id = f"{comp.label}_Z"                                            # component costs
            if Z_id in Exe_Eco_Costs:
                comp.set_Z_costs(Exe_Eco_Costs[f"{Z_id}"])
            else:
                if not comp.component() in["source", "sink"]:
                    comp.set_Z_costs_standard()

        # SETTING UP THE MATRIX
        exergy_cost_matrix = np.zeros([num_variables,num_variables])
        exergy_cost_vector = np.zeros(num_variables)
        counter = 0

        # ADD KNOWN BUS COSTS TO MATRIX
        for bus in self.E_F:
            # bus.Ex_C_col is the position of the variable going over the system boarder (external)
            if len(bus.Ex_C_col)>0 and len(bus.C_external)>0:
                exergy_cost_matrix[counter][bus.Ex_C_col["external_bus"]] = 1
                exergy_cost_vector[counter] = bus.C_external["bus"]
                counter += 1

        # ADD AUXILIARY EQUATION FOR BUS
        for bus in self.E_F + self.E_P + self.internal_busses:
            # all outputs have same c
            if not any([comp.component() in ["source", "sink"] for comp in bus.comps.index]):
                supplied_components = []
                for comp, comp_data in bus.comps.iterrows():
                    if comp_data["base"] == "bus":      # goes out of bus
                        supplied_components.append(comp)
                # if connected to external
                if len(bus.Ex_C_col)>0:
                    if bus.E_external["bus"] < 0:       # bus gives exergy to external
                        #print("bus gives exergy to external ", bus.label)
                        for comp in supplied_components:
                            exergy_cost_matrix[counter][bus.Ex_C_col["external_bus"]] = 1 / abs(bus.E_external["bus"])
                            exergy_cost_matrix[counter][comp.Ex_C_col[f"{bus.label}_bus"]] = -1 / (comp.E_bus["massless"] ) #/ comp.calc_bus_efficiency(bus)
                            # is comp.E_bus the component value of the bus
                            exergy_cost_vector[counter] = 0
                            counter += 1
                    elif bus.E_external["bus"] > 0:    # bus needs external power
                        #print("bus needs external power ", bus.label)
                        for i in range(len(supplied_components) -1):
                            exergy_cost_matrix[counter][supplied_components[i].Ex_C_col[f"{bus.label}_bus"]] = 1 / supplied_components[i].E_bus["massless"]
                            exergy_cost_matrix[counter][supplied_components[i+1].Ex_C_col[f"{bus.label}_bus"]] = -1 / (comp.E_bus["massless"] )     #/ comp.calc_bus_efficiency(bus)
                            # is comp.E_bus the component value of the bus
                            exergy_cost_vector[counter] = 0
                            counter += 1
                # if not connected to external --> same as if needs external power
                else:
                    #print("bus not connected to external ", bus.label)
                    for i in range(len(supplied_components) -1):
                        exergy_cost_matrix[counter][supplied_components[i].Ex_C_col[f"{bus.label}_bus"]] = 1 / (supplied_components[i].E_bus["massless"] )  #/ comp.calc_bus_efficiency(bus)
                        # is comp.E_bus the component value of the bus
                        exergy_cost_matrix[counter][supplied_components[i+1].Ex_C_col[f"{bus.label}_bus"]] = -1 / (supplied_components[i+1].E_bus["massless"] ) #/ comp.calc_bus_efficiency(bus)
                        # is comp.E_bus the component value of the bus
                        exergy_cost_vector[counter] = 0
                        counter += 1

            # for each component C_bus = C_component (also for internal busses) as long as there are no bus costs
                for comp in bus.comps.index:
                    exergy_cost_matrix[counter][comp.Ex_C_col[f"{bus.label}_bus"]] = 1
                    exergy_cost_matrix[counter][comp.Ex_C_col[f"{bus.label}_component"]] = -1
                    exergy_cost_vector[counter] = 0
                    counter += 1


        # BUS BALANCE
        for bus in self.E_F + self.E_P + self.internal_busses:
            if not any([comp.component() in ["source", "sink"] for comp in bus.comps.index]):
                for comp, comp_data in bus.comps.iterrows():
                    if comp_data["base"] == "bus":              # seen from position of bus
                        exergy_cost_matrix[counter][comp.Ex_C_col[f"{bus.label}_bus"]] = -1
                        # goes out of the bus and into the component
                    else:
                        exergy_cost_matrix[counter][comp.Ex_C_col[f"{bus.label}_bus"]] = +1
                        # goes out of the component and into the bus

                if len(bus.Ex_C_col) > 0:                    # bus is connected to external source/sink
                    if bus.E_external["bus"] < 0:      # bus gives exergy to external     same as    if bus in self.E_P:
                        exergy_cost_matrix[counter][bus.Ex_C_col["external_bus"]] = -1
                    elif bus.E_external["bus"] > 0:       # bus needs exergy from external same as    elif bus in self.E_F:
                        exergy_cost_matrix[counter][bus.Ex_C_col["external_bus"]] = +1

                counter += 1

        # productive components
        for comp in self.nw.comps['object']:
            comp.exergy_cost_line = counter
            if not comp.dissipative.val:
                #print(comp.component())

                # ADD KNOWN SOURCE COSTS TO MATRIX
                if comp.component() == "source":
                    # thermal exergy related source costs
                    exergy_cost_matrix[counter][comp.outl[0].Ex_C_col["therm"]] = +1
                    exergy_cost_vector[counter] = comp.outl[0].C_therm
                    counter +=1

                    # mechanical exergy related source costs
                    exergy_cost_matrix[counter][comp.outl[0].Ex_C_col["mech"]] = +1
                    exergy_cost_vector[counter] = comp.outl[0].C_mech
                    counter +=1

                    # chemical exergy related source costs
                    exergy_cost_matrix[counter][comp.outl[0].Ex_C_col["chemical"]] = +1
                    exergy_cost_vector[counter] = comp.outl[0].C_chemical
                    counter +=1

                # ADD COST BALANCE FOR EACH COMPONENT that isn't a source or sink
                # and has an associated cost value or is a cycle closer (-> cost value = 0)
                # TO DO: error message if Z_cost is nan
                # make a check for all component Z cost simultaneously
                if comp.component() == "cycle closer":
                    exergy_cost_matrix[counter][comp.inl[0].Ex_C_col["therm"]] = +1
                    exergy_cost_matrix[counter][comp.outl[0].Ex_C_col["therm"]] = -1
                    exergy_cost_vector[counter] = 0
                    counter +=1
                    exergy_cost_matrix[counter][comp.inl[0].Ex_C_col["mech"]] = +1
                    exergy_cost_matrix[counter][comp.outl[0].Ex_C_col["mech"]] = -1
                    exergy_cost_vector[counter] = 0
                    counter +=1
                    exergy_cost_matrix[counter][comp.inl[0].Ex_C_col["chemical"]] = +1
                    # bei Kreisprozessen muss ein Wert vorgegeben werden damit das System lösbar ist
                    #exergy_cost_matrix[counter][comp.outl[0].Ex_C_col["chemical"]] = -1
                    exergy_cost_vector[counter] = 0
                    counter +=1
                if not (comp.component() in ["source", "sink", "cycle closer"] or np.isnan(comp.Z_costs)):
                    for comp_inl in comp.inl:
                        exergy_cost_matrix[counter][comp_inl.Ex_C_col["therm"]] = +1         # thermal exergy related costs
                        exergy_cost_matrix[counter][comp_inl.Ex_C_col["mech"]] = +1          # mechanical exergy related costs
                        exergy_cost_matrix[counter][comp_inl.Ex_C_col["chemical"]] = +1      # chemical exergy related costs
                    for comp_outl in comp.outl:
                        exergy_cost_matrix[counter][comp_outl.Ex_C_col["therm"]] = -1        # thermal exergy related costs
                        exergy_cost_matrix[counter][comp_outl.Ex_C_col["mech"]] = -1         # mechanical exergy related costs
                        exergy_cost_matrix[counter][comp_outl.Ex_C_col["chemical"]] = -1     # chemical exergy related costs

                    for bus in self.E_F + self.E_P + self.E_L + self.internal_busses:
                        if comp in bus.comps.index:             # if this bus is connected to this component
                            if bus.comps.loc[comp, 'base'] == 'bus':
                                exergy_cost_matrix[counter][comp.Ex_C_col[f"{bus.label}_component"]] = +1   # goes into component
                                #print("goes in")
                            else:
                                exergy_cost_matrix[counter][comp.Ex_C_col[f"{bus.label}_component"]] = -1   # goes out of component
                                #print("goes out")

                    exergy_cost_vector[counter] = -comp.Z_costs
                    counter +=1

                # ADD AUXILIARY EQUATIONS (F AND P RULES)
                aux_eqs = comp.aux_eqs(exergy_cost_matrix, exergy_cost_vector, counter, Tamb_SI)
                # print(aux_eqs)
                exergy_cost_matrix = aux_eqs[0]
                exergy_cost_vector = aux_eqs[1]
                counter = aux_eqs[2]

        # dissipative components (need to be handled last)
        for comp in self.nw.comps['object']:
            if comp.dissipative.val:
                print(comp.component(), "(dissipative)")
                print("serving: ", comp.serving_components)
                if comp.serving_components is None:
                    # assign dissipative costs to the whole system
                    numComps = sum(1 for comp in self.nw.comps['object'] if comp.component() not in ["source", "sink"])
                    if  numComps ==1:
                        print("Network consists of dissipative component only")
                    comp.serving_components = []
                    for c in self.nw.comps['object']:
                        if not c == comp and not c.component() in ["source", "sink", "cycle closer"] and not c.dissipative.val:
                            comp.serving_components.append(c)
                dis_eqs = comp.dissipative_balance(exergy_cost_matrix, exergy_cost_vector, counter, Tamb_SI)     # changes Z_costs of serving component
                # print(dis_eqs)
                exergy_cost_matrix = dis_eqs[0]
                exergy_cost_vector = dis_eqs[1]
                counter = dis_eqs[2]

        if counter!=num_variables:
            print("not enough equations: ", num_variables-counter, " missing")

        print(np.round(exergy_cost_matrix,3))
        print(np.round(exergy_cost_vector,3))

        try:
            C_sol = np.linalg.solve (exergy_cost_matrix, exergy_cost_vector)
            # C_sol = np.linalg.lstsq(A, b)[0]        # not solve, bc solve only works for well-determined systems, lstsq also for over-determined
            # carefull: works also for under-determined systems: need to add exception if this happens
        except np.linalg.LinAlgError:
            msg = f"System is not solvable\nA = \n{np.round(exergy_cost_matrix,3)} \n b =\n{np.round(exergy_cost_vector,3)}"
            logger.error(msg)
            raise hlp.TESPyNetworkError(msg)
        if np.isnan(exergy_cost_matrix).any():
            msg = f"System contains nan value\nA = \n{np.round(exergy_cost_matrix,3)} \n b =\n{np.round(exergy_cost_vector,3)}"
            logger.error(msg)
            raise hlp.TESPyNetworkError(msg)
        #print(np.round(C_sol,3))

        for conn in self.nw.conns["object"]:
            conn.C_therm = C_sol[conn.Ex_C_col["therm"]]
            conn.C_mech = C_sol[conn.Ex_C_col["mech"]]
            conn.C_physical = conn.C_therm + conn.C_mech
            conn.C_chemical = C_sol[conn.Ex_C_col["chemical"]]
            conn.C_tot = conn.C_physical + conn.C_chemical

            conn.c_therm = conn.C_therm / conn.Ex_therm if conn.Ex_therm != 0 else 0
            conn.c_mech = conn.C_mech / conn.Ex_mech if conn.Ex_mech != 0 else 0
            conn.c_physical = conn.C_physical / conn.Ex_physical if conn.Ex_physical != 0 else 0
            conn.c_chemical = conn.C_chemical / conn.Ex_chemical if conn.Ex_chemical != 0 else 0
            conn.c_tot = conn.C_tot / (conn.Ex_physical + conn.Ex_chemical) if (conn.Ex_physical + conn.Ex_chemical) != 0 else 0

        for bus in self.E_F + self.E_P:
            for comp, comp_data in bus.comps.iterrows():
                if not comp.component() in["source", "sink"]:
                    comp.C_bus = C_sol[comp.Ex_C_col[f"{bus.label}_bus"]]
                    comp.c_bus["bus"] = comp.C_bus / comp.E_bus["massless"] if comp.E_bus["massless"] != 0 else 0
                    if comp_data["base"] == "bus":
                        comp.c_bus["component"] = comp.c_bus["bus"] / comp.calc_bus_efficiency(bus)
                    else:
                        comp.c_bus["component"] = comp.c_bus["bus"] * comp.calc_bus_efficiency(bus)

                #print(comp.label, "c_bus bus based: ", comp.c_bus["bus"])
                #print(comp.label, "c_bus component based: ", comp.c_bus["component"])

        # exergoeconomic balance (C_F, C_P, C_L) of components
        for cp in self.nw.comps['object']:
            cp.exergoeconomic_balance(Tamb_SI)
            if not cp.component() in["source", "sink", "cycle closer"]:
                cp.calc_C_L()

        # exergoeconomic analysis of network
        self.network_data_exergo = pd.Series(
            index=['C_F', 'C_P', 'C_D', 'C_L'], dtype='float64'
        )
        self.network_data_exergo[:] = 0
        """
        # evaluate busses exergoeconomically
        for b in self.E_F + self.E_P + self.internal_busses + self.E_L:
            if cp in b.comps.index:
                if b.comps.loc[cp, 'base'] == 'bus':
                    C_bus = cp.C_bus
                    bus_efficiency = cp.calc_bus_efficiency(b)
                    C_F = C_bus / bus_efficiency

                    if b in self.E_F:
                        self.network_data_exergo.loc['C_F'] += C_F
                    elif b in self.E_P:
                        self.network_data_exergo.loc['C_P'] -= C_F
                    elif b in self.E_L:
                        self.network_data_exergo.loc['C_L'] -= C_F

                else:
                    C_bus = cp.C_bus
                    bus_efficiency = cp.calc_bus_efficiency(b)
                    C_P = C_bus * bus_efficiency

                    if b in self.E_F:
                        self.network_data_exergo.loc['C_F'] -= C_P
                    elif b in self.E_P:
                        self.network_data_exergo.loc['C_P'] += C_P
                    elif b in self.E_L:
                        self.network_data_exergo.loc['C_L'] += C_P

        self.network_data_exergo.loc['C_F'] = abs(self.network_data_exergo.loc['C_F'])
        self.network_data_exergo.loc['C_P'] = abs(self.network_data_exergo.loc['C_P'])
        """
        """+F+F+F+F++++END++++F+F+F+F+"""

    def evaluate_busses(self, cp):
        """Evaluate the exergy balances of busses.

        Parameters
        ----------
        cp : tespy.components.component.Component
            Component to analyse the bus exergy balance of.
        """
        cp_on_num_busses = 0
        for b in self.E_F + self.E_P + self.internal_busses + self.E_L:
            if cp in b.comps.index:
                if cp_on_num_busses > 0:
                    msg = (
                        'The component ' + cp.label + ' is on multiple '
                        'busses in the exergy analysis. Make sure that no '
                        'component is connected to more than one of the '
                        'busses passed to the exergy_analysis method.')
                    logger.error(msg)
                    raise hlp.TESPyNetworkError(msg)
                # todo: E_bus als dict mit den versch. werten
                if b.comps.loc[cp, 'base'] == 'bus':
                    E_bus = sum(e for e in cp.E_bus.values() if e)
                    self.bus_data.loc[cp.label, 'E_P'] = E_bus
                    bus_efficiency = cp.calc_bus_efficiency(b)
                    E_F = E_bus / bus_efficiency
                    self.bus_data.loc[cp.label, 'E_F'] = E_F

                    if b in self.E_F:
                        self.network_data.loc['E_F'] += E_F
                    elif b in self.E_P:
                        self.network_data.loc['E_P'] -= E_F
                    elif b in self.E_L:
                        self.network_data.loc['E_L'] -= E_F

                    for key, value in cp.E_bus.items():
                        if value == 0:
                            continue

                        if key != "massless":
                            # this should be a source
                            category = categorize_fluids(cp.outl[0])
                        else:
                            category = "work"

                        if (cp.fkt_group, category) not in self.sankey_data[b.label].index:
                            self.sankey_data[b.label].loc[(cp.fkt_group, category), :] = 0

                        self.sankey_data[b.label].loc[(cp.fkt_group, category), key] += value / bus_efficiency
                else:
                    E_bus = sum(e for e in cp.E_bus.values() if e)
                    bus_efficiency = cp.calc_bus_efficiency(b)
                    E_P = E_bus * bus_efficiency
                    self.bus_data.loc[cp.label, 'E_P'] = E_P
                    self.bus_data.loc[cp.label, 'E_F'] = E_bus

                    if b in self.E_F:
                        self.network_data.loc['E_F'] -= E_P
                    elif b in self.E_P:
                        self.network_data.loc['E_P'] += E_P
                    elif b in self.E_L:
                        self.network_data.loc['E_L'] += E_P

                    for key, value in cp.E_bus.items():
                        if value == 0:
                            continue

                        if key != "massless":
                            # this should be a sink
                            category = categorize_fluids(cp.inl[0])
                        else:
                            category = "work"

                        if (b.label, category) not in self.sankey_data[cp.fkt_group].index:
                            self.sankey_data[cp.fkt_group].loc[(b.label, category), :] = 0

                        self.sankey_data[cp.fkt_group].loc[(b.label, category), key] += value * bus_efficiency

                self.bus_data.loc[cp.label, 'base'] = b.comps.loc[cp, 'base']
                self.bus_data.loc[cp.label, 'group'] = cp.fkt_group

                cp_on_num_busses += 1

        self.bus_data['E_D'] = self.bus_data['E_F'] - self.bus_data['E_P']

    def create_group_data(self):
        """Collect the component group exergy data."""
        for group in self.sankey_data.keys():
            E_D = self.aggregation_data[
                self.aggregation_data['group'] == group
            ]['E_D'].sum()
            self.sankey_data[group].loc[('E_D', "E_D"), :] = [0., 0., E_D]
        # establish connections for fuel exergy via bus balance
        for b in self.E_F:
            input_value = self.calculate_group_input_value(b.label)
            self.sankey_data['E_F'].loc[(b.label, "E_F"), self.exergy_cats] = (
                self.sankey_data[b.label].loc[:, self.exergy_cats].sum()
                - input_value.sum()
            )

        # establish connections for product exergy via bus balance
        for b in self.E_P:
            input_value = self.calculate_group_input_value(b.label)
            self.sankey_data[b.label].loc[('E_P', "E_P"), self.exergy_cats] = (
                input_value.sum()
                - self.sankey_data[b.label].loc[:, self.exergy_cats].sum()
            )

        # establish connections for exergy loss via bus balance
        for b in self.E_L:
            input_value = self.calculate_group_input_value(b.label)
            self.sankey_data[b.label].loc[('E_L', 'E_L'), self.exergy_cats] = (
                input_value.sum()
                - self.sankey_data[b.label].loc[:, self.exergy_cats].sum()
            )

        for fkt_group, data in self.sankey_data.items():
            mask = self.component_data['group'] == fkt_group
            comps = self.component_data.loc[mask].index
            for comp in comps:
                comp_obj = self.nw.get_comp(comp)
                sources = self.nw.conns[self.nw.conns['source'] == comp_obj]
                for conn in sources['object']:
                    if conn.target.label not in comps:
                        target_group = self.component_data.loc[conn.target.label, 'group']
                        target_value_chemical = (
                            conn.Ex_chemical
                            if hasattr(conn, "Ex_chemical") else 0.
                        )
                        target_value_physical = conn.Ex_physical
                        cat = categorize_fluids(conn)
                        if (target_group, cat) in data.index:
                            self.sankey_data[fkt_group].loc[
                                (target_group, cat), 'physical'] += target_value_physical
                            self.sankey_data[fkt_group].loc[
                                (target_group, cat), 'chemical'] += target_value_chemical
                        else:
                            self.sankey_data[fkt_group].loc[(target_group, cat), :] = [
                                target_value_chemical, target_value_physical, 0
                            ]

        # create overview of component groups
        self.group_data = pd.DataFrame(
            columns=['E_in', 'E_out', 'E_D'], dtype='float64'
        )

        for fkt_group in self.component_data['group'].unique():
            self.group_data.loc[fkt_group, 'E_in'] = (
                self.calculate_group_input_value(fkt_group).sum().sum()
            )
            self.group_data.loc[fkt_group, 'E_D'] = (
                self.sankey_data[fkt_group].loc[idx['E_D', :], self.exergy_cats].sum().sum()
            )

        # calculate missing values
        self.group_data['E_out'] = (
            self.group_data['E_in'] - self.group_data['E_D'])
        self.group_data['y_Dk'] = (
            self.group_data['E_D'] / self.network_data.loc['E_F'])
        self.group_data['y*_Dk'] = (
            self.group_data['E_D'] / self.network_data.loc['E_D'])

        # ToDo: Transform this into a test
        # assert self.group_data['E_D'].sum() == self.network_data["E_D"]

        return

    def calculate_group_input_value(self, group_label):
        """Calculate the total exergy input of a component group."""
        value = pd.DataFrame(
            columns=self.exergy_cats, index=[0], data=[[0., 0., 0.]]
        )
        for fkt_group, data in self.sankey_data.items():
            if group_label in data.index:
                value += data.loc[idx[group_label, :], self.exergy_cats].values
        return value

    def single_group_input(self, group_label, group_data):
        """Calculate the total exergy input of a component group."""
        inputs = []
        for fkt_group, data in group_data.items():
            if group_label in data.index.get_level_values("target_group") and fkt_group != group_label:
                inputs += [fkt_group]

        return inputs

    def remove_transit_groups(self, group_data):
        """Remove transit only component groups from sankey display.

        Method is recursively called if a group was removed from display to
        catch cases, where multiple groups are attached in line without
        any streams leaving the line.

        Parameters
        ----------
        group_data : dict
            Dictionary containing the modified component group data.
        """
        for fkt_group in group_data.copy().keys():
            if fkt_group in self.reserved_fkt_groups:
                continue
            source_groups = self.single_group_input(fkt_group, group_data)
            if len(source_groups) == 1 and len(group_data[fkt_group].index.get_level_values("target_group").unique()) == 1:

                source_group = source_groups[0]

                group_data[source_group] = group_data[source_group].add(group_data[fkt_group], fill_value=0)
                to_drop = group_data[source_group].index.get_level_values("target_group") == fkt_group
                group_data[source_group] = group_data[source_group].loc[~to_drop]
                del group_data[fkt_group]
                # recursive call in case multiple components are attached in
                # a line without any conversion
                self.remove_transit_groups(group_data)
                # exit to stop further iterations inside the groups
                return
            # remove groups without any connection
            elif len(source_groups) == 0 and len(group_data[fkt_group]) == 0:
                del group_data[fkt_group]

    def generate_plotly_sankey_input(
            self, node_order=[], colors={}, display_thresold=1e-3,
            disaggregate_flows=False):
        """Generate input data for sankey plots.

        Only exergy flow above the display threshold is included. All
        component groups with transit only (one input and one output) are cut
        out of the display.

        Parameters
        ----------
        node_order : list
            Order for the nodes in the sankey diagram (optional).
            In case no order is passed, a generic order will be used.

        colors : dict
            Dictionary containing a color for every stream type (optional).
            Stream type is the key, the color is the corresponding value.
            Stream types are:

            - :code:`E_F`, :code:`E_P`, :code:`E_L`, :code:`E_D`
            - names of the pure fluids of the tespy Network
            - :code:`mix` (used in case of any gas mixture)
            - labels of internal busses

        display_threshold : float
            Minimum value of flows to be displayed in the diagram, defaults to
            1e-3.

        disaggregate_flows : boolean
            Separate every flow by chemical, physical and massless exergy,
            defaults to False.

        Returns
        -------
        tuple
            Tuple containing the links and node_order for the plotly sankey
            diagram.
        """
        group_data = self.sankey_data.copy()
        cols = self.exergy_cats
        for fkt_group, data in self.sankey_data.items():
            mask = data.loc[:, cols].abs().sum(axis=1) >= display_thresold
            group_data[fkt_group] = group_data[fkt_group].loc[mask]

        self.remove_transit_groups(group_data)

        if len(node_order) == 0:
            node_order = (
                ['E_F'] + [b.label for b in self.E_F] +
                [fkt_group for fkt_group in self.group_data.index] +
                [b.label for b in self.internal_busses + self.E_P + self.E_L] +
                ['E_P', 'E_L', 'E_D']
            )
        else:
            missing = []
            for node in group_data.keys():
                if node not in node_order:
                    missing += [node]

            if len(missing) > 0:
                msg = (
                    'The list of nodes passed is missing the following '
                    'nodes: "' + '", "'.join(missing) + '".')
                logger.error(msg)
                raise ValueError(msg)

        colordict = {
            "E_F": "rgba(242, 142, 43, 0.90)",
            "E_P": "rgba(118, 183, 178, 0.90)",
            "E_D": "rgba(176, 122, 161, 0.90)",
            "E_L": "rgba(156, 117, 95, 0.90)",
            "combustion-gas": "rgba(237, 201, 72, 0.90)",
            "non-combustion-gas": "rgba(186, 176, 172, 0.90)",
            "two-phase-fluid": "rgba(89, 161, 79, 0.90)",
            "incompressible": "rgba(255, 157, 167, 0.90)",
            "work": "rgba(78, 121, 167, 0.90)",
            "heat": "rgba(225, 87, 89, 0.90)",
            np.nan: "rgba(100, 100, 100, 1.00)"
        }
        colordict.update(colors)

        links = {
            'source': [],
            'target': [],
            'value': [],
            'color': []
        }
        for fkt_group, data in group_data.items():
            source_id = node_order.index(fkt_group)
            for target in data.index:
                for col in cols:
                    # how to aggregate here?
                    if data.loc[target, col] > 0.:
                        links['source'] += [source_id]
                        links['target'] += [node_order.index(target[0])]
                        links['value'] += [data.loc[target, col]]
                        links['color'].append(colordict[target[1]])

        return links, node_order

    def print_results(
            self, sort_desc=True,
            busses=True, components=True, connections=True, groups=True,
            network=True, aggregation=True, Exe_Eco_An = False):
        r"""Print the results of the exergy analysis to prompt.

        - The results are sorted beginning with the component having the
          biggest exergy destruction by default.
        - Components with an exergy destruction smaller than 1000 W is not
          printed to prompt by default.

        Parameters
        ----------
        sort_des : boolean
            Sort the component results descending by exergy destruction.

        busses : boolean
            Print bus results, default value :code:`True`.

        components : boolean
            Print component results, default value :code:`True`.

        connections : boolean
            Print connection results, default value :code:`True`.

        network : boolean
            Print network results, default value :code:`True`.

        aggregation : boolean
            Print aggregated component results, default value :code:`True`.
        """
        if connections:
            print('##### RESULTS: Connection specific physical exergy and ' +
                  'chemical exergy #####')
            print(tabulate(
                self.connection_data, headers='keys',
                tablefmt='psql', floatfmt='.3e'))

        if components:
            df = self.component_data.copy()
            df = df.loc[:, df.columns != 'group']
            if sort_desc:
                df.sort_values(by=['E_D'], ascending=False, inplace=True)

            print('##### RESULTS: Component exergy analysis #####')
            print(tabulate(
                df, headers='keys', tablefmt='psql', floatfmt='.3e'))

        if busses:
            df = self.bus_data.copy()
            df = df.loc[:, (df.columns != 'group') & (df.columns != 'base')]
            if sort_desc:
                df.sort_values(by=['E_D'], ascending=False, inplace=True)

            print('##### RESULTS: Bus exergy analysis #####')
            print(tabulate(
                df, headers='keys', tablefmt='psql', floatfmt='.3e'))

        if aggregation:
            df = self.aggregation_data.copy()
            df = df.loc[:, df.columns != 'group']
            if sort_desc:
                df.sort_values(by=['E_D'], ascending=False, inplace=True)

            print('##### RESULTS: Aggregation of components and busses #####')
            print(tabulate(
                df, headers='keys', tablefmt='psql', floatfmt='.3e'))

        if network:
            print('##### RESULTS: Network exergy analysis #####')
            print(tabulate(
                self.network_data.to_frame().transpose(),
                headers='keys', tablefmt='psql', floatfmt='.3e',
                showindex=False))

        if groups:
            df = self.group_data.copy()
            if sort_desc:
                df.sort_values(by=['E_D'], ascending=False, inplace=True)

            print('##### RESULTS: Functional groups exergy flows #####')
            print(tabulate(
                df, headers='keys', tablefmt='psql', floatfmt='.3e'))

        """+F+F+F+F++++START++++F+F+F+F+"""
        # Exergoeconomic Results for Connections
        # creating data frame here bc this is after analysis where c and C values have been calculated already

        """
        conn_exergoec_data_cols = ['c_TOT', 'C_TOT']
        self.connection_exergoec_data = pd.DataFrame(
            columns=conn_exergoec_data_cols,
            dtype='float64'
        )
        for conn in self.nw.conns['object']:
            conn_exergoec_data = [
                conn.c_tot, conn.C_tot
            ]
            self.connection_exergoec_data.loc[conn.label] = conn_exergoec_data

        """

        # long version with subdivision of c and C
        conn_exergoec_data_cols = ['c_T', 'c_M', 'c_CH', 'c_tot', 'C_T', 'C_M', 'C_CH', 'C_tot']
        self.connection_exergoec_data = pd.DataFrame(
            columns=conn_exergoec_data_cols,
            dtype='float64'
        )
        for conn in self.nw.conns['object']:
            conn_exergoec_data = [
                conn.c_therm, conn.c_mech, conn.c_chemical, conn.c_tot,
                conn.C_therm, conn.C_mech, conn.C_chemical, conn.C_tot
            ]
            self.connection_exergoec_data.loc[conn.label] = conn_exergoec_data

        if connections and Exe_Eco_An:
            print('##### RESULTS: Connection exergoeconomic analysis #####')
            print(tabulate(
                self.connection_exergoec_data, headers='keys',
                tablefmt='psql', floatfmt='.3e'))


        # Exergoeconomic Results for Components
        # creating data frame here bc this is after analysis where c and C values have been calculated already
        comp_exergoec_data_cols = ['C_F', 'C_P', 'C_D', 'C_L', 'C_bus', 'Z', 'r', 'f']
        self.component_exergoec_data = pd.DataFrame(
            columns=comp_exergoec_data_cols,
            dtype='float64'
        )
        for comp in self.nw.comps['object']:
            comp_exergoec_data = [
                comp.C_F, comp.C_P, comp.C_D, comp.C_L, comp.C_bus, comp.Z_costs, comp.r, comp.f
            ]
            self.component_exergoec_data.loc[comp.label] = comp_exergoec_data

        if components and Exe_Eco_An:
            print('##### RESULTS: Component exergoeconomic analysis #####')
            print(tabulate(
                self.component_exergoec_data, headers='keys',
                tablefmt='psql', floatfmt='.3e'))

        if network and Exe_Eco_An:
            print('##### RESULTS: Network exergoeconomic analysis #####')
            print(tabulate(
                self.network_data_exergo.to_frame().transpose(),
                headers='keys', tablefmt='psql', floatfmt='.3e',
                showindex=False))

        """+F+F+F+F++++END++++F+F+F+F+"""
