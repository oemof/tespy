# -*- coding: utf-8

"""Module for connections between components.

Components in this module:

    - :func:`tespy.connections.connection` (mass flow)
    - :func:`tespy.connections.bus` (energy flow)
    - :func:`tespy.connections.ref` (referenced fluid states container)


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/connections.py

SPDX-License-Identifier: MIT
"""

import logging
import warnings

import numpy as np
import pandas as pd

from tespy.components.components import component
from tespy.tools.characteristics import char_line
from tespy.tools.data_containers import dc_cp
from tespy.tools.data_containers import dc_flu
from tespy.tools.data_containers import dc_prop
from tespy.tools.data_containers import dc_simple
from tespy.tools.helpers import TESPyConnectionError

# pass the warning messages to the logger
logging.captureWarnings(True)


class connection:
    r"""
    Class connection is the container for fluid properties between components.

    Parameters
    ----------
    m : float/tespy.connections.ref/tespy.tools.helpers.dc_prop
        Mass flow specification.

    m0 : float
        Starting value specification for mass flow.

    p : float/tespy.connections.ref/tespy.tools.helpers.dc_prop
        Pressure specification.

    p0 : float
        Starting value specification for pressure.

    h : float/tespy.connections.ref/tespy.tools.helpers.dc_prop
        Enthalpy specification.

    h0 : float
        Starting value specification for enthalpy.

    fluid : dict/tespy.tools.helpers.dc_flu
        Fluid compostition specification.

    fluid0 : dict
        Starting value specification for fluid compostition.

    fluid_balance : boolean
        Fluid balance equation specification.

    x : float/tespy.tools.helpers.dc_prop
        Gas phase mass fraction specification.

    T : float/tespy.connections.ref/tespy.tools.helpers.dc_prop
        Temperature specification.

    Td_bp : float/tespy.tools.helpers.dc_prop
        Temperature difference to boiling point at pressure corresponding
        pressure of this connection in K.

    v : float/tespy.tools.helpers.dc_prop
        Volumetric flow specification.

    state : str
        State of the pure fluid on this connection: liquid ('l') or gaseous
        ('g').

    design : list
        List containing design parameters (stated as string).

    offdesign : list
        List containing offdesign parameters (stated as string).

    design_path : str
        Path to individual design case for this connection.

    local_offdesign : boolean
        Treat this connection in offdesign mode in a design calculation.

    local_design : boolean
        Treat this connection in design mode in an offdesign calculation.

    printout: boolean
        Include this connection in the network's results printout.

    label : str
        Label of the connection. The default value is:
        :code:`'source:source_id_target:target_id'`.

    Note
    ----
    - The fluid balance parameter applies a balancing of the fluid vector on
      the specified conntion to 100 %. For example, you have four fluid
      components (a, b, c and d) in your vector, you set two of them
      (a and b) and want the other two (components c and d) to be a result of
      your calculation. If you set this parameter to True, the equation
      (0 = 1 - a - b - c - d) will be applied.

    - The specification of values for design and/or offdesign is used for
      automatic switch from design to offdesign calculation: All parameters
      given in 'design', e.g. :code:`design=['T', 'p']`, are unset in any
      offdesign calculation, parameters given in 'offdesign' are set for
      offdesign calculation.

    Example
    -------
    This example shows how to create connections and specify parameters. First
    create the required components and connect them in the next step. After
    that, it is possible specify parameters with the :code:`set_attr` method.

    >>> from tespy.components import sink, source
    >>> from tespy.connections import connection, ref
    >>> from tespy.tools.data_containers import dc_flu, dc_prop
    >>> import numpy as np
    >>> so1 = source('source1')
    >>> so2 = source('source2')
    >>> si1 = sink('sink1')
    >>> si2 = sink('sink2')
    >>> so_si1 = connection(so1, 'out1', si1, 'in1')
    >>> so_si2 = connection(so2, 'out1', si2, 'in1')

    There are different ways of setting parameters on connections: Specify
        - a numeric value  (for attributes mass flow, pressure and enthalpy)
        - a numeric starting value (for attributes mass flow, pressure and
          enthalpy)
        - a dictionary (for attributes fluid and fluid0)
        - a boolean value (for attributes fluid_balance, local_design,
          local_offdesign).
        - a referenced value (mass flow, pressure, temperature, enthalpy).
        - a data_container (for attributes fluid and fluid0 dc_flu, for other
          fluid propertie attributes dc_prop).
        - numpy.nan or None (unsetting a value).
        - a string (for attributes design_paht and state).
        - a list (for attributes design and offdesign).

    >>> so_si1.set_attr(v=0.012, m0=10, p=5, h=400, fluid={'H2O': 1, 'N2': 0})

    Individual specification of pressure with a data container. Setting the
    unit individually will overwrite the network's unit specification for this
    connection.

    >>> p = dc_prop(val=50, val_set=True, unit='bar', unit_set=True)
    >>> so_si2.set_attr(m=ref(so_si1, 2, -5), p=p, h0=700, T=200,
    ... fluid={'N2': 1}, fluid_balance=True)

    The set_attr method automatically converts your input in data_container
    information.

    >>> type(so_si1.v)
    <class 'tespy.tools.data_containers.dc_prop'>
    >>> type(so_si1.fluid)
    <class 'tespy.tools.data_containers.dc_flu'>

    If you want get a spcific value use the logic: connection.property.*.
    Aditionally, it is possible to use the :code:`get_attr` method.

    >>> so_si1.m.val0
    10
    >>> so_si1.m.get_attr('val_set')
    False
    >>> type(so_si2.m.ref)
    <class 'tespy.connections.ref'>
    >>> so_si2.fluid.get_attr('balance')
    True
    >>> so_si2.m.ref.get_attr('d')
    -5
    >>> so_si2.m.ref_set
    True
    >>> type(so_si2.m.ref.get_attr('obj'))
    <class 'tespy.connections.connection'>

    Unset the specified temperature and specify temperature difference to
    boiling point instead.

    >>> so_si2.T.val_set
    True
    >>> so_si2.set_attr(Td_bp=5, T=np.nan)
    >>> so_si2.T.val_set
    False
    >>> so_si2.Td_bp.val
    5
    >>> so_si2.set_attr(Td_bp=None)
    >>> so_si2.Td_bp.val_set
    False

    Specify the state keyword: The fluid will be forced to liquid or gaseous
    state in this case.

    >>> so_si2.set_attr(state='l')
    >>> so_si2.state.is_set
    True
    >>> so_si2.set_attr(state=np.nan)
    >>> so_si2.state.is_set
    False
    >>> so_si2.set_attr(state='g')
    >>> so_si2.state.is_set
    True
    >>> so_si2.set_attr(state=None)
    >>> so_si2.state.is_set
    False
    >>> so_si2.set_attr(label='myconnection')
    >>> so_si2.label
    'myconnection'
    """

    def __init__(self, comp1, outlet_id, comp2, inlet_id, **kwargs):

        # check input parameters
        if not (isinstance(comp1, component) and
                isinstance(comp2, component)):
            msg = ('Error creating connection. Check if comp1, comp2 are of '
                   'type component.')
            logging.error(msg)
            raise TypeError(msg)

        if comp1 == comp2:
            msg = ('Error creating connection. Cannot connect component ' +
                   comp1.label + ' to itself.')
            logging.error(msg)
            raise TESPyConnectionError(msg)

        if outlet_id not in comp1.outlets():
            msg = ('Error creating connection. Specified oulet_id (' +
                   outlet_id + ') is not valid for component ' +
                   comp1.component() + '. Valid ids are: ' +
                   str(comp1.outlets()) + '.')
            logging.error(msg)
            raise ValueError(msg)

        if inlet_id not in comp2.inlets():
            msg = (
                'Error creating connection. Specified inlet_id (' + inlet_id +
                ') is not valid for component ' + comp2.component() +
                '. Valid ids are: ' + str(comp2.inlets()) + '.')
            logging.error(msg)
            raise ValueError(msg)

        self.label = (
            comp1.label + ':' + outlet_id + '_' + comp2.label + ':' + inlet_id)

        # set specified values
        self.source = comp1
        self.source_id = outlet_id
        self.target = comp2
        self.target_id = inlet_id

        # defaults
        self.new_design = True
        self.design_path = None
        self.design = []
        self.offdesign = []
        self.local_design = False
        self.local_offdesign = False
        self.printout = True

        # set default values for kwargs
        self.variables = self.attr()
        self.variables0 = [x + '0' for x in self.variables.keys()]
        self.__dict__.update(self.variables)
        self.set_attr(**kwargs)

        msg = (
            'Created connection ' + self.source.label + ' (' + self.source_id +
            ') -> ' + self.target.label + ' (' + self.target_id + ').')
        logging.debug(msg)

    def set_attr(self, **kwargs):
        r"""
        Set, reset or unset attributes of a connection.

        Parameters
        ----------
        m : float/tespy.connections.ref/tespy.tools.helpers.dc_prop
            Mass flow specification.

        m0 : float
            Starting value specification for mass flow.

        p : float/tespy.connections.ref/tespy.tools.helpers.dc_prop
            Pressure specification.

        p0 : float
            Starting value specification for pressure.

        h : float/tespy.connections.ref/tespy.tools.helpers.dc_prop
            Enthalpy specification.

        h0 : float
            Starting value specification for enthalpy.

        fluid : dict/tespy.tools.helpers.dc_flu
            Fluid compostition specification.

        fluid0 : dict
            Starting value specification for fluid compostition.

        fluid_balance : boolean
            Fluid balance equation specification.

        x : float/tespy.tools.helpers.dc_prop
            Gas phase mass fraction specification.

        T : float/tespy.connections.ref/tespy.tools.helpers.dc_prop
            Temperature specification.

        Td_bp : float/tespy.tools.helpers.dc_prop
            Temperature difference to boiling point at pressure corresponding
            pressure of this connection in K.

        v : float/tespy.tools.helpers.dc_prop
            Volumetric flow specification.

        state : str
            State of the pure fluid on this connection: liquid ('l') or gaseous
            ('g').

        design : list
            List containing design parameters (stated as string).

        offdesign : list
            List containing offdesign parameters (stated as string).

        design_path : str
            Path to individual design case for this connection.

        local_offdesign : boolean
            Treat this connection in offdesign mode in a design calculation.

        local_design : boolean
            Treat this connection in design mode in an offdesign calculation.

        printout: boolean
            Include this connection in the network's results printout.

        label : str
            Label of the connection. The default value is:
            :code:`'source:source_id_target:target_id'`.

        Note
        ----
        - The fluid balance parameter applies a balancing of the fluid vector
          on the specified conntion to 100 %. For example, you have four fluid
          components (a, b, c and d) in your vector, you set two of them
          (a and b) and want the other two (components c and d) to be a result
          of your calculation. If you set this parameter to True, the equation
          (0 = 1 - a - b - c - d) will be applied.

        - The specification of values for design and/or offdesign is used for
          automatic switch from design to offdesign calculation: All parameters
          given in 'design', e.g. :code:`design=['T', 'p']`, are unset in any
          offdesign calculation, parameters given in 'offdesign' are set for
          offdesign calculation.

        - The property state is applied on pure fluids only. If you specify the
          desired state of the fluid at a connection the convergence check will
          adjust the enthalpy values of that connection for the first
          iterations in order to meet the state requirement.
        """
        # set specified values
        for key in kwargs:
            if key in self.variables.keys() or key in self.variables0:
                # fluid specification
                if 'fluid' in key:
                    if isinstance(kwargs[key], dict):
                        # starting values
                        if key in self.variables0:
                            self.fluid.set_attr(val0=kwargs[key].copy())
                        # specified parameters
                        else:
                            self.fluid.set_attr(val=kwargs[key].copy())
                            self.fluid.set_attr(
                                val_set={f: True for f in kwargs[key].keys()})

                    elif isinstance(kwargs[key], dc_flu):
                        # data container for fluids
                        self.__dict__.update({key: kwargs[key]})

                    else:
                        # bad datatype
                        msg = (
                            'Datatype for fluid vector specification must be '
                            'tespy.tools.data_containers.dc_flu or dict.')
                        logging.error(msg)
                        raise TypeError(msg)

                elif key == 'state':
                    if kwargs[key] in ['l', 'g']:
                        self.state.set_attr(val=kwargs[key], is_set=True)
                    elif isinstance(kwargs[key], dc_simple):
                        self.state = kwargs[key]
                    elif kwargs[key] is None:
                        self.get_attr(key).set_attr(is_set=False)
                    elif (isinstance(kwargs[key], float) or
                          isinstance(kwargs[key], np.float64) or
                          isinstance(kwargs[key], np.int64) or
                          isinstance(kwargs[key], int)):
                        if np.isnan(kwargs[key]):
                            self.get_attr(key).set_attr(is_set=False)
                        else:
                            msg = (
                                'To unset the state specification either use '
                                'np.nan or None.')
                            logging.error(msg)
                            raise ValueError(msg)
                    else:
                        msg = (
                            'Keyword argument "state" must either be '
                            '"l" or "g" or be None or np.nan.')
                        logging.error(msg)
                        raise TypeError(msg)

                elif kwargs[key] is None:
                    self.get_attr(key).set_attr(val_set=False)
                    self.get_attr(key).set_attr(ref_set=False)

                elif (isinstance(kwargs[key], float) or
                      isinstance(kwargs[key], np.float64) or
                      isinstance(kwargs[key], np.int64) or
                      isinstance(kwargs[key], int)):
                    if np.isnan(kwargs[key]):
                        self.get_attr(key).set_attr(val_set=False)
                        self.get_attr(key).set_attr(ref_set=False)
                    else:
                        # value specification
                        if key in self.variables:
                            self.get_attr(key).set_attr(
                                val_set=True,
                                val=kwargs[key])
                        # starting value specification
                        else:
                            self.get_attr(key.replace('0', '')).set_attr(
                                val0=kwargs[key])

                # reference object
                elif isinstance(kwargs[key], ref):
                    if key in ['x', 'v', 'Td_bp']:
                        msg = (
                            'References for volumetric flow, vapor mass '
                            'fraction and subcooling/superheating are not '
                            'implemented.')
                        logging.error(msg)
                        raise NotImplementedError(msg)
                    else:
                        self.get_attr(key).set_attr(ref=kwargs[key])
                        self.get_attr(key).set_attr(ref_set=True)

                # data container specification
                elif isinstance(kwargs[key], dc_prop):
                    self.__dict__.update({key: kwargs[key]})

                # invalid datatype for keyword
                else:
                    msg = 'Bad datatype for keyword argument ' + key + '.'
                    logging.error(msg)
                    raise TypeError(msg)

            # fluid balance
            elif key == 'fluid_balance':
                if isinstance(kwargs[key], bool):
                    self.get_attr('fluid').set_attr(balance=kwargs[key])
                else:
                    msg = (
                        'Datatype for keyword argument fluid_balance must be '
                        'boolean.')
                    logging.error(msg)
                    raise TypeError(msg)

            # design/offdesign parameter list
            elif key == 'design' or key == 'offdesign':
                if not isinstance(kwargs[key], list):
                    msg = 'Please provide the ' + key + ' parameters as list!'
                    logging.error(msg)
                    raise TypeError(msg)
                elif set(kwargs[key]).issubset(self.variables.keys()):
                    self.__dict__.update({key: kwargs[key]})
                else:
                    msg = (
                        'Available parameters for (off-)design specification '
                        'are: ' + str(self.variables.keys()) + '.')
                    logging.error(msg)
                    raise ValueError(msg)

            # design path
            elif key == 'design_path':
                if isinstance(kwargs[key], str):
                    self.__dict__.update({key: kwargs[key]})
                elif np.isnan(kwargs[key]):
                    self.design_path = None
                else:
                    msg = (
                        'Please provide the design_path parameter as string '
                        'or as nan.')
                    logging.error(msg)
                    raise TypeError(msg)

                self.new_design = True

            # other boolean keywords
            elif key in ['printout', 'local_design', 'local_offdesign']:
                if not isinstance(kwargs[key], bool):
                    msg = ('Please provide the ' + key + ' as boolean.')
                    logging.error(msg)
                    raise TypeError(msg)
                else:
                    self.__dict__.update({key: kwargs[key]})

            # label
            elif key == 'label':
                if isinstance(kwargs[key], str):
                    self.__dict__.update({key: kwargs[key]})
                else:
                    msg = 'Please provide the label as string.'
                    logging.error(msg)
                    raise TypeError(msg)

            # invalid keyword
            else:
                msg = 'Connection has no attribute ' + key + '.'
                logging.error(msg)
                raise KeyError(msg)

    def get_attr(self, key):
        r"""
        Get the value of a connection's attribute.

        Parameters
        ----------
        key : str
            The attribute you want to retrieve.

        Returns
        -------
        out :
            Specified attribute.
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            msg = 'Connection has no attribute \"' + key + '\".'
            logging.error(msg)
            raise KeyError(msg)

    @staticmethod
    def attr():
        r"""
        Return available attributes of a connection.

        Returns
        -------
        out : list
            List of available attributes of a connection.
        """
        return {'m': dc_prop(), 'p': dc_prop(), 'h': dc_prop(), 'T': dc_prop(),
                'x': dc_prop(), 'v': dc_prop(), 'vol': dc_prop(),
                's': dc_prop(),
                'fluid': dc_flu(), 'Td_bp': dc_prop(), 'state': dc_simple()}

    def to_flow(self):
        r"""
        Return the SI-values for the network variables.

        Returns
        -------
        out : list
            List of mass flow and fluid property information.
        """
        return [self.m.val_SI, self.p.val_SI, self.h.val_SI, self.fluid.val]

    def to_flow_design(self):
        r"""
        Return the SI-values for the network variables at design point.

        Returns
        -------
        out : list
            List of mass flow and fluid property information.
        """
        return [self.m.design, self.p.design, self.h.design, self.fluid.design]


class bus:
    r"""
    A bus is used to connect different energy flows.

    Parameters
    ----------
    label : str
        Label for the bus.

    P : float
        Total power/heat flow specification for bus, :math:`P \text{/W}`.

    printout : bool
        Print the results of this bus to prompt with the
        :func:`tespy.networks.networks.network.print_results` method. Standard
        value is :code:`True`.

    Example
    -------
    Busses are used to connect energy flow of different components. They can
    also be used to introduce efficiencies of energy conversion, e.g. in
    motors, generator or boilers. This example takes the combustion engine
    example at :func:`tespy.components.combustion.combustion_engine` and adds
    a flue gas cooler and a circulation pump for the cooling water. Then
    busses for heat output, thermal input and electricity output are
    implemented.

    >>> from tespy.components import (sink, source, combustion_engine,
    ... heat_exchanger, merge, splitter, pump)
    >>> from tespy.connections import connection, ref, bus
    >>> from tespy.networks import network
    >>> import numpy as np
    >>> import shutil
    >>> fluid_list = ['Ar', 'N2', 'O2', 'CO2', 'CH4', 'H2O']
    >>> nw = network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ... p_range=[0.5, 10], T_range=[150, 1200], iterinfo=False)
    >>> amb = source('ambient')
    >>> sf = source('fuel')
    >>> fg = sink('flue gas outlet')
    >>> cw_in = source('cooling water inlet')
    >>> sp = splitter('cooling water splitter', num_out=2)
    >>> me = merge('cooling water merge', num_in=2)
    >>> cw_out = sink('cooling water outlet')
    >>> fgc = heat_exchanger('flue gas cooler')
    >>> pu = pump('cooling water pump')
    >>> chp = combustion_engine(label='internal combustion engine')
    >>> amb_comb = connection(amb, 'out1', chp, 'in3')
    >>> sf_comb = connection(sf, 'out1', chp, 'in4')
    >>> comb_fgc = connection(chp, 'out3', fgc, 'in1')
    >>> fgc_fg = connection(fgc, 'out1', fg, 'in1')
    >>> nw.add_conns(sf_comb, amb_comb, comb_fgc, fgc_fg)
    >>> cw_pu = connection(cw_in, 'out1', pu, 'in1')
    >>> pu_sp = connection(pu, 'out1', sp, 'in1')
    >>> sp_chp1 = connection(sp, 'out1', chp, 'in1')
    >>> sp_chp2 = connection(sp, 'out2', chp, 'in2')
    >>> chp1_me = connection(chp, 'out1', me, 'in1')
    >>> chp2_me = connection(chp, 'out2', me, 'in2')
    >>> me_fgc = connection(me, 'out1', fgc, 'in2')
    >>> fgc_cw = connection(fgc, 'out2', cw_out, 'in1')
    >>> nw.add_conns(cw_pu, pu_sp, sp_chp1, sp_chp2, chp1_me, chp2_me, me_fgc,
    ... fgc_cw)
    >>> chp.set_attr(pr1=0.99, lamb=1.0,
    ... design=['pr1'], offdesign=['zeta1'])
    >>> fgc.set_attr(pr1=0.999, pr2=0.98, design=['pr1', 'pr2'],
    ... offdesign=['zeta1', 'zeta2', 'kA_char'])
    >>> pu.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char'])
    >>> amb_comb.set_attr(p=5, T=30, fluid={'Ar': 0.0129, 'N2': 0.7553,
    ... 'H2O': 0, 'CH4': 0, 'CO2': 0.0004, 'O2': 0.2314})
    >>> sf_comb.set_attr(m0=0.1, T=30, fluid={'CO2': 0, 'Ar': 0, 'N2': 0,
    ... 'O2': 0, 'H2O': 0, 'CH4': 1})
    >>> cw_pu.set_attr(p=3, T=60, fluid={'CO2': 0, 'Ar': 0, 'N2': 0,
    ... 'O2': 0, 'H2O': 1, 'CH4': 0})
    >>> sp_chp2.set_attr(m=ref(sp_chp1, 1, 0))

    Cooling water mass flow is calculated given the feed water temperature
    (90 °C). The pressure at the cooling water outlet should be identical to
    pressure before pump. The flue gases of the combustion engine leave the
    flue gas cooler at 120 °C.

    >>> fgc_cw.set_attr(p=ref(cw_pu, 1, 0), T=90)
    >>> fgc_fg.set_attr(T=120, design=['T'])

    Now add the busses, pump and combustion engine generator will get a
    characteristic function for conversion efficiency. In case of the
    combustion engine we have to individually choose the bus parameter as there
    are several available (P, TI, Q1, Q2, Q). For heat exchangers or
    turbomachinery the bus parameter is unique. The characteristic function for
    the flue gas cooler is set to -1, as the heat transferred is always
    negative in class heat_exchanger by definition. Instead of specifying the
    combustion engine power output we will define the total electrical power
    output at the bus.

    >>> load = np.array([0.2, 0.4, 0.6, 0.8, 1, 1.2])
    >>> gen_efficiency = np.array([0.9, 0.94, 0.97, 0.99, 1, 0.99]) * 0.98
    >>> mot_efficiency = 1 / (
    ... np.array([0.9, 0.94, 0.97, 0.99, 1, 0.99]) * 0.98)
    >>> gen = char_line(x=load, y=gen_efficiency)
    >>> mot = char_line(x=load, y=mot_efficiency)
    >>> power_bus = bus('total power output', P=10e6)
    >>> heat_bus = bus('total heat input')
    >>> fuel_bus = bus('thermal input')

    You can check, if the bus value is set or not identically to connections.
    Unsetting is possible using :code:`np.nan` or :code:`None`. For
    demonstration we will specify a value for the heat bus and unset it again.

    >>> heat_bus.P.is_set
    False
    >>> heat_bus.set_attr(P=1e5)
    >>> heat_bus.P.is_set
    True
    >>> heat_bus.set_attr(P=None)
    >>> heat_bus.P.is_set
    False

    >>> power_bus.add_comps({'comp': chp, 'char': gen, 'param': 'P'},
    ... {'comp': pu, 'char': mot})
    >>> heat_bus.add_comps({'comp': chp, 'param': 'Q'},
    ... {'comp': fgc, 'char': -1})
    >>> fuel_bus.add_comps({'comp': chp, 'param': 'TI'},
    ... {'comp': pu, 'char': mot})
    >>> nw.add_busses(power_bus, heat_bus, fuel_bus)
    >>> mode = 'design'
    >>> nw.solve(mode=mode)
    >>> nw.save('tmp')

    The heat bus characteristic for the flue gas cooler has automatically been
    transformed into an array. The total heat output can be seen in the
    enthalpy rise of the cooling water in the combustion engine and the flue
    gas cooler.

    >>> heat_bus.comps.loc[fgc]['char'].x
    array([0., 3.])
    >>> heat_bus.comps.loc[fgc]['char'].y
    array([-1., -1.])
    >>> round(chp.ti.val)
    25813247.0
    >>> round(chp.Q1.val + chp.Q2.val, 0)
    8896898.0
    >>> round(fgc_cw.m.val_SI * (fgc_cw.h.val_SI - pu_sp.h.val_SI), 0)
    12553049.0
    >>> round(heat_bus.P.val, 0)
    12553049.0
    >>> round(1 / pu.calc_bus_efficiency(power_bus), 2)
    0.98
    >>> power_bus.set_attr(P=7.5e6)
    >>> mode = 'offdesign'
    >>> nw.solve(mode=mode, design_path='tmp', init_path='tmp')
    >>> round(chp.ti.val)
    21187528.0
    >>> round(chp.P.val / chp.P.design, 3)
    0.761
    >>> round(1 / pu.calc_bus_efficiency(power_bus), 3)
    0.967
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def __init__(self, label, **kwargs):

        self.comps = pd.DataFrame(
            columns=['param', 'P_ref', 'char', 'efficiency', 'base'])

        self.label = label
        self.P = dc_cp(val=np.nan, is_set=False)
        self.char = char_line(x=np.array([0, 3]), y=np.array([1, 1]))
        self.printout = True

        self.set_attr(**kwargs)

        msg = 'Created bus ' + self.label + '.'
        logging.debug(msg)

    def set_attr(self, **kwargs):
        r"""
        Set, reset or unset attributes of a bus object.

        Parameters
        ----------
        label : str
            Label for the bus.

        P : float
            Total power/heat flow specification for bus, :math:`P \text{/W}`.

        printout : bool
            Print the results of this bus to prompt with the
            :func:`tespy.networks.networks.network.print_results` method.
            Standard value is :code:`True`.

        Note
        ----
        Specify :math:`P=\text{nan}`, if you want to unset the value of P.
        """
        for key in kwargs:
            if key == 'P':
                if (isinstance(kwargs[key], float) or
                        isinstance(kwargs[key], np.float64) or
                        isinstance(kwargs[key], np.int64) or
                        isinstance(kwargs[key], int)):
                    if np.isnan(kwargs[key]):
                        self.P.set_attr(is_set=False)
                    else:
                        self.P.set_attr(val=kwargs[key], is_set=True)
                elif kwargs[key] is None:
                    self.P.set_attr(is_set=False)
                else:
                    msg = ('Keyword argument ' + key + ' must be numeric.')
                    logging.error(msg)
                    raise TypeError(msg)

            elif key == 'printout':
                if not isinstance(kwargs[key], bool):
                    msg = ('Please provide the ' + key + ' as boolean.')
                    logging.error(msg)
                    raise TypeError(msg)
                else:
                    self.__dict__.update({key: kwargs[key]})

            # invalid keyword
            else:
                msg = 'A bus has no attribute ' + key + '.'
                logging.error(msg)
                raise KeyError(msg)

    def get_attr(self, key):
        r"""
        Get the value of a busses attribute.

        Parameters
        ----------
        key : str
            The attribute you want to retrieve.

        Returns
        -------
        out :
            Specified attribute.
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            msg = 'Bus ' + self.label + ' has no attribute ' + key + '.'
            logging.error(msg)
            raise KeyError(msg)

    def add_comps(self, *args):
        r"""
        Add components to a bus.

        Parameters
        ----------
        *args : dict
            Dictionaries containing the component information to be added to
            the bus. The information are described below.

        Note
        ----
        **Required Key**

        - comp (tespy.components.components.component): Component you want to
          add to the bus.

        **Optional Keys**

        - param (str): Bus parameter, optional.

            - You do not need to provide a parameter, if the component only has
              one option for the bus (turbomachines, heat exchangers,
              combustion chamber).
            - For instance, you do neet do provide a parameter, if you want to
              add a combustion engine ('Q', 'Q1', 'Q2', 'TI', 'P', 'Qloss').

        - char (float/tespy.components.characteristics.characteristics):
          Characteristic function for this components share to the bus value,
          optional.

            - If you do not provide a characteristic line at all, TESPy assumes
              a constant factor of 1.
            - If you provide a numeric value instead of a characteristic line,
              TESPy takes this numeric value as a constant factor.
            - Provide a TESPy.characteristic (cmp_char), if you want the factor
              to follow a characteristic line.

        - P_ref (float): Energy flow specification for reference case,
          :math:`P \text{/W}`, optional.
        - base (str): Base value for characteristic line and efficiency
          calculation. The base can either be :code:`'component'` (default) or
          :code:`'bus'`.

            - In case you choose :code:`'component'`, the characteristic line
              input will follow the value of the component's bus function and
              the efficiency definition is
              :math:`\eta=\frac{P_\mathrm{bus}}{P_\mathrm{component}}`.
            - In case you choose :code:`'bus'`, the characteristic line
              input will follow the bus value of the component and the
              efficiency definition is
              :math:`\eta=\frac{P_\mathrm{component}}{P_\mathrm{bus}}`.
        """
        for c in args:
            if isinstance(c, dict):
                if 'comp' in c.keys() or 'c' in c.keys():
                    if 'c' in c.keys():
                        comp = c['c']
                        msg = (
                            'The dictionary keyword "c" will be removed in '
                            'TESPy version 0.4.0. Please use "comp" instead.')
                        warnings.warn(msg, FutureWarning, stacklevel=2)
                    else:
                        comp = c['comp']
                    # default values
                    if isinstance(comp, component):
                        self.comps.loc[comp] = [
                            None, np.nan, self.char, np.nan, 'component']
                    else:
                        msg = 'Keyword "comp" must hold a TESPy component.'
                        logging.error(msg)
                        raise TypeError(msg)
                else:
                    msg = 'You must provide the component "comp".'
                    logging.error(msg)
                    raise TypeError(msg)

                for k, v in c.items():
                    if k == 'param' or k == 'p':
                        if k == 'p':
                            msg = (
                                'The dictionary keyword "p" will be removed '
                                'TESPy version 0.4.0. Please use "param" '
                                'instead.')
                            warnings.warn(msg, FutureWarning, stacklevel=2)
                        if isinstance(v, str) or v is None:
                            self.comps.loc[comp, 'param'] = v
                        else:
                            msg = (
                                'The bus parameter selection must be a '
                                'string (at bus ' + self.label + ').')
                            logging.error(msg)
                            raise TypeError(msg)

                    elif k == 'char':
                        if isinstance(v, char_line):
                            self.comps.loc[comp, 'char'] = v
                        elif (isinstance(v, float) or
                              isinstance(v, np.float64) or
                              isinstance(v, np.int64) or
                              isinstance(v, int)):
                            x = np.array([0, 3])
                            y = np.array([1, 1]) * v
                            self.comps.loc[comp, 'char'] = (
                                    char_line(x=x, y=y))
                        else:
                            msg = (
                                'Char must be a number or a TESPy '
                                'characteristics char line.')
                            logging.error(msg)
                            raise TypeError(msg)

                    elif k == 'P_ref':
                        if (v is None or isinstance(v, float) or
                                isinstance(v, np.float64) or
                                isinstance(v, np.int64) or
                                isinstance(v, int)):
                            self.comps.loc[comp, 'P_ref'] = v
                        else:
                            msg = 'Reference value must be numeric.'
                            logging.error(msg)
                            raise TypeError(msg)

                    elif k == 'base':
                        if v in ['bus', 'component']:
                            self.comps.loc[comp, 'base'] = v
                        else:
                            msg = (
                                'The base value must be "bus" or "component".')
                            logging.error(msg)
                            raise ValueError(msg)

            else:
                msg = (
                    'Provide arguments as dictionaries. See the documentation '
                    'of bus.add_comps() for more information.')
                logging.error(msg)
                raise TypeError(msg)

            msg = (
                'Added component ' + comp.label + ' to bus ' +
                self.label + '.')
            logging.debug(msg)


class ref:
    r"""
    Reference fluid properties from one connection to another connection.

    Parameters
    ----------
    obj : tespy.connections.connection
        Connection to be referenced.

    f : float
        Factor to multiply specified property with.

    d : float
        Delta to add after multiplication.

    Note
    ----
    Reference the mass flow of one connection :math:`\dot{m}` to another mass
    flow :math:`\dot{m}_{ref}`

    .. math::

        \dot{m} = \dot{m}_\mathrm{ref} \cdot f + d

    """

    def __init__(self, ref_obj, factor, delta):

        if not isinstance(ref_obj, connection):
            msg = 'First parameter must be object of type connection.'
            logging.error(msg)
            raise TypeError(msg)

        if not (isinstance(factor, int) or isinstance(factor, float)):
            msg = 'Second parameter must be of type int or float.'
            logging.error(msg)
            raise TypeError(msg)

        if not (isinstance(delta, int) or isinstance(delta, float)):
            msg = 'Thrid parameter must be of type int or float.'
            logging.error(msg)
            raise TypeError(msg)

        self.obj = ref_obj
        self.f = factor
        self.d = delta

        msg = ('Created reference object with factor ' + str(self.f) +
               ' and delta ' + str(self.d) + ' referring to connection ' +
               ref_obj.source.label + ' (' + ref_obj.source_id + ') -> ' +
               ref_obj.target.label + ' (' + ref_obj.target_id + ').')
        logging.debug(msg)

    def get_attr(self, key):
        r"""
        Get the value of a reference attribute.

        Parameters
        ----------
        key : str
            The attribute you want to retrieve.

        Returns
        -------
        out :
            Specified attribute.
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            msg = 'Reference has no attribute \"' + key + '\".'
            logging.error(msg)
            raise KeyError(msg)
