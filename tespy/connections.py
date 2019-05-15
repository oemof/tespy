# -*- coding: utf-8

"""
.. module:: connections
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import numpy as np
import pandas as pd

import logging

from tespy.tools.helpers import (TESPyConnectionError, data_container, dc_prop, dc_flu, dc_cp, dc_simple)
from tespy.components import components as cmp
from tespy.components import characteristics as cmp_char


class connection:
    r"""
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
        State of the pure fluid on this connection: liquid ('l') or gaseous ('g').

    design : list
        List containing design parameters (stated as string).

    offdesign : list
        List containing offdesign parameters (stated as string).

    Note
    ----
    - The fluid balance parameter applies a balancing of the fluid vector on
      the specified conntion to 100 %. For example, you have four fluid
      components (a, b, c and d) in your vector, you set two of them
      (a and b) and want the other two (components c and d) to be a result of
      your calculation. If you set this parameter to True, the equation
      (0 = 1 - a - b - c - d) will be applied.

    - The specification of values for design and/or offdesign is used for
      automatic switch from design to offdesign calculation: All parameters given
      in 'design', e. g. :code:`design=['T', 'p']`, are unset in any offdesign
      calculation, parameters given in 'offdesign' are set for offdesign
      calculation.

    Example
    -------
    >>> from tespy import con, cmp, hlp
    >>> import numpy as np
    >>> so1 = cmp.source('source1')
    >>> so2 = cmp.source('source2')
    >>> si1 = cmp.sink('sink1')
    >>> si2 = cmp.sink('sink2')
    >>> so_si1 = con.connection(so1, 'out1', si1, 'in1')
    >>> so_si2 = con.connection(so2, 'out1', si2, 'in1')
    >>> so_si1.set_attr(v=0.012, m0=10, p=5, h=400, fluid={'H2O': 1, 'N2': 0})
    >>> p = hlp.dc_prop(val=50, val_set=True, unit='bar', unit_set=True)
    >>> so_si2.set_attr(m=con.ref(so_si1, 2, -5), p=p, h0=700, T=200, fluid={'N2': 1}, fluid_balance=True)
    >>> type(so_si1.p)
    <class 'tespy.tools.helpers.dc_prop'>
    >>> type(so_si1.fluid)
    <class 'tespy.tools.helpers.dc_flu'>
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
    >>> type(so_si2.m.ref.get_attr('obj'))
    <class 'tespy.connections.connection'>
    >>> so_si2.set_attr(Td_bp=5, T=np.nan)
    >>> type(so_si2.Td_bp)
    <class 'tespy.tools.helpers.dc_prop'>
    >>> so_si2.Td_bp.val
    5
    """

    def __init__(self, comp1, outlet_id, comp2, inlet_id, **kwargs):

        # check input parameters
        if not (isinstance(comp1, cmp.component) and
                isinstance(comp2, cmp.component)):
            msg = ('Error creating connection. Check if comp1, comp2 are of type component.')
            logging.error(msg)
            raise TypeError(msg)

        if comp1 == comp2:
            msg = ('Error creating connection. Can\'t connect component ' + comp1.label + ' to itself.')
            logging.error(msg)
            raise TESPyConnectionError(msg)

        if outlet_id not in comp1.outlets():
            msg = ('Error creating connection. Specified oulet_id (' + outlet_id + ') is not valid for component ' +
                   comp1.component() + '. Valid ids are: ' + str(comp1.outlets()) + '.')
            logging.error(msg)
            raise ValueError(msg)

        if inlet_id not in comp2.inlets():
            msg = ('Error creating connection. Specified inlet_id (' + inlet_id + ') is not valid for component ' +
                   comp2.component() + '. Valid ids are: ' + str(comp2.inlets()) + '.')
            logging.error(msg)
            raise ValueError(msg)

        # set specified values
        self.s = comp1
        self.s_id = outlet_id
        self.t = comp2
        self.t_id = inlet_id

        self.design = []
        self.offdesign = []

        # set default values for kwargs
        var = self.attr()

        for key in self.attr().keys():
            self.__dict__.update({key: var[key]})

        self.set_attr(**kwargs)

        msg = 'Created connection ' + self.s.label + ' (' + self.s_id + ') -> ' + self.t.label + ' (' + self.t_id + ').'
        logging.debug(msg)

    def set_attr(self, **kwargs):
        r"""
        Sets, resets or unsets attributes of a component for provided keyword arguments.

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

        x : float/tespy.connections.ref/tespy.tools.helpers.dc_prop
            Gas phase mass fraction specification.

        T : float/tespy.connections.ref/tespy.tools.helpers.dc_prop
            Temperature specification.

        Td_bp : float/tespy.tools.helpers.dc_prop
            Temperature difference to boiling point at pressure corresponding
            pressure of this connection in K.

        v : float/tespy.tools.helpers.dc_prop
            Volumetric flow specification.

        state : str
            State of the pure fluid on this connection: liquid ('l') or gaseous ('g').

        design : list
            List containing design parameters (stated as String).

        offdesign : list
            List containing offdesign parameters (stated as String).

        Note
        ----
        - The fluid balance parameter applies a balancing of the fluid vector on
          the specified conntion to 100 %. For example, you have four fluid
          components (a, b, c and d) in your vector, you set two of them
          (a and b) and want the other two (components c and d) to be a result of
          your calculation. If you set this parameter to True, the equation
          (0 = 1 - a - b - c - d) will be applied.

        - The specification of values for design and/or offdesign is used for
          automatic switch from design to offdesign calculation: All parameters given
          in 'design', e. g. :code:`design=['T', 'p']`, are unset in any offdesign
          calculation, parameters given in 'offdesign' are set for offdesign
          calculation.

        - The property state is applied on pure fluids only. If you specify the
          desired state of the fluid at a connection the convergence check will
          adjust the enthalpy values of that connection for the first iterations
          in order to meet the state requirement.
        """
        var = self.attr()
        var0 = [x + '0' for x in var.keys()]

        # set specified values
        for key in kwargs:
            if key in var.keys() or key in var0:
                if 'fluid' in key:
                    # fluid specification
                    if isinstance(kwargs[key], dict):
                        # starting values
                        if key in var0:
                            self.get_attr(key.replace('0', '')).set_attr(val0=kwargs[key])
                        # specified parameters
                        else:
                            self.get_attr(key).set_attr(val=kwargs[key].copy())
                            for f in kwargs[key]:
                                kwargs[key][f] = True
                            self.get_attr(key).set_attr(val_set=kwargs[key])

                    elif isinstance(kwargs[key], dc_flu):
                        # data container for fluids
                        self.__dict__.update({key: kwargs[key]})

                    else:
                        # bad datatype
                        msg = 'Bad datatype for connection keyword ' + key + '.'
                        logging.error(msg)
                        raise TypeError(msg)

                elif key == 'state':
                    if kwargs[key] in ['l', 'g']:
                        self.state.set_attr(val=kwargs[key], val_set=True)
                    elif isinstance(kwargs[key], dc_simple):
                        self.state = kwargs[key]
                    else:
                        if (isinstance(kwargs[key], float) or
                                isinstance(kwargs[key], np.float64) or
                                isinstance(kwargs[key], np.int64) or
                                isinstance(kwargs[key], int)):
                            if np.isnan(kwargs[key]):
                                self.state.set_attr(val=kwargs[key], val_set=False)
                            else:
                                msg = 'Datatype for keyword argument ' + str(key) + ' must be str.'
                                logging.error(msg)
                                raise TypeError(msg)
                        else:
                            msg = 'Keyword argument ' + str(key) + ' must be \'l\' or \'g\'.'
                            logging.error(msg)
                            raise ValueError(msg)


                elif (isinstance(kwargs[key], float) or
                        isinstance(kwargs[key], np.float64) or
                        isinstance(kwargs[key], np.int64) or
                        isinstance(kwargs[key], int)):
                    # unset
                    if np.isnan(kwargs[key]) and key not in var0:
                        self.get_attr(key).set_attr(val_set=False, ref_set=False)
                    # starting value
                    elif key in var0:
                        self.get_attr(key.replace('0', '')).set_attr(val0=kwargs[key])
                    # set/reset
                    else:
                        self.get_attr(key).set_attr(val_set=True, val=kwargs[key], val0=kwargs[key])

                # reference object
                elif isinstance(kwargs[key], ref):
                    if key == 'x' or key == 'v' or key == 'Td_bp':
                        msg = 'References for volumetric flow, vapour mass fraction and subcooling/superheating not implemented.'
                        logging.warning(msg)
                    else:
                        self.get_attr(key).set_attr(ref=kwargs[key])
                        self.get_attr(key).set_attr(ref_set=True)

                # data container specification
                elif isinstance(kwargs[key], dc_prop):
                    self.__dict__.update({key: kwargs[key]})

                # invalid datatype for keyword
                else:
                    msg = 'Bad datatype for keyword argument ' + str(key) + '.'
                    logging.error(msg)
                    raise TypeError(msg)

            # fluid balance
            elif key == 'fluid_balance':
                if isinstance(kwargs[key], bool):
                    self.get_attr('fluid').set_attr(balance=kwargs[key])
                else:
                    msg = ('Datatype for keyword argument ' + str(key) + ' must be boolean.')
                    logging.error(msg)
                    raise TypeError(msg)

            elif key == 'design' or key == 'offdesign':
                if not isinstance(kwargs[key], list):
                    msg = 'Please provide the ' + key + ' parameters as list!'
                    logging.error(msg)
                    raise TypeError(msg)
                if set(kwargs[key]).issubset(var.keys()):
                    self.__dict__.update({key: kwargs[key]})
                else:
                    msg = ('Available parameters for (off-)design specification are: ' + str(var.keys()) + '.')
                    logging.error(msg)
                    raise ValueError(msg)

            # invalid keyword
            else:
                msg = 'Connection has no attribute ' + str(key) + '.'
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

    def attr(self):
        r"""
        Returns the list of available attributes of a connection.

        Returns
        -------
        out : list
            List of available attributes of a connection.
        """
        return {'m': dc_prop(), 'p': dc_prop(), 'h': dc_prop(), 'T': dc_prop(),
                'x': dc_prop(), 'v': dc_prop(),
                'fluid': dc_flu(), 'Td_bp': dc_prop(), 'state': dc_simple()}

    def to_flow(self):
        r"""
        Returns a list with the SI-values for the network variables (m, p, h, fluid vector) of a connection.

        Returns
        -------
        out : list
            List of mass flow and fluid property information.
        """
        return [self.m.val_SI, self.p.val_SI, self.h.val_SI, self.fluid.val]

    def to_flow_design(self):
        r"""
        Returns a list with the SI-values for the network variables (m, p, h, fluid vector) of a connection at design point.

        Returns
        -------
        out : list
            List of mass flow and fluid property information.
        """
        return [self.m.design, self.p.design, self.h.design, self.fluid.design]


class bus:
    r"""
    A bus is used to connect different energy flows (power, heat flow, thermal input).

    Parameters
    ----------
    label : str
        Label for the bus.

    P : float
        Total power/heat flow specification for bus, :math:`P \text{/W}`.

    Example
    -------
    >>> from tespy import cmp, con, nwk, cmp_char
    >>> import shutil
    >>> import numpy as np
    >>> fluid_list = ['Ar', 'N2', 'O2', 'CO2', 'CH4', 'H2O']
    >>> nw = nwk.network(fluids=fluid_list, p_unit='bar', T_unit='C', p_range=[0.5, 10], T_range=[10, 1200])
    >>> nw.set_printoptions(print_level='none')

    >>> amb = cmp.source('ambient')
    >>> sf = cmp.source('fuel')
    >>> fg = cmp.sink('flue gas outlet')
    >>> cw_in1 = cmp.source('cooling water inlet1')
    >>> cw_in2 = cmp.source('cooling water inlet2')
    >>> cw_out1 = cmp.sink('cooling water outlet1')
    >>> cw_out2 = cmp.sink('cooling water outlet2')
    >>> split = cmp.splitter('splitter')
    >>> merge = cmp.merge('merge')
    >>> chp = cmp.cogeneration_unit(label='cogeneration unit')
    >>> amb_comb = con.connection(amb, 'out1', chp, 'in3')
    >>> sf_comb = con.connection(sf, 'out1', chp, 'in4')
    >>> comb_fg = con.connection(chp, 'out3', fg, 'in1')
    >>> nw.add_conns(sf_comb, amb_comb, comb_fg)
    >>> cw1_chp1 = con.connection(cw_in1, 'out1', chp, 'in1')
    >>> cw2_chp2 = con.connection(cw_in2, 'out1', chp, 'in2')
    >>> nw.add_conns(cw1_chp1, cw2_chp2)
    >>> chp1_cw = con.connection(chp, 'out1', cw_out1, 'in1')
    >>> chp2_cw = con.connection(chp, 'out2', cw_out2, 'in1')
    >>> nw.add_conns(chp1_cw, chp2_cw)
    >>> chp.set_attr(fuel='CH4', pr1=0.99, pr2=0.99, P=1e6, lamb=1.2)
    >>> amb_comb.set_attr(p=5, T=30, fluid={'Ar': 0.0129, 'N2': 0.7553, 'H2O': 0, 'CH4': 0, 'CO2': 0.0004, 'O2': 0.2314})
    >>> sf_comb.set_attr(T=30, fluid={'CO2': 0, 'Ar': 0, 'N2': 0, 'O2': 0, 'H2O': 0, 'CH4': 1})
    >>> cw1_chp1.set_attr(p=3, T=60, m=50, fluid={'CO2': 0, 'Ar': 0, 'N2': 0, 'O2': 0, 'H2O': 1, 'CH4': 0})
    >>> cw2_chp2.set_attr(p=3, T=80, m=50, fluid={'CO2': 0, 'Ar': 0, 'N2': 0, 'O2': 0, 'H2O': 1, 'CH4': 0})

    >>> si = cmp.sink('sink')
    >>> so = cmp.source('source')
    >>> t = cmp.turbine('turbine')
    >>> c = cmp.heat_exchanger_simple('condenser')
    >>> inc = con.connection(so, 'out1', t, 'in1')
    >>> ws = con.connection(t, 'out1', c, 'in1')
    >>> outg = con.connection(c, 'out1', si, 'in1')
    >>> nw.add_conns(inc, ws, outg)
    >>> c.set_attr(pr=1)
    >>> t.set_attr(eta_s=0.8, design=['eta_s'], offdesign=['eta_s_char', 'cone'])
    >>> inc.set_attr(fluid={'CO2': 0, 'Ar': 0, 'N2': 0, 'O2': 0, 'H2O': 1, 'CH4': 0}, T=600, p=50, design=['p'])
    >>> outg.set_attr(p=0.5, x=0)

    >>> power_bus = con.bus('power')
    >>> heat_bus = con.bus('heat')
    >>> power_bus.set_attr(P=-2e6)
    >>> load = np.array([0.2, 0.4, 0.6, 0.8, 1, 1.2])
    >>> efficiency = np.array([0.8, 0.9, 0.93, 0.95, 0.955, 0.94])
    >>> t_gen = cmp_char.characteristics(x=load, y=efficiency)
    >>> cog_gen = cmp_char.characteristics(x=load, y=-efficiency)
    >>> power_bus.add_comps({'c': t, 'char': t_gen}, {'c': chp, 'char': cog_gen, 'p': 'P'})
    >>> heat_bus.add_comps({'c': c, 'char': -1}, {'c': chp, 'p': 'Q'})
    >>> nw.add_busses(power_bus, heat_bus)
    >>> type(heat_bus)
    <class 'tespy.connections.bus'>
    >>> type(t_gen)
    <class 'tespy.components.characteristics.characteristics'>
    >>> all(power_bus.comps.loc[chp]['char'].x == load)
    True
    >>> all(power_bus.comps.loc[chp]['char'].y == -efficiency)
    True
    >>> heat_bus.comps.loc[c]['char'].y
    array([-1, -1, -1, -1])

    >>> nw.solve('design')
    >>> nw.save('tmp')
    >>> chp.set_attr(P=np.nan)
    >>> t.set_attr(P=-2e6)
    >>> power_bus.set_attr(P=-2.5e6)
    >>> nw.solve('offdesign', design_path='tmp')
    >>> round(chp.P.val, 0)
    662235.0
    >>> round(heat_bus.P.val, 0)
    5366642.0
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def __init__(self, label, **kwargs):

        self.comps = pd.DataFrame(columns=['param', 'P_ref', 'char'])

        self.label = label
        self.P = dc_cp(val=np.nan, val_set=False)
        self.char = cmp_char.characteristics(x=np.array([0, 1, 2, 3]),
                                             y=np.array([1, 1, 1, 1]))

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

        Note
        ----
        Specify :math:`P=\text{nan}`, if you want to unset the value.
        """
        self.label = kwargs.get('label', self.label)
        self.P.val = kwargs.get('P', self.P.val)

        if np.isnan(self.P.val):
            self.P.val_set = False
        else:
            self.P.val_set = True

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
        c : dict
            Dictionary containing the component information to be added to the bus.
            These information are described in the notes!

        Note
        ----
        Keys for the dictionary c:

        - c (tespy.components.components.component): Component you want to add to the bus.
        - p (str): Bus parameter, optional.

            - You do not need to provide a parameter, if the component only has one
              option for the bus (turbomachines, heat exchangers, combustion
              chamber).
            - For instance, you do neet do provide a parameter, if you want to add
              a cogeneration unit ('Q', 'Q1', 'Q2', 'TI', 'P', 'Qloss').

        - char (float/tespy.components.characteristics.characteristics):
          Characteristic function for this components share to the bus value, optional.

            - If you do not provide a characteristic line at all, TESPy assumes a
              constant factor of 1.
            - If you provide a numeric value instead of a characteristic line,
              TESPy takes this numeric value as a constant factor.
            - Provide a TESPy.characteristic (cmp_char), if you want the factor
              to follow a characteristic line.

        - P_ref (float): Energy flow specification for reference case, :math:`P \text{/W}`, optional.
        """
        for c in args:
            if isinstance(c, dict):
                if 'c' in c.keys():
                    if isinstance(c['c'], cmp.component):
                        self.comps.loc[c['c']] = [None, np.nan, self.char]
                    else:
                        msg = ('Keyword c must hold a TESPy component.')
                        logging.error(msg)
                        raise TypeError(msg)
                else:
                        msg = ('You must provide the component c.')
                        logging.error(msg)
                        raise TypeError(msg)

                for k, v in c.items():
                    if k == 'p':
                        if isinstance(v, str) or v is None:
                            self.comps.loc[c['c']]['param'] = v
                        else:
                            msg = ('Parameter p must be a string.')
                            logging.error(msg)
                            raise TypeError(msg)

                    elif k == 'char':
                        if isinstance(v, cmp_char.characteristics):
                            self.comps.loc[c['c']]['char'] = v
                        elif (isinstance(v, float) or
                              isinstance(v, np.float64) or
                              isinstance(v, np.int64) or
                              isinstance(v, int)):
                            x = np.array([0, 1, 2, 3])
                            y = np.array([1, 1, 1, 1]) * v
                            self.comps.loc[c['c']]['char'] = (cmp_char.characteristics(x=x, y=y))
                        else:
                            msg = ('Char must be a number or a TESPy characteristics.')
                            logging.error(msg)
                            raise TypeError(msg)

                    elif k == 'P_ref':
                        if (v is None or isinstance(v, float) or
                                isinstance(v, np.float64) or
                                isinstance(v, np.int64) or
                                isinstance(v, int)):
                            self.comps.loc[c['c']]['P_ref'] = v
                        else:
                            msg = ('Reference value must be numeric.')
                            logging.error(msg)
                            raise TypeError(msg)
            else:
                msg = ('Provide arguments as dicts. See the documentation of bus.add_comps() for more information.')
                logging.error(msg)
                raise TESPyConnectionError(msg)

            msg = 'Added component ' + c['c'].label + ' to bus ' + self.label + '.'
            logging.debug(msg)


class ref:
    r"""
    Reference class to reference fluid properties from one connection to another connection.

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
    Reference the mass flow of one connection :math:`\dot{m}` to another mass flow
    :math:`\dot{m}_{ref}`

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

        msg = ('Created reference object with factor ' + str(self.f) + ' and delta ' + str(self.d) + ' referring to connection ' +
               ref_obj.s.label + ' (' + ref_obj.s_id + ') -> ' + ref_obj.t.label + ' (' + ref_obj.t_id + ').')
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
