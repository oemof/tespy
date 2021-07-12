# -*- coding: utf-8

"""Module of class Bus.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/connections/bus.py

SPDX-License-Identifier: MIT
"""

import logging

import numpy as np
import pandas as pd

from tespy.components.component import Component
from tespy.tools.characteristics import CharLine
from tespy.tools.data_containers import DataContainerSimple as dc_simple

# pass the warning messages to the logger
logging.captureWarnings(True)


class Bus:
    r"""
    A bus is used to connect different energy flows.

    Parameters
    ----------
    label : str
        Label for the bus.

    P : float
        Total power/heat flow specification for bus, :math:`P \text{/W}`.

    printout : boolean
        Print the results of this bus to prompt with the
        :py:meth:`tespy.networks.network.Network.print_results` method.
        Standard value is :code:`True`.

    Example
    -------
    Busses are used to connect energy flow of different components. They can
    also be used to introduce efficiencies of energy conversion, e.g. in
    motors, generator or boilers. This example takes the combustion engine
    example at
    :py:class:`tespy.components.combustion.combustion_engine.CombustionEngine`
    and adds a flue gas cooler and a circulation pump for the cooling water.
    Then busses for heat output, thermal input and electricity output are
    implemented.

    >>> from tespy.components import (Sink, Source, CombustionEngine,
    ... HeatExchanger, Merge, Splitter, Pump)
    >>> from tespy.connections import Connection, Ref, Bus
    >>> from tespy.networks import Network
    >>> from tespy.tools import CharLine
    >>> import numpy as np
    >>> import shutil
    >>> fluid_list = ['Ar', 'N2', 'O2', 'CO2', 'CH4', 'H2O']
    >>> nw = Network(fluids=fluid_list, p_unit='bar', T_unit='C',
    ... p_range=[0.5, 10], iterinfo=False)
    >>> amb = Source('ambient')
    >>> sf = Source('fuel')
    >>> fg = Sink('flue gas outlet')
    >>> cw_in = Source('cooling water inlet')
    >>> sp = Splitter('cooling water splitter', num_out=2)
    >>> me = Merge('cooling water merge', num_in=2)
    >>> cw_out = Sink('cooling water outlet')
    >>> fgc = HeatExchanger('flue gas cooler')
    >>> pu = Pump('cooling water pump')
    >>> chp = CombustionEngine(label='internal combustion engine')
    >>> amb_comb = Connection(amb, 'out1', chp, 'in3')
    >>> sf_comb = Connection(sf, 'out1', chp, 'in4')
    >>> comb_fgc = Connection(chp, 'out3', fgc, 'in1')
    >>> fgc_fg = Connection(fgc, 'out1', fg, 'in1')
    >>> nw.add_conns(sf_comb, amb_comb, comb_fgc, fgc_fg)
    >>> cw_pu = Connection(cw_in, 'out1', pu, 'in1')
    >>> pu_sp = Connection(pu, 'out1', sp, 'in1')
    >>> sp_chp1 = Connection(sp, 'out1', chp, 'in1')
    >>> sp_chp2 = Connection(sp, 'out2', chp, 'in2')
    >>> chp1_me = Connection(chp, 'out1', me, 'in1')
    >>> chp2_me = Connection(chp, 'out2', me, 'in2')
    >>> me_fgc = Connection(me, 'out1', fgc, 'in2')
    >>> fgc_cw = Connection(fgc, 'out2', cw_out, 'in1')
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
    >>> sp_chp2.set_attr(m=Ref(sp_chp1, 1, 0))

    Cooling water mass flow is calculated given the feed water temperature
    (90 °C). The pressure at the cooling water outlet should be identical to
    pressure before pump. The flue gases of the combustion engine leave the
    flue gas cooler at 120 °C.

    >>> fgc_cw.set_attr(p=Ref(cw_pu, 1, 0), T=90)
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
    >>> eff = np.array([0.9, 0.94, 0.97, 0.99, 1, 0.99]) * 0.98
    >>> gen = CharLine(x=load, y=eff)
    >>> mot = CharLine(x=load, y=eff)
    >>> power_bus = Bus('total power output', P=-10e6)
    >>> heat_bus = Bus('total heat input')
    >>> fuel_bus = Bus('thermal input')

    You can check, if the bus value is set or not identically to connections.
    Unsetting is possible using :code:`np.nan` or :code:`None`. For
    demonstration we will specify a value for the heat bus and unset it again.

    >>> heat_bus.P.is_set
    False
    >>> heat_bus.set_attr(P=-1e5)
    >>> heat_bus.P.is_set
    True
    >>> heat_bus.set_attr(P=None)
    >>> heat_bus.P.is_set
    False

    >>> power_bus.add_comps({'comp': chp, 'char': gen, 'param': 'P'},
    ... {'comp': pu, 'char': mot, 'base': 'bus'})
    >>> heat_bus.add_comps({'comp': chp, 'param': 'Q', 'char': -1},
    ... {'comp': fgc, 'char': -1})
    >>> fuel_bus.add_comps({'comp': chp, 'param': 'TI'},
    ... {'comp': pu, 'char': mot})
    >>> nw.add_busses(power_bus, heat_bus, fuel_bus)
    >>> mode = 'design'
    >>> nw.solve(mode=mode)
    >>> nw.save('tmp')

    The heat bus characteristic for the combustion engine and the flue gas
    cooler have automatically been transformed into an array. The total heat
    output can be seen on both, the heat bus and in the enthalpy rise of the
    cooling water in the combustion engine and the flue gas cooler. Therefore
    these values must be identical.

    >>> heat_bus.comps.loc[fgc]['char'].x
    array([0., 3.])
    >>> heat_bus.comps.loc[fgc]['char'].y
    array([-1., -1.])
    >>> round(chp.ti.val, 0)
    25819387.0
    >>> round(chp.Q1.val + chp.Q2.val, 0)
    -8899014.0
    >>> round(fgc_cw.m.val_SI * (fgc_cw.h.val_SI - pu_sp.h.val_SI), 0)
    12476723.0
    >>> round(heat_bus.P.val, 0)
    12476723.0
    >>> round(pu.calc_bus_efficiency(power_bus), 2)
    0.98
    >>> power_bus.set_attr(P=-7.5e6)
    >>> mode = 'offdesign'
    >>> nw.solve(mode=mode, design_path='tmp', init_path='tmp')
    >>> round(chp.ti.val, 0)
    21192700.0
    >>> round(chp.P.val / chp.P.design, 3)
    0.761
    >>> round(pu.calc_bus_efficiency(power_bus), 3)
    0.968
    >>> shutil.rmtree('./tmp', ignore_errors=True)
    """

    def __init__(self, label, **kwargs):

        self.comps = pd.DataFrame(
            columns=['param', 'P_ref', 'char', 'efficiency', 'base'],
            dtype='object')

        self.label = label
        self.P = dc_simple(val=np.nan, is_set=False)
        self.char = CharLine(x=np.array([0, 3]), y=np.array([1, 1]))
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

        printout : boolean
            Print the results of this bus to prompt with the
            :py:meth:`tespy.networks.network.Network.print_results` method.
            Standard value is :code:`True`.

        Note
        ----
        Specify :math:`P=\text{nan}`, if you want to unset the value of P.
        """
        for key in kwargs:
            try:
                float(kwargs[key])
                is_numeric = True
            except (TypeError, ValueError):
                is_numeric = False
            if key == 'P':
                if is_numeric:
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

        - comp (tespy.components.component.Component): Component you want to
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
            - Provide a :py:class:`tespy.tools.characteristics.CharLine`, if
              you want the factor to follow a characteristic line.

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
                if 'comp' in c:
                    comp = c['comp']
                    # default values
                    if isinstance(comp, Component):
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
                    if k == 'param':
                        if isinstance(v, str) or v is None:
                            self.comps.loc[comp, 'param'] = v
                        else:
                            msg = (
                                'The bus parameter selection must be a '
                                'string (at bus ' + self.label + ').')
                            logging.error(msg)
                            raise TypeError(msg)

                    elif k == 'char':
                        try:
                            float(v)
                            is_numeric = True
                        except (TypeError, ValueError):
                            is_numeric = False
                        if isinstance(v, CharLine):
                            self.comps.loc[comp, 'char'] = v
                        elif is_numeric:
                            x = np.array([0, 3])
                            y = np.array([1, 1]) * v
                            self.comps.loc[comp, 'char'] = (
                                    CharLine(x=x, y=y))
                        else:
                            msg = (
                                'Char must be a number or a TESPy '
                                'characteristics char line.')
                            logging.error(msg)
                            raise TypeError(msg)

                    elif k == 'P_ref':
                        try:
                            float(v)
                            is_numeric = True
                        except (TypeError, ValueError):
                            is_numeric = False
                        if v is None or is_numeric:
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
