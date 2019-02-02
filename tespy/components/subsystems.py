# -*- coding: utf-8

"""
.. module:: subsystems
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import numpy as np
import logging

from tespy import networks as nwk
from tespy.connections import connection, ref
from tespy.components import components as cmp


class subsystem:
    r"""
    Class subsystem is the base class of all TESPy subsystems.

    Parameters
    ----------
    label : str
        The label of the subsystem.

    **kwargs :
        See the class documentation of desired subsystem for available keywords.
        If you want to define your own subsystem, provide the keywords in the :func:`tespy.components.subsystems.subystem.attr` method.

    Note
    ----
    The initialisation method (__init__), setter method (set_attr) and getter method (get_attr)
    are used for instances of class subsystem and its children.

    Keywords for keyword arguments depend on the type of subsystem you want to create/use.

    Example
    -------
    Basic example for a setting up a tespy.components.subsystems.subsystem object.
    This example does not run a tespy calculation!

    >>> from tespy import subsys
    >>> mysub = subsys.subsystem('mySubsystem')
    >>> type(mysub)
    <class 'tespy.components.subsystems.subsystem'>
    >>> mysub.get_attr('label')
    'mySubsystem'
    >>> type(mysub.get_network())
    <class 'tespy.networks.network'>
    """
    def __init__(self, label, **kwargs):

        if not isinstance(label, str):
            msg = 'Subsystem label must be of type str!'
            logging.error(msg)
            raise ValueError(msg)

        elif len([x for x in [';', ', ', '.'] if x in label]) > 0:
            msg = 'Can\'t use ' + str([';', ', ', '.']) + ' in label.'
            logging.error(msg)
            raise ValueError(msg)
        else:
            self.label = label

        # set default values
        for key in self.attr():
            self.__dict__.update({key: np.nan})
            self.__dict__.update({key + '_set': False})

        self.subsys_init()
        self.set_attr(**kwargs)

    def set_attr(self, **kwargs):
        r"""
        Sets, resets or unsets attributes of a component for provided keyword arguments.

        Parameters
        ----------
        label : str
            The label of the subsystem.

        **kwargs :
            See the class documentation of desired subsystem for available keywords.
            If you want to define your own subsystem, provide the keywords in the :func:`tespy.components.subsystems.subystem.attr` method.

        Note
        ----
        Allowed keywords in kwargs are obtained from class documentation as all
        subsystems share the :func:`tespy.components.subsystems.subsystem.set_attr` method.
        """
        # set provided values,  check for invalid keys
        for key in kwargs:
            if key in self.attr():
                if (isinstance(kwargs[key], float) or
                        isinstance(kwargs[key], np.float64) or
                        isinstance(kwargs[key], int)):
                    if np.isnan(kwargs[key]):
                        self.__dict__.update({key: np.nan})
                        self.__dict__.update({key + '_set': False})
                    else:
                        self.__dict__.update({key: kwargs[key]})
                        self.__dict__.update({key + '_set': True})

                elif isinstance(kwargs[key], str):
                    self.__dict__.update({key: kwargs[key]})
            else:
                msg = ('Component ' + self.label + ' has no attribute ' + str(key))
                logging.error(msg)
                raise KeyError(msg)

        self.set_comps()
        self.set_conns()

        msg = 'Created subsystem ' + self.label + '.'
        logging.debug(msg)

    def get_attr(self, key):
        r"""
        Get the value of a subsystem's attribute.

        Parameters
        ----------
        key : str
            The attribute you want to retrieve.

        Returns
        -------
        out :
            Value of specified attribute.
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            msg = 'Subsystem ' + self.label + ' has no attribute ' + key + '.'
            logging.error(msg)
            raise KeyError(msg)

    def attr(self):
        return ['num_i', 'num_o']

    def subsys_init(self):
        r"""
        Creates network of the subsystem on subsystem creation.
        """
        self.create_comps()
        self.create_conns()

    def create_comps(self):
        return

    def create_conns(self):
        self.conns = []

    def set_comps(self):
        return

    def set_conns(self):
        if not hasattr(self, 'nw'):
            self.create_network()
        return

    def create_network(self):
        self.nw = nwk.network(fluids=[])
        for c in self.conns:
            self.nw.add_conns(c)

    def get_network(self):
        return self.nw

# %%


class dr_eva_forced(subsystem):
    r"""
    Drum with evaporator in forced flow.

    Inlets/Outlets

        - in1, in2
        - out1, out2

    Image

        .. image:: _images/subsys_dr_eva_forced.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the subsystem.

    pr1_eva : float
        Pressure drop at hot side of evaporator.

    pr2_eva : float
        Pressure drop at cold side of evaporator.

    eta_s : float
        Circulation pump isentropic efficiency.

    PP : float
        Pinch point of evaporator.

    circ_num : float
        Circulation number, ratio of mass flow through cold side of
        evaporator to mass flow at the drum's inlet.

    Example
    -------
    >>> from tespy import subsys
    >>> import numpy as np
    >>> mysub = subsys.dr_eva_forced('drum with evaporator')
    >>> mysub.set_attr(pr1_eva=0.99, pr2_eva=0.98, PP=20, circ_num=3.4)
    >>> mysub.get_attr('pr1_eva')
    0.99
    >>> mysub.evaporator.label
    'drum with evaporator_evaporator'
    >>> mysub.inlet.num_inter.val
    2
    >>> mysub.pump.eta_s.is_set
    False
    >>> mysub.evaporator.ttd_l.val == mysub.PP
    True
    >>> mysub.conns[3].m.ref_set
    True
    >>> mysub.set_attr(circ_num=np.nan)
    >>> mysub.conns[3].m.ref_set
    False
    """

    def attr(self):
        return ([n for n in subsystem.attr(self) if
                 n != 'num_i' and n != 'num_o'] +
                ['pr1_eva', 'pr2_eva', 'eta_s', 'PP', 'circ_num'])

    def create_comps(self):

        self.num_i = 2
        self.num_o = 2
        self.inlet = cmp.subsys_interface(label=self.label + '_inlet',
                                          num_inter=self.num_i)
        self.outlet = cmp.subsys_interface(label=self.label + '_outlet',
                                           num_inter=self.num_o)
        self.pump = cmp.pump(label=self.label + '_pump')
        self.drum = cmp.drum(label=self.label + '_drum')
        self.evaporator = cmp.heat_exchanger(label=self.label + '_evaporator')

    def set_comps(self):

        self.evaporator.set_attr(ttd_l=self.PP)
        self.evaporator.set_attr(pr1=self.pr1_eva)
        self.evaporator.set_attr(pr2=self.pr2_eva)
        self.pump.set_attr(eta_s=self.eta_s)

    def create_conns(self):

        self.conns = []

        self.conns += [connection(self.inlet, 'out1', self.evaporator, 'in1')]
        self.conns += [connection(self.evaporator, 'out1', self.outlet, 'in1')]
        self.conns += [connection(self.inlet, 'out2', self.drum, 'in1')]
        self.conns += [connection(self.drum, 'out1', self.pump, 'in1')]
        self.conns += [connection(self.pump, 'out1', self.evaporator, 'in2')]
        self.conns += [connection(self.evaporator, 'out2', self.drum, 'in2')]
        self.conns += [connection(self.drum, 'out2', self.outlet, 'in2')]

    def set_conns(self):

        if self.circ_num_set:
            self.conns[3].set_attr(m=ref(self.conns[-1], self.circ_num, 0))
        else:
            self.conns[3].set_attr(m=np.nan)

# %%


class dr_eva_natural(subsystem):
    r"""
    Drum with evaporator in natural flow.

    Inlets/Outlets

        - in1, in2
        - out1, out2

    Image

        .. image:: _images/subsys_dr_eva_natural.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the subsystem.

    pr1_eva : float
        Pressure drop at hot side of evaporator.

    PP : float
        Pinch point of evaporator.

    circ_num : float
        Circulation number, ratio of mass flow through cold side of
        evaporator to mass flow at the drum's inlet.

    Example
    -------
    >>> from tespy import subsys
    >>> mysub = subsys.dr_eva_natural('drum with evaporator')
    >>> mysub.set_attr(pr1_eva=0.99, PP=20, circ_num=3.4)
    >>> mysub.get_attr('PP')
    20
    >>> mysub.drum.label
    'drum with evaporator_drum'
    """

    def attr(self):
        return ([n for n in subsystem.attr(self) if
                 n != 'num_i' and n != 'num_o'] +
                ['pr1_eva', 'PP', 'circ_num'])

    def create_comps(self):

        self.num_i = 2
        self.num_o = 2
        self.inlet = cmp.subsys_interface(label=self.label + '_inlet',
                                          num_inter=self.num_i)
        self.outlet = cmp.subsys_interface(label=self.label + '_outlet',
                                           num_inter=self.num_o)
        self.drum = cmp.drum(label=self.label + '_drum')
        self.evaporator = cmp.heat_exchanger(label=self.label + '_evaporator',
                                             mode='man')

    def set_comps(self):

        self.evaporator.set_attr(ttd_l=self.PP)
        self.evaporator.set_attr(pr1=self.pr1_eva)

    def create_conns(self):

        self.conns = []

        self.conns += [connection(self.inlet, 'out1', self.evaporator, 'in1')]
        self.conns += [connection(self.evaporator, 'out1', self.outlet, 'in1')]
        self.conns += [connection(self.inlet, 'out2', self.drum, 'in1')]
        self.conns += [connection(self.drum, 'out1', self.evaporator, 'in2')]
        self.conns += [connection(self.evaporator, 'out2', self.drum, 'in2')]
        self.conns += [connection(self.drum, 'out2', self.outlet, 'in2')]

    def set_conns(self):

        if self.circ_num_set:
            self.conns[3].set_attr(m=ref(self.conns[-1], self.circ_num, 0))
        else:
            self.conns[3].set_attr(m=np.nan)

# %%


class ph_desup_cond(subsystem):
    r"""
    Preheater with desuperheater, condesate at hot side outlet.

    Inlets/Outlets

        - in1, in2
        - out1, out2

    Image

        .. image:: _images/subsys_ph_desup_cond.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the subsystem.

    pr1_desup : float
        Pressure drop at hot side of desuperheater.

    pr2_desup : float
        Pressure drop at cold side of desuperheater.

    pr1_cond : float
        Pressure drop at hot side of condensator.

    pr2_cond : float
        Pressure drop at cold side of condensator.

    ttd : float
        Upper terminal temperature difference at condenser.

    Example
    -------
    >>> from tespy import subsys
    >>> mysub = subsys.ph_desup_cond('preheater')
    """

    def attr(self):
        return ([n for n in subsystem.attr(self) if
                 n != 'num_i' and n != 'num_o'] +
                ['ttd', 'pr1_desup', 'pr2_desup', 'pr1_cond', 'pr2_cond'])

    def create_comps(self):

        self.num_i = 2
        self.num_o = 2
        self.inlet = cmp.subsys_interface(label=self.label + '_inlet',
                                          num_inter=self.num_i)
        self.outlet = cmp.subsys_interface(label=self.label + '_outlet',
                                           num_inter=self.num_o)
        self.desup = cmp.desuperheater(label=self.label + '_desup')
        self.condenser = cmp.condenser(label=self.label + '_condenser')
        self.valve = cmp.valve(label=self.label + '_valve', mode='man')

    def set_comps(self):

        self.condenser.set_attr(ttd_u=self.ttd)
        self.desup.set_attr(pr1=self.pr1_desup)
        self.desup.set_attr(pr2=self.pr2_desup)
        self.condenser.set_attr(pr1=self.pr1_cond)
        self.condenser.set_attr(pr2=self.pr2_cond)

    def create_conns(self):

        self.conns = []

        self.conns += [connection(self.inlet, 'out1', self.desup, 'in1')]
        self.conns += [connection(self.desup, 'out1', self.condenser, 'in1')]
        self.conns += [connection(self.condenser, 'out1', self.valve, 'in1')]
        self.conns += [connection(self.valve, 'out1', self.outlet, 'in1')]
        self.conns += [connection(self.inlet, 'out2', self.condenser, 'in2')]
        self.conns += [connection(self.condenser, 'out2', self.desup, 'in2')]
        self.conns += [connection(self.desup, 'out2', self.outlet, 'in2')]

# %%


class ph_desup_cond_subc(subsystem):
    r"""
    Preheater with desuperheater and subcooler, subcooled liquid at hot side outlet.

    Inlets/Outlets

        - in1, in2
        - out1, out2

    Image

        .. image:: _images/subsys_ph_desup_cond_subc.svg
           :scale: 100 %
           :alt: alternative text
           :align: center

    Parameters
    ----------
    label : str
        The label of the subsystem.

    pr1_desup : float
        Pressure drop at hot side of desuperheater.

    pr2_desup : float
        Pressure drop at cold side of desuperheater.

    pr1_cond : float
        Pressure drop at hot side of condensator.

    pr2_cond : float
        Pressure drop at cold side of condensator.

    pr1_subc : float
        Pressure drop at hot side of subcooler.

    pr2_subc : float
        Pressure drop at cold side of subcooler.

    ttd : float
        Upper terminal temperature difference at condenser.

    pr_v : float
        Pressure ratio of the valve at hot side outlet.

    Example
    -------
    >>> from tespy import subsys
    >>> mysub = subsys.ph_desup_cond_subc('preheater')
    """

    def attr(self):
        return ([n for n in subsystem.attr(self) if
                 n != 'num_i' and n != 'num_o'] +
                ['ttd', 'pr1_desup', 'pr2_desup',
                 'pr1_cond', 'pr2_cond', 'pr1_subc', 'pr2_subc', 'pr_v'])

    def create_comps(self):

        self.num_i = 2
        self.num_o = 2
        self.inlet = cmp.subsys_interface(label=self.label + '_inlet',
                                          num_inter=self.num_i)
        self.outlet = cmp.subsys_interface(label=self.label + '_outlet',
                                           num_inter=self.num_o)
        self.desup = cmp.desuperheater(label=self.label + '_desup')
        self.condenser = cmp.condenser(label=self.label + '_condenser')
        self.subcooler = cmp.heat_exchanger(label=self.label + '_subcooler')
        self.valve = cmp.valve(label=self.label + '_valve', mode='man')

    def set_comps(self):

        self.condenser.set_attr(ttd_u=self.ttd)
        self.desup.set_attr(pr1=self.pr1_desup)
        self.desup.set_attr(pr2=self.pr2_desup)
        self.condenser.set_attr(pr1=self.pr1_cond)
        self.condenser.set_attr(pr2=self.pr2_cond)
        self.subcooler.set_attr(pr1=self.pr1_subc)
        self.subcooler.set_attr(pr2=self.pr2_subc)
        self.valve.set_attr(pr=self.pr_v)

    def create_conns(self):

        self.conns = []

        self.conns += [connection(self.inlet, 'out1', self.desup, 'in1')]
        self.conns += [connection(self.desup, 'out1', self.condenser, 'in1')]
        self.conns += [connection(self.condenser, 'out1',
                                  self.subcooler, 'in1')]
        self.conns += [connection(self.subcooler, 'out1', self.valve, 'in1')]
        self.conns += [connection(self.valve, 'out1', self.outlet, 'in1')]
        self.conns += [connection(self.inlet, 'out2', self.subcooler, 'in2')]
        self.conns += [connection(self.subcooler, 'out2',
                                  self.condenser, 'in2')]
        self.conns += [connection(self.condenser, 'out2', self.desup, 'in2')]
        self.conns += [connection(self.desup, 'out2', self.outlet, 'in2')]
