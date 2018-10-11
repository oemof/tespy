"""
.. module:: components
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import numpy as np

from tespy.helpers import (
    MyComponentError, MyNetworkError
)

from tespy import networks as nwk
from tespy.connections import connection, ref
from tespy.components import components as cmp


class subsystem:
    r"""

    :param label: label for subsystem
    :type label: str
    :param kwargs: for the keyword arguments see :code:`subsystem.attr()`
    :returns: no return value
    :raises: - :code:`TypeError`, if label is not of type str
               components
             - :code:`ValueError`, if label contains forbidden characters
               (';', ',', '.')

    **example**

    .. code-block:: python

        ph_desup = ph_desup_cond('preheater 1', ttd=5)

    creates a subsystem of a desuperheater and a condenser with a vessel at
    the hot outlet labeled "preheater 1" and sets the
    terminal temperature difference (hot side inlet at condenser to
    cold side outlet after condenser) to 5 K

    initialisation method is used for instances of class component and its
    children

    allowed keywords in kwargs are 'mode' and additional keywords depending
    on the type of subsystem you want to create
    """
    def __init__(self, label, **kwargs):

        if not isinstance(label, str):
            raise MyComponentError('Subsystem label must be of type str!')

        elif len([x for x in [';', ', ', '.'] if x in label]) > 0:
            raise MyComponentError('Can\'t use ' + str([';', ', ', '.']) + ' ',
                                   'in label.')
        else:
            self.label = label

        # set default values
        for key in self.attr():
            self.__dict__.update({key: np.nan})
            self.__dict__.update({key + '_set': False})

        self.subsys_init()
        self.set_attr(**kwargs)

    def set_attr(self, **kwargs):
        """
        set attributes of subsystem
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
                msg = ('Component ' + self.label + ' has no attribute ' +
                       str(key))
                raise ValueError(msg)

        self.set_comps()
        self.set_conns()

    def get_attr(self, key):
        """
        get attribute of subsystem
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            print(self, ' \"', self.label, '\" has no attribute \"',
                  key, '\"')
            return None

    def display_information(self):
        """
        print all attributes of subsystem
        """
        print(self.__dict__)

    def attr(self):
        return ['num_i', 'num_o']

    def from_network(self, nw):
        """
        create a subsystem from a network

        .. note::
            Does not work by now!
        """
        num_inter = 0
        for c in nw.conns.index:
            self.nw.add_conns(c)
            if isinstance(c, cmp.subsys_interface):
                num_inter += 1

        if num_inter != 2:
            raise MyNetworkError('A subsystem must exactly have one inlet '
                                 'and one outlet interface!')

    def subsys_init(self):
        """
        subsystem initialisation:

        - create components and connections
        - set components and connections parameters
        - create network
        """
        self.create_comps()
        self.create_conns()

    def create_comps(self):
        return

    def create_conns(self):
        return

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
    """
    **available parameters**

    - pr1_eva: pressure drop at hot side of evaporator
    - pr2_eva: pressure drop at cold side of evaporator
    - eta_s: isentropic efficiency of the pump
    - PP: pinch point
    - circ_num: circulation number, ratio of mass flow through cold side of
      evaporator to mass flow at the drum inlet

    **inlets and outlets**

    - in1, in2
    - out1, out2

    .. image:: _images/subsys_dr_eva_forced.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
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
    """
    **available parameters**

    - pr1_eva: pressure drop at hot side of evaporator
    - PP: pinch point
    - circ_num: circulation number, ratio of mass flow through cold side of
      evaporator to mass flow at the drum inlet

    **inlets and outlets**

    - in1, in2
    - out1, out2

    .. image:: _images/subsys_dr_eva_natural.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
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
    """
    **available parameters**

    - ttd: upper terminal temperature difference of condenser
    - pr1_desup: pressure drop at hot side of desuperheater
    - pr2_desup: pressure drop at cold side of desuperheater
    - pr1_cond: pressure drop at hot side of condenser
    - pr2_cond: pressure drop at cold side of condenser

    **inlets and outlets**

    - in1, in2
    - out1, out2

    .. image:: _images/subsys_ph_desup_cond.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
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
        self.vessel = cmp.vessel(label=self.label + '_vessel', mode='man')

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
        self.conns += [connection(self.condenser, 'out1', self.vessel, 'in1')]
        self.conns += [connection(self.vessel, 'out1', self.outlet, 'in1')]
        self.conns += [connection(self.inlet, 'out2', self.condenser, 'in2')]
        self.conns += [connection(self.condenser, 'out2', self.desup, 'in2')]
        self.conns += [connection(self.desup, 'out2', self.outlet, 'in2')]

# %%


class ph_desup_cond_subc(subsystem):
    """
    **available parameters**

    - ttd: upper terminal temperature difference of condenser
    - pr1_desup: pressure drop at hot side of desuperheater
    - pr2_desup: pressure drop at cold side of desuperheater
    - pr1_cond: pressure drop at hot side of condenser
    - pr2_cond: pressure drop at cold side of condenser
    - pr1_subc: pressure drop at hot side of subcooler
    - pr2_subc: pressure drop at cold side of subcooler

    **inlets and outlets**

    - in1, in2
    - out1, out2

    .. image:: _images/subsys_ph_desup_cond_subc.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
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
        self.vessel = cmp.vessel(label=self.label + '_vessel', mode='man')

    def set_comps(self):

        self.condenser.set_attr(ttd_u=self.ttd)
        self.desup.set_attr(pr1=self.pr1_desup)
        self.desup.set_attr(pr2=self.pr2_desup)
        self.condenser.set_attr(pr1=self.pr1_cond)
        self.condenser.set_attr(pr2=self.pr2_cond)
        self.subcooler.set_attr(pr1=self.pr1_subc)
        self.subcooler.set_attr(pr2=self.pr2_subc)
        self.vessel.set_attr(pr=self.pr_v)

    def create_conns(self):

        self.conns = []

        self.conns += [connection(self.inlet, 'out1', self.desup, 'in1')]
        self.conns += [connection(self.desup, 'out1', self.condenser, 'in1')]
        self.conns += [connection(self.condenser, 'out1',
                                  self.subcooler, 'in1')]
        self.conns += [connection(self.subcooler, 'out1', self.vessel, 'in1')]
        self.conns += [connection(self.vessel, 'out1', self.outlet, 'in1')]
        self.conns += [connection(self.inlet, 'out2', self.subcooler, 'in2')]
        self.conns += [connection(self.subcooler, 'out2',
                                  self.condenser, 'in2')]
        self.conns += [connection(self.condenser, 'out2', self.desup, 'in2')]
        self.conns += [connection(self.desup, 'out2', self.outlet, 'in2')]

# %%


class ph_desup_inl_cond_subc(subsystem):
    """
    **available parameters**

    - ttd: upper terminal temperature difference of condenser
    - pr1_desup: pressure drop at hot side of desuperheater
    - pr2_desup: pressure drop at cold side of desuperheater
    - pr1_cond: pressure drop at hot side of condenser
    - pr2_cond: pressure drop at cold side of condenser
    - pr1_subc: pressure drop at hot side of subcooler
    - pr2_subc: pressure drop at cold side of subcooler

    **inlets and outlets**

    - in1, in2, in3
    - out1, out2

    .. image:: _images/subsys_ph_desup_inl_cond_subc.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
    """

    def attr(self):
        return ([n for n in subsystem.attr(self) if
                 n != 'num_i' and n != 'num_o'] +
                ['ttd', 'pr1_desup', 'pr2_desup',
                 'pr1_cond', 'pr2_cond', 'pr1_subc', 'pr2_subc'])

    def create_comps(self):

        self.num_i = 3
        self.num_o = 2
        self.inlet = cmp.subsys_interface(label=self.label + '_inlet',
                                          num_inter=self.num_i)
        self.outlet = cmp.subsys_interface(label=self.label + '_outlet',
                                           num_inter=self.num_o)
        self.desup = cmp.desuperheater(label=self.label + '_desup')
        self.condenser = cmp.condenser(label=self.label + '_condenser')
        self.subcooler = cmp.heat_exchanger(label=self.label + '_subcooler')
        self.vessel = cmp.vessel(label=self.label + '_vessel', mode='man')
        self.merge = cmp.merge(label=self.label + '_merge')

    def set_comps(self):

        self.condenser.set_attr(ttd_u=self.ttd)
        self.desup.set_attr(pr1=self.pr1_desup)
        self.desup.set_attr(pr2=self.pr2_desup)
        self.condenser.set_attr(pr1=self.pr1_cond)
        self.condenser.set_attr(pr2=self.pr2_cond)
        self.subcooler.set_attr(pr1=self.pr1_cond)
        self.subcooler.set_attr(pr2=self.pr2_cond)

    def create_conns(self):

        self.conns = []

        self.conns += [connection(self.inlet, 'out1', self.desup, 'in1')]
        self.conns += [connection(self.desup, 'out1', self.merge, 'in1')]
        self.conns += [connection(self.inlet, 'out3', self.merge, 'in2')]
        self.conns += [connection(self.merge, 'out1', self.condenser, 'in1')]
        self.conns += [connection(self.condenser, 'out1',
                                  self.subcooler, 'in1')]
        self.conns += [connection(self.subcooler, 'out1', self.vessel, 'in1')]
        self.conns += [connection(self.vessel, 'out1', self.outlet, 'in1')]
        self.conns += [connection(self.inlet, 'out2', self.subcooler, 'in2')]
        self.conns += [connection(self.subcooler, 'out2',
                                  self.condenser, 'in2')]
        self.conns += [connection(self.condenser, 'out2', self.desup, 'in2')]
        self.conns += [connection(self.desup, 'out2', self.outlet, 'in2')]
