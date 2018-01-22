"""
.. module:: components
    :platforms: all
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import numpy as np

from tespy.helpers import (
    MyComponentError, MyNetworkError
)

from tespy import networks as nwk
from tespy.connections import connection, ref
from tespy.components import components as comp


class subsystem:
    r"""

    :param label: label for subsystem
    :type label: str
    :param **kwargs: for the keyword arguments see :code:`subsystem.attr()`
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

    **initialisation method is used for instances of class component and its
    children**

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

# set provided values,  check for invalid keys
        invalid_keys = np.array([])
        for key in kwargs:
            if key not in self.attr():
                invalid_keys = np.append(invalid_keys, key)
            else:
                if type(kwargs[key]) == float or type(kwargs[key]) == int:
                    self.__dict__.update({key: kwargs[key]})
                    self.__dict__.update({key + '_set': True})
                else:
                    raise MyComponentError('Specified value does not match '
                                           'requirements. Only numbers are '
                                           'allowed as parameters.')

# print invalid keywords
        if len(invalid_keys) > 0:
            print('\'', invalid_keys, '\' are invalid attributes.',
                  'Available attributes for object \'', self,
                  '\' are:',  self.attr())

        self.subsys_init()

    def set_attr(self, **kwargs):
        """
        set attribute of subsystem
        """

        invalid_keys = np.array([])
        for key in kwargs:
            if key not in self.attr():
                invalid_keys = np.append(invalid_keys, key)
            else:
                if np.isnan(kwargs[key]):
                    self.__dict__.update({key: np.nan})
                    self.__dict__.update({key + '_set': False})
                elif (type(kwargs[key]) == float or
                      type(kwargs[key]) == np.float64 or
                      type(kwargs[key]) == int or
                      key == 'fuel'):
                    self.__dict__.update({key: kwargs[key]})
                    self.__dict__.update({key + '_set': True})
                else:
                    msg = ('Specified value does not match requirements. '
                           'Only numeric parameters are allowed.')
                    raise TypeError(msg)

        if len(invalid_keys) > 0:
            print('\'', invalid_keys, '\' are invalid attributes. '
                  'Available attributes for object \'',
                  self.__class__.__name__, '\' are:', self.attr())

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
            if isinstance(c, comp.subsys_interface):
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
        self.set_comps()

        self.create_conns()
        self.set_conns()

        self.create_network()

    def create_comps(self):
        return

    def create_conns(self):
        return

    def set_comps(self):
        return

    def set_conns(self):
        return

    def create_network(self):
        self.nw = nwk.network()
        for c in self.conns:
            self.nw.add_conns(c)

    def get_network(self):
        return self.nw

# %%


class dr_eva_forced(subsystem):
    """
    **available parameters**

    - dp1_eva: pressure drop at hot side of evaporator
    - dp2_eva: pressure drop at cold side of evaporator
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
                ['dp1_eva', 'dp2_eva', 'eta_s', 'PP', 'circ_num'])

    def create_comps(self):

        self.num_i = 2
        self.num_o = 2
        self.inlet = comp.subsys_interface(label=self.label + '_inlet',
                                           num_inter=self.num_i)
        self.outlet = comp.subsys_interface(label=self.label + '_outlet',
                                            num_inter=self.num_o)
        self.pump = comp.pump(label=self.label + '_pump')
        self.drum = comp.drum(label=self.label + '_drum')
        self.evaporator = comp.heat_exchanger(label=self.label + '_evaporator')

    def set_comps(self):

        self.evaporator.set_attr(ttd_l=self.PP)
        self.evaporator.set_attr(dp1=self.dp1_eva)
        self.evaporator.set_attr(dp2=self.dp2_eva)
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

    - dp1_eva: pressure drop at hot side of evaporator
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
                ['dp1_eva', 'PP', 'circ_num'])

    def create_comps(self):

        self.num_i = 2
        self.num_o = 2
        self.inlet = comp.subsys_interface(label=self.label + '_inlet',
                                           num_inter=self.num_i)
        self.outlet = comp.subsys_interface(label=self.label + '_outlet',
                                            num_inter=self.num_o)
        self.drum = comp.drum(label=self.label + '_drum')
        self.evaporator = comp.heat_exchanger(label=self.label + '_evaporator',
                                              mode='man')

    def set_comps(self):

        self.evaporator.set_attr(ttd_l=self.PP)
        self.evaporator.set_attr(dp1=self.dp1_eva)

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
    - dp1_desup: pressure drop at hot side of desuperheater
    - dp2_desup: pressure drop at cold side of desuperheater
    - dp1_cond: pressure drop at hot side of condenser
    - dp2_cond: pressure drop at cold side of condenser

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
                ['ttd', 'dp1_desup', 'dp2_desup', 'dp1_cond', 'dp2_cond'])

    def create_comps(self):

        self.num_i = 2
        self.num_o = 2
        self.inlet = comp.subsys_interface(label=self.label + '_inlet',
                                           num_inter=self.num_i)
        self.outlet = comp.subsys_interface(label=self.label + '_outlet',
                                            num_inter=self.num_o)
        self.desup = comp.desuperheater(label=self.label + '_desup')
        self.condenser = comp.condenser(label=self.label + '_condenser')
        self.vessel = comp.vessel(label=self.label + '_vessel', mode='man')

    def set_comps(self):

        self.condenser.set_attr(ttd_u=self.ttd)
        self.desup.set_attr(dp1=self.dp1_desup)
        self.desup.set_attr(dp2=self.dp2_desup)
        self.condenser.set_attr(dp1=self.dp1_cond)
        self.condenser.set_attr(dp2=self.dp2_cond)

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

    def attr(self):
        return ([n for n in subsystem.attr(self) if
                 n != 'num_i' and n != 'num_o'] +
                ['ttd', 'dp1_desup', 'dp2_desup',
                 'dp1_cond', 'dp2_cond', 'dp1_subc', 'dp2_subc'])

    def create_comps(self):

        self.num_i = 2
        self.num_o = 2
        self.inlet = comp.subsys_interface(label=self.label + '_inlet',
                                           num_inter=self.num_i)
        self.outlet = comp.subsys_interface(label=self.label + '_outlet',
                                            num_inter=self.num_o)
        self.desup = comp.desuperheater(label=self.label + '_desup')
        self.condenser = comp.condenser(label=self.label + '_condenser')
        self.subcooler = comp.heat_exchanger(label=self.label + '_subcooler')
        self.vessel = comp.vessel(label=self.label + '_vessel', mode='man')

    def set_comps(self):

        self.condenser.set_attr(ttd_u=self.ttd)
        self.desup.set_attr(dp1=self.dp1_desup)
        self.desup.set_attr(dp2=self.dp2_desup)
        self.condenser.set_attr(dp1=self.dp1_cond)
        self.condenser.set_attr(dp2=self.dp2_cond)
        self.subcooler.set_attr(dp1=self.dp1_cond)
        self.subcooler.set_attr(dp2=self.dp2_cond)

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

    def attr(self):
        return ([n for n in subsystem.attr(self) if
                 n != 'num_i' and n != 'num_o'] +
                ['ttd', 'dp1_desup', 'dp2_desup',
                 'dp1_cond', 'dp2_cond', 'dp1_subc', 'dp2_subc'])

    def create_comps(self):

        self.num_i = 3
        self.num_o = 2
        self.inlet = comp.subsys_interface(label=self.label + '_inlet',
                                           num_inter=self.num_i)
        self.outlet = comp.subsys_interface(label=self.label + '_outlet',
                                            num_inter=self.num_o)
        self.desup = comp.desuperheater(label=self.label + '_desup')
        self.condenser = comp.condenser(label=self.label + '_condenser')
        self.subcooler = comp.heat_exchanger(label=self.label + '_subcooler')
        self.vessel = comp.vessel(label=self.label + '_vessel', mode='man')
        self.merge = comp.merge(label=self.label + '_merge')

    def set_comps(self):

        self.condenser.set_attr(ttd_u=self.ttd)
        self.desup.set_attr(dp1=self.dp1_desup)
        self.desup.set_attr(dp2=self.dp2_desup)
        self.condenser.set_attr(dp1=self.dp1_cond)
        self.condenser.set_attr(dp2=self.dp2_cond)
        self.subcooler.set_attr(dp1=self.dp1_cond)
        self.subcooler.set_attr(dp2=self.dp2_cond)

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
