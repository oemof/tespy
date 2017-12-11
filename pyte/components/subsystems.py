"""
.. module:: components
    :platforms: all
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import numpy as np
import math

from pyte.helpers import (
    MyComponentError, MyNetworkError
)

from pyte import networks as nwk
from pyte.connections import connection, ref
from pyte.components import components as comp


class subsystem:
    """

    """

    def __init__(self, **kwargs):

        if 'label' not in kwargs:
            raise MyComponentError('Subsystems must have a label!')

        if not isinstance(kwargs['label'], str):
            raise MyComponentError('Subsystem label must be of type str!')

        if len([x for x in [';', ', ', '.'] if x in kwargs['label']]) > 0:
            raise MyComponentError('Can\'t use ' + str([';', ', ', '.']) + ' ',
                                   'in label.')

# set default values
        for key in self.attr():
            if key not in subsystem.attr(self):
                self.__dict__.update({key: 0})
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
                elif key == 'label':
                    self.label = kwargs['label']
                else:
                    raise MyComponentError('Specified value does not match '
                                           'requirements. Only numbers are '
                                           'allowed as parameters.')

# print invalid keywords
        if len(invalid_keys) > 0:
            print('\'', invalid_keys, '\' are invalid attributes.',
                  'Available attributes for object \'', self,
                  '\' are:',  self.attr())

        self.nw = nwk.network()
        self.comp_init()

    def get_attr(self, key):
        """
        Print specific attribute of component
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            print(self, ' \"', self.label, '\" has no attribute \"',
                  key, '\"')
            return None

    def display_information(self):
        """
        Print all attributes of component
        """
        print(self.__dict__)

    def attr(self):
        return ['label', 'num_i', 'num_o']

    def from_network(self, nw):
        num_inter = 0
        for c in nw.conns.index:
            self.nw.add_conns(c)
            if isinstance(c, comp.subsys_interface):
                num_inter += 1

        if num_inter != 2:
            raise MyNetworkError('A subsystem must exactly have one inlet '
                                 'and one outlet interface!')

    def comp_init(self):
        return

    def comps(self):
        return

    def conns(self):
        return

    def create_network(self):
        for c in self.conns():
            self.nw.add_conns(c)

    def get_network(self):
        return self.nw

# %%


class dr_eva_forced(subsystem):

    def attr(self):
        return ([n for n in subsystem.attr(self) if
                 n != 'num_i' and n != 'num_o'] +
                ['dp1_eva', 'dp2_eva', 'eta_s', 'PP', 'circ_num'])

    def comp_init(self):
        self.num_i = 2
        self.num_o = 2
        self.comps()
        self.conns()
        self.create_network()

    def comps(self):
        self.inlet = comp.subsys_interface(label=self.label + 'inlet',
                                           num_inter=self.num_i)
        self.outlet = comp.subsys_interface(label=self.label + 'outlet',
                                            num_inter=self.num_o)
        self.pump = comp.pump(label=self.label + 'pump')
        self.drum = comp.drum(label=self.label + 'drum')
        self.evaporator = comp.heat_exchanger(label=self.label + 'evaporator')
        if self.PP_set:
            self.evaporator.set_attr(ttd_l=self.PP)
        if self.dp1_eva_set:
            self.evaporator.set_attr(dp1=self.dp1_eva)
        if self.dp2_eva_set:
            self.evaporator.set_attr(dp2=self.dp2_eva)
        if self.eta_s_set:
            self.pump.set_attr(eta_s=self.eta_s)

    def conns(self):
        conns = [connection(self.inlet, 'out1', self.evaporator, 'in1'),
                 connection(self.evaporator, 'out1', self.outlet, 'in1'),
                 connection(self.inlet, 'out2', self.drum, 'in1'),
                 connection(self.drum, 'out1', self.pump, 'in1'),
                 connection(self.pump, 'out1', self.evaporator, 'in2'),
                 connection(self.evaporator, 'out2', self.drum, 'in2'),
                 connection(self.drum, 'out2', self.outlet, 'in2')]
        if self.circ_num_set:
            conns[3].set_attr(m=ref(conns[-1], self.circ_num, 0))
        return conns

# %%


class dr_eva_natural(subsystem):

    def attr(self):
        return ([n for n in subsystem.attr(self) if
                 n != 'num_i' and n != 'num_o'] +
                ['dp1_eva', 'PP', 'circ_num'])

    def comp_init(self):
        self.num_i = 2
        self.num_o = 2
        self.comps()
        self.conns()
        self.create_network()

    def comps(self):
        self.inlet = comp.subsys_interface(label=self.label + 'inlet',
                                           num_inter=self.num_i)
        self.outlet = comp.subsys_interface(label=self.label + 'outlet',
                                            num_inter=self.num_o)
        self.drum = comp.drum(label=self.label + 'drum')
        self.evaporator = comp.heat_exchanger(label=self.label + 'evaporator', mode='man')
        if self.PP_set:
            self.evaporator.set_attr(ttd_l=self.PP)
        if self.dp1_eva_set:
            self.evaporator.set_attr(dp1=self.dp1_eva)

    def conns(self):
        conns = [connection(self.inlet, 'out1', self.evaporator, 'in1'),
                 connection(self.evaporator, 'out1', self.outlet, 'in1'),
                 connection(self.inlet, 'out2', self.drum, 'in1'),
                 connection(self.drum, 'out1', self.evaporator, 'in2'),
                 connection(self.evaporator, 'out2', self.drum, 'in2'),
                 connection(self.drum, 'out2', self.outlet, 'in2')]
        if self.circ_num_set:
            conns[3].set_attr(m=ref(conns[-1], self.circ_num, 0))
        return conns

# %%


class ph_desup_cond(subsystem):

    def attr(self):
        return ([n for n in subsystem.attr(self) if
                 n != 'num_i' and n != 'num_o'] +
                ['ttd', 'dp1_desup', 'dp2_desup', 'dp1_cond', 'dp2_cond'])

    def comp_init(self):
        self.num_i = 2
        self.num_o = 2
        self.comps()
        self.conns()
        self.create_network()

    def comps(self):
        self.inlet = comp.subsys_interface(label=self.label + 'inlet',
                                           num_inter=self.num_i)
        self.outlet = comp.subsys_interface(label=self.label + 'outlet',
                                            num_inter=self.num_o)
        self.desup = comp.desuperheater(label=self.label + 'desup')
        self.condenser = comp.condenser(label=self.label + 'condenser')
        self.vessel = comp.vessel(label=self.label + 'vessel', mode='man')
        if self.ttd_set:
            self.condenser.set_attr(ttd_u=self.ttd)
        if self.dp1_desup_set:
            self.desup.set_attr(dp1=self.dp1_desup)
        if self.dp2_desup_set:
            self.desup.set_attr(dp2=self.dp2_desup)
        if self.dp1_cond_set:
            self.condenser.set_attr(dp1=self.dp1_cond)
        if self.dp2_cond_set:
            self.condenser.set_attr(dp2=self.dp2_cond)

    def conns(self):
        conns = [connection(self.inlet, 'out1', self.desup, 'in1'),
                 connection(self.desup, 'out1', self.condenser, 'in1'),
                 connection(self.condenser, 'out1', self.vessel, 'in1'),
                 connection(self.vessel, 'out1', self.outlet, 'in1'),
                 connection(self.inlet, 'out2', self.condenser, 'in2'),
                 connection(self.condenser, 'out2', self.desup, 'in2'),
                 connection(self.desup, 'out2', self.outlet, 'in2')]
        return conns

# %%


class ph_desup_cond_subc(subsystem):

    def attr(self):
        return ([n for n in subsystem.attr(self) if
                 n != 'num_i' and n != 'num_o'] +
                ['ttd', 'dp1_desup', 'dp2_desup',
                 'dp1_cond', 'dp2_cond', 'dp1_subc', 'dp2_subc'])

    def comp_init(self):
        self.num_i = 2
        self.num_o = 2
        self.comps()
        self.conns()
        self.create_network()

    def comps(self):
        self.inlet = comp.subsys_interface(label=self.label + 'inlet',
                                           num_inter=self.num_i)
        self.outlet = comp.subsys_interface(label=self.label + 'outlet',
                                            num_inter=self.num_o)
        self.desup = comp.desuperheater(label=self.label + 'desup')
        self.condenser = comp.condenser(label=self.label + 'condenser')
        self.subcooler = comp.heat_exchanger(label=self.label + 'subcooler')
        self.vessel = comp.vessel(label=self.label + 'vessel', mode='man')
        if self.ttd_set:
            self.condenser.set_attr(ttd_u=self.ttd)
        if self.dp1_desup_set:
            self.desup.set_attr(dp1=self.dp1_desup)
        if self.dp2_desup_set:
            self.desup.set_attr(dp2=self.dp2_desup)
        if self.dp1_cond_set:
            self.condenser.set_attr(dp1=self.dp1_cond)
        if self.dp2_cond_set:
            self.condenser.set_attr(dp2=self.dp2_cond)
        if self.dp1_subc_set:
            self.subcooler.set_attr(dp1=self.dp1_cond)
        if self.dp2_subc_set:
            self.subcooler.set_attr(dp2=self.dp2_cond)

    def conns(self):
        conns = [connection(self.inlet, 'out1', self.desup, 'in1'),
                 connection(self.desup, 'out1', self.condenser, 'in1'),
                 connection(self.condenser, 'out1', self.subcooler, 'in1'),
                 connection(self.subcooler, 'out1', self.vessel, 'in1'),
                 connection(self.vessel, 'out1', self.outlet, 'in1'),
                 connection(self.inlet, 'out2', self.subcooler, 'in2'),
                 connection(self.subcooler, 'out2', self.condenser, 'in2'),
                 connection(self.condenser, 'out2', self.desup, 'in2'),
                 connection(self.desup, 'out2', self.outlet, 'in2')]
        return conns

# %%


class ph_desup_inl_cond_subc(subsystem):

    def attr(self):
        return ([n for n in subsystem.attr(self) if
                 n != 'num_i' and n != 'num_o'] +
                ['ttd', 'dp1_desup', 'dp2_desup',
                 'dp1_cond', 'dp2_cond', 'dp1_subc', 'dp2_subc'])

    def comp_init(self):
        self.num_i = 3
        self.num_o = 2
        self.comps()
        self.conns()
        self.create_network()

    def comps(self):
        self.inlet = comp.subsys_interface(label=self.label + 'inlet',
                                           num_inter=self.num_i)
        self.outlet = comp.subsys_interface(label=self.label + 'outlet',
                                            num_inter=self.num_o)
        self.desup = comp.desuperheater(label=self.label + 'desup')
        self.condenser = comp.condenser(label=self.label + 'condenser')
        self.subcooler = comp.heat_exchanger(label=self.label + 'subcooler')
        self.vessel = comp.vessel(label=self.label + 'vessel', mode='man')
        self.merge = comp.merge(label=self.label + 'merge')
        if self.ttd_set:
            self.condenser.set_attr(ttd_u=self.ttd)
        if self.dp1_desup_set:
            self.desup.set_attr(dp1=self.dp1_desup)
        if self.dp2_desup_set:
            self.desup.set_attr(dp2=self.dp2_desup)
        if self.dp1_cond_set:
            self.condenser.set_attr(dp1=self.dp1_cond)
        if self.dp2_cond_set:
            self.condenser.set_attr(dp2=self.dp2_cond)
        if self.dp1_subc_set:
            self.subcooler.set_attr(dp1=self.dp1_cond)
        if self.dp2_subc_set:
            self.subcooler.set_attr(dp2=self.dp2_cond)

    def conns(self):
        conns = [connection(self.inlet, 'out1', self.desup, 'in1'),
                 connection(self.desup, 'out1', self.merge, 'in1'),
                 connection(self.inlet, 'out2', self.merge, 'in2'),
                 connection(self.merge, 'out1', self.condenser, 'in1'),
                 connection(self.condenser, 'out1', self.subcooler, 'in1'),
                 connection(self.subcooler, 'out1', self.vessel, 'in1'),
                 connection(self.vessel, 'out1', self.outlet, 'in1'),
                 connection(self.inlet, 'out2', self.subcooler, 'in2'),
                 connection(self.subcooler, 'out2', self.condenser, 'in2'),
                 connection(self.condenser, 'out2', self.desup, 'in2'),
                 connection(self.desup, 'out2', self.outlet, 'in2')]
        return conns
