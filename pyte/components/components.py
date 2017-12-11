"""
.. module:: components
    :platforms: all
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import numpy as np
import math

import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI as CPPSI

from pyte.helpers import (
    num_fluids, fluid_structure, MyComponentError,
    v_mix_ph, h_mix_pT, s_mix_pT, s_mix_ph, T_mix_ph, T_mix_ps,
    dT_mix_dph, dT_mix_pdh, dT_mix_ph_dfluid, h_mix_pQ, dh_mix_dpQ,
    molar_masses, err
)

from pyte.components import characteristics as cmp_char


def init_target(nw, c, start):
    """
    propagates the fluids towards connections target,
    ends when reaching sink, merge or combustion chamber

    :param nw: network to operate on
    :type nw: pyte.networks.network
    :param c: connection to initialise
    :type c: pyte.connections.connection
    :param start: fluid propagation startingpoint, in some cases needed
        to exit the recursion
    :type start: pyte.connections.connection
    :returns: no return value

    .. note::
        This method is the same as the method in the network class of the
        networks module. This is necessary as the combustion chambers
        convergence check requires the method while the networks module
        requires the components module. Check, if the cicular imports can be
        avoided in a more elegant way.
    """
    if (len(c.t.inlets()) == 1 and len(c.t.outlets()) == 1 or
        isinstance(c.t, heat_exchanger) or
        isinstance(c.t, subsys_interface)):

        inconn = [x for x in nw.comps.loc[c.s].o if
                  x in nw.comps.loc[c.t].i]
        inconn_id = nw.comps.loc[c.t].i.tolist().index(inconn[0])
        outconn = nw.comps.loc[c.t].o.tolist()[inconn_id]
        for fluid, x in c.fluid.items():
            if not outconn.fluid_set[fluid]:
                outconn.fluid[fluid] = x

        init_target(nw, outconn, start)

    if isinstance(c.t, splitter):
        for outconn in nw.comps.loc[c.t].o:
            for fluid, x in c.fluid.items():
                if not outconn.fluid_set[fluid]:
                    outconn.fluid[fluid] = x

            init_target(nw, outconn, start)

    if isinstance(c.t, drum) and c.t != start:
        start = c.t
        for outconn in nw.comps.loc[c.t].o:
            for fluid, x in c.fluid.items():
                if not outconn.fluid_set[fluid]:
                    outconn.fluid[fluid] = x

            init_target(nw, outconn, start)


class component:
    """
    class component
    """
    printOnCreate = False
    printOnChange = False

    def __init__(self, label, **kwargs):
        r"""
        object initialisation

        allowed keywords in kwargs are 'mode' and additional keywords depending
        on the type of component you want to specify

        :param label: label for component
        :type label: str
        :param **kwargs: for the keyword arguments see :code:`component.attr()`
        :returns: no return value
        :raises: - :code:`TypeError`, if label is not of type str
                   components
                 - :code:`ValueError`, if label contains forbidden characters
                   (';', ',', '.')
                 - :code:`ValueError`, if outlet_id or inlet_id are not allowed
                   for ids for comp1 or comp2

        **example**

        .. code-block:: python

            cond = condenser('main condenser', ttd_u=5)

        creates component condenser labeled 'main condenser' and sets the
        terminal temperature difference at the upper side (hot side inlet to
        cold side outlet) to 5 K
        """

        if not isinstance(label, str):
            msg = 'Component label must be of type str!'
            raise TypeError(msg)
        elif len([x for x in [';', ',', '.'] if x in label]) > 0:
            msg = ('Can\'t use ' + str([';', ',', '.']) + ' ',
                   'in label (' + str(self.component()) + ').')
            raise ValueError(msg)
        else:
            self.label = label

        self.mode = kwargs.get('mode', 'auto')

        # set default values
        for key in self.attr():
            if key != 'mode':
                self.__dict__.update({key: 0})
                self.__dict__.update({key + '_set': False})

        # set provided values, check for invalid keys
        invalid_keys = np.array([])
        for key in kwargs:
            if key != 'mode':
                if key not in self.attr():
                    invalid_keys = np.append(invalid_keys, key)
                else:
                    if (type(kwargs[key]) == float or
                        type(kwargs[key]) == np.float64 or
                        type(kwargs[key]) == int or
                        key == 'fuel'):
                        self.__dict__.update({key: kwargs[key]})
                        self.__dict__.update({key + '_set': True})
                    else:
                        msg = ('Specified value does not match requirements. '
                               'Only numeric parameters are allowed.')
                        raise TypeError(msg)

        # print invalid keywords
        if len(invalid_keys) > 0:
            print('\'', invalid_keys, '\' are invalid attributes.',
                  'Available attributes for object \'', self,
                  '\' are:', self.attr())

        if component.printOnCreate:
            print('Created ', self, '.')
            print(self.__dict__)

    def set_attr(self, **kwargs):
        """
        sets, resets or unsets attributes of a connection, for the keyword
        arguments, return values and errors see object initialisation
        """
        for key in kwargs:
            print(type(kwargs[key]))
        invalid_keys = np.array([])
        for key in kwargs:
            if key not in self.attr():
                invalid_keys = np.append(invalid_keys, key)
            else:
                if np.isnan(kwargs[key]):
                    self.__dict__.update({key + '_set': False})
                elif (type(kwargs[key]) == float or
                      type(kwargs[key]) == np.float64 or
                      type(kwargs[key]) == int or
                      key == 'fuel'):
                    self.__dict__.update({key: kwargs[key]})
                    self.__dict__.update({key + '_set': True})
                elif type(kwargs[key]) == str and key == 'mode':
                    if kwargs[key] != '':
                        self.__dict__.update({key: kwargs[key]})
                    else:
                        msg = 'mode must be string and cannot be empty.'
                        raise TypeError(msg)

                else:
                    msg = ('Specified value does not match requirements. '
                           'Only numeric parameters are allowed.')
                    raise TypeError(msg)

        if len(invalid_keys) > 0:
            print('\'', invalid_keys, '\' are invalid attributes. '
                  'Available attributes for object \'', self.component(),
                  '\' are:', self.attr())

        if component.printOnChange:
            print('Updated ', self, '.')
            print(self.__dict__)

    def get_attr(self, key):
        """
        get the value of a components attribute

        :param key: attribute to return its value
        :type key: str
        :returns:
            - :code:`self.__dict__[key]` if object has attribute key
            - :code:`None` if object has no attribute key
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            print(self.component(), ' \"', self.label, '\" '
                  'has no attribute \"', key, '\"')
            return None

    def comp_init(self, nw):
        return

    def attr(self):
        return ['mode']

    def inlets(self):
        return []

    def outlets(self):
        return []

    def design(self):
        return []

    def offdesign(self):
        return []

    def equations(self, nw):
        return []

    def derivatives(self, nw):
        return []

    def initialise_source_p(self, c):
        return 0

    def initialise_target_p(self, c):
        return 0

    def initialise_source_h(self, c):
        return 0

    def initialise_target_h(self, c):
        return 0

    def calc_parameters(self, inlets, outlets, mode):
        return

    def print_parameters(self, inlets, outlets):
        return

    def initialise_fluids(self, nw):
        return

    def convergence_check(self, nw):
        return

# %%

    def fluid_res(self, inlets, outlets):
        r"""
        returns residual values for fluid equations

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: vec_res (*list*) - a list containing the residual values

        **components with one inlet and one outlet**

        .. math:: 0 = fluid_{i,in} - fluid_{i,out} \;
            \forall i \in \mathrm{fluid}

        **component heat exchanger or subsystem interface**

        .. math:: 0 = fluid_{i,in_{j}} - fluid_{i,out_{j}} \;
            \forall i \in \mathrm{fluid}, \; \forall j \in inlets/outlets

        **component splitter**

        .. math:: 0 = fluid_{i,in} - fluid_{i,out_{j}} \;
            \forall i \in \mathrm{fluid}, \; \forall j \in outlets

        **component merge**

        .. math::
            0 = \dot{m}_{in_{j}} \cdot fluid_{i,in_{j}} -
                \dot {m}_{out} \cdot fluid_{i,out} \\
            \forall i \in \mathrm{fluid}, \; \forall j \in inlets

        **component drum**

        .. math::
            0 = fluid_{i,in_1} - fluid_{i,out_{j}} \;
            \forall i \in \mathrm{fluid}, \; \forall j \in inlets

        **component separator**

        .. math::
            0 = \dot{m}_{in} \cdot fluid_{i,in} -
                \dot {m}_{out_{j}} \cdot fluid_{i,out_{j}} \;
            \forall i \in \mathrm{fluid}, \; \forall j \in outlets

        """
        vec_res = []

        if len(self.inlets()) == 1 and len(self.outlets()) == 1:
            for fluid, x in inlets[0].fluid.items():
                vec_res += [x - outlets[0].fluid[fluid]]
            return vec_res

        if (isinstance(self, subsys_interface) or
            isinstance(self, heat_exchanger)):
            for i in range(len(inlets)):
                for fluid, x in inlets[i].fluid.items():
                    vec_res += [x - outlets[i].fluid[fluid]]
            return vec_res

        if isinstance(self, splitter):
            for o in outlets:
                for fluid, x in inlets[0].fluid.items():
                    vec_res += [x - o.fluid[fluid]]
            return vec_res

        if isinstance(self, merge):
            res = 0
            for fluid, x in outlets[0].fluid.items():
                res = -x * outlets[0].m
                for i in inlets:
                    res += i.fluid[fluid] * i.m
                vec_res += [res]
            return vec_res

        if isinstance(self, drum):
            for o in outlets:
                for fluid, x in inlets[0].fluid.items():
                    vec_res += [x - o.fluid[fluid]]
            return vec_res

        if isinstance(self, separator):

            for fluid, x in inlets[0].fluid.items():
                res = x * inlets[0].m
                for o in outlets:
                    res -= o.fluid[fluid] * o.m
                vec_res += [res]
            return vec_res

        if isinstance(self, source) or isinstance(self, sink):
            return None

    def fluid_deriv(self, inlets, outlets):
        r"""
        returns derivatives for fluid equations

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: mat_deriv (*list*) - a list containing the derivatives

        **example:**

        component with one inlet and one outlet and 3 fluids in fluid vector

        .. math::
            \left(
            \begin{array}{cccccc}
                0 & 0 & 0 & 1 & 0 & 0\\
                0 & 0 & 0 & -1 & 0 & 0\\
            \end{array}
            \right)

            \left(
            \begin{array}{cccccc}
                0 & 0 & 0 & 0 & 1 & 0\\
                0 & 0 & 0 & 0 & -1 & 0\\
            \end{array}
            \right)

            \left(
            \begin{array}{cccccc}
                0 & 0 & 0 & 0 & 0 & 1\\
                0 & 0 & 0 & 0 & 0 & -1\\
            \end{array}
            \right)
        """
        num_i = len(inlets)
        num_o = len(outlets)
        num_fl = len(inlets[0].fluid)

        if len(self.inlets()) == 1 and len(self.outlets()) == 1:
            mat_deriv = np.zeros((num_fl, num_i + num_o, 3 + num_fl))
            i = 0
            for fluid, x in inlets[0].fluid.items():
                mat_deriv[i, 0, i + 3] = 1
                mat_deriv[i, 1, i + 3] = -1
                i += 1
            return mat_deriv.tolist()

        if isinstance(self, heat_exchanger):
            mat_deriv = np.zeros((num_fl * 2, num_i + num_o, 3 + num_fl))
            i = 0
            for fluid in inlets[0].fluid.keys():
                mat_deriv[i, 0, i + 3] = 1
                mat_deriv[i, 2, i + 3] = -1
                i += 1
            j = 0
            for fluid in inlets[1].fluid.keys():
                mat_deriv[i + j, 1, j + 3] = 1
                mat_deriv[i + j, 3, j + 3] = -1
                j += 1
            return mat_deriv.tolist()

        if isinstance(self, splitter):
            mat_deriv = np.zeros((num_fl * num_i * num_o,
                                  num_i + num_o, 3 + num_fl))
            k = 0
            for o in outlets:
                i = 0
                for fluid, x in inlets[0].fluid.items():
                    mat_deriv[i + k * num_fl, 0, i + 3] = 1
                    mat_deriv[i + k * num_fl, k + 1, i + 3] = -1
                    i += 1
                k += 1
            return mat_deriv.tolist()

        if isinstance(self, merge):
            mat_deriv = np.zeros((num_fl, num_i + num_o, 3 + num_fl))
            j = 0
            for fluid, x in outlets[0].fluid.items():
                k = 0
                for i in inlets:
                    mat_deriv[j, k, 0] = i.fluid[fluid]
                    mat_deriv[j, k, j + 3] = i.m
                    k += 1
                mat_deriv[j, k, 0] = -x
                mat_deriv[j, k, j + 3] = -outlets[0].m
                j += 1
            return mat_deriv.tolist()

        if isinstance(self, drum):
            mat_deriv = np.zeros((num_o * num_fl, num_i + num_o, 3 + num_fl))
            k = 0
            for o in outlets:
                i = 0
                for fluid, x in inlets[0].fluid.items():
                    mat_deriv[i + k * num_fl, 0, i + 3] = 1
                    mat_deriv[i + k * num_fl, k + 2, i + 3] = -1
                    i += 1
                k += 1
            return mat_deriv.tolist()

        if isinstance(self, separator):
            mat_deriv = np.zeros((num_fl, num_i + num_o, 3 + num_fl))
            j = 0
            for fluid, x in inlets[0].fluid.items():
                k = 0
                for o in outlets:
                    mat_deriv[j, k, 0] = -o.fluid[fluid]
                    mat_deriv[j, k, j + 3] = -o.m
                    k += 1
                mat_deriv[j, 0, 0] = x
                mat_deriv[j, 0, j + 3] = inlets[0].m
                j += 1
            return mat_deriv.tolist()

        if isinstance(self, source) or isinstance(self, sink):
            return None

        if isinstance(self, subsys_interface):
            mat_deriv = np.zeros((num_fl * num_i, num_i + num_o, 3 + num_fl))
            for i in range(num_i):
                j = 0
                for fluid in inlets[i].fluid.keys():
                    mat_deriv[i * num_fl + j, i, j + 3] = 1
                    mat_deriv[i * num_fl + j, num_i + i, j + 3] = -1
                    j += 1
            return mat_deriv.tolist()

# %%

    def mass_flow_res(self, inlets, outlets):
        r"""
        returns residual values for mass flow equations

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: vec_res (*list*) - a list containing the residual values

        **all components but heat exchanger and subsystem interface**

        .. math:: 0 = \sum \dot{m}_{in,i} - \sum \dot{m}_{out,j} \;
            \forall i \in inlets, \forall j \in outlets

        **heat exchanger and subsystem interface (same number of inlets and
        outlets)**

        .. math:: 0 = \dot{m}_{in,i} - \dot{m}_{out,i} \;
            \forall i \in inlets/outlets
        """

        if (isinstance(self, split) or
            isinstance(self, merge) or
            isinstance(self, combustion_chamber) or
            isinstance(self, drum) or
            (len(self.inlets()) == 1 and len(self.outlets()) == 1)):
            res = 0
            for i in inlets:
                res += i.m
            for o in outlets:
                res -= o.m
            return [res]

        if (isinstance(self, subsys_interface) or
            isinstance(self, heat_exchanger)):
            vec_res = []
            for i in range(len(inlets)):
                vec_res += [inlets[i].m - outlets[i].m]
            return vec_res

        if isinstance(self, source) or isinstance(self, sink):
            return None



    def mass_flow_deriv(self, inlets, outlets):
        r"""
        returns derivatives for mass flow equations

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: mat_deriv (*list*) - a list containing the derivatives

        **example**

        merge with three inlets and one outlet (one fluid in fluid vector)

        .. math::
            \left(
            \begin{array}{cccc}
                1 & 0 & 0 & 0\\
                1 & 0 & 0 & 0\\
                1 & 0 & 0 & 0\\
                -1 & 0 & 0 & 0\\
            \end{array}
            \right)
        """
        num_i = len(inlets)
        num_o = len(outlets)
        num_fl = len(inlets[0].fluid)

        if (isinstance(self, split) or
            isinstance(self, merge) or
            isinstance(self, combustion_chamber) or
            isinstance(self, drum) or
            (len(self.inlets()) == 1 and len(self.outlets()) == 1)):
            mat_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            j = 0
            for i in inlets:
                mat_deriv[0, j, 0] = 1
                j += 1
            k = 0
            for o in outlets:
                mat_deriv[0, k + j, 0] = -1
                k += 1
            return mat_deriv.tolist()

        if (isinstance(self, subsys_interface) or
            isinstance(self, heat_exchanger)):
            mat_deriv = np.zeros((num_i, num_i + num_o, num_fl + 3))
            for i in range(num_i):
                mat_deriv[i, i, 0] = 1
            for j in range(num_o):
                mat_deriv[j, j + i + 1, 0] = -1
            return mat_deriv.tolist()

        if isinstance(self, source) or isinstance(self, sink):
            return None
# %%

    def ddx_func(self, inlets, outlets, func, dx, pos):
        r"""
        calculates derivative of the function func to dx at components inlet or
        outlet in position pos

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :param func: function to calculate derivative
        :type func: function
        :param dx: dx
        :type dx: str
        :param pos: position of inlet or outlet, logic: ['in1', 'in2', ...,
                    'out1', ...] -> 0, 1, ..., n, n + 1, ...
        :type pos: int
        :returns: deriv (list or float) - partial derivative of the function
                  func to dx

        .. math::

            \frac{\partial f}{\partial x} = \frac{f(x + d) + f(x - d)}
            {2 \cdot d}
        """

        dm, dp, dh, df = 0, 0, 0, 0
        if dx == 'm':
            dm = 1e-3
        elif dx == 'p':
            dp = 1
        elif dx == 'h':
            dh = 1
        else:
            df = 1e-4

        if dx == 'fluid':
            deriv = []
            for f in inlets[0].fluid.keys():
                val = (inlets + outlets)[pos].fluid[f]
                exp = 0
                if (inlets + outlets)[pos].fluid[f] + df <= 1:
                    (inlets + outlets)[pos].fluid[f] += df
                else:
                    (inlets + outlets)[pos].fluid[f] = 1
                exp += func(inlets, outlets)
                if (inlets + outlets)[pos].fluid[f] - 2 * df >= 0:
                    (inlets + outlets)[pos].fluid[f] -= 2 * df
                else:
                    (inlets + outlets)[pos].fluid[f] = 0
                exp -= func(inlets, outlets)
                (inlets + outlets)[pos].fluid[f] = val

                deriv += [exp / (2 * (dm + dp + dh + df))]

        else:
            exp = 0
            (inlets + outlets)[pos].m += dm
            (inlets + outlets)[pos].p += dp
            (inlets + outlets)[pos].h += dh
            exp += func(inlets, outlets)

            (inlets + outlets)[pos].m -= 2 * dm
            (inlets + outlets)[pos].p -= 2 * dp
            (inlets + outlets)[pos].h -= 2 * dh
            exp -= func(inlets, outlets)
            deriv = exp / (2 * (dm + dp + dh + df))

            (inlets + outlets)[pos].m += dm
            (inlets + outlets)[pos].p += dp
            (inlets + outlets)[pos].h += dh

        return deriv

# %%

    def zeta_func(self, inlets, outlets):
        r"""
        calculates pressure drop from zeta (zeta1 for heat exchangers)

        :param inlets: components connections at inlets
        :type inlets: list
        :param outlets: components connections at outlets
        :type outlets: list
        :returns: residual value for the pressure drop

        .. math::

            \zeta = \frac{\Delta p \cdot v \cdot 2}{c^2}\\
            c = \frac{\dot{m} \cdot v}{A}

        As the surface area A will not change from design to offdesign
        calculation, it is possible to handle this the following way:

        .. math::
            0 = \zeta - \frac{(p_{in} - p_{out}) \cdot \pi^2}{8 \cdot
            \dot{m}_{in}^2 \cdot \frac{v_{in} + v_{out}}{2}}
        """
        i = inlets[0].as_list()
        o = outlets[0].as_list()
        if hasattr(self, 'zeta'):
            val = self.zeta
        else:
            val = self.zeta1
        return (val - (i[1] - o[1]) * math.pi ** 2 /
                (8 * i[0] ** 2 * (v_mix_ph(i) + v_mix_ph(o)) / 2))

# %%


class source(component):
    """
    component source
        - a flow originates from this component
    """

    def outlets(self):
        return ['out1']

    def component(self):
        return 'source'

# %%


class sink(component):
    """
    component sink
        - a flow drains in this component
    """

    def inlets(self):
        return ['in1']

    def component(self):
        return 'sink'

# %%


class turbomachine(component):
    """
    component turbomachine
        - can be subdivided in pump, compressor and turbine

    **available parameters**
        - P: power
        - eta_s: isentropic efficiency
        - dp: inlet to outlet pressure ratio
        - char: characteristic curve to use, characteristics are generated in
                preprocessing of offdesign calculations

    **design parameters**
        - dp, eta_s

    **offdesign parameters**
        - char

    **inlets and outlets**
        - in1
        - out1
    """

    def attr(self):
        return (component.attr(self) + ['P', 'eta_s', 'dp', 'char'])

    def design(self):
        return ['dp', 'eta_s']

    def offdesign(self):
        return ['char']

    def inlets(self):
        return ['in1']

    def outlets(self):
        return ['out1']

    def component(self):
        return 'turbomachine'

    def equations(self, nw):
        r"""
        returns vector vec_res with result of equations for this component

        :param nw: network using this component object
        :type nw: pyte.networks.network object

        **equations**

        .. math::
            0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in} \right) - P\\
            0 = dp \cdot p_{in} - p_{out}

        for isentropic efficiency, characteristics and additional equations
        see the subclasses:
            - turbine
            - compressor
            - pump
        """

        vec_res = []
        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())

        vec_res += self.fluid_res(inlets, outlets)
        vec_res += self.mass_flow_res(inlets, outlets)

        if self.P_set:
            vec_res += [inlets[0].m * (outlets[0].h - inlets[0].h) - self.P]

        if self.dp_set:
            vec_res += [self.dp * inlets[0].p - outlets[0].p]

        if self.eta_s_set:
            self.eta_s_res = self.eta_s_func(inlets, outlets)
            vec_res += [self.eta_s_res]

        if self.char_set:
            vec_res += self.char_func(inlets, outlets).tolist()

        vec_res += self.additional_equations(nw)

        return vec_res

    def additional_equations(self, nw):
        """
        returns vector vec_res with result of additional equations for this
        component

        :param nw: network using this component object
        :type nw: pyte.networks.network
        :returns: vec_res (*list*) - vector of residual values
        """
        return []

    def derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition

        :param nw: network using this component object
        :type nw: pyte.networks.network
        :returns: mat_deriv (*numpy array*) - matrix of partial derivatives

        **example**

        matrix of partial derivatives for a turbine with one fluid in fluid
        vector and specified value for power P

        .. math::
            \left(
            \begin{array}{cccc}
                0 & 0 & 0 & 1\\
                0 & 0 & 0 & -1\\
                1 & 0 & 0 & 0\\
                -1 & 0 & 0 & 0\\
                h_{out} - h_{in} & 0 & -\dot{m}_{in} & 0\\
                0 & 0 & \dot{m}_{in} & 0
            \end{array}
            \right)

        .. note::
            in this example you can see, that there is no equation regarding
            pressure change, thus the pressure at the inlet and the outlet must
            be defined externally through other components or by connection
            parametrisation
        """
        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inlets), len(outlets)
        num_fl = len(nw.fluids)
        mat_deriv = []

        mat_deriv += self.fluid_deriv(inlets, outlets)
        mat_deriv += self.mass_flow_deriv(inlets, outlets)

        if self.P_set:
            P_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
            for k in range(num_i + num_o - 1):
                P_deriv[k, 0, 0] = outlets[0].h - inlets[0].h
                P_deriv[k, 0, 2] = -inlets[0].m
                P_deriv[k, k + 1, 2] = inlets[0].m
            mat_deriv += P_deriv.tolist()

        if self.dp_set:
            dp_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
            for k in range(num_i + num_o - 1):
                dp_deriv[k, 0, 1] = self.dp
                dp_deriv[k, k + 1, 1] = -1
            mat_deriv += dp_deriv.tolist()

        if self.eta_s_set:
            mat_deriv += self.eta_s_deriv(inlets, outlets)

        if self.char_set:
            mat_deriv += self.char_deriv(inlets, outlets)

        mat_deriv += self.additional_derivatives(nw)

        return np.asarray(mat_deriv)

    def additional_derivatives(self, nw):
        """
        returns matrix mat_deriv with partial derivatives for additional
        equations of this component

        :param nw: network using this component object
        :type nw: pyte.networks.network
        :returns: mat_deriv (*numpy array*) - matrix of partial derivatives
        """
        return []

    def eta_s_func(self, inlets, outlets):
        """
        see subclasses
        """
        msg = ('If you want to use eta_s as parameter, '
               'please specify which type of turbomachine you are using.')
        raise MyComponentError(msg)

    def h_os(self, inlets, outlets):
        """
        calculates the enthalpy at the outlet if compression or expansion is
        isentropic

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: h (*float*) - enthalpy after isentropic state change
        """
        if isinstance(inlets[0], float) or isinstance(inlets[0], int):
            i = inlets
            o = outlets
        else:
            i = inlets[0].as_list()
            o = outlets[0].as_list()
        if num_fluids(i[3]) == 1:
            for fluid, x in i[3].items():
                if x > err:
                    return CPPSI('H', 'P', o[1], 'S',
                                 CPPSI('S', 'P', i[1], 'H', i[2], fluid),
                                 fluid)
        else:
            T_mix = T_mix_ph(i)
            s_mix = s_mix_pT(i, T_mix)
            T_mix_o = T_mix_ps(o, s_mix)
            return h_mix_pT(o, T_mix_o)

    def char_func(self, inlets, outlets):
        raise MyComponentError('Function not available for this component.')

    def eta_s_deriv(self, inlets, outlets):
        """
        calculates partial derivatives of the isentropic efficiency function
            - if the residual value for this equation is lower than the square
              value of the global error tolerance skip calculation
            - calculates the partial derivatives for enthalpy and pressure at
              inlet and for pressure at outlet numerically
            - partial derivative to enthalpy at outlet can be calculated
              analytically, :code:`-1` for expansion and :code:`-self.eta_s`
              for compression
        """

        num_i, num_o = len(inlets), len(outlets)
        num_fl = len(inlets[0].fluid)
        mat_deriv = np.zeros((1, num_i + num_o, num_fl + 3))

        if abs(self.eta_s_res) > err ** (2):

            for i in range(num_i + num_o):
                mat_deriv[0, i, 1] = self.ddx_func(inlets, outlets,
                                                   self.eta_s_func, 'p', i)
                if i == 0:
                    mat_deriv[0, i, 2] = self.ddx_func(inlets, outlets,
                                                       self.eta_s_func, 'h', i)
                else:
                    if isinstance(self, turbine):
                        mat_deriv[0, i, 2] = -1
                    else:
                        mat_deriv[0, i, 2] = -self.eta_s

        else:
            for i in range(num_i + num_o):
                mat_deriv[0, i, 1] = -1
                mat_deriv[0, i, 2] = -1

        return mat_deriv.tolist()

    def char_deriv(self, inlets, outlets):
        raise MyComponentError('Function not available.')

    def convergence_check(self, nw):
        """
        performs a convergence check (class turbine overwrites this method)
            - check if isentropic efficiency or characteristic is set
            - manipulate enthalpies at inlet and outlet if not specified by
              user, if function for isentropic efficiency cannot be calculated
        """

        i, o = nw.comps.loc[self].i, nw.comps.loc[self].o

        if not self.eta_s_set and not self.char_set:
            return

        expr = False
        while not expr:
            try:
                self.eta_s_func(i, o)
                expr = True
            except:
                if not i[0].h_set:
                    i[0].h = i[0].h * 0.8 + 2e4
                if not o[0].h_set:
                    o[0].o = i[0].h + 5e4

    def calc_parameters(self, inlets, outlets, mode):
        """
        parameter calculation pre- or postprocessing

        **postprocessing**
            - calculate power P
            - calculate isentropic efficiency
            - calculate pressure ratio

        **preprocessing**
            - set references for inlet :code:`self.i0` and outlet
              :code:`self.o0` flows
            - set attribute for isentropic enthalpy difference
              :code:`self.dh_s0` at reference
            - generate characteristics for components

        TODO:
            - check if component specific parts are well located in subclasses
        """

        if mode == 'post':

            self.P = inlets[0].m * (outlets[0].h - inlets[0].h)

            if isinstance(self, pump) or isinstance(self, compressor):
                self.eta_s = ((self.h_os(inlets, outlets) - inlets[0].h) /
                              (outlets[0].h - inlets[0].h))
            elif isinstance(self, turbine):
                self.eta_s = ((outlets[0].h - inlets[0].h) /
                              (self.h_os(inlets, outlets) - inlets[0].h))
            else:
                raise MyComponentError('If you want to use eta_s as parameter,'
                                       ' please specify which type of '
                                       'turbomachine you are using.')

            self.dp = outlets[0].p / inlets[0].p

        if mode == 'pre':

            self.i0 = inlets[0].as_list()
            self.o0 = outlets[0].as_list()
            self.dh_s0 = (self.h_os(self.i0, self.o0) - self.i0[2])

            print('Creating characteristics for component ', self)

            if isinstance(self, pump):
                self.char = cmp_char.pump(self.i0[0] * v_mix_ph(self.i0),
                                          self.eta_s)

            if isinstance(self, turbine):
                self.char = cmp_char.turbine(self.eta_s)
                nu_new = np.linspace(self.char.nu[0],
                                     self.char.nu[-1], 1001)
                self.nu0 = nu_new[np.argmax(self.char.eta(nu_new))]

            if isinstance(self, compressor):
                self.char = cmp_char.compressor()

    def print_parameters(self, inlets, outlets):
        i1 = inlets[0].as_list()
        o1 = outlets[0].as_list()
        print('##### ', self.label, ' #####')
        if self.eta_s > 1:
            print('!!!!! Error in parametrisation of the model, '
                  'eta_s higher than 1 !!!!!')
        print('P = ', self.P, 'W; '
              'eta_s = ', self.eta_s, '; '
              'dp = ', self.dp, '; '
              'm = ', inlets[0].m, 'kg / s; '
              'Sirr = ', inlets[0].m * (s_mix_ph(o1) - s_mix_ph(i1)), 'W / K')

        if isinstance(self, compressor) and not isinstance(self.char, int):
            n = math.sqrt(T_mix_ph(self.i0)) / math.sqrt(T_mix_ph(i1))
            m = (
                (i1[0] * math.sqrt(T_mix_ph(i1)) / i1[1]) /
                (self.i0[0] * math.sqrt(T_mix_ph(self.i0)) / self.i0[1])
                )
            vigv = self.char.get_vigv(n, m, (o1[1] *
                                      self.i0[1]) / (i1[1] * self.o0[1]))
            if abs(self.vigv - vigv) > err:
                print('!!!!! Selected inlet guide vane angle is not feasible '
                      '!!!!!')
                if self.vigv > vigv:
                    print('calculated maximum angle:', vigv,
                          'selected:', self.vigv)
                else:
                    print('calculated minimum angle:', vigv,
                          'selected:', self.vigv)
            else:
                print('vigv = ', self.vigv)

# %%


class pump(turbomachine):
    """
    component pump
        - for parametetrisation see parent class turbomachine
    """
    def component(self):
        return 'pump'

    def eta_s_func(self, inlets, outlets):
        r"""
        equation for isentropic efficiency of a pump

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*float*) - residual value of equation

        .. math::
            0 = -\left( h_{out} - h_{in} \right) \cdot \eta_{s,c} +
            \left( h_{out,s} -  h_{in} \right)
        """
        return (-(outlets[0].h - inlets[0].h) * self.eta_s +
                (self.h_os(inlets, outlets) - inlets[0].h))

    def char_func(self, inlets, outlets):
        r"""
        equation for characteristics of a pump

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*numpy array*) - residual value of equation

        .. math::
            0 = -\left( h_{out} - h_{in} \right) \cdot char\left( \dot{m}_{in}
            \cdot v_{in} \right) + \left( h_{out,s} - h_{in} \right)
        """
        i = inlets[0].as_list()
        o = outlets[0].as_list()
        return np.array([(-(o[2] - i[2]) * self.char.eta(
                         i[0] * v_mix_ph(i)) + (self.h_os(i, o) - i[2]))])

    def char_deriv(self, inlets, outlets):
        r"""
        calculates the derivatives for the characteristics

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: mat_deriv (*list*) - matrix of derivatives

        **example**
        one fluid in fluid vector

        .. math::
            \left(
            \begin{array}{cccc}
                \frac{\partial char}{\partial \dot{m}_{in}} &
                \frac{\partial char}{\partial p_{in} &
                \frac{\partial char}{\partial h_{in} & 0\\
                0 & \frac{\partial char}{\partial p_{out} &
                \frac{\partial char}{\partial h_{out} & 0\\
            \end{array}
            \right)
        """
        num_i, num_o = len(inlets), len(outlets)
        num_fl = len(inlets[0].fluid)
        mat_deriv = np.zeros((1, num_i + num_o, num_fl + 3))

        mat_deriv[0, 0, 0] = (
            self.ddx_func(inlets, outlets, self.char_func, 'm', 0))
        for i in range(2):
            mat_deriv[0, i, 1] = (
                self.ddx_func(inlets, outlets, self.char_func, 'p', i))
            mat_deriv[0, i, 2] = (
                self.ddx_func(inlets, outlets, self.char_func, 'h', i))

        return mat_deriv.tolist()

    def initialise_source_p(self, c):
        return 10e5

    def initialise_target_p(self, c):
        return 1e5

    def initialise_source_h(self, c):
        return 3e5

    def initialise_target_h(self, c):
        return 2.9e5

# %%


class compressor(turbomachine):
    """
    component compressor
        - for parametetrisation see parent class turbomachine

    **additional parameter**
        - vigv: variable inlet guide vane angle
    """
    def component(self):
        return 'compressor'

    def attr(self):
        return turbomachine.attr(self) + ['vigv']

    def eta_s_func(self, inlets, outlets):
        r"""
        equation for isentropic efficiency of a compressor is the same as for a
        pump, see documentation of component pump
        """
        return (-(outlets[0].h - inlets[0].h) * self.eta_s +
                (self.h_os(inlets, outlets) - inlets[0].h))

    def char_func(self, inlets, outlets):
        r"""
        equation(s) for characteristics of compressor
            - returns one value, if vigv is not set
            - returns two values, if vigv is set

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*numpy array*) - residual value(s) of equation(s):

            - :code:`np.array([val1, val2])` if vigv_set
            - :code:`np.array([val2])` else

        .. math::

            n = \sqrt{\frac{T_{in,ref}}{T_{in}}}\\
            m = \frac{\dot{m}_{in} \cdot \sqrt{T_{in}} \cdot p_{in,ref}}
            {\dot{m}_{in,ref} \cdot \sqrt{T_{in,ref}} \cdot p_{in}}\\
            val_1 = \frac{p_{out} \cdot p_{in,ref}}{p_{in} \cdot p_{out,ref}} -
            pr_{c}(char(m))\\
            val_2 = \frac{\eta_{s,c}}{\eta_{s,c,ref}} - \eta_{s,c}(char(m))

        **parameters**
            - n: speedline index (rotational speed is constant)
            - m: nondimensional mass flow
            - val1: change ratio to reference case in mass flow and pressure
            - val2: change of isentropic efficiency to reference case

        **logic**
            - calculate n
            - calculate m
            - calculate dn for convergence stability reasons (move speedline
              inside of feasible range of compressor map)
        **if vigv is set**
            - calculate possible vigv range and adjust user specified vigv
              angle, if not inside feasible range (throws warning in post-
              processing)
            - create new speedline to look up values for val1 and val2
            - val1, val2 are relative factors for pressure ratio and isentropic
              efficiency
        **else**
            - set vigv (from compressor map with pressure ratio)
            - calculate relative factor for isentropic efficiency

        TODO:
            - parameter memorisation for improved performance
        """
        if isinstance(inlets[0], float):
            i = inlets
            o = outlets
        else:
            i = inlets[0].as_list()
            o = outlets[0].as_list()
        n = math.sqrt(T_mix_ph(self.i0)) / math.sqrt(T_mix_ph(i))
        m = (i[0] * math.sqrt(T_mix_ph(i)) * self.i0[1] /
             (self.i0[0] * math.sqrt(T_mix_ph(self.i0)) * i[1]))

        dn = 0

        if n < min(self.char.pr.keys()):
            dn = min(self.char.pr.keys()) - n
        if n > max(self.char.pr.keys()):
            dn = max(self.char.pr.keys()) - n

        if self.vigv_set:

            vigv_range = self.char.get_vigv_range(n + dn, m)

            dvigv = 0
            if self.vigv < vigv_range[0]:
                dvigv = vigv_range[0] - self.vigv + 1
            if self.vigv > vigv_range[1]:
                dvigv = vigv_range[1] - self.vigv - 1

            speedline = self.char.get_speedline(n + dn, self.vigv + dvigv)

            return np.array([
                    o[1] * self.i0[1] / (i[1] * self.o0[1]) - speedline[0](m),
                    ((self.h_os(i, o) - i[2]) / (o[2] - i[2])) /
                    (self.dh_s0 / (self.o0[2] - self.i0[2])) -
                    speedline[1](m)
                ])

        else:

            self.vigv = self.char.get_vigv(
                n + dn, m, o[1] / i[1] / (self.o0[1] / self.i0[1]))

            return np.array([
                    ((self.h_os(i, o) - i[2]) / (o[2] - i[2])) /
                    (self.dh_s0 / (self.o0[2] - self.i0[2])) -
                    self.char.get_eta(n + dn, m, self.vigv)
                ])

    def char_deriv(self, inlets, outlets):
        r"""
        calculates the derivatives for the characteristics
            - if vigv is set two sets of equations are used

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: mat_deriv (*list*) - matrix of derivatives

        **example**

        see method char_deriv of class pump for an example

        TODO:
            - improve asthetics, this part of code looks horrible
        """
        num_i = len(inlets)
        num_o = len(outlets)
        num_fl = len(inlets[0].fluid.keys())

        m11 = self.ddx_func(inlets, outlets, self.char_func, 'm', 0)
        p11 = self.ddx_func(inlets, outlets, self.char_func, 'p', 0)
        h11 = self.ddx_func(inlets, outlets, self.char_func, 'h', 0)

        p21 = self.ddx_func(inlets, outlets, self.char_func, 'p', 1)
        h21 = self.ddx_func(inlets, outlets, self.char_func, 'h', 1)

        if self.vigv_set:
            deriv = np.zeros((2, num_i + num_o, num_fl + 3))
            deriv[0, 0, 0] = m11[0]
            deriv[0, 0, 1] = p11[0]
            deriv[0, 0, 2] = h11[0]
            deriv[0, 1, 1] = p21[0]
            deriv[0, 1, 2] = h21[0]
            deriv[1, 0, 0] = m11[1]
            deriv[1, 0, 1] = p11[1]
            deriv[1, 0, 2] = h11[1]
            deriv[1, 1, 1] = p21[1]
            deriv[1, 1, 2] = h21[1]
            return deriv.tolist()
        else:
            deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            deriv[0, 0, 0] = m11[0]
            deriv[0, 0, 1] = p11[0]
            deriv[0, 0, 2] = h11[0]
            deriv[0, 1, 1] = p21[0]
            deriv[0, 1, 2] = h21[0]
            return deriv.tolist()

    def initialise_source_p(self, c):
        return 10e5

    def initialise_target_p(self, c):
        return 1e5

    def initialise_source_h(self, c):
        return 6e5

    def initialise_target_h(self, c):
        return 4e5


# %%

class turbine(turbomachine):
    def component(self):
        return 'turbine'

    def attr(self):
        return turbomachine.attr(self) + ['cone']

    def offdesign(self):
        return turbomachine.offdesign(self) + ['cone']

    def additional_equations(self, nw):
        vec_res = []
        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())

        if self.cone_set:
            vec_res += [self.cone_func(inlets, outlets)]

        return vec_res

    def additional_derivatives(self, nw):

        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inlets), len(outlets)
        num_fl = len(nw.fluids)
        mat_deriv = []

        if self.cone_set:

            cone_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            cone_deriv[0, 0, 0] = -1
            cone_deriv[0, 0, 1] = (
                self.ddx_func(inlets, outlets, self.cone_func, 'p', 0))
            cone_deriv[0, 0, 2] = (
                self.ddx_func(inlets, outlets, self.cone_func, 'h', 0))
            cone_deriv[0, 1, 2] = (
                self.ddx_func(inlets, outlets, self.cone_func, 'p', 1))
            mat_deriv += cone_deriv.tolist()

        return mat_deriv

    def eta_s_func(self, inlets, outlets):
        return (-(outlets[0].h - inlets[0].h) +
                (self.h_os(inlets, outlets) - inlets[0].h) * self.eta_s)

    def cone_func(self, inlets, outlets):
        """
        Formulation with (p0 * v0) / (p * v)
        """
        i = inlets[0].as_list()
        o = outlets[0].as_list()
        n = 1
        return (self.i0[0] * i[1] / self.i0[1] * math.sqrt(
                self.i0[1] * v_mix_ph(self.i0) / (i[1] * v_mix_ph(i))) *
                math.sqrt(abs((1 - (o[1] / i[1]) ** ((n + 1) / n)) /
                (1 - (self.o0[1] / self.i0[1]) ** ((n + 1) / n)))) - i[0])

        """
        Formulation with T0 / T
        """
#        n = 1
#        return (
#            self.i1_0[0] * i1[1] / self.i1_0[1] * math.sqrt(
#                T_mix_ph(self.i1_0) / T_mix_ph(i1)) *
#            math.sqrt((1 - (o1[1] / i1[1]) ** ((n + 1) / n)) /
#                (1 - (self.o1_0[1] / self.i1_0[1]) ** ((n + 1) / n))) - i1[0])

    def char_func(self, inlets, outlets):
        i = inlets[0].as_list()
        o = outlets[0].as_list()
        expr = (math.sqrt(abs(self.dh_s0)) /
                math.sqrt(abs(i[2] - self.h_os(i, o))) + self.nu0 - 1)
        if expr < 0:
            expr = 0.2
        if expr > 1:
            expr = 0.8
        return (
            np.array([(-(o[2] - i[2]) +
                      (self.o0[2] - self.i0[2]) / self.dh_s0 *
                      self.char.eta(expr) *
                      (self.h_os(i, o) - i[2]))]))

    def char_deriv(self, inlets, outlets):

        num_i, num_o = len(inlets), len(outlets)
        num_fl = len(inlets[0].fluid)

        mat_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
        for i in range(num_i + num_o):
            mat_deriv[0, i, 1] = (
                self.ddx_func(inlets, outlets, self.char_func, 'p', i))
            mat_deriv[0, i, 2] = (
                self.ddx_func(inlets, outlets, self.char_func, 'h', i))

        return mat_deriv.tolist()

    def convergence_check(self, nw):
        i, o = nw.comps.loc[self].i, nw.comps.loc[self].o

        if i[0].p < o[0].p:
            o[0].p = i[0].p / 2
        if i[0].h <= o[0].h:
            o[0].h = i[0].h * 0.9

    def initialise_source_p(self, c):
        return 0.5e5

    def initialise_target_p(self, c):
        return 25e5

    def initialise_source_h(self, c):
        return 15e5

    def initialise_target_h(self, c):
        return 20e5

# %%


class split(component):
    """
    component split
    has one inlet and two outlets
    TODO: demixture (separator) requires power and cooling, add equations
    """

    def inlets(self):
        return ['in1']

    def outlets(self):
        return ['out1', 'out2']

    def component(self):
        return 'split'

    def equations(self, nw):

        """
        equations for splitter
        """
        vec_res = []
        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())

        if isinstance(self, splitter):
            vec_res += self.fluid_res(inlets, outlets)

            for o in outlets:
                vec_res += [inlets[0].h - o.h]

        """
        equations for seperator
        TODO: fluid separation requires power!
        """
        if isinstance(self, separator):
            vec_res += self.fluid_res(inlets, outlets)

            if num_fluids(inlets[0].fluid) <= 1:
                for o in outlets:
                    vec_res += [inlets[0].h - o.h]
            else:
                for o in outlets:
                    vec_res += [T_mix_ph(inlets[0].as_list()) -
                                T_mix_ph(o.as_list())]

        vec_res += self.mass_flow_res(inlets, outlets)

        for o in outlets:
            vec_res += [inlets[0].p - o.p]

        return vec_res

    def derivatives(self, nw):

        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inlets), len(outlets)
        num_fl = len(nw.fluids)
        mat_deriv = []

        """
        derivatives for splitter
        """
        if isinstance(self, splitter):
            mat_deriv += self.fluid_deriv(inlets, outlets)

            h_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
            k = 0
            for o in outlets:
                h_deriv[k, 0, 2] = 1
                h_deriv[k, k + 1, 2] = -1
                k += 1

            mat_deriv += h_deriv.tolist()
        """
        derivatives for seperator
        """
        if isinstance(self, separator):
            mat_deriv = np.vstack([mat_deriv,
                                   self.fluid_deriv(inlets, outlets)])

            if num_fluids(inlets[0].fluid) <= 1:

                h_deriv = np.zeros((num_i + num_o - 1, num_i + num_o,
                                    num_fl + 3))
                k = 0
                for o in outlets:
                    h_deriv[k, 0, 2] = 1
                    h_deriv[k, k + 1, 2] = -1
                    k += 1

                mat_deriv += h_deriv.tolist()
            else:

                T_deriv = np.zeros((num_i + num_o - 1, num_i + num_o,
                                    num_fl + 3))
                k = 0
                for o in outlets:
                    T_deriv[k, 0, 1] = dT_mix_dph(inlets[0].as_list())
                    T_deriv[k, 0, 2] = dT_mix_pdh(inlets[0].as_list())
                    T_deriv[k, 0, 3:] = dT_mix_ph_dfluid(inlets[0].as_list())
                    T_deriv[k, k + 1, 1] = -dT_mix_dph(o.as_list())
                    T_deriv[k, k + 1, 2] = -dT_mix_pdh(o.as_list())
                    T_deriv[k, k + 1, 3:] = -1 * dT_mix_ph_dfluid(o.as_list())
                    k += 1

                mat_deriv += T_deriv.tolist()

        mat_deriv += self.mass_flow_deriv(inlets, outlets)

        p_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
        k = 0
        for o in outlets:
            p_deriv[k, 0, 1] = 1
            p_deriv[k, k + 1, 1] = -1
            k += 1
        mat_deriv += p_deriv.tolist()

        return np.asarray(mat_deriv)

    def initialise_source_p(self, c):
        return 1e5

    def initialise_target_p(self, c):
        return 1e5

    def initialise_source_h(self, c):
        return 5e5

    def initialise_target_h(self, c):
        return 5e5

    def print_parameters(self, inlets, outlets):

        print('##### ', self.label, ' #####')
        print('m_in = ', inlets[0].m, 'kg / s; ')
        i = 1
        for o in outlets:
            print('m_out' + str(i) + ' = ', o.m, 'kg / s; ')
            i += 1
        if isinstance(self, separator):
            print('; fluid_in:', inlets[0].fluid, '; ')
            i = 1
            for o in outlets:
                print('fluid_out' + str(i) + ' = ', o.fluid, 'kg / s; ')
                i += 1


class splitter(split):

    def component(self):
        return 'splitter'


class separator(split):

    def component(self):
        return 'separator'

# %%


class merge(component):
    """
    component merge
    has two inlets and one outlet
    """

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1']

    def component(self):
        return 'merge'

    def equations(self, nw):
        vec_res = []
        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())

        vec_res += self.fluid_res(inlets, outlets)
        vec_res += self.mass_flow_res(inlets, outlets)

        h_res = -outlets[0].m * outlets[0].h
        for i in inlets:
            h_res += i.m * i.h
        vec_res += [h_res]

        for i in inlets:
            vec_res += [outlets[0].p - i.p]

        return vec_res

    def derivatives(self, nw):

        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inlets), len(outlets)
        num_fl = len(nw.fluids)
        mat_deriv = []

        mat_deriv += self.fluid_deriv(inlets, outlets)
        mat_deriv += self.mass_flow_deriv(inlets, outlets)

        h_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
        h_deriv[0, num_i, 0] = -outlets[0].h
        h_deriv[0, num_i, 2] = -outlets[0].m
        k = 0
        for i in inlets:
            h_deriv[0, k, 0] = i.h
            h_deriv[0, k, 2] = i.m
            k += 1
        mat_deriv += h_deriv.tolist()

        p_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
        k = 0
        for i in inlets:
            p_deriv[k, k, 1] = -1
            p_deriv[k, num_i, 1] = 1
            k += 1
        mat_deriv += p_deriv.tolist()

        return np.asarray(mat_deriv)

    def initialise_fluids(self, nw):

        for outconn in nw.comps.loc[self].o:
            inlets = nw.comps.loc[self].i.tolist()
            for fluid in nw.fluids:
                if not outconn.fluid_set[fluid]:
                    x = 0
                    m = 0
                    for inlet in inlets:
                        m += inlet.m
                        x += inlet.fluid[fluid] * inlet.m

                    outconn.fluid[fluid] = x / m

    def initialise_source_p(self, c):
        return 1e5

    def initialise_target_p(self, c):
        return 1e5

    def initialise_source_h(self, c):
        return 5e5

    def initialise_target_h(self, c):
        return 5e5

    def print_parameters(self, inlets, outlets):

        print('##### ', self.label, ' #####')
        j = 1
        for i in inlets:
            print('m_in' + str(j) + ' = ', i.m, 'kg / s; ')
            j += 1
        print('m_out = ', outlets[0].m, 'kg / s; ')

# %%


class combustion_chamber(component):
    """
    Component combustion chamber.
    TODO:
    """
    def comp_init(self, nw):
        self.o2 = [x for x in nw.fluids if x in
                   [a.replace(' ', '') for a in CP.get_aliases('O2')]][0]
        self.co2 = [x for x in nw.fluids if x in
                    [a.replace(' ', '') for a in CP.get_aliases('CO2')]][0]
        self.h2o = [x for x in nw.fluids if x in
                    [a.replace(' ', '') for a in CP.get_aliases('H2O')]][0]
        self.n2 = [x for x in nw.fluids if x in
                   [a.replace(' ', '') for a in CP.get_aliases('N2')]][0]

        self.n_c = fluid_structure(self.fuel)['C']
        self.n_h = fluid_structure(self.fuel)['H']

        self.hi = self.lhv()

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1']

    def attr(self):
        return component.attr(self) + ['Q', 'fuel', 'lamb']

    def fuels(self):
        return ['methane']

    def component(self):
        return 'combustion chamber'

    def equations(self, nw):

        vec_res = []
        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())

        if not self.fuel_set:
            raise MyComponentError('Must specify fuel for combustion chamber.')

        for fluid in inlets[0].fluid.keys():
            vec_res += [self.reaction_balance(inlets, outlets, fluid)]

        vec_res += self.mass_flow_res(inlets, outlets)

        for i in inlets:
            vec_res += [outlets[0].p - i.p]

        self.eb_res = self.energy_balance(inlets, outlets)
        vec_res += [self.eb_res]

        if self.Q_set:
            val = 0
            for i in inlets:
                val += i.m * i.fluid[self.fuel] * self.hi
            vec_res += [val - self.Q]

        return vec_res

    def derivatives(self, nw):

        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inlets), len(outlets)
        num_fl = len(nw.fluids)
        mat_deriv = []

        j = 0
        fl_deriv = np.zeros((num_fl, num_i + num_o, num_fl + 3))
        for fluid in inlets[0].fluid.keys():
            for i in range(num_i + num_o):
                fl_deriv[j, i, 0] = self.drb_dx(inlets, outlets, 'm', i, fluid)
                fl_deriv[j, i, 3:] = (
                    self.drb_dx(inlets, outlets, 'fluid', i, fluid))

            j += 1

        mat_deriv += fl_deriv.tolist()

        mat_deriv += self.mass_flow_deriv(inlets, outlets)

        p_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
        for k in range(num_i + num_o - 1):
            p_deriv[k][num_i][1] = 1
            p_deriv[k][k][1] = -1

        mat_deriv += p_deriv.tolist()

        eb_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
        for i in range(num_i + num_o):
            eb_deriv[0, i, 0] = (
                self.ddx_func(inlets, outlets, self.energy_balance, 'm', i))
            eb_deriv[0, i, 1] = (
                self.ddx_func(inlets, outlets, self.energy_balance, 'p', i))
            if i >= num_i:
                eb_deriv[0, i, 2] = -(inlets + outlets)[i].m
            else:
                eb_deriv[0, i, 2] = (inlets + outlets)[i].m

        mat_deriv += eb_deriv.tolist()

        if self.Q_set:
            Q_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            k = 0
            for i in inlets:
                Q_deriv[0, k, 0] = i.fluid[self.fuel] * self.hi
                k += 1

            mat_deriv += Q_deriv.tolist()

        return np.asarray(mat_deriv)

    def lhv(self):
        val = 50.015e6
        return val

    def reaction_balance(self, inlets, outlets, fluid):

        n_fuel = 0
        for i in inlets:
            n_fuel += (i.m * i.fluid[self.fuel] /
                       molar_masses[self.fuel])

        n_oxygen = 0
        for i in inlets:
            n_oxygen += (i.m * i.fluid[self.o2] /
                         molar_masses[self.o2])

        if not self.lamb_set:
            self.lamb = n_oxygen / (n_fuel * (self.n_c + self.n_h / 4))

        n_fuel_exc = 0
        if self.lamb < 1:
            n_fuel_exc = n_fuel - n_oxygen / (self.n_c + self.n_h / 4)

        if fluid == self.co2:
            dm = ((n_fuel - n_fuel_exc) *
                  self.n_c * molar_masses[self.co2])
        elif fluid == self.h2o:
            dm = ((n_fuel - n_fuel_exc) *
                  self.n_h / 2 * molar_masses[self.h2o])
        elif fluid == self.o2:
            if self.lamb < 1:
                dm = -n_oxygen * molar_masses[self.o2]
            else:
                dm = -n_oxygen / self.lamb * molar_masses[self.o2]
        elif fluid == self.fuel:
            dm = -(n_fuel - n_fuel_exc) * molar_masses[self.fuel]
        else:
            dm = 0

        res = dm

        for i in inlets:
            res += i.fluid[fluid] * i.m
        for o in outlets:
            res -= o.fluid[fluid] * o.m
        return res

    def energy_balance(self, inlets, outlets):
        T_ref = 500
        p_ref = 1e5

        res = 0
        for i in inlets:
            res += i.m * (i.h - h_mix_pT([i.m, p_ref, i.h, i.fluid], T_ref))
            res += i.m * i.fluid[self.fuel] * self.hi
        for o in outlets:
            res -= o.m * (o.h - h_mix_pT([o.m, p_ref, o.h, o.fluid], T_ref))
            res -= o.m * o.fluid[self.fuel] * self.hi

        return res

    def drb_dx(self, inlets, outlets, dx, pos, fluid):

        dm, dp, dh, df = 0, 0, 0, 0
        if dx == 'm':
            dm = 1e-5
        elif dx == 'p':
            dp = 1e-5
        elif dx == 'h':
            dh = 1e-5
        else:
            df = 1e-5

        if dx == 'fluid':
            deriv = []
            for f in inlets[0].fluid.keys():
                val = (inlets + outlets)[pos].fluid[f]
                exp = 0
                if (inlets + outlets)[pos].fluid[f] + df <= 1:
                    (inlets + outlets)[pos].fluid[f] += df
                else:
                    (inlets + outlets)[pos].fluid[f] = 1
                exp += self.reaction_balance(inlets, outlets, fluid)
                if (inlets + outlets)[pos].fluid[f] - 2 * df >= 0:
                    (inlets + outlets)[pos].fluid[f] -= 2 * df
                else:
                    (inlets + outlets)[pos].fluid[f] = 0
                exp -= self.reaction_balance(inlets, outlets, fluid)
                (inlets + outlets)[pos].fluid[f] = val

                deriv += [exp / (2 * (dm + dp + dh + df))]

        else:
            exp = 0
            (inlets + outlets)[pos].m += dm
            (inlets + outlets)[pos].p += dp
            (inlets + outlets)[pos].h += dh
            exp += self.reaction_balance(inlets, outlets, fluid)

            (inlets + outlets)[pos].m -= 2 * dm
            (inlets + outlets)[pos].p -= 2 * dp
            (inlets + outlets)[pos].h -= 2 * dh
            exp -= self.reaction_balance(inlets, outlets, fluid)
            deriv = exp / (2 * (dm + dp + dh + df))

            (inlets + outlets)[pos].m += dm
            (inlets + outlets)[pos].p += dp
            (inlets + outlets)[pos].h += dh

        return deriv

    def initialise_fluids(self, nw):

        N2 = 0.7655
        O2 = 0.2345

        n_fuel = 1
        lamb = 2
        m_o2 = (n_fuel * (self.n_c + self.n_h / 4) *
                molar_masses[self.o2] * lamb)
        m_co2 = n_fuel * self.n_c * molar_masses[self.co2]
        m_h2o = n_fuel * self.n_h / 2 * molar_masses[self.h2o]
        m_n2 = m_o2 / O2 * N2
        m_o2 /= lamb
        m = m_n2 + m_h2o + m_co2 + m_o2

        fg = {
            self.n2: m_n2 / m,
            self.co2: m_co2 / m,
            self.o2: m_o2 / m,
            self.h2o: m_h2o / m
        }

        for c in nw.comps.loc[self].o:
            for fluid, x in c.fluid.items():
                if not c.fluid_set[fluid] and fluid in fg.keys():
                    c.fluid[fluid] = fg[fluid]

    def convergence_check(self, nw):
        if nw.init_file is not None or nw.iter > 5:
            return

        for o in nw.comps.loc[self].o:
            fluids = [f for f in o.fluid.keys() if not o.fluid_set[f]]
            for f in fluids:
                if f == self.o2:
                    if o.fluid[f] > 0.25:
                        o.fluid[f] = 0.2
                elif f == self.n2:
                    if o.fluid[f] < 0.6 or o.fluid[f] > 0.8:
                        o.fluid[f] = 0.65
                elif f == self.co2:
                    if o.fluid[f] > 0.15:
                        o.fluid[f] = 0.10
                elif f == self.h2o:
                    if o.fluid[f] > 0.15:
                        o.fluid[f] = 0.10
                elif f == self.fuel:
                    if o.fluid[f] > 0.1:
                        o.fluid[f] = 0.05
                else:
                    if o.fluid[f] > 0.05:
                        o.fluid[f] = 0.02

        for c in nw.comps.loc[self].o:
            init_target(nw, c, c.t)

    def calc_parameters(self, inlets, outlets, mode):
        self.Q = 0
        for i in inlets:
            self.Q += i.m * i.fluid[self.fuel] * self.hi

    def print_parameters(self, inlets, outlets):
        print('##### ', self.label, ' #####')
        print('Q = ', self.Q,
              'lambda = ', self.lamb)
        j = 1
        for i in inlets:
            print('m_in' + str(j) + ' = ', i.m, 'kg / s; ')
            j += 1
        print('m_out = ', outlets[0].m, 'kg / s; ')

    def initialise_source_p(self, c):
        return 5e5

    def initialise_target_p(self, c):
        return 5e5

    def initialise_source_h(self, c):
        return 1e6

    def initialise_target_h(self, c):
        return 5e5

# %%


class vessel(component):
    """
    component Vessel
    state difference - pressure dp: 0 < dp <= 1
    TODO: add zeta
    """

    def attr(self):
        return component.attr(self) + ['dp', 'zeta']

    def design(self):
        return ['dp']

    def offdesign(self):
        return ['zeta']

    def inlets(self):
        return ['in1']

    def outlets(self):
        return ['out1']

    def component(self):
        return 'vessel'

    def equations(self, nw):

        vec_res = []
        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())

        vec_res += self.fluid_res(inlets, outlets)
        vec_res += self.mass_flow_res(inlets, outlets)

        vec_res += [inlets[0].h - outlets[0].h]

        if self.dp_set:
            vec_res += [inlets[0].p * self.dp - outlets[0].p]

        if self.zeta_set:
            vec_res += [self.zeta_func(inlets, outlets)]

        return vec_res

    def derivatives(self, nw):

        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inlets), len(outlets)
        num_fl = len(nw.fluids)
        mat_deriv = []

        mat_deriv += self.fluid_deriv(inlets, outlets)
        mat_deriv += self.mass_flow_deriv(inlets, outlets)

        h_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
        for k in range(num_i + num_o - 1):
            h_deriv[k, 0, 2] = 1
            h_deriv[k, k + 1, 2] = -1
        mat_deriv += h_deriv.tolist()

        if self.dp_set:
            dp_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            dp_deriv[0, 0, 1] = self.dp
            dp_deriv[0, 1, 1] = -1
            mat_deriv += dp_deriv.tolist()

        if self.zeta_set:
            zeta_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            for i in range(2):
                if i == 0:
                    zeta_deriv[0, i, 0] = (
                        self.ddx_func(inlets, outlets, self.zeta_func, 'm', i))
                zeta_deriv[0, i, 1] = (
                    self.ddx_func(inlets, outlets, self.zeta_func, 'p', i))
                zeta_deriv[0, i, 2] = (
                    self.ddx_func(inlets, outlets, self.zeta_func, 'h', i))
            mat_deriv += zeta_deriv.tolist()

        return np.asarray(mat_deriv)

    def calc_parameters(self, inlets, outlets, mode):
        self.dp = outlets[0].p / inlets[0].p
        self.zeta = ((inlets[0].p - outlets[0].p) * math.pi ** 2 /
                     (8 * inlets[0].m ** 2 *
                     (v_mix_ph(inlets[0].as_list()) +
                      v_mix_ph(outlets[0].as_list())) / 2))

    def print_parameters(self, inlets, outlets):
        print('##### ', self.label, ' #####')
        print('dp = ', self.dp,
              'zeta = ', self.zeta, 'kg / m^4 * s '
              'm = ', inlets[0].m, 'kg / s; '
              'Sirr = ', inlets[0].m * (s_mix_ph(outlets[0].as_list()) -
                                        s_mix_ph(inlets[0].as_list())), 'W / K'
              )

    def initialise_source_p(self, c):
        return 4e5

    def initialise_target_p(self, c):
        return 5e5

    def initialise_source_h(self, c):
        return 5e5

    def initialise_target_h(self, c):
        return 5e5

# %%


class heat_exchanger_simple(component):
    """
    simplified heat exchanger.
    """

    def attr(self):
        return component.attr(self) + ['Q', 'dp', 'zeta', 'D', 'L', 'ks']

    def design(self):
        return ['dp']

    def offdesign(self):
        return ['zeta']

    def inlets(self):
        return ['in1']

    def outlets(self):
        return ['out1']

    def component(self):
        return 'simplified heat exchanger'

    def equations(self, nw):

        vec_res = []
        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())

        vec_res += self.fluid_res(inlets, outlets)
        vec_res += self.mass_flow_res(inlets, outlets)

        if self.Q_set:
            vec_res += [inlets[0].m * (outlets[0].h - inlets[0].h) - self.Q]

        if self.dp_set:
            vec_res += [inlets[0].p * self.dp - outlets[0].p]

        if self.zeta_set:
            vec_res += [self.zeta_func(inlets, outlets)]

        if self.ks_set and self.D_set and self.L_set:
            vec_res += [self.lamb_func(inlets, outlets)]

        return vec_res

    def derivatives(self, nw):

        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inlets), len(outlets)
        num_fl = len(nw.fluids)
        mat_deriv = []

        mat_deriv += self.fluid_deriv(inlets, outlets)
        mat_deriv += self.mass_flow_deriv(inlets, outlets)

        if self.Q_set:
            Q_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            Q_deriv[0, 0, 0] = outlets[0].h - inlets[0].h
            Q_deriv[0, 0, 2] = -inlets[0].m
            Q_deriv[0, 1, 2] = inlets[0].m
            mat_deriv += Q_deriv.tolist()

        if self.dp_set:
            dp_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            dp_deriv[0, 0, 1] = self.dp
            dp_deriv[0, 1, 1] = -1
            mat_deriv += dp_deriv.tolist()

        if self.zeta_set:
            zeta_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            for i in range(2):
                if i == 0:
                    zeta_deriv[0, i, 0] = (
                        self.ddx_func(inlets, outlets, self.zeta_func, 'm', i))
                zeta_deriv[0, i, 1] = (
                    self.ddx_func(inlets, outlets, self.zeta_func, 'p', i))
                zeta_deriv[0, i, 2] = (
                    self.ddx_func(inlets, outlets, self.zeta_func, 'h', i))
            mat_deriv += zeta_deriv.tolist()

        if self.ks_set and self.D_set and self.L_set:
            lamb_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            for i in range(2):
                if i == 0:
                    lamb_deriv[0, i, 0] = (
                        self.ddx_func(inlets, outlets, self.lamb_func, 'm', i))
                lamb_deriv[0, i, 1] = (
                    self.ddx_func(inlets, outlets, self.lamb_func, 'p', i))
                lamb_deriv[0, i, 2] = (
                    self.ddx_func(inlets, outlets, self.lamb_func, 'h', i))
            mat_deriv += lamb_deriv.tolist()

        return np.asarray(mat_deriv)

    def lamb_func(self, inlets, outlets):
        """
        add calculation for viscosity of ideal gas mixtures
        """
        return (1 / (v_mix_ph(inlets[0].as_list()) +
                     v_mix_ph(outlets[0].as_list())) *
                lamb(re, self.D, self.ks) * self.L / self.D)

    def calc_parameters(self, inlets, outlets, mode):

        self.Q = inlets[0].m * (outlets[0].h - inlets[0].h)
        self.dp = outlets[0].p / inlets[0].p
        self.zeta = ((inlets[0].p - outlets[0].p) * math.pi ** 2 /
                     (8 * inlets[0].m ** 2 *
                     (v_mix_ph(inlets[0].as_list()) +
                      v_mix_ph(outlets[0].as_list())) / 2))

    def print_parameters(self, inlets, outlets):

        print('##### ', self.label, ' #####')
        print('Q = ', self.Q, 'W; '
              'dp = ', self.dp, '; '
              'zeta = ', self.zeta, 'kg / m^4 * s '
              'm = ', inlets[0].m, 'kg / s;'
              'Sq = ', inlets[0].m * (s_mix_ph(outlets[0].as_list()) -
                                      s_mix_ph(inlets[0].as_list())), 'W / K'
              )

    def initialise_source_p(self, c):
        return 1e5

    def initialise_target_p(self, c):
        return 1e5

    def initialise_source_h(self, c):
        if self.Q < 0:
            return 5e5
        elif self.Q > 0:
            return 15e5
        else:
            return 10e5

    def initialise_target_h(self, c):
        if self.Q < 0:
            return 15e5
        elif self.Q > 0:
            return 5e5
        else:
            return 10e5

# %%


class pipe(heat_exchanger_simple):
    """
    Pipe with heatloss and pressure drop.
    """
    def component(self):
        return 'pipe'

# %%
"""
TODO:
    - add components economiser, steam generator and superheater
    - add pressure loss in part load if fluid is condensing or evaporating
    -
"""


class heat_exchanger(component):
    """
    component heat exchanger
    TODO:   adjust calculation of ttd and kA with phase shifts
            check dT_G derivatives - >  fluid derivative necessary?
    """
    def attr(self):
        return (component.attr(self) +
                ['Q', 'kA', 'dT_log', 'ttd_u', 'ttd_l',
                 'dp1', 'dp2', 'zeta1', 'zeta2'])

    """
    hot side: index 1
    cold side: index 2
    """

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1', 'out2']

    def design(self):
        return ['ttd_u', 'ttd_l', 'dp1', 'dp2']

    def offdesign(self):
        return ['kA', 'zeta1', 'zeta2']

    def component(self):
        return 'heat_exchanger'

    def equations(self, nw):

        vec_res = []
        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())

        vec_res += self.fluid_res(inlets, outlets)
        vec_res += self.mass_flow_res(inlets, outlets)

        vec_res += [inlets[0].m * (outlets[0].h - inlets[0].h) +
                    inlets[1].m * (outlets[1].h - inlets[1].h)]

        if self.Q_set:
            vec_res += [inlets[0].m * (outlets[0].h - inlets[0].h) - self.Q]

        if self.kA_set:
            vec_res += [self.kA_func(inlets, outlets)]

        if self.dT_log_set:
            vec_res += [self.dT_log_func(inlets, outlets)]

        if self.ttd_u_set:
            vec_res += [self.ttd_u_func(inlets, outlets)]

        if self.ttd_l_set:
            vec_res += [self.ttd_l_func(inlets, outlets)]

        if self.dp1_set:
            vec_res += [self.dp1 * inlets[0].p - outlets[0].p]

        if self.dp2_set:
            vec_res += [self.dp2 * inlets[1].p - outlets[1].p]

        if self.zeta1_set:
            vec_res += [self.zeta_func(inlets, outlets)]

        if self.zeta2_set:
            vec_res += [self.zeta2_func(inlets, outlets)]

        vec_res += self.additional_equations(nw)

        return vec_res

    def additional_equations(self, nw):
        """
        returns vector for additional equations (specific to components)
        """
        return []

    def derivatives(self, nw):

        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inlets), len(outlets)
        num_fl = len(nw.fluids)
        mat_deriv = []

        mat_deriv += self.fluid_deriv(inlets, outlets)
        mat_deriv += self.mass_flow_deriv(inlets, outlets)

        q_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
        for k in range(num_i):
            q_deriv[0, k, 0] = outlets[k].h - inlets[k].h

            q_deriv[0, k, 2] = -inlets[k].m
        q_deriv[0, 2, 2] = inlets[0].m
        q_deriv[0, 3, 2] = inlets[1].m
        mat_deriv += q_deriv.tolist()

        if self.Q_set:
            Q_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            Q_deriv[0, 0, 0] = outlets[0].h - inlets[0].h
            Q_deriv[0, 0, 2] = -inlets[0].m
            Q_deriv[0, 1, 2] = inlets[0].m
            mat_deriv += Q_deriv.tolist()

        if self.kA_set:
            kA_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            kA_deriv[0, 0, 0] = outlets[0].h - inlets[0].h
            for i in range(num_i + num_o):
                kA_deriv[0, i, 1] = (
                    self.ddx_func(inlets, outlets, self.kA_func, 'p', i))
                kA_deriv[0, i, 2] = (
                    self.ddx_func(inlets, outlets, self.kA_func, 'h', i))
            mat_deriv += kA_deriv.tolist()

        if self.dT_log_set:
            mat_deriv += [[[0,
                       self.ddx_func(i1, i2, o1, o2, self.dT_log_func, 'p11'),
                       self.ddx_func(i1, i2, o1, o2, self.dT_log_func, 'h11')] + z,
                      [0,
                       self.ddx_func(i1, i2, o1, o2, self.dT_log_func, 'p12'),
                       self.ddx_func(i1, i2, o1, o2, self.dT_log_func, 'h12')] + z,
                      [0,
                       self.ddx_func(i1, i2, o1, o2, self.dT_log_func, 'p21'),
                       self.ddx_func(i1, i2, o1, o2, self.dT_log_func, 'h21') + i1[0]] + z,
                      [0,
                       self.ddx_func(i1, i2, o1, o2, self.dT_log_func, 'p22'),
                       self.ddx_func(i1, i2, o1, o2, self.dT_log_func, 'h22')] + z]]

        if self.ttd_u_set:
            mat_deriv += self.ttd_u_deriv(inlets, outlets)

        if self.ttd_l_set:
            mat_deriv += self.ttd_l_deriv(inlets, outlets)

        if self.dp1_set:
            dp1_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            dp1_deriv[0, 0, 1] = self.dp1
            dp1_deriv[0, 2, 1] = -1
            mat_deriv += dp1_deriv.tolist()

        if self.dp2_set:
            dp2_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            dp2_deriv[0, 1, 1] = self.dp2
            dp2_deriv[0, 3, 1] = -1
            mat_deriv += dp2_deriv.tolist()

        if self.zeta1_set:
            zeta1_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            for i in range(2):
                if i == 0:
                    zeta1_deriv[0, i * 2, 0] = (
                        self.ddx_func(inlets, outlets,
                                      self.zeta_func, 'm', i * 2))
                zeta1_deriv[0, i * 2, 1] = (
                    self.ddx_func(inlets, outlets,
                                  self.zeta_func, 'p', i * 2))
                zeta1_deriv[0, i * 2, 2] = (
                    self.ddx_func(inlets, outlets,
                                  self.zeta_func, 'h', i * 2))
            mat_deriv += zeta1_deriv.tolist()

        if self.zeta2_set:
            zeta2_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            for i in range(2):
                if i == 0:
                    zeta2_deriv[0, i * 2 + 1, 0] = (
                        self.ddx_func(inlets, outlets,
                                      self.zeta2_func, 'm', i * 2 + 1))
                zeta2_deriv[0, i * 2 + 1, 1] = (
                    self.ddx_func(inlets, outlets,
                                  self.zeta2_func, 'p', i * 2 + 1))
                zeta2_deriv[0, i * 2 + 1, 2] = (
                    self.ddx_func(inlets, outlets,
                                  self.zeta2_func, 'h', i * 2 + 1))
            mat_deriv += zeta2_deriv.tolist()

        mat_deriv += self.additional_derivatives(nw)

        return np.asarray(mat_deriv)

    def additional_derivatives(self, nw):
        return []

    def zeta2_func(self, inlets, outlets):
        """
        pressure loss in part load
        """
        i = inlets[1].as_list()
        o = outlets[1].as_list()
        return (self.zeta2 - (i[1] - o[1]) * math.pi ** 2 /
                (8 * i[0] ** 2 * (v_mix_ph(i) + v_mix_ph(o)) / 2))

    def kA_func(self, inlets, outlets):

        i1 = inlets[0].as_list()
        i2 = inlets[1].as_list()
        o1 = outlets[0].as_list()
        o2 = outlets[1].as_list()

        T_i1 = T_mix_ph(i1)
        T_i2 = T_mix_ph(i2)
        T_o1 = T_mix_ph(o1)
        T_o2 = T_mix_ph(o2)

        ii1 = 0
        io2 = 0
        while T_i1 <= T_o2:
            try:
                T_o2 = T_mix_ph([o2[0], o2[1], o2[2] - io2 * 10000, o2[3]])
                io2 += 1
            except:
                T_i1 = T_mix_ph([i1[0], i1[1], i1[2] + ii1 * 10000, i1[3]])
                ii1 += 1

        i = 0
        while T_o1 <= T_i2:
            i += 1
            T_o1 = T_mix_ph([o1[0], o1[1], o1[2] + i * 10000, o1[3]])

        return (i1[0] * (o1[2] - i1[2]) + self.kA *
                (T_o1 - T_i2 - T_i1 + T_o2) /
                math.log((T_o1 - T_i2) / (T_i1 - T_o2)))

    def dT_log_func(self, inlets, outlets):

        i1 = inlets[0].as_list()
        i2 = inlets[1].as_list()
        o1 = outlets[0].as_list()
        o2 = outlets[1].as_list()

        T_i1 = T_mix_ph(i1)
        T_i2 = T_mix_ph(i2)
        T_o1 = T_mix_ph(o1)
        T_o2 = T_mix_ph(o2)

        ii1 = 0
        io2 = 0
        while T_i1 <= T_o2:
            try:
                T_o2 = T_mix_ph([o2[0], o2[1], o2[2] - io2 * 10000, o2[3]])
                io2 += 1
            except:
                T_i1 = T_mix_ph([i1[0], i1[1], i1[2] + ii1 * 10000, i1[3]])
                ii1 += 1

        i = 0
        while T_o1 <= T_i2:
            i += 1
            T_o1 = T_mix_ph([o1[0], o1[1], o1[2] + i * 10000, o1[3]])

        return (self.dT_log *
                math.log((T_o1 - T_i2) / (T_i1 - T_o2)) -
                T_o1 + T_i2 + T_i1 - T_o2)

    def ttd_u_func(self, inlets, outlets):
        i1 = inlets[0].as_list()
        o2 = outlets[1].as_list()
        return self.ttd_u - T_mix_ph(i1) + T_mix_ph(o2)

    def ttd_l_func(self, inlets, outlets):
        i2 = inlets[1].as_list()
        o1 = outlets[0].as_list()
        return self.ttd_l - T_mix_ph(o1) + T_mix_ph(i2)

    def ttd_u_deriv(self, inlets, outlets):
        deriv = np.zeros((1, 4, len(inlets[0].fluid) + 3))
        for i in range(2):
            deriv[0, i * 3, 1] = (
                self.ddx_func(inlets, outlets, self.ttd_u_func, 'p', i * 3))
            deriv[0, i * 3, 2] = (
                self.ddx_func(inlets, outlets, self.ttd_u_func, 'h', i * 3))
        return deriv.tolist()

    def ttd_l_deriv(self, inlets, outlets):
        deriv = np.zeros((1, 4, len(inlets[0].fluid) + 3))
        for i in range(2):
            deriv[0, i + 1, 1] = (
                self.ddx_func(inlets, outlets, self.ttd_l_func, 'p', i + 1))
            deriv[0, i + 1, 2] = (
                self.ddx_func(inlets, outlets, self.ttd_l_func, 'h', i + 1))
        return deriv.tolist()

    def convergence_check(self, nw):
        i, o = nw.comps.loc[self].i.tolist(), nw.comps.loc[self].o.tolist()

        if i[0].h < o[0].h and not o[0].h_set:
            o[0].h = i[0].h / 2
        if i[1].h > o[1].h and not i[1].h_set:
            i[1].h = o[1].h / 2

        if self.ttd_u_set:
            expr = False
            while not expr:
                try:
                    self.ttd_u_func(i, o)
                    expr = True
                except:
                    if not i[0].h_set:
                        i[0].h *= 1.05
                    if not o[1].h_set:
                        o[1].h *= 1.05

        if self.ttd_l_set:
            expr = False
            while not expr:
                try:
                    self.ttd_l_func(i, o)
                    expr = True
                except:
                    if not i[1].h_set:
                        i[1].h *= 1.05
                    if not o[0].h_set:
                        o[0].h *= 1.05

    def initialise_source_p(self, c):
        return 50e5

    def initialise_target_p(self, c):
        return 50e5

    def initialise_source_h(self, c):
        flow = [c.m0, c.p0, c.h, c.fluid]
        if c.s_id == 'out1':
            T = 200 + 273.15
            return h_mix_pT(flow, T)
        else:
            T = 250 + 273.15
            return h_mix_pT(flow, T)

    def initialise_target_h(self, c):
        flow = [c.m0, c.p0, c.h, c.fluid]
        if c.t_id == 'in1':
            T = 300 + 273.15
            return h_mix_pT(flow, T)
        else:
            T = 220 + 273.15
            return h_mix_pT(flow, T)

    def calc_parameters(self, inlets, outlets, mode):

        T_i2 = T_mix_ph(inlets[1].as_list())
        T_o1 = T_mix_ph(outlets[0].as_list())

        if isinstance(self, condenser):
            i1 = inlets[0].as_list()
            T_i1 = T_mix_ph([i1[0], i1[1], h_mix_pQ(i1, 1), i1[3]])
        else:
            T_i1 = T_mix_ph(inlets[0].as_list())
        T_o2 = T_mix_ph(outlets[1].as_list())

        self.Q = inlets[0].m * (outlets[0].h - inlets[0].h)
        self.ttd_u = T_i1 - T_o2
        self.ttd_l = T_o1 - T_i2
        if T_i1 <= T_o2 or T_o1 <= T_i2:
            self.dT_log = np.nan
            self.kA = np.nan
        else:
            self.dT_log = ((T_o1 - T_i2 - T_i1 + T_o2) /
                           math.log((T_o1 - T_i2) / (T_i1 - T_o2)))
            self.kA = -self.Q / self.dT_log

        self.dp1 = outlets[0].p / inlets[0].p
        self.dp2 = outlets[1].p / inlets[1].p
        self.zeta1 = ((inlets[0].p - outlets[0].p) * math.pi ** 2 /
                      (8 * inlets[0].m ** 2 *
                      (v_mix_ph(inlets[0].as_list()) +
                       v_mix_ph(outlets[0].as_list())) / 2))
        self.zeta2 = ((inlets[1].p - outlets[1].p) * math.pi ** 2 /
                      (8 * inlets[1].m ** 2 *
                      (v_mix_ph(inlets[1].as_list()) +
                       v_mix_ph(outlets[1].as_list())) / 2))

    def print_parameters(self, inlets, outlets):

        print('##### ', self.label, ' #####')
        if self.ttd_u < 0 and (self.dT_log_set or self.kA_set):
            print('!!!!! ERROR calculating condenser: !!!!!\n'
                  'Negative value for TTD at given logarithmic temperature '
                  'difference or kA, result may be wrong.')
        print('Q = ', self.Q, 'W; '
              'dT_G (upper) = ', self.ttd_u, 'K; '
              'dT_G (lower) = ', self.ttd_l, 'K; '
              'dT_log = ', self.dT_log, 'K; '
              'kA = ', self.kA, 'W / K; '
              'dp1 = ', self.dp1, '; '
              'dp2 = ', self.dp2, '; '
              'zeta1 = ', self.zeta1, '; '
              'zeta2 = ', self.zeta2, '; '
              'm1 = ', inlets[0].m, 'kg / s; '
              'm2 = ', inlets[1].m, 'kg / s; '
              'Sirr = ', inlets[1].m * (s_mix_ph(outlets[1].as_list()) -
                                        s_mix_ph(inlets[1].as_list())) +
              inlets[0].m * (s_mix_ph(outlets[0].as_list()) -
                             s_mix_ph(inlets[0].as_list())), 'W / K'
              )

# %%


class condenser(heat_exchanger):
    """
    component condenser inhertis from component heat exchanger
    """

    def component(self):
        return 'condenser'

    def design(self):
        return [n for n in heat_exchanger.design(self) if n != 'dp1']

    def offdesign(self):
        return [n for n in heat_exchanger.offdesign(self) if n != 'zeta1']

    def additional_equations(self, nw):

        vec_res = []
        outlets = nw.comps.loc[self].o.tolist()

        o1 = outlets[0].as_list()
        vec_res += [o1[2] - h_mix_pQ(o1, 0)]

        return vec_res

    def additional_derivatives(self, nw):

        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inlets), len(outlets)
        num_fl = len(nw.fluids)
        mat_deriv = []

        o1 = outlets[0].as_list()
        x_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
        x_deriv[0, 2, 1] = -dh_mix_dpQ(o1, 0)
        x_deriv[0, 2, 2] = 1
        mat_deriv += x_deriv.tolist()

        return mat_deriv

    def kA_func(self, inlets, outlets):

        i1 = inlets[0].as_list()
        i2 = inlets[1].as_list()
        o1 = outlets[0].as_list()
        o2 = outlets[1].as_list()

        T_i1 = T_mix_ph([i1[0], i1[1], h_mix_pQ(i1, 1), i1[3]])
        T_i2 = T_mix_ph(i2)
        T_o1 = T_mix_ph(o1)
        T_o2 = T_mix_ph(o2)

        io2 = 0
        while T_i1 <= T_o2:
            try:
                T_o2 = T_mix_ph([o2[0], o2[1], o2[2] - io2 * 10000, o2[3]])
                io2 += 1
            except:
                None

        i = 0
        while T_o1 <= T_i2:
            i += 1
            T_o1 = T_mix_ph([o1[0], o1[1], o1[2] + i * 10000, o1[3]])

        return (i1[0] * (o1[2] - i1[2]) + self.kA *
                (T_o1 - T_i2 - T_i1 + T_o2) /
                math.log((T_o1 - T_i2) / (T_i1 - T_o2)))

    def dT_log_func(self, inlets, outlets):

        T_i1 = T_mix_ph([i1[0], i1[1], h_mix_pQ(i1, 1), i1[3]])
        T_i2 = T_mix_ph(i2)
        T_o1 = T_mix_ph(o1)
        T_o2 = T_mix_ph(o2)

        io2 = 0
        while T_i1 <= T_o2:
            try:
                T_o2 = T_mix_ph([o2[0], o2[1], o2[2] - io2 * 10000, o2[3]])
                io2 += 1
            except:
                None

        i = 0
        while T_o1 <= T_i2:
            i += 1
            T_o1 = T_mix_ph([o1[0], o1[1], o1[2] + i * 10000, o1[3]])

        return (self.dT_log *
                math.log((T_o1 - T_i2) / (T_i1 - T_o2)) -
                T_o1 + T_i2 + T_i1 - T_o2)

    def ttd_u_func(self, inlets, outlets):
        i1 = inlets[0].as_list()
        o2 = outlets[1].as_list()
        return (self.ttd_u -
                T_mix_ph([i1[0], i1[1], h_mix_pQ(i1, 1), i1[3]]) +
                T_mix_ph(o2))

# %%


class desuperheater(heat_exchanger):
    """
    component desuperheater inhertis from component heat exchanger
    TODO:
    """

    def component(self):
        return 'desuperheater'

    def design(self):
        return heat_exchanger.design(self)

    def offdesign(self):
        return heat_exchanger.offdesign(self)

    def additional_equations(self, nw):

        vec_res = []
        outlets = nw.comps.loc[self].o.tolist()

        o1 = outlets[0].as_list()
        vec_res += [o1[2] - h_mix_pQ(o1, 1)]

        return vec_res

    def additional_derivatives(self, nw):

        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inlets), len(outlets)
        num_fl = len(nw.fluids)
        mat_deriv = []

        o1 = outlets[0].as_list()
        x_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
        x_deriv[0, 2, 1] = -dh_mix_dpQ(o1, 1)
        x_deriv[0, 2, 2] = 1
        mat_deriv += x_deriv.tolist()

        return mat_deriv


# %%


class drum(component):
    """
    this component is working with pure fluids only
    TODO:

    index 1
         - in: from eco
         - out: to evaporator
    index 2
         - in: from evaporator
         - out: to superheater
    """
    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1', 'out2']

    def component(self):
        return 'drum'

    def equations(self, nw):

        vec_res = []
        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())

        vec_res += self.fluid_res(inlets, outlets)
        vec_res += self.mass_flow_res(inlets, outlets)

        E_res = 0
        for i in inlets:
            E_res += i.m * i.h
        for o in outlets:
            E_res -= o.m * o.h
        vec_res += [E_res]

        p = inlets[0].p
        for c in [inlets[1]] + outlets:
            vec_res += [p - c.p]

        vec_res += [h_mix_pQ(outlets[0].as_list(), 0) - outlets[0].h]
        vec_res += [h_mix_pQ(outlets[1].as_list(), 1) - outlets[1].h]

        return vec_res

    def derivatives(self, nw):

        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inlets), len(outlets)
        num_fl = len(nw.fluids)
        mat_deriv = []

        mat_deriv += self.fluid_deriv(inlets, outlets)
        mat_deriv += self.mass_flow_deriv(inlets, outlets)

        E_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
        k = 0
        for i in inlets:
            E_deriv[0, k, 0] = i.h
            E_deriv[0, k, 2] = i.m
            k += 1
        j = 0
        for o in outlets:
            E_deriv[0, j + k, 0] = -o.h
            E_deriv[0, j + k, 2] = -o.m
            j += 1
        mat_deriv += E_deriv.tolist()

        p_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
        for k in range(num_i + num_o - 1):
            p_deriv[k, 0, 1] = 1
            p_deriv[k, k + 1, 1] = -1
        mat_deriv += p_deriv.tolist()

        o1 = outlets[0].as_list()
        o2 = outlets[1].as_list()

        x_deriv = np.zeros((num_o, num_i + num_o, num_fl + 3))
        x_deriv[0, 2, 1] = dh_mix_dpQ(o1, 0)
        x_deriv[0, 2, 2] = -1
        x_deriv[1, 3, 1] = dh_mix_dpQ(o2, 1)
        x_deriv[1, 3, 2] = -1
        mat_deriv += x_deriv.tolist()

        return np.asarray(mat_deriv)

    def initialise_source_p(self, c):
        return 10e5

    def initialise_target_p(self, c):
        return 10e5

    def initialise_source_h(self, c):
        if c.s_id == 'out1':
            return h_mix_pQ(c.as_list(), 0)
        else:
            return h_mix_pQ(c.as_list(), 1)

    def initialise_target_h(self, c):
        if c.t_id == 'in1':
            return h_mix_pQ(c.as_list(), 0)
        else:
            return h_mix_pQ(c.as_list(), 0.7)
# %%


class subsys_interface(component):
    """
    simplified heat exchanger.
    """

    def attr(self):
        return component.attr(self) + ['num_inter']

    def inlets(self):
        if self.num_inter_set:
            return ['in' + str(i + 1) for i in range(self.num_inter)]
        else:
            return ['in1']

    def outlets(self):
        if self.num_inter_set:
            return ['out' + str(i + 1) for i in range(self.num_inter)]
        else:
            return ['out1']

    def component(self):
        return 'subsystem interface'

    def equations(self, nw):

        vec_res = []
        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())

        vec_res += self.fluid_res(inlets, outlets)
        vec_res += self.mass_flow_res(inlets, outlets)
        for j in range(len(inlets)):
            i = inlets[j]
            o = outlets[j]
            vec_res += [i.p - o.p]
        for j in range(len(inlets)):
            i = inlets[j]
            o = outlets[j]
            vec_res += [i.h - o.h]

        return vec_res

    def derivatives(self, nw):

        inlets, outlets = (nw.comps.loc[self].i.tolist(),
                           nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inlets), len(outlets)
        num_fl = len(nw.fluids)
        mat_deriv = []

        mat_deriv += self.fluid_deriv(inlets, outlets)
        mat_deriv += self.mass_flow_deriv(inlets, outlets)

        p_deriv = np.zeros((num_i, num_i + num_o, num_fl + 3))
        for i in range(num_i):
            p_deriv[i, i, 1] = 1
        for j in range(num_o):
            p_deriv[j, j + i + 1, 1] = -1
        mat_deriv += p_deriv.tolist()

        h_deriv = np.zeros((num_i, num_i + num_o, num_fl + 3))
        for i in range(num_i):
            h_deriv[i, i, 2] = 1
        for j in range(num_o):
            h_deriv[j, j + i + 1, 2] = -1
        mat_deriv += h_deriv.tolist()

        return np.asarray(mat_deriv)
