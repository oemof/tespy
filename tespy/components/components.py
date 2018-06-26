"""
.. module:: components
    :synopsis:

.. moduleauthor:: Francesco Witte <francesco.witte@hs-flensburg.de>
"""

import numpy as np
import math

import CoolProp.CoolProp as CP

from tespy.helpers import (
    num_fluids, fluid_structure, MyComponentError, tespy_fluid,
    v_mix_ph, h_mix_pT, h_mix_ps, s_mix_pT, s_mix_ph, T_mix_ph, visc_mix_ph,
    dT_mix_dph, dT_mix_pdh, dT_mix_ph_dfluid, h_mix_pQ, dh_mix_dpQ,
    h_ps, s_ph,
    molar_massflow, lamb,
    molar_masses, err,
    dc_cp, dc_cc
)

from tespy.components import characteristics as cmp_char


def init_target(nw, c, start):
    """
    propagates the fluids towards connections target,
    ends when reaching sink, merge or combustion chamber

    :param nw: network to operate on
    :type nw: tespy.networks.network
    :param c: connection to initialise
    :type c: tespy.connections.connection
    :param start: fluid propagation startingpoint, in some cases needed
        to exit the recursion
    :type start: tespy.connections.connection
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
        for fluid, x in c.fluid.val.items():
            if not outconn.fluid.val_set[fluid]:
                outconn.fluid.val[fluid] = x

        init_target(nw, outconn, start)

    if isinstance(c.t, splitter):
        for outconn in nw.comps.loc[c.t].o:
            for fluid, x in c.fluid.val.items():
                if not outconn.fluid.val_set[fluid]:
                    outconn.fluid.val[fluid] = x

            init_target(nw, outconn, start)

    if isinstance(c.t, drum) and c.t != start:
        start = c.t
        for outconn in nw.comps.loc[c.t].o:
            for fluid, x in c.fluid.val.items():
                if not outconn.fluid.val_set[fluid]:
                    outconn.fluid.val[fluid] = x

            init_target(nw, outconn, start)


class component:
    r"""

    :param label: label for component
    :type label: str
    :param kwargs: for the keyword arguments see :code:`component.attr()`
    :returns: no return value
    :raises: - :code:`TypeError`, if label is not of type str
               components
             - :code:`ValueError`, if label contains forbidden characters
               (';', ',', '.')

    **example**

    .. code-block:: python

        cond = condenser('main condenser', ttd_u=5)

    creates component condenser labeled 'main condenser' and sets the
    terminal temperature difference at the upper side (hot side inlet to
    cold side outlet) to 5 K

    initialisation method is used for instances of class component and
    its children`

    allowed keywords in kwargs are 'mode' and additional keywords depending
    on the type of component you want to create
    """

    def __init__(self, label, **kwargs):

        # check if components label is of type str and for prohibited chars
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

        # check calculation mode declaration
        if self.mode not in ['man', 'auto']:
            msg = 'Mode must be \'man\' or \'auto\'.'
            raise TypeError(msg)

        # set default design and offdesign parameters
        self.design = self.default_design()
        self.offdesign = self.default_offdesign()

        # add container for components attributes
        var = self.attr_prop()

        for key in var.keys():
            self.__dict__.update({key: var[key]})

        self.set_attr(**kwargs)

#        print('Created ', self, '.')
#        print(self.__dict__)

    def set_attr(self, **kwargs):
        """
        sets, resets or unsets attributes of a connection, for the keyword
        arguments, return values and errors see object initialisation
        """
        var = self.attr()

        # set specified values
        for key in kwargs:
            if key in var:

                # data container specification
                if (isinstance(kwargs[key], dc_cp) or
                        isinstance(kwargs[key], dc_cc)):
                    self.__dict__.update({key: kwargs[key]})

                elif isinstance(self.get_attr(key), dc_cp):
                    # value specification for component properties
                    if (isinstance(kwargs[key], float) or
                            isinstance(kwargs[key], np.float64) or
                            isinstance(kwargs[key], int)):
                        if np.isnan(kwargs[key]):
                            self.get_attr(key).set_attr(is_set=False)
                        else:
                            self.get_attr(key).set_attr(val=kwargs[key])
                            self.get_attr(key).set_attr(is_set=True)

                    elif isinstance(kwargs[key], str):
                        self.get_attr(key).set_attr(val=kwargs[key])
                        self.get_attr(key).set_attr(is_set=True)

                    elif isinstance(kwargs[key], dict):
                        self.get_attr(key).set_attr(val=kwargs[key])
                        self.get_attr(key).set_attr(is_set=True)

                    elif kwargs[key] == 'var':
                        self.get_attr(key).set_attr(val=1)
                        self.get_attr(key).set_attr(is_set=True)
                        self.get_attr(key).set_attr(is_var=True)

                    # invalid datatype for keyword
                    else:
                        msg = 'Bad datatype for keyword argument ' + str(key)
                        raise TypeError(msg)

                elif isinstance(self.get_attr(key), dc_cc):
                    # value specification for component characteristics
                    if isinstance(kwargs[key], str):
                        self.get_attr(key).set_attr(method=kwargs[key])

                # invalid datatype for keyword
                else:
                    msg = 'Bad datatype for keyword argument ' + str(key)
                    raise TypeError(msg)

            elif key == 'design' or key == 'offdesign':
                if not isinstance(kwargs[key], list):
                    msg = 'Please provide the design parameters as list!'
                    raise ValueError(msg)
                if set(kwargs[key]).issubset(list(self.attr())):
                    self.__dict__.update({key: kwargs[key]})
                else:
                    msg = ('Available parameters for (off-)design'
                           'specification are: ' +
                           str(self.attr()) + '.')
                    raise ValueError(msg)

            elif key == 'mode':
                if kwargs[key] in ['man', 'auto']:
                    self.__dict__.update({key: kwargs[key]})
                else:
                    msg = 'mode must be \'man\' or \'auto\'.'
                    raise TypeError(msg)

            # invalid keyword
            else:
                msg = ('Component ' + self.label + ' has no attribute ' +
                       str(key))
                raise ValueError(msg)

#        print('Updated ', self, '.')
#        print(self.__dict__)

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
        return ['mode', 'design', 'offdesign']

    def attr_prop(self):
        return {}

    def inlets(self):
        return []

    def outlets(self):
        return []

    def default_design(self):
        return []

    def default_offdesign(self):
        return []

    def equations(self, nw):
        return []

    def derivatives(self, nw):
        return []

    def initialise_source(self, c, key):
        r"""
        returns a starting value for fluid properties at components outlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: val (*float*) - starting value for pressure at components
                  outlet in corresponding unit system, :math:`val = 0`
        """
        return 0

    def initialise_target(self, c, key):
        r"""
        returns a starting value for fluid properties at components inlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: val (*float*) - starting value for property at components
                  inlet in corresponding unit system, :math:`val = 0`
        """
        return 0

    def calc_parameters(self, nw, mode):
        return

    def print_parameters(self, nw):
        return

    def initialise_fluids(self, nw):
        return

    def convergence_check(self, nw):
        return

# %%

    def fluid_res(self, inl, outl):
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
            for fluid, x in inl[0].fluid.val.items():
                vec_res += [x - outl[0].fluid.val[fluid]]
            return vec_res

        if (isinstance(self, subsys_interface) or
                isinstance(self, heat_exchanger)):
            for i in range(len(inl)):
                for fluid, x in inl[i].fluid.val.items():
                    vec_res += [x - outl[i].fluid.val[fluid]]
            return vec_res

        if isinstance(self, splitter):
            for o in outl:
                for fluid, x in inl[0].fluid.val.items():
                    vec_res += [x - o.fluid.val[fluid]]
            return vec_res

        if isinstance(self, merge):
            res = 0
            for fluid, x in outl[0].fluid.val.items():
                res = -x * outl[0].m.val_SI
                for i in inl:
                    res += i.fluid.val[fluid] * i.m.val_SI
                vec_res += [res]
            return vec_res

        if isinstance(self, drum):
            for o in outl:
                for fluid, x in inl[0].fluid.val.items():
                    vec_res += [x - o.fluid.val[fluid]]
            return vec_res

        if isinstance(self, separator):

            for fluid, x in inl[0].fluid.val.items():
                res = x * inl[0].m.val_SI
                for o in outl:
                    res -= o.fluid.val[fluid] * o.m.val_SI
                vec_res += [res]
            return vec_res

        if isinstance(self, source) or isinstance(self, sink):
            return None

    def fluid_deriv(self, inl, outl):
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
        num_i = len(inl)
        num_o = len(outl)
        num_fl = len(inl[0].fluid.val)

        if len(self.inlets()) == 1 and len(self.outlets()) == 1:
            mat_deriv = np.zeros((num_fl, num_i + num_o, 3 + num_fl))
            i = 0
            for fluid, x in inl[0].fluid.val.items():
                mat_deriv[i, 0, i + 3] = 1
                mat_deriv[i, 1, i + 3] = -1
                i += 1
            return mat_deriv.tolist()

        if isinstance(self, heat_exchanger):
            mat_deriv = np.zeros((num_fl * 2, num_i + num_o, 3 + num_fl))
            i = 0
            for fluid in inl[0].fluid.val.keys():
                mat_deriv[i, 0, i + 3] = 1
                mat_deriv[i, 2, i + 3] = -1
                i += 1
            j = 0
            for fluid in inl[1].fluid.val.keys():
                mat_deriv[i + j, 1, j + 3] = 1
                mat_deriv[i + j, 3, j + 3] = -1
                j += 1
            return mat_deriv.tolist()

        if isinstance(self, splitter):
            mat_deriv = np.zeros((num_fl * num_i * num_o,
                                  num_i + num_o, 3 + num_fl))
            k = 0
            for o in outl:
                i = 0
                for fluid, x in inl[0].fluid.val.items():
                    mat_deriv[i + k * num_fl, 0, i + 3] = 1
                    mat_deriv[i + k * num_fl, k + 1, i + 3] = -1
                    i += 1
                k += 1
            return mat_deriv.tolist()

        if isinstance(self, merge):
            mat_deriv = np.zeros((num_fl, num_i + num_o, 3 + num_fl))
            j = 0
            for fluid, x in outl[0].fluid.val.items():
                k = 0
                for i in inl:
                    mat_deriv[j, k, 0] = i.fluid.val[fluid]
                    mat_deriv[j, k, j + 3] = i.m.val_SI
                    k += 1
                mat_deriv[j, k, 0] = -x
                mat_deriv[j, k, j + 3] = -outl[0].m.val_SI
                j += 1
            return mat_deriv.tolist()

        if isinstance(self, drum):
            mat_deriv = np.zeros((num_o * num_fl, num_i + num_o, 3 + num_fl))
            k = 0
            for o in outl:
                i = 0
                for fluid, x in inl[0].fluid.val.items():
                    mat_deriv[i + k * num_fl, 0, i + 3] = 1
                    mat_deriv[i + k * num_fl, k + 2, i + 3] = -1
                    i += 1
                k += 1
            return mat_deriv.tolist()

        if isinstance(self, separator):
            mat_deriv = np.zeros((num_fl, num_i + num_o, 3 + num_fl))
            j = 0
            for fluid, x in inl[0].fluid.val.items():
                k = 0
                for o in outl:
                    mat_deriv[j, k, 0] = -o.fluid.val[fluid]
                    mat_deriv[j, k, j + 3] = -o.m.val_SI
                    k += 1
                mat_deriv[j, 0, 0] = x
                mat_deriv[j, 0, j + 3] = inl[0].m.val_SI
                j += 1
            return mat_deriv.tolist()

        if isinstance(self, source) or isinstance(self, sink):
            return None

        if isinstance(self, subsys_interface):
            mat_deriv = np.zeros((num_fl * num_i, num_i + num_o, 3 + num_fl))
            for i in range(num_i):
                j = 0
                for fluid in inl[i].fluid.val.keys():
                    mat_deriv[i * num_fl + j, i, j + 3] = 1
                    mat_deriv[i * num_fl + j, num_i + i, j + 3] = -1
                    j += 1
            return mat_deriv.tolist()

# %%

    def mass_flow_res(self, inl, outl):
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

        heat exchanger and subsystem interface (same number of inlets and
        outlets

        .. math:: 0 = \dot{m}_{in,i} - \dot{m}_{out,i} \;
            \forall i \in inlets/outlets
        """

        if (isinstance(self, split) or
                isinstance(self, merge) or
                isinstance(self, combustion_chamber) or
                isinstance(self, combustion_chamber_stoich) or
                isinstance(self, drum) or
                (len(self.inlets()) == 1 and len(self.outlets()) == 1)):
            res = 0
            for i in inl:
                res += i.m.val_SI
            for o in outl:
                res -= o.m.val_SI
            return [res]

        if (isinstance(self, subsys_interface) or
                isinstance(self, heat_exchanger)):
            vec_res = []
            for i in range(len(inl)):
                vec_res += [inl[i].m.val_SI - outl[i].m.val_SI]
            return vec_res

        if isinstance(self, source) or isinstance(self, sink):
            return None

    def mass_flow_deriv(self, inl, outl):
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
        num_i = len(inl)
        num_o = len(outl)
        num_fl = len(inl[0].fluid.val)

        if (isinstance(self, split) or
                isinstance(self, merge) or
                isinstance(self, combustion_chamber) or
                isinstance(self, combustion_chamber_stoich) or
                isinstance(self, drum) or
                (len(self.inlets()) == 1 and len(self.outlets()) == 1)):
            mat_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            j = 0
            for i in inl:
                mat_deriv[0, j, 0] = 1
                j += 1
            k = 0
            for o in outl:
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

    def ddx_func(self, inl, outl, func, dx, pos):
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
            dm = 1e-4
        elif dx == 'p':
            dp = 1
        elif dx == 'h':
            dh = 1
        else:
            df = 1e-5

        if dx == 'fluid':
            deriv = []
            for f in inl[0].fluid.val.keys():
                val = (inl + outl)[pos].fluid.val[f]
                exp = 0
                if (inl + outl)[pos].fluid.val[f] + df <= 1:
                    (inl + outl)[pos].fluid.val[f] += df
                else:
                    (inl + outl)[pos].fluid.val[f] = 1
                exp += func(inl, outl)
                if (inl + outl)[pos].fluid.val[f] - 2 * df >= 0:
                    (inl + outl)[pos].fluid.val[f] -= 2 * df
                else:
                    (inl + outl)[pos].fluid.val[f] = 0
                exp -= func(inl, outl)
                (inl + outl)[pos].fluid.val[f] = val

                deriv += [exp / (2 * (dm + dp + dh + df))]

        else:
            exp = 0
            (inl + outl)[pos].m.val_SI += dm
            (inl + outl)[pos].p.val_SI += dp
            (inl + outl)[pos].h.val_SI += dh
            exp += func(inl, outl)

            (inl + outl)[pos].m.val_SI -= 2 * dm
            (inl + outl)[pos].p.val_SI -= 2 * dp
            (inl + outl)[pos].h.val_SI -= 2 * dh
            exp -= func(inl, outl)
            deriv = exp / (2 * (dm + dp + dh + df))

            (inl + outl)[pos].m.val_SI += dm
            (inl + outl)[pos].p.val_SI += dp
            (inl + outl)[pos].h.val_SI += dh

        return deriv

# %%

    def zeta_func(self, inl, outl):
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

        As the cross sectional area A will not change from design to offdesign
        calculation, it is possible to handle this the following way:

        .. math::
            0 = \zeta - \frac{(p_{in} - p_{out}) \cdot \pi^2}{8 \cdot
            \dot{m}_{in}^2 \cdot \frac{v_{in} + v_{out}}{2}}
        """
        i = inl[0].to_flow()
        o = outl[0].to_flow()
        if hasattr(self, 'zeta'):
            val = self.zeta.val
        else:
            val = self.zeta1.val
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
    component turbomachine can be subdivided in pump, compressor and turbine

    **available parameters**

    - P: power
    - eta_s: isentropic efficiency
    - pr: outlet to inlet pressure ratio
    - char: characteristic curve to use, characteristics are generated in
      preprocessing of offdesign calculations

    **equations**

    see tespy.components.components.turbomachine.equations

    **default design parameters**

    - pr, eta_s

    **default offdesign parameters**

    - char

    **inlets and outlets**

    - in1
    - out1
    """
    def attr(self):
        return ['P', 'eta_s', 'pr', 'eta_s_char']

    def attr_prop(self):
        return {'P': dc_cp(), 'eta_s': dc_cp(), 'pr': dc_cp(),
                'eta_s_char': dc_cc()}

    def default_design(self):
        return ['pr', 'eta_s']

    def default_offdesign(self):
        return ['eta_s_char']

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
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - vector of residual values

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_res`
        - :func:`tespy.components.components.component.mass_flow_res`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in} \right) - P\\
            0 = pr \cdot p_{in} - p_{out}

        isentropic efficiency

        - :func:`tespy.components.components.pump.eta_s_func`
        - :func:`tespy.components.components.compressor.eta_s_func`
        - :func:`tespy.components.components.turbine.eta_s_func`

        characteristics

        - :func:`tespy.components.components.pump.char_func`
        - :func:`tespy.components.components.compressor.char_func`
        - :func:`tespy.components.components.turbine.char_func`

        **additional equations**

        - :func:`tespy.components.components.turbomachine.additional_equations`
        - :func:`tespy.components.components.turbine.additional_equations`
        """

        vec_res = []
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        vec_res += self.fluid_res(inl, outl)
        vec_res += self.mass_flow_res(inl, outl)

        if self.P.is_set:
            vec_res += [inl[0].m.val_SI *
                        (outl[0].h.val_SI - inl[0].h.val_SI) -
                        self.P.val]

        if self.pr.is_set:
            vec_res += [self.pr.val * inl[0].p.val_SI - outl[0].p.val_SI]

        if self.eta_s.is_set:
            self.eta_s_res = self.eta_s_func(inl, outl)
            vec_res += [self.eta_s_res]

        vec_res += self.additional_equations(nw)

        return vec_res

    def additional_equations(self, nw):
        """
        returns vector vec_res with result of additional equations for this
        component

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - vector of residual values
        """
        return []

    def derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition

        :param nw: network using this component object
        :type nw: tespy.networks.network
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
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inl), len(outl)
        num_fl = len(nw.fluids)
        mat_deriv = []

        mat_deriv += self.fluid_deriv(inl, outl)
        mat_deriv += self.mass_flow_deriv(inl, outl)

        if self.P.is_set:
            P_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
            for k in range(num_i + num_o - 1):
                P_deriv[k, 0, 0] = outl[0].h.val_SI - inl[0].h.val_SI
                P_deriv[k, 0, 2] = -inl[0].m.val_SI
                P_deriv[k, k + 1, 2] = inl[0].m.val_SI
            mat_deriv += P_deriv.tolist()

        if self.pr.is_set:
            pr_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
            for k in range(num_i + num_o - 1):
                pr_deriv[k, 0, 1] = self.pr.val
                pr_deriv[k, k + 1, 1] = -1
            mat_deriv += pr_deriv.tolist()

        if self.eta_s.is_set:
            mat_deriv += self.eta_s_deriv(inl, outl)

        mat_deriv += self.additional_derivatives(nw)

        return np.asarray(mat_deriv)

    def additional_derivatives(self, nw):
        """
        returns matrix mat_deriv with partial derivatives for additional
        equations of this component

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*numpy array*) - matrix of partial derivatives
        """
        return []

    def eta_s_func(self, inl, outl):
        """
        see subclasses
        """
        msg = ('If you want to use eta_s as parameter, '
               'please specify which type of turbomachine you are using.')
        raise MyComponentError(msg)

    def eta_s_deriv(self, inl, outl):
        """
        see subclasses
        """
        msg = ('If you want to use eta_s as parameter, '
               'please specify which type of turbomachine you are using.')
        raise MyComponentError(msg)

    def h_os(self, inl, outl):
        """
        calculates the enthalpy at the outlet if compression or expansion is
        isentropic

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: h (*float*) - enthalpy after isentropic state change
        """
        if isinstance(inl[0], float) or isinstance(inl[0], int):
            i = inl
            o = outl
        else:
            i = inl[0].to_flow()
            o = outl[0].to_flow()

        if num_fluids(i[3]) == 1:
            for fluid, x in i[3].items():
                if x > err:
                    return h_ps(o[1], s_ph(i[1], i[2], fluid), fluid)
        else:
            T_mix = T_mix_ph(i)
            s_mix = s_mix_pT(i, T_mix)
            return h_mix_ps(o, s_mix)

    def char_func(self, inl, outl):
        raise MyComponentError('Function not available for this component.')

    def char_deriv(self, inl, outl):
        raise MyComponentError('Function not available.')

    def calc_parameters(self, nw, mode):
        """
        parameter calculation pre- or postprocessing

        **postprocessing**

        - calculate power P
        - calculate pressure ratio

        **preprocessing**

        - set references for inlet :code:`self.i0` and outlet
          :code:`self.o0` flows
        - set attribute for isentropic enthalpy difference
          :code:`self.dh_s0` at reference

        """
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        if (mode == 'pre' and 'P' in self.offdesign) or mode == 'post':
            self.P.val = inl[0].m.val_SI * (outl[0].h.val_SI - inl[0].h.val_SI)

        if (mode == 'pre' and 'pr' in self.offdesign) or mode == 'post':
            self.pr.val = outl[0].p.val_SI / inl[0].p.val_SI

        if mode == 'pre':

            self.i0 = inl[0].to_flow()
            self.o0 = outl[0].to_flow()
            self.i0[3] = self.i0[3].copy()
            self.o0[3] = self.i0[3].copy()
            self.dh_s0 = (self.h_os(self.i0, self.o0) - self.i0[2])

    def print_parameters(self, nw):

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())
        i1 = inl[0].to_flow()
        o1 = outl[0].to_flow()
        print('##### ', self.label, ' #####')
        if self.eta_s.val > 1:
            print('!!!!! Error in parametrisation of the model, '
                  'eta_s higher than 1 !!!!!')
        print('P = ', self.P.val, 'W; '
              'eta_s = ', self.eta_s.val, '; '
              'pr = ', self.pr.val, '; '
              'm = ', inl[0].m.val_SI, 'kg / s; '
              'Sirr = ', inl[0].m.val_SI *
              (s_mix_ph(o1) - s_mix_ph(i1)), 'W / K')

# %%


class pump(turbomachine):
    """
    **available parameters**

    - P: power
    - eta_s: isentropic efficiency
    - pr: outlet to inlet pressure ratio
    - char: characteristic curve to use, characteristics are generated in
      preprocessing of offdesign calculations

    **equations**

    see tespy.components.components.turbomachine.equations

    **default design parameters**

    - pr, eta_s

    **default offdesign parameters**

    - eta_s_char (method: None, parameter: v)

    .. note::

        Using the characteristic function for isentropic efficiency of the pump
        (char) is partly leading to unstable calculations, it is recommended
        to use a constant values for now.

    **inlets and outlets**

    - in1
    - out1

    .. image:: _images/pump.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
    """

    def comp_init(self, nw):

        if self.flow_char.func is None:
            method = self.flow_char.method
            x = self.flow_char.x
            y = self.flow_char.y
            self.flow_char.func = cmp_char.characteristics(method=method,
                                                           x=x, y=y)

    def component(self):
        return 'pump'

    def attr(self):
        return ['P', 'eta_s', 'pr', 'eta_s_char', 'flow_char']

    def attr_prop(self):
        return {'P': dc_cp(), 'eta_s': dc_cp(), 'pr': dc_cp(),
                'eta_s_char': dc_cc(), 'flow_char': dc_cc()}

    def additional_equations(self, nw):
        r"""
        additional equations for pumps

        - applies characteristic function for isentropic efficiency

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - residual value vector

        **optional equations**

        - :func:`tespy.components.components.turbine.char_func`
        """
        vec_res = []
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        if self.eta_s_char.is_set:
            vec_res += self.char_func(inl, outl).tolist()

        if self.flow_char.is_set:
            vec_res += self.flow_char_func(inl, outl).tolist()

        return vec_res

    def additional_derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition for the additional equations

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*list*) - matrix of partial derivatives
        """
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())
        mat_deriv = []

        if self.eta_s_char.is_set:
            mat_deriv += self.char_deriv(inl, outl)

        if self.flow_char.is_set:
            mat_deriv += self.flow_char_deriv(inl, outl)

        return mat_deriv

    def eta_s_func(self, inl, outl):
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
        return (-(outl[0].h.val_SI - inl[0].h.val_SI) * self.eta_s.val +
                (self.h_os(inl, outl) - inl[0].h.val_SI))

    def eta_s_deriv(self, inl, outl):
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

        num_i, num_o = len(inl), len(outl)
        num_fl = len(inl[0].fluid.val)
        mat_deriv = np.zeros((1, num_i + num_o, num_fl + 3))

        if abs(self.eta_s_res) > err ** (2):

            for i in range(num_i + num_o):
                mat_deriv[0, i, 1] = self.ddx_func(inl, outl,
                                                   self.eta_s_func, 'p', i)
                if i == 0:
                    mat_deriv[0, i, 2] = self.ddx_func(inl, outl,
                                                       self.eta_s_func, 'h', i)
                else:
                    mat_deriv[0, i, 2] = -self.eta_s.val

        else:
            for i in range(num_i + num_o):
                mat_deriv[0, i, 1] = -1
                mat_deriv[0, i, 2] = -1

        return mat_deriv.tolist()

    def char_func(self, inl, outl):
        r"""
        isentropic efficiency characteristic of a pump

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*numpy array*) - residual value of equation

        .. math::
            0 = -\left( h_{out} - h_{in} \right) \cdot char\left( \dot{m}_{in}
            \cdot v_{in} \right) + \left( h_{out,s} - h_{in} \right)
        """
        i = inl[0].to_flow()
        o = outl[0].to_flow()
        return np.array([((o[2] - i[2]) * self.dh_s0 /
                          (self.o0[2] - self.i0[2]) *
                          self.eta_s_char.func.f_x(i[0] * v_mix_ph(i)) -
                          (self.h_os(i, o) - i[2]))])

    def char_deriv(self, inl, outl):
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
                \frac{\partial char}{\partial p_{in}} &
                \frac{\partial char}{\partial h_{in}} & 0\\
                0 & \frac{\partial char}{\partial p_{out}} &
                \frac{\partial char}{\partial h_{out}} & 0\\
            \end{array}
            \right)

        """
        num_i, num_o = len(inl), len(outl)
        num_fl = len(inl[0].fluid.val)
        mat_deriv = np.zeros((1, num_i + num_o, num_fl + 3))

        mat_deriv[0, 0, 0] = (
            self.ddx_func(inl, outl, self.char_func, 'm', 0))
        for i in range(2):
            mat_deriv[0, i, 1] = (
                self.ddx_func(inl, outl, self.char_func, 'p', i))
            mat_deriv[0, i, 2] = (
                self.ddx_func(inl, outl, self.char_func, 'h', i))

        return mat_deriv.tolist()

    def flow_char_func(self, inl, outl):
        r"""
        equation for characteristics of a pump

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*numpy array*) - residual value of equation

        .. math::
            0 = p_{out} - p_{in} - char\left( \dot{m}_{in} \cdot v_{in} \right)
        """
        i = inl[0].to_flow()
        o = outl[0].to_flow()

        expr = i[0] * v_mix_ph(i)

        if expr > self.flow_char.func.x[-1]:
            expr = self.flow_char.func.x[-1]
        elif expr < self.flow_char.func.x[0]:
            expr = self.flow_char.func.x[0]

        return np.array([o[1] - i[1] - self.flow_char.func.f_x(expr)])

    def flow_char_deriv(self, inl, outl):
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
                \frac{\partial char}{\partial p_{in}} &
                \frac{\partial char}{\partial h_{in}} & 0\\
                0 & \frac{\partial char}{\partial p_{out}} &
                \frac{\partial char}{\partial h_{out}} & 0\\
            \end{array}
            \right)

        """
        num_i, num_o = len(inl), len(outl)
        num_fl = len(inl[0].fluid.val)
        mat_deriv = np.zeros((1, num_i + num_o, num_fl + 3))

        i = inl[0].to_flow()
        expr = i[0] * v_mix_ph(i)

        if expr > self.flow_char.func.x[-1] and not inl[0].m.val_set:
            inl[0].m.val_SI = self.flow_char.func.x[-2] / v_mix_ph(i)
        elif expr < self.flow_char.func.x[1] and not inl[0].m.val_set:
            inl[0].m.val_SI = self.flow_char.func.x[1] / v_mix_ph(i)

        mat_deriv[0, 0, 0] = (
            self.ddx_func(inl, outl, self.flow_char_func, 'm', 0))
        mat_deriv[0, 0, 2] = (
            self.ddx_func(inl, outl, self.flow_char_func, 'h', 0))
        for i in range(2):
            mat_deriv[0, i, 1] = (
                self.ddx_func(inl, outl, self.flow_char_func, 'p', i))

        return mat_deriv.tolist()

    def convergence_check(self, nw):
        """
        performs a convergence check

        - check if isentropic efficiency or characteristic is set
        - manipulate enthalpies at inlet and outlet if not specified by
          user, if function for isentropic efficiency cannot be calculated

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: no return value

         **Improvements**

         - work on this convergence check as there is no guarantee for
           successful performance
        """

        i, o = nw.comps.loc[self].i, nw.comps.loc[self].o

        if not o[0].p.val_set and o[0].p.val_SI < i[0].p.val_SI:
                o[0].p.val_SI = o[0].p.val_SI * 2
        if not i[0].p.val_set and o[0].p.val_SI < i[0].p.val_SI:
                i[0].p.val_SI = o[0].p.val_SI * 0.5

        if not o[0].h.val_set and o[0].h.val_SI < i[0].h.val_SI:
                o[0].h.val_SI = o[0].h.val_SI * 1.1
        if not i[0].h.val_set and o[0].h.val_SI < i[0].h.val_SI:
                i[0].h.val_SI = o[0].h.val_SI * 0.9

    def initialise_source(self, c, key):
        r"""
        returns a starting value for fluid properties at components outlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    outlet, :math:`val = 10^6 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    outlet,
                    :math:`val = 3 \cdot 10^5 \; \frac{\text{J}}{\text{kg}}`
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            return 3e5
        else:
            return 0

    def initialise_target(self, c, key):
        r"""
        returns a starting value for fluid properties at components inlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    inlet, :math:`val = 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    inlet,
                    :math:`val = 2,9 \cdot 10^5 \; \frac{\text{J}}{\text{kg}}`
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            return 2.9e5
        else:
            return 0

    def calc_parameters(self, nw, mode):
        """
        component specific parameter calculation pre- or postprocessing

        **postprocessing**

        - calculate isentropic efficiency

        **preprocessing**

        - generate characteristics for component
        """

        turbomachine.calc_parameters(self, nw, mode)

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        if (mode == 'pre' and 'eta_s' in self.offdesign) or mode == 'post':
            self.eta_s.val = ((self.h_os(inl, outl) - inl[0].h.val_SI) /
                              (outl[0].h.val_SI - inl[0].h.val_SI))
            if self.eta_s.val > 1 or self.eta_s.val <= 0:
                msg = ('Invalid value for isentropic efficiency.\n'
                       'eta_s =', self.eta_s.val)
                print(msg)
                nw.errors += [self]

        if (mode == 'pre' and 'eta_s_char' in self.offdesign):
            print('Creating characteristics for component ', self)
            v_opt = (self.i0[0] *
                     (v_mix_ph(self.i0) + v_mix_ph(self.o0)) / 2)
            H_opt = ((self.o0[1] - self.i0[1]) /
                     (9.81 * 2 / (v_mix_ph(self.i0) + v_mix_ph(self.o0))))
            self.eta_s_char.func = cmp_char.pump(v_opt, H_opt)

# %%


class compressor(turbomachine):
    """
    **available parameters**

    - P: power
    - eta_s: isentropic efficiency
    - pr: outlet to inlet pressure ratio
    - char: characteristic curve to use, characteristics are generated in
      preprocessing of offdesign calculations

    **additional parameter**

    - vigv: variable inlet guide vane angle

    **equations**

    see tespy.components.components.turbomachine.equations

    **default design parameters**

    - pr, eta_s

    **default offdesign parameters**

    - char_map

    **inlets and outlets**

    - in1
    - out1

    .. image:: _images/compressor.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
    """

    def component(self):
        return 'compressor'

    def attr(self):
        return ['P', 'eta_s', 'pr', 'vigv', 'char_map']

    def attr_prop(self):
        return {'P': dc_cp(), 'eta_s': dc_cp(), 'pr': dc_cp(), 'vigv': dc_cp(),
                'char_map': dc_cc(func=cmp_char.compressor())}

    def default_offdesign(self):
        return ['char_map']

    def additional_equations(self, nw):
        r"""
        additional equations for compressor

        - applies characteristic compressor map

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - residual value vector

        **optional equations**

        - :func:`tespy.components.components.compressor.char_func`
        """
        vec_res = []
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        if self.char_map.is_set:
            vec_res += self.char_func(inl, outl).tolist()

        return vec_res

    def additional_derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition for the additional equations

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*list*) - matrix of partial derivatives
        """
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())
        mat_deriv = []

        if self.char_map.is_set:
            mat_deriv += self.char_deriv(inl, outl)

        return mat_deriv

    def eta_s_func(self, inl, outl):
        r"""
        equation for isentropic efficiency of a compressor

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*float*) - residual value of equation

        .. math::
            0 = -\left( h_{out} - h_{in} \right) \cdot \eta_{s,c} +
            \left( h_{out,s} -  h_{in} \right)
        """
        return (-(outl[0].h.val_SI - inl[0].h.val_SI) * self.eta_s.val +
                (self.h_os(inl, outl) - inl[0].h.val_SI))

    def eta_s_deriv(self, inl, outl):
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

        num_i, num_o = len(inl), len(outl)
        num_fl = len(inl[0].fluid.val)
        mat_deriv = np.zeros((1, num_i + num_o, num_fl + 3))

        if abs(self.eta_s_res) > err ** (2):

            for i in range(num_i + num_o):
                mat_deriv[0, i, 1] = self.ddx_func(inl, outl,
                                                   self.eta_s_func, 'p', i)
                if i == 0:
                    mat_deriv[0, i, 2] = self.ddx_func(inl, outl,
                                                       self.eta_s_func, 'h', i)
                else:
                    mat_deriv[0, i, 2] = -self.eta_s.val

        else:
            for i in range(num_i + num_o):
                mat_deriv[0, i, 1] = -1
                mat_deriv[0, i, 2] = -1

        return mat_deriv.tolist()

    def char_func(self, inl, outl):
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
        """
        if isinstance(inl[0], float):
            i = inl
            o = outl
        else:
            i = inl[0].to_flow()
            o = outl[0].to_flow()
        n = math.sqrt(T_mix_ph(self.i0)) / math.sqrt(T_mix_ph(i))
        m = (i[0] * math.sqrt(T_mix_ph(i)) * self.i0[1] /
             (self.i0[0] * math.sqrt(T_mix_ph(self.i0)) * i[1]))

        dn = 0

        if n < min(self.char_map.func.pr.keys()):
            dn = min(self.char_map.func.pr.keys()) - n
        if n > max(self.char_map.func.pr.keys()):
            dn = max(self.char_map.func.pr.keys()) - n

        if self.vigv.is_set:

            vigv_range = self.char_map.func.get_vigv_range(n + dn, m)

            dvigv = 0
            if self.vigv.val < vigv_range[0]:
                dvigv = vigv_range[0] - self.vigv.val + 0.01
            if self.vigv.val > vigv_range[1]:
                dvigv = vigv_range[1] - self.vigv.val - 0.01

            speedline = self.char_map.func.get_speedline(n + dn,
                                                         self.vigv.val + dvigv)

            return np.array([
                    o[1] * self.i0[1] / (i[1] * self.o0[1]) - speedline[0](m),
                    ((self.h_os(i, o) - i[2]) / (o[2] - i[2])) /
                    (self.dh_s0 / (self.o0[2] - self.i0[2])) -
                    speedline[1](m)
                ])

        else:

            self.vigv.val = self.char_map.func.get_vigv(
                n + dn, m, o[1] / i[1] / (self.o0[1] / self.i0[1]))

            return np.array([
                    ((self.h_os(i, o) - i[2]) / (o[2] - i[2])) /
                    (self.dh_s0 / (self.o0[2] - self.i0[2])) -
                    self.char_map.func.get_eta(n + dn, m, self.vigv.val)
                ])

    def char_deriv(self, inl, outl):
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

        **Improvements**

        - improve asthetics, this part of code looks horrible
        """
        num_i = len(inl)
        num_o = len(outl)
        num_fl = len(inl[0].fluid.val)

        m11 = self.ddx_func(inl, outl, self.char_func, 'm', 0)
        p11 = self.ddx_func(inl, outl, self.char_func, 'p', 0)
        h11 = self.ddx_func(inl, outl, self.char_func, 'h', 0)

        p21 = self.ddx_func(inl, outl, self.char_func, 'p', 1)
        h21 = self.ddx_func(inl, outl, self.char_func, 'h', 1)

        if self.vigv.is_set:
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

    def convergence_check(self, nw):
        """
        performs a convergence check

        - check if isentropic efficiency or characteristic is set
        - manipulate enthalpies at inlet and outlet if not specified by
          user, if function for isentropic efficiency cannot be calculated

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: no return value

         **Improvements:**

         - work on this convergence check as there is no guarantee for
           successful performance
        """

        i, o = nw.comps.loc[self].i, nw.comps.loc[self].o

        if not o[0].p.val_set and o[0].p.val_SI < i[0].p.val_SI:
            o[0].p.val_SI = o[0].p.val_SI * 2
        if not i[0].p.val_set and o[0].p.val_SI < i[0].p.val_SI:
            i[0].p.val_SI = o[0].p.val_SI * 0.5

        if not o[0].h.val_set and o[0].h.val_SI < i[0].h.val_SI:
            o[0].h.val_SI = o[0].h.val_SI * 1.1
        if not i[0].h.val_set and o[0].h.val_SI < i[0].h.val_SI:
            i[0].h.val_SI = o[0].h.val_SI * 0.9

    def initialise_source(self, c, key):
        r"""
        returns a starting value for fluid properties at components outlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    outlet, :math:`val = 10^6 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    outlet,
                    :math:`val = 6 \cdot 10^5 \; \frac{\text{J}}{\text{kg}}`
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            return 6e5
        else:
            return 0

    def initialise_target(self, c, key):
        r"""
        returns a starting value for fluid properties at components inlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    inlet, :math:`val = 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    inlet,
                    :math:`val = 4 \cdot 10^5 \; \frac{\text{J}}{\text{kg}}`
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            return 4e5
        else:
            return 0

    def calc_parameters(self, nw, mode):
        """
        component specific parameter calculation pre- or postprocessing

        **postprocessing**

        - calculate isentropic efficiency

        **preprocessing**

        - generate characteristics for component
        """

        turbomachine.calc_parameters(self, nw, mode)

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        if (mode == 'pre' and 'eta_s' in self.offdesign) or mode == 'post':
            self.eta_s.val = ((self.h_os(inl, outl) -
                               inl[0].h.val_SI) /
                              (outl[0].h.val_SI - inl[0].h.val_SI))
            if self.eta_s.val > 1 or self.eta_s.val <= 0:
                msg = ('Invalid value for isentropic efficiency.\n'
                       'eta_s =', self.eta_s.val)
                print(msg)
                nw.errors += [self]

#        if (mode == 'pre' and 'char_map' in self.offdesign):
#            print('Creating characteristics for component ', self)
#            self.char_map.func = cmp_char.compressor()

    def print_parameters(self, nw):

        turbomachine.print_parameters(self, nw)

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        i1 = inl[0].to_flow()
        o1 = outl[0].to_flow()

        if self.char_map.is_set:
            n = math.sqrt(T_mix_ph(self.i0)) / math.sqrt(T_mix_ph(i1))
            m = (
                (i1[0] * math.sqrt(T_mix_ph(i1)) / i1[1]) /
                (self.i0[0] * math.sqrt(T_mix_ph(self.i0)) / self.i0[1])
                )
            vigv = self.char_map.func.get_vigv(n, m, (o1[1] * self.i0[1]) /
                                                     (i1[1] * self.o0[1]))
            if abs(self.vigv.val - vigv) > err:
                print('!!!!! Selected inlet guide vane angle is not feasible '
                      '!!!!!')
                if self.vigv.val > vigv:
                    print('calculated maximum angle:', vigv,
                          'selected:', self.vigv.val)
                else:
                    print('calculated minimum angle:', vigv,
                          'selected:', self.vigv.val)
            else:
                print('vigv = ', self.vigv.val)

# %%


class turbine(turbomachine):
    """
    **available parameters**

    - P: power
    - eta_s: isentropic efficiency
    - pr: outlet to inlet pressure ratio
    - char: characteristic curve to use, characteristics are generated in
      preprocessing of offdesign calculations

    **additional parameter**

    - cone: cone law to apply in offdesign calculation

    **equations**

    see tespy.components.components.turbomachine.equations

    **default design parameters**

    - pr, eta_s

    **default offdesign parameters**

    - eta_s_char (method: TRAUPEL, parameter: dh_s)

    **inlets and outlets**

    - in1
    - out1

    .. image:: _images/turbine.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
    """

    def comp_init(self, nw):

        if self.eta_s_char.func is None:
            method = self.eta_s_char.method
            x = self.eta_s_char.x
            y = self.eta_s_char.y
            self.eta_s_char.func = cmp_char.turbine(method=method, x=x, y=y)

        if self.cone.func is None:
            method = self.cone.method
            x = self.cone.x
            y = self.cone.y
            self.cone.func = cmp_char.characteristics(method=method, x=x, y=y)

    def component(self):
        return 'turbine'

    def attr(self):
        return ['P', 'eta_s', 'pr',  'eta_s_char', 'cone']

    def attr_prop(self):
        return {'P': dc_cp(), 'eta_s': dc_cp(), 'pr': dc_cp(),
                'eta_s_char': dc_cc(method='TRAUPEL', param='dh_s'),
                'cone': dc_cc(method='default')}

    def default_offdesign(self):
        return turbomachine.default_offdesign(self) + ['cone']

    def additional_equations(self, nw):
        r"""
        additional equations for turbines

        - applies characteristic function for isentropic efficiency
        - applies stodolas law in offdesign calculation

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - residual value vector

        **optional equations**

        - :func:`tespy.components.components.turbine.char_func`
        - :func:`tespy.components.components.turbine.cone_func`
        """
        vec_res = []
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        if self.eta_s_char.is_set:
            vec_res += self.char_func(inl, outl).tolist()

        if self.cone.is_set:
            vec_res += [self.cone_func(inl, outl)]

        return vec_res

    def additional_derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition for the additional equations

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*list*) - matrix of partial derivatives
        """
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inl), len(outl)
        num_fl = len(nw.fluids)
        mat_deriv = []

        if self.eta_s_char.is_set:
            mat_deriv += self.char_deriv(inl, outl)

        if self.cone.is_set:
            cone_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            cone_deriv[0, 0, 0] = -1
            cone_deriv[0, 0, 1] = (
                self.ddx_func(inl, outl, self.cone_func, 'p', 0))
            cone_deriv[0, 0, 2] = (
                self.ddx_func(inl, outl, self.cone_func, 'h', 0))
            cone_deriv[0, 1, 2] = (
                self.ddx_func(inl, outl, self.cone_func, 'p', 1))
            mat_deriv += cone_deriv.tolist()

        return mat_deriv

    def eta_s_func(self, inl, outl):
        r"""
        equation for isentropic efficiency of a turbine

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*float*) - residual value of equation

        .. math::
            0 = -\left( h_{out} - h_{in} \right) +
            \left( h_{out,s} -  h_{in} \right) \cdot \eta_{s,e}
        """
        return (-(outl[0].h.val_SI - inl[0].h.val_SI) +
                (self.h_os(inl, outl) - inl[0].h.val_SI) *
                self.eta_s.val)

    def eta_s_deriv(self, inl, outl):
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

        num_i, num_o = len(inl), len(outl)
        num_fl = len(inl[0].fluid.val)
        mat_deriv = np.zeros((1, num_i + num_o, num_fl + 3))

        if abs(self.eta_s_res) > err ** (2):

            for i in range(num_i + num_o):
                mat_deriv[0, i, 1] = self.ddx_func(inl, outl,
                                                   self.eta_s_func, 'p', i)
                if i == 0:
                    mat_deriv[0, i, 2] = self.ddx_func(inl, outl,
                                                       self.eta_s_func, 'h', i)
                else:
                    mat_deriv[0, i, 2] = -1

        else:
            for i in range(num_i + num_o):
                mat_deriv[0, i, 1] = -1
                mat_deriv[0, i, 2] = -1

        return mat_deriv.tolist()

    def cone_func(self, inl, outl):
        r"""
        equation for stodolas cone law

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*float*) - residual value of equation

        .. math::
            0 = \frac{\dot{m}_{in,0} \cdot p_{in}}{p_{in,0}} \cdot
            \sqrt{\frac{p_{in,0} \cdot v_{in}}{p_{in} \cdot v_{in,0}}} \cdot
            \sqrt{\frac{1 - \left(\frac{p_{out}}{p_{in}} \right)^{2}}
            {1 - \left(\frac{p_{out,0}}{p_{in,0}} \right)^{2}}} - \dot{m}_{in}
        """
        i = inl[0].to_flow()
        o = outl[0].to_flow()
        n = 1
        return (self.i0[0] * i[1] / self.i0[1] * math.sqrt(
                    self.i0[1] * v_mix_ph(self.i0) / (i[1] * v_mix_ph(i))) *
                math.sqrt(abs((1 - (o[1] / i[1]) ** ((n + 1) / n)) /
                          (1 - (self.o0[1] / self.i0[1]) ** ((n + 1) / n)))) -
                i[0])

        # Formulation with T0 / T
#        n = 1
#        return (
#            self.i1_0[0] * i1[1] / self.i1_0[1] * math.sqrt(
#                T_mix_ph(self.i1_0) / T_mix_ph(i1)) *
#            math.sqrt((1 - (o1[1] / i1[1]) ** ((n + 1) / n)) /
#                (1 - (self.o1_0[1] / self.i1_0[1]) ** ((n + 1) / n))) - i1[0])

    def char_func(self, inl, outl):
        r"""
        equation for turbine characteristics

        - calculate the isentropic efficiency as function of characteristic
          line
        - default method is TRAUPEL, see
          tespy.components.characteristics.turbine for more information on
          available methods

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*numpy array*) - residual value of equation

        .. math::
            0 = - \left( h_{out} - h_{in} \right) + \eta_{s,e,0} \cdot f\left(
            expr \right) \cdot
            \Delta h_{s}
        """
        i = inl[0].to_flow()
        o = outl[0].to_flow()

        if self.eta_s_char.param == 'dh_s':
            expr = math.sqrt(self.dh_s0 / (self.h_os(i, o) - i[2]))
        elif self.eta_s_char.param == 'm':
            expr = i[0] / self.i0[0]
        elif self.eta_s_char.param == 'v':
            expr = i[0] * v_mix_ph(i) / (self.i0[0] * v_mix_ph(self.i0))
        elif self.eta_s_char.param == 'pr':
            expr = (o[1] * self.i0[1]) / (i[1] * self.o0[1])
        else:
            msg = ('Please choose the parameter, you want to link the '
                   'isentropic efficiency to.')
            raise MyComponentError(msg)

        if expr > self.eta_s_char.func.x[-1]:
            expr = self.eta_s_char.func.x[-1]
        if expr < self.eta_s_char.func.x[0]:
            expr = self.eta_s_char.func.x[0]

        return np.array([(-(o[2] - i[2]) + (self.o0[2] - self.i0[2]) /
                          self.dh_s0 * self.eta_s_char.func.f_x(expr) *
                          (self.h_os(i, o) - i[2]))])

    def char_deriv(self, inl, outl):
        r"""
        partial derivatives for turbine characteristics

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: mat_deriv (*list*) - matrix of partial derivatives
        """
        num_i, num_o = len(inl), len(outl)
        num_fl = len(inl[0].fluid.val)

        mat_deriv = np.zeros((1, num_i + num_o, num_fl + 3))

        mat_deriv[0, 0, 0] = (
            self.ddx_func(inl, outl, self.char_func, 'm', 0))
        for i in range(num_i + num_o):
            mat_deriv[0, i, 1] = (
                self.ddx_func(inl, outl, self.char_func, 'p', i))
            mat_deriv[0, i, 2] = (
                self.ddx_func(inl, outl, self.char_func, 'h', i))

        return mat_deriv.tolist()

    def convergence_check(self, nw):
        r"""
        prevent impossible fluid properties in calculation

        - set :math:`p_{out} = \frac{p_{in}}{2}` if :math:`p_{out}>p_{in}`
        - set :math:`h_{out} = 0,9 \cdot h_{in}` if :math:`h_{out}>h_{in}`
        """
        i, o = nw.comps.loc[self].i, nw.comps.loc[self].o

        if i[0].p.val_SI < o[0].p.val_SI and not o[0].p.val_set:
            o[0].p.val_SI = i[0].p.val_SI / 2

        if i[0].h.val_SI < 10e5 and not i[0].h.val_set:
            i[0].h.val_SI = 10e5

        if o[0].h.val_SI < 8e5 and not o[0].h.val_set:
            o[0].h.val_SI = 8e5

        if i[0].h.val_SI <= o[0].h.val_SI and not o[0].h.val_set:
            o[0].h.val_SI = i[0].h.val_SI * 0.75

    def initialise_source(self, c, key):
        r"""
        returns a starting value for fluid properties at components outlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    outlet, :math:`val = 0.5 \cdot 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    outlet,
                    :math:`val = 15 \cdot 10^5 \; \frac{\text{J}}{\text{kg}}`
        """
        if key == 'p':
            return 0.5e5
        elif key == 'h':
            return 15e5
        else:
            return 0

    def initialise_target(self, c, key):
        r"""
        returns a starting value for fluid properties at components inlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    inlet, :math:`val = 25 \cdot 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    inlet,
                    :math:`val = 20 \cdot 10^5 \; \frac{\text{J}}{\text{kg}}`
        """
        if key == 'p':
            return 25e5
        elif key == 'h':
            return 20e5
        else:
            return 0

    def calc_parameters(self, nw, mode):
        """
        component specific parameter calculation pre- or postprocessing

        **postprocessing**

        - calculate isentropic efficiency

        **preprocessing**

        - generate characteristics for component
        """

        turbomachine.calc_parameters(self, nw, mode)

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        if (mode == 'pre' and 'eta_s' in self.offdesign) or mode == 'post':
            self.eta_s.val = ((outl[0].h.val_SI - inl[0].h.val_SI) /
                              (self.h_os(inl, outl) - inl[0].h.val_SI))
            if self.eta_s.val > 1 or self.eta_s.val <= 0:
                msg = ('Invalid value for isentropic efficiency.\n'
                       'eta_s =', self.eta_s)
                print(msg)
                nw.errors += [self]

# %%


class split(component):
    """
    component split

    - can be subdivided in splitter and separator

    **available parameters**

    - num_out: number of outlets (default value: 2)

    **equations**

    see tespy.components.components.split.equations

    **inlets and outlets**

    - in1
    - specify number of outlets with :code:`num_out`

    .. image:: _images/split.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
    """
    def attr(self):
        return ['num_out']

    def attr_prop(self):
        return {'num_out': dc_cp()}

    def inlets(self):
        return ['in1']

    def outlets(self):
        if self.num_out.is_set:
            return ['out' + str(i + 1) for i in range(self.num_out.val)]
        else:
            self.set_attr(num_out=2)
            return self.outlets()

    def component(self):
        return 'split'

    def equations(self, nw):
        r"""
        returns vector vec_res with result of equations for this component

        - equations are different for splitter and separator

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - vector of residual values

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_res`
        - :func:`tespy.components.components.component.mass_flow_res`

        .. math::
            0 = p_{in} - p_{out,i} \;
            \forall i \in \mathrm{outlets}

        **equations for splitter**

        .. math::
            0 = h_{in} - h_{out,i} \;
            \forall i \in \mathrm{outlets}\\

        **equations for separator**

        .. math::
            0 = T_{in} - T_{out,i} \;
            \forall i \in \mathrm{outlets}\\

        **TODO**

        - fluid separation requires power and cooling, equations have not
          been implemented!
        """
        vec_res = []
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        vec_res += self.fluid_res(inl, outl)
        vec_res += self.mass_flow_res(inl, outl)

        # pressure is the same at all connections
        for o in outl:
            vec_res += [inl[0].p.val_SI - o.p.val_SI]

        # different equations for splitter and separator
        if isinstance(self, splitter):
            for o in outl:
                vec_res += [inl[0].h.val_SI - o.h.val_SI]

        # different equations for splitter and separator
        if isinstance(self, separator):
            if num_fluids(inl[0].fluid.val) <= 1:
                for o in outl:
                    vec_res += [inl[0].h.val_SI - o.h.val_SI]
            else:
                for o in outl:
                    vec_res += [T_mix_ph(inl[0].to_flow()) -
                                T_mix_ph(o.to_flow())]

        return vec_res

    def derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*numpy array*) - matrix of partial derivatives

        """

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inl), len(outl)
        num_fl = len(nw.fluids)
        mat_deriv = []

        mat_deriv += self.fluid_deriv(inl, outl)
        mat_deriv += self.mass_flow_deriv(inl, outl)

        p_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
        k = 0
        for o in outl:
            p_deriv[k, 0, 1] = 1
            p_deriv[k, k + 1, 1] = -1
            k += 1
        mat_deriv += p_deriv.tolist()

        if isinstance(self, splitter):

            h_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
            k = 0
            for o in outl:
                h_deriv[k, 0, 2] = 1
                h_deriv[k, k + 1, 2] = -1
                k += 1

            mat_deriv += h_deriv.tolist()

        if isinstance(self, separator):
            if num_fluids(inl[0].fluid.val) <= 1:

                h_deriv = np.zeros((num_i + num_o - 1, num_i + num_o,
                                    num_fl + 3))
                k = 0
                for o in outl:
                    h_deriv[k, 0, 2] = 1
                    h_deriv[k, k + 1, 2] = -1
                    k += 1

                mat_deriv += h_deriv.tolist()
            else:

                T_deriv = np.zeros((num_i + num_o - 1, num_i + num_o,
                                    num_fl + 3))
                i = inl[0].to_flow()
                k = 0
                for o in outl:
                    o = o.to_flow()
                    T_deriv[k, 0, 1] = dT_mix_dph(i)
                    T_deriv[k, 0, 2] = dT_mix_pdh(i)
                    T_deriv[k, 0, 3:] = dT_mix_ph_dfluid(i)
                    T_deriv[k, k + 1, 1] = -dT_mix_dph(o)
                    T_deriv[k, k + 1, 2] = -dT_mix_pdh(o)
                    T_deriv[k, k + 1, 3:] = -1 * dT_mix_ph_dfluid(o)
                    k += 1

                mat_deriv += T_deriv.tolist()

        return np.asarray(mat_deriv)

    def initialise_source(self, c, key):
        r"""
        returns a starting value for fluid properties at components outlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    outlet, :math:`val = 1 \cdot 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    outlet,
                    :math:`val = 5 \cdot 10^5 \; \frac{\text{J}}{\text{kg}}`
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            return 5e5
        else:
            return 0

    def initialise_target(self, c, key):
        r"""
        returns a starting value for fluid properties at components inlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    inlet, :math:`val = 1 \cdot 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    inlet,
                    :math:`val = 5 \cdot 10^5 \; \frac{\text{J}}{\text{kg}}`
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            return 5e5
        else:
            return 0

    def print_parameters(self, nw):

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        print('##### ', self.label, ' #####')
        print('m_in = ', inl[0].m.val_SI, 'kg / s; ')
        i = 1
        for o in outl:
            print('m_out' + str(i) + ' = ', o.m.val_SI, 'kg / s; ')
            i += 1
        if isinstance(self, separator):
            print('; fluid_in:', inl[0].fluid, '; ')
            i = 1
            for o in outl:
                print('fluid_out' + str(i) + ' = ', o.fluid, 'kg / s; ')
                i += 1


class splitter(split):
    """
    **available parameters**

    - num_out: number of outlets (default value: 2)

    **equations**

    see tespy.components.components.split.equations

    **inlets and outlets**

    - in1
    - specify number of outlets with :code:`num_out`

    .. image:: _images/split.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
    """

    def component(self):
        return 'splitter'


class separator(split):
    """
    **available parameters**

    - num_out: number of outlets (default value: 2)

    **equations**

    see tespy.components.components.split.equations

    **inlets and outlets**

    - in1
    - specify number of outlets with :code:`num_out`

    .. image:: _images/split.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
    """

    def component(self):
        return 'separator'

# %%


class merge(component):
    """
    **available parameters**

    - num_in: number of inlets (default value: 2)

    **equations**

    see tespy.components.components.merge.equations

    **inlets and outlets**

    - specify number of inlets with :code:`num_in`
    - out1

    .. image:: _images/merge.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
    """
    def attr(self):
        return ['num_in']

    def attr_prop(self):
        return {'num_in': dc_cp()}

    def inlets(self):
        if self.num_in.is_set:
            return ['in' + str(i + 1) for i in range(self.num_in.val)]
        else:
            self.set_attr(num_in=2)
            return self.inlets()

    def outlets(self):
        return ['out1']

    def component(self):
        return 'merge'

    def equations(self, nw):
        r"""
        returns vector vec_res with result of equations for this component

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - vector of residual values

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_res`
        - :func:`tespy.components.components.component.mass_flow_res`

        .. math::
            0 = - \dot{m}_{out} \cdot h_{out} + \sum_{i} \dot{m}_{in,i} \cdot
            h_{in,i} \; \forall i \in \mathrm{inlets}

            0 = p_{in,i} - p_{out} \;
            \forall i \in \mathrm{inlets}
        """
        vec_res = []
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        vec_res += self.fluid_res(inl, outl)
        vec_res += self.mass_flow_res(inl, outl)

        h_res = -outl[0].m.val_SI * outl[0].h.val_SI
        for i in inl:
            h_res += i.m.val_SI * i.h.val_SI
        vec_res += [h_res]

        for i in inl:
            vec_res += [outl[0].p.val_SI - i.p.val_SI]

        return vec_res

    def derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*numpy array*) - matrix of partial derivatives
        """

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inl), len(outl)
        num_fl = len(nw.fluids)
        mat_deriv = []

        mat_deriv += self.fluid_deriv(inl, outl)
        mat_deriv += self.mass_flow_deriv(inl, outl)

        h_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
        h_deriv[0, num_i, 0] = -outl[0].h.val_SI
        h_deriv[0, num_i, 2] = -outl[0].m.val_SI
        k = 0
        for i in inl:
            h_deriv[0, k, 0] = i.h.val_SI
            h_deriv[0, k, 2] = i.m.val_SI
            k += 1
        mat_deriv += h_deriv.tolist()

        p_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
        k = 0
        for i in inl:
            p_deriv[k, k, 1] = -1
            p_deriv[k, num_i, 1] = 1
            k += 1
        mat_deriv += p_deriv.tolist()

        return np.asarray(mat_deriv)

    def initialise_fluids(self, nw):
        r"""
        fluid initialisation for fluid mixture at outlet of the merge

        - it is recommended to specify starting values for mass flows at merges

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: no return value
        """
        for outconn in nw.comps.loc[self].o:
            inl = nw.comps.loc[self].i.tolist()
            for fluid in nw.fluids:
                if not outconn.fluid.val_set[fluid]:
                    x = 0
                    m = 0
                    for i in inl:
                        m += i.m.val_SI
                        x += i.fluid.val[fluid] * i.m.val_SI

                    outconn.fluid.val[fluid] = x / m

    def initialise_source(self, c, key):
        r"""
        returns a starting value for fluid properties at components outlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    outlet, :math:`val = 1 \cdot 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    outlet,
                    :math:`val = 5 \cdot 10^5 \; \frac{\text{J}}{\text{kg}}`
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            return 5e5
        else:
            return 0

    def initialise_target(self, c, key):
        r"""
        returns a starting value for fluid properties at components inlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    inlet, :math:`val = 1 \cdot 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    inlet,
                    :math:`val = 5 \cdot 10^5 \; \frac{\text{J}}{\text{kg}}`
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            return 5e5
        else:
            return 0

    def print_parameters(self, nw):

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        print('##### ', self.label, ' #####')
        j = 1
        for i in inl:
            print('m_in' + str(j) + ' = ', i.m.val_SI, 'kg / s; ')
            j += 1
        print('m_out = ', outl[0].m.val_SI, 'kg / s; ')

# %%


class combustion_chamber(component):
    r"""
    **available parameters**

    - fuel: fuel for combustion chamber
    - lamb: air to stoichiometric air ratio
    - ti: thermal input (:math:`{LHV \cdot \dot{m}_f}`)

    **equations**

    see tespy.components.components.combustion_chamber.equations

    **available fuels**

    - methane
    - ethane
    - propane
    - butane
    - hydrogen

    **inlets and outlets**

    - in1, in2
    - out1

    .. image:: _images/combustion_chamber.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
    """

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1']

    def attr(self):
        return ['fuel', 'lamb', 'ti']

    def attr_prop(self):
        return {'fuel': dc_cp(), 'lamb': dc_cp(), 'ti': dc_cp()}

    def fuels(self):
        return ['methane', 'ethane', 'propane', 'butane',
                'hydrogen']

    def component(self):
        return 'combustion chamber'

    def comp_init(self, nw):

        if not self.fuel.is_set:
            msg = 'Must specify fuel for combustion chamber.'
            raise MyComponentError(msg)

        if (len([x for x in nw.fluids if x in [a.replace(' ', '') for a in
                 CP.get_aliases(self.fuel.val)]]) == 0):
            msg = ('The fuel you specified does not match the fuels available'
                   ' within the network.')
            raise MyComponentError(msg)

        if (len([x for x in self.fuels() if x in [a.replace(' ', '') for a in
                 CP.get_aliases(self.fuel.val)]])) == 0:
            msg = ('The fuel you specified is not available. Available fuels '
                   'are: ' + str(self.fuels()) + '.')
            raise MyComponentError(msg)

        self.fuel.val = [x for x in nw.fluids if x in [
                a.replace(' ', '') for a in CP.get_aliases(self.fuel.val)]][0]

        self.o2 = [x for x in nw.fluids if x in
                   [a.replace(' ', '') for a in CP.get_aliases('O2')]][0]
        self.co2 = [x for x in nw.fluids if x in
                    [a.replace(' ', '') for a in CP.get_aliases('CO2')]][0]
        self.h2o = [x for x in nw.fluids if x in
                    [a.replace(' ', '') for a in CP.get_aliases('H2O')]][0]
        self.n2 = [x for x in nw.fluids if x in
                   [a.replace(' ', '') for a in CP.get_aliases('N2')]][0]

        structure = fluid_structure(self.fuel.val)

        self.n = {}
        for el in ['C', 'H', 'O']:
            if el in structure.keys():
                self.n[el] = structure[el]
            else:
                self.n[el] = 0

        self.lhv = self.calc_lhv()

    def calc_lhv(self):
        r"""
        calculates the lower heating value of the combustion chambers fuel

        :returns: val (*float*) - lhv of the specified fuel

        **equation**

        .. math::
            LHV = -\frac{\sum_i {\Delta H_f^0}_i -
            \sum_j {\Delta H_f^0}_j }
            {M_{fuel}}\\
            \forall i \in \text{reation products},\\
            \forall j \in \text{reation educts},\\
            \Delta H_f^0: \text{molar formation enthalpy}

        =============== =====================================
         substance       :math:`\frac{\Delta H_f^0}{kJ/mol}`
        =============== =====================================
         hydrogen        0
         methane         -74.85
         ethane          -84.68
         propane         -103.8
         butane          -124.51
        --------------- -------------------------------------
         oxygen          0
         carbondioxide   -393.5
         water (g)       -241.8
        =============== =====================================

        """

        hf = {}
        hf['hydrogen'] = 0
        hf['methane'] = -74.85
        hf['ethane'] = -84.68
        hf['propane'] = -103.8
        hf['butane'] = -124.51
        hf[self.o2] = 0
        hf[self.co2] = -393.5
        # water (gaseous)
        hf[self.h2o] = -241.8

        key = set(list(hf.keys())).intersection(
                set([a.replace(' ', '')
                     for a in CP.get_aliases(self.fuel.val)]))

        val = (-(self.n['H'] / 2 * hf[self.h2o] + self.n['C'] * hf[self.co2] -
                 ((self.n['C'] + self.n['H'] / 4) * hf[self.o2] +
                  hf[list(key)[0]])) /
               molar_masses[self.fuel.val] * 1000)

        return val

    def equations(self, nw):
        r"""
        returns vector vec_res with result of equations for this component

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - vector of residual values

        **mandatory equations**

        - :func:`tespy.components.components.combustion_chamber.reaction_balance`
        - :func:`tespy.components.components.component.mass_flow_res`

        .. math::

            0 = p_{in,i} - p_{out} \;
            \forall i \in \mathrm{inlets}

        - :func:`tespy.components.components.combustion_chamber.energy_balance`

        **optional equations**

        - :func:`tespy.components.components.combustion_chamber.lambda_func`
        - :func:`tespy.components.components.combustion_chamber.ti_func`

        """

        vec_res = []
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        for fluid in nw.fluids:
            vec_res += [self.reaction_balance(inl, outl, fluid)]

        vec_res += self.mass_flow_res(inl, outl)

        for i in inl:
            vec_res += [outl[0].p.val_SI - i.p.val_SI]

        vec_res += [self.energy_balance(inl, outl)]

        if self.lamb.is_set:
            vec_res += [self.lambda_func(inl, outl)]

        if self.ti.is_set:
            vec_res += [self.ti_func(inl, outl)]

        return vec_res

    def derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*numpy array*) - matrix of partial derivatives
        """

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inl), len(outl)
        num_fl = len(nw.fluids)
        mat_deriv = []

        # derivatives for reaction balance
        j = 0
        fl_deriv = np.zeros((num_fl, num_i + num_o, num_fl + 3))
        for fluid in nw.fluids:
            for i in range(num_i + num_o):
                fl_deriv[j, i, 0] = self.drb_dx(inl, outl, 'm', i, fluid)
                fl_deriv[j, i, 3:] = (
                    self.drb_dx(inl, outl, 'fluid', i, fluid))

            j += 1
        mat_deriv += fl_deriv.tolist()

        # derivatives for mass balance
        mat_deriv += self.mass_flow_deriv(inl, outl)

        # derivatives for pressure equations
        p_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
        for k in range(num_i + num_o - 1):
            p_deriv[k][num_i][1] = 1
            p_deriv[k][k][1] = -1
        mat_deriv += p_deriv.tolist()

        # derivatives for energy balance
        eb_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
        for i in range(num_i + num_o):
            eb_deriv[0, i, 0] = (
                self.ddx_func(inl, outl, self.energy_balance, 'm', i))
            eb_deriv[0, i, 1] = (
                self.ddx_func(inl, outl, self.energy_balance, 'p', i))
            if i >= num_i:
                eb_deriv[0, i, 2] = -(inl + outl)[i].m.val_SI
            else:
                eb_deriv[0, i, 2] = (inl + outl)[i].m.val_SI
        mat_deriv += eb_deriv.tolist()

        if self.lamb.is_set:
            # derivatives for specified lambda
            lamb_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            for i in range(num_i):
                lamb_deriv[0, i, 0] = (
                    self.ddx_func(inl, outl, self.lambda_func, 'm', i))
                lamb_deriv[0, i, 3:] = (
                    self.ddx_func(inl, outl, self.lambda_func,
                                  'fluid', i))
            mat_deriv += lamb_deriv.tolist()

        if self.ti.is_set:
            # derivatives for specified lambda
            ti_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            for i in range(num_i):
                ti_deriv[0, i, 0] = (
                    self.ddx_func(inl, outl, self.ti_func, 'm', i))
                ti_deriv[0, i, 3:] = (
                    self.ddx_func(inl, outl, self.ti_func,
                                  'fluid', i))
            for o in range(num_o):
                ti_deriv[0, o + i, 0] = (
                    self.ddx_func(inl, outl, self.ti_func, 'm', o + i))
                ti_deriv[0, o + i, 3:] = (
                    self.ddx_func(inl, outl, self.ti_func,
                                  'fluid', o + i))
            mat_deriv += ti_deriv.tolist()

        return np.asarray(mat_deriv)

    def reaction_balance(self, inl, outl, fluid):
        r"""
        calculates the reactions mass balance for one fluid

        - determine molar mass flows of fuel and oxygen
        - calculate excess fuel
        - calculate residual value of the fluids balance

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :param fluid: fluid to calculate the reaction balance for
        :type fluid: str
        :returns: res (*float*) - residual value of mass balance

        **reaction balance equations**

        .. math::
            res = \sum_i \left(x_{fluid,i} \cdot \dot{m}_{i}\right) -
            \sum_j \left(x_{fluid,j} \cdot \dot{m}_{j}\right) \;
            \forall i \in inlets, \; \forall j \in outlets

            \dot{m}_{fluid,m} = \sum_i \frac{x_{fluid,i} \cdot \dot{m}_{i}}
            {M_{fluid}} \; \forall i \in inlets\\

            \lambda = \frac{\dot{m}_{f,m}}{\dot{m}_{O_2,m} \cdot
            \left(n_{C,fuel} + 0.25 \cdot n_{H,fuel}\right)}

        *fuel*

        .. math::
            0 = res - \left(\dot{m}_{f,m} - \dot{m}_{f,exc,m}\right)
            \cdot M_{fuel}\\

            \dot{m}_{f,exc,m} = \begin{cases}
            0 & \lambda \geq 1\\
            \dot{m}_{f,m} - \frac{\dot{m}_{O_2,m}}
            {n_{C,fuel} + 0.25 \cdot n_{H,fuel}} & \lambda < 1
            \end{cases}

        *oxygen*

        .. math::
            0 = res - \begin{cases}
            -\frac{\dot{m}_{O_2,m} \cdot M_{O_2}}{\lambda} & \lambda \geq 1\\
            - \dot{m}_{O_2,m} \cdot M_{O_2} & \lambda < 1
            \end{cases}

        *water*

        .. math::
            0 = res + \left( \dot{m}_{f,m} - \dot{m}_{f,exc,m} \right)
            \cdot 0.5 \cdot n_{H,fuel} \cdot M_{H_2O}

        *carbondioxide*

        .. math::
            0 = res + \left( \dot{m}_{f,m} - \dot{m}_{f,exc,m} \right)
            \cdot n_{C,fuel} \cdot M_{CO_2}

        *other*

        .. math::
            0 = res

        """

        n_fuel = 0
        for i in inl:
            n_fuel += (i.m.val_SI * i.fluid.val[self.fuel.val] /
                       molar_masses[self.fuel.val])

        n_oxygen = 0
        for i in inl:
            n_oxygen += (i.m.val_SI * i.fluid.val[self.o2] /
                         molar_masses[self.o2])

        if not self.lamb.is_set:
            self.lamb.val = n_oxygen / (
                    n_fuel * (self.n['C'] + self.n['H'] / 4))

        n_fuel_exc = 0
        if self.lamb.val < 1:
            n_fuel_exc = n_fuel - n_oxygen / (self.n['C'] + self.n['H'] / 4)

        if fluid == self.co2:
            dm = ((n_fuel - n_fuel_exc) *
                  self.n['C'] * molar_masses[self.co2])
        elif fluid == self.h2o:
            dm = ((n_fuel - n_fuel_exc) *
                  self.n['H'] / 2 * molar_masses[self.h2o])
        elif fluid == self.o2:
            if self.lamb.val < 1:
                dm = -n_oxygen * molar_masses[self.o2]
            else:
                dm = -n_oxygen / self.lamb.val * molar_masses[self.o2]
        elif fluid == self.fuel.val:
            dm = -(n_fuel - n_fuel_exc) * molar_masses[self.fuel.val]
        else:
            dm = 0

        res = dm

        for i in inl:
            res += i.fluid.val[fluid] * i.m.val_SI
        for o in outl:
            res -= o.fluid.val[fluid] * o.m.val_SI
        return res

    def energy_balance(self, inl, outl):
        r"""
        calculates the energy balance of the adiabatic combustion chamber

        - reference temperature: 500 K
        - reference pressure: 1 bar

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: res (*float*) - residual value of energy balance

        .. math::
            0 = \dot{m}_{in,i} \cdot \left( h_{in,i} - h_{in,i,ref} \right) -
            \dot{m}_{out,j} \cdot \left( h_{out,j} - h_{out,j,ref} \right) +
            H_{I,f} \cdot \left( \dot{m}_{in,i} \cdot x_{f,i} -
            \dot{m}_{out,j} \cdot x_{f,j} \right)

        """
        T_ref = 500
        p_ref = 1e5

        res = 0
        for i in inl:
            res += i.m.val_SI * (i.h.val_SI -
                                 h_mix_pT([i.m.val_SI, p_ref, i.h.val_SI,
                                           i.fluid.val], T_ref))
            res += i.m.val_SI * i.fluid.val[self.fuel.val] * self.lhv
        for o in outl:
            res -= o.m.val_SI * (o.h.val_SI -
                                 h_mix_pT([o.m.val_SI, p_ref, o.h.val_SI,
                                           o.fluid.val], T_ref))
            res -= o.m.val_SI * o.fluid.val[self.fuel.val] * self.lhv

        return res

    def lambda_func(self, inl, outl):
        r"""
        calculates the residual for specified lambda

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: res (*float*) - residual value of equation

        .. math::

            \dot{m}_{fluid,m} = \sum_i \frac{x_{fluid,i} \cdot \dot{m}_{i}}
            {M_{fluid}} \; \forall i \in inlets\\

            0 = \frac{\dot{m}_{f,m}}{\dot{m}_{O_2,m} \cdot
            \left(n_{C,fuel} + 0.25 \cdot n_{H,fuel}\right)} - \lambda
        """
        n_fuel = 0
        for i in inl:
            n_fuel += (i.m.val_SI * i.fluid.val[self.fuel.val] /
                       molar_masses[self.fuel.val])

        n_oxygen = 0
        for i in inl:
            n_oxygen += (i.m.val_SI * i.fluid.val[self.o2] /
                         molar_masses[self.o2])

        return (n_oxygen / (n_fuel * (self.n['C'] + self.n['H'] / 4)) -
                self.lamb.val)

    def ti_func(self, inl, outl):
        r"""
        calculates the residual for specified thermal input

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: res (*float*) - residual value of equation

        .. math::

            0 = ti - \dot{m}_f \cdot LHV
        """
        m_fuel = 0
        for i in inl:
            m_fuel += (i.m.val_SI * i.fluid.val[self.fuel.val])

        for o in outl:
            m_fuel -= (o.m.val_SI * o.fluid.val[self.fuel.val])

        return (self.ti.val - m_fuel * self.lhv)

    def drb_dx(self, inl, outl, dx, pos, fluid):
        r"""
        calculates derivative of the reaction balance to dx at components inlet
        or outlet in position pos for the fluid fluid

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :param dx: dx
        :type dx: str
        :param pos: position of inlet or outlet, logic: ['in1', 'in2', ...,
                    'out1', ...] -> 0, 1, ..., n, n + 1, ...
        :type pos: int
        :param fluid: calculate reaction balance for this fluid
        :type fluid: str
        :returns: deriv (*list* or *float*) - partial derivative of the
                  function reaction balance to dx

        .. math::

            \frac{\partial f}{\partial x} = \frac{f(x + d) + f(x - d)}
            {2 \cdot d}
        """

        dm, dp, dh, df = 0, 0, 0, 0
        if dx == 'm':
            dm = 1e-4
        elif dx == 'p':
            dp = 1
        elif dx == 'h':
            dh = 1
        else:
            df = 1e-5

        if dx == 'fluid':
            deriv = []
            for f in inl[0].fluid.val.keys():
                val = (inl + outl)[pos].fluid.val[f]
                exp = 0
                if (inl + outl)[pos].fluid.val[f] + df <= 1:
                    (inl + outl)[pos].fluid.val[f] += df
                else:
                    (inl + outl)[pos].fluid.val[f] = 1
                exp += self.reaction_balance(inl, outl, fluid)
                if (inl + outl)[pos].fluid.val[f] - 2 * df >= 0:
                    (inl + outl)[pos].fluid.val[f] -= 2 * df
                else:
                    (inl + outl)[pos].fluid.val[f] = 0
                exp -= self.reaction_balance(inl, outl, fluid)
                (inl + outl)[pos].fluid.val[f] = val

                deriv += [exp / (2 * (dm + dp + dh + df))]

        else:
            exp = 0
            (inl + outl)[pos].m.val_SI += dm
            (inl + outl)[pos].p.val_SI += dp
            (inl + outl)[pos].h.val_SI += dh
            exp += self.reaction_balance(inl, outl, fluid)

            (inl + outl)[pos].m.val_SI -= 2 * dm
            (inl + outl)[pos].p.val_SI -= 2 * dp
            (inl + outl)[pos].h.val_SI -= 2 * dh
            exp -= self.reaction_balance(inl, outl, fluid)
            deriv = exp / (2 * (dm + dp + dh + df))

            (inl + outl)[pos].m.val_SI += dm
            (inl + outl)[pos].p.val_SI += dp
            (inl + outl)[pos].h.val_SI += dh

        return deriv

    def initialise_fluids(self, nw):
        r"""
        calculates reaction balance with given lambda for good generic
        starting values

        - sets the fluid composition at the combustion chambers outlet

         for the reaction balance equations see
         :func:`tespy.components.components.combustion_chamber.reaction_balance`

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: no return value
        """
        N_2 = 0.7655
        O_2 = 0.2345

        n_fuel = 1
        lamb = 3
        m_co2 = n_fuel * self.n['C'] * molar_masses[self.co2]
        m_h2o = n_fuel * self.n['H'] / 2 * molar_masses[self.h2o]

        n_o2 = (m_co2 / molar_masses[self.co2] +
                0.5 * m_h2o / molar_masses[self.h2o]) * lamb

        m_air = n_o2 * molar_masses[self.o2] / O_2
        m_fuel = n_fuel * molar_masses[self.fuel.val]
        m_fg = m_air + m_fuel

        m_o2 = n_o2 * molar_masses[self.o2] * (1 - 1 / lamb)
        m_n2 = N_2 * m_air

        fg = {
            self.n2: m_n2 / m_fg,
            self.co2: m_co2 / m_fg,
            self.o2: m_o2 / m_fg,
            self.h2o: m_h2o / m_fg
        }

        for c in nw.comps.loc[self].o:
            for fluid, x in c.fluid.val.items():
                if not c.fluid.val_set[fluid] and fluid in fg.keys():
                    c.fluid.val[fluid] = fg[fluid]

    def convergence_check(self, nw):
        r"""
        prevent impossible fluid properties in calculation

        - check if mass fractions of fluid components at combustion chambers
          outlet are within typical range
        - propagate the corrected fluid composition towards target

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: no return value
        """
        m = 0
        for i in nw.comps.loc[self].i:
            if i.m.val_SI < 0 and not i.m.val_set:
                i.m.val_SI = 0.01
            m += i.m.val_SI

        for o in nw.comps.loc[self].o:
            fluids = [f for f in o.fluid.val.keys() if not o.fluid.val_set[f]]
            for f in fluids:
                if f not in [self.o2, self.co2, self.h2o, self.fuel]:
                    m_f = 0
                    for i in nw.comps.loc[self].i:
                        m_f += i.fluid.val[f] * i.m.val_SI

                    if abs(o.fluid.val[f] - m_f / m) > 0.03:
                        o.fluid.val[f] = m_f / m

                elif f == self.o2:
                    if o.fluid.val[f] > 0.25:
                        o.fluid.val[f] = 0.2
                    if o.fluid.val[f] < 0.05:
                        o.fluid.val[f] = 0.05

                elif f == self.co2:
                    if o.fluid.val[f] > 0.10:
                        o.fluid.val[f] = 0.10
                    if o.fluid.val[f] < 0.02:
                        o.fluid.val[f] = 0.02

                elif f == self.h2o:
                    if o.fluid.val[f] > 0.10:
                        o.fluid.val[f] = 0.10
                    if o.fluid.val[f] < 0.02:
                        o.fluid.val[f] = 0.02

                elif f == self.fuel:
                    if o.fluid.val[f] > 0:
                        o.fluid.val[f] = 0

                else:
                    continue

        for c in nw.comps.loc[self].o:
            if c.m.val_SI < 0 and not c.m.val_set:
                c.m.val_SI = 10
            init_target(nw, c, c.t)

            if c.h.val_SI < 7.5e5 and not c.h.val_set:
                c.h.val_SI = 1e6

        if self.lamb.val < 1 and not self.lamb.is_set:
            self.lamb.val = 3

    def initialise_source(self, c, key):
        r"""
        returns a starting value for fluid properties at components outlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    outlet, :math:`val = 5 \cdot 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    outlet,
                    :math:`val = 10 \cdot 10^5 \; \frac{\text{J}}{\text{kg}}`
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            return 10e5
        else:
            return 0

    def initialise_target(self, c, key):
        r"""
        returns a starting value for fluid properties at components inlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    inlet, :math:`val = 5 \cdot 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    inlet,
                    :math:`val = 5 \cdot 10^5 \; \frac{\text{J}}{\text{kg}}`
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            return 5e5
        else:
            return 0

    def calc_parameters(self, nw, mode):

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        self.ti.val = 0
        for i in inl:
            self.ti.val += i.m.val_SI * i.fluid.val[self.fuel.val] * self.lhv

        n_fuel = 0
        for i in inl:
            n_fuel += (i.m.val_SI * i.fluid.val[self.fuel.val] /
                       molar_masses[self.fuel.val])

        n_oxygen = 0
        for i in inl:
            n_oxygen += (i.m.val_SI * i.fluid.val[self.o2] /
                         molar_masses[self.o2])

        if mode == 'post':
            if not self.lamb.is_set:
                self.lamb.val = n_oxygen / (
                        n_fuel * (self.n['C'] + self.n['H'] / 4))

        if mode == 'pre':
            if 'lamb' in self.offdesign:
                self.lamb.val = n_oxygen / (
                        n_fuel * (self.n['C'] + self.n['H'] / 4))

    def print_parameters(self, nw):

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        print('##### ', self.label, ' #####')
        print('Thermal Input = ', self.ti.val,
              'lambda = ', self.lamb.val)
        j = 1
        for i in inl:
            print('m_in' + str(j) + ' = ', i.m.val_SI, 'kg / s; ')
            j += 1
        print('m_out = ', outl[0].m.val_SI, 'kg / s; ')

# %%


class combustion_chamber_stoich(combustion_chamber):
    r"""

    .. note::
        This combustion chamber uses fresh air and its fuel as the only
        reactive gas components. Therefore note the following restrictions. You
        are to

        - specify the fluid composition of the fresh air,
        - fully define the fuel's fluid components,
        - provide the aliases of the fresh air and the fuel and
        - make sure, both of the aliases are part of the network fluid vector.

        If you choose 'Air' or 'air' as alias for the fresh air, TESPy will use
        the fluid properties from CoolProp's air. Else, a custom fluid
        'TESPy::yourairalias' will be created.

        The name of the flue gas will be: 'TESPy::yourfuelalias_fg'. It is also
        possible to use fluid mixtures for the fuel, e. g.
        :code:`fuel={CH4: 0.9, 'CO2': 0.1}`. If you specify a fluid mixture for
        the fuel, TESPy will automatically create a custom fluid called:
        'TESPy::yourfuelalias'. For more information see the examples section.

    **available parameters**

    - fuel: fuel composition
    - fuel_alias: alias for fuel
    - air: air composition
    - air_alias: alias for air
    - path: path to existing fluid property table
    - lamb: air to stoichiometric air ratio
    - ti: thermal input (:math:`{LHV \cdot \dot{m}_f}`)

    **equations**

    see :func:`tespy.components.components.combustion_chamber_stoich.equations`

    **available fuel gases**

    - methane
    - ethane
    - propane
    - butane
    - hydrogen

    **inlets and outlets**

    - in1, in2
    - out1

    .. image:: _images/combustion_chamber.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
    """

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1']

    def attr(self):
        return ['fuel', 'fuel_alias', 'air', 'air_alias', 'path', 'lamb', 'ti']

    def attr_prop(self):
        return {'fuel': dc_cp(), 'fuel_alias': dc_cp(),
                'air': dc_cp(), 'air_alias': dc_cp(),
                'path': dc_cp(),
                'lamb': dc_cp(), 'ti': dc_cp()}

    def fuels(self):
        return ['methane', 'ethane', 'propane', 'butane',
                'hydrogen']

    def component(self):
        return 'combustion chamber stoichiometric flue gas'

    def comp_init(self, nw):

        if not self.fuel.is_set or not isinstance(self.fuel.val, dict):
            msg = 'Must specify fuel composition for combustion chamber.'
            raise MyComponentError(msg)

        if not self.fuel_alias.is_set:
            msg = 'Must specify fuel alias for combustion chamber.'
            raise MyComponentError(msg)
        if 'TESPy::' in self.fuel_alias.val:
            msg = 'Can not use \'TESPy::\' at this point.'
            raise MyComponentError(msg)

        if not self.air.is_set or not isinstance(self.air.val, dict):
            msg = ('Must specify air composition for combustion chamber.')
            raise MyComponentError(msg)

        if not self.air_alias.is_set:
            msg = 'Must specify air alias for combustion chamber.'
            raise MyComponentError(msg)
        if 'TESPy::' in self.air_alias.val:
            msg = 'Can not use \'TESPy::\' at this point.'
            raise MyComponentError(msg)

        for f in self.air.val.keys():
            alias = [x for x in nw.fluids if x in
                     [a.replace(' ', '') for a in CP.get_aliases(f)]]
            if len(alias) > 0:
                self.air.val[alias[0]] = self.air.val.pop(f)

        for f in self.fuel.val.keys():
            alias = [x for x in nw.fluids if x in
                     [a.replace(' ', '') for a in CP.get_aliases(f)]]
            if len(alias) > 0:
                self.fuel.val[alias[0]] = self.fuel.val.pop(f)

        for f in self.fuel.val.keys():
            alias = [x for x in self.air.val.keys() if x in
                     [a.replace(' ', '') for a in CP.get_aliases(f)]]
            if len(alias) > 0:
                self.fuel.val[alias[0]] = self.fuel.val.pop(f)

        fluids = list(self.air.val.keys()) + list(self.fuel.val.keys())

        alias = [x for x in fluids if x in
                 [a.replace(' ', '') for a in CP.get_aliases('O2')]]
        if len(alias) == 0:
            msg = 'Oxygen missing in input fluids.'
            raise MyComponentError(msg)
        else:
            self.o2 = alias[0]

        self.co2 = [x for x in nw.fluids if x in
                    [a.replace(' ', '') for a in CP.get_aliases('CO2')]]
        if len(self.co2) == 0:
            self.co2 = 'CO2'
        else:
            self.co2 = self.co2[0]

        self.h2o = [x for x in nw.fluids if x in
                    [a.replace(' ', '') for a in CP.get_aliases('H2O')]]
        if len(self.h2o) == 0:
            self.h2o = 'H2O'
        else:
            self.h2o = self.h2o[0]

        self.lhv = self.calc_lhv()
        self.stoich_flue_gas(nw)

    def calc_lhv(self):
        r"""
        calculates the lower heating value of the combustion chambers fuel

        :returns: val (*float*) - lhv of the specified fuel

        **equation**

        .. math::
            LHV = \sum_{fuels} \left(-\frac{\sum_i {\Delta H_f^0}_i -
            \sum_j {\Delta H_f^0}_j }
            {M_{fuel}} \cdot x_{fuel} \right)\\
            \forall i \in \text{reation products},\\
            \forall j \in \text{reation educts},\\
            \forall fuel \in \text{fuels},\\
            \Delta H_f^0: \text{molar formation enthalpy},\\
            x_{fuel}: \text{mass fraction of fuel in fuel mixture}

        =============== =====================================
         substance       :math:`\frac{\Delta H_f^0}{kJ/mol}`
        =============== =====================================
         hydrogen        0
         methane         -74.85
         ethane          -84.68
         propane         -103.8
         butane          -124.51
        --------------- -------------------------------------
         oxygen          0
         carbondioxide   -393.5
         water (g)       -241.8
        =============== =====================================

        """

        hf = {}
        hf['hydrogen'] = 0
        hf['methane'] = -74.85
        hf['ethane'] = -84.68
        hf['propane'] = -103.8
        hf['butane'] = -124.51
        hf['O2'] = 0
        hf['CO2'] = -393.5
        # water (gaseous)
        hf['H2O'] = -241.8

        lhv = 0

        for f, x in self.fuel.val.items():
            molar_masses[f] = CP.PropsSI('M', f)
            fl = set(list(hf.keys())).intersection(
                    set([a.replace(' ', '')
                         for a in CP.get_aliases(f)]))
            if len(fl) == 0:
                continue

            if list(fl)[0] in self.fuels():
                structure = fluid_structure(f)

                n = {}
                for el in ['C', 'H', 'O']:
                    if el in structure.keys():
                        n[el] = structure[el]
                    else:
                        n[el] = 0

                lhv += (-(n['H'] / 2 * hf['H2O'] + n['C'] * hf['CO2'] -
                          ((n['C'] + n['H'] / 4) * hf['O2'] +
                           hf[list(fl)[0]])) / molar_masses[f] * 1000) * x

        return lhv

    def stoich_flue_gas(self, nw):
        r"""
        calculates the fluid composition of the stoichiometric flue gas and
        creates a custom fluid

        - uses one mole of fuel as reference quantity and :math:`\lambda=1`
          for stoichiometric flue gas calculation (no oxygen in flue gas)
        - calculate molar quantities of (reactive) fuel components to determine
          water and carbondioxide mass fraction in flue gas
        - calculate required molar quantity for oxygen and required fresh
          air mass
        - calculate residual mass fractions for non reactive components of
          fresh air in the flue gas
        - calculate flue gas fluid composition
        - generate custom fluid porperties

        :returns: - no return value

        **reactive components in fuel**

        .. math::

            m_{fuel} = \frac{1}{M_{fuel}}\\
            m_{CO_2} = \sum_{i} \frac{x_{i} \cdot m_{fuel} \cdot num_{C,i}
            \cdot M_{CO_{2}}}{M_{i}}\\
            m_{H_{2}O} = \sum_{i} \frac{x_{i} \cdot m_{fuel} \cdot num_{H,i}
            \cdot M_{H_{2}O}}{2 \cdot M_{i}}\\
            \forall i \in \text{fuels in fuel vector},\\
            num = \text{number of atoms in molecule}

        **other components of fuel vector**

        .. math::

            m_{fg,j} = x_{j} \cdot m_{fuel}\\
            \forall j \in \text{non fuels in fuel vecotr, e. g. } CO_2,\\
            m_{fg,j} = \text{mass of fluid component j in flue gas}

        **non reactive components in air**

        .. math::

            n_{O_2} = \left( \frac{m_{CO_2}}{M_{CO_2}} +
            \frac{m_{H_{2}O}}{0,5 \cdot M_{H_{2}O}} \right) \cdot \lambda,\\
            n_{O_2} = \text{mol of oxygen required}\\
            m_{air} = \frac{n_{O_2} \cdot M_{O_2}}{x_{O_{2}, air}},\\
            m_{air} = \text{required total air mass}\\
            m_{fg,j} = x_{j, air} \cdot m_{air}\\
            m_{fg, O_2} = 0,\\
            m_{fg,j} = \text{mass of fluid component j in flue gas}

        **flue gas composition**

        .. math::

            x_{fg,j} = \frac{m_{fg, j}}{m_{air} + m_{fuel}}
        """

        lamb = 1
        n_fuel = 1
        m_fuel = 1 / molar_massflow(self.fuel.val) * n_fuel
        m_fuel_fg = m_fuel
        m_co2 = 0
        m_h2o = 0
        molar_masses[self.h2o] = CP.PropsSI('M', self.h2o)
        molar_masses[self.co2] = CP.PropsSI('M', self.co2)
        molar_masses[self.o2] = CP.PropsSI('M', self.o2)

        fg = {}

        for f, x in self.fuel.val.items():
            fl = set(list(self.fuels())).intersection(
                    set([a.replace(' ', '')
                         for a in CP.get_aliases(f)]))

            if len(fl) == 0:
                fg[f] = x * m_fuel
                continue
            else:
                n_fluid = x * m_fuel / molar_masses[f]
                m_fuel_fg -= n_fluid * molar_masses[f]
                structure = fluid_structure(f)
                n = {}
                for el in ['C', 'H', 'O']:
                    if el in structure.keys():
                        n[el] = structure[el]
                    else:
                        n[el] = 0

                m_co2 += n_fluid * n['C'] * molar_masses[self.co2]
                m_h2o += n_fluid * n['H'] / 2 * molar_masses[self.h2o]

        fg[self.co2] = m_co2
        fg[self.h2o] = m_h2o

        n_o2 = (m_co2 / molar_masses[self.co2] +
                0.5 * m_h2o / molar_masses[self.h2o]) * lamb
        m_air = n_o2 * molar_masses[self.o2] / self.air.val[self.o2]

        self.air_min = m_air / m_fuel

        for f, x in self.air.val.items():
            if f != self.o2:
                if f in fg.keys():
                    fg[f] += m_air * x
                else:
                    fg[f] = m_air * x

        m_fg = m_fuel + m_air

        for f in fg.keys():
            fg[f] /= m_fg

        tespy_fluid(self.fuel_alias.val, self.fuel.val,
                    [1000, nw.p_range_SI[1]], nw.T_range_SI, path=self.path)

        tespy_fluid(self.fuel_alias.val + '_fg', fg,
                    [1000, nw.p_range_SI[1]], nw.T_range_SI, path=self.path)

        if self.air_alias.val not in ['Air', 'air']:
            tespy_fluid(self.air_alias.val, self.air.val,
                        [1000, nw.p_range_SI[1]], nw.T_range_SI,
                        path=self.path)

    def reaction_balance(self, inl, outl, fluid):
        r"""
        calculates the reactions mass balance for one fluid

        - determine molar mass flows of fuel and fresh air
        - calculate excess fuel
        - calculate residual value of the fluids balance

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :param fluid: fluid to calculate the reaction balance for
        :type fluid: str
        :returns: res (*float*) - residual value of mass balance

        **reaction balance equations**

        .. math::
            res = \sum_i \left(x_{fluid,i} \cdot \dot{m}_{i}\right) -
            \sum_j \left(x_{fluid,j} \cdot \dot{m}_{j}\right) \;
            \forall i \in inlets, \; \forall j \in outlets

            \dot{m}_{fluid,m} = \sum_i \frac{x_{fluid,i} \cdot \dot{m}_{i}}
            {M_{fluid}} \; \forall i \in inlets\\

            \lambda = \frac{\dot{m}_{f,m}}{\dot{m}_{O_2,m} \cdot
            \left(n_{C,fuel} + 0.25 \cdot n_{H,fuel}\right)}

        *fuel*

        .. math::
            0 = res - \left(\dot{m}_{f,m} - \dot{m}_{f,exc,m}\right)
            \cdot M_{fuel}\\

            \dot{m}_{f,exc,m} = \begin{cases}
            0 & \lambda \geq 1\\
            \dot{m}_{f,m} - \frac{\dot{m}_{O_2,m}}
            {n_{C,fuel} + 0.25 \cdot n_{H,fuel}} & \lambda < 1
            \end{cases}

        *oxygen*

        .. math::
            0 = res - \begin{cases}
            -\frac{\dot{m}_{O_2,m} \cdot M_{O_2}}{\lambda} & \lambda \geq 1\\
            - \dot{m}_{O_2,m} \cdot M_{O_2} & \lambda < 1
            \end{cases}

        *water*

        .. math::
            0 = res + \left( \dot{m}_{f,m} - \dot{m}_{f,exc,m} \right)
            \cdot 0.5 \cdot n_{H,fuel} \cdot M_{H_2O}

        *carbondioxide*

        .. math::
            0 = res + \left( \dot{m}_{f,m} - \dot{m}_{f,exc,m} \right)
            \cdot n_{C,fuel} \cdot M_{CO_2}

        *other*

        .. math::
            0 = res

        """
        if self.air_alias.val in ['air', 'Air']:
            air = self.air_alias.val
        else:
            air = 'TESPy::' + self.air_alias.val
        fuel = 'TESPy::' + self.fuel_alias.val
        flue_gas = 'TESPy::' + self.fuel_alias.val + '_fg'

        m_fuel = 0
        for i in inl:
            m_fuel += i.m.val_SI * i.fluid.val[fuel]

        m_air = 0
        for i in inl:
            m_air += i.m.val_SI * i.fluid.val[air]

        if not self.lamb.is_set:
            self.lamb.val = (m_air / m_fuel) / self.air_min

        m_air_min = self.air_min * m_fuel

        m_fuel_exc = 0
        if self.lamb.val < 1:
            m_fuel_exc = m_fuel - m_air / (self.lamb.val * self.air_min)

        if fluid == air:
            if self.lamb.val >= 1:
                dm = -m_air_min
            else:
                dm = -m_air
        elif fluid == fuel:
            dm = -(m_fuel - m_fuel_exc)
        elif fluid == flue_gas:
            dm = m_air_min + m_fuel
        else:
            dm = 0

        res = dm

        for i in inl:
            res += i.fluid.val[fluid] * i.m.val_SI
        for o in outl:
            res -= o.fluid.val[fluid] * o.m.val_SI
        return res

    def energy_balance(self, inl, outl):
        r"""
        calculates the energy balance of the adiabatic combustion chamber

        - reference temperature: 500 K
        - reference pressure: 1 bar

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: res (*float*) - residual value of energy balance

        .. math::
            0 = \dot{m}_{in,i} \cdot \left( h_{in,i} - h_{in,i,ref} \right) -
            \dot{m}_{out,j} \cdot \left( h_{out,j} - h_{out,j,ref} \right) +
            H_{I,f} \cdot \left( \dot{m}_{in,i} \cdot x_{f,i} -
            \dot{m}_{out,j} \cdot x_{f,j} \right)

        """
        fuel = 'TESPy::' + self.fuel_alias.val

        T_ref = 500
        p_ref = 1e5

        res = 0
        for i in inl:
            res += i.m.val_SI * (i.h.val_SI -
                                 h_mix_pT([i.m.val_SI, p_ref, i.h.val_SI,
                                           i.fluid.val], T_ref))
            res += i.m.val_SI * i.fluid.val[fuel] * self.lhv
        for o in outl:
            res -= o.m.val_SI * (o.h.val_SI -
                                 h_mix_pT([o.m.val_SI, p_ref, o.h.val_SI,
                                           o.fluid.val], T_ref))
            res -= o.m.val_SI * o.fluid.val[fuel] * self.lhv

        return res

    def lambda_func(self, inl, outl):
        r"""
        calculates the residual for specified thermal input

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: res (*float*) - residual value of equation

        .. math::

            0 = ti - \dot{m}_f \cdot LHV
        """
        if self.air_alias.val in ['air', 'Air']:
            air = self.air_alias.val
        else:
            air = 'TESPy::' + self.air_alias.val
        fuel = 'TESPy::' + self.fuel_alias.val

        m_air = 0
        m_fuel = 0

        for i in inl:
            m_air += (i.m.val_SI * i.fluid.val[air])
            m_fuel += (i.m.val_SI * i.fluid.val[fuel])

        return (self.lamb.val - (m_air / m_fuel) / self.air_min)

    def ti_func(self, inl, outl):
        r"""
        calculates the residual for specified thermal input

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: res (*float*) - residual value of equation

        .. math::

            0 = ti - \dot{m}_f \cdot LHV
        """
        fuel = 'TESPy::' + self.fuel_alias.val

        m_fuel = 0
        for i in inl:
            m_fuel += (i.m.val_SI * i.fluid.val[fuel])

        for o in outl:
            m_fuel -= (o.m.val_SI * o.fluid.val[fuel])

        return (self.ti.val - m_fuel * self.lhv)

    def initialise_fluids(self, nw):
        r"""
        calculates reaction balance with given lambda for good generic
        starting values

        - sets the fluid composition at the combustion chambers outlet

         for the reaction balance equations see
         :func:`tespy.components.components.combustion_chamber.reaction_balance`

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: no return value
        """

        if self.air_alias.val in ['air', 'Air']:
            air = self.air_alias.val
        else:
            air = 'TESPy::' + self.air_alias.val
        flue_gas = 'TESPy::' + self.fuel_alias.val + "_fg"

        for c in nw.comps.loc[self].o:
            if not c.fluid.val_set[air]:
                c.fluid.val[air] = 0.8
            if not c.fluid.val_set[flue_gas]:
                c.fluid.val[flue_gas] = 0.2

    def convergence_check(self, nw):
        r"""
        prevent impossible fluid properties in calculation

        - check if mass fractions of fluid components at combustion chambers
          outlet are within typical range
        - propagate the corrected fluid composition towards target

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: no return value
        """

        if self.air_alias.val in ['air', 'Air']:
            air = self.air_alias.val
        else:
            air = 'TESPy::' + self.air_alias.val
        flue_gas = 'TESPy::' + self.fuel_alias.val + "_fg"
        fuel = 'TESPy::' + self.fuel_alias.val

        for c in nw.comps.loc[self].o:
            if not c.fluid.val_set[air]:
                if c.fluid.val[air] > 0.95:
                    c.fluid.val[air] = 0.95
                if c.fluid.val[air] < 0.5:
                    c.fluid.val[air] = 0.5

            if not c.fluid.val_set[flue_gas]:
                if c.fluid.val[flue_gas] > 0.5:
                    c.fluid.val[flue_gas] = 0.5
                if c.fluid.val[flue_gas] < 0.05:
                    c.fluid.val[flue_gas] = 0.05

            if not c.fluid.val_set[fuel]:
                if c.fluid.val[fuel] > 0:
                    c.fluid.val[fuel] = 0

            init_target(nw, c, c.t)

        for i in nw.comps.loc[self].i:
            if i.m.val_SI < 0 and not i.m.val_set:
                i.m.val_SI = 0.01

        for c in nw.comps.loc[self].o:
            if c.m.val_SI < 0 and not c.m.val_set:
                c.m.val_SI = 10
            init_target(nw, c, c.t)

        if self.lamb.val < 1 and not self.lamb.is_set:
            self.lamb.val = 2

    def calc_parameters(self, nw, mode):

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        if self.air_alias.val in ['air', 'Air']:
            air = self.air_alias.val
        else:
            air = 'TESPy::' + self.air_alias.val
        fuel = 'TESPy::' + self.fuel_alias.val

        m_fuel = 0
        for i in inl:
            m_fuel += i.m.val_SI * i.fluid.val[fuel]

        m_air = 0
        for i in inl:
            m_air += i.m.val_SI * i.fluid.val[air]

        if mode == 'post':
            if not self.lamb.is_set:
                self.lamb.val = (m_air / m_fuel) / self.air_min

        if mode == 'pre':
            if 'lamb' in self.offdesign:
                self.lamb.val = (m_air / m_fuel) / self.air_min

        self.ti.val = 0
        for i in inl:
            self.ti.val += i.m.val_SI * i.fluid.val[fuel] * self.lhv

# %%


class vessel(component):
    r"""
    **available parameters**

    - pr: outlet to inlet pressure ratio
    - zeta: geometry independent friction coefficient
      :math:`[\zeta]=\frac{\text{Pa}}{\text{m}^4}`, also see
      :func:`tespy.components.components.component.zeta_func`

    **equations**

    see tespy.components.components.vessel.equations

    **default design parameters**

    - pr

    **default offdesign parameters**

    - zeta

    **inlets and outlets**

    - in1
    - out1

    .. image:: _images/vessel.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
    """

    def attr(self):
        return ['pr', 'zeta']

    def attr_prop(self):
        return {'pr': dc_cp(), 'zeta': dc_cp()}

    def default_design(self):
        return ['pr']

    def default_offdesign(self):
        return ['zeta']

    def inlets(self):
        return ['in1']

    def outlets(self):
        return ['out1']

    def component(self):
        return 'vessel'

    def equations(self, nw):
        r"""
        returns vector vec_res with result of equations for this component

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - vector of residual values

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_res`
        - :func:`tespy.components.components.component.mass_flow_res`

        .. math::

            0 = h_{in} - h_{out}

        **optional equations**

        .. math::

            0 = p_{in} \cdot pr - p_{out}

        - :func:`tespy.components.components.component.zeta_func`

        """
        vec_res = []
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        vec_res += self.fluid_res(inl, outl)
        vec_res += self.mass_flow_res(inl, outl)

        vec_res += [inl[0].h.val_SI - outl[0].h.val_SI]

        if self.pr.is_set:
            vec_res += [inl[0].p.val_SI * self.pr.val - outl[0].p.val_SI]

        if self.zeta.is_set:
            vec_res += [self.zeta_func(inl, outl)]

        return vec_res

    def derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*numpy array*) - matrix of partial derivatives
        """

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inl), len(outl)
        num_fl = len(nw.fluids)
        mat_deriv = []

        mat_deriv += self.fluid_deriv(inl, outl)
        mat_deriv += self.mass_flow_deriv(inl, outl)

        h_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
        for k in range(num_i + num_o - 1):
            h_deriv[k, 0, 2] = 1
            h_deriv[k, k + 1, 2] = -1
        mat_deriv += h_deriv.tolist()

        if self.pr.is_set:
            pr_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            pr_deriv[0, 0, 1] = self.pr.val
            pr_deriv[0, 1, 1] = -1
            mat_deriv += pr_deriv.tolist()

        if self.zeta.is_set:
            zeta_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            for i in range(2):
                if i == 0:
                    zeta_deriv[0, i, 0] = (
                        self.ddx_func(inl, outl, self.zeta_func, 'm', i))
                zeta_deriv[0, i, 1] = (
                    self.ddx_func(inl, outl, self.zeta_func, 'p', i))
                zeta_deriv[0, i, 2] = (
                    self.ddx_func(inl, outl, self.zeta_func, 'h', i))
            mat_deriv += zeta_deriv.tolist()

        return np.asarray(mat_deriv)

    def initialise_source(self, c, key):
        r"""
        returns a starting value for fluid properties at components outlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    outlet, :math:`val = 4 \cdot 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    outlet,
                    :math:`val = 5 \cdot 10^5 \; \frac{\text{J}}{\text{kg}}`
        """
        if key == 'p':
            return 4e5
        elif key == 'h':
            return 5e5
        else:
            return 0

    def initialise_target(self, c, key):
        r"""
        returns a starting value for fluid properties at components inlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    inlet, :math:`val = 5 \cdot 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    inlet,
                    :math:`val = 5 \cdot 10^5 \; \frac{\text{J}}{\text{kg}}`
        """
        if key == 'p':
            return 5e5
        elif key == 'h':
            return 5e5
        else:
            return 0

    def calc_parameters(self, nw, mode):

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        if mode == 'post':
            self.pr.val = outl[0].p.val_SI / inl[0].p.val_SI
            self.zeta.val = ((inl[0].p.val_SI - outl[0].p.val_SI) *
                             math.pi ** 2 /
                             (8 * inl[0].m.val_SI ** 2 *
                             (v_mix_ph(inl[0].to_flow()) +
                              v_mix_ph(outl[0].to_flow())) / 2))

        if mode == 'pre':
            if 'pr' in self.offdesign:
                self.pr.val = outl[0].p.val_SI / inl[0].p.val_SI
            if 'zeta' in self.offdesign:
                self.zeta.val = ((inl[0].p.val_SI - outl[0].p.val_SI) *
                                 math.pi ** 2 /
                                 (8 * inl[0].m.val_SI ** 2 *
                                 (v_mix_ph(inl[0].to_flow()) +
                                  v_mix_ph(outl[0].to_flow())) / 2))

    def print_parameters(self, nw):

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        print('##### ', self.label, ' #####')
        print('pr = ', self.pr.val, '; '
              'zeta = ', self.zeta.val, 'kg / m^4 * s ; '
              'm = ', inl[0].m.val_SI, 'kg / s ; '
              'Sirr = ', inl[0].m.val_SI * (s_mix_ph(outl[0].to_flow()) -
                                            s_mix_ph(inl[0].to_flow())),
              'W / K')

# %%


class heat_exchanger_simple(component):
    r"""
    **available parameters**

    - Q: heat flux
    - pr: outlet to inlet pressure ratio
    - zeta: geometry independent friction coefficient
      :math:`[\zeta]=\frac{\text{Pa}}{\text{m}^4}`, also see
      :func:`tespy.components.components.component.zeta_func`
    - D: diameter of the pipes
    - L: length of the pipes
    - ks: pipes roughness
    - kA: area independent heat transition coefficient,
      :math:`kA=\frac{\text{W}}{\text{K}}`
    - t_a: ambient temperature, provide parameter in K

    .. note::
        for now, it is not possible to make these parameters part of the
        variable space. Thus you need to provide

        - D, L and ks, if you want to calculate pressure drop from darcy
          friction factor and
        - kA and t_a, if you want to calculate the heat flux on basis of the
          ambient conditions

    **equations**

    see tespy.components.components.heat_exchager_simple.equations

    **default design parameters**

    - pr

    **default offdesign parameters**

    - kA (method: HE_COLD, param: m): *be aware that you must provide t_a and
      t_a_design, if you want the heat flux calculated by this method*

    **inlets and outlets**

    - in1
    - out1

    .. image:: _images/pipe.svg
       :scale: 100 %
       :alt: alternative text
       :align: center

    **Improvements**

    - check design and default offdesign parameters
    """

    def comp_init(self, nw):

        if self.kA_char.func is None:
            method = self.kA_char.method
            x = self.kA_char.x
            y = self.kA_char.y
            self.kA_char.func = cmp_char.heat_ex(method=method, x=x, y=y)

        self.t_a.val_SI = ((self.t_a.val + nw.T[nw.T_unit][0]) *
                           nw.T[nw.T_unit][1])
        self.t_a_design.val_SI = ((self.t_a_design.val + nw.T[nw.T_unit][0]) *
                                  nw.T[nw.T_unit][1])

    def attr(self):
        return ['Q', 'pr', 'zeta', 'D', 'L', 'ks',
                'kA', 't_a', 't_a_design', 'kA_char']

    def attr_prop(self):
        return {'Q': dc_cp(), 'pr': dc_cp(), 'zeta': dc_cp(),
                'D': dc_cp(), 'L': dc_cp(), 'ks': dc_cp(),
                'kA': dc_cp(), 't_a': dc_cp(), 't_a_design': dc_cp(),
                'kA_char': dc_cc(method='HE_HOT', param='m')}

    def default_design(self):
        return ['pr']

    def default_offdesign(self):
        return ['kA']

    def inlets(self):
        return ['in1']

    def outlets(self):
        return ['out1']

    def component(self):
        return 'simplified heat exchanger'

    def equations(self, nw):
        r"""
        returns vector vec_res with result of equations for this component

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - vector of residual values

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_res`
        - :func:`tespy.components.components.component.mass_flow_res`

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) - \dot{Q}

        .. math::

            0 = p_{in} \cdot pr - p_{out}

        - :func:`tespy.components.components.component.zeta_func`
        - :func:`tespy.components.components.component.lamb_func`
        - :func:`tespy.components.components.component.kA_func`

        """

        vec_res = []
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        vec_res += self.fluid_res(inl, outl)
        vec_res += self.mass_flow_res(inl, outl)

        if self.Q.is_set:
            vec_res += [inl[0].m.val_SI *
                        (outl[0].h.val_SI - inl[0].h.val_SI) - self.Q.val]

        if self.pr.is_set:
            vec_res += [inl[0].p.val_SI * self.pr.val - outl[0].p.val_SI]

        if self.zeta.is_set:
            vec_res += [self.zeta_func(inl, outl)]

        if self.ks.is_set and self.D.is_set and self.L.is_set:
            vec_res += [self.lamb_func(inl, outl)]

        if self.kA.is_set and self.t_a.is_set:
            vec_res += [self.kA_func(inl, outl)]

        return vec_res

    def derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*numpy array*) - matrix of partial derivatives
        """

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inl), len(outl)
        num_fl = len(nw.fluids)
        mat_deriv = []

        mat_deriv += self.fluid_deriv(inl, outl)
        mat_deriv += self.mass_flow_deriv(inl, outl)

        if self.Q.is_set:
            Q_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            Q_deriv[0, 0, 0] = outl[0].h.val_SI - inl[0].h.val_SI
            Q_deriv[0, 0, 2] = -inl[0].m.val_SI
            Q_deriv[0, 1, 2] = inl[0].m.val_SI
            mat_deriv += Q_deriv.tolist()

        if self.pr.is_set:
            pr_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            pr_deriv[0, 0, 1] = self.pr.val
            pr_deriv[0, 1, 1] = -1
            mat_deriv += pr_deriv.tolist()

        if self.zeta.is_set:
            zeta_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            for i in range(2):
                if i == 0:
                    zeta_deriv[0, i, 0] = (
                        self.ddx_func(inl, outl, self.zeta_func, 'm', i))
                zeta_deriv[0, i, 1] = (
                    self.ddx_func(inl, outl, self.zeta_func, 'p', i))
                zeta_deriv[0, i, 2] = (
                    self.ddx_func(inl, outl, self.zeta_func, 'h', i))
            mat_deriv += zeta_deriv.tolist()

        if self.ks.is_set and self.D.is_set and self.L.is_set:
            lamb_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            for i in range(2):
                if i == 0:
                    lamb_deriv[0, i, 0] = (
                        self.ddx_func(inl, outl, self.lamb_func, 'm', i))
                lamb_deriv[0, i, 1] = (
                    self.ddx_func(inl, outl, self.lamb_func, 'p', i))
                lamb_deriv[0, i, 2] = (
                    self.ddx_func(inl, outl, self.lamb_func, 'h', i))
            mat_deriv += lamb_deriv.tolist()

        if self.kA.is_set and self.t_a.is_set:
            kA_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            kA_deriv[0, 0, 0] = self.ddx_func(inl, outl, self.kA_func, 'm', 0)
            for i in range(2):
                kA_deriv[0, i, 1] = (
                    self.ddx_func(inl, outl, self.kA_func, 'p', i))
                kA_deriv[0, i, 2] = (
                    self.ddx_func(inl, outl, self.kA_func, 'h', i))
            mat_deriv += kA_deriv.tolist()

        return np.asarray(mat_deriv)

    def lamb_func(self, inl, outl):
        r"""
        equation for pressure drop from darcy friction factor

        - calculate reynolds and darcy friction factor
        - calculate pressure drop

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*float*) - residual value of equation

        .. math::

            Re = \frac{4 \cdot \dot{m}_{in}}{\pi \cdot D \cdot
            \frac{\eta_{in}+\eta_{out}}{2}}\\

            0 = p_{in} - p_{out} - \frac{8 \cdot \dot{m}_{in}^2 \cdot
            \frac{v_{in}+v_{out}}{2} \cdot L \cdot \lambda\left(
            Re, ks, D\right)}{\pi^2 \cdot D^5}\\

            \eta: \text{dynamic viscosity}\\
            v: \text{specific volume}\\
            \lambda: \text{darcy friction factor}
        """
        i, o = inl[0].to_flow(), outl[0].to_flow()
        visc_i, visc_o = visc_mix_ph(i), visc_mix_ph(o)
        v_i, v_o = v_mix_ph(i), v_mix_ph(o)
        re = 4 * inl[0].m.val_SI / (math.pi * self.D.val *
                                    (visc_i + visc_o) / 2)

        return ((inl[0].p.val_SI - outl[0].p.val_SI) -
                8 * inl[0].m.val_SI ** 2 * (v_i + v_o) / 2 * self.L.val *
                lamb(re, self.ks.val, self.D.val) /
                (math.pi ** 2 * self.D.val ** 5))

    def kA_func(self, inl, outl):
        r"""
        equation for heat flux from ambient conditions

        - determine hot side and cold side of the heat exchanger
        - calculate heat flux

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*float*) - residual value of equation

        .. math::

            ttd_u = \begin{cases}
            t_a - T_{out} & t_a > T_{in}\\
            T_{in} - t_a & t_a \leq T_{in}
            \end{cases}

            ttd_l = \begin{cases}
            t_a - T_{in} & t_a > T_{in}\\
            T_{out} - t_a & t_a \leq T_{in}
            \end{cases}

            0 = \dot{m}_{in} \cdot \left( h_{out} - h_{in}\right) +
            kA \cdot \frac{ttd_u - ttd_l}{\ln{\frac{ttd_u}{ttd_l}}}

            t_a: \text{ambient temperature}
        """

        i, o = inl[0], outl[0]
        T_i = T_mix_ph(i.to_flow())
        T_o = T_mix_ph(o.to_flow())

        if self.t_a.val_SI > T_i:
            ttd_u = self.t_a.val_SI - T_o
            ttd_l = self.t_a.val_SI - T_i
        else:
            ttd_u = T_i - self.t_a.val_SI
            ttd_l = T_o - self.t_a.val_SI

        if self.kA_char.param == 'm':
            expr = i.m.val_SI / self.i0[0]
        else:
            expr = 1

        if expr > self.kA_char.func.x[-1]:
            expr = self.kA_char.func.x[-1]
        if expr < self.kA_char.func.x[0]:
            expr = self.kA_char.func.x[0]

        fkA = self.kA_char.func.f_x(expr)

        return (i.m.val_SI * (o.h.val_SI - i.h.val_SI) + self.kA.val * fkA * (
                (ttd_u - ttd_l) / math.log(ttd_u / ttd_l)))

    def initialise_source(self, c, key):
        r"""
        returns a starting value for fluid properties at components outlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    outlet, :math:`val = 1 \cdot 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    outlet,

        .. math::
            h = \begin{cases}
            1 \cdot 10^5 \; \frac{\text{J}}{\text{kg}} & Q < 0\\
            3 \cdot 10^5 \; \frac{\text{J}}{\text{kg}} & Q = 0\\
            5 \cdot 10^5 \; \frac{\text{J}}{\text{kg}} & Q > 0
            \end{cases}`
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            if self.Q.val < 0:
                return 1e5
            elif self.Q.val > 0:
                return 5e5
            else:
                return 3e5
        else:
            return 0

    def initialise_target(self, c, key):
        r"""
        returns a starting value for fluid properties at components inlet

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    inlet, :math:`val = 1 \cdot 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    inlet,
                    :math:`val = 5 \cdot 10^5 \; \frac{\text{J}}{\text{kg}}`

        .. math::
            h = \begin{cases}
            5 \cdot 10^5 \; \frac{\text{J}}{\text{kg}} & Q < 0\\
            3 \cdot 10^5 \; \frac{\text{J}}{\text{kg}} & Q = 0\\
            1 \cdot 10^5 \; \frac{\text{J}}{\text{kg}} & Q > 0
            \end{cases}
        """
        if key == 'p':
            return 1e5
        elif key == 'h':
            if self.Q.val < 0:
                return 5e5
            elif self.Q.val > 0:
                return 1e5
            else:
                return 3e5
        else:
            return 0

    def calc_parameters(self, nw, mode):

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        if mode == 'pre':

            self.i0 = inl[0].to_flow()
            self.o0 = outl[0].to_flow()
            self.i0[3] = self.i0[3].copy()
            self.o0[3] = self.o0[3].copy()

        if nw.mode == 'design':
            if self.t_a_design.is_set:
                t_a = self.t_a_design.val
            else:
                t_a = np.nan

        if nw.mode == 'offdesign':
            if mode == 'pre':
                if self.t_a_design.is_set:
                    t_a = self.t_a_design.val_SI
                else:
                    t_a = np.nan
            else:
                if self.t_a.is_set:
                    t_a = self.t_a.val_SI
                else:
                    t_a = np.nan

        if t_a != np.nan:

            T_i = T_mix_ph(inl[0].to_flow())
            T_o = T_mix_ph(outl[0].to_flow())

            if t_a > T_i:
                ttd_u = t_a - T_o
                ttd_l = t_a - T_i
            else:
                ttd_u = T_i - t_a
                ttd_l = T_o - t_a

            if ttd_u < 0 or ttd_l < 0:
                msg = ('Invalid value for terminal temperature '
                       'difference.'
                       'ttd_u =', ttd_u,
                       'ttd_l =', ttd_l)
                print(msg)
                nw.errors += [self]

            self.kA.val = inl[0].m.val_SI * (
                            outl[0].h.val_SI - inl[0].h.val_SI) / (
                            (ttd_u - ttd_l) / math.log(ttd_l / ttd_u))

        if (mode == 'pre' and 'Q' in self.offdesign) or mode == 'post':
            self.Q.val = inl[0].m.val_SI * (outl[0].h.val_SI -
                                            inl[0].h.val_SI)
        if (mode == 'pre' and 'pr' in self.offdesign) or mode == 'post':
            self.pr.val = outl[0].p.val_SI / inl[0].p.val_SI
        if (mode == 'pre' and 'zeta' in self.offdesign) or mode == 'post':
            self.zeta.val = ((inl[0].p.val_SI - outl[0].p.val_SI) *
                             math.pi ** 2 /
                             (8 * inl[0].m.val_SI ** 2 *
                             (v_mix_ph(inl[0].to_flow()) +
                              v_mix_ph(outl[0].to_flow())) / 2))

        # improve this part (for heat exchangers only atm)
        if self.kA.is_set:
            expr = inl[0].m.val_SI / self.i0[0]
            if (expr > self.kA_char.func.x[-1] or
                    expr < self.kA_char.func.x[0]):
                msg = ('Warning: Expression for characteristics out of bounds:'
                       ' value is ' + str(expr))
                print(msg)
                nw.errors += [self]

    def print_parameters(self, nw):

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        print('##### ', self.label, ' #####')
        print('Q = ', self.Q.val, 'W; '
              'pr = ', self.pr.val, '; '
              'zeta = ', self.zeta.val, 'kg / m^4 * s; '
              'm = ', inl[0].m.val_SI, 'kg / s; '
              'Sq = ', inl[0].m.val_SI * (s_mix_ph(outl[0].to_flow()) -
                                          s_mix_ph(inl[0].to_flow())),
              'W / K; ')
        if self.t_a.is_set or self.t_a_design.is_set:
            print('kA = ', self.kA.val, 'W / (m^2 * K)')

# %%


class pipe(heat_exchanger_simple):
    r"""

    class pipe is an alias of class heat_exchanger_simple

    **available parameters**

    - Q: heat flux
    - pr: outlet to inlet pressure ratio
    - zeta: geometry independent friction coefficient
      :math:`[\zeta]=\frac{\text{Pa}}{\text{m}^4}`, also see
      :func:`tespy.components.components.component.zeta_func`
    - D: diameter of the pipes
    - L: length of the pipes
    - ks: pipes roughness
    - kA: area independent heat transition coefficient,
      :math:`kA=\frac{\text{W}}{\text{K}}`
    - t_a: ambient temperature, provide parameter in network's temperature unit
    - t_a_design: ambient temperature design case, provide parameter in
      network's temperature unit

    **equations**

    see tespy.components.components.heat_exchager_simple.equations

    **default design parameters**

    - pr

    **default offdesign parameters**

    - kA (method: HE_COLD, param: m): *be aware that you must provide t_a and
      t_a_design, if you want the heat flux calculated by this method*

    **inlets and outlets**

    - in1
    - out1

    .. image:: _images/pipe.svg
       :scale: 100 %
       :alt: alternative text
       :align: center

    **Improvements**

    - check design and default offdesign parameters
    """
    def component(self):
        return 'pipe'

# %%


class heat_exchanger(component):
    r"""
    - all components of class heat exchanger are counter flow heat exchangers

    **available parameters**

    - Q: heat flux
    - kA: area independent heat transition coefficient,
      :math:`kA=\frac{\text{W}}{\text{K}}`
    - ttd_u: upper terminal temperature difference
    - ttd_l: lower terminal temperature difference
    - pr1: outlet to inlet pressure ratio at hot side
    - pr2: outlet to inlet pressure ratio at cold side
    - zeta1: geometry independent friction coefficient hot side
      :math:`[\zeta]=\frac{\text{Pa}}{\text{m}^4}`, also see
      :func:`tespy.components.components.component.zeta_func`
    - zeta2: geometry independent friction coefficient cold side
      :math:`[\zeta]=\frac{\text{Pa}}{\text{m}^4}`, also see
      :func:`tespy.components.components.heat_exchanger.zeta2_func`

    **equations**

    see tespy.components.components.heat_exchager.equations

    **default design parameters**

    - pr1, pr2, ttd_u, ttd_l

    **default offdesign parameters**

    - zeta1, zeta2, kA (using kA_char1/kA_char2, method: HE_HOT/HE_COLD,
      param: m/m)

    **inlets and outlets**

    - in1, in2 (index 1: hot side, index 2: cold side)
    - out1, out2 (index 1: hot side, index 2: cold side)

    .. image:: _images/heat_exchanger.svg
       :scale: 100 %
       :alt: alternative text
       :align: center

    **Improvements**

    - add the partial derivatives for specified logarithmic temperature
      difference
    - add direct current heat exchangers
    """

    def comp_init(self, nw):

        if self.kA_char1.func is None:
            method = self.kA_char1.method
            x = self.kA_char1.x
            y = self.kA_char1.y
            self.kA_char1.func = cmp_char.heat_ex(method=method, x=x, y=y)

        if self.kA_char2.func is None:
            method = self.kA_char2.method
            x = self.kA_char2.x
            y = self.kA_char2.y
            self.kA_char2.func = cmp_char.heat_ex(method=method, x=x, y=y)

    def attr(self):
        return ['Q', 'kA', 'td_log', 'kA_char1', 'kA_char2',
                'ttd_u', 'ttd_l',
                'pr1', 'pr2', 'zeta1', 'zeta2']

    def attr_prop(self):
        return {'Q': dc_cp(), 'kA': dc_cp(), 'td_log': dc_cp(),
                'kA_char1': dc_cc(method='HE_HOT', param='m'),
                'kA_char2': dc_cc(method='HE_COLD', param='m'),
                'ttd_u': dc_cp(), 'ttd_l': dc_cp(),
                'pr1': dc_cp(), 'pr2': dc_cp(),
                'zeta1': dc_cp(), 'zeta2': dc_cp()}
        # derivatives for logarithmic temperature difference not implemented
#        return (component.attr(self) +
#                ['Q', 'kA', 'td_log', 'ttd_u', 'ttd_l',
#                 'pr1', 'pr2', 'zeta1', 'zeta2'])

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1', 'out2']

    def default_design(self):
        return ['ttd_u', 'ttd_l', 'pr1', 'pr2']

    def default_offdesign(self):
        return ['kA', 'zeta1', 'zeta2']

    def component(self):
        return 'heat_exchanger'

    def equations(self, nw):
        r"""
        returns vector vec_res with result of equations for this component

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - vector of residual values

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_res`
        - :func:`tespy.components.components.component.mass_flow_res`

        .. math::

            0 = \dot{m}_{1,in} \cdot \left(h_{1,out} - h_{1,in} \right) +
            \dot{m}_{2,in} \cdot \left(h_{2,out} - h_{2,in} \right)

        **optional equations**

        .. math::

            0 = \dot{m}_{in} \cdot \left(h_{out} - h_{in} \right) - \dot{Q}

        - :func:`tespy.components.components.component.kA_func`
        - :func:`tespy.components.components.component.ttd_u_func`
        - :func:`tespy.components.components.component.ttd_l_func`

        .. math::

            0 = p_{1,in} \cdot pr1 - p_{1,out}\\
            0 = p_{2,in} \cdot pr2 - p_{2,out}

        - :func:`tespy.components.components.component.zeta_func`
        - :func:`tespy.components.components.component.zeta2_func`

        **additional equations**

        - :func:`tespy.components.components.heat_exchanger.additional_equations`
        - :func:`tespy.components.components.condenser.additional_equations`
        - :func:`tespy.components.components.desuperheater.additional_equations`

        """
        vec_res = []
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        vec_res += self.fluid_res(inl, outl)
        vec_res += self.mass_flow_res(inl, outl)

        vec_res += [inl[0].m.val_SI * (outl[0].h.val_SI -
                                       inl[0].h.val_SI) +
                    inl[1].m.val_SI * (outl[1].h.val_SI -
                                       inl[1].h.val_SI)]

        if self.Q.is_set:
            vec_res += [inl[0].m.val_SI * (outl[0].h.val_SI -
                                           inl[0].h.val_SI) - self.Q.val]

        if self.kA.is_set:
            vec_res += [self.kA_func(inl, outl)]

        # derivatives for logarithmic temperature difference not implemented
#        if self.td_log_set:
#            vec_res += [self.td_log_func(inl, outl)]

        if self.ttd_u.is_set:
            vec_res += [self.ttd_u_func(inl, outl)]

        if self.ttd_l.is_set:
            vec_res += [self.ttd_l_func(inl, outl)]

        if self.pr1.is_set:
            vec_res += [self.pr1.val * inl[0].p.val_SI - outl[0].p.val_SI]

        if self.pr2.is_set:
            vec_res += [self.pr2.val * inl[1].p.val_SI - outl[1].p.val_SI]

        if self.zeta1.is_set:
            vec_res += [self.zeta_func(inl, outl)]

        if self.zeta2.is_set:
            vec_res += [self.zeta2_func(inl, outl)]

        vec_res += self.additional_equations(nw)

        return vec_res

    def additional_equations(self, nw):
        """
        returns vector vec_res with result of additional equations for this
        component

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - vector of residual values
        """
        return []

    def derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*numpy array*) - matrix of partial derivatives
        """

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inl), len(outl)
        num_fl = len(nw.fluids)
        mat_deriv = []

        mat_deriv += self.fluid_deriv(inl, outl)
        mat_deriv += self.mass_flow_deriv(inl, outl)

        q_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
        for k in range(num_i):
            q_deriv[0, k, 0] = outl[k].h.val_SI - inl[k].h.val_SI

            q_deriv[0, k, 2] = -inl[k].m.val_SI
        q_deriv[0, 2, 2] = inl[0].m.val_SI
        q_deriv[0, 3, 2] = inl[1].m.val_SI
        mat_deriv += q_deriv.tolist()

        if self.Q.is_set:
            Q_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            Q_deriv[0, 0, 0] = outl[0].h.val_SI - inl[0].h.val_SI
            Q_deriv[0, 0, 2] = -inl[0].m.val_SI
            Q_deriv[0, 1, 2] = inl[0].m.val_SI
            mat_deriv += Q_deriv.tolist()

        if self.kA.is_set:
            kA_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            kA_deriv[0, 0, 0] = self.ddx_func(inl, outl, self.kA_func, 'm', 0)
            kA_deriv[0, 1, 0] = self.ddx_func(inl, outl, self.kA_func, 'm', 1)
            for i in range(num_i + num_o):
                kA_deriv[0, i, 1] = (
                    self.ddx_func(inl, outl, self.kA_func, 'p', i))
                kA_deriv[0, i, 2] = (
                    self.ddx_func(inl, outl, self.kA_func, 'h', i))
            mat_deriv += kA_deriv.tolist()

        # derivatives for logarithmic temperature difference not implemented
#        if self.td_log_set:
#            mat_deriv += [[[0,
#                       self.ddx_func(i1, i2, o1, o2, self.td_log_func, 'p11'),
#                       self.ddx_func(i1, i2, o1, o2, self.td_log_func, 'h11')] + z,
#                      [0,
#                       self.ddx_func(i1, i2, o1, o2, self.td_log_func, 'p12'),
#                       self.ddx_func(i1, i2, o1, o2, self.td_log_func, 'h12')] + z,
#                      [0,
#                       self.ddx_func(i1, i2, o1, o2, self.td_log_func, 'p21'),
#                       self.ddx_func(i1, i2, o1, o2, self.td_log_func, 'h21') + i1[0]] + z,
#                      [0,
#                       self.ddx_func(i1, i2, o1, o2, self.td_log_func, 'p22'),
#                       self.ddx_func(i1, i2, o1, o2, self.td_log_func, 'h22')] + z]]

        if self.ttd_u.is_set:
            mat_deriv += self.ttd_u_deriv(inl, outl)

        if self.ttd_l.is_set:
            mat_deriv += self.ttd_l_deriv(inl, outl)

        if self.pr1.is_set:
            pr1_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            pr1_deriv[0, 0, 1] = self.pr1.val
            pr1_deriv[0, 2, 1] = -1
            mat_deriv += pr1_deriv.tolist()

        if self.pr2.is_set:
            pr2_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            pr2_deriv[0, 1, 1] = self.pr2.val
            pr2_deriv[0, 3, 1] = -1
            mat_deriv += pr2_deriv.tolist()

        if self.zeta1.is_set:
            zeta1_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            for i in range(2):
                if i == 0:
                    zeta1_deriv[0, i * 2, 0] = (
                        self.ddx_func(inl, outl,
                                      self.zeta_func, 'm', i * 2))
                zeta1_deriv[0, i * 2, 1] = (
                    self.ddx_func(inl, outl,
                                  self.zeta_func, 'p', i * 2))
                zeta1_deriv[0, i * 2, 2] = (
                    self.ddx_func(inl, outl,
                                  self.zeta_func, 'h', i * 2))
            mat_deriv += zeta1_deriv.tolist()

        if self.zeta2.is_set:
            zeta2_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
            for i in range(2):
                if i == 0:
                    zeta2_deriv[0, i * 2 + 1, 0] = (
                        self.ddx_func(inl, outl,
                                      self.zeta2_func, 'm', i * 2 + 1))
                zeta2_deriv[0, i * 2 + 1, 1] = (
                    self.ddx_func(inl, outl,
                                  self.zeta2_func, 'p', i * 2 + 1))
                zeta2_deriv[0, i * 2 + 1, 2] = (
                    self.ddx_func(inl, outl,
                                  self.zeta2_func, 'h', i * 2 + 1))
            mat_deriv += zeta2_deriv.tolist()

        mat_deriv += self.additional_derivatives(nw)

        return np.asarray(mat_deriv)

    def additional_derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition for additional equations

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*list*) - matrix of partial derivatives
        """
        return []

    def zeta2_func(self, inl, outl):
        r"""
        calculates pressure drop from zeta2

        :param inlets: components connections at inlets
        :type inlets: list
        :param outlets: components connections at outlets
        :type outlets: list
        :returns: residual value for the pressure drop

        .. math::

            \zeta_2 = \frac{\Delta p_2 \cdot v_2 \cdot 2}{c_2^2}\\
            c_2 = \frac{\dot{m}_2 \cdot v_2}{A_2}

        As the cross sectional area A will not change from design to offdesign
        calculation, it is possible to handle this the following way:

        .. math::
            0 = \zeta_2 - \frac{(p_{2,in} - p_{2,out}) \cdot \pi^2}{8 \cdot
            \dot{m}_{2,in}^2 \cdot \frac{v_{2,in} + v_{2,out}}{2}}
        """
        i = inl[1].to_flow()
        o = outl[1].to_flow()
        return (self.zeta2.val - (i[1] - o[1]) * math.pi ** 2 /
                (8 * i[0] ** 2 * (v_mix_ph(i) + v_mix_ph(o)) / 2))

    def kA_func(self, inl, outl):
        r"""
        equation for heat flux from conditions on both sides of heat exchanger

        - calculate temperatures at inlets and outlets
        - perform convergence correction, if temperature levels do not
          match logic:

              * :math:`T_{1,in} > T_{2,out}`?
              * :math:`T_{1,out} < T_{2,in}`?

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*float*) - residual value of equation

        .. math::

            0 = \dot{m}_{1,in} \cdot \left( h_{1,out} - h_{1,in}\right) +
            kA \cdot \frac{T_{1,out} - T_{2,in} - T_{1,in} + T_{2,out}}
            {\ln{\frac{T_{1,out} - T_{2,in}}{T_{1,in} - T_{2,out}}}}
        """

        i1 = inl[0].to_flow()
        i2 = inl[1].to_flow()
        o1 = outl[0].to_flow()
        o2 = outl[1].to_flow()

        T_i1 = T_mix_ph(i1)
        T_i2 = T_mix_ph(i2)
        T_o1 = T_mix_ph(o1)
        T_o2 = T_mix_ph(o2)

        if T_i1 <= T_o2 and not inl[0].T.val_set:
            T_i1 = T_o2 + 2
        if T_i1 <= T_o2 and not outl[1].T.val_set:
            T_o2 = T_i1 - 1
        if T_i1 <= T_o2 and inl[0].T.val_set and outl[1].T.val_set:
            msg = ('Infeasibility at ' + str(self.label) + ': Upper '
                   'temperature difference is negative!')
            raise MyComponentError(msg)

        if T_o1 <= T_i2 and not outl[0].T.val_set:
            T_o1 = T_i2 + 1
        if T_o1 <= T_i2 and not inl[1].T.val_set:
            T_i2 = T_o1 - 1
        if T_o1 <= T_i2 and inl[1].T.val_set and outl[0].T.val_set:
            msg = ('Infeasibility at ' + str(self.label) + ': Lower '
                   'temperature difference is negative!')
            raise MyComponentError(msg)

        if self.kA_char1.param == 'm':
            expr = i1[0] / self.i10[0]
        else:
            expr = 1

        if expr > self.kA_char1.func.x[-1]:
            expr = self.kA_char1.func.x[-1]
        if expr < self.kA_char1.func.x[0]:
            expr = self.kA_char1.func.x[0]

        fkA1 = self.kA_char1.func.f_x(expr)

        if self.kA_char2.param == 'm':
            expr = i2[0] / self.i20[0]
        else:
            expr = 1

        if expr > self.kA_char2.func.x[-1]:
            expr = self.kA_char2.func.x[-1]
        if expr < self.kA_char2.func.x[0]:
            expr = self.kA_char2.func.x[0]

        fkA2 = self.kA_char2.func.f_x(expr)

        return (i1[0] * (o1[2] - i1[2]) + self.kA.val * fkA1 * fkA2 *
                (T_o1 - T_i2 - T_i1 + T_o2) /
                math.log((T_o1 - T_i2) / (T_i1 - T_o2)))

    def td_log_func(self, inl, outl):
        r"""
        equation for logarithmic temperature difference

        - calculate temperatures at inlets and outlets
        - perform convergence correction, if temperature levels do not
          match logic:

              * :math:`T_{1,in} > T_{2,out}`?
              * :math:`T_{1,out} < T_{2,in}`?

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*float*) - residual value of equation

        .. math::

            0 = td_{log} \cdot
            \frac{\ln{\frac{T_{1,out} - T_{2,in}}{T_{1,in} - T_{2,out}}}}
            {T_{1,out} - T_{2,in} - T_{1,in} + T_{2,out}}
        """

        i1 = inl[0].to_flow()
        i2 = inl[1].to_flow()
        o1 = outl[0].to_flow()
        o2 = outl[1].to_flow()

        T_i1 = T_mix_ph(i1)
        T_i2 = T_mix_ph(i2)
        T_o1 = T_mix_ph(o1)
        T_o2 = T_mix_ph(o2)

        if T_i1 <= T_o2 and not inl[0].T.val_set:
            T_i1 = T_o2 + 1
        if T_i1 <= T_o2 and not outl[1].T.val_set:
            T_o2 = T_i1 - 1
        if T_i1 <= T_o2 and inl[0].T.val_set and outl[1].T.val_set:
            msg = ('Infeasibility at ' + str(self.label) + ': Upper '
                   'temperature difference is negative!')
            raise MyComponentError(msg)

        if T_o1 <= T_i2 and not outl[0].T.val_set:
            T_o1 = T_i2 + 1
        if T_o1 <= T_i2 and not inl[1].T.val_set:
            T_i2 = T_o1 - 1
        if T_o1 <= T_i2 and inl[1].T.val_set and outl[0].T.val_set:
            msg = ('Infeasibility at ' + str(self.label) + ': Lower '
                   'temperature difference is negative!')
            raise MyComponentError(msg)

        return (self.td_log.val *
                math.log((T_o1 - T_i2) / (T_i1 - T_o2)) -
                T_o1 + T_i2 + T_i1 - T_o2)

    def ttd_u_func(self, inl, outl):
        r"""
        equation for upper terminal temperature difference

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*float*) - residual value of equation

        .. math::

            0 = ttd_{u} - T_{1,in} + T_{2,out}
        """
        i1 = inl[0].to_flow()
        o2 = outl[1].to_flow()
        return self.ttd_u.val - T_mix_ph(i1) + T_mix_ph(o2)

    def ttd_l_func(self, inl, outl):
        r"""
        equation for lower terminal temperature difference

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*float*) - residual value of equation

        .. math::

            0 = ttd_{l} - T_{1,out} + T_{2,in}
        """
        i2 = inl[1].to_flow()
        o1 = outl[0].to_flow()
        return self.ttd_l.val - T_mix_ph(o1) + T_mix_ph(i2)

    def ttd_u_deriv(self, inl, outl):
        r"""
        calculate matrix of partial derivatives towards pressure and
        enthalpy for upper terminal temperature equation

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*list*) - matrix of partial derivatives
        """
        deriv = np.zeros((1, 4, len(inl[0].fluid.val) + 3))
        for i in range(2):
            deriv[0, i * 3, 1] = (
                self.ddx_func(inl, outl, self.ttd_u_func, 'p', i * 3))
            deriv[0, i * 3, 2] = (
                self.ddx_func(inl, outl, self.ttd_u_func, 'h', i * 3))
        return deriv.tolist()

    def ttd_l_deriv(self, inl, outl):
        r"""
        calculate matrix of partial derivatives towards pressure and
        enthalpy for lower terminal temperature equation

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*list*) - matrix of partial derivatives
        """
        deriv = np.zeros((1, 4, len(inl[0].fluid.val) + 3))
        for i in range(2):
            deriv[0, i + 1, 1] = (
                self.ddx_func(inl, outl, self.ttd_l_func, 'p', i + 1))
            deriv[0, i + 1, 2] = (
                self.ddx_func(inl, outl, self.ttd_l_func, 'h', i + 1))
        return deriv.tolist()

    def convergence_check(self, nw):
        r"""
        prevent bad values for fluid properties in calculation

        - :math:`h_{1,in} > h_{1,out}`?
        - :math:`h_{2,in} < h_{2,out}`?

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: no return value
        """

        i, o = nw.comps.loc[self].i.tolist(), nw.comps.loc[self].o.tolist()

        if i[0].h.val_SI < o[0].h.val_SI and not o[0].h.val_set:
            o[0].h.val_SI = i[0].h.val_SI / 2
        if i[1].h.val_SI > o[1].h.val_SI and not i[1].h.val_set:
            i[1].h.val_SI = o[1].h.val_SI / 2

        if self.ttd_l.is_set:
            h_min_o1 = h_mix_pT(o[0].to_flow(), nw.T_range_SI[0])
            h_min_i2 = h_mix_pT(i[1].to_flow(), nw.T_range_SI[0])
            if not o[0].h.val_set and o[0].h.val_SI < h_min_o1 * 2:
                o[0].h.val_SI = h_min_o1 * 2
            if not i[1].h.val_set and i[1].h.val_SI < h_min_i2:
                i[1].h.val_SI = h_min_i2 * 1.1

        if self.ttd_u.is_set:
            h_min_i1 = h_mix_pT(i[0].to_flow(), nw.T_range_SI[0])
            h_min_o2 = h_mix_pT(o[1].to_flow(), nw.T_range_SI[0])
            if not i[0].h.val_set and i[0].h.val_SI < h_min_i1 * 2:
                i[0].h.val_SI = h_min_i1 * 2
            if not o[1].h.val_set and o[1].h.val_SI < h_min_o2:
                o[1].h.val_SI = h_min_o2 * 1.1

    def initialise_source(self, c, key):
        r"""
        returns a starting value for fluid properties at components outlets

        - set starting temperatures in a way, that they match required logic

              * :math:`T_{1,in} > T_{2,out}`?
              * :math:`T_{1,out} < T_{2,in}`?

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    outlets, :math:`val = 5 \cdot 10^6 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    outlets,

                      - outlet 1:
                        :math:`h = h(p,\;T=473.15 \text{K})`

                      - outlet 2:
                        :math:`h = h(p,\;T=523.15 \text{K})`
        """
        if key == 'p':
            return 50e5
        elif key == 'h':
            flow = [c.m.val0, c.p.val_SI, c.h.val_SI, c.fluid.val]
            if c.s_id == 'out1':
                T = 200 + 273.15
                return h_mix_pT(flow, T)
            else:
                T = 250 + 273.15
                return h_mix_pT(flow, T)
        else:
            return 0

    def initialise_target(self, c, key):
        r"""
        returns a starting value for fluid properties at components inlets

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - val (*float*) - starting value for pressure at components
                    inlets, :math:`val = 5 \cdot 10^6 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    inlets,

                      - inlet 1:
                        :math:`h = h(p,\;T=573.15 \text{K})`

                      - inlet 2:
                        :math:`h = h(p,\;T=493.15 \text{K})`
        """
        if key == 'p':
            return 50e5
        elif key == 'h':
            flow = [c.m.val0, c.p.val_SI, c.h.val_SI, c.fluid.val]
            if c.t_id == 'in1':
                T = 300 + 273.15
                return h_mix_pT(flow, T)
            else:
                T = 220 + 273.15
                return h_mix_pT(flow, T)
        else:
            return 0

    def calc_parameters(self, nw, mode):

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        if mode == 'pre':

            self.i10 = inl[0].to_flow()
            self.i20 = inl[1].to_flow()
            self.o10 = outl[0].to_flow()
            self.o20 = outl[1].to_flow()
            self.i10[3] = self.i10[3].copy()
            self.i20[3] = self.i20[3].copy()
            self.o10[3] = self.o10[3].copy()
            self.o20[3] = self.o20[3].copy()

        T_i2 = T_mix_ph(inl[1].to_flow())
        T_o1 = T_mix_ph(outl[0].to_flow())

        if isinstance(self, condenser):
            i1 = inl[0].to_flow()
            T_i1 = T_mix_ph([i1[0], i1[1], h_mix_pQ(i1, 1), i1[3]])
        else:
            T_i1 = T_mix_ph(inl[0].to_flow())
        T_o2 = T_mix_ph(outl[1].to_flow())
        if (mode == 'pre' and 'ttd_u' in self.offdesign) or mode == 'post':
            self.ttd_u.val = T_i1 - T_o2
        if (mode == 'pre' and 'ttd_l' in self.offdesign) or mode == 'post':
            self.ttd_l.val = T_o1 - T_i2

        if self.ttd_u.val < 0 or self.ttd_l.val < 0:
            msg = ('Invalid value for terminal temperature '
                   'difference.'
                   'ttd_u =', self.ttd_u.val,
                   'ttd_l =', self.ttd_l.val)
            print(msg)
            nw.errors += [self]

        if (mode == 'pre' and 'Q' in self.offdesign) or mode == 'post':
            self.Q.val = inl[0].m.val_SI * (outl[0].h.val_SI -
                                            inl[0].h.val_SI)

        if (mode == 'pre' and 'kA' in self.offdesign) or mode == 'post':
            if T_i1 <= T_o2 or T_o1 <= T_i2:
                self.td_log.val = np.nan
                self.kA.val = np.nan
            else:
                self.td_log.val = ((T_o1 - T_i2 - T_i1 + T_o2) /
                                   math.log((T_o1 - T_i2) / (T_i1 - T_o2)))
                self.kA.val = -(inl[0].m.val_SI * (
                                outl[0].h.val_SI - inl[0].h.val_SI) /
                                self.td_log.val)

        if (mode == 'pre' and 'pr1' in self.offdesign) or mode == 'post':
            self.pr1.val = outl[0].p.val_SI / inl[0].p.val_SI
        if (mode == 'pre' and 'pr2' in self.offdesign) or mode == 'post':
            self.pr2.val = outl[1].p.val_SI / inl[1].p.val_SI
        if (mode == 'pre' and 'zeta1' in self.offdesign) or mode == 'post':
            self.zeta1.val = ((inl[0].p.val_SI - outl[0].p.val_SI) *
                              math.pi ** 2 /
                              (8 * inl[0].m.val_SI ** 2 *
                              (v_mix_ph(inl[0].to_flow()) +
                               v_mix_ph(outl[0].to_flow())) / 2))
        if (mode == 'pre' and 'zeta2' in self.offdesign) or mode == 'post':
            self.zeta2.val = ((inl[1].p.val_SI - outl[1].p.val_SI) *
                              math.pi ** 2 /
                              (8 * inl[1].m.val_SI ** 2 *
                              (v_mix_ph(inl[1].to_flow()) +
                               v_mix_ph(outl[1].to_flow())) / 2))

        # improve this part (for heat exchangers only atm)
        if self.kA.is_set:
            expr = inl[0].m.val_SI / self.i10[0]
            if (expr > self.kA_char1.func.x[-1] or
                    expr < self.kA_char1.func.x[0]):
                msg = ('Warning: Expression for characteristics out of bounds:'
                       ' value is ' + str(expr))
                print(msg)
                nw.errors += [self]

            expr = inl[1].m.val_SI / self.i20[0]
            if (expr > self.kA_char2.func.x[-1] or
                    expr < self.kA_char2.func.x[0]):
                msg = ('Warning: Expression for characteristics out of bounds:'
                       ' value is ' + str(expr))
                print(msg)
                nw.errors += [self]

    def print_parameters(self, nw):

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        print('##### ', self.label, ' #####')
        if self.ttd_u.val < 0 and self.kA.is_set:
            print('!!!!! ERROR calculating heat exchanger: !!!!!\n'
                  'Negative value for TTD at given logarithmic temperature '
                  'difference or kA, result may be wrong.')
        print('Q = ', self.Q.val, 'W; '
              'ttd_u = ', self.ttd_u.val, 'K; '
              'ttd_l = ', self.ttd_l.val, 'K; '
              'td_log = ', self.td_log.val, 'K; '
              'kA = ', self.kA.val, 'W / K; '
              'pr1 = ', self.pr1.val, '; '
              'pr2 = ', self.pr2.val, '; '
              'zeta1 = ', self.zeta1.val, '; '
              'zeta2 = ', self.zeta2.val, '; '
              'm1 = ', inl[0].m.val_SI, 'kg / s; '
              'm2 = ', inl[1].m.val_SI, 'kg / s; '
              'Sirr = ', inl[1].m.val_SI * (s_mix_ph(outl[1].to_flow()) -
                                            s_mix_ph(inl[1].to_flow())) +
              inl[0].m.val_SI * (s_mix_ph(outl[0].to_flow()) -
                                 s_mix_ph(inl[0].to_flow())), 'W / K'
              )

# %%


class condenser(heat_exchanger):
    r"""

    - has additional equation for enthalpy at hot side outlet
    - pressure drop via zeta at hot side is not an offdesign parameter
    - has different calculation method for given heat transition coefficient
      and upper terminal temperature difference compared to parent class

    **available parameters**

    - Q: heat flux
    - kA: area independent heat transition coefficient,
      :math:`kA=\frac{\text{W}}{\text{K}}`
    - ttd_u: upper terminal temperature difference
    - ttd_l: lower terminal temperature difference
    - pr1: outlet to inlet pressure ratio at hot side
    - pr2: outlet to inlet pressure ratio at cold side
    - zeta1: geometry independent friction coefficient hot side
      :math:`[\zeta]=\frac{\text{Pa}}{\text{m}^4}`, also see
      :func:`tespy.components.components.component.zeta_func`
    - zeta2: geometry independent friction coefficient cold side
      :math:`[\zeta]=\frac{\text{Pa}}{\text{m}^4}`, also see
      :func:`tespy.components.components.heat_exchanger.zeta2_func`

    **equations**

    see tespy.components.components.heat_exchager.equations

    **default design parameters**

    - pr2, ttd_u, ttd_l

    **default offdesign parameters**

    - zeta2, kA (using kA_char1/kA_char2, method: COND_HOT/COND_COLD,
      param: m/m)

    **inlets and outlets**

    - in1, in2 (index 1: hot side, index 2: cold side)
    - out1, out2 (index 1: hot side, index 2: cold side)

    .. image:: _images/condenser.svg
       :scale: 100 %
       :alt: alternative text
       :align: center

    **Improvements**

    - see parent class
    """

    def component(self):
        return 'condenser'

    def attr_prop(self):
        return {'Q': dc_cp(), 'kA': dc_cp(), 'td_log': dc_cp(),
                'kA_char1': dc_cc(method='COND_HOT', param='m'),
                'kA_char2': dc_cc(method='COND_COLD', param='m'),
                'ttd_u': dc_cp(), 'ttd_l': dc_cp(),
                'pr1': dc_cp(), 'pr2': dc_cp(),
                'zeta1': dc_cp(), 'zeta2': dc_cp()}

    def default_design(self):
        return [n for n in heat_exchanger.default_design(self) if n != 'pr1']

    def default_offdesign(self):
        return [n for n in heat_exchanger.default_offdesign(self) if
                n != 'zeta1']

    def additional_equations(self, nw):
        r"""
        returns vector vec_res with result of additional equations for this
        component

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - vector of residual values

        **mandatory equations**

        .. math::

            0 = h_{1,out} - h\left(p, x=0 \right)\\
            x: \text{vapour mass fraction}
        """
        vec_res = []
        outl = nw.comps.loc[self].o.tolist()

        o1 = outl[0].to_flow()
        vec_res += [o1[2] - h_mix_pQ(o1, 0)]

        return vec_res

    def additional_derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition for additional equations

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*list*) - matrix of partial derivatives
        """

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inl), len(outl)
        num_fl = len(nw.fluids)
        mat_deriv = []

        o1 = outl[0].to_flow()
        x_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
        x_deriv[0, 2, 1] = -dh_mix_dpQ(o1, 0)
        x_deriv[0, 2, 2] = 1
        mat_deriv += x_deriv.tolist()

        return mat_deriv

    def kA_func(self, inl, outl):
        r"""
        equation for heat flux from conditions on both sides of heat exchanger

        - calculate temperatures at inlets and outlets
        - perform convergence correction, if temperature levels do not
          match logic:

              * :math:`T_{1,in} > T_{2,out}`?
              * :math:`T_{1,out} < T_{2,in}`?

        - kA refers to boiling temperature at hot side inlet

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*float*) - residual value of equation

        .. math::

            0 = ttd_{u} - T_s \left(p_{1,in}\right) + T_{2,out}

        .. math::

            0 = \dot{m}_{1,in} \cdot \left( h_{1,out} - h_{1,in}\right) +
            kA \cdot \frac{T_{1,out} - T_{2,in} - T_s \left(p_{1,in}\right) +
            T_{2,out}}
            {\ln{\frac{T_{1,out} - T_{2,in}}
            {T_s \left(p_{1,in}\right) - T_{2,out}}}}
        """

        i1 = inl[0].to_flow()
        i2 = inl[1].to_flow()
        o1 = outl[0].to_flow()
        o2 = outl[1].to_flow()

        T_i1 = T_mix_ph([i1[0], i1[1], h_mix_pQ(i1, 1), i1[3]])
        T_i2 = T_mix_ph(i2)
        T_o1 = T_mix_ph(o1)
        T_o2 = T_mix_ph(o2)

        if T_i1 <= T_o2 and not inl[0].T.val_set:
            T_i1 = T_o2 + 0.5
        if T_i1 <= T_o2 and not outl[1].T.val_set:
            T_o2 = T_i1 - 0.5

        if T_o1 <= T_i2 and not outl[0].T.val_set:
            T_o1 = T_i2 + 1
        if T_o1 <= T_i2 and not inl[1].T.val_set:
            T_i2 = T_o1 - 1

        if self.kA_char1.param == 'm':
            expr = i1[0] / self.i10[0]
        else:
            expr = 1

        if expr > self.kA_char1.func.x[-1]:
            expr = self.kA_char1.func.x[-1]
        if expr < self.kA_char1.func.x[0]:
            expr = self.kA_char1.func.x[0]

        fkA1 = self.kA_char1.func.f_x(expr)

        if self.kA_char2.param == 'm':
            expr = i2[0] / self.i20[0]
        else:
            expr = 1

        if expr > self.kA_char2.func.x[-1]:
            expr = self.kA_char2.func.x[-1]
        if expr < self.kA_char2.func.x[0]:
            expr = self.kA_char2.func.x[0]

        fkA2 = self.kA_char2.func.f_x(expr)

        return (i1[0] * (o1[2] - i1[2]) + self.kA.val * fkA1 * fkA2 *
                (T_o1 - T_i2 - T_i1 + T_o2) /
                math.log((T_o1 - T_i2) / (T_i1 - T_o2)))

    # function for logarithmic temperature difference not implemented
#    def td_log_func(self, inl, outl):
#
#        T_i1 = T_mix_ph([i1[0], i1[1], h_mix_pQ(i1, 1), i1[3]])
#        T_i2 = T_mix_ph(i2)
#        T_o1 = T_mix_ph(o1)
#        T_o2 = T_mix_ph(o2)
#
#        io2 = 0
#        while T_i1 <= T_o2:
#            try:
#                T_o2 = T_mix_ph([o2[0], o2[1], o2[2] - io2 * 10000, o2[3]])
#                io2 += 1
#            except:
#                None
#
#        i = 0
#        while T_o1 <= T_i2:
#            i += 1
#            T_o1 = T_mix_ph([o1[0], o1[1], o1[2] + i * 10000, o1[3]])
#
#        return (self.td_log *
#                math.log((T_o1 - T_i2) / (T_i1 - T_o2)) -
#                T_o1 + T_i2 + T_i1 - T_o2)

    def ttd_u_func(self, inl, outl):
        r"""
        equation for upper terminal temperature difference

        - ttd_u refers to boiling temperature at hot side inlet

        :param inlets: the components connections at the inlets
        :type inlets: list
        :param outlets: the components connections at the outlets
        :type outlets: list
        :returns: val (*float*) - residual value of equation

        .. math::

            0 = ttd_{u} - T_s \left(p_{1,in}\right) + T_{2,out}
        """
        i1 = inl[0].to_flow()
        o2 = outl[1].to_flow()
        return (self.ttd_u.val -
                T_mix_ph([i1[0], i1[1], h_mix_pQ(i1, 1), i1[3]]) +
                T_mix_ph(o2))

# %%


class desuperheater(heat_exchanger):
    r"""

    - has additional equation for enthalpy at hot side outlet

    **available parameters**

    - Q: heat flux
    - kA: area independent heat transition coefficient,
      :math:`kA=\frac{\text{W}}{\text{K}}`
    - ttd_u: upper terminal temperature difference
    - ttd_l: lower terminal temperature difference
    - pr1: outlet to inlet pressure ratio at hot side
    - pr2: outlet to inlet pressure ratio at cold side
    - zeta1: geometry independent friction coefficient hot side
      :math:`[\zeta]=\frac{\text{Pa}}{\text{m}^4}`, also see
      :func:`tespy.components.components.component.zeta_func`
    - zeta2: geometry independent friction coefficient cold side
      :math:`[\zeta]=\frac{\text{Pa}}{\text{m}^4}`, also see
      :func:`tespy.components.components.heat_exchanger.zeta2_func`

    **equations**

    see tespy.components.components.heat_exchager.equations

    **default design parameters**

    - pr1, pr2, ttd_u, ttd_l

    **default offdesign parameters**

    - zeta1, zeta2, kA

    **inlets and outlets**

    - in1, in2 (index 1: hot side, index 2: cold side)
    - out1, out2 (index 1: hot side, index 2: cold side)

    .. image:: _images/heat_exchanger.svg
       :scale: 100 %
       :alt: alternative text
       :align: center

    **Improvements**

    - see parent class
    """

    def component(self):
        return 'desuperheater'

    def default_design(self):
        return heat_exchanger.default_design(self)

    def default_offdesign(self):
        return heat_exchanger.default_offdesign(self)

    def additional_equations(self, nw):
        r"""
        returns vector vec_res with result of additional equations for this
        component

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - vector of residual values

        **mandatory equations**

        .. math::

            0 = h_{1,out} - h\left(p, x=1 \right)\\
            x: \text{vapour mass fraction}
        """

        vec_res = []
        outl = nw.comps.loc[self].o.tolist()

        o1 = outl[0].to_flow()
        vec_res += [o1[2] - h_mix_pQ(o1, 1)]

        return vec_res

    def additional_derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition for additional equations

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*list*) - matrix of partial derivatives
        """

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inl), len(outl)
        num_fl = len(nw.fluids)
        mat_deriv = []

        o1 = outl[0].to_flow()
        x_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
        x_deriv[0, 2, 1] = -dh_mix_dpQ(o1, 1)
        x_deriv[0, 2, 2] = 1
        mat_deriv += x_deriv.tolist()

        return mat_deriv


# %%


class drum(component):
    r"""

    - assumes, that the fluid composition between outlet 1 and inlet 2 does
      not change!

    **no specification parameters available**

    **equations**

    see tespy.components.components.drum.equations

    **inlets and outlets**

    - in1, in2 (index 1: from economiser, index 2: from evaporator)
    - out1, out2 (index 1: to evaporator, index 2: to superheater)

    .. image:: _images/drum.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
    """

    def inlets(self):
        return ['in1', 'in2']

    def outlets(self):
        return ['out1', 'out2']

    def component(self):
        return 'drum'

    def equations(self, nw):
        r"""
        returns vector vec_res with result of equations for this component

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - vector of residual values

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_res`
        - :func:`tespy.components.components.component.mass_flow_res`

        .. math::

            0 = \sum_i \left(\dot{m}_{i,in} \cdot h_{i,in} \right) -
            \sum_j \left(\dot{m}_{j,out} \cdot h_{j,out} \right)\;
            \forall i \in inlets, \; \forall j \in outlets\\
            0 = h_{1,out} - h\left(p, x=0 \right)\\
            0 = h_{2,out} - h\left(p, x=1 \right)\\
            x: \text{vapour mass fraction}

        """

        vec_res = []
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        vec_res += self.fluid_res(inl, outl)
        vec_res += self.mass_flow_res(inl, outl)

        E_res = 0
        for i in inl:
            E_res += i.m.val_SI * i.h.val_SI
        for o in outl:
            E_res -= o.m.val_SI * o.h.val_SI
        vec_res += [E_res]

        p = inl[0].p.val_SI
        for c in [inl[1]] + outl:
            vec_res += [p - c.p.val_SI]

        vec_res += [h_mix_pQ(outl[0].to_flow(), 0) - outl[0].h.val_SI]
        vec_res += [h_mix_pQ(outl[1].to_flow(), 1) - outl[1].h.val_SI]

        return vec_res

    def derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*numpy array*) - matrix of partial derivatives
        """

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inl), len(outl)
        num_fl = len(nw.fluids)
        mat_deriv = []

        mat_deriv += self.fluid_deriv(inl, outl)
        mat_deriv += self.mass_flow_deriv(inl, outl)

        E_deriv = np.zeros((1, num_i + num_o, num_fl + 3))
        k = 0
        for i in inl:
            E_deriv[0, k, 0] = i.h.val_SI
            E_deriv[0, k, 2] = i.m.val_SI
            k += 1
        j = 0
        for o in outl:
            E_deriv[0, j + k, 0] = -o.h.val_SI
            E_deriv[0, j + k, 2] = -o.m.val_SI
            j += 1
        mat_deriv += E_deriv.tolist()

        p_deriv = np.zeros((num_i + num_o - 1, num_i + num_o, num_fl + 3))
        for k in range(num_i + num_o - 1):
            p_deriv[k, 0, 1] = 1
            p_deriv[k, k + 1, 1] = -1
        mat_deriv += p_deriv.tolist()

        o1 = outl[0].to_flow()
        o2 = outl[1].to_flow()

        x_deriv = np.zeros((num_o, num_i + num_o, num_fl + 3))
        x_deriv[0, 2, 1] = dh_mix_dpQ(o1, 0)
        x_deriv[0, 2, 2] = -1
        x_deriv[1, 3, 1] = dh_mix_dpQ(o2, 1)
        x_deriv[1, 3, 2] = -1
        mat_deriv += x_deriv.tolist()

        return np.asarray(mat_deriv)

    def initialise_source(self, c, key):
        r"""
        returns a starting value for fluid properties at components outlets

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    outlets, :math:`val = 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    outlets,

                      - outlet 1:
                        :math:`h = h(p,\;x=0)`

                      - outlet 2:
                        :math:`h = h(p,\;x=1)`
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            if c.s_id == 'out1':
                return h_mix_pQ(c.to_flow(), 0)
            else:
                return h_mix_pQ(c.to_flow(), 1)
        else:
            return 0

    def initialise_target(self, c, key):
        r"""
        returns a starting value for fluid properties at components inlets

        :param c: connection to apply initialisation
        :type c: tespy.connections.connection
        :param key: property
        :type key: str
        :returns: - p (*float*) - starting value for pressure at components
                    inlets, :math:`val = 10^5 \; \text{Pa}`
                  - h (*float*) - starting value for enthalpy at components
                    inlets,

                      - inlet 1:
                        :math:`h = h(p,\;x=0)`

                      - inlet 2:
                        :math:`h = h(p,\;x=0.7)`
        """
        if key == 'p':
            return 10e5
        elif key == 'h':
            if c.t_id == 'in1':
                return h_mix_pQ(c.to_flow(), 0)
            else:
                return h_mix_pQ(c.to_flow(), 0.7)
        else:
            return 0

# %%


class subsys_interface(component):
    r"""
    interface for subsystems

    - passes fluid properties/flow information at inlet i to outlet i
    - no transformation of any fluid properties

    **available parameters**

    - num_inter: number of connections for the interface

    **equations**

    see tespy.components.components.subsys_interface.equations

    **inlets and outlets**

    - specify number of inlets and outlets with :code:`num_inter`
    - predefined value: 1

    .. image:: _images/subsys_interface.svg
       :scale: 100 %
       :alt: alternative text
       :align: center
    """

    def attr(self):
        return ['num_inter']

    def attr_prop(self):
        return {'num_inter': dc_cp()}

    def inlets(self):
        if self.num_inter.is_set:
            return ['in' + str(i + 1) for i in range(self.num_inter.val)]
        else:
            return ['in1']

    def outlets(self):
        if self.num_inter.is_set:
            return ['out' + str(i + 1) for i in range(self.num_inter.val)]
        else:
            return ['out1']

    def component(self):
        return 'subsystem interface'

    def equations(self, nw):
        r"""
        returns vector vec_res with result of equations for this component

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: vec_res (*list*) - vector of residual values

        **mandatory equations**

        - :func:`tespy.components.components.component.fluid_res`
        - :func:`tespy.components.components.component.mass_flow_res`

        .. math::

            0 = p_{i,in} - p_{i,out}\\
            0 = h_{i,in} - h_{i,out}\\
            \forall i \in inlets/outlets

        """

        vec_res = []
        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())

        vec_res += self.fluid_res(inl, outl)
        vec_res += self.mass_flow_res(inl, outl)
        for j in range(len(inl)):
            i = inl[j]
            o = outl[j]
            vec_res += [i.p.val_SI - o.p.val_SI]
        for j in range(len(inl)):
            i = inl[j]
            o = outl[j]
            vec_res += [i.h.val_SI - o.h.val_SI]

        return vec_res

    def derivatives(self, nw):
        r"""
        calculate matrix of partial derivatives towards mass flow, pressure,
        enthalpy and fluid composition

        :param nw: network using this component object
        :type nw: tespy.networks.network
        :returns: mat_deriv (*numpy array*) - matrix of partial derivatives
        """

        inl, outl = (nw.comps.loc[self].i.tolist(),
                     nw.comps.loc[self].o.tolist())
        num_i, num_o = len(inl), len(outl)
        num_fl = len(nw.fluids)
        mat_deriv = []

        mat_deriv += self.fluid_deriv(inl, outl)
        mat_deriv += self.mass_flow_deriv(inl, outl)

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
