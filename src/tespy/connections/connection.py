# -*- coding: utf-8
"""Module of class Connection and class Ref.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/connections/connection.py

SPDX-License-Identifier: MIT
"""

from collections import OrderedDict

import numpy as np

from tespy.components.component import Component
from tespy.tools import fluid_properties as fp
from tespy.tools import logger
from tespy.tools.data_containers import FluidComposition as dc_flu
from tespy.tools.data_containers import FluidProperties as dc_prop
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import CoolPropWrapper
from tespy.tools.fluid_properties import Q_mix_ph
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import T_sat_p
from tespy.tools.fluid_properties import dT_mix_dph
from tespy.tools.fluid_properties import dT_mix_pdh
from tespy.tools.fluid_properties import dT_sat_dp
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import s_mix_ph
from tespy.tools.fluid_properties import v_mix_ph
from tespy.tools.fluid_properties.helpers import get_number_of_fluids
from tespy.tools.global_vars import ERR
from tespy.tools.global_vars import fluid_property_data as fpd
from tespy.tools.helpers import TESPyConnectionError
from tespy.tools.helpers import convert_from_SI


class Connection:
    r"""
    Class connection is the container for fluid properties between components.

    Parameters
    ----------
    m : float, tespy.connections.connection.Ref
        Mass flow specification.

    m0 : float
        Starting value specification for mass flow.

    p : float, tespy.connections.connection.Ref
        Pressure specification.

    p0 : float
        Starting value specification for pressure.

    h : float, tespy.connections.connection.Ref
        Enthalpy specification.

    h0 : float
        Starting value specification for enthalpy.

    fluid : dict
        Fluid compostition specification.

    fluid0 : dict
        Starting value specification for fluid compostition.

    fluid_balance : boolean
        Fluid balance equation specification.

    x : float
        Gas phase mass fraction specification.

    T : float, tespy.connections.connection.Ref
        Temperature specification.

    Td_bp : float
        Temperature difference to boiling point at pressure corresponding
        pressure of this connection in K.

    v : float
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

    printout : boolean
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

    >>> from tespy.components import Sink, Source
    >>> from tespy.connections import Connection, Ref
    >>> from tespy.tools import FluidComposition as dc_flu
    >>> from tespy.tools import FluidProperties as dc_prop
    >>> import numpy as np
    >>> so1 = Source('source1')
    >>> so2 = Source('source2')
    >>> si1 = Sink('sink1')
    >>> si2 = Sink('sink2')
    >>> so_si1 = Connection(so1, 'out1', si1, 'in1', label='connection 1')
    >>> so_si2 = Connection(so2, 'out1', si2, 'in1')
    >>> so_si1.label
    'connection 1'
    >>> so_si2.label
    'source2:out1_sink2:in1'

    There are different ways of setting parameters on connections: Specify

    - a numeric value  (for attributes mass flow, pressure and enthalpy)
    - a numeric starting value (for attributes mass flow, pressure and
      enthalpy)
    - a dictionary (for attributes fluid and fluid0)
    - a boolean value (for attributes fluid_balance, local_design,
      local_offdesign).
    - a referenced value (mass flow, pressure, temperature, enthalpy).
    - numpy.nan or None (unsetting a value).
    - a string (for attributes design_path and state).
    - a list (for attributes design and offdesign).

    >>> so_si1.set_attr(v=0.012, m0=10, p=5, h=400, fluid={'H2O': 1, 'N2': 0})
    >>> so_si2.set_attr(m=Ref(so_si1, 2, -5), h0=700, T=200,
    ... fluid={'N2': 1}, fluid_balance=True,
    ... design=['T'], offdesign=['m', 'v'])

    The set_attr method automatically converts your input in data_container
    information.

    >>> type(so_si1.v)
    <class 'tespy.tools.data_containers.FluidProperties'>
    >>> type(so_si1.fluid)
    <class 'tespy.tools.data_containers.FluidComposition'>

    If you want get a spcific value use the logic: connection.property.*.
    Aditionally, it is possible to use the :code:`get_attr` method.

    >>> so_si1.m.val0
    10
    >>> so_si1.m.val_set
    False
    >>> so_si1.m.get_attr('val_set')
    False
    >>> type(so_si2.m.ref)
    <class 'tespy.connections.connection.Ref'>
    >>> so_si2.fluid.get_attr('balance')
    True
    >>> so_si2.m.ref.get_attr('delta')
    -5
    >>> so_si2.m.ref_set
    True
    >>> type(so_si2.m.ref.get_attr('obj'))
    <class 'tespy.connections.connection.Connection'>

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
    """

    def __init__(self, comp1, outlet_id, comp2, inlet_id,
                 label=None, **kwargs):

        # check input parameters
        if not (isinstance(comp1, Component) and
                isinstance(comp2, Component)):
            msg = ('Error creating connection. Check if comp1, comp2 are of '
                   'type component.')
            logger.error(msg)
            raise TypeError(msg)

        if comp1 == comp2:
            msg = ('Error creating connection. Cannot connect component ' +
                   comp1.label + ' to itself.')
            logger.error(msg)
            raise TESPyConnectionError(msg)

        if outlet_id not in comp1.outlets():
            msg = ('Error creating connection. Specified oulet_id (' +
                   outlet_id + ') is not valid for component ' +
                   comp1.component() + '. Valid ids are: ' +
                   str(comp1.outlets()) + '.')
            logger.error(msg)
            raise ValueError(msg)

        if inlet_id not in comp2.inlets():
            msg = (
                'Error creating connection. Specified inlet_id (' + inlet_id +
                ') is not valid for component ' + comp2.component() +
                '. Valid ids are: ' + str(comp2.inlets()) + '.')
            logger.error(msg)
            raise ValueError(msg)

        if label is None:
            self.label = (
                comp1.label + ':' + outlet_id + '_' +
                comp2.label + ':' + inlet_id)
        else:
            self.label = label

        if not isinstance(self.label, str):
            msg = 'Please provide the label as string.'
            logger.error(msg)
            raise TypeError(msg)

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
        self.property_data = self.get_parameters()
        self.parameters = {
            k: v for k, v in self.get_parameters().items()
            if isinstance(v, dc_prop) and v.func is not None
        }
        self.property_data0 = [x + '0' for x in self.property_data.keys()]
        self.__dict__.update(self.property_data)
        self.mixing_rule = None
        self.set_attr(**kwargs)

        msg = (
            'Created connection ' + self.source.label + ' (' + self.source_id +
            ') -> ' + self.target.label + ' (' + self.target_id + ').')
        logger.debug(msg)

    def set_attr(self, **kwargs):
        r"""
        Set, reset or unset attributes of a connection.

        Parameters
        ----------
        m : float, tespy.connections.connection.Ref
            Mass flow specification.

        m0 : float
            Starting value specification for mass flow.

        p : float, tespy.connections.connection.Ref
            Pressure specification.

        p0 : float
            Starting value specification for pressure.

        h : float, tespy.connections.connection.Ref
            Enthalpy specification.

        h0 : float
            Starting value specification for enthalpy.

        fluid : dict
            Fluid compostition specification.

        fluid0 : dict
            Starting value specification for fluid compostition.

        fluid_balance : boolean
            Fluid balance equation specification.

        x : float
            Gas phase mass fraction specification.

        T : float, tespy.connections.connection.Ref
            Temperature specification.

        Td_bp : float
            Temperature difference to boiling point at pressure corresponding
            pressure of this connection in K.

        v : float
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

        printout : boolean
            Include this connection in the network's results printout.

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
            if key == 'label':
                # bad datatype
                msg = 'Label can only be specified on instance creation.'
                logger.error(msg)
                raise TESPyConnectionError(msg)
            elif key in self.property_data or key in self.property_data0:
                # fluid specification
                try:
                    float(kwargs[key])
                    is_numeric = True
                except (TypeError, ValueError):
                    is_numeric = False
                if 'fluid' in key and key != 'fluid_balance':
                    if isinstance(kwargs[key], dict):
                        # starting values
                        if key in self.property_data0:
                            self.fluid.set_attr(val0=kwargs[key].copy())
                        # specified parameters
                        else:
                            self.fluid.set_attr(val=kwargs[key].copy())
                            self.fluid.set_attr(
                                val_set={f: True for f in kwargs[key].keys()})

                    else:
                        # bad datatype
                        msg = (
                            'Datatype for fluid vector specification must be '
                            'dict.')
                        logger.error(msg)
                        raise TypeError(msg)

                elif key == 'state':
                    if kwargs[key] in ['l', 'g']:
                        self.state.set_attr(val=kwargs[key], is_set=True)
                    elif kwargs[key] is None:
                        self.state.set_attr(is_set=False)
                    elif is_numeric:
                        if np.isnan(kwargs[key]):
                            self.get_attr(key).set_attr(is_set=False)
                        else:
                            msg = (
                                'To unset the state specification either use '
                                'np.nan or None.')
                            logger.error(msg)
                            raise ValueError(msg)
                    else:
                        msg = (
                            'Keyword argument "state" must either be '
                            '"l" or "g" or be None or np.nan.')
                        logger.error(msg)
                        raise TypeError(msg)

                elif kwargs[key] is None:
                    self.get_attr(key).set_attr(val_set=False)
                    if f"{key}_ref" in self.property_data:
                        self.get_attr(key).set_attr(val_set=False)
                    if key in ["m", "p", "h"]:
                        self.get_attr(key).is_var = True

                elif is_numeric:
                    if np.isnan(kwargs[key]):
                        self.get_attr(key).set_attr(val_set=False)
                        if f"{key}_ref" in self.property_data:
                            self.get_attr(key).set_attr(val_set=False)
                        if key in ["m", "p", "h"]:
                            self.get_attr(key).is_var = True
                    else:
                        # value specification
                        if key in self.property_data:
                            self.get_attr(key).set_attr(
                                val_set=True,
                                val=kwargs[key])
                            if key in ["m", "p", "h"]:
                                self.get_attr(key).is_var = False
                        # starting value specification
                        else:
                            self.get_attr(key.replace('0', '')).set_attr(
                                val0=kwargs[key])

                # reference object
                elif isinstance(kwargs[key], Ref):
                    if f"{key}_ref" not in self.property_data:
                        msg = f"Referencing {key} is not implemented."
                        logger.error(msg)
                        raise NotImplementedError(msg)
                    else:
                        self.get_attr(f"{key}_ref").set_attr(val=kwargs[key])
                        self.get_attr(f"{key}_ref").set_attr(val_set=True)

                # invalid datatype for keyword
                else:
                    msg = 'Bad datatype for keyword argument ' + key + '.'
                    logger.error(msg)
                    raise TypeError(msg)

            # fluid balance
            elif key == 'fluid_balance':
                if isinstance(kwargs[key], bool):
                    self.get_attr('fluid').set_attr(balance=kwargs[key])
                else:
                    msg = (
                        'Datatype for keyword argument fluid_balance must be '
                        'boolean.')
                    logger.error(msg)
                    raise TypeError(msg)

            # fluid variable
            elif key == 'fluid_variable':
                if isinstance(kwargs[key], bool):
                    self.get_attr('fluid').set_attr(is_var=kwargs[key])
                else:
                    msg = (
                        'Datatype for keyword argument fluid_variable must be '
                        'boolean.')
                    logger.error(msg)
                    raise TypeError(msg)

            # design/offdesign parameter list
            elif key == 'design' or key == 'offdesign':
                if not isinstance(kwargs[key], list):
                    msg = 'Please provide the ' + key + ' parameters as list!'
                    logger.error(msg)
                    raise TypeError(msg)
                elif set(kwargs[key]).issubset(self.property_data.keys()):
                    self.__dict__.update({key: kwargs[key]})
                else:
                    params = ', '.join(self.property_data.keys())
                    msg = (
                        'Available parameters for (off-)design specification '
                        'are: ' + params + '.')
                    logger.error(msg)
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
                    logger.error(msg)
                    raise TypeError(msg)

                self.new_design = True

            # other boolean keywords
            elif key in ['printout', 'local_design', 'local_offdesign']:
                if not isinstance(kwargs[key], bool):
                    msg = ('Please provide the ' + key + ' as boolean.')
                    logger.error(msg)
                    raise TypeError(msg)
                else:
                    self.__dict__.update({key: kwargs[key]})

            elif key == "mixing_rule":
                self.mixing_rule = kwargs[key]

            # invalid keyword
            else:
                msg = 'Connection has no attribute ' + key + '.'
                logger.error(msg)
                raise KeyError(msg)

        if "fluid" in kwargs:
            for fluid in kwargs["fluid"]:
                engine = CoolPropWrapper
                back_end = "HEOS"
                self._create_fluid_wrapper(fluid, engine, back_end)

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
            logger.error(msg)
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
        return {
            'm': dc_prop(), 'p': dc_prop(), 'h': dc_prop(), 'T': dc_prop(),
            'x': dc_prop(), 'v': dc_prop(), 'vol': dc_prop(),
            's': dc_prop(),
            'fluid': dc_flu(), 'Td_bp': dc_prop(), 'state': dc_simple()
        }

    def _create_fluid_wrapper(self, fluid, engine, back_end):
        self.fluid.back_end[fluid] = back_end
        self.fluid.engine[fluid] = engine
        self.fluid.wrapper[fluid] = engine(fluid, back_end)

    def preprocess(self):
        self.num_eq = 0
        self.it = 0
        for parameter in self.parameters:
            if self.get_attr(parameter).val_set:
                self.num_eq += self.parameters[parameter].num_eq

        self.residual = np.zeros(self.num_eq)
        self.jacobian = OrderedDict()

    def get_parameters(self):
        return {
            'm': dc_prop(is_var=True),
            'p': dc_prop(is_var=True),
            'h': dc_prop(is_var=True),
            'vol': dc_prop(),
            's': dc_prop(),
            'fluid': dc_flu(),
            'state': dc_simple(),
            "T": dc_prop(**{
                "func": self.T_func, "deriv": self.T_deriv, "num_eq": 1
            }),
            "v": dc_prop(**{
                "func": self.v_func, "deriv": self.v_deriv, "num_eq": 1
            }),
            "x": dc_prop(**{
                "func": self.x_func, "deriv": self.x_deriv, "num_eq": 1
            }),
            "Td_bp": dc_prop(**{
                "func": self.Td_bp_func, "deriv": self.Td_bp_deriv, "num_eq": 1
            }),
            "m_ref": dc_prop(**{
                "func": self.primary_ref_func, "deriv": self.primary_ref_deriv,
                "num_eq": 1, "func_params": {"variable": "m"}
            }),
            "p_ref": dc_prop(**{
                "func": self.primary_ref_func, "deriv": self.primary_ref_deriv,
                "num_eq": 1, "func_params": {"variable": "p"}
            }),
            "h_ref": dc_prop(**{
                "func": self.primary_ref_func, "deriv": self.primary_ref_deriv,
                "num_eq": 1, "func_params": {"variable": "h"}
            }),
        }

    def build_fluid_data(self):
        self.fluid_data = {
            fluid: {
                "wrapper": self.fluid.wrapper[fluid],
                "mass_fraction": self.fluid.val[fluid]
            } for fluid in self.fluid.val
        }

    def primary_ref_func(self, k, **kwargs):
        variable = kwargs['variable']
        self.get_attr(variable)
        ref = self.get_attr(f"{variable}_ref").val
        self.residual[k] = (
            self.get_attr(variable).val_SI
            - ref.obj.get_attr(variable).val_SI * ref.factor + ref.delta
        )

    def primary_ref_deriv(self, k, **kwargs):
        variable = kwargs['variable']
        ref = self.get_attr(f"{variable}_ref").val
        if self.get_attr(variable).is_var:
            self.jacobian[k, self.get_attr(variable).J_col] = 1

        if ref.obj.get_attr(variable).is_var:
            self.jacobian[k, ref.obj.get_attr(variable).J_col] = -ref.factor

    def calc_T(self):
        return T_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data, self.mixing_rule)

    def T_func(self, k, **kwargs):
        self.residual[k] = self.T.val_SI - self.calc_T()

    def T_deriv(self, k, **kwargs):
        if self.p.is_var:
            self.jacobian[k, self.p.J_col] = (
                -dT_mix_dph(self.p.val_SI, self.h.val_SI, self.fluid_data, self.mixing_rule)
            )
        if self.h.is_var:
            self.jacobian[k, self.h.J_col] = (
                -dT_mix_pdh(self.p.val_SI, self.h.val_SI, self.fluid_data, self.mixing_rule)
            )
        # if len(self.fluid.val) > 1:
        #     self.jacobian[0, 0, 3:] = dT_mix_ph_dfluid(
        #         c.p.val_SI, c.h.val_SI, self.fluid_data, T0=c.T.val_SI
        # )

    def calc_vol(self):
        try:
            vol = v_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data, self.mixing_rule, T0=self.T.val_SI)
        except NotImplementedError:
            vol = np.nan
        return vol

    def v_func(self, k, **kwargs):
        self.residual[k] = (
            v_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data)
            * self.m.val_SI - self.v.val_SI
        )

    def v_deriv(self, k, **kwargs):
        if self.m.is_var:
            self.jacobian[k, self.m.J_col] = v_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data)
        if self.p.is_var:
            self.jacobian[k, self.p.J_col] = dv_mix_dph(self.p.val_SI, self.h.val_SI, self.fluid_data) * self.m.val_SI
        if self.h.is_var:
            self.jacobian[k, self.h.J_col] = dv_mix_pdh(self.p.val_SI, self.h.val_SI, self.fluid_data) * self.m.val_SI

    def calc_x(self):
        try:
            x = Q_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data)
        except NotImplementedError:
            x = np.nan
        return x

    def x_func(self, k, **kwargs):
        # saturated steam fraction
        # how to check if condition?
        # if (np.absolute(self.residual[k]) > ERR ** 2 or
        #         self.iter % 2 == 0 or self.always_all_equations):
        self.residual[k] = self.h.val_SI - h_mix_pQ(self.p.val_SI, self.x.val_SI, self.fluid_data)

    def x_deriv(self, k, **kwargs):
        # if not self.increment_filter[col + 1]:
        if self.p.is_var:
            self.jacobian[k, self.p.J_col] = -dh_mix_dpQ(self.p.val_SI, self.x.val_SI, self.fluid_data)
        if self.h.is_var: # and self.it == 0
            self.jacobian[k, self.h.J_col] = 1

    def calc_Td_bp(self):
        try:
            Td_bp = self.calc_T() - T_sat_p(self.p.val_SI, self.fluid_data)
        except NotImplementedError:
            Td_bp = np.nan
        return Td_bp

    def Td_bp_func(self, k, **kwargs):

        # temperature difference to boiling point
            # if (np.absolute(self.residual[k]) > ERR ** 2 or
            #         self.iter % 2 == 0 or self.always_all_equations):
        self.residual[k] = self.calc_Td_bp() - self.Td_bp.val_SI

    def Td_bp_deriv(self, k, **kwargs):
        # if not self.increment_filter[col + 1]:
        if self.p.is_var:
            self.jacobian[k, self.p.J_col] = (
                dT_mix_dph(self.p.val_SI, self.h.val_SI, self.fluid_data)
                - dT_sat_dp(self.p.val_SI, self.fluid_data)
            )
        if self.h.is_var:
            # if not self.increment_filter[col + 2]:
            self.jacobian[k, self.h.J_col] = dT_mix_pdh(
                self.p.val_SI, self.h.val_SI, self.fluid_data
            )

    def calc_s(self):
        return s_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data, self.mixing_rule, T0=self.T.val_SI)

    def solve(self, increment_filter):
        k = 0
        for parameter, data in self.parameters.items():
            if self.get_attr(parameter).val_set:
                data.func(k, **data.func_params)
                data.deriv(k, **data.func_params)

                k += 1

    def calc_results(self):
        self.T.val_SI = self.calc_T()
        number_fluids = get_number_of_fluids(self.fluid_data)
        if (
                number_fluids > 1
                and abs(
                    h_mix_pT(self.p.val_SI, self.T.val_SI, self.fluid_data, self.mixing_rule)
                    - self.h.val_SI
                ) > ERR ** .5
        ):
            self.T.val_SI = np.nan
            self.vol.val_SI = np.nan
            self.v.val_SI = np.nan
            self.s.val_SI = np.nan
            msg = (
                'Could not find a feasible value for mixture temperature '
                'at connection ' + self.label + '. The values for '
                'temperature, specific volume, volumetric flow and '
                'entropy are set to nan.')
            logger.error(msg)
        else:
            self.vol.val_SI = self.calc_vol()
            self.v.val_SI = self.vol.val_SI * self.m.val_SI
            self.s.val_SI = self.calc_s()

        if number_fluids == 1:
            if not self.x.val_set:
                self.x.val_SI = self.calc_x()
            if not self.Td_bp.val_set:
                self.Td_bp.val_SI = self.calc_Td_bp()

            for prop in fpd.keys():
                self.get_attr(prop).val = convert_from_SI(
                    prop, self.get_attr(prop).val_SI, self.get_attr(prop).unit
                )

        self.m.val0 = self.m.val
        self.p.val0 = self.p.val
        self.h.val0 = self.h.val
        self.fluid.val0 = self.fluid.val.copy()

    def check_pressure_bounds(self, fluid):
        if self.p.val_SI < self.fluid.wrapper[fluid]._p_min:
            c.p.val_SI = self.fluid.wrapper[fluid]._p_min * 1.01
            logger.debug(self._property_range_message('p'))
        elif self.p.val_SI > self.fluid.wrapper[fluid]._p_max:
            self.p.val_SI = self.fluid.wrapper[fluid]._p_max
            logger.debug(self._property_range_message('p'))

    def check_enthalpy_bounds(self, fluid):
        # enthalpy
        try:
            hmin = self.fluid.wrapper[fluid].h_pT(
                self.p.val_SI, self.fluid.wrapper[fluid]._T_min
            )
        except ValueError:
            f = 1.05
            hmin = self.fluid.wrapper[fluid].h_pT(
                self.p.val_SI, self.fluid.wrapper[fluid]._T_min * f
            )

        T = self.fluid.wrapper[fluid]._T_max
        while True:
            try:
                hmax = self.fluid.wrapper[fluid].h_pT(self.p.val_SI, T)
                break
            except ValueError as e:
                T *= 0.99
                if T < self.fluid.wrapper[fluid]._T_min:
                    raise ValueError(e) from e

        if self.h.val_SI < hmin:
            if hmin < 0:
                self.h.val_SI = hmin * 0.9999
            else:
                self.h.val_SI = hmin * 1.0001
            logger.debug(self._property_range_message('h'))

        elif self.h.val_SI > hmax:
            self.h.val_SI = hmax * 0.9999
            logger.debug(self._property_range_message('h'))

    def check_two_phase_bounds(self, fluid):

        if (self.Td_bp.val_SI > 0 or (self.state.val == 'g' and self.state.is_set)):
            h = self.fluid.wrapper[fluid].h_pQ(self.p.val_SI, 1)
            if self.h.val_SI < h:
                self.h.val_SI = h * 1.01
                logger.debug(self._property_range_message('h'))
        elif (self.Td_bp.val_SI < 0 or (self.state.val == 'l' and self.state.is_set)):
            h = self.fluid.wrapper[fluid].h_pQ(self.p.val_SI, 0)
            if self.h.val_SI > h:
                self.h.val_SI = h * 0.99
                logger.debug(self._property_range_message('h'))

    def check_temperature_bounds(self):
        r"""
        Check if temperature is within user specified limits.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to check fluid properties.
        """
        Tmin = max(
            [w._T_min for f, w in self.fluid.wrapper.items() if self.fluid.val[f] > ERR]
        ) * 1.01
        Tmax = min(
            [w._T_max for f, w in self.fluid.wrapper.items() if self.fluid.val[f] > ERR]
        ) * 0.99
        hmin = h_mix_pT(self.p.val_SI, Tmin, self.fluid_data, self.mixing_rule)
        hmax = h_mix_pT(self.p.val_SI, Tmax, self.fluid_data, self.mixing_rule)

        if self.h.val_SI < hmin:
            self.h.val_SI = hmin
            logger.debug(self._property_range_message('h'))

        if self.h.val_SI > hmax:
            self.h.val_SI = hmax
            logger.debug(self._property_range_message('h'))

    def _property_range_message(self, prop):
        r"""
        Return debugging message for fluid property range adjustments.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to check fluid properties.

        prop : str
            Fluid property.

        Returns
        -------
        msg : str
            Debugging message.
        """
        msg = (
            f"{fpd[prop]['text'][0].upper()} {fpd[prop]['text'][1:]} out of "
            f"fluid property range at connection {self.label}, adjusting value "
            f"to {self.get_attr(prop).val_SI} {fpd[prop]['SI_unit']}."
        )
        return msg

    def get_physical_exergy(self, p0, T0):
        r"""
        Get the value of a connection's specific physical exergy.

        Parameters
        ----------
        p0 : float
            Ambient pressure p0 / Pa.

        T0 : float
            Ambient temperature T0 / K.

        Note
        ----
            .. math::

                e^\mathrm{PH} = e^\mathrm{T} + e^\mathrm{M}\\
                E^\mathrm{T} = \dot{m} \cdot e^\mathrm{T}\\
                E^\mathrm{M} = \dot{m} \cdot e^\mathrm{M}\\
                E^\mathrm{PH} = \dot{m} \cdot e^\mathrm{PH}
        """
        self.ex_therm, self.ex_mech = fp.calc_physical_exergy(self, p0, T0)
        self.Ex_therm = self.ex_therm * self.m.val_SI
        self.Ex_mech = self.ex_mech * self.m.val_SI

        self.ex_physical = self.ex_therm + self.ex_mech
        self.Ex_physical = self.m.val_SI * self.ex_physical

    def get_chemical_exergy(self, p0, T0, Chem_Ex):
        r"""
        Get the value of a connection's specific chemical exergy.

        Parameters
        ----------
        p0 : float
            Ambient pressure p0 / Pa.

        T0 : float
            Ambient temperature T0 / K.

        Chem_Ex : dict
            Lookup table for standard specific chemical exergy.

        Note
        ----
            .. math::

                E^\mathrm{CH} = \dot{m} \cdot e^\mathrm{CH}
        """
        if Chem_Ex is None:
            self.ex_chemical = 0
        else:
            self.ex_chemical = fp.calc_chemical_exergy(self, p0, T0, Chem_Ex)

        self.Ex_chemical = self.m.val_SI * self.ex_chemical


class Ref:
    r"""
    A reference object is used to reference (unknown) properties of connections
    to other connections.

    For example, reference the mass flow of one connection :math:`\dot{m}` to
    another mass flow :math:`\dot{m}_{ref}`:

    .. math::

        \dot{m} = \dot{m}_\mathrm{ref} \cdot \mathrm{factor} + \mathrm{delta}

    Parameters
    ----------
    obj : tespy.connections.connection.Connection
        Connection to be referenced.

    factor : float
        Factor to multiply specified property with.

    delta : float
        Delta to add after multiplication.
    """

    def __init__(self, ref_obj, factor, delta):

        if not isinstance(ref_obj, Connection):
            msg = 'First parameter must be object of type connection.'
            logger.error(msg)
            raise TypeError(msg)

        if not (isinstance(factor, int) or isinstance(factor, float)):
            msg = 'Second parameter must be of type int or float.'
            logger.error(msg)
            raise TypeError(msg)

        if not (isinstance(delta, int) or isinstance(delta, float)):
            msg = 'Thrid parameter must be of type int or float.'
            logger.error(msg)
            raise TypeError(msg)

        self.obj = ref_obj
        self.factor = factor
        self.delta = delta
        self.delta_SI = None

        msg = ('Created reference object with factor ' + str(self.factor) +
               ' and delta ' + str(self.delta) + ' referring to connection ' +
               ref_obj.source.label + ' (' + ref_obj.source_id + ') -> ' +
               ref_obj.target.label + ' (' + ref_obj.target_id + ').')
        logger.debug(msg)

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
            logger.error(msg)
            raise KeyError(msg)
