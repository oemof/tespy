# -*- coding: utf-8
"""Module of class Connection and class Ref.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/connections/connection.py

SPDX-License-Identifier: MIT
"""

import warnings

import numpy as np
import pint

from tespy.components import Subsystem
from tespy.components.component import Component
from tespy.tools import fluid_properties as fp
from tespy.tools import logger
from tespy.tools.data_containers import DataContainer as dc
from tespy.tools.data_containers import FluidComposition as dc_flu
from tespy.tools.data_containers import FluidProperties as dc_prop
from tespy.tools.data_containers import ReferencedFluidProperties as dc_ref
from tespy.tools.data_containers import SimpleDataContainer as dc_simple
from tespy.tools.fluid_properties import CoolPropWrapper
from tespy.tools.fluid_properties import Q_mix_ph
from tespy.tools.fluid_properties import T_mix_ph
from tespy.tools.fluid_properties import T_sat_p
from tespy.tools.fluid_properties import dh_mix_dpQ
from tespy.tools.fluid_properties import dT_mix_dph
from tespy.tools.fluid_properties import dT_mix_pdh
from tespy.tools.fluid_properties import dv_mix_dph
from tespy.tools.fluid_properties import dv_mix_pdh
from tespy.tools.fluid_properties import h_mix_pQ
from tespy.tools.fluid_properties import h_mix_pT
from tespy.tools.fluid_properties import phase_mix_ph
from tespy.tools.fluid_properties import s_mix_ph
from tespy.tools.fluid_properties import v_mix_ph
from tespy.tools.fluid_properties import viscosity_mix_ph
from tespy.tools.fluid_properties.functions import T_bubble_p
from tespy.tools.fluid_properties.functions import T_dew_p
from tespy.tools.fluid_properties.functions import p_bubble_T
from tespy.tools.fluid_properties.functions import p_dew_T
from tespy.tools.fluid_properties.functions import p_sat_T
from tespy.tools.fluid_properties.helpers import get_mixture_temperature_range
from tespy.tools.fluid_properties.helpers import single_fluid
from tespy.tools.fluid_properties.wrappers import wrapper_registry
from tespy.tools.global_vars import ERR
from tespy.tools.global_vars import fluid_property_data as fpd
from tespy.tools.helpers import TESPyConnectionError
from tespy.tools.helpers import TESPyNetworkError
from tespy.tools.helpers import _get_dependents
from tespy.tools.helpers import _get_vector_dependents
from tespy.tools.helpers import _is_variable
from tespy.tools.helpers import _partial_derivative
from tespy.tools.helpers import _partial_derivative_vecvar
from tespy.tools.units import SI_UNITS


def connection_registry(type):
    connection_registry.items[type.__name__] = type
    return type


connection_registry.items = {}


class ConnectionBase:

    def __init__(self):
        pass

    def _remap_if_subsystem(self, source, target):
        # If the connected source or target is a subsystem we must
        # remap the source and target to its outlet/inlet
        if isinstance(source, Subsystem):
            source = source.outlet

        if isinstance(target, Subsystem):
            target = target.inlet

        return source, target

    def _check_types(self, source, target):
        # check input parameters
        if not (isinstance(source, Component) and
                isinstance(target, Component)):
            msg = (
                "Error creating connection. Check if source and target are "
                "tespy.components."
            )
            logger.error(msg)
            raise TypeError(msg)

    def _check_self_connect(self, source, target):
        if source == target:
            msg = (
                "Error creating connection. Cannot connect component "
                f"{source.label} to itself."
            )
            logger.error(msg)
            raise TESPyConnectionError(msg)

    def _check_connector_id(self, component, connector_id, connecter_locations):
        if connector_id not in connecter_locations:
            msg = (
                "Error creating connection. Specified connector for "
                f"{component.label} of class {component.__class__.__name__} "
                f"({connector_id})  is not available. Select one of the "
                f"following connectors {', '.join(connecter_locations)}."
            )
            logger.error(msg)
            raise ValueError(msg)

    def _parameter_specification(self, key, value):

        is_numeric = False
        is_quantity = False

        if isinstance(value, pint.Quantity):
            is_quantity = True
        else:
            try:
                float(value)
                is_numeric = True
            except (TypeError, ValueError):
                pass

        if key == "Td_bp":
            msg = (
                "The parameter 'Td_bp' is depreciated and will be removed in "
                "the next major release of tespy. Please use 'td_bubble' or "
                "'td_dew' instead. In contrast to 'Td_bp' not the following: "
                "A positive value for 'td_bubble' indicates a liquid state, "
                "meaning the temperature of the fluid will be lower than the "
                "associated bubble temperature by the specified value."
                "A positive value for 'td_dew' indicates a gaseous state, "
                "meaning the temperature of the fluid will be higher than the "
                "associated dew temperature by the specified value."
            )
            warnings.warn(msg, FutureWarning)

        if value is None:
            self.get_attr(key).set_attr(is_set=False)

            if f"{key}_ref" in self.property_data:
                self.get_attr(f"{key}_ref").set_attr(is_set=False)

        elif is_numeric or is_quantity:
            # value specification
            if key in self.property_data:
                self.get_attr(key).set_attr(is_set=True, val=value)
            else:
                self.get_attr(key.replace('0', '')).set_attr(val0=value)

        # reference object
        elif isinstance(value, Ref):
            if f"{key}_ref" not in self.property_data:
                msg = f"Referencing {key} is not implemented."
                logger.error(msg)
                raise NotImplementedError(msg)
            else:
                self.get_attr(f"{key}_ref").set_attr(ref=value)
                self.get_attr(f"{key}_ref").set_attr(is_set=True)

        # invalid datatype for keyword
        else:
            msg = f"Wrong datatype for keyword argument {key}."
            logger.error(msg)
            raise TypeError(msg)

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

    def _serialize(self):
        export = {}
        export.update({"source": self.source.label})
        export.update({"target": self.target.label})
        for k in self._serializable():
            export.update({k: self.get_attr(k)})
        for k in self.property_data:
            data = self.get_attr(k)
            export.update({k: data._serialize()})

        return {self.label: export}

    @staticmethod
    def _serializable():
        return [
            "source_id", "target_id",
            "design_path", "design", "offdesign",
            "local_design", "local_design",
            "printout"
        ]

    def get_variables(self):
        return {}

    def _preprocess(self, row_idx):
        self.num_eq = 0

        self._structure_matrix = {}
        self._rhs = {}
        self._equation_set_lookup = {}

        for parameter in self.parameters:
            container = self.get_attr(parameter)
            if container.is_set and container.func is not None:
                num_eq = self.parameters[parameter].num_eq
                # the row index matches the location in the network's rhs
                # and matrix
                for i in range(self.num_eq, self.num_eq + num_eq):
                    self._equation_set_lookup[i + row_idx] = parameter
                    self._rhs[i + row_idx] = 0
                # the structure matrix function also computes the rhs
                if container.structure_matrix is not None:
                    container.structure_matrix(
                        row_idx + self.num_eq, **container.func_params
                    )

                self.num_eq += num_eq

    def _prepare_for_solver(self, system_dependencies, eq_counter):
        self.num_eq = 0
        self.it = 0
        self.equations = {}
        self._equation_lookup = {}
        self._equation_scalar_dependents_lookup = {}
        self._equation_vector_dependents_lookup = {}

        for eq_num, value in self._equation_set_lookup.items():
            if eq_num in system_dependencies:
                continue

            if value not in self.equations:
                data = self.parameters[value]
                self.equations.update({value: data})
                self._assign_dependents_and_eq_mapping(
                    value, data, self.equations, eq_counter
                )
                self.num_eq += data.num_eq
                eq_counter += data.num_eq

        self.residual = {}
        self.jacobian = {}

        return eq_counter

    def _assign_dependents_and_eq_mapping(self, value, data, eq_dict, eq_counter):
        if data.dependents is None:
            scalar_dependents = [[] for _ in range(data.num_eq)]
            vector_dependents = [{} for _ in range(data.num_eq)]
        else:
            dependents = data.dependents(**data.func_params)
            if type(dependents) == list:
                scalar_dependents = _get_dependents(dependents)
                vector_dependents = [{} for _ in range(data.num_eq)]
            else:
                scalar_dependents = _get_dependents(dependents["scalars"])
                vector_dependents = _get_vector_dependents(dependents["vectors"])

                # this is a temporary fix
                if len(vector_dependents) < data.num_eq:
                    vector_dependents = [{} for _ in range(data.num_eq)]

        eq_dict[value]._scalar_dependents = scalar_dependents
        eq_dict[value]._vector_dependents = vector_dependents
        eq_dict[value]._first_eq_index = eq_counter

        for i in range(data.num_eq):
            self._equation_lookup[eq_counter + i] = (value, i)
            self._equation_scalar_dependents_lookup[eq_counter + i] = scalar_dependents[i]
            self._equation_vector_dependents_lookup[eq_counter + i] = vector_dependents[i]

    def _partial_derivative(self, var, eq_num, value, increment_filter=None, **kwargs):
        result = _partial_derivative(var, value, increment_filter, **kwargs)
        if result is not None:
            self.jacobian[eq_num, var.J_col] = result

    def _adjust_to_property_limits(self, nw):
        pass

    @classmethod
    def _print_attributes(cls):
        return []

    @classmethod
    def _result_attributes(cls):
        return []

    @classmethod
    def _get_result_cols(cls, all_fluids):
        return []

    def calc_results(self):
        return True

    def collect_results(self, all_fluids):
        return None


@connection_registry
class Connection(ConnectionBase):
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

    >>> so_si1.set_attr(v=0.012, m0=10, p=5, h=400, fluid={'H2O': 1})
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
    >>> so_si1.m.is_set
    False
    >>> so_si1.m.get_attr('is_set')
    False
    >>> type(so_si2.m_ref.ref)
    <class 'tespy.connections.connection.Ref'>
    >>> so_si2.fluid_balance.is_set
    True
    >>> so_si2.m_ref.ref.get_attr('delta')
    -5
    >>> so_si2.m_ref.is_set
    True
    >>> type(so_si2.m_ref.ref.get_attr('obj'))
    <class 'tespy.connections.connection.Connection'>

    Unset the specified temperature and specify temperature difference to
    boiling point (deprecated) instead.

    >>> so_si2.T.is_set
    True
    >>> so_si2.set_attr(Td_bp=5, T=None)
    >>> so_si2.T.is_set
    False
    >>> so_si2.Td_bp.val
    5.0
    >>> so_si2.set_attr(Td_bp=None)
    >>> so_si2.Td_bp.is_set
    False

    Bubble line or dew line temperature difference:

    >>> so_si2.set_attr(td_bubble=5)
    >>> so_si2.td_bubble.is_set
    True
    >>> so_si2.td_bubble.val
    5.0
    >>> so_si2.set_attr(td_bubble=None)
    >>> so_si2.td_bubble.is_set
    False
    >>> so_si2.set_attr(td_dew=5)
    >>> so_si2.td_dew.is_set
    True
    >>> so_si2.td_dew.val
    5.0
    >>> so_si2.set_attr(td_dew=None)
    >>> so_si2.td_dew.is_set
    False

    Specify the state keyword: The fluid will be forced to liquid or gaseous
    state in this case.

    >>> so_si2.set_attr(state='l')
    >>> so_si2.state.is_set
    True
    >>> so_si2.set_attr(state=None)
    >>> so_si2.state.is_set
    False
    >>> so_si2.set_attr(state='g')
    >>> so_si2.state.is_set
    True
    >>> so_si2.set_attr(state=None)
    >>> so_si2.state.is_set
    False
    """

    def __init__(self, source, outlet_id, target, inlet_id,
                 label=None, **kwargs):

        source, target = self._remap_if_subsystem(source, target)
        self._check_types(source, target)
        self._check_self_connect(source, target)
        self._check_connector_id(source, outlet_id, source.outlets())
        self._check_connector_id(target, inlet_id, target.inlets())

        self.label = f"{source.label}:{outlet_id}_{target.label}:{inlet_id}"
        if label is not None:
            self.label = label
            if not isinstance(label, str):
                msg = "Please provide the label as string."
                logger.error(msg)
                raise TypeError(msg)

        # set specified values
        self.source = source
        self.source_id = outlet_id
        self.target = target
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
            if hasattr(v, "func") and v.func is not None
        }
        self.state = dc_simple()
        self.phase = dc_simple()
        self.property_data0 = [x + '0' for x in self.property_data.keys()]
        self.__dict__.update(self.property_data)
        self.mixing_rule = None
        msg = (
            f"Created connection from {self.source.label} ({self.source_id}) "
            f"to {self.target.label} ({self.target_id})."
        )
        logger.debug(msg)

        self.set_attr(**kwargs)

    def _reset_design(self, redesign):
        for value in self.get_variables().values():
            value.design = np.nan

        self.fluid.design = {}

        self.new_design = True

        # switch connections to design mode
        if redesign:
            for var in self.design:
                self.get_attr(var).is_set = True

            for var in self.offdesign:
                self.get_attr(var).is_set = False

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
            Fluid composition specification.

        fluid0 : dict
            Starting value specification for fluid composition.

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
          on the specified connection to 100 %. For example, you have four
          fluid components (a, b, c and d) in your vector, you set two of them
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
                msg = 'Label can only be specified on instance creation.'
                logger.error(msg)
                raise TESPyConnectionError(msg)
            elif 'fluid' in key:
                self._fluid_specification(key, kwargs[key])

            elif key in self.property_data or key in self.property_data0:
                self._parameter_specification(key, kwargs[key])

            elif key == 'state':
                if kwargs[key] in ['l', 'g']:
                    self.state.set_attr(_val=kwargs[key], is_set=True)
                elif kwargs[key] is None:
                    self.state.set_attr(is_set=False)
                else:
                    msg = (
                        'Keyword argument "state" must either be '
                        '"l" or "g" or be None.'
                    )
                    logger.error(msg)
                    raise TypeError(msg)

            # design/offdesign parameter list
            elif key in ['design', 'offdesign']:
                if not isinstance(kwargs[key], list):
                    msg = f"Please provide the {key} parameters as list!"
                    logger.error(msg)
                    raise TypeError(msg)
                elif set(kwargs[key]).issubset(self.property_data.keys()):
                    self.__dict__.update({key: kwargs[key]})
                else:
                    params = ', '.join(self.property_data.keys())
                    msg = (
                        "Available parameters for (off-)design specification "
                        f"are: {params}."
                    )
                    logger.error(msg)
                    raise ValueError(msg)

            # design path
            elif key == 'design_path':
                self.__dict__.update({key: kwargs[key]})
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
                msg = f"Connection has no attribute {key}."
                logger.error(msg)
                raise KeyError(msg)

    def _fluid_specification(self, key, value):

        self._check_fluid_datatypes(key, value)

        if key == "fluid":
            for fluid, fraction in value.items():
                if "::" in fluid:
                    back_end, fluid = fluid.split("::")
                else:
                    back_end = None

                if fraction is None:
                    if fluid in self.fluid.is_set:
                        self.fluid.is_set.remove(fluid)

                else:
                    self.fluid.val[fluid] = fraction
                    self.fluid.is_set.add(fluid)
                    self.fluid.back_end[fluid] = back_end

        elif key == "fluid0":
            self.fluid.val0.update(value)

        elif key == "fluid_engines":
            self.fluid.engine = value

        elif key == "fluid_balance":
            self.fluid_balance.is_set = value

        else:
            msg = f"Connections do not have an attribute named {key}"
            logger.error(msg)
            raise KeyError(msg)

    def _check_fluid_datatypes(self, key, value):
        if key == "fluid_balance":
            if not isinstance(value, bool):
                msg = "Datatype for 'fluid_balance' must be boolean."
                logger.error(msg)
                raise TypeError(msg)
        else:
            if not isinstance(value, dict):
                msg = "Datatype for fluid vector specification must be dict."
                logger.error(msg)
                raise TypeError(msg)

    def _serialize(self):
        export = super()._serialize()
        export[self.label].update({"state": self.state._serialize()})
        return export

    def _deserialize(self, data, all_connections):
        arglist = [
            _ for _ in data
            if _ not in ["source", "source_id", "target", "target_id", "label", "fluid"]
            and "ref" not in _
        ]

        for arg in arglist:
            container = self.get_attr(arg)
            if isinstance(container, dc):
                container.set_attr(**data[arg])
            else:
                self.set_attr(**{arg: data[arg]})

        for f, engine in data["fluid"]["engine"].items():
            data["fluid"]["engine"][f] = wrapper_registry.items[engine]

        self.fluid.set_attr(**data["fluid"])
        self._create_fluid_wrapper()

        arglist_ref = [_ for _ in data if "ref" in _]

        for arg in arglist_ref:
            if len(data[arg]) > 0:
                param = arg.replace("_ref", "")
                ref = Ref(
                    all_connections[data[arg]["conn"]],
                    data[arg]["factor"],
                    data[arg]["delta"]
                )
                self.set_attr(**{param: ref})

    def _serializable(self):
        return super()._serializable() + ["mixing_rule"]

    def _create_fluid_wrapper(self):
        for fluid in self.fluid.val:
            if fluid in self.fluid.wrapper:
                continue
            if fluid not in self.fluid.engine:
                self.fluid.engine[fluid] = CoolPropWrapper

            back_end = None
            if fluid in self.fluid.back_end:
                back_end = self.fluid.back_end[fluid]
            else:
                self.fluid.back_end[fluid] = None

            self.fluid.wrapper[fluid] = self.fluid.engine[fluid](fluid, back_end)

    def _precalc_guess_values(self):
        """
        Precalculate enthalpy values for connections.

        Precalculation is performed only if temperature, vapor mass fraction,
        temperature difference to boiling point or phase is specified.

        Parameters
        ----------
        c : tespy.connections.connection.Connection
            Connection to precalculate values for.
        """
        # starting values for specified vapour content or temperature
        if not self.h.is_var:
            return

        if not self.good_starting_values:
            if self.x.is_set:
                fluid = fp.single_fluid(self.fluid_data)
                if self.p.is_var and self.p.val_SI > self.fluid.wrapper[fluid]._p_crit:
                    self.p.set_reference_val_SI(self.fluid.wrapper[fluid]._p_crit * 0.9)
                self.h.set_reference_val_SI(
                    fp.h_mix_pQ(self.p.val_SI, self.x.val_SI, self.fluid_data, self.mixing_rule)
                )
            if self.T.is_set:
                try:
                    self.h.set_reference_val_SI(
                        fp.h_mix_pT(self.p.val_SI, self.T.val_SI, self.fluid_data, self.mixing_rule)
                    )
                except ValueError:
                    pass

        # starting values for specified quality, specified subcooling/overheating
        # and state specification. These should be recalculated even with
        # good starting values, for example, when one exchanges enthalpy
        # with boiling point temperature difference.
        if (self.Td_bp.is_set or self.state.is_set or self.td_dew.is_set or self.td_bubble.is_set):
            fluid = fp.single_fluid(self.fluid_data)
            if self.p.is_var and self.p.val_SI > self.fluid.wrapper[fluid]._p_crit:
                self.p.set_reference_val_SI(self.fluid.wrapper[fluid]._p_crit * 0.9)
            if (
                    (self.Td_bp.val_SI > 0 and self.Td_bp.is_set)
                    or (self.state.val == 'g' and self.state.is_set)
                    or (self.td_dew.val_SI >= 0 and self.td_dew.is_set)
                    or (self.td_bubble.val_SI < 0 and self.td_bubble.is_set)
                ):
                h = fp.h_mix_pQ(self.p.val_SI, 1, self.fluid_data)
                if self.h.val_SI < h:
                    self.h.set_reference_val_SI(h + 1e3)

            elif (
                    (self.Td_bp.val_SI < 0 and self.Td_bp.is_set)
                    or (self.state.val == 'l' and self.state.is_set)
                    or (self.td_bubble.val_SI >= 0 and self.td_bubble.is_set)
                    or (self.td_dew.val_SI < 0 and self.td_dew.is_set)
                ):
                h = fp.h_mix_pQ(self.p.val_SI, 0, self.fluid_data)
                if self.h.val_SI > h:
                    self.h.set_reference_val_SI(h - 1e3)

    def _presolve(self):
        if len(self.fluid.is_var) > 0:
            return []

        specifications = []
        for name, container in self.property_data.items():
            if name in ["p", "h", "T", "x", "Td_bp", "td_bubble", "td_dew"]:
                if container.is_set:
                    specifications += [name]

        num_specs = len(specifications)

        if num_specs > 2:
            msg = (
                "You have specified more than 2 parameters for the connection "
                f"{self.label} with a known fluid compoistion: "
                f"{', '.join(specifications)}. This overdetermines the state "
                "of the fluid."
            )
            raise TESPyNetworkError(msg)

        presolved_equations = []
        if self.h.is_var and not self.p.is_var:
            if self.T.is_set:
                self.h.set_reference_val_SI(h_mix_pT(self.p.val_SI, self.T.val_SI, self.fluid_data, self.mixing_rule))
                self.h._potential_var = False
                if "T" in self._equation_set_lookup.values():
                    presolved_equations += ["T"]
                msg = f"Determined h by known p and T at {self.label}."
                logger.info(msg)

            elif self.Td_bp.is_set:
                T_sat = T_sat_p(self.p.val_SI, self.fluid_data)
                self.h.set_reference_val_SI(h_mix_pT(self.p.val_SI, T_sat + self.Td_bp.val_SI, self.fluid_data))
                self.h._potential_var = False
                if "Td_bp" in self._equation_set_lookup.values():
                    presolved_equations += ["Td_bp"]
                msg = f"Determined h by known p and Td_bp at {self.label}."
                logger.info(msg)

            elif self.td_bubble.is_set:
                T_bubble = T_bubble_p(self.p.val_SI, self.fluid_data)
                # fix for pure fluids: T cannot be too close to saturation
                if abs(self.td_bubble.val_SI) <= 1e-3:
                    # at saturation we can use h_mix_pQ
                    self.h.set_reference_val_SI(h_mix_pQ(self.p.val_SI, 0, self.fluid_data))
                else:
                    self.h.set_reference_val_SI(h_mix_pT(self.p.val_SI, T_bubble - self.td_bubble.val_SI, self.fluid_data))
                self.h._potential_var = False
                if "td_bubble" in self._equation_set_lookup.values():
                    presolved_equations += ["td_bubble"]
                msg = f"Determined h by known p and td_bubble at {self.label}."
                logger.info(msg)

            elif self.td_dew.is_set:
                T_dew = T_dew_p(self.p.val_SI, self.fluid_data)
                # fix for pure fluids: T cannot be too close to saturation
                if abs(self.td_dew.val_SI) <= 1e-3:
                    # at saturation we can use h_mix_pQ
                    self.h.set_reference_val_SI(h_mix_pQ(self.p.val_SI, 1, self.fluid_data))
                else:
                    self.h.set_reference_val_SI(h_mix_pT(self.p.val_SI, T_dew + self.td_dew.val_SI, self.fluid_data))
                self.h._potential_var = False
                if "td_dew" in self._equation_set_lookup.values():
                    presolved_equations += ["td_dew"]
                msg = f"Determined h by known p and td_bubble at {self.label}."
                logger.info(msg)

            elif self.x.is_set:
                self.h.set_reference_val_SI(h_mix_pQ(self.p.val_SI, self.x.val_SI, self.fluid_data))
                self.h._potential_var = False
                if "x" in self._equation_set_lookup.values():
                    presolved_equations += ["x"]
                msg = f"Determined h by known p and x at {self.label}."
                logger.info(msg)

        elif self.h.is_var and self.p.is_var:
            if self.T.is_set and self.x.is_set:
                self.p.set_reference_val_SI(p_sat_T(self.T.val_SI, self.fluid_data))
                self.p._potential_var = False
                self.h.set_reference_val_SI(h_mix_pQ(self.p.val_SI, self.x.val_SI, self.fluid_data))
                self.h._potential_var = False
                if "T" in self._equation_set_lookup.values():
                    presolved_equations += ["T"]
                if "x" in self._equation_set_lookup.values():
                    presolved_equations += ["x"]
                msg = f"Determined h and p by known T and x at {self.label}."
                logger.info(msg)

            elif self.T.is_set and self.Td_bp.is_set:
                self.p.set_reference_val_SI(p_sat_T(self.T.val_SI - self.Td_bp.val_SI, self.fluid_data))
                self.p._potential_var = False
                self.h.set_reference_val_SI(h_mix_pT(self.p.val_SI, self.T.val_SI, self.fluid_data))
                self.h._potential_var = False
                if "T" in self._equation_set_lookup.values():
                    presolved_equations += ["T"]
                if "Td_bp" in self._equation_set_lookup.values():
                    presolved_equations += ["Td_bp"]
                msg = f"Determined h and p by known T and Td_bp at {self.label}."
                logger.info(msg)

            elif self.T.is_set and self.td_bubble.is_set:
                self.p.set_reference_val_SI(p_bubble_T(self.T.val_SI + self.td_bubble.val_SI, self.fluid_data))
                self.p._potential_var = False
                if round(self.td_bubble.val_SI, 6) == 0:
                    self.h.set_reference_val_SI(h_mix_pQ(self.p.val_SI, 0, self.fluid_data))
                else:
                    self.h.set_reference_val_SI(h_mix_pT(self.p.val_SI, self.T.val_SI, self.fluid_data))
                self.h._potential_var = False
                if "T" in self._equation_set_lookup.values():
                    presolved_equations += ["T"]
                if "td_bubble" in self._equation_set_lookup.values():
                    presolved_equations += ["td_bubble"]
                msg = f"Determined h and p by known T and td_bubble at {self.label}."
                logger.info(msg)

            elif self.T.is_set and self.td_dew.is_set:
                self.p.set_reference_val_SI(p_dew_T(self.T.val_SI - self.td_dew.val_SI, self.fluid_data))
                self.p._potential_var = False
                if round(self.td_dew.val_SI, 6) == 0:
                    self.h.set_reference_val_SI(h_mix_pQ(self.p.val_SI, 1, self.fluid_data))
                else:
                    self.h.set_reference_val_SI(h_mix_pT(self.p.val_SI, self.T.val_SI, self.fluid_data))
                self.h._potential_var = False
                if "T" in self._equation_set_lookup.values():
                    presolved_equations += ["T"]
                if "td_dew" in self._equation_set_lookup.values():
                    presolved_equations += ["td_dew"]
                msg = f"Determined h and p by known T and td_dew at {self.label}."
                logger.info(msg)

        presolved_equations = [
            key for parameter in presolved_equations
            for key, value in self._equation_set_lookup.items()
            if value == parameter
        ]
        return presolved_equations

    def _partial_derivative_fluid(self, var, eq_num, value, dx, increment_filter=None, **kwargs):
        result = _partial_derivative_vecvar(var, value, dx, increment_filter, **kwargs)
        if result is not None:
            self.jacobian[eq_num, var.J_col[dx]] = result

    def reset_fluid_vector(self):
        self.fluid = dc_flu()

    def get_variables(self):
        return {"m": self.m, "p": self.p, "h": self.h}

    def get_parameters(self):
        return {
            "m": dc_prop(d=1e-4, quantity="mass_flow"),
            "p": dc_prop(d=1e-3, quantity="pressure"),
            "h": dc_prop(d=1e-3, quantity="enthalpy"),
            "vol": dc_prop(quantity="specific_volume"),
            "s": dc_prop(quantity="entropy"),
            "fluid": dc_flu(d=1e-5),
            "fluid_balance": dc_simple(
                func=self.fluid_balance_func,
                deriv=self.fluid_balance_deriv,
                _val=False, num_eq_sets=1,
                dependents=self.fluid_balance_dependents
            ),
            "T": dc_prop(
                func=self.T_func, deriv=self.T_deriv,
                dependents=self.T_dependents, num_eq=1,
                quantity="temperature"
            ),
            "v": dc_prop(
                func=self.v_func, deriv=self.v_deriv,
                dependents=self.v_dependents, num_eq=1,
                quantity="volumetric_flow"
            ),
            "x": dc_prop(
                func=self.x_func, deriv=self.x_deriv,
                dependents=self.x_dependents, num_eq=1,
                quantity="quality"
            ),
            "Td_bp": dc_prop(
                func=self.Td_bp_func, deriv=self.Td_bp_deriv,
                dependents=self.Td_bp_dependents, num_eq=1,
                quantity="temperature_difference"
            ),
            "td_dew": dc_prop(
                func=self.td_dew_func,
                dependents=self.td_dew_dependents, num_eq=1,
                quantity="temperature_difference"
            ),
            "td_bubble": dc_prop(
                func=self.td_bubble_func, #deriv=self.td_bubble_deriv,
                dependents=self.td_bubble_dependents, num_eq=1,
                quantity="temperature_difference"
            ),
            "m_ref": dc_ref(
                func=self.primary_ref_func,
                num_eq=1, func_params={"variable": "m"},
                structure_matrix=self.primary_ref_structure_matrix,
                quantity="mass_flow"
            ),
            "p_ref": dc_ref(
                func=self.primary_ref_func,
                num_eq=1, func_params={"variable": "p"},
                structure_matrix=self.primary_ref_structure_matrix,
                quantity="pressure"
            ),
            "h_ref": dc_ref(
                func=self.primary_ref_func,
                num_eq=1, func_params={"variable": "h"},
                structure_matrix=self.primary_ref_structure_matrix,
                quantity="enthalpy"
            ),
            "T_ref": dc_ref(
                func=self.T_ref_func, deriv=self.T_ref_deriv,
                dependents=self.T_ref_dependents, num_eq=1,
                quantity="temperature_difference"  # reference has delta T
            ),
            "v_ref": dc_ref(
                func=self.v_ref_func, deriv=self.v_ref_deriv,
                dependents=self.v_ref_dependents, num_eq=1,
                quantity="volumetric_flow"
            ),

        }

    def get_fluid_data(self):
        return {
            fluid: {
                "wrapper": self.fluid.wrapper[fluid],
                "mass_fraction": self.fluid.val[fluid]
            } for fluid in self.fluid.val
        }

    fluid_data = property(get_fluid_data)

    def primary_ref_func(self, **kwargs):
        variable = kwargs["variable"]
        self.get_attr(variable)
        ref = self.get_attr(f"{variable}_ref").ref
        return (
            self.get_attr(variable).val_SI
            - (ref.obj.get_attr(variable).val_SI * ref.factor + ref.delta_SI)
        )

    def primary_ref_structure_matrix(self, k, **kwargs):
        variable = kwargs["variable"]
        ref = self.get_attr(f"{variable}_ref").ref
        self._structure_matrix[k, self.get_attr(variable).sm_col] = 1
        self._structure_matrix[k, ref.obj.get_attr(variable).sm_col] = -ref.factor
        self._rhs[k] = ref.delta_SI

    def calc_T(self, T0=None):
        if T0 is None:
            T0 = self.T.val_SI
        return T_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data, self.mixing_rule, T0=T0)

    def T_func(self, **kwargs):
        return self.calc_T() - self.T.val_SI

    def T_deriv(self, increment_filter, k, **kwargs):
        if _is_variable(self.p):
            self.jacobian[k, self.p.J_col] = (
                dT_mix_dph(self.p.val_SI, self.h.val_SI, self.fluid_data, self.mixing_rule, self.T.val_SI)
            )
        if _is_variable(self.h):
            self.jacobian[k, self.h.J_col] = (
                dT_mix_pdh(self.p.val_SI, self.h.val_SI, self.fluid_data, self.mixing_rule, self.T.val_SI)
            )

    def T_dependents(self):
        return [self.p, self.h]

    def T_ref_func(self, **kwargs):
        ref = self.T_ref.ref
        return self.calc_T() - (ref.obj.calc_T() * ref.factor + ref.delta_SI)

    def T_ref_deriv(self, increment_filter, k, **kwargs):
        # first part of sum is identical to direct temperature specification
        self.T_deriv(increment_filter, k, **kwargs)
        ref = self.T_ref.ref
        if _is_variable(ref.obj.p):
            self.jacobian[k, ref.obj.p.J_col] = -(
                dT_mix_dph(ref.obj.p.val_SI, ref.obj.h.val_SI, ref.obj.fluid_data, ref.obj.mixing_rule)
            ) * ref.factor
        if _is_variable(ref.obj.h):
            self.jacobian[k, ref.obj.h.J_col] = -(
                dT_mix_pdh(ref.obj.p.val_SI, ref.obj.h.val_SI, ref.obj.fluid_data, ref.obj.mixing_rule)
            ) * ref.factor

    def T_ref_dependents(self):
        ref = self.T_ref.ref
        return self.T_dependents() + ref.obj.T_dependents()

    def calc_viscosity(self, T0=None):
        try:
            return viscosity_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data, self.mixing_rule, T0=T0)
        except NotImplementedError:
            return np.nan

    def calc_vol(self, T0=None):
        try:
            return v_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data, self.mixing_rule, T0=T0)
        except NotImplementedError:
            return np.nan

    def v_func(self, **kwargs):
        return self.calc_vol(T0=self.T.val_SI) * self.m.val_SI - self.v.val_SI

    def v_deriv(self, increment_filter, k, **kwargs):
        if _is_variable(self.m):
            self._partial_derivative(self.m, k, self.calc_vol(T0=self.T.val_SI))
        if _is_variable(self.p):
            self._partial_derivative(
                self.p, k,
                dv_mix_dph(self.p.val_SI, self.h.val_SI, self.fluid_data)
                * self.m.val_SI
            )
        if _is_variable(self.h):
            self._partial_derivative(
                self.h, k,
                dv_mix_pdh(self.p.val_SI, self.h.val_SI, self.fluid_data)
                * self.m.val_SI
            )

    def v_dependents(self):
        return [self.m, self.p, self.h]

    def v_ref_func(self, **kwargs):
        ref = self.v_ref.ref
        return (
            self.calc_vol(T0=self.T.val_SI) * self.m.val_SI
            - (
                ref.obj.calc_vol(T0=ref.obj.T.val_SI) * ref.obj.m.val_SI
                * ref.factor + ref.delta_SI
            )
        )

    def v_ref_deriv(self, increment_filter, k, **kwargs):
        # first part of sum is identical to direct volumetric flow specification
        self.v_deriv(increment_filter, k, **kwargs)

        ref = self.v_ref.ref
        if ref.obj.m.is_var:
            self.jacobian[k, ref.obj.m.J_col] = -(
                ref.obj.calc_vol(T0=ref.obj.T.val_SI) * ref.factor
            )
        if ref.obj.p.is_var:
            self.jacobian[k, ref.obj.p.J_col] = -(
                dv_mix_dph(ref.obj.p.val_SI, ref.obj.h.val_SI, ref.obj.fluid_data)
                * ref.obj.m.val_SI * ref.factor
            )
        if ref.obj.h.is_var:
            self.jacobian[k, ref.obj.h.J_col] = -(
                dv_mix_pdh(ref.obj.p.val_SI, ref.obj.h.val_SI, ref.obj.fluid_data)
                * ref.obj.m.val_SI * ref.factor
            )

    def v_ref_dependents(self):
        ref = self.v_ref.ref
        return self.v_dependents() + ref.obj.v_dependents()

    def calc_x(self):
        try:
            return Q_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data)
        except NotImplementedError:
            return np.nan

    def x_func(self, **kwargs):
        # saturated steam fraction
        return (
            self.h.val_SI
            - h_mix_pQ(self.p.val_SI, self.x.val_SI, self.fluid_data)
        )

    def x_deriv(self, increment_filter, k, **kwargs):
        if self.p.is_var:
            self.jacobian[k, self.p.J_col] = -dh_mix_dpQ(self.p.val_SI, self.x.val_SI, self.fluid_data)
        if self.h.is_var:
            self.jacobian[k, self.h.J_col] = 1

    def x_dependents(self):
        return [self.p, self.h]

    def calc_T_sat(self):
        try:
            return T_sat_p(self.p.val_SI, self.fluid_data)
        except NotImplementedError:
            return np.nan

    def calc_Td_bp(self):
        try:
            return self.calc_T() - T_sat_p(self.p.val_SI, self.fluid_data)
        except NotImplementedError:
            return np.nan

    def calc_td_dew(self):
        try:
            return self.calc_T() - T_dew_p(self.p.val_SI, self.fluid_data)
        except NotImplementedError:
            return np.nan

    def calc_td_bubble(self):
        try:
            return T_bubble_p(self.p.val_SI, self.fluid_data) - self.calc_T()
        except NotImplementedError:
            return np.nan

    def Td_bp_func(self, **kwargs):
        # temperature difference to boiling point
        return self.calc_Td_bp() - self.Td_bp.val_SI

    def Td_bp_deriv(self, increment_filter, k, **kwargs):
        f = self.Td_bp_func
        self._partial_derivative(self.p, k, f)
        self._partial_derivative(self.h, k, f)

    def Td_bp_dependents(self):
        return [self.p, self.h]

    def td_dew_func(self, **kwargs):
        return self.calc_td_dew() - self.td_dew.val_SI

    def td_dew_dependents(self):
        return [self.p, self.h]

    def td_bubble_func(self, **kwargs):
        return self.calc_td_bubble() - self.td_bubble.val_SI

    def td_bubble_dependents(self):
        return [self.p, self.h]

    def fluid_balance_func(self, **kwargs):
        residual = 1 - sum(self.fluid.val[f] for f in self.fluid.is_set)
        residual -= sum(self.fluid.val[f] for f in self.fluid.is_var)
        return residual

    def fluid_balance_deriv(self, increment_filter, k, **kwargs):
        for f in self.fluid.is_var:
            self.jacobian[k, self.fluid.J_col[f]] = -self.fluid.val[f]

    def fluid_balance_dependents(self):
        return {
            "scalars": [[]],
            "vectors": [{self.fluid: self.fluid.is_var}]
        }

    def calc_s(self):
        try:
            return s_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data, self.mixing_rule, T0=self.T.val_SI)
        except NotImplementedError:
            return np.nan

    def calc_Q(self):
        return Q_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data)

    def calc_phase(self):
        try:
            return phase_mix_ph(self.p.val_SI, self.h.val_SI, self.fluid_data)
        except NotImplementedError:
            return np.nan

    def calc_results(self, units):
        self.T.val_SI = self.calc_T()
        fluid = single_fluid(self.fluid_data)
        _converged = True
        if fluid is None:
            # this is a mixture
            h_from_T = h_mix_pT(self.p.val_SI, self.T.val_SI, self.fluid_data, self.mixing_rule)
            if (
                abs(h_from_T - self.h.val_SI) > ERR ** .5 and
                abs((h_from_T - self.h.val_SI) / self.h.val_SI) > ERR ** .5
            ):
                self.T.val_SI = np.nan
                self.vol.val_SI = np.nan
                self.v.val_SI = np.nan
                self.s.val_SI = np.nan
                msg = (
                    "Could not find a feasible value for mixture temperature "
                    f"at connection {self.label}. The values of temperature, "
                    "specific volume and entropy are set to nan."
                )
                logger.error(msg)
                _converged = False
            else:
                _, Tmax = get_mixture_temperature_range(self.fluid_data)
                if self.T.val_SI > Tmax:
                    msg = (
                        "The temperature value of the mixture is above the "
                        "upper temperature limit of a mixture component. The "
                        "resulting temperature may have larger deviations "
                        "compared to the tolerance specified in the "
                        "corresponding substance property library."
                    )
                    logger.warning(msg)
        else:
            # these are pure fluids
            # two-phase properties are calculated based on pressure
            if self.p.val_SI < self.fluid.wrapper[fluid]._p_crit:
                if not self.x.is_set:
                    try:
                        self.x.val_SI = self.calc_x()
                    except ValueError:
                        self.x.val_SI = np.nan
                if not self.Td_bp.is_set:
                    try:
                        self.Td_bp.val_SI = self.calc_Td_bp()
                    except ValueError:
                        self.Td_bp.val_SI = np.nan
                try:
                    self.phase.val = self.calc_phase()
                except ValueError:
                    self.phase.val = "phase not recognized"
            else:
                self.x.val_SI = np.nan
                self.Td_bp.val_SI = np.nan
                self.phase.val = "phase not recognized"

        if _converged:
            self.vol.val_SI = self.calc_vol()
            self.v.val_SI = self.vol.val_SI * self.m.val_SI
            self.s.val_SI = self.calc_s()

        for prop in self._result_attributes():
            param = self.get_attr(prop)
            result = param._get_val_from_SI(units)
            if param.is_set and round(result.magnitude, 3) != round(param.val, 3):
                _converged = False
                msg = (
                    "The simulation converged but the calculated result "
                    f"{result} for the fixed input parameter {prop} of "
                    f"connection {self.label} is not equal to the originally "
                    f"specified value of {param.val}. Usually, this can "
                    "happen, when a method internally manipulates the "
                    "associated equation during iteration in order to allow "
                    "progress in situations, when the equation is otherwise "
                    "not well defined for the current values of the "
                    "variables, e.g. in case a negative root would need to be "
                    "evaluated. Often, this can happen during the first "
                    "iterations and then will resolve itself as convergence "
                    "progresses. In this case it did not, meaning convergence "
                    "was not actually achieved."
                )
                logger.warning(msg)
            else:
                if not param.is_set:
                    param.set_val_from_SI(units)

        self.m.set_val0_from_SI(units)
        self.p.set_val0_from_SI(units)
        self.h.set_val0_from_SI(units)
        self.fluid.val0 = self.fluid.val.copy()

        return _converged

    def _set_design_params(self, data, units):
        for var in self._result_attributes():
            unit = data[f"{var}_unit"]
            if unit == "C":
                if var == "T":
                    unit = "degC"
                elif var == "Td_bp":
                    unit = "delta_degC"
            elif "kgK" in unit:
                unit = unit.replace("kgK", "kg/K")
            elif unit == "-":
                unit = "1"
            param = self.get_attr(var)
            param.design = units.ureg.Quantity(
                float(data[var]),
                unit
            ).to(SI_UNITS[param.quantity]).magnitude

        for fluid in self.fluid.val:
            self.fluid.design[fluid] = float(data[fluid])

    def _set_starting_values(self, data, units):
        for prop in self.get_variables():
            var = self.get_attr(prop)
            var.val0 = units.ureg.Quantity(
                float(data[prop]),
                data[f"{prop}_unit"]
            )

        for fluid in self.fluid.is_var:
            self.fluid.val[fluid] = float(data[fluid])
            self.fluid.val0[fluid] = float(self.fluid.val[fluid])

    @classmethod
    def _result_attributes(cls):
        return ["m", "p", "h", "T", "v", "s", "vol", "x", "Td_bp"]

    @classmethod
    def _get_result_cols(cls, all_fluids):
        return [
            col for prop in cls._result_attributes()
            for col in [prop, f"{prop}_unit"]
        ] + list(all_fluids) + ['phase']

    @classmethod
    def _print_attributes(cls):
        return ["m", "p", "h", "T", "x", "phase"]

    def collect_results(self, all_fluids):
        return [
            _ for key in self._result_attributes()
            for _ in [self.get_attr(key).val, self.get_attr(key).unit]
        ] + [
            self.fluid.val[fluid] if fluid in self.fluid.val else np.nan
            for fluid in all_fluids
        ] + [
            self.phase.val
        ]

    def _adjust_to_property_limits(self, nw):
        r"""
        Check for invalid fluid property values.
        TODO: The network passed to this method should be putting the value
        limits to the connections in the preprocessing, then it can be
        omitted here.
        """
        fl = fp.single_fluid(self.fluid_data)

        # pure fluid
        if fl is not None:
            # pressure
            if self.p.is_var:
                self._adjust_pressure(fl)

            # enthalpy
            if self.h.is_var:
                self._adjust_enthalpy(fl)

                # two-phase related
                if (self.Td_bp.is_set or self.state.is_set or self.x.is_set or self.td_bubble.is_set or self.td_dew.is_set) and self.it < 10:
                    self._adjust_to_two_phase(fl)

        # mixture
        elif self.it < 5 and not self.good_starting_values:
            # pressure
            if self.p.is_var:
                if self.p.val_SI <= nw.p_range_SI[0]:
                    self.p.set_reference_val_SI(nw.p_range_SI[0])
                    logger.debug(self._property_range_message('p'))

                elif self.p.val_SI >= nw.p_range_SI[1]:
                    self.p.set_reference_val_SI(nw.p_range_SI[1])
                    logger.debug(self._property_range_message('p'))

            # enthalpy
            if self.h.is_var:
                if self.h.val_SI < nw.h_range_SI[0]:
                    self.h.set_reference_val_SI(nw.h_range_SI[0])
                    logger.debug(self._property_range_message('h'))

                elif self.h.val_SI > nw.h_range_SI[1]:
                    self.h.set_reference_val_SI(nw.h_range_SI[1])
                    logger.debug(self._property_range_message('h'))

                # temperature
                if self.T.is_set:
                    self._adjust_to_temperature_limits()

        # mass flow
        if self.m.is_var:
            if self.m.val_SI <= nw.m_range_SI[0]:
                self.m.set_reference_val_SI(nw.m_range_SI[0])
                logger.debug(self._property_range_message('m'))

            elif self.m.val_SI >= nw.m_range_SI[1]:
                self.m.set_reference_val_SI(nw.m_range_SI[1])
                logger.debug(self._property_range_message('m'))

    def _adjust_pressure(self, fluid):
        if self.p.val_SI > self.fluid.wrapper[fluid]._p_max:
            self.p.set_reference_val_SI(self.fluid.wrapper[fluid]._p_max)
            logger.debug(self._property_range_message('p'))

        elif self.p.val_SI < self.fluid.wrapper[fluid]._p_min:
            try:
                # if this works, the temperature is higher than the minimum
                # temperature, we can access pressure values below minimum
                # pressure
                self.fluid.wrapper[fluid].T_ph(self.p.val_SI, self.h.val_SI)
            except ValueError:
                self.p.set_reference_val_SI(self.fluid.wrapper[fluid]._p_min + 1e1)
                logger.debug(self._property_range_message('p'))

    def _adjust_enthalpy(self, fluid):
        # enthalpy
        try:
            hmin = self.fluid.wrapper[fluid].h_pT(
                self.p.val_SI, self.fluid.wrapper[fluid]._T_min + 1e-1
            )
        except ValueError:
            f = 1.05
            hmin = self.fluid.wrapper[fluid].h_pT(
                self.p.val_SI, self.fluid.wrapper[fluid]._T_min * f
            )
        if self.h.val_SI < hmin:
            if hmin < 0:
                self.h.set_reference_val_SI(hmin * 0.9999)
            else:
                self.h.set_reference_val_SI(hmin * 1.0001)
            logger.debug(self._property_range_message('h'))
        else:

            T = self.fluid.wrapper[fluid]._T_max
            # T_max depends on pressure for incompressibles
            while True:
                try:
                    hmax = self.fluid.wrapper[fluid].h_pT(self.p.val_SI, T)
                    break
                except ValueError as e:
                    T *= 0.99
                    if T < self.fluid.wrapper[fluid]._T_min:
                        raise ValueError(e) from e

            if self.h.val_SI > hmax:
                self.h.set_reference_val_SI(hmax * 0.9999)
                logger.debug(self._property_range_message('h'))

    def _adjust_to_two_phase(self, fluid):

        if self.p.val_SI > self.fluid.wrapper[fluid]._p_crit:
            self.p.set_reference_val_SI(self.fluid.wrapper[fluid]._p_crit * 0.9)
        # this is supposed to never be accessed with INCOMP backend but it is
        # not enforced. With INCOMP backend this causes a crash
        if self.Td_bp.is_set or self.state.is_set:
            if self.Td_bp.val_SI > 0 or self.state.val == 'g':
                h = self.fluid.wrapper[fluid].h_pQ(self.p.val_SI, 1)
                if self.h.val_SI < h:
                    self.h.set_reference_val_SI(h + 1e3)
                    logger.debug(self._property_range_message('h'))
            else:
                h = self.fluid.wrapper[fluid].h_pQ(self.p.val_SI, 0)
                if self.h.val_SI > h:
                    self.h.set_reference_val_SI(h - 1e3)
                    logger.debug(self._property_range_message('h'))

        elif self.td_bubble.is_set:
            h = self.fluid.wrapper[fluid].h_pQ(self.p.val_SI, 0)
            if self.td_bubble.val_SI >= 0:
                if self.h.val_SI > h:
                    self.h.set_reference_val_SI(h - 1e3)
            else:
                if self.h.val_SI < h:
                    self.h.set_reference_val_SI(h + 1e3)

        elif self.td_dew.is_set:
            h = self.fluid.wrapper[fluid].h_pQ(self.p.val_SI, 1)
            if self.td_dew.val_SI >= 0:
                if self.h.val_SI < h:
                    self.h.set_reference_val_SI(h + 1e3)
            else:
                if self.h.val_SI > h:
                    self.h.set_reference_val_SI(h - 1e3)

        elif self.x.is_set:
            h = self.fluid.wrapper[fluid].h_pQ(self.p.val_SI, self.x.val_SI)
            self.h.set_reference_val_SI(h)


    def _adjust_to_temperature_limits(self):
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
            f"{fpd[prop]['text'][0].upper()}{fpd[prop]['text'][1:]} out of "
            f"fluid property range at connection {self.label}, adjusting value "
            f"to {self.get_attr(prop).val_SI} {fpd[prop]['SI_unit']}."
        )
        return msg

    def _to_exerpy(self, pamb, Tamb):
        connection_json = {}

        self._get_physical_exergy(pamb, Tamb)

        connection_json[self.label] = {
            "source_component": self.source.label,
            "source_connector": int(self.source_id.removeprefix("out")) - 1,
            "target_component": self.target.label,
            "target_connector": int(self.target_id.removeprefix("in")) - 1
        }
        connection_json[self.label].update({f"mass_composition": self.fluid.val})
        connection_json[self.label].update({"kind": "material"})
        for param in ["m", "T", "p", "h", "s", "v"]:
            connection_json[self.label].update({
                param: self.get_attr(param).val_SI
            })
        connection_json[self.label].update(
            {"e_T": self.ex_therm, "e_M": self.ex_mech, "e_PH": self.ex_physical}
        )

        return connection_json

    def _get_physical_exergy(self, pamb, Tamb):
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
        self.ex_therm, self.ex_mech = fp.functions.calc_physical_exergy(
            self.h.val_SI, self.s.val_SI, self.p.val_SI, pamb, Tamb,
            self.fluid_data, self.mixing_rule, self.T.val_SI
        )
        self.Ex_therm = self.ex_therm * self.m.val_SI
        self.Ex_mech = self.ex_mech * self.m.val_SI

        self.ex_physical = self.ex_therm + self.ex_mech
        self.Ex_physical = self.m.val_SI * self.ex_physical

    def _get_chemical_exergy(self, pamb, Tamb, Chem_Ex):
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
            self.ex_chemical = fp.functions.calc_chemical_exergy(
                pamb, Tamb, self.fluid_data, Chem_Ex, self.mixing_rule,
                self.T.val_SI
            )

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

        msg = (
            f"Created reference object with factor {self.factor} and delta "
            f"{self.delta} referring to connection {ref_obj.label}"
        )
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
