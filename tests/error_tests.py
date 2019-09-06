# -*- coding: utf-8

from nose.tools import eq_, raises

from tespy import nwk, cmp, con, hlp, subsys, cmp_char


# %% bulk tests


class specification_error_tests:

    def setup(self):
        self.nw = nwk.network(['water', 'air'])
        self.comp = cmp.cogeneration_unit('cogeneration unit')
        self.pipe = cmp.pipe('pipe')
        self.conn = con.connection(self.comp, 'out1', self.pipe, 'in1')
        self.bus = con.bus('mybus')
        self.sub = subsys.subsystem('MySub')

    @raises(ValueError)
    def cmp_instanciation_ValueError(self, label):
        cmp.cogeneration_unit(label)

    @raises(ValueError)
    def set_attr_ValueError(self, instance, **kwargs):
        instance.set_attr(**kwargs)

    @raises(ValueError)
    def create_connection_ValueError(self, pos):
        if pos == 'source':
            con.connection(self.comp, 'out6', self.pipe, 'in1')
        elif pos == 'target':
            con.connection(self.comp, 'out1', self.pipe, 'in5')

    @raises(TypeError)
    def set_attr_TypeError(self, instance, **kwargs):
        instance.set_attr(**kwargs)

    @raises(TypeError)
    def bus_add_comps_TypeError(self, c):
        self.bus.add_comps(c)

    @raises(TypeError)
    def test_create_conn_TypeError(self):
        con.connection(self.comp, 'out1', self.conn, 'in1')

    @raises(TypeError)
    def create_ref_TypeError(self, params):
        con.ref(params[0], params[1], params[2])

    @raises(KeyError)
    def set_attr_KeyError(self, instance, **kwargs):
        instance.set_attr(**kwargs)

    @raises(KeyError)
    def get_attr_KeyError(self, instance, key):
        instance.get_attr(key)

    @raises(hlp.TESPyComponentError)
    def test_TESPyComponentError(self):
        self.comp.set_attr(interface=True)

    @raises(hlp.TESPyConnectionError)
    def test_con_TESPyConnectionError(self):
        con.connection(self.comp, 'out1', self.comp, 'in1')

    @raises(hlp.TESPyConnectionError)
    def test_bus_TESPyConnectionError(self):
        self.bus.add_comps(self.comp)

    @raises(hlp.TESPyNetworkError)
    def test_network_bus_duplicate(self):
        self.nw.add_busses(self.bus, self.bus)

    @raises(hlp.TESPyNetworkError)
    def test_network_buslabel_duplicate(self):
        bus = con.bus('mybus')
        self.nw.add_busses(self.bus)
        self.nw.add_busses(bus)

    @raises(TypeError)
    def test_network_bus_type(self):
        self.nw.add_busses(self.conn)

    def test_set_attr_errors(self):
        #
        labels = [5, 'Label,', 'Labe;l', 'Label.']
        for l in labels:
            self.cmp_instanciation_ValueError(l)

        # ValueErrors
        self.set_attr_ValueError(self.comp, offdesign=['Q'])

        self.set_attr_ValueError(self.conn, offdesign=['f'])
        self.set_attr_ValueError(self.conn, state='f')

        self.set_attr_ValueError(self.nw, m_unit='kg')
        self.set_attr_ValueError(self.nw, h_unit='kg')
        self.set_attr_ValueError(self.nw, p_unit='kg')
        self.set_attr_ValueError(self.nw, T_unit='kg')
        self.set_attr_ValueError(self.nw, v_unit='kg')

        self.create_connection_ValueError('source')
        self.create_connection_ValueError('target')

        # TypeErrors
        self.set_attr_TypeError(self.comp, P=[5])
        self.set_attr_TypeError(self.comp, tiP_char=None)
        self.set_attr_TypeError(self.comp, design='f')
        self.set_attr_TypeError(self.comp, fuel=hlp.dc_cp(val='CH4'))

        self.set_attr_TypeError(self.conn, design='h')
        self.set_attr_TypeError(self.conn, fluid_balance=1)
        self.set_attr_TypeError(self.conn, h0=[4])
        self.set_attr_TypeError(self.conn, fluid=5)
        self.set_attr_TypeError(self.conn, state=5)

        self.set_attr_TypeError(self.nw, m_range=5)
        self.set_attr_TypeError(self.nw, p_range=5)
        self.set_attr_TypeError(self.nw, h_range=5)
        self.set_attr_TypeError(self.nw, T_range=5)

        self.bus_add_comps_TypeError({'c': self.conn})
        self.bus_add_comps_TypeError({'f': self.comp})
        self.bus_add_comps_TypeError({'c': self.comp, 'char': 'Hi'})
        self.bus_add_comps_TypeError({'c': self.comp, 'p': 5})
        self.bus_add_comps_TypeError({'c': self.comp, 'P_ref': 'what'})

        self.create_ref_TypeError([self.conn, 7, 'hi'])
        self.create_ref_TypeError([self.conn, 'hi', 0])
        self.create_ref_TypeError([self.comp, 1, 0])

        # KeyErrors
        self.set_attr_KeyError(self.comp, wow=5)
        self.set_attr_KeyError(self.conn, jey=5)
        self.set_attr_KeyError(self.sub, a=7)

        self.get_attr_KeyError(self.comp, 'wow')
        self.get_attr_KeyError(self.conn, 'key')
        self.get_attr_KeyError(self.bus, 'components')
        self.get_attr_KeyError(con.ref(self.conn, 1, 0), 'comp')
        self.get_attr_KeyError(self.nw, 'test')
        self.get_attr_KeyError(self.sub, 'test')
        self.get_attr_KeyError(cmp_char.characteristics(), 'test')
        self.get_attr_KeyError(hlp.data_container(), 'somekey')


# %% Single tests

@raises(ValueError)
def test_interface_ValueError():
    # interface specification
    cmp.sink('sink', interface=5)

##############################################################################
# networks


@raises(hlp.TESPyNetworkError)
def test_network_instanciation_no_fluids():
    nw = nwk.network([])
    so = cmp.source('source')
    si = cmp.sink('sink')
    conn = con.connection(so, 'out1', si, 'in1')
    nw.add_conns(conn)
    nw.solve('design', init_only=True)


@raises(ValueError)
def test_network_print_level():
    nwk.network(['INCOMP::DowQ']).set_printoptions(print_level='error')


@raises(TypeError)
def test_network_instanciation_single_fluid():
    nwk.network('water')


@raises(TypeError)
def test_network_add_conns():
    nwk.network(['water']).add_conns(cmp.component('test'))


@raises(hlp.TESPyNetworkError)
def test_network_connection_error_source():
    nw = nwk.network(['water'])
    source = cmp.source('source')
    sink1 = cmp.sink('sink1')
    sink2 = cmp.sink('sink2')
    a = con.connection(source, 'out1', sink1, 'in1')
    b = con.connection(source, 'out1', sink2, 'in1')
    nw.add_conns(a, b)
    nw.check_network()


@raises(hlp.TESPyNetworkError)
def test_network_connection_error_target():
    nw = nwk.network(['water'])
    source1 = cmp.source('source1')
    source2 = cmp.source('source2')
    sink = cmp.sink('sink')
    a = con.connection(source1, 'out1', sink, 'in1')
    b = con.connection(source2, 'out1', sink, 'in1')
    nw.add_conns(a, b)
    nw.check_network()


@raises(hlp.TESPyNetworkError)
def test_network_network_consistency_inlets():
    nw = nwk.network(['water'])
    merge = cmp.merge('merge')
    sink = cmp.sink('label')
    a = con.connection(merge, 'out1', sink, 'in1')
    nw.add_conns(a)
    nw.check_network()


@raises(hlp.TESPyNetworkError)
def test_network_network_consistency_outlets():
    nw = nwk.network(['water', 'air'])
    source = cmp.source('source')
    splitter = cmp.splitter('splitter')
    a = con.connection(source, 'out1', splitter, 'in1')
    nw.add_conns(a)
    nw.check_network()


@raises(hlp.TESPyNetworkError)
def test_network_component_labels():
    nw = nwk.network(['water'])
    source = cmp.source('label')
    sink = cmp.sink('label')
    a = con.connection(source, 'out1', sink, 'in1')
    nw.add_conns(a)
    nw.check_network()


@raises(hlp.TESPyNetworkError)
def test_network_offdesign_path():
    nw = nwk.network(['water'])
    source = cmp.source('source')
    sink = cmp.sink('sink')
    a = con.connection(source, 'out1', sink, 'in1')
    nw.add_conns(a)
    nw.solve('offdesign')


@raises(ValueError)
def test_network_mode():
    nw = nwk.network(['water'])
    source = cmp.source('source')
    sink = cmp.sink('sink')
    a = con.connection(source, 'out1', sink, 'in1')
    nw.add_conns(a)
    nw.solve('ofdesign')


@raises(hlp.TESPyNetworkError)
def test_network_underdetermination():
    nw = nwk.network(['water'])
    source = cmp.source('source')
    sink = cmp.sink('sink')
    a = con.connection(source, 'out1', sink, 'in1', m=1)
    nw.add_conns(a)
    nw.solve('design')


@raises(hlp.TESPyNetworkError)
def test_network_overdetermination():
    nw = nwk.network(['water'])
    source = cmp.source('source')
    sink = cmp.sink('sink')
    a = con.connection(source, 'out1', sink, 'in1', m=1, p=1e5, x=1, h=1e6, fluid={'water': 1}, fluid_balance=True)
    nw.add_conns(a)
    nw.solve('design')


def test_network_linear_dependency():
    nw = nwk.network(['water'])
    source = cmp.source('source')
    sink = cmp.sink('sink')
    a = con.connection(source, 'out1', sink, 'in1', m=1, p=1e5, h=1e6, x=1)
    nw.add_conns(a)
    nw.solve('design')
    eq_(nw.lin_dep, True, 'This test must result in a linear dependency of the jacobian matrix.')


def test_network_no_progress():
    nw = nwk.network(['water'])
    source = cmp.source('source')
    pipe = cmp.pipe('pipe', pr=1, Q=-100e3)
    sink = cmp.sink('sink')
    a = con.connection(source, 'out1', pipe, 'in1', m=1, p=1e5, T=280, fluid={'water': 1})
    b = con.connection(pipe, 'out1', sink, 'in1')
    nw.add_conns(a, b)
    nw.solve('design')
    eq_(nw.progress, False, 'This test must result in a calculation making no progress, as the pipe\'s outlet enthalpy is below fluid property range.')


def test_network_max_iter():
    nw = nwk.network(['water'])
    source = cmp.source('source')
    pipe = cmp.pipe('pipe', pr=1, Q=100e3)
    sink = cmp.sink('sink')
    a = con.connection(source, 'out1', pipe, 'in1', m=1, p=1e5, T=280, fluid={'water': 1})
    b = con.connection(pipe, 'out1', sink, 'in1')
    nw.add_conns(a, b)
    nw.solve('design', max_iter=2)
    eq_(nw.max_iter, nw.iter + 1, 'This test must result in the itercount being equal to the max iter statement.')

##############################################################################
# subsystems


@raises(ValueError)
def test_subsys_label_str():
    subsys.subsystem(5)


@raises(ValueError)
def test_subsys_label_forbidden():
    subsys.subsystem('label;')

##############################################################################
# characteristics


@raises(KeyError)
def test_char_missing_key():
    cmp_char.characteristics(a=6)


@raises(KeyError)
def test_char_missing_key():
    cmp_char.characteristics(a=6)


@raises(ValueError)
def test_char_number_of_points():
    cmp_char.characteristics(x=[0, 1, 2], y=[1, 2, 3, 4])


@raises(KeyError)
def test_char_map_missing_key():
    cmp_char.char_map(a=6)


@raises(ValueError)
def test_char_map_number_of_points():
    cmp_char.char_map(x=[0, 1, 2], y=[[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3]])


@raises(ValueError)
def test_char_map_number_of_dimensions():
    cmp_char.char_map(x=[0, 1, 2], y=[[1, 2, 3, 4], [1, 2, 3, 4]])

##############################################################################
# helpers


@raises(TypeError)
def test_tespy_fluid_alias_type():
    hlp.tespy_fluid(5, {'water': 1}, [0, 1], [0, 1])


@raises(ValueError)
def test_tespy_fluid_alias_value():
    hlp.tespy_fluid('IDGAS::water', {'water': 1}, [0, 1], [0, 1])


##############################################################################
# components


class combustion_chamber_error_tests:

    @raises(hlp.TESPyComponentError)
    def test_combustion_chamber_missing_fuel(self):
        """
        Test no fuel in network.
        """
        nw = nwk.network(['H2O', 'N2', 'O2', 'Ar', 'CO2'])
        comb = cmp.combustion_chamber('combustion chamber')
        c1 = con.connection(cmp.source('air'), 'out1', comb, 'in1')
        c2 = con.connection(cmp.source('fuel'), 'out1', comb, 'in2')
        c3 = con.connection(comb, 'out1', cmp.sink('flue gas'), 'in1')
        nw.add_conns(c1, c2, c3)
        nw.solve('design', init_only=True)


class combustion_chamber_stoich_error_tests:

    def setup(self):
        self.nw = nwk.network(['TESPy::fuel', 'TESPy::fuel_fg', 'Air'])
        self.comb = cmp.combustion_chamber_stoich('combustion chamber')
        c1 = con.connection(cmp.source('air'), 'out1', self.comb, 'in1')
        c2 = con.connection(cmp.source('fuel'), 'out1', self.comb, 'in2')
        c3 = con.connection(self.comb, 'out1', cmp.sink('flue gas'), 'in1')
        self.nw.add_conns(c1, c2, c3)

    @raises(hlp.TESPyComponentError)
    def test_cc_stoich_missing_fuel(self):
        """
        Test missing fuel composition.
        """
        self.nw.solve('design', init_only=True)

    @raises(hlp.TESPyComponentError)
    def test_cc_stoich_missing_fuel_alias(self):
        """
        Test missing fuel alias.
        """
        self.comb.set_attr(fuel={'CH4': 1})
        self.nw.solve('design', init_only=True)

    @raises(hlp.TESPyComponentError)
    def test_cc_stoich_bad_fuel_alias(self):
        """
        Test bad name for fuel alias.
        """
        self.comb.set_attr(fuel={'CH4': 1}, fuel_alias='TESPy::fuel')
        self.nw.solve('design', init_only=True)

    @raises(hlp.TESPyComponentError)
    def test_cc_stoich_missing_air(self):
        """
        Test missing air composition.
        """
        self.comb.set_attr(fuel={'CH4': 1}, fuel_alias='fuel')
        self.nw.solve('design', init_only=True)

    @raises(hlp.TESPyComponentError)
    def test_cc_stoich_missing_air_alias(self):
        """
        Test missing air alias.
        """
        self.comb.set_attr(fuel={'CH4': 1}, fuel_alias='fuel', air={'N2': 0.76, 'O2': 0.24})
        self.nw.solve('design', init_only=True)

    @raises(hlp.TESPyComponentError)
    def test_cc_stoich_bad_air_alias(self):
        """
        Test bad name for air alias.
        """
        self.comb.set_attr(fuel={'CH4': 1}, fuel_alias='fuel', air={'N2': 0.76, 'O2': 0.24}, air_alias='TESPy::air')
        self.nw.solve('design', init_only=True)

    @raises(hlp.TESPyComponentError)

    def test_cc_stoich_missing_oxygen(self):
        """
        Test bad name for air alias.
        """
        self.comb.set_attr(fuel={'CH4': 1}, fuel_alias='fuel', air={'N2': 1}, air_alias='myair')
        self.nw.solve('design', init_only=True)


class cogeneration_unit_bus_error_tests:

    def setup(self):
        self.nw = nwk.network(['water', 'air'])
        self.comp = cmp.cogeneration_unit('cogeneration unit')
        self.bus = con.bus('power')
        self.bus.add_comps({'c': self.comp, 'p': 'Param'})

    @raises(ValueError)
    def test_missing_bus_param_func(self):
        """
        Test missing bus parameter in bus function for cogeneration unit.
        """
        self.comp.bus_func(self.bus.comps.loc[self.comp])

    @raises(ValueError)
    def test_missing_bus_param_deriv(self):
        """
        Test missing bus parameter in bus derivatives for cogeneration unit.
        """
        self.comp.num_vars = 2
        self.comp.inl = [con.connection(self.comp, 'out1', cmp.sink('sink'), 'in1')]
        self.comp.inl[0].fluid = hlp.dc_flu(val={'water': 1})
        self.comp.bus_deriv(self.bus.comps.loc[self.comp])

