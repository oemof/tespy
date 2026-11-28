import json

from pytest import approx
from pytest import raises

from tespy.components import Motor
from tespy.components import Pipe
from tespy.components import PowerSource
from tespy.components import Pump
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Splitter
from tespy.components import Subsystem
from tespy.components import Valve
from tespy.connections import Connection
from tespy.connections import PowerConnection
from tespy.networks import Network
from tespy.tools.helpers import TESPyComponentError
from tespy.tools.helpers import TESPyNetworkError


class SplitterSubsystem(Subsystem):
    def __init__(self, label):
        self.num_in = 1
        self.num_out = 2
        super().__init__(label)

    def create_network(self):
        splitter = Splitter("splitter", num_out=2)
        self.add_conns(
            Connection(self.inlet, "out1", splitter, "in1", label="feed"),
            Connection(splitter, "out1", self.outlet, "in1", label="first"),
            Connection(splitter, "out2", self.outlet, "in2", label="second"),
        )


def _network_with_splitter_subsystem():
    nw = Network(iterinfo=False)
    sub = SplitterSubsystem("sub")
    c_in = Connection(Source("source"), "out1", sub, "in1", label="c_in")
    c_out1 = Connection(sub, "out1", Sink("sink1"), "in1", label="c_out1")
    c_out2 = Connection(sub, "out2", Sink("sink2"), "in1", label="c_out2")
    nw.add_subsystems(sub)
    nw.add_conns(c_in, c_out1, c_out2)
    return nw, sub


def test_get_conn_miss_raises():
    sub = SplitterSubsystem("sub")
    with raises(KeyError, match="sub_feed"):
        sub.get_conn("sub_feed")


def test_get_comp_miss_raises():
    sub = SplitterSubsystem("sub")
    with raises(KeyError, match="sub_splitter"):
        sub.get_comp("sub_splitter")


def test_del_subsystems_with_boundary_leaves_valid_network():
    nw, sub = _network_with_splitter_subsystem()
    source = Source("other source")
    sink = Sink("other sink")
    nw.add_conns(Connection(source, "out1", sink, "in1", label="other"))

    nw.del_subsystems(sub, with_boundary=True)

    assert list(nw.conns.index) == ["other"]
    assert list(nw.comps.index) == ["other sink", "other source"]
    nw.check_topology()


def test_subsystem_label_with_separator_raises():
    with raises(ValueError, match="must not contain '/'"):
        SplitterSubsystem("a/b")



class Drive(Subsystem):
    def __init__(self, label):
        self.num_in = 1
        self.num_out = 1
        self.num_power_in = 1
        super().__init__(label)

    def create_network(self):
        pump = Pump("pump")
        motor = Motor("motor")
        self.add_conns(
            Connection(self.inlet, "out1", pump, "in1", label="suction"),
            Connection(pump, "out1", self.outlet, "in1", label="discharge"),
            PowerConnection(self.inlet, "power_out1", motor, "power_in", label="electric"),
            PowerConnection(motor, "power_out", pump, "power", label="shaft"),
        )


class Station(Subsystem):
    def __init__(self, label):
        self.num_in = 1
        self.num_out = 1
        self.num_power_in = 1
        super().__init__(label)

    def create_network(self):
        drive = Drive("drive")
        valve = Valve("valve")
        self.add_conns(
            Connection(self.inlet, "out1", drive, "in1", label="feed"),
            PowerConnection(self.inlet, "power_out1", drive, "power_in1", label="power"),
            Connection(drive, "out1", valve, "in1", label="discharge"),
            Connection(valve, "out1", self.outlet, "in1", label="throttled"),
        )


class Plant(Subsystem):
    def __init__(self, label):
        self.num_in = 1
        self.num_out = 1
        self.num_power_in = 1
        super().__init__(label)

    def create_network(self):
        station = Station("station")
        self.add_conns(
            Connection(self.inlet, "out1", station, "in1", label="feed"),
            PowerConnection(self.inlet, "power_out1", station, "power_in1", label="power"),
            Connection(station, "out1", self.outlet, "in1", label="product"),
        )


def _network_with(sub, solve=False):
    nw = Network(iterinfo=False)
    nw.units.set_defaults(temperature="degC", pressure="bar")
    c1 = Connection(Source("feed"), "out1", sub, "in1", label="1")
    c2 = Connection(sub, "out1", Sink("product"), "in1", label="2")
    e1 = PowerConnection(PowerSource("grid"), "power", sub, "power_in1", label="e1")
    nw.add_subsystems(sub)
    nw.add_conns(c1, c2, e1)
    if solve:
        c1.set_attr(fluid={"water": 1}, m=1, p=1, T=20)
        c2.set_attr(p=4)
        for s in sub._walk():
            for comp in s.comps.values():
                if isinstance(comp, Pump):
                    comp.set_attr(eta_s=0.8)
                elif isinstance(comp, Motor):
                    comp.set_attr(eta=0.95)
                elif isinstance(comp, Valve):
                    comp.set_attr(dp=1)
        nw.solve("design")
        nw.assert_convergence()
    return nw


def test_nested_labels():
    station = Station("station")
    drive = station.get_subsystem("drive")
    assert drive.parent is station
    assert drive.path == "station/drive"
    assert drive.get_comp("pump").label == "station_drive_pump"
    assert drive.get_conn("shaft").label == "station_drive_shaft"
    assert station.get_conn("discharge").label == "station_discharge"
    assert sorted(station.comps) == ["inlet", "outlet", "valve"]
    assert sorted(drive.comps) == ["inlet", "motor", "outlet", "pump"]


def test_three_levels_solve():
    plant = Plant("plant")
    nw = _network_with(plant, solve=True)
    assert nw.get_subsystem("plant/station/drive").path == "plant/station/drive"
    pump = nw.get_comp("plant_station_drive_pump")
    assert nw.subsystem_of(pump) is nw.get_subsystem("plant/station/drive")
    assert len(nw.conns) == 3 + 3 + 4 + 4


def test_membership_of_connections():
    station = Station("station")
    nw = _network_with(station)
    drive = station.get_subsystem("drive")

    expected = {
        "1": None,
        "2": None,
        "e1": None,
        "station_feed": station,
        "station_power": station,
        "station_discharge": station,
        "station_throttled": station,
        "station_drive_suction": drive,
        "station_drive_discharge": drive,
        "station_drive_electric": drive,
        "station_drive_shaft": drive,
    }
    assert {
        label: nw.subsystem_of(nw.get_conn(label)) for label in expected
    } == expected

    assert sorted(c.label for c in station.boundary_connections) == ["1", "2", "e1"]
    assert sorted(c.label for c in drive.boundary_connections) == [
        "station_discharge", "station_feed", "station_power"
    ]
    assert sorted(c.label for c in station.get_connections()) == [
        "station_discharge", "station_feed", "station_power", "station_throttled"
    ]
    assert len(station.get_connections(recursive=True)) == 8
    assert len(station.get_components(recursive=True)) == 7
    assert nw.subsystem_of(nw.get_comp("feed")) is None


def test_same_class_twice():
    nw = Network(iterinfo=False)
    a = Station("a")
    b = Station("b")
    nw.add_subsystems(a, b)
    nw.add_conns(
        Connection(Source("feed"), "out1", a, "in1", label="1"),
        Connection(a, "out1", b, "in1", label="2"),
        Connection(b, "out1", Sink("product"), "in1", label="3"),
        PowerConnection(PowerSource("grid1"), "power", a, "power_in1", label="e1"),
        PowerConnection(PowerSource("grid2"), "power", b, "power_in1", label="e2"),
    )
    nw.check_topology()
    assert nw.subsystem_of(nw.get_comp("b_drive_pump")) is b.get_subsystem("drive")


class TwoDrives(Subsystem):
    def __init__(self, label):
        self.num_in = 1
        self.num_out = 1
        self.num_power_in = 2
        super().__init__(label)

    def create_network(self):
        first = Drive("first")
        second = Drive("second")
        self.add_conns(
            Connection(self.inlet, "out1", first, "in1", label="feed"),
            Connection(first, "out1", second, "in1", label="between"),
            Connection(second, "out1", self.outlet, "in1", label="product"),
            PowerConnection(self.inlet, "power_out1", first, "power_in1", label="power1"),
            PowerConnection(self.inlet, "power_out2", second, "power_in1", label="power2"),
        )


def test_same_child_class_twice_in_parent():
    sub = TwoDrives("pumps")
    assert sorted(sub.subsystems) == ["first", "second"]
    assert sub.get_subsystem("second").get_comp("pump").label == "pumps_second_pump"


def test_single_component_subsystem():
    class Single(Subsystem):
        def __init__(self, label):
            self.num_in = 0
            self.num_out = 0
            super().__init__(label)

        def create_network(self):
            self.add_comps(Source("source"))

    sub = Single("single")
    assert sub.get_comp("source").label == "single_source"


def test_label_collision_names_both_paths():
    class Colliding(Subsystem):
        def __init__(self, label):
            self.num_in = 1
            self.num_out = 1
            self.num_power_in = 1
            super().__init__(label)

        def create_network(self):
            drive = Drive("drive")
            self.add_comps(Pump("drive_pump"))
            self.add_conns(
                Connection(self.inlet, "out1", drive, "in1", label="feed"),
                PowerConnection(self.inlet, "power_out1", drive, "power_in1", label="power"),
                Connection(drive, "out1", self.outlet, "in1", label="product"),
            )

    with raises(TESPyComponentError, match="'pump' of the subsystem a/drive"):
        Colliding("a")


def test_connecting_to_inner_component_raises():
    class Inner(Subsystem):
        def __init__(self, label):
            self.num_in = 1
            self.num_out = 1
            super().__init__(label)

        def create_network(self):
            drive = Drive("drive")
            self.add_conns(
                Connection(self.inlet, "out1", drive.get_comp("pump"), "in1", label="feed"),
            )

    with raises(TESPyComponentError, match="part of the subsystem drive"):
        Inner("inner")


def test_rename_and_reparent_rules():
    station = Station("station")
    drive = station.get_subsystem("drive")
    with raises(TESPyComponentError, match="cannot be renamed"):
        drive.label = "other"
    with raises(TESPyComponentError, match="already part of the subsystem station"):
        Station("other").add_subsystems(drive)
    with raises(TESPyComponentError, match="contains it"):
        drive.add_subsystems(station)
    with raises(ValueError, match="top level subsystem"):
        Network().add_subsystems(drive)

    drive = Drive("drive")
    drive.label = "renamed"
    assert drive.get_comp("pump").label == "renamed_pump"

    nw = _network_with(station)
    assert station._network() is nw
    with raises(TESPyComponentError, match="before adding it to a network"):
        station.add_subsystems(Drive("late"))
    with raises(TESPyComponentError, match="part of a network"):
        Station("other").add_subsystems(station)


def test_del_nested_subsystem():
    station = Station("station")
    nw = _network_with(station)
    nw.del_subsystems(station, with_boundary=True)
    assert len(nw.conns) == 0
    assert len(nw.comps) == 0
    assert station._network() is None
    nw2 = _network_with(station)
    nw2.check_topology()


def _line_network():
    nw = Network(iterinfo=False)
    nw.units.set_defaults(temperature="degC", pressure="bar", pressure_difference="bar")
    so = Source("source")
    pipe1 = Pipe("pipe 1")
    pipe2 = Pipe("pipe 2")
    pipe3 = Pipe("pipe 3")
    si = Sink("sink")
    c1 = Connection(so, "out1", pipe1, "in1", label="1")
    nw.add_conns(
        c1,
        Connection(pipe1, "out1", pipe2, "in1", label="2"),
        Connection(pipe2, "out1", pipe3, "in1", label="3"),
        Connection(pipe3, "out1", si, "in1", label="4"),
    )
    c1.set_attr(fluid={"water": 1}, m=1, p=5, T=20)
    for pipe in [pipe1, pipe2, pipe3]:
        pipe.set_attr(dp=0.5, Q=1000)
    return nw


def test_group_keeps_labels_and_results():
    nw = _line_network()
    nw.solve("design")
    before = nw.results["Connection"].copy()

    pipes = [nw.get_comp(f"pipe {i}") for i in [1, 2, 3]]
    line = nw.group("line", *pipes[:2])
    assert [c.label for c in line.get_components()] == ["pipe 1", "pipe 2"]
    assert nw.subsystem_of(nw.get_conn("2")) is line
    assert nw.subsystem_of(nw.get_conn("3")) is None
    assert sorted(c.label for c in line.boundary_connections) == ["1", "3"]
    assert list(line.conns) == ["2"]

    nw.solve("design")
    nw.assert_convergence()
    after = nw.results["Connection"]
    assert list(after.index) == list(before.index)
    columns = before.select_dtypes("number").columns
    assert after[columns].to_numpy() == approx(before[columns].to_numpy(), rel=1e-6)


def test_nested_groups_and_ungroup():
    nw = _line_network()
    pipes = [nw.get_comp(f"pipe {i}") for i in [1, 2, 3]]
    line = nw.group("line", *pipes)
    first = nw.group("first", pipes[0], pipes[1], parent=line)
    assert first.path == "line/first"
    assert nw.get_subsystem("line/first") is first
    assert nw.subsystem_of(nw.get_conn("2")) is first
    assert nw.subsystem_of(nw.get_conn("3")) is line
    assert list(first.conns) == ["2"]
    assert list(line.conns) == ["3"]

    nw.ungroup(line)
    assert list(nw.subsystems) == ["first"]
    assert first.parent is None
    assert nw.subsystem_of(pipes[2]) is None

    nw.ungroup(first)
    assert nw.subsystems == {}
    assert all(nw.subsystem_of(pipe) is None for pipe in pipes)


def test_group_of_subsystem_keeps_labels():
    station = Station("station")
    nw = _network_with(station)
    feed = nw.get_comp("feed")
    group = nw.group("supply", feed, station)
    assert station.path == "supply/station"
    assert station.get_comp("valve").label == "station_valve"
    assert nw.get_subsystem("supply/station/drive") is station.get_subsystem("drive")
    assert nw.subsystem_of(nw.get_conn("1")) is group
    assert sorted(c.label for c in group.boundary_connections) == ["2", "e1"]


def test_group_rules():
    nw = _line_network()
    pipe1 = nw.get_comp("pipe 1")
    line = nw.group("line", pipe1)
    with raises(TESPyNetworkError, match="not on the network level"):
        nw.group("other", pipe1)
    with raises(TESPyNetworkError, match="not part of the network"):
        nw.group("other", Pipe("pipe 1"))
    with raises(ValueError, match="already a subsystem"):
        nw.group("line", nw.get_comp("pipe 2"))

    station = Station("station")
    nw2 = _network_with(station)
    with raises(TESPyNetworkError, match="must be a group"):
        nw2.group("inner", parent=station)
    with raises(TESPyNetworkError, match="not on the network level"):
        nw2.group("inner", station.get_comp("valve"))
    with raises(TESPyNetworkError, match="Only groups"):
        nw2.ungroup(station)
    with raises(TESPyNetworkError, match="Only groups"):
        nw2.ungroup(line)


def test_group_ports():
    nw = Network(iterinfo=False)
    pump = Pump("pump")
    motor = Motor("motor")
    nw.add_conns(PowerConnection(motor, "power_out", pump, "power", label="shaft"))
    drive = nw.group("drive", pump, motor)
    drive.add_port("suction", pump, "in1")
    drive.add_port("discharge", pump, "out1")
    drive.add_port("electric", motor, "power_in")

    c1 = Connection(Source("feed"), "out1", drive, "suction", label="1")
    c2 = Connection(drive, "discharge", Sink("product"), "in1", label="2")
    e1 = PowerConnection(PowerSource("grid"), "power", drive, "electric", label="e1")
    assert (c1.target, c1.target_id) == (pump, "in1")
    assert (c2.source, c2.source_id) == (pump, "out1")
    assert (e1.target, e1.target_id) == (motor, "power_in")
    nw.add_conns(c1, c2, e1)
    nw.check_topology()

    with raises(ValueError, match="has no port inlet"):
        Connection(Source("other"), "out1", drive, "inlet")
    with raises(TESPyComponentError, match="already has a port"):
        drive.add_port("suction", pump, "in1")
    with raises(TESPyComponentError, match="has no port"):
        drive.add_port("other", pump, "in2")
    with raises(TESPyComponentError, match="not part of the group"):
        drive.add_port("other", Pump("pump 2"), "in1")
    with raises(TESPyComponentError, match="given by its interfaces"):
        Station("station").add_port("other", pump, "in1")


def test_interface_ports():
    station = Station("station")
    assert {key: (comp.label, port) for key, (comp, port) in station.ports.items()} == {
        "in1": ("station_inlet", "in1"),
        "power_in1": ("station_inlet", "power_in1"),
        "out1": ("station_outlet", "out1"),
    }
    with raises(ValueError, match="has no port in2"):
        Connection(Source("source"), "out1", station, "in2")


def _round_trip(nw):
    data = json.loads(json.dumps(nw.export()))
    return Network.from_dict(data)


def _membership(nw):
    def path(obj):
        sub = nw.subsystem_of(obj)
        return sub.path if sub is not None else ""

    return (
        {label: path(comp) for label, comp in nw.comps["object"].items()},
        {label: path(conn) for label, conn in nw.conns["object"].items()},
    )


def test_export_contains_tree():
    station = Station("station")
    nw = _network_with(station)
    data = nw.export()["Subsystem"]
    assert list(data) == ["station"]
    assert data["station"]["class"] == "Station"
    assert data["station"]["components"] == {
        "inlet": "station_inlet", "valve": "station_valve", "outlet": "station_outlet"
    }
    assert data["station"]["ports"]["power_in1"] == {
        "component": "station_inlet", "port": "power_in1"
    }
    drive = data["station"]["subsystems"]["drive"]
    assert drive["components"]["pump"] == "station_drive_pump"
    assert drive["connections"]["shaft"] == "station_drive_shaft"
    assert "subsystems" not in drive


def test_round_trip_keeps_tree_membership_and_results():
    nw = _network_with(Plant("plant"), solve=True)
    loaded = _round_trip(nw)
    assert _membership(loaded) == _membership(nw)
    assert loaded.export()["Subsystem"] == nw.export()["Subsystem"]

    drive = loaded.get_subsystem("plant/station/drive")
    assert drive.get_comp("pump") is loaded.get_comp("plant_station_drive_pump")
    assert drive.get_conn("shaft") is loaded.get_conn("plant_station_drive_shaft")
    assert type(drive) is Subsystem
    assert not drive.is_group

    loaded.solve("design")
    loaded.assert_convergence()
    columns = ["m", "p", "h", "T"]
    assert loaded.results["Connection"].loc[:, columns].to_numpy() == approx(
        nw.results["Connection"].loc[:, columns].to_numpy(), rel=1e-6
    )


def test_class_name_survives_two_round_trips():
    nw = _network_with(Station("station"))
    twice = _round_trip(_round_trip(nw))
    data = twice.export()["Subsystem"]["station"]
    assert data["class"] == "Station"
    assert data["subsystems"]["drive"]["class"] == "Drive"


def test_round_trip_of_groups_with_ports():
    nw = _line_network()
    pipes = [nw.get_comp(f"pipe {i}") for i in [1, 2, 3]]
    line = nw.group("line", *pipes)
    first = nw.group("first", pipes[0], pipes[1], parent=line)
    line.add_port("entry", pipes[0], "in1")

    loaded = _round_trip(nw)
    assert _membership(loaded) == _membership(nw)
    group = loaded.get_subsystem("line")
    assert group.is_group
    assert loaded.get_subsystem("line/first").is_group
    assert loaded.get_comp("pipe 1").label == "pipe 1"
    conn = Connection(Source("other"), "out1", group, "entry")
    assert conn.target is loaded.get_comp("pipe 1")
    assert list(first.conns) == ["2"]


def test_file_without_subsystem_section_loads():
    nw = _network_with(Station("station"))
    data = json.loads(json.dumps(nw.export()))
    del data["Subsystem"]
    loaded = Network.from_dict(data)
    assert loaded.subsystems == {}
    assert len(loaded.conns) == len(nw.conns)


def _tampered(change):
    nw = _network_with(Station("station"))
    data = json.loads(json.dumps(nw.export()))
    change(data["Subsystem"]["station"])
    return data


def test_reading_checks():
    def missing_member(sub):
        sub["components"]["other"] = "station_other"

    def two_subsystems(sub):
        sub["subsystems"]["drive"]["components"]["valve"] = "station_drive_valve"
        sub["components"]["drive_valve"] = "station_drive_valve"

    def wrong_level(sub):
        sub["connections"]["shaft"] = "station_drive_shaft"

    def port_outside(sub):
        sub["ports"]["other"] = {"component": "feed", "port": "out1"}

    def inconsistent_labels(sub):
        sub["components"]["valve"] = "valve"

    cases = [
        (missing_member, "component station_other of the subsystem station does not exist"),
        (wrong_level, "connection station_drive_shaft is listed in the subsystem station"),
        (port_outside, "port other of the subsystem station points to the component feed"),
        (inconsistent_labels, "do not follow from their local labels"),
    ]
    for change, message in cases:
        with raises(TESPyNetworkError, match=message):
            Network.from_dict(_tampered(change))


def test_reading_check_member_of_two_subsystems():
    nw = _line_network()
    nw.group("a", nw.get_comp("pipe 1"))
    nw.group("b", nw.get_comp("pipe 2"))
    data = json.loads(json.dumps(nw.export()))
    data["Subsystem"]["b"]["components"]["pipe 1"] = "pipe 1"
    with raises(TESPyNetworkError, match="member of the subsystems a and b"):
        Network.from_dict(data)


def test_from_network_keeps_nesting_and_places_twice():
    inner = _network_with(Station("station"))
    first = Subsystem.from_network("a", inner)
    second = Subsystem.from_network("b", inner)
    assert first.get_subsystem("station/drive").get_comp("pump").label == "a_station_drive_pump"

    nw = Network(iterinfo=False)
    nw.units.set_defaults(temperature="degC", pressure="bar")
    c1 = Connection(Source("feed"), "out1", first, "in1", label="1")
    c2 = Connection(first, "out1", second, "in1", label="2")
    c3 = Connection(second, "out1", Sink("product"), "in1", label="3")
    nw.add_subsystems(first, second)
    nw.add_conns(
        c1, c2, c3,
        PowerConnection(PowerSource("grid 1"), "power", first, "power_in1", label="e1"),
        PowerConnection(PowerSource("grid 2"), "power", second, "power_in1", label="e2"),
    )
    nw.check_topology()
    assert nw.subsystem_of(nw.get_comp("b_station_drive_pump")) is second.get_subsystem(
        "station/drive"
    )
    assert nw.subsystem_of(nw.get_conn("a_station_drive_shaft")).path == "a/station/drive"


def test_from_network_drops_boundary_members_of_groups():
    nw = _line_network()
    group = nw.group("supply", nw.get_comp("source"), nw.get_comp("pipe 1"))
    group.add_port("entry", nw.get_comp("source"), "out1")
    sub = Subsystem.from_network("line", nw)
    supply = sub.get_subsystem("supply")
    assert [c.label for c in supply.get_components()] == ["line_pipe 1"]
    assert supply.ports == {}
    assert supply.is_group


def test_results_path_column():
    nw = _network_with(Station("station"), solve=True)
    connections = nw.results["Connection"]["subsystem"]
    assert connections["1"] == ""
    assert connections["station_feed"] == "station"
    assert connections["station_drive_suction"] == "station/drive"
    assert nw.results["PowerConnection"].loc["station_drive_shaft", "subsystem"] == "station/drive"
    assert nw.results["Pump"].loc["station_drive_pump", "subsystem"] == "station/drive"
    assert nw.results["Valve"].loc["station_valve", "subsystem"] == "station"


def test_results_path_column_without_subsystems():
    nw = _line_network()
    nw.solve("design")
    assert (nw.results["Connection"]["subsystem"] == "").all()
    assert (nw.results["Pipe"]["subsystem"] == "").all()


def test_get_results_of_nested_subsystem_and_group():
    station = Station("station")
    nw = _network_with(station, solve=True)
    drive = station.get_subsystem("drive")

    results = nw.get_results(station, recursive=False)
    assert list(results["Connection"].index) == [
        "station_discharge", "station_feed", "station_throttled"
    ]
    assert list(results["Pump"].index) == []
    assert list(results["Valve"].index) == ["station_valve"]

    results = nw.get_results(station)
    assert len(results["Connection"]) == 5
    assert list(results["Pump"].index) == ["station_drive_pump"]

    results = nw.get_results(drive)
    assert list(results["PowerConnection"].index) == [
        "station_drive_electric", "station_drive_shaft"
    ]
    assert list(results["Pump"].index) == ["station_drive_pump"]
    assert list(results["Valve"].index) == []

    group = nw.group("supply", nw.get_comp("feed"), nw.get_comp("grid"), station)
    results = nw.get_results(group, recursive=False)
    assert list(results["Connection"].index) == ["1"]
    assert nw.results["Connection"].loc["1", "subsystem"] == "supply"
    assert nw.results["Pump"].loc["station_drive_pump", "subsystem"] == "supply/station/drive"
    nw.ungroup(group)
    assert nw.results["Connection"].loc["1", "subsystem"] == ""

    assert len(nw.get_results()["Connection"]) == len(nw.results["Connection"])


def test_print_results_of_nested_subsystem():
    station = Station("station")
    nw = _network_with(station, solve=True)
    nw.print_results(subsystem=station, print_results=False)
    nw.print_results(subsystem=station.get_subsystem("drive"), print_results=False)


def _line_subsystem_network():
    inner = _line_network()
    sub = Subsystem.from_network("plant", inner)
    nw = Network(iterinfo=False)
    nw.add_subsystems(sub)
    nw.add_conns(
        Connection(Source("source"), "out1", sub, "in1", label="in"),
        Connection(sub, "out1", Sink("sink"), "in1", label="out"),
    )
    return nw, sub


def test_export_after_deleting_subsystem_connection():
    nw, sub = _line_subsystem_network()
    conn = sub.get_conn("2")
    nw.del_conns(conn)
    assert conn not in sub.conns.values()
    assert "plant_2" not in json.dumps(nw.export()["Subsystem"])
    nw.add_conns(
        Connection(nw.get_comp("plant_pipe 1"), "out1", nw.get_comp("plant_pipe 2"), "in1", label="new")
    )
    loaded = _round_trip(nw)
    assert loaded.subsystem_of(loaded.get_conn("new")).path == "plant"


def test_export_after_deleting_group_connection_and_component():
    nw = _line_network()
    pipe1 = nw.get_comp("pipe 1")
    group = nw.group("line", pipe1, nw.get_comp("pipe 2"))
    group.add_port("entry", pipe1, "in1")
    nw.del_conns(nw.get_conn("2"))
    assert list(group.conns) == []
    assert json.loads(json.dumps(nw.export()))["Subsystem"]["line"]["connections"] == {}

    nw.del_conns(nw.get_conn("1"))
    assert pipe1 not in group.get_components()
    assert group.ports == {}
    assert nw.subsystem_of(pipe1) is None
    data = json.loads(json.dumps(nw.export()))
    assert data["Subsystem"]["line"]["components"] == {"pipe 2": "pipe 2"}


def test_group_inside_prefixing_subsystem():
    inner = _line_network()
    inner.group("line", *[inner.get_comp(f"pipe {i}") for i in [1, 2]])
    sub = Subsystem.from_network("plant", inner)
    nw = Network(iterinfo=False)
    nw.add_subsystems(sub)
    nw.add_conns(
        Connection(Source("source"), "out1", sub, "in1", label="in"),
        Connection(sub, "out1", Sink("sink"), "in1", label="out"),
    )
    line = nw.get_subsystem("plant/line")
    pipe1 = nw.get_comp("plant_pipe 1")

    first = nw.group("first", pipe1, parent=line)
    assert first.comps == {"pipe 1": pipe1}
    assert pipe1.label == "plant_pipe 1"

    nw.ungroup(first)
    nw.ungroup(line)
    assert pipe1.label == "plant_pipe 1"
    assert nw.get_comp("plant_pipe 1") is pipe1
    assert sub.get_comp("pipe 1") is pipe1
    assert sub.get_conn("2").label == "plant_2"
    _round_trip(nw)
