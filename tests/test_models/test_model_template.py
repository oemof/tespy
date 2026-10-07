"""Unit tests for ModelTemplate base class."""
from unittest.mock import MagicMock
from unittest.mock import call
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from tespy.models.template import ModelTemplate

# ---------------------------------------------------------------------------
# Minimal concrete subclass - no real network operations needed for most tests
# ---------------------------------------------------------------------------

class SimpleModel(ModelTemplate):
    """Concrete subclass that exposes all four parameter forms."""

    def _parameter_lookup(self):
        return {
            "pressure": ["Connections", "c1", "p"],
            "efficiency": ["Components", "comp1", "eta_s"],
            "cop": {"get": self._calc_cop},
            "ude_target": {"set": self._set_ude},
            "custom": {"get": self._get_custom, "set": self._set_custom},
        }

    # helpers used in lookup callables
    def _calc_cop(self):
        return 3.5

    def _set_ude(self, value):
        self._ude_value = value

    def _get_custom(self):
        return getattr(self, "_custom_value", 0.0)

    def _set_custom(self, value):
        self._custom_value = value

    def _create_network(self):
        super()._create_network()
        # replace the real network with a mock so tests stay unit-level
        self.nw = MagicMock()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mock_conn(value):
    """Return a mock connection whose .get_attr(x).val == value."""
    conn = MagicMock()
    conn.get_attr.return_value.val = value
    return conn


def _mock_comp(value):
    """Return a mock component whose .get_attr(x).val == value."""
    comp = MagicMock()
    comp.get_attr.return_value.val = value
    return comp


# ---------------------------------------------------------------------------
# get_parameter
# ---------------------------------------------------------------------------

class TestGetParameter:

    def setup_method(self):
        self.model = SimpleModel()
        self.model.nw.get_conn.return_value = _mock_conn(5.0)
        self.model.nw.get_comp.return_value = _mock_comp(0.9)

    def test_connection_attribute(self):
        val = self.model.get_parameter("pressure")
        self.model.nw.get_conn.assert_called_with("c1")
        assert val == 5.0

    def test_component_attribute(self):
        val = self.model.get_parameter("efficiency")
        self.model.nw.get_comp.assert_called_with("comp1")
        assert val == 0.9

    def test_get_only_callable(self):
        assert self.model.get_parameter("cop") == pytest.approx(3.5)

    def test_get_set_callable_reads_via_get(self):
        self.model._custom_value = 42.0
        assert self.model.get_parameter("custom") == pytest.approx(42.0)

    def test_write_only_raises_attribute_error(self):
        with pytest.raises(AttributeError, match="write-only"):
            self.model.get_parameter("ude_target")

    def test_missing_key_raises_key_error(self):
        with pytest.raises(KeyError, match="nonexistent"):
            self.model.get_parameter("nonexistent")


# ---------------------------------------------------------------------------
# set_parameters
# ---------------------------------------------------------------------------

class TestSetParameters:

    def setup_method(self):
        self.model = SimpleModel()
        mock_conn = MagicMock()
        mock_comp = MagicMock()
        self.model.nw.get_conn.return_value = mock_conn
        self.model.nw.get_comp.return_value = mock_comp
        self.mock_conn = mock_conn
        self.mock_comp = mock_comp

    def test_sets_connection_attribute(self):
        self.model.set_parameters(pressure=10.0)
        self.model.nw.get_conn.assert_called_with("c1")
        self.mock_conn.set_attr.assert_called_with(p=10.0)

    def test_sets_component_attribute(self):
        self.model.set_parameters(efficiency=0.85)
        self.model.nw.get_comp.assert_called_with("comp1")
        self.mock_comp.set_attr.assert_called_with(eta_s=0.85)

    def test_set_only_callable(self):
        self.model.set_parameters(ude_target=99.0)
        assert self.model._ude_value == pytest.approx(99.0)

    def test_get_set_callable_writes_via_set(self):
        self.model.set_parameters(custom=7.0)
        assert self.model._custom_value == pytest.approx(7.0)

    def test_read_only_raises_attribute_error(self):
        with pytest.raises(AttributeError, match="read-only"):
            self.model.set_parameters(cop=1.0)

    def test_missing_key_raises_key_error(self):
        with pytest.raises(KeyError, match="nonexistent"):
            self.model.set_parameters(nonexistent=1.0)

    def test_multiple_params_on_same_connection_batched(self):
        """Two parameters on the same connection must result in one set_attr call."""

        class TwoParamModel(ModelTemplate):
            def _parameter_lookup(self):
                return {
                    "pressure": ["Connections", "c1", "p"],
                    "temperature": ["Connections", "c1", "T"],
                }

            def _create_network(self):
                super()._create_network()
                self.nw = MagicMock()

        m = TwoParamModel()
        mock_conn = MagicMock()
        m.nw.get_conn.return_value = mock_conn

        m.set_parameters(pressure=10.0, temperature=200.0)

        # get_conn should be called once per unique label
        m.nw.get_conn.assert_called_once_with("c1")
        # set_attr should receive both attributes in one call
        mock_conn.set_attr.assert_called_once_with(p=10.0, T=200.0)


# ---------------------------------------------------------------------------
# _order_min_change
# ---------------------------------------------------------------------------

class TestOrderMinChange:

    def test_trivial_single_point(self):
        pts = np.array([[1.0, 2.0]])
        order = ModelTemplate._order_min_change(pts, start_idx=0)
        assert list(order) == [0]

    def test_two_points_start_first(self):
        pts = np.array([[0.0], [10.0]])
        order = ModelTemplate._order_min_change(pts, start_idx=0)
        assert list(order) == [0, 1]

    def test_two_points_start_second(self):
        pts = np.array([[0.0], [10.0]])
        order = ModelTemplate._order_min_change(pts, start_idx=1)
        assert list(order) == [1, 0]

    def test_visits_all_points(self):
        rng = np.random.default_rng(0)
        pts = rng.random((10, 3))
        order = ModelTemplate._order_min_change(pts, start_idx=3)
        assert sorted(order) == list(range(10))

    def test_starts_at_requested_index(self):
        pts = np.array([[0.0], [5.0], [10.0]])
        order = ModelTemplate._order_min_change(pts, start_idx=2)
        assert order[0] == 2

    def test_greedy_nearest_neighbour(self):
        # points on a line; greedy from idx 0 should go 0->1->2
        pts = np.array([[0.0], [1.0], [100.0]])
        order = ModelTemplate._order_min_change(pts, start_idx=0)
        assert list(order) == [0, 1, 2]

# ---------------------------------------------------------------------------
# get_objectives
# ---------------------------------------------------------------------------

class TestGetObjectives:

    def test_returns_list_of_parameter_values(self):
        model = SimpleModel()
        model.nw = MagicMock()
        model.nw.get_conn.return_value = _mock_conn(5.0)
        model.nw.get_comp.return_value = _mock_comp(0.9)

        result = model.get_objectives(["pressure", "efficiency"])
        assert result == pytest.approx([5.0, 0.9])

    def test_callable_objective(self):
        model = SimpleModel()
        model.nw = MagicMock()
        result = model.get_objectives(["cop"])
        assert result == pytest.approx([3.5])


# ---------------------------------------------------------------------------
# sensitivity_analysis
# ---------------------------------------------------------------------------

class TestSensitivityAnalysis:

    def _make_model(self, side_effects):
        """Return a SimpleModel whose solve_model_design and get_parameter are mocked."""
        model = SimpleModel()
        model.nw = MagicMock()

        call_counter = {"n": 0}

        def fake_solve(**kwargs):
            model.nw.status = 0

        def fake_get_parameter(name):
            # return current state values for initial distance computation
            return side_effects.get(name, 0.0)

        model.solve_model_design = MagicMock(side_effect=fake_solve)
        model.get_parameter = MagicMock(side_effect=fake_get_parameter)
        model.get_results = MagicMock(return_value={"result": 1.0})
        return model

    def test_returns_dataframe_with_correct_columns(self):
        model = self._make_model({"pressure": 5.0, "efficiency": 0.8})
        result = model.sensitivity_analysis(
            param_dict={"pressure": [5.0, 10.0, 15.0], "efficiency": [0.7, 0.8, 0.9]},
            result_param_list=["result"],
        )
        assert isinstance(result, pd.DataFrame)
        assert "pressure" in result.columns
        assert "efficiency" in result.columns
        assert "result" in result.columns

    def test_returns_correct_number_of_rows(self):
        model = self._make_model({"pressure": 5.0, "efficiency": 0.8})
        result = model.sensitivity_analysis(
            param_dict={"pressure": [5.0, 10.0, 15.0], "efficiency": [0.7, 0.8, 0.9]},
            result_param_list=[],
        )
        assert len(result) == 3

    def test_calls_solve_once_per_row(self):
        model = self._make_model({"pressure": 5.0})
        model.sensitivity_analysis(
            param_dict={"pressure": [5.0, 10.0, 15.0, 20.0]},
            result_param_list=[],
        )
        assert model.solve_model_design.call_count == 4

    def test_raises_on_unequal_parameter_lengths(self):
        model = SimpleModel()
        model.nw = MagicMock()
        with pytest.raises(ValueError, match="same number"):
            model.sensitivity_analysis(
                param_dict={"pressure": [1.0, 2.0], "efficiency": [0.7]},
                result_param_list=[],
            )

    def test_raises_without_param_dict(self):
        model = SimpleModel()
        with pytest.raises(ValueError):
            model.sensitivity_analysis(result_param_list=[])

    def test_output_sorted_by_original_index(self):
        """Results must be reordered back to the original param_dict order."""
        model = self._make_model({"pressure": 15.0})  # current state near last value
        result = model.sensitivity_analysis(
            param_dict={"pressure": [5.0, 10.0, 15.0]},
            result_param_list=[],
        )
        assert list(result["pressure"]) == pytest.approx([5.0, 10.0, 15.0])

    def test_none_result_on_failed_solve(self):
        model = SimpleModel()
        model.nw = MagicMock()

        def failing_solve(**kwargs):
            model.nw.status = 3

        model.solve_model_design = MagicMock(side_effect=failing_solve)
        model.get_parameter = MagicMock(return_value=5.0)
        model.get_results = MagicMock(return_value={"result": 1.0})

        result = model.sensitivity_analysis(
            param_dict={"pressure": [5.0]},
            result_param_list=["result"],
        )
        assert result["result"].iloc[0] is None


# ---------------------------------------------------------------------------
# diagram cache isolation
# ---------------------------------------------------------------------------

def test_diagram_cache_is_per_instance():
    """Each ModelTemplate instance must have its own cache dict."""
    m1 = SimpleModel()
    m2 = SimpleModel()
    m1._diagram_cache["fluid"] = "something"
    assert "fluid" not in m2._diagram_cache


# ---------------------------------------------------------------------------
# diagram creation, registration and graceful degradation
# ---------------------------------------------------------------------------

class TestDiagramHandling:

    def setup_method(self):
        self.model = SimpleModel()
        self.model.nw.conns.index = ["c1"]

    def test_create_diagram_is_used_and_cached(self):
        diagram = MagicMock()
        self.model.create_diagram = MagicMock(return_value=diagram)
        assert self.model._get_diagram("water") is diagram
        assert self.model._get_diagram("water") is diagram
        self.model.create_diagram.assert_called_once_with("water")

    def test_register_diagram_takes_precedence(self):
        diagram = MagicMock()
        self.model.create_diagram = MagicMock()
        self.model.register_diagram("water", diagram)
        assert self.model._get_diagram("water") is diagram
        self.model.create_diagram.assert_not_called()

    def test_failing_creation_warns_and_caches_none(self, caplog):
        import logging
        self.model.create_diagram = MagicMock(side_effect=ValueError("nope"))
        with caplog.at_level(logging.WARNING):
            assert self.model._get_diagram("water") is None
        assert "water" in caplog.text
        # second call must not retry the creation
        assert self.model._get_diagram("water") is None
        self.model.create_diagram.assert_called_once()

    def test_failing_creation_raises_when_strict(self):
        self.model.create_diagram = MagicMock(side_effect=ValueError("nope"))
        with pytest.raises(ValueError, match="nope"):
            self.model._get_diagram("water", strict=True)

    def test_cached_failure_raises_when_strict(self):
        self.model._diagram_cache["water"] = None
        with pytest.raises(ValueError, match="water"):
            self.model._get_diagram("water", strict=True)

    def test_supports_fluid_diagram_unknown_connection(self):
        assert not self.model.supports_fluid_diagram("nonexistent")

    def test_supports_fluid_diagram_mixture(self):
        with patch("tespy.models.template.single_fluid", return_value=None):
            assert not self.model.supports_fluid_diagram("c1")

    def test_supports_fluid_diagram_pure_fluid(self):
        self.model.register_diagram("water", MagicMock())
        with patch("tespy.models.template.single_fluid", return_value="water"):
            assert self.model.supports_fluid_diagram("c1")

    def test_supports_fluid_diagram_failing_creation(self):
        self.model.create_diagram = MagicMock(side_effect=ValueError("nope"))
        with patch("tespy.models.template.single_fluid", return_value="water"):
            assert not self.model.supports_fluid_diagram("c1")

    def test_mixture_prepare_raises_when_strict(self):
        with patch("tespy.models.template.single_fluid", return_value=None):
            with pytest.raises(ValueError, match="mixture"):
                self.model._prepare_diagram_and_process_data("c1", strict=True)

    def test_explicit_diagram_is_passed_through(self):
        diagram = MagicMock()
        diagram.calc_individual_isoline = MagicMock(return_value="isoline")
        plotting_data = ({"comp": {"some": "data"}}, {"a": {"s": 1.0}})
        with (
            patch("tespy.models.template.single_fluid", return_value="water"),
            patch("tespy.models.template.get_plotting_data", return_value=plotting_data)
        ):
            processes, points, used = (
                self.model._prepare_diagram_and_process_data("c1", diagram=diagram)
            )
        assert used is diagram
        assert processes == {"comp": "isoline"}

    def _plot_points_only(self, method):
        import matplotlib
        matplotlib.use("Agg")
        points = {
            "a": {"p": 1.0, "h": 100.0, "T": 25.0, "s": 1.0, "vol": 0.5},
            "b": {"p": 5.0, "h": 200.0, "T": 50.0, "s": 2.0, "vol": 0.1},
        }
        self.model.nw.units.get_default = MagicMock(return_value="unit")
        with (
            patch("tespy.models.template.single_fluid", return_value=None),
            patch("tespy.models.template.get_plotting_data", return_value=({}, points))
        ):
            fig, ax = method("c1")
        return fig, ax

    def test_ts_plot_falls_back_to_points_for_mixture(self, caplog):
        import logging
        with caplog.at_level(logging.WARNING):
            fig, ax = self._plot_points_only(
                self.model.plot_Ts_diagram_matplotlib
            )
        assert "mixture" in caplog.text
        # the two state points must be scattered onto the axes
        assert len(ax.collections) == 2

    def test_logph_plot_falls_back_to_points_for_mixture(self):
        fig, ax = self._plot_points_only(
            self.model.plot_logph_diagram_matplotlib
        )
        assert ax.get_yscale() == "log"
        assert len(ax.collections) == 2
