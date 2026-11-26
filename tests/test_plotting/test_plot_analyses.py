import math

import pandas as pd
import pytest

from tespy.plotting import operating_map_data
from tespy.plotting import plot_operating_map
from tespy.plotting import plot_operating_map_plotly
from tespy.plotting import plot_sweep
from tespy.plotting import plot_sweep_plotly
from tespy.plotting import sweep_data

SWEEP = {
    "T_cond": [35, 45, 55, 35, 45, 55],
    "COP": [6.0, 4.6, 3.7, 6.1, 4.8, 3.9],
    "fluid": ["R290", "R290", "R290", "R600a", "R600a", "R600a"],
}
MAP = {
    "T_cond": [35, 45, 35, 45, 55],
    "T_evap": [0, 0, 10, 10, 10],
    "COP": [5.0, 4.0, 6.0, 5.0, 4.0],
}


def test_sweep_one_series_per_group():
    description = sweep_data(SWEEP, "T_cond", "COP", group="fluid", units={"T_cond": "°C"})
    assert [s.name for s in description.series] == ["R290", "R600a"]
    assert [s.slot for s in description.series] == [0, 1]
    assert list(description.series[1].y) == [6.1, 4.8, 3.9]
    assert description.x_axis.unit == "°C"


def test_sweep_from_dataframe_without_group():
    description = sweep_data(pd.DataFrame(SWEEP), "T_cond", "COP")
    assert len(description.series) == 1


def test_operating_map_grid():
    description = operating_map_data(MAP, "T_cond", "T_evap", "COP")
    grid = description.series[0]
    assert list(grid.x) == [35, 45, 55]
    assert list(grid.y) == [0, 10]
    assert grid.z[1] == [6.0, 5.0, 4.0]
    assert math.isnan(grid.z[0][2])
    assert description.z_axis.label == "COP"


def test_table_plots_with_both_renderers():
    assert plot_sweep(SWEEP, "T_cond", "COP", group="fluid").has_data()
    assert plot_operating_map(MAP, "T_cond", "T_evap", "COP").has_data()
    assert plot_sweep_plotly(SWEEP, "T_cond", "COP").data[0].type == "scatter"
    assert plot_operating_map_plotly(MAP, "T_cond", "T_evap", "COP").data[0].type == "heatmap"
