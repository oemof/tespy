import importlib.metadata
import os

import matplotlib
import pytest

matplotlib.use("Agg")

from matplotlib import pyplot as plt  # noqa: E402
from matplotlib.offsetbox import OffsetImage  # noqa: E402
from matplotlib.text import Text  # noqa: E402

from tespy.plotting.branding import add_logo  # noqa: E402
from tespy.plotting.branding import info_lines  # noqa: E402
from tespy.plotting.branding import logo_path  # noqa: E402


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.mark.parametrize("variant", ["light", "dark"])
@pytest.mark.parametrize("fmt", ["svg", "png"])
def test_logo_files_are_shipped(variant, fmt):
    path = logo_path(variant, fmt)
    assert os.path.isfile(path)
    assert path.endswith(f".{fmt}")


def test_unknown_logo_variant_raises():
    with pytest.raises(ValueError, match="light"):
        logo_path("blue")


def test_unknown_logo_format_raises():
    with pytest.raises(ValueError, match="png"):
        logo_path("light", "jpg")


def test_info_lines_name_tespy_version_first():
    version = importlib.metadata.version("tespy")
    assert info_lines() == [f"TESPy {version}"]
    assert info_lines(["heat pump", "COP 3.7"])[1:] == ["heat pump", "COP 3.7"]


def test_info_lines_accept_single_string():
    assert info_lines("heat pump")[1:] == ["heat pump"]


def test_add_logo_to_matplotlib_axes():
    _, ax = plt.subplots()
    add_logo(ax, lines=["heat pump"])
    assert len(ax.findobj(OffsetImage)) == 1
    texts = [t.get_text() for t in ax.findobj(Text)]
    assert "heat pump" in "\n".join(texts)


def test_add_logo_rejects_unknown_location():
    _, ax = plt.subplots()
    with pytest.raises(ValueError, match="lower right"):
        add_logo(ax, location="middle")
