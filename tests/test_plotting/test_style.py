import pytest

from tespy.plotting import settings
from tespy.plotting import style as style_module
from tespy.plotting.style import Style
from tespy.plotting.style import available_styles
from tespy.plotting.style import get_style
from tespy.plotting.style import register_style


@pytest.fixture(autouse=True)
def _isolated_registry(monkeypatch):
    monkeypatch.setattr(style_module, "_REGISTRY", dict(style_module._REGISTRY))
    settings.reset_defaults()
    yield
    settings.reset_defaults()


def _luminance(color):
    channels = [int(color[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    linear = [
        c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
        for c in channels
    ]
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]


def test_shipped_styles_are_registered():
    assert {"tespy", "tespy-paper", "tespy-dark"} <= set(available_styles())


def test_default_style_is_tespy():
    assert get_style().name == "tespy"
    assert get_style(None) is get_style("tespy")


def test_style_object_is_passed_through():
    custom = get_style("tespy").derive("mine", cycle="#000000")
    assert get_style(custom) is custom


def test_unknown_style_name_raises():
    with pytest.raises(KeyError, match="tespy-dark"):
        get_style("nonexistent")


def test_dark_style_uses_dark_logo_and_surface():
    dark = get_style("tespy-dark")
    assert dark.logo == "dark"
    assert _luminance(dark.surface) < 0.05
    assert _luminance(get_style("tespy").surface) > 0.9


@pytest.mark.parametrize("name", ["tespy", "tespy-paper", "tespy-dark"])
def test_sequential_colormap_runs_away_from_the_surface(name):
    style = get_style(name)
    surface = _luminance(style.surface)
    distances = [abs(_luminance(c) - surface) for c in style.sequential]
    assert distances == sorted(distances)


@pytest.mark.parametrize("name", ["tespy", "tespy-paper", "tespy-dark"])
def test_diverging_colormap_has_neutral_midpoint(name):
    colors = get_style(name).diverging
    assert len(colors) % 2 == 1
    middle = colors[len(colors) // 2]
    red, green, blue = (int(middle[i:i + 2], 16) for i in (1, 3, 5))
    assert max(red, green, blue) - min(red, green, blue) <= 12


def test_slots_cycle_through_colors_markers_and_dashes():
    style = get_style("tespy")
    n = len(style.categorical)
    assert style.color(n) == style.color(0)
    assert style.marker(len(style.markers)) == style.marker(0)
    assert style.dash(len(style.dashes)) == style.dash(0)


def test_paper_style_distinguishes_slots_by_dash():
    style = get_style("tespy-paper")
    assert len(set(style.dashes)) >= 4


def test_derive_keeps_other_roles():
    base = get_style("tespy")
    derived = base.derive("mine", hot="#ff0000")
    assert derived.name == "mine"
    assert derived.hot == "#ff0000"
    assert derived.cold == base.cold
    assert base.hot != "#ff0000"


def test_colors_are_normalized_to_hex():
    style = get_style("tespy").derive("named", hot="red", cold="#00F")
    assert style.hot == "#ff0000"
    assert style.cold == "#0000ff"


def test_invalid_color_raises():
    with pytest.raises(ValueError):
        get_style("tespy").derive("broken", hot="not-a-color")


def test_empty_categorical_raises():
    with pytest.raises(ValueError, match="categorical"):
        get_style("tespy").derive("broken", categorical=())


def test_unknown_marker_raises():
    with pytest.raises(ValueError, match="marker"):
        get_style("tespy").derive("broken", markers=("star-of-david",))


def test_register_and_use_custom_style():
    custom = get_style("tespy").derive("company", cycle="#00395b")
    register_style(custom)
    assert get_style("company") is custom
    settings.set_defaults(style="company")
    assert get_style().name == "company"


def test_register_existing_name_needs_overwrite():
    custom = get_style("tespy").derive("tespy")
    with pytest.raises(ValueError, match="overwrite"):
        register_style(custom)
    register_style(custom, overwrite=True)
    assert get_style("tespy") is custom


def test_dict_round_trip():
    style = get_style("tespy-dark")
    assert Style.from_dict(style.to_dict()) == style


def test_from_dict_based_on_shipped_style():
    style = Style.from_dict({
        "name": "company", "based_on": "tespy",
        "categorical": ["#00395B", "#B54036", "#74ADC0", "#EC6707"]
    })
    assert style.categorical[0] == "#00395b"
    assert style.surface == get_style("tespy").surface


def test_from_toml(tmp_path):
    path = tmp_path / "style.toml"
    path.write_text(
        'name = "company"\n'
        'based_on = "tespy-dark"\n'
        'hot = "#B54036"\n'
    )
    style = Style.from_toml(path)
    assert style.name == "company"
    assert style.hot == "#b54036"
    assert style.surface == get_style("tespy-dark").surface


def test_rc_params_follow_style():
    style = get_style("tespy-dark")
    params = style.rc_params()
    assert params["axes.facecolor"] == style.surface
    assert params["figure.facecolor"] == style.surface
    assert params["text.color"] == style.ink


def test_matplotlib_colormap():
    from matplotlib.colors import Colormap
    style = get_style("tespy")
    cmap = style.colormap()
    assert isinstance(cmap, Colormap)
    assert cmap.name == "tespy"
    assert style.colormap("diverging").name == "tespy_diverging"
    with pytest.raises(ValueError, match="sequential"):
        style.colormap("rainbow")


def test_plotly_colorscale():
    scale = get_style("tespy").colorscale("diverging")
    assert scale[0][0] == 0
    assert scale[-1][0] == 1
    assert [color for _, color in scale] == list(get_style("tespy").diverging)


def test_defaults_can_be_changed_and_reset():
    settings.set_defaults(logo=False, logo_location="upper left")
    assert settings.get_defaults()["logo"] is False
    assert settings.get_defaults()["logo_location"] == "upper left"
    settings.reset_defaults()
    assert settings.get_defaults() == {
        "style": "tespy", "logo": True, "logo_location": "lower right"
    }


def test_invalid_default_location_raises():
    with pytest.raises(ValueError, match="lower right"):
        settings.set_defaults(logo_location="middle")


def test_default_style_must_be_registered():
    with pytest.raises(KeyError):
        settings.set_defaults(style="nonexistent")
