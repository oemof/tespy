# -*- coding: utf-8

"""Session defaults of the plotting module.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/plotting/settings.py

SPDX-License-Identifier: MIT
"""

#: Corners for the logo and the information block.
LOGO_LOCATIONS = ("upper left", "upper right", "lower left", "lower right")

_FACTORY_DEFAULTS = {
    "style": "tespy", "logo": True, "logo_location": "lower right"
}
_DEFAULTS = dict(_FACTORY_DEFAULTS)


def check_location(location):
    if location not in LOGO_LOCATIONS:
        msg = (
            f"The logo location {location} is unknown, use one of "
            f"{', '.join(LOGO_LOCATIONS)}."
        )
        raise ValueError(msg)


def set_defaults(style=None, logo=None, logo_location=None):
    """Change the defaults of all plots of this session.

    Parameters
    ----------
    style : str, optional
        Name of a registered style.

    logo : bool, optional
        Show the TESPy logo and the information block.

    logo_location : str, optional
        Corner of the logo, see :py:data:`LOGO_LOCATIONS`.

    Example
    -------
    >>> from tespy.plotting.settings import get_defaults, set_defaults
    >>> set_defaults(style="tespy-dark", logo=False)
    >>> get_defaults()["style"]
    'tespy-dark'
    >>> set_defaults(style="tespy", logo=True)
    """
    if style is not None:
        from tespy.plotting.style import get_style
        _DEFAULTS["style"] = get_style(style).name
    if logo is not None:
        _DEFAULTS["logo"] = bool(logo)
    if logo_location is not None:
        check_location(logo_location)
        _DEFAULTS["logo_location"] = logo_location


def get_defaults():
    """Return a copy of the current defaults."""
    return dict(_DEFAULTS)


def reset_defaults():
    """Restore the defaults TESPy ships with."""
    _DEFAULTS.clear()
    _DEFAULTS.update(_FACTORY_DEFAULTS)


def resolve(logo=None, logo_location=None):
    """Logo switch and location of a call, the defaults where not given."""
    logo = _DEFAULTS["logo"] if logo is None else logo
    location = _DEFAULTS["logo_location"] if logo_location is None else logo_location
    check_location(location)
    return logo, location
