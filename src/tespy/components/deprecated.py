# -*- coding: utf-8
import logging

from tespy.tools.helpers import TESPyComponentError

# %%


class PlaceHolderForError:
    r"""Throw errors for deprecated classes."""

    def __init__(self, *args, **kwargs):
        msg = (
            'Version 0.4.x introduced PEP 8 conform class names for all '
            'classes used in TESPy. Please refer to '
            'https://tespy.readthedocs.io/en/main/whats_new.html '
            'learn about the changes necessary to adapt your script to the '
            'latest API.'
        )
        logging.error(msg)
        raise TESPyComponentError(msg)


class cycle_closer(PlaceHolderForError):
    pass


class sink(PlaceHolderForError):
    pass


class source(PlaceHolderForError):
    pass


class subsystem_interface(PlaceHolderForError):
    pass


class combustion_chamber(PlaceHolderForError):
    pass


class combustion_chamber_stoich(PlaceHolderForError):
    pass


class combustion_engine(PlaceHolderForError):
    pass


class orc_evaporator(PlaceHolderForError):
    pass


class condenser(PlaceHolderForError):
    pass


class desuperheater(PlaceHolderForError):
    pass


class heat_exchanger(PlaceHolderForError):
    pass


class heat_exchanger_simple(PlaceHolderForError):
    pass


class parabolic_trough(PlaceHolderForError):
    pass


class solar_collector(PlaceHolderForError):
    pass


class droplet_separator(PlaceHolderForError):
    pass


class drum(PlaceHolderForError):
    pass


class merge(PlaceHolderForError):
    pass


class node(PlaceHolderForError):
    pass


class separator(PlaceHolderForError):
    pass


class splitter(PlaceHolderForError):
    pass


class pipe(PlaceHolderForError):
    pass


class valve(PlaceHolderForError):
    pass


class water_electrolyzer(PlaceHolderForError):
    pass


class compressor(PlaceHolderForError):
    pass


class pump(PlaceHolderForError):
    pass


class turbine(PlaceHolderForError):
    pass


class subsystem(PlaceHolderForError):
    pass
