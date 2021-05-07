# -*- coding: utf-8
import logging

from tespy.tools.helpers import TESPyConnectionError

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
        raise TESPyConnectionError(msg)


class bus(PlaceHolderForError):
    pass


class connection(PlaceHolderForError):
    pass


class ref(PlaceHolderForError):
    pass
