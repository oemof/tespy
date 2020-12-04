# -*- coding: utf-8

"""Module of class Ref.


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/connections/ref.py

SPDX-License-Identifier: MIT
"""

import logging
import warnings

import numpy as np
import pandas as pd

from tespy.components.component import Component
from tespy.tools.characteristics import CharLine
from tespy.tools.data_containers import dc_cp
from tespy.tools.data_containers import dc_flu
from tespy.tools.data_containers import dc_prop
from tespy.tools.data_containers import dc_simple
from tespy.tools.helpers import TESPyConnectionError

# pass the warning messages to the logger
logging.captureWarnings(True)


class Ref:
    r"""
    Reference fluid properties from one connection to another connection.

    Parameters
    ----------
    obj : tespy.connections.connection
        Connection to be referenced.

    f : float
        Factor to multiply specified property with.

    d : float
        Delta to add after multiplication.

    Note
    ----
    Reference the mass flow of one connection :math:`\dot{m}` to another mass
    flow :math:`\dot{m}_{ref}`

    .. math::

        \dot{m} = \dot{m}_\mathrm{ref} \cdot f + d

    """

    def __init__(self, ref_obj, factor, delta):

        if not isinstance(ref_obj, Connection):
            msg = 'First parameter must be object of type connection.'
            logging.error(msg)
            raise TypeError(msg)

        if not (isinstance(factor, int) or isinstance(factor, float)):
            msg = 'Second parameter must be of type int or float.'
            logging.error(msg)
            raise TypeError(msg)

        if not (isinstance(delta, int) or isinstance(delta, float)):
            msg = 'Thrid parameter must be of type int or float.'
            logging.error(msg)
            raise TypeError(msg)

        self.obj = ref_obj
        self.f = factor
        self.d = delta

        msg = ('Created reference object with factor ' + str(self.f) +
               ' and delta ' + str(self.d) + ' referring to connection ' +
               ref_obj.source.label + ' (' + ref_obj.source_id + ') -> ' +
               ref_obj.target.label + ' (' + ref_obj.target_id + ').')
        logging.debug(msg)

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
            logging.error(msg)
            raise KeyError(msg)
