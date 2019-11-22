# -*- coding: utf-8

"""This module contains the component: subsystem


This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location tespy/components/subsystems.py

SPDX-License-Identifier: MIT
"""

import numpy as np
import logging

from tespy.networks import networks

# %%


class subsystem:
    r"""
    Class subsystem is the base class of all TESPy subsystems.

    Parameters
    ----------
    label : str
        The label of the subsystem.

    **kwargs :
        See the class documentation of desired subsystem for available
        keywords. If you want to define your own subsystem, provide the
        keywords in the :func:`tespy.components.subsystems.subystem.attr`
        method.

    Note
    ----
    The initialisation method (__init__), setter method (set_attr) and getter
    method (get_attr) are used for instances of class subsystem and its
    children.

    Keywords for keyword arguments depend on the type of subsystem you want to
    create/use.

    Example
    -------
    Basic example for a setting up a tespy.components.subsystems.subsystem
    object.
    This example does not run a tespy calculation!

    >>> from tespy import subsys
    >>> mysub = subsys.subsystem('mySubsystem')
    >>> type(mysub)
    <class 'tespy.components.subsystems.subsystem'>
    >>> mysub.get_attr('label')
    'mySubsystem'
    >>> type(mysub.get_network())
    <class 'tespy.networks.network'>
    """
    def __init__(self, label, **kwargs):

        if not isinstance(label, str):
            msg = 'Subsystem label must be of type str!'
            logging.error(msg)
            raise ValueError(msg)

        elif len([x for x in [';', ', ', '.'] if x in label]) > 0:
            msg = 'Can\'t use ' + str([';', ', ', '.']) + ' in label.'
            logging.error(msg)
            raise ValueError(msg)
        else:
            self.label = label

        # set default values
        for key in self.attr():
            self.__dict__.update({key: np.nan})
            self.__dict__.update({key + '_set': False})

        self.subsys_init()
        self.set_attr(**kwargs)

    def set_attr(self, **kwargs):
        r"""
        Sets, resets or unsets attributes of a component for provided keyword
        arguments.

        Parameters
        ----------
        label : str
            The label of the subsystem.

        **kwargs :
            See the class documentation of desired subsystem for available
            keywords. If you want to define your own subsystem, provide the
            keywords in the :func:`tespy.components.subsystems.subystem.attr`
            method.

        Note
        ----
        Allowed keywords in kwargs are obtained from class documentation as all
        subsystems share the
        :func:`tespy.components.subsystems.subsystem.set_attr` method.
        """
        # set provided values,  check for invalid keys
        for key in kwargs:
            if key in self.attr():
                if (isinstance(kwargs[key], float) or
                        isinstance(kwargs[key], np.float64) or
                        isinstance(kwargs[key], int)):
                    if np.isnan(kwargs[key]):
                        self.__dict__.update({key: np.nan})
                        self.__dict__.update({key + '_set': False})
                    else:
                        self.__dict__.update({key: kwargs[key]})
                        self.__dict__.update({key + '_set': True})

                elif isinstance(kwargs[key], str):
                    self.__dict__.update({key: kwargs[key]})
            else:
                msg = ('Component ' + self.label + ' has no attribute ' +
                       str(key) + '.')
                logging.error(msg)
                raise KeyError(msg)

        self.set_comps()
        self.set_conns()

        msg = 'Created subsystem ' + self.label + '.'
        logging.debug(msg)

    def get_attr(self, key):
        r"""
        Get the value of a subsystem's attribute.

        Parameters
        ----------
        key : str
            The attribute you want to retrieve.

        Returns
        -------
        out :
            Value of specified attribute.
        """
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            msg = 'Subsystem ' + self.label + ' has no attribute ' + key + '.'
            logging.error(msg)
            raise KeyError(msg)

    def attr(self):
        return ['num_i', 'num_o']

    def subsys_init(self):
        r"""
        Creates network of the subsystem on subsystem creation.
        """
        self.create_comps()
        self.create_conns()

    def create_comps(self):
        return

    def create_conns(self):
        self.conns = []

    def set_comps(self):
        return

    def set_conns(self):
        if not hasattr(self, 'nw'):
            self.create_network()
        return

    def create_network(self):
        self.nw = networks.network(fluids=[])
        for c in self.conns:
            self.nw.add_conns(c)

    def get_network(self):
        return self.nw
