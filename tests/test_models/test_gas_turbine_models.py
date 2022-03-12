# -*- coding: utf-8 -*-

"""Module for testing two tespy simulation against each other.

This file is part of project TESPy (github.com/oemof/tespy). It's copyrighted
by the contributors recorded in the version control history of the file,
available from its original location
tests/test_models/test_gasturbine_model.py

SPDX-License-Identifier: MIT
"""
import shutil

from tespy.components import CombustionChamber
from tespy.components import Compressor
from tespy.components import Sink
from tespy.components import Source
from tespy.components import Turbine
from tespy.connections import Connection
from tespy.networks import Network


class TestGasturbine:

    def __init__(self) -> None:
        pass
