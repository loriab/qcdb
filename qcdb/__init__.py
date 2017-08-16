#
# @BEGIN LICENSE
#
# QCDB: quantum chemistry common driver and databases
#
# Copyright (c) 2011-2017 The QCDB Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of QCDB.
#
# QCDB is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# QCDB is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with QCDB; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

"""Module to facilitate quantum chemical computations on chemical
databases. Contains Molecule class and physical constants from psi4 suite.

"""
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
__version__ = '0.4'
__author__ = 'Lori A. Burns'

# Load Python modules
import sys
from .molecule import Molecule
from .dbproc import *
from .options import *
from .qcformat import *
from . import cfour
from . import jajo
from . import orca
from .orient import OrientMols
from .dbwrap import Database, DB4 #DatabaseWrapper  #ReactionDatum, Reagent, Reaction
from .libmintspointgrp import SymmetryOperation, PointGroup
from .libmintsbasisset import BasisSet

# Load items that are useful to access from an input file
from .psiutil import *
from .physconst import *
