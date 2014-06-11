#!/usr/bin/env python

#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

import os
#import re
import sys
import glob
#import math
import string
import importlib
import collections

qcdbpkg_path = os.path.dirname(__file__)
sys.path.append(qcdbpkg_path + '/../')
import qcdb
import qcdb.basislist
from qcdb.exceptions import *
sys.path.append(qcdbpkg_path + '/../databases')


# instructions
print """
 Welcome to herd-db.
#    Just execute in a directory of output files.
"""

actionable_data = {
    'SCF TOTAL ENERGY': '.usemeraw',
    'MP2 CORRELATION ENERGY': '.mp2.usemecorl',
    'MP3 CORRELATION ENERGY': '.mp3.usemecorl',
    'MP4 CORRELATION ENERGY': '.mp4.usemecorl',
    'CCSD CORRELATION ENERGY': '.ccsd.usemecorl',
    'CCSD(T) CORRELATION ENERGY': '.ccsdt.usemecorl',
    'CCSDT CORRELATION ENERGY': '.ccsdfullt.usemecorl',
    'CCSDT(Q) CORRELATION ENERGY': '.ccsdtq.usemecorl',
    }


def identify_qcprog(filename):
    """

    """
    qcprogs = {
        'PSI4: An Open-Source Ab Initio Electronic Structure Package': 'psi4',
        }

    with open(sample, 'r') as handle:
        contents = handle.readlines()
    for line in contents:
        for target in qcprogs.keys():
            if target in line:
                return qcprogs[target]

# query database, qcprog, and directory name
sample = glob.glob('*.out')[0]
db_name = sample.split('-')[0]
qcprog = identify_qcprog(sample)
dirprefix = os.path.split(os.getcwd())[-1]
if len(glob.glob('*-CP.out')) > 0:
    mode = 'ACTV_CP'
elif len(glob.glob('*-unCP.out')) > 0:
    mode = 'ACTV'
else:
    mode = 'ACTV_SA'

# Load module for QC program
try:
    qcmod = importlib.import_module('qcdb.' + qcprog)
except ImportError:
    print('\nPython module for QC program %s failed to load\n\n' % (qcprog))
    print('\nSearch path that was tried:\n')
    print(", ".join(map(str, sys.path)))
    raise ValidationError("Python module loading problem for QC program " + str(qcprog))
else:
    print qcmod
    qcmtdIN = qcmod.qcmtdIN

# Load module for requested database
try:
    database = __import__(db_name)
except ImportError:
    print('\nPython module for database %s failed to load\n\n' % (db_name))
    print('\nSearch path that was tried:\n')
    print(", ".join(map(str, sys.path)))
    raise ValidationError("Python module loading problem for database " + str(db_name))
else:
    dbse = database.dbse
    HRXN = database.HRXN
    #ACTV = database.ACTV
    ACTV = getattr(database, mode)
    RXNM = database.RXNM  # TODO alter RXNM too?
    BIND = database.BIND
    TAGL = database.TAGL
    GEOS = database.GEOS
    try:
        DATA = database.DATA
    except AttributeError:
        DATA = {}


print """
        <<< SCANNED SETTINGS  SCANNED SETTINGS  SCANNED SETTINGS  SCANNED SETTINGS >>>

                          dbse = %s
                        qcprog = %s
                     dirprefix = %s

        <<< SCANNED SETTINGS  DISREGARD RESULTS IF INAPPROPRIATE  SCANNED SETTINGS >>>

""" % (dbse, qcprog, dirprefix)

# commence iteration through reactions
psivar = collections.defaultdict(dict)
for rxn in HRXN:
    DASHCORR = {}
    index = dbse + '-' + str(rxn)
    textline = '%33s ' % (index)

    # collect interesting results for rgts and rxn
    for rgt in ACTV[index]:
        complete = False
        rxnm_wt = RXNM[index][ACTV[index][ACTV[index].index(rgt)]]
        try:
            with open(rgt + '.out', 'r') as handle:
                contents = handle.read()
        except IOError:
            textline += """|  %14s %3d """ % ('<<< MIA >>>', rxnm_wt)
            continue
        pv, _tmp, _tmp2 = qcmod.harvest_output(contents)
        for key, val in pv.items():
            psivar[key][rgt] = val
            if key == 'SUCCESS':
                complete = True
        textline += """|  %14s %3d """ % ('Completed' if complete else 'Running', rxnm_wt)
    print textline

# prepare useme footer
footer = '%-23s ' % ('#__elecE_in_hartree')
for i in range(max([len(ACTV[rgt]) for rgt in ACTV])):
    #footer += '%16s %4s ' % ('Reagent' + string.uppercase[i], 'Wt')
    footer += '%16s ' % ('Reagent' + string.uppercase[i])
footer += '\n'

# print main results to useme
for datum in psivar.keys():
    usemecontents = ''
    try:
        usemeext = actionable_data[datum]
    except KeyError:
        continue

    for rxn in HRXN:
        index = dbse + '-' + str(rxn)
        textline = '%-23s ' % (index)

        for rgt in ACTV[index]:
            rxnm_wt = RXNM[index][ACTV[index][ACTV[index].index(rgt)]]
            try:
                #textline += '%16.8f %4d ' % (psivar[datum][rgt], rxnm_wt)
                textline += '%16.8f ' % (psivar[datum][rgt])
            except KeyError:
                break
        else:
            usemecontents += textline + '\n'
            continue

    if len(usemecontents) > 0:
        with open(dirprefix + usemeext, 'w') as handle:
            handle.write(usemecontents)
            handle.write(footer)
