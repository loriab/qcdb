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
import argparse
import importlib
import collections

qcdbpkg_path = os.path.dirname(__file__)
sys.path.append(qcdbpkg_path + '/../')
import qcdb
import qcdb.basislist
from qcdb.exceptions import *
sys.path.append(qcdbpkg_path + '/../databases')


# instructions
parser = argparse.ArgumentParser(description='Process quantum chemical results from output files.')
parser.add_argument('-d', '--dbmodule', help='force choice of database module')
parser.add_argument('-p', '--prefix', help='force prefix for usemefiles, defaults to dir name')
parser.add_argument('-q', '--qcprog', help='force choice of QC program parser')
parser.add_argument('-s', '--style', help='stype of usemefile (3col or Wt)')
parser.add_argument('-a', '--actv', help='force ACTV mode (cp, uncp, sapt), defaults to sample')
args = parser.parse_args()

actionable_data = {
    #'SCF TOTAL ENERGY': 'usemeraw'
    'HF TOTAL ENERGY': 'usemeraw',
    'MP2 CORRELATION ENERGY': 'mp2.usemecorl',
    'MP2 SAME-SPIN CORRELATION ENERGY': 'mp2.usemetrip',
    'MP3 CORRELATION ENERGY': 'mp3.usemecorl',
    'MP4 CORRELATION ENERGY': 'mp4.usemecorl',
    'CCSD CORRELATION ENERGY': 'ccsd.usemecorl',
    'CCSD SAME-SPIN CORRELATION ENERGY': 'ccsd.usemetrip',
    'CCSD(T) CORRELATION ENERGY': 'ccsdt.usemecorl',
    'CCSDT CORRELATION ENERGY': 'ccsdfullt.usemecorl',
    'CCSDT(Q) CORRELATION ENERGY': 'ccsdtq.usemecorl',
    
    'HF-CABS TOTAL ENERGY': 'f12.usemeraw',
    'MP2-F12 CORRELATION ENERGY': 'mp2f12.usemecorl',
    'MP2-F12 SAME-SPIN CORRELATION ENERGY': 'mp2f12.usemetrip',

    'CCSD-F12A CORRELATION ENERGY': 'ccsdaf12.usemecorl',
    'CCSD-F12A SAME-SPIN CORRELATION ENERGY': 'ccsdaf12.usemetrip',
    'CCSD-F12B CORRELATION ENERGY': 'ccsdbf12.usemecorl',
    'CCSD-F12B SAME-SPIN CORRELATION ENERGY': 'ccsdbf12.usemetrip',
    '(T)-F12AB CORRECTION ENERGY': 'ccsdnstabf12.usemecrct',
    'CCSD-F12C CORRELATION ENERGY': 'ccsdcf12.usemecorl',
    'CCSD-F12C SAME-SPIN CORRELATION ENERGY': 'ccsdcf12.usemetrip',
    '(T)-F12C CORRECTION ENERGY': 'ccsdnstcf12.usemecrct',

    'DFT FUNCTIONAL TOTAL ENERGY': 'DFT.usemeraw',
    'DISPERSION CORRECTION ENERGY': '-nobas.DFTdX.usemedash',
    'DOUBLE-HYBRID CORRECTION ENERGY': 'DHDFT.usemeraw',  # violation of conventions to get plain dhdft E!
}

def identify_qcprog(filename):
    """

    """
    qcprogs = {
        'PSI4: An Open-Source Ab Initio Electronic Structure Package': 'psi4',
        '***  PROGRAM SYSTEM MOLPRO  ***': 'molpro2',
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
if db_name == 'NBC1':  # TODO
    db_name = 'NBC10'
if db_name == 'HBC1':
    db_name = 'HBC6'
db_name = db_name if args.dbmodule is None else args.dbmodule
qcprog = identify_qcprog(sample) if args.qcprog is None else args.qcprog
dirprefix = os.path.split(os.getcwd())[-1] if args.prefix is None else args.prefix
usemeold = False if args.prefix is 'new' else True
if args.actv is None:
    if len(glob.glob('*-CP.out')) > 0:
        mode = 'ACTV_CP'
    elif len(glob.glob('*-unCP.out')) > 0:
        mode = 'ACTV'
    elif len(glob.glob('*reagent.out')) > 0:
        mode = 'ACTV'
    else:
        mode = 'ACTV_SA'
elif args.actv.upper() == 'UNCP':
    mode = 'ACTV'
else:
    mode = 'ACTV_' + args.actv.upper()

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
                      dbmodule = %s
                        qcprog = %s
                          actv = %s
                     dirprefix = %s

        <<< SCANNED SETTINGS  DISREGARD RESULTS IF INAPPROPRIATE  SCANNED SETTINGS >>>

""" % (dbse, db_name, qcprog, mode, dirprefix)

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
    if usemeold:
        footer += '%16s ' % ('Reagent' + string.uppercase[i])
    else:
        footer += '%16s %4s ' % ('Reagent' + string.uppercase[i], 'Wt')
footer += '\n'

isDHDFT = True if 'DOUBLE-HYBRID CORRECTION ENERGY' in psivar.keys() else False

# print main results to useme
print ''
for datum in psivar.keys():
    usemecontents = ''
    tally = 0
    try:
        usemeext = actionable_data[datum]
    except KeyError:
        continue

    if isDHDFT and usemeext in ['DFT.usemeraw', 'mp2.usemecorl']:  # sad hack
        continue

    for rxn in HRXN:
        index = dbse + '-' + str(rxn)
        textline = ''
        complete = True

        for rgt in ACTV[index]:
            rxnm_wt = RXNM[index][ACTV[index][ACTV[index].index(rgt)]]
            try:
                value = psivar[datum][rgt]
                if datum == 'DOUBLE-HYBRID CORRECTION ENERGY':
                    value += psivar['DFT FUNCTIONAL TOTAL ENERGY'][rgt]

                if usemeold:
                    textline += '%16.8f ' % (value)
                    if rxnm_wt == -2:
                        textline += '%16.8f ' % (value)
                else:
                    textline += '%16.8f %4d ' % (value, rxnm_wt)
            except KeyError:
                textline += '%16s ' % ('')
                complete = False

        if len(textline.strip()) > 0:
            if complete:
                usemecontents += '%-23s %s\n' % (index, textline)
                tally += 1
            else:
                usemecontents += '#%-22s %s\n' % (index, textline)

    if len(usemecontents) > 0:
        with open(dirprefix + '.' + usemeext, 'w') as handle:
            handle.write(usemecontents)
            handle.write(footer)
            print '        writing %4d entries to %s' % (tally, dirprefix + '.' + usemeext)

