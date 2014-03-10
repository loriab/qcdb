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

import sys
import math
import os
import re
import importlib
import collections

qcdbpkg_path = os.path.dirname(__file__)
sys.path.append(qcdbpkg_path + '/../')
import qcdb
import qcdb.basislist
from qcdb.exceptions import *
sys.path.append(qcdbpkg_path + '/../databases')

# load docstring info from database files (doesn't actually import database modules)
DBdocstrings = qcdb.dictify_database_docstrings()

# instructions
print """
 Welcome to imake-db.
    Just fill in the variables when prompted.
    Hit ENTER to accept default.
    Strings should not be in quotes.
    Elements in arrays should be space-delimited.
    Nothing is case sensitive.
"""

# query database name
module_choices = dict(zip([x.upper() for x in DBdocstrings.keys()], DBdocstrings.keys()))

print '\n Choose your database.'
for item in module_choices.keys():
    print '    %-12s   %s' % ('[' + module_choices[item] + ']', DBdocstrings[module_choices[item]]['general'][0].lstrip(" |"))
print '\n'

user_obedient = False
while not user_obedient:
    temp = raw_input('    dbse = ').strip()
    if temp.upper() in module_choices.keys():
        db_name = module_choices[temp.upper()]
        user_obedient = True

# query database subset
subset_choices = dict(zip([x.upper() for x in DBdocstrings[db_name]['subset'].keys()], DBdocstrings[db_name]['subset'].keys()))

print '\n Choose your subset (multiple allowed).'
for key, val in DBdocstrings[db_name]['subset'].items():
    print '    %-12s   %s' % ('[' + key + ']', val)
print '\n'

subset = []
user_obedient = False
while not user_obedient:
    temp = raw_input('    subset [all] = ').strip()
    ltemp = temp.split()
    if temp == "":
        user_obedient = True
    for item in ltemp:
        if item.upper() in subset_choices.keys():
            subset.append(subset_choices[item.upper()])
            user_obedient = True
        else:
            user_obedient = False
            subset = []
            break

# query qc program
print """
 Choose your quantum chemistry program.
    #[qchem]
    [molpro]       writes Molpro input files
    [molpro2]      writes Molpro input files
    [psi4]         writes Psi4 input files
    #[nwchem]
    #[xyz]          writes basic xyz files only
"""

user_obedient = False
while not user_obedient:
    temp = raw_input('    qcprog = ').strip()
    if temp.lower() in ['molpro', 'psi4', 'molpro2']:
        qcprog = temp.lower()
        user_obedient = True

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

# query quantum chemical method(s)
method_choices = dict(zip([x.upper() for x in qcmtdIN.keys()], qcmtdIN.keys()))

print '\n Choose your quantum chemical methods (multiple allowed).'
for key, val in qcmtdIN.items():
    print '    %-12s' % ('[' + key + ']')
print '\n'

methods = []
user_obedient = False
while not user_obedient:
    temp = raw_input('    methods = ').strip()
    ltemp = temp.split()
    for item in ltemp:
        if item.upper() in method_choices:
            methods.append(method_choices[item.upper()])
            user_obedient = True
        else:
            user_obedient = False
            methods = []
            break

# set up options dict
options = collections.defaultdict(lambda: collections.defaultdict(dict))

# query basis set(s)
print """
 Choose your basis set (multiple allowed).
    e.g., aug-cc-pvdz or 6-31+G* or cc-pvtz may-cc-pvtz aug-cc-pvtz
"""

bases = []
user_obedient = False
while not user_obedient:
    temp = raw_input('    bases = ').strip()
    ltemp = temp.split()
    for item in ltemp:
        btemp = qcdb.basislist.corresponding_orbital(item)
        if btemp:
            bases.append(btemp)
            user_obedient = True
        else:
            print '    Basis set %s not recognized.' % (item)
            proceed = qcdb.query_yes_no('    Proceed anyway? =', False)
            if proceed:
                bases.append(item)
                user_obedient = True
            else:
                bases = []
                user_obedient = False
                break

# query castup preference
print """
 Do cast up from smaller basis set?
"""

user_obedient = False
while not user_obedient:
    castup = qcdb.query_yes_no('    castup [F] = ', False)
    user_obedient = True
options['SCF']['BASIS_GUESS']['value'] = castup

# query directory prefix
print """
 State your destination directory prefix.
"""

user_obedient = False
while not user_obedient:
    temp = raw_input('    dirprefix [try] = ').strip()
    if temp == "":
        dirprefix = 'try'
        user_obedient = True
    if temp.isalnum():
        dirprefix = temp
        user_obedient = True

# query memory
print """
 Choose your memory usage in MB.
"""

user_obedient = False
while not user_obedient:
    temp = raw_input('    memory [1600] = ').strip()
    if temp == "":
        memory = 1600
        user_obedient = True
    if temp.isdigit():
        memory = int(temp)
        user_obedient = True

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
    ACTV = database.ACTV
    RXNM = database.RXNM
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
                        subset = %s
                        qcprog = %s
                       methods = %s
                         bases = %s
                     dirprefix = %s
                   memory [MB] = %d
                       usesymm =
                       cast up = %s

        <<< SCANNED SETTINGS  DISREGARD RESULTS IF INAPPROPRIATE  SCANNED SETTINGS >>>

""" % (dbse, subset, qcprog, methods, bases, dirprefix, memory, castup)


# establish multiplicity hash table
mult = {
   1: "singlet",
   2: "doublet",
   3: "triplet",
   4: "quartet",
   5: "quintet",
   6: "sextet",
   7: "septet",
   8: "octet"}


# file extension
fext = 'xyz' if qcprog == 'xyz' else 'in'

# merge and condense HRXN from subset
if len(subset) == 0:
    pass
else:
    temp = []
    for item in subset:
        if item == 'small':
            temp.append(database.HRXN_SM)
        elif item == 'large':
            temp.append(database.HRXN_LG)
        elif item == 'equilibrium':
            temp.append(database.HRXN_EQ)
        else:
            try:
                temp.append(getattr(database, item))
            except AttributeError:
                try:
                    temp.append(getattr(database, 'HRXN_' + item))
                except AttributeError:
                    raise ValidationError('Special subset \'%s\' not available for database %s.' % (item, db_name))
    HRXN = qcdb.drop_duplicates(temp)

# assemble reagent list from reaction list
temp = []
for rxn in HRXN:
    temp.append(database.ACTV['%s-%s' % (dbse, rxn)])
    temp.append(database.ACTV_CP['%s-%s' % (dbse, rxn)])
HSYS = qcdb.drop_duplicates(temp)

# commence the file-writing loop
tdir = '-'.join([dirprefix, dbse, qcprog])
try:
    os.mkdir(tdir)
except OSError:
    print 'Warning: directory %s already present.' % (tdir)

for basis in bases:
    options['GLOBALS']['BASIS']['value'] = basis
    basdir = qcdb.basislist.sanitize_basisname(basis)
    basdir = re.sub('-', '', basdir)

    for method in methods:
        mtddir = qcdb.basislist.sanitize_basisname(method).upper()
        mtddir = re.sub('-', '', mtddir)

        subdir = '-'.join([basdir, mtddir])
        try:
            os.mkdir(tdir + '/' + subdir)
        except OSError:
            print 'Warning: directory %s/%s already present.' % (tdir, subdir)

        # TODO: forcing c1 symm skipped - still needed for xdm and molpro

        for system in HSYS:
            # QC program may reorient but at least input file geometry will match database
            GEOS[system].fix_orientation(True)
            GEOS[system].PYmove_to_com = False
            GEOS[system].tagline = 'index %s label %s' % (system, TAGL[system])
            GEOS[system].update_geometry()

            dertype = 0

            try:
                if qcprog == 'molpro':
                    infile = qcmod.MolproIn(memory, method, basis, GEOS[system], system, castup).format_infile_string()

                elif qcprog == 'psi4' or qcprog == 'molpro2':
                    infile = qcmod.Infile(memory, GEOS[system], method, dertype, options).format_infile_string()

            except FragmentCountError:
                pass
                # We're passing ACTV rgt list for SAPT methods so this error is to be expected

            else:
                sfile = tdir + '/' + subdir + '/' + system + '.' + fext
                with open(sfile, 'w') as handle:
                    handle.write(infile)
