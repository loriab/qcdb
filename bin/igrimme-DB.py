#!/usr/bin/python

import sys
import math
import os
import re
import string
qcdbpkg_path = os.path.expanduser('~loriab/linux/qcdb')
sys.path.append(qcdbpkg_path)
import qcdb
import qcdb.basislist
sys.path.append(qcdbpkg_path + '/databases')
import qcdb.exceptions


# load docstring info from database files (doesn't actually import database modules)
DBdocstrings = qcdb.dictify_database_docstrings()

# instructions
print """
 Welcome to igrimme-db.
    Just fill in the variables when prompted.
    Hit ENTER to accept default.
    Strings should not be in quotes.
    Elements in arrays should be space-delimited.
    Nothing is case sensitive.
"""

# query database name
module_choices = dict(zip([x.upper() for x in DBdocstrings.keys()], DBdocstrings.keys()))

print 'n Choose your database.'
for item in module_choices.keys():
    print '    %-12s   %s' % ('[' + module_choices[item] + ']', DBdocstrings[module_choices[item]]['general'][0].lstrip(" |"))
print '\n',

user_obedient = False
while not user_obedient:
    temp = raw_input('    dbse = ').strip()
    if temp.upper() in module_choices.keys():
        db_name = module_choices[temp.upper()]
        user_obedient = True

# query functionals
print '\n State your functional with Psi4 names (multiple allowed).\n'

functionals = []
user_obedient = False
while not user_obedient:
    temp = raw_input('    functionals [B3LYP PBE PBE0 BP86 B2PLYP B97-D] = ').strip()
    ltemp = temp.split()
    if temp == "":
        functionals = ['B3LYP', 'PBE', 'PBE0', 'BP86', 'B2PLYP', 'B97-D']
        user_obedient = True
    for item in ltemp:
        functionals.append(item.upper())
        user_obedient = True

# query dash level
print """
 Choose your -D correction level (multiple allowed).
    [d2]
    [d3zero]
    [d3bj]
"""

variants = []
user_obedient = False
while not user_obedient:
    temp = raw_input('    variants [d2 d3zero d3bj] = ').strip()
    ltemp = temp.split()
    if temp == "":
        variants = ['d2', 'd3zero', 'd3bj']
        user_obedient = True
    for item in ltemp:
        if item.lower() in ['d2', 'd3zero', 'd3bj']:
            variants.append(item.lower())
            user_obedient = True
        else:
            user_obedient = False
            variants = []
            break


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
                   functionals = %s
                      variants = %s

        <<< SCANNED SETTINGS  DISREGARD RESULTS IF INAPPROPRIATE  SCANNED SETTINGS >>>

""" % (dbse, functionals, variants)


# commence main computation loop
maxrgt = max([len(ACTV[rgt]) for rgt in ACTV])

h2 = qcdb.Molecule("""\nH\nH 1 0.74\n""")
h2.update_geometry()

for func in functionals:
    for variant in variants:
        print """\n  ==> %s-%s <==""" % (func.upper(), variant.upper())

        # catch if this is an invalid functional-dashlvl before run database
        try:
            h2.grimme_dash(func=func, dashlvl=variant)
        except qcdb.ValidationError:
            print 'No output or files will be generated.'
            continue

        #ddir = 'output-dftd3_$set-$functional-$variant'
        #try:
        #    os.mkdir(ddir)
        #except OSError:
        #    print 'Warning: directory %s already present.' % (ddir)

        dfile = dbse + '-' + func + '-' + variant + '.usemedash'
        usemedash = open(dfile, 'w')

        # print table header
        textline = '\n                                 '
        for i in range(maxrgt):
            textline += '             REAGENT_%s' % (string.uppercase[i])
        textline += '            REACTION\n                              SYS'
        for i in range(maxrgt):
            textline += '          DASHCORR  WT'
        textline += '       DISP      REF\n'
        print textline

        # commence iteration through reactions
        for rxn in HRXN:
            DASHCORR = {}
            index = dbse + '-' + str(rxn)

            # compute correction for rgts and rxn
            for rgt in ACTV[index]:
                GEOS[rgt].update_geometry()
                DASHCORR[rgt], temp = GEOS[rgt].grimme_dash(func=func, dashlvl=variant)
            IEDASH = sum([RXNM[index][rgt] * DASHCORR[rgt] for rgt in ACTV[index]]) * qcdb.psi_hartree2kcalmol

            # print main results
            textline = '%33s ' % (index)
            for i in range(maxrgt):
                if i < len(ACTV[index]):
                    textline += '|  %14.8f %3d ' % (DASHCORR[ACTV[index][i]], RXNM[index][ACTV[index][i]])
                else:
                    textline += '|                     '
            textline += '| %8.3f %8.3f' % (IEDASH, BIND[index])
            print textline

            # print main results to useme
            textline = '%-23s ' % (index)
            for rgt in ACTV[index]:
                textline += '%16.8f %4d ' % (DASHCORR[rgt], RXNM[index][rgt])
            textline += '\n'
            usemedash.write(textline)

        # print useme footer
        textline = '%-23s ' % ('#__elecE_in_hartree')
        for i in range(maxrgt):
            textline += '%16s %4s ' % ('Reagent' + string.uppercase[i], 'Wt')
        textline += '\n'
        usemedash.write(textline)
