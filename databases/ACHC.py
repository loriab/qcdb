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

"""
| Database of Adenine-Cytosine curves along all possible axes.
| Geometries from Parker et al. JCTC 11 4197 (2015).
| Reference interaction energies from Parker et al. JCTC 11 4197 (2015).


- **cp**  ``'off'`` <erase this comment and after unless on is a valid option> || ``'on'``

- **rlxd** ``'off'`` <erase this comment and after unless on is valid option> || ``'on'``


- **benchmark**

  - ``'<benchmark_name>'`` <Reference>.
  - |dl| ``'<default_benchmark_name>'`` |dr| <Reference>.

- **subset**

  - ``'small'`` <members_description>
  - ``'large'`` <members_description>
  - ``'<subset>'`` <members_description>

"""
import re
import qcdb

# <<< ACHC Database Module >>>
dbse = 'ACHC'

# <<< Database Members >>>
AXIS_Rrat = {}
AXIS_Rang = {}

Rise = ['AC-3.0_____', 'AC-3.2_____', 'AC-3.4_____', 'AC-3.6_____',
        'AC-3.8_____', 'AC-4.0_____', 'AC-4.4_____', 'AC-5.0_____']
dist = [3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.4, 5.0]
AXIS_Rang.update(dict(zip(Rise, dist)))
AXIS_Rrat.update(dict(zip(Rise, [d / 3.4 for d in dist])))

Shift = ['AC-3.4_n2.0____', 'AC-3.4_n1.6____', 'AC-3.4_n1.2____', 'AC-3.4_n0.8____',
         'AC-3.4_n0.4____', 'AC-3.4_____', 'AC-3.4_0.4____', 'AC-3.4_0.8____',
         'AC-3.4_1.2____', 'AC-3.4_1.6____', 'AC-3.4_2.0____']
dist = [-2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0]
AXIS_Rang.update(dict(zip(Shift, dist)))
AXIS_Rrat.update(dict(zip(Shift, [d / -0.4 for d in dist])))

Slide = ['AC-3.4__n2.0___', 'AC-3.4__n1.6___', 'AC-3.4__n1.2___', 'AC-3.4__n0.8___',
         'AC-3.4__n0.4___', 'AC-3.4_____', 'AC-3.4__0.4___', 'AC-3.4__0.8___',
         'AC-3.4__1.2___', 'AC-3.4__1.6___', 'AC-3.4__2.0___']
dist = [-2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0]
AXIS_Rang.update(dict(zip(Slide, dist)))
AXIS_Rrat.update(dict(zip(Slide, [d / 1.2 for d in dist])))

Twist = ['AC-3.4_____', 'AC-3.4___30__', 'AC-3.4___60__', 'AC-3.4___90__',
         'AC-3.4___120__', 'AC-3.4___150__', 'AC-3.4___180__']
dist = [0., 30., 60., 90., 120., 150., 180.]
AXIS_Rang.update(dict(zip(Twist, dist)))
AXIS_Rrat.update(dict(zip(Twist, [d / 150. for d in dist])))

Roll = ['AC-3.6____n20_', 'AC-3.6____n16_', 'AC-3.6____n12_',
        'AC-3.6____n8_', 'AC-3.6____n4_', 'AC-3.6_____', 'AC-3.6____4_',
        'AC-3.6____8_', 'AC-3.6____12_', 'AC-3.6____16_', 'AC-3.6____20_']
dist = [-20., -16., -12., -8., -4., 0., 4., 8., 12., 16., 20.]
AXIS_Rang.update(dict(zip(Roll, dist)))
AXIS_Rrat.update(dict(zip(Roll, [d / 4. for d in dist])))

Tilt = ['AC-3.6_____n20', 'AC-3.6_____n16', 'AC-3.6_____n12', 'AC-3.6_____n8',
        'AC-3.6_____n4', 'AC-3.6_____', 'AC-3.6_____4', 'AC-3.6_____8',
        'AC-3.6_____12', 'AC-3.6_____16', 'AC-3.6_____20']
dist = [-20., -16., -12., -8., -4., 0, 4., 8., 12., 16., 20.]
AXIS_Rang.update(dict(zip(Tilt, dist)))
AXIS_Rrat.update(dict(zip(Tilt, [d / -4. for d in dist])))

HRXN = []
temp = sum([Rise, Twist, Shift, Slide, Roll, Tilt], [])
[HRXN.append(i) for i in temp if not HRXN.count(i)]
#HRXN_SM = []
#HRXN_LG = []
HRXN_EQ = ['AC-3.4_____', 'AC-3.4_1.2____', 'AC-3.4__n0.4___',
           'AC-3.4___150__','AC-3.6____4_', 'AC-3.6_____n4']

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV_CP = {}  # order of active reagents per counterpoise-corrected reaction
ACTV_SA = {}  # order of active reagents for non-supermolecular calculations
for rxn in HRXN:

    RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                      '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                      '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                      '%s-%s-monoA-unCP' % (dbse, rxn) : -1,
                                      '%s-%s-monoB-unCP' % (dbse, rxn) : -1 }

    ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

    ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                      '%s-%s-monoA-CP'   % (dbse, rxn),
                                      '%s-%s-monoB-CP'   % (dbse, rxn) ]

    ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                      '%s-%s-monoA-unCP' % (dbse, rxn),
                                      '%s-%s-monoB-unCP' % (dbse, rxn) ]

# <<< Reference Values [kcal/mol] >>>
BIND_ACHC0 = {}
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.0_____'           )] = -3.141
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.2_____'           )] = -5.741
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4_____'           )] = -6.153  # rise eq
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6_____'           )] = -5.620
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.8_____'           )] = -4.768
BIND_ACHC0['%s-%s'            % (dbse, 'AC-4.0_____'           )] = -3.888
BIND_ACHC0['%s-%s'            % (dbse, 'AC-4.4_____'           )] = -2.445
BIND_ACHC0['%s-%s'            % (dbse, 'AC-5.0_____'           )] = -1.170

BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4_____'           )] = -6.153
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4___30__'         )] = -5.609
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4___60__'         )] = -4.834
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4___90__'         )] = -6.696
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4___120__'        )] = -6.967
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4___150__'        )] = -8.009  # twist eq
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4___180__'        )] = -6.533

BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4_n2.0____'       )] = -6.179
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4_n1.6____'       )] = -6.090
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4_n1.2____'       )] = -6.069
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4_n0.8____'       )] = -6.245
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4_n0.4____'       )] = -6.396  # shift eq
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4_____'           )] = -6.153
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4_0.4____'        )] = -5.470
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4_0.8____'        )] = -4.733
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4_1.2____'        )] = -4.331
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4_1.6____'        )] = -4.250
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4_2.0____'        )] = -4.186

BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4__n2.0___'       )] = -3.645
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4__n1.6___'       )] = -4.417
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4__n1.2___'       )] = -5.170
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4__n0.8___'       )] = -5.673
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4__n0.4___'       )] = -5.915
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4_____'           )] = -6.153
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4__0.4___'        )] = -6.608
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4__0.8___'        )] = -7.136
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4__1.2___'        )] = -7.345  # slide eq
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4__1.6___'        )] = -7.041
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.4__2.0___'        )] = -6.427

BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6____n20_'        )] = -4.360
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6____n16_'        )] = -5.075
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6____n12_'        )] = -5.420
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6____n8_'         )] = -5.564
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6____n4_'         )] = -5.611
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6_____'           )] = -5.620
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6____4_'          )] = -5.621  # roll eq
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6____8_'          )] = -5.617
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6____12_'         )] = -5.585
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6____16_'         )] = -5.464
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6____20_'         )] = -5.137

BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6_____n20'        )] = -4.373
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6_____n16'        )] = -5.029
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6_____n12'        )] = -5.382
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6_____n8'         )] = -5.555
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6_____n4'         )] = -5.620  # tilt eq
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6_____'           )] = -5.620
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6_____4'          )] = -5.572
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6_____8'          )] = -5.478
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6_____12'         )] = -5.330
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6_____16'         )] = -5.107
BIND_ACHC0['%s-%s'            % (dbse, 'AC-3.6_____20'         )] = -4.776
# Set default
BIND = BIND_ACHC0
# Reference information
BINDINFO_ACHC0 = {}
for rxn in HRXN:
    BINDINFO_ACHC0['%s-%s' % (dbse, rxn)] = {'citation': 'achc', 'method': 'DWCCSDTF12', 'mode': 'CP', 'basis': 'adz'}

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, 'AC-3.0_____'           )] = """a-c nucleobase complex at 3.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.0_____'           )] = """Dimer from a-c nucleobase complex at 3.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.0_____'           )] = """Monomer A from a-c nucleobase complex at 3.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.0_____'           )] = """Monomer B from a-c nucleobase complex at 3.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.0_____'           )] = """Monomer A from a-c nucleobase complex at 3.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.0_____'           )] = """Monomer B from a-c nucleobase complex at 3.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.2_____'           )] = """a-c nucleobase complex at 3.2 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.2_____'           )] = """Dimer from a-c nucleobase complex at 3.2 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.2_____'           )] = """Monomer A from a-c nucleobase complex at 3.2 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.2_____'           )] = """Monomer B from a-c nucleobase complex at 3.2 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.2_____'           )] = """Monomer A from a-c nucleobase complex at 3.2 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.2_____'           )] = """Monomer B from a-c nucleobase complex at 3.2 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.4_1.2____'        )] = """a-c nucleobase complex at 1.2 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4_1.2____'        )] = """Dimer from a-c nucleobase complex at 1.2 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4_1.2____'        )] = """Monomer A from a-c nucleobase complex at 1.2 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4_1.2____'        )] = """Monomer B from a-c nucleobase complex at 1.2 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4_1.2____'        )] = """Monomer A from a-c nucleobase complex at 1.2 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4_1.2____'        )] = """Monomer B from a-c nucleobase complex at 1.2 slide (A) and 0.0 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4_1.6____'        )] = """a-c nucleobase complex at 1.6 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4_1.6____'        )] = """Dimer from a-c nucleobase complex at 1.6 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4_1.6____'        )] = """Monomer A from a-c nucleobase complex at 1.6 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4_1.6____'        )] = """Monomer B from a-c nucleobase complex at 1.6 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4_1.6____'        )] = """Monomer A from a-c nucleobase complex at 1.6 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4_1.6____'        )] = """Monomer B from a-c nucleobase complex at 1.6 slide (A) and 0.0 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4_2.0____'        )] = """a-c nucleobase complex at 2.0 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4_2.0____'        )] = """Dimer from a-c nucleobase complex at 2.0 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4_2.0____'        )] = """Monomer A from a-c nucleobase complex at 2.0 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4_2.0____'        )] = """Monomer B from a-c nucleobase complex at 2.0 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4_2.0____'        )] = """Monomer A from a-c nucleobase complex at 2.0 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4_2.0____'        )] = """Monomer B from a-c nucleobase complex at 2.0 slide (A) and 0.0 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4_0.4____'        )] = """a-c nucleobase complex at 0.4 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4_0.4____'        )] = """Dimer from a-c nucleobase complex at 0.4 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4_0.4____'        )] = """Monomer A from a-c nucleobase complex at 0.4 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4_0.4____'        )] = """Monomer B from a-c nucleobase complex at 0.4 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4_0.4____'        )] = """Monomer A from a-c nucleobase complex at 0.4 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4_0.4____'        )] = """Monomer B from a-c nucleobase complex at 0.4 slide (A) and 0.0 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4_0.8____'        )] = """a-c nucleobase complex at 0.8 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4_0.8____'        )] = """Dimer from a-c nucleobase complex at 0.8 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4_0.8____'        )] = """Monomer A from a-c nucleobase complex at 0.8 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4_0.8____'        )] = """Monomer B from a-c nucleobase complex at 0.8 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4_0.8____'        )] = """Monomer A from a-c nucleobase complex at 0.8 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4_0.8____'        )] = """Monomer B from a-c nucleobase complex at 0.8 slide (A) and 0.0 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4__1.2___'        )] = """a-c nucleobase complex at 0.0 slide (A) and 1.2 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4__1.2___'        )] = """Dimer from a-c nucleobase complex at 0.0 slide (A) and 1.2 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4__1.2___'        )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and 1.2 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4__1.2___'        )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and 1.2 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4__1.2___'        )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and 1.2 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4__1.2___'        )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and 1.2 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4__1.6___'        )] = """a-c nucleobase complex at 0.0 slide (A) and 1.6 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4__1.6___'        )] = """Dimer from a-c nucleobase complex at 0.0 slide (A) and 1.6 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4__1.6___'        )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and 1.6 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4__1.6___'        )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and 1.6 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4__1.6___'        )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and 1.6 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4__1.6___'        )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and 1.6 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4__2.0___'        )] = """a-c nucleobase complex at 0.0 slide (A) and 2.0 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4__2.0___'        )] = """Dimer from a-c nucleobase complex at 0.0 slide (A) and 2.0 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4__2.0___'        )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and 2.0 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4__2.0___'        )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and 2.0 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4__2.0___'        )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and 2.0 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4__2.0___'        )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and 2.0 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4__0.4___'        )] = """a-c nucleobase complex at 0.0 slide (A) and 0.4 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4__0.4___'        )] = """Dimer from a-c nucleobase complex at 0.0 slide (A) and 0.4 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4__0.4___'        )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and 0.4 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4__0.4___'        )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and 0.4 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4__0.4___'        )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and 0.4 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4__0.4___'        )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and 0.4 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4__0.8___'        )] = """a-c nucleobase complex at 0.0 slide (A) and 0.8 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4__0.8___'        )] = """Dimer from a-c nucleobase complex at 0.0 slide (A) and 0.8 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4__0.8___'        )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and 0.8 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4__0.8___'        )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and 0.8 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4__0.8___'        )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and 0.8 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4__0.8___'        )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and 0.8 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4___120__'        )] = """a-c nucleobase complex at 3.4 rise (A) and 120.0 twist (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4___120__'        )] = """Dimer from a-c nucleobase complex at 3.4 rise (A) and 120.0 twist (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4___120__'        )] = """Monomer A from a-c nucleobase complex at 3.4 rise (A) and 120.0 twist (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4___120__'        )] = """Monomer B from a-c nucleobase complex at 3.4 rise (A) and 120.0 twist (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4___120__'        )] = """Monomer A from a-c nucleobase complex at 3.4 rise (A) and 120.0 twist (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4___120__'        )] = """Monomer B from a-c nucleobase complex at 3.4 rise (A) and 120.0 twist (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.4___150__'        )] = """a-c nucleobase complex at 3.4 rise (A) and 150.0 twist (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4___150__'        )] = """Dimer from a-c nucleobase complex at 3.4 rise (A) and 150.0 twist (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4___150__'        )] = """Monomer A from a-c nucleobase complex at 3.4 rise (A) and 150.0 twist (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4___150__'        )] = """Monomer B from a-c nucleobase complex at 3.4 rise (A) and 150.0 twist (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4___150__'        )] = """Monomer A from a-c nucleobase complex at 3.4 rise (A) and 150.0 twist (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4___150__'        )] = """Monomer B from a-c nucleobase complex at 3.4 rise (A) and 150.0 twist (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.4___180__'        )] = """a-c nucleobase complex at 3.4 rise (A) and 180.0 twist (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4___180__'        )] = """Dimer from a-c nucleobase complex at 3.4 rise (A) and 180.0 twist (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4___180__'        )] = """Monomer A from a-c nucleobase complex at 3.4 rise (A) and 180.0 twist (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4___180__'        )] = """Monomer B from a-c nucleobase complex at 3.4 rise (A) and 180.0 twist (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4___180__'        )] = """Monomer A from a-c nucleobase complex at 3.4 rise (A) and 180.0 twist (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4___180__'        )] = """Monomer B from a-c nucleobase complex at 3.4 rise (A) and 180.0 twist (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.4___30__'         )] = """a-c nucleobase complex at 3.4 rise (A) and 30.0 twist (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4___30__'         )] = """Dimer from a-c nucleobase complex at 3.4 rise (A) and 30.0 twist (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4___30__'         )] = """Monomer A from a-c nucleobase complex at 3.4 rise (A) and 30.0 twist (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4___30__'         )] = """Monomer B from a-c nucleobase complex at 3.4 rise (A) and 30.0 twist (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4___30__'         )] = """Monomer A from a-c nucleobase complex at 3.4 rise (A) and 30.0 twist (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4___30__'         )] = """Monomer B from a-c nucleobase complex at 3.4 rise (A) and 30.0 twist (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.4___60__'         )] = """a-c nucleobase complex at 3.4 rise (A) and 60.0 twist (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4___60__'         )] = """Dimer from a-c nucleobase complex at 3.4 rise (A) and 60.0 twist (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4___60__'         )] = """Monomer A from a-c nucleobase complex at 3.4 rise (A) and 60.0 twist (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4___60__'         )] = """Monomer B from a-c nucleobase complex at 3.4 rise (A) and 60.0 twist (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4___60__'         )] = """Monomer A from a-c nucleobase complex at 3.4 rise (A) and 60.0 twist (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4___60__'         )] = """Monomer B from a-c nucleobase complex at 3.4 rise (A) and 60.0 twist (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.4___90__'         )] = """a-c nucleobase complex at 3.4 rise (A) and 90.0 twist (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4___90__'         )] = """Dimer from a-c nucleobase complex at 3.4 rise (A) and 90.0 twist (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4___90__'         )] = """Monomer A from a-c nucleobase complex at 3.4 rise (A) and 90.0 twist (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4___90__'         )] = """Monomer B from a-c nucleobase complex at 3.4 rise (A) and 90.0 twist (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4___90__'         )] = """Monomer A from a-c nucleobase complex at 3.4 rise (A) and 90.0 twist (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4___90__'         )] = """Monomer B from a-c nucleobase complex at 3.4 rise (A) and 90.0 twist (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.4_____'           )] = """a-c nucleobase complex at 3.4 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4_____'           )] = """Dimer from a-c nucleobase complex at 3.4 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4_____'           )] = """Monomer A from a-c nucleobase complex at 3.4 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4_____'           )] = """Monomer B from a-c nucleobase complex at 3.4 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4_____'           )] = """Monomer A from a-c nucleobase complex at 3.4 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4_____'           )] = """Monomer B from a-c nucleobase complex at 3.4 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.4__n1.2___'       )] = """a-c nucleobase complex at 0.0 slide (A) and -1.2 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4__n1.2___'       )] = """Dimer from a-c nucleobase complex at 0.0 slide (A) and -1.2 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4__n1.2___'       )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and -1.2 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4__n1.2___'       )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and -1.2 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4__n1.2___'       )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and -1.2 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4__n1.2___'       )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and -1.2 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4__n1.6___'       )] = """a-c nucleobase complex at 0.0 slide (A) and -1.6 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4__n1.6___'       )] = """Dimer from a-c nucleobase complex at 0.0 slide (A) and -1.6 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4__n1.6___'       )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and -1.6 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4__n1.6___'       )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and -1.6 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4__n1.6___'       )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and -1.6 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4__n1.6___'       )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and -1.6 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4__n2.0___'       )] = """a-c nucleobase complex at 0.0 slide (A) and -2.0 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4__n2.0___'       )] = """Dimer from a-c nucleobase complex at 0.0 slide (A) and -2.0 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4__n2.0___'       )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and -2.0 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4__n2.0___'       )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and -2.0 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4__n2.0___'       )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and -2.0 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4__n2.0___'       )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and -2.0 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4__n0.4___'       )] = """a-c nucleobase complex at 0.0 slide (A) and -0.4 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4__n0.4___'       )] = """Dimer from a-c nucleobase complex at 0.0 slide (A) and -0.4 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4__n0.4___'       )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and -0.4 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4__n0.4___'       )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and -0.4 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4__n0.4___'       )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and -0.4 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4__n0.4___'       )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and -0.4 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4__n0.8___'       )] = """a-c nucleobase complex at 0.0 slide (A) and -0.8 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4__n0.8___'       )] = """Dimer from a-c nucleobase complex at 0.0 slide (A) and -0.8 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4__n0.8___'       )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and -0.8 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4__n0.8___'       )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and -0.8 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4__n0.8___'       )] = """Monomer A from a-c nucleobase complex at 0.0 slide (A) and -0.8 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4__n0.8___'       )] = """Monomer B from a-c nucleobase complex at 0.0 slide (A) and -0.8 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4_n1.2____'       )] = """a-c nucleobase complex at -1.2 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4_n1.2____'       )] = """Dimer from a-c nucleobase complex at -1.2 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4_n1.2____'       )] = """Monomer A from a-c nucleobase complex at -1.2 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4_n1.2____'       )] = """Monomer B from a-c nucleobase complex at -1.2 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4_n1.2____'       )] = """Monomer A from a-c nucleobase complex at -1.2 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4_n1.2____'       )] = """Monomer B from a-c nucleobase complex at -1.2 slide (A) and 0.0 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4_n1.6____'       )] = """a-c nucleobase complex at -1.6 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4_n1.6____'       )] = """Dimer from a-c nucleobase complex at -1.6 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4_n1.6____'       )] = """Monomer A from a-c nucleobase complex at -1.6 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4_n1.6____'       )] = """Monomer B from a-c nucleobase complex at -1.6 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4_n1.6____'       )] = """Monomer A from a-c nucleobase complex at -1.6 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4_n1.6____'       )] = """Monomer B from a-c nucleobase complex at -1.6 slide (A) and 0.0 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4_n2.0____'       )] = """a-c nucleobase complex at -2.0 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4_n2.0____'       )] = """Dimer from a-c nucleobase complex at -2.0 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4_n2.0____'       )] = """Monomer A from a-c nucleobase complex at -2.0 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4_n2.0____'       )] = """Monomer B from a-c nucleobase complex at -2.0 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4_n2.0____'       )] = """Monomer A from a-c nucleobase complex at -2.0 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4_n2.0____'       )] = """Monomer B from a-c nucleobase complex at -2.0 slide (A) and 0.0 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4_n0.4____'       )] = """a-c nucleobase complex at -0.4 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4_n0.4____'       )] = """Dimer from a-c nucleobase complex at -0.4 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4_n0.4____'       )] = """Monomer A from a-c nucleobase complex at -0.4 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4_n0.4____'       )] = """Monomer B from a-c nucleobase complex at -0.4 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4_n0.4____'       )] = """Monomer A from a-c nucleobase complex at -0.4 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4_n0.4____'       )] = """Monomer B from a-c nucleobase complex at -0.4 slide (A) and 0.0 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.4_n0.8____'       )] = """a-c nucleobase complex at -0.8 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.4_n0.8____'       )] = """Dimer from a-c nucleobase complex at -0.8 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.4_n0.8____'       )] = """Monomer A from a-c nucleobase complex at -0.8 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.4_n0.8____'       )] = """Monomer B from a-c nucleobase complex at -0.8 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.4_n0.8____'       )] = """Monomer A from a-c nucleobase complex at -0.8 slide (A) and 0.0 shift (A) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.4_n0.8____'       )] = """Monomer B from a-c nucleobase complex at -0.8 slide (A) and 0.0 shift (A) """
TAGL['%s-%s'            % (dbse, 'AC-3.6____12_'        )] = """a-c nucleobase complex at 0.0 tilt (deg) and 12.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6____12_'        )] = """Dimer from a-c nucleobase complex at 0.0 tilt (deg) and 12.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6____12_'        )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and 12.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6____12_'        )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and 12.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6____12_'        )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and 12.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6____12_'        )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and 12.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6____16_'        )] = """a-c nucleobase complex at 0.0 tilt (deg) and 16.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6____16_'        )] = """Dimer from a-c nucleobase complex at 0.0 tilt (deg) and 16.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6____16_'        )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and 16.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6____16_'        )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and 16.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6____16_'        )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and 16.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6____16_'        )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and 16.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6____20_'        )] = """a-c nucleobase complex at 0.0 tilt (deg) and 20.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6____20_'        )] = """Dimer from a-c nucleobase complex at 0.0 tilt (deg) and 20.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6____20_'        )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and 20.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6____20_'        )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and 20.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6____20_'        )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and 20.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6____20_'        )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and 20.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6____4_'         )] = """a-c nucleobase complex at 0.0 tilt (deg) and 4.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6____4_'         )] = """Dimer from a-c nucleobase complex at 0.0 tilt (deg) and 4.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6____4_'         )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and 4.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6____4_'         )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and 4.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6____4_'         )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and 4.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6____4_'         )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and 4.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6____8_'         )] = """a-c nucleobase complex at 0.0 tilt (deg) and 8.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6____8_'         )] = """Dimer from a-c nucleobase complex at 0.0 tilt (deg) and 8.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6____8_'         )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and 8.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6____8_'         )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and 8.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6____8_'         )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and 8.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6____8_'         )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and 8.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6_____'           )] = """a-c nucleobase complex at 3.6 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6_____'           )] = """Dimer from a-c nucleobase complex at 3.6 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6_____'           )] = """Monomer A from a-c nucleobase complex at 3.6 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6_____'           )] = """Monomer B from a-c nucleobase complex at 3.6 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6_____'           )] = """Monomer A from a-c nucleobase complex at 3.6 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6_____'           )] = """Monomer B from a-c nucleobase complex at 3.6 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6_____12'        )] = """a-c nucleobase complex at 12.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6_____12'        )] = """Dimer from a-c nucleobase complex at 12.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6_____12'        )] = """Monomer A from a-c nucleobase complex at 12.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6_____12'        )] = """Monomer B from a-c nucleobase complex at 12.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6_____12'        )] = """Monomer A from a-c nucleobase complex at 12.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6_____12'        )] = """Monomer B from a-c nucleobase complex at 12.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6_____16'        )] = """a-c nucleobase complex at 16.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6_____16'        )] = """Dimer from a-c nucleobase complex at 16.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6_____16'        )] = """Monomer A from a-c nucleobase complex at 16.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6_____16'        )] = """Monomer B from a-c nucleobase complex at 16.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6_____16'        )] = """Monomer A from a-c nucleobase complex at 16.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6_____16'        )] = """Monomer B from a-c nucleobase complex at 16.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6_____20'        )] = """a-c nucleobase complex at 20.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6_____20'        )] = """Dimer from a-c nucleobase complex at 20.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6_____20'        )] = """Monomer A from a-c nucleobase complex at 20.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6_____20'        )] = """Monomer B from a-c nucleobase complex at 20.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6_____20'        )] = """Monomer A from a-c nucleobase complex at 20.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6_____20'        )] = """Monomer B from a-c nucleobase complex at 20.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6_____4'         )] = """a-c nucleobase complex at 4.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6_____4'         )] = """Dimer from a-c nucleobase complex at 4.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6_____4'         )] = """Monomer A from a-c nucleobase complex at 4.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6_____4'         )] = """Monomer B from a-c nucleobase complex at 4.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6_____4'         )] = """Monomer A from a-c nucleobase complex at 4.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6_____4'         )] = """Monomer B from a-c nucleobase complex at 4.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6_____8'         )] = """a-c nucleobase complex at 8.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6_____8'         )] = """Dimer from a-c nucleobase complex at 8.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6_____8'         )] = """Monomer A from a-c nucleobase complex at 8.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6_____8'         )] = """Monomer B from a-c nucleobase complex at 8.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6_____8'         )] = """Monomer A from a-c nucleobase complex at 8.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6_____8'         )] = """Monomer B from a-c nucleobase complex at 8.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6_____n12'       )] = """a-c nucleobase complex at -12.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6_____n12'       )] = """Dimer from a-c nucleobase complex at -12.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6_____n12'       )] = """Monomer A from a-c nucleobase complex at -12.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6_____n12'       )] = """Monomer B from a-c nucleobase complex at -12.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6_____n12'       )] = """Monomer A from a-c nucleobase complex at -12.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6_____n12'       )] = """Monomer B from a-c nucleobase complex at -12.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6_____n16'       )] = """a-c nucleobase complex at -16.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6_____n16'       )] = """Dimer from a-c nucleobase complex at -16.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6_____n16'       )] = """Monomer A from a-c nucleobase complex at -16.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6_____n16'       )] = """Monomer B from a-c nucleobase complex at -16.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6_____n16'       )] = """Monomer A from a-c nucleobase complex at -16.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6_____n16'       )] = """Monomer B from a-c nucleobase complex at -16.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6_____n20'       )] = """a-c nucleobase complex at -20.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6_____n20'       )] = """Dimer from a-c nucleobase complex at -20.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6_____n20'       )] = """Monomer A from a-c nucleobase complex at -20.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6_____n20'       )] = """Monomer B from a-c nucleobase complex at -20.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6_____n20'       )] = """Monomer A from a-c nucleobase complex at -20.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6_____n20'       )] = """Monomer B from a-c nucleobase complex at -20.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6_____n4'        )] = """a-c nucleobase complex at -4.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6_____n4'        )] = """Dimer from a-c nucleobase complex at -4.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6_____n4'        )] = """Monomer A from a-c nucleobase complex at -4.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6_____n4'        )] = """Monomer B from a-c nucleobase complex at -4.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6_____n4'        )] = """Monomer A from a-c nucleobase complex at -4.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6_____n4'        )] = """Monomer B from a-c nucleobase complex at -4.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6_____n8'        )] = """a-c nucleobase complex at -8.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6_____n8'        )] = """Dimer from a-c nucleobase complex at -8.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6_____n8'        )] = """Monomer A from a-c nucleobase complex at -8.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6_____n8'        )] = """Monomer B from a-c nucleobase complex at -8.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6_____n8'        )] = """Monomer A from a-c nucleobase complex at -8.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6_____n8'        )] = """Monomer B from a-c nucleobase complex at -8.0 tilt (deg) and 0.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6____n12_'       )] = """a-c nucleobase complex at 0.0 tilt (deg) and -12.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6____n12_'       )] = """Dimer from a-c nucleobase complex at 0.0 tilt (deg) and -12.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6____n12_'       )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and -12.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6____n12_'       )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and -12.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6____n12_'       )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and -12.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6____n12_'       )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and -12.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6____n16_'       )] = """a-c nucleobase complex at 0.0 tilt (deg) and -16.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6____n16_'       )] = """Dimer from a-c nucleobase complex at 0.0 tilt (deg) and -16.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6____n16_'       )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and -16.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6____n16_'       )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and -16.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6____n16_'       )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and -16.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6____n16_'       )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and -16.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6____n20_'       )] = """a-c nucleobase complex at 0.0 tilt (deg) and -20.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6____n20_'       )] = """Dimer from a-c nucleobase complex at 0.0 tilt (deg) and -20.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6____n20_'       )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and -20.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6____n20_'       )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and -20.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6____n20_'       )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and -20.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6____n20_'       )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and -20.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6____n4_'        )] = """a-c nucleobase complex at 0.0 tilt (deg) and -4.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6____n4_'        )] = """Dimer from a-c nucleobase complex at 0.0 tilt (deg) and -4.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6____n4_'        )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and -4.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6____n4_'        )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and -4.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6____n4_'        )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and -4.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6____n4_'        )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and -4.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.6____n8_'        )] = """a-c nucleobase complex at 0.0 tilt (deg) and -8.0 roll (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.6____n8_'        )] = """Dimer from a-c nucleobase complex at 0.0 tilt (deg) and -8.0 roll (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.6____n8_'        )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and -8.0 roll (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.6____n8_'        )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and -8.0 roll (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.6____n8_'        )] = """Monomer A from a-c nucleobase complex at 0.0 tilt (deg) and -8.0 roll (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.6____n8_'        )] = """Monomer B from a-c nucleobase complex at 0.0 tilt (deg) and -8.0 roll (deg) """
TAGL['%s-%s'            % (dbse, 'AC-3.8_____'           )] = """a-c nucleobase complex at 3.8 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-3.8_____'           )] = """Dimer from a-c nucleobase complex at 3.8 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-3.8_____'           )] = """Monomer A from a-c nucleobase complex at 3.8 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-3.8_____'           )] = """Monomer B from a-c nucleobase complex at 3.8 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-3.8_____'           )] = """Monomer A from a-c nucleobase complex at 3.8 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-3.8_____'           )] = """Monomer B from a-c nucleobase complex at 3.8 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s'            % (dbse, 'AC-4.0_____'           )] = """a-c nucleobase complex at 4.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-4.0_____'           )] = """Dimer from a-c nucleobase complex at 4.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-4.0_____'           )] = """Monomer A from a-c nucleobase complex at 4.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-4.0_____'           )] = """Monomer B from a-c nucleobase complex at 4.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-4.0_____'           )] = """Monomer A from a-c nucleobase complex at 4.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-4.0_____'           )] = """Monomer B from a-c nucleobase complex at 4.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s'            % (dbse, 'AC-4.4_____'           )] = """a-c nucleobase complex at 4.4 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-4.4_____'           )] = """Dimer from a-c nucleobase complex at 4.4 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-4.4_____'           )] = """Monomer A from a-c nucleobase complex at 4.4 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-4.4_____'           )] = """Monomer B from a-c nucleobase complex at 4.4 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-4.4_____'           )] = """Monomer A from a-c nucleobase complex at 4.4 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-4.4_____'           )] = """Monomer B from a-c nucleobase complex at 4.4 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s'            % (dbse, 'AC-5.0_____'           )] = """a-c nucleobase complex at 5.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-dimer'      % (dbse, 'AC-5.0_____'           )] = """Dimer from a-c nucleobase complex at 5.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'AC-5.0_____'           )] = """Monomer A from a-c nucleobase complex at 5.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'AC-5.0_____'           )] = """Monomer B from a-c nucleobase complex at 5.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'AC-5.0_____'           )] = """Monomer A from a-c nucleobase complex at 5.0 rise (A) and 0.0 twist (deg) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'AC-5.0_____'           )] = """Monomer B from a-c nucleobase complex at 5.0 rise (A) and 0.0 twist (deg) """

TAGL['dbse'] = 'interaction energy curves for adenine-cytosine stacked nucleobases through six translations and rotations'
TAGL['rise'] = 'translation curve in +Z, up the helix'
TAGL['twist'] = 'rotation curve ccw in +Z, up the helix'
TAGL['shift'] = 'translation curve in +X, toward the major groove'
TAGL['slide'] = "translation curve in +Y, toward the 5'-3' strand"
TAGL['roll'] = "rotation curve ccw in +Y, toward the 5'-3' strand"
TAGL['tilt'] = 'rotation curve ccw in +X, toward the major groove'
TAGL['equilibrium'] = 'minimum-energy systems on translation and rotation curves'
TAGL['default'] = 'entire database'

# <<< Geometry Specification Strings >>>
GEOS = {}

GEOS['%s-%s-%s' % (dbse, 'AC-3.0_____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.197258000000     3.000058000000
    C                1.217664000000    -0.193037000000     3.000004000000
    N                0.132391000000    -1.008972000000     2.999965000000
    C               -1.110962000000    -0.503844000000     2.999971000000
    C               -1.352449000000     0.921627000000     3.000023000000
    C               -0.263556000000     1.733781000000     3.000063000000
    H                1.815200000000     1.780954000000     3.000092000000
    O                2.388816000000    -0.594965000000     2.999991000000
    N               -2.130425000000    -1.368706000000     2.999929000000
    H               -1.956075000000    -2.396574000000     2.999871000000
    H               -3.078136000000    -1.025139000000     2.999928000000
    H               -2.359817000000     1.329225000000     3.000032000000
    H               -0.332189000000     2.820853000000     3.000106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.2_____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.197258000000     3.200058000000
    C                1.217664000000    -0.193037000000     3.200004000000
    N                0.132391000000    -1.008972000000     3.199965000000
    C               -1.110962000000    -0.503844000000     3.199971000000
    C               -1.352449000000     0.921627000000     3.200023000000
    C               -0.263556000000     1.733781000000     3.200063000000
    H                1.815200000000     1.780954000000     3.200092000000
    O                2.388816000000    -0.594965000000     3.199991000000
    N               -2.130425000000    -1.368706000000     3.199929000000
    H               -1.956075000000    -2.396574000000     3.199871000000
    H               -3.078136000000    -1.025139000000     3.199928000000
    H               -2.359817000000     1.329225000000     3.200032000000
    H               -0.332189000000     2.820853000000     3.200106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4_1.2____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     2.397258000000     3.400058000000
    C                1.217664000000     1.006963000000     3.400004000000
    N                0.132391000000     0.191028000000     3.399965000000
    C               -1.110962000000     0.696156000000     3.399971000000
    C               -1.352449000000     2.121627000000     3.400023000000
    C               -0.263556000000     2.933781000000     3.400063000000
    H                1.815200000000     2.980954000000     3.400092000000
    O                2.388816000000     0.605035000000     3.399991000000
    N               -2.130425000000    -0.168706000000     3.399929000000
    H               -1.956075000000    -1.196574000000     3.399871000000
    H               -3.078136000000     0.174861000000     3.399928000000
    H               -2.359817000000     2.529225000000     3.400032000000
    H               -0.332189000000     4.020853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4_1.6____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     2.797258000000     3.400058000000
    C                1.217664000000     1.406963000000     3.400004000000
    N                0.132391000000     0.591028000000     3.399965000000
    C               -1.110962000000     1.096156000000     3.399971000000
    C               -1.352449000000     2.521627000000     3.400023000000
    C               -0.263556000000     3.333781000000     3.400063000000
    H                1.815200000000     3.380954000000     3.400092000000
    O                2.388816000000     1.005035000000     3.399991000000
    N               -2.130425000000     0.231294000000     3.399929000000
    H               -1.956075000000    -0.796574000000     3.399871000000
    H               -3.078136000000     0.574861000000     3.399928000000
    H               -2.359817000000     2.929225000000     3.400032000000
    H               -0.332189000000     4.420853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4_2.0____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     3.197258000000     3.400058000000
    C                1.217664000000     1.806963000000     3.400004000000
    N                0.132391000000     0.991028000000     3.399965000000
    C               -1.110962000000     1.496156000000     3.399971000000
    C               -1.352449000000     2.921627000000     3.400023000000
    C               -0.263556000000     3.733781000000     3.400063000000
    H                1.815200000000     3.780954000000     3.400092000000
    O                2.388816000000     1.405035000000     3.399991000000
    N               -2.130425000000     0.631294000000     3.399929000000
    H               -1.956075000000    -0.396574000000     3.399871000000
    H               -3.078136000000     0.974861000000     3.399928000000
    H               -2.359817000000     3.329225000000     3.400032000000
    H               -0.332189000000     4.820853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4_0.4____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.597258000000     3.400058000000
    C                1.217664000000     0.206963000000     3.400004000000
    N                0.132391000000    -0.608972000000     3.399965000000
    C               -1.110962000000    -0.103844000000     3.399971000000
    C               -1.352449000000     1.321627000000     3.400023000000
    C               -0.263556000000     2.133781000000     3.400063000000
    H                1.815200000000     2.180954000000     3.400092000000
    O                2.388816000000    -0.194965000000     3.399991000000
    N               -2.130425000000    -0.968706000000     3.399929000000
    H               -1.956075000000    -1.996574000000     3.399871000000
    H               -3.078136000000    -0.625139000000     3.399928000000
    H               -2.359817000000     1.729225000000     3.400032000000
    H               -0.332189000000     3.220853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4_0.8____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.997258000000     3.400058000000
    C                1.217664000000     0.606963000000     3.400004000000
    N                0.132391000000    -0.208972000000     3.399965000000
    C               -1.110962000000     0.296156000000     3.399971000000
    C               -1.352449000000     1.721627000000     3.400023000000
    C               -0.263556000000     2.533781000000     3.400063000000
    H                1.815200000000     2.580954000000     3.400092000000
    O                2.388816000000     0.205035000000     3.399991000000
    N               -2.130425000000    -0.568706000000     3.399929000000
    H               -1.956075000000    -1.596574000000     3.399871000000
    H               -3.078136000000    -0.225139000000     3.399928000000
    H               -2.359817000000     2.129225000000     3.400032000000
    H               -0.332189000000     3.620853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4__1.2___', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                2.189061000000     1.197258000000     3.400058000000
    C                2.417664000000    -0.193037000000     3.400004000000
    N                1.332391000000    -1.008972000000     3.399965000000
    C                0.089038000000    -0.503844000000     3.399971000000
    C               -0.152449000000     0.921627000000     3.400023000000
    C                0.936444000000     1.733781000000     3.400063000000
    H                3.015200000000     1.780954000000     3.400092000000
    O                3.588816000000    -0.594965000000     3.399991000000
    N               -0.930425000000    -1.368706000000     3.399929000000
    H               -0.756075000000    -2.396574000000     3.399871000000
    H               -1.878136000000    -1.025139000000     3.399928000000
    H               -1.159817000000     1.329225000000     3.400032000000
    H                0.867811000000     2.820853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4__1.6___', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                2.589061000000     1.197258000000     3.400058000000
    C                2.817664000000    -0.193037000000     3.400004000000
    N                1.732391000000    -1.008972000000     3.399965000000
    C                0.489038000000    -0.503844000000     3.399971000000
    C                0.247551000000     0.921627000000     3.400023000000
    C                1.336444000000     1.733781000000     3.400063000000
    H                3.415200000000     1.780954000000     3.400092000000
    O                3.988816000000    -0.594965000000     3.399991000000
    N               -0.530425000000    -1.368706000000     3.399929000000
    H               -0.356075000000    -2.396574000000     3.399871000000
    H               -1.478136000000    -1.025139000000     3.399928000000
    H               -0.759817000000     1.329225000000     3.400032000000
    H                1.267811000000     2.820853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4__2.0___', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                2.989061000000     1.197258000000     3.400058000000
    C                3.217664000000    -0.193037000000     3.400004000000
    N                2.132391000000    -1.008972000000     3.399965000000
    C                0.889038000000    -0.503844000000     3.399971000000
    C                0.647551000000     0.921627000000     3.400023000000
    C                1.736444000000     1.733781000000     3.400063000000
    H                3.815200000000     1.780954000000     3.400092000000
    O                4.388816000000    -0.594965000000     3.399991000000
    N               -0.130425000000    -1.368706000000     3.399929000000
    H                0.043925000000    -2.396574000000     3.399871000000
    H               -1.078136000000    -1.025139000000     3.399928000000
    H               -0.359817000000     1.329225000000     3.400032000000
    H                1.667811000000     2.820853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4__0.4___', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                1.389061000000     1.197258000000     3.400058000000
    C                1.617664000000    -0.193037000000     3.400004000000
    N                0.532391000000    -1.008972000000     3.399965000000
    C               -0.710962000000    -0.503844000000     3.399971000000
    C               -0.952449000000     0.921627000000     3.400023000000
    C                0.136444000000     1.733781000000     3.400063000000
    H                2.215200000000     1.780954000000     3.400092000000
    O                2.788816000000    -0.594965000000     3.399991000000
    N               -1.730425000000    -1.368706000000     3.399929000000
    H               -1.556075000000    -2.396574000000     3.399871000000
    H               -2.678136000000    -1.025139000000     3.399928000000
    H               -1.959817000000     1.329225000000     3.400032000000
    H                0.067811000000     2.820853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4__0.8___', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                1.789061000000     1.197258000000     3.400058000000
    C                2.017664000000    -0.193037000000     3.400004000000
    N                0.932391000000    -1.008972000000     3.399965000000
    C               -0.310962000000    -0.503844000000     3.399971000000
    C               -0.552449000000     0.921627000000     3.400023000000
    C                0.536444000000     1.733781000000     3.400063000000
    H                2.615200000000     1.780954000000     3.400092000000
    O                3.188816000000    -0.594965000000     3.399991000000
    N               -1.330425000000    -1.368706000000     3.399929000000
    H               -1.156075000000    -2.396574000000     3.399871000000
    H               -2.278136000000    -1.025139000000     3.399928000000
    H               -1.559817000000     1.329225000000     3.400032000000
    H                0.467811000000     2.820853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4___120__', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N               -1.531386000000     0.257923000000     3.400058000000
    C               -0.441657000000     1.151046000000     3.400004000000
    N                0.807600000000     0.619140000000     3.399965000000
    C                0.991823000000    -0.710199000000     3.399971000000
    C               -0.121928000000    -1.632069000000     3.400023000000
    C               -1.369720000000    -1.095137000000     3.400063000000
    H               -2.449951000000     0.681532000000     3.400092000000
    O               -0.679153000000     2.366258000000     3.399991000000
    N                2.250547000000    -1.160649000000     3.399929000000
    H                3.053531000000    -0.495724000000     3.399871000000
    H                2.426864000000    -2.153174000000     3.399928000000
    H                0.028766000000    -2.708274000000     3.400032000000
    H               -2.276836000000    -1.698111000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4___150__', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N               -1.455181000000    -0.542325000000     3.400058000000
    C               -0.958009000000     0.776007000000     3.400004000000
    N                0.389832000000     0.939991000000     3.399965000000
    C                1.214043000000    -0.119139000000     3.399971000000
    C                0.710442000000    -1.474377000000     3.400023000000
    C               -0.638644000000    -1.633276000000     3.400063000000
    H               -2.462486000000    -0.634751000000     3.400092000000
    O               -1.771293000000     1.709663000000     3.399991000000
    N                2.529355000000     0.120122000000     3.399929000000
    H                2.892298000000     1.097456000000     3.399871000000
    H                3.178313000000    -0.651272000000     3.399928000000
    H                1.379049000000    -2.331051000000     3.400032000000
    H               -1.122742000000    -2.609025000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4___180__', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N               -0.989061000000    -1.197258000000     3.400058000000
    C               -1.217664000000     0.193037000000     3.400004000000
    N               -0.132391000000     1.008972000000     3.399965000000
    C                1.110962000000     0.503844000000     3.399971000000
    C                1.352449000000    -0.921627000000     3.400023000000
    C                0.263556000000    -1.733781000000     3.400063000000
    H               -1.815200000000    -1.780954000000     3.400092000000
    O               -2.388816000000     0.594965000000     3.399991000000
    N                2.130425000000     1.368706000000     3.399929000000
    H                1.956075000000     2.396574000000     3.399871000000
    H                3.078136000000     1.025139000000     3.399928000000
    H                2.359817000000    -1.329225000000     3.400032000000
    H                0.332189000000    -2.820853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4___30__', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.257923000000     1.531386000000     3.400058000000
    C                1.151046000000     0.441657000000     3.400004000000
    N                0.619140000000    -0.807600000000     3.399965000000
    C               -0.710199000000    -0.991823000000     3.399971000000
    C               -1.632069000000     0.121928000000     3.400023000000
    C               -1.095137000000     1.369720000000     3.400063000000
    H                0.681532000000     2.449951000000     3.400092000000
    O                2.366258000000     0.679153000000     3.399991000000
    N               -1.160649000000    -2.250547000000     3.399929000000
    H               -0.495724000000    -3.053531000000     3.399871000000
    H               -2.153174000000    -2.426864000000     3.399928000000
    H               -2.708274000000    -0.028766000000     3.400032000000
    H               -1.698111000000     2.276836000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4___60__', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N               -0.542325000000     1.455181000000     3.400058000000
    C                0.776007000000     0.958009000000     3.400004000000
    N                0.939991000000    -0.389832000000     3.399965000000
    C               -0.119139000000    -1.214043000000     3.399971000000
    C               -1.474377000000    -0.710442000000     3.400023000000
    C               -1.633276000000     0.638644000000     3.400063000000
    H               -0.634751000000     2.462486000000     3.400092000000
    O                1.709663000000     1.771293000000     3.399991000000
    N                0.120122000000    -2.529355000000     3.399929000000
    H                1.097456000000    -2.892298000000     3.399871000000
    H               -0.651272000000    -3.178313000000     3.399928000000
    H               -2.331051000000    -1.379049000000     3.400032000000
    H               -2.609025000000     1.122742000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4___90__', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N               -1.197258000000     0.989061000000     3.400058000000
    C                0.193037000000     1.217664000000     3.400004000000
    N                1.008972000000     0.132391000000     3.399965000000
    C                0.503844000000    -1.110962000000     3.399971000000
    C               -0.921627000000    -1.352449000000     3.400023000000
    C               -1.733781000000    -0.263556000000     3.400063000000
    H               -1.780954000000     1.815200000000     3.400092000000
    O                0.594965000000     2.388816000000     3.399991000000
    N                1.368706000000    -2.130425000000     3.399929000000
    H                2.396574000000    -1.956075000000     3.399871000000
    H                1.025139000000    -3.078136000000     3.399928000000
    H               -1.329225000000    -2.359817000000     3.400032000000
    H               -2.820853000000    -0.332189000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4_____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.197258000000     3.400058000000
    C                1.217664000000    -0.193037000000     3.400004000000
    N                0.132391000000    -1.008972000000     3.399965000000
    C               -1.110962000000    -0.503844000000     3.399971000000
    C               -1.352449000000     0.921627000000     3.400023000000
    C               -0.263556000000     1.733781000000     3.400063000000
    H                1.815200000000     1.780954000000     3.400092000000
    O                2.388816000000    -0.594965000000     3.399991000000
    N               -2.130425000000    -1.368706000000     3.399929000000
    H               -1.956075000000    -2.396574000000     3.399871000000
    H               -3.078136000000    -1.025139000000     3.399928000000
    H               -2.359817000000     1.329225000000     3.400032000000
    H               -0.332189000000     2.820853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4__n1.2___', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N               -0.210939000000     1.197258000000     3.400058000000
    C                0.017664000000    -0.193037000000     3.400004000000
    N               -1.067609000000    -1.008972000000     3.399965000000
    C               -2.310962000000    -0.503844000000     3.399971000000
    C               -2.552449000000     0.921627000000     3.400023000000
    C               -1.463556000000     1.733781000000     3.400063000000
    H                0.615200000000     1.780954000000     3.400092000000
    O                1.188816000000    -0.594965000000     3.399991000000
    N               -3.330425000000    -1.368706000000     3.399929000000
    H               -3.156075000000    -2.396574000000     3.399871000000
    H               -4.278136000000    -1.025139000000     3.399928000000
    H               -3.559817000000     1.329225000000     3.400032000000
    H               -1.532189000000     2.820853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4__n1.6___', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N               -0.610939000000     1.197258000000     3.400058000000
    C               -0.382336000000    -0.193037000000     3.400004000000
    N               -1.467609000000    -1.008972000000     3.399965000000
    C               -2.710962000000    -0.503844000000     3.399971000000
    C               -2.952449000000     0.921627000000     3.400023000000
    C               -1.863556000000     1.733781000000     3.400063000000
    H                0.215200000000     1.780954000000     3.400092000000
    O                0.788816000000    -0.594965000000     3.399991000000
    N               -3.730425000000    -1.368706000000     3.399929000000
    H               -3.556075000000    -2.396574000000     3.399871000000
    H               -4.678136000000    -1.025139000000     3.399928000000
    H               -3.959817000000     1.329225000000     3.400032000000
    H               -1.932189000000     2.820853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4__n2.0___', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N               -1.010939000000     1.197258000000     3.400058000000
    C               -0.782336000000    -0.193037000000     3.400004000000
    N               -1.867609000000    -1.008972000000     3.399965000000
    C               -3.110962000000    -0.503844000000     3.399971000000
    C               -3.352449000000     0.921627000000     3.400023000000
    C               -2.263556000000     1.733781000000     3.400063000000
    H               -0.184800000000     1.780954000000     3.400092000000
    O                0.388816000000    -0.594965000000     3.399991000000
    N               -4.130425000000    -1.368706000000     3.399929000000
    H               -3.956075000000    -2.396574000000     3.399871000000
    H               -5.078136000000    -1.025139000000     3.399928000000
    H               -4.359817000000     1.329225000000     3.400032000000
    H               -2.332189000000     2.820853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4__n0.4___', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.589061000000     1.197258000000     3.400058000000
    C                0.817664000000    -0.193037000000     3.400004000000
    N               -0.267609000000    -1.008972000000     3.399965000000
    C               -1.510962000000    -0.503844000000     3.399971000000
    C               -1.752449000000     0.921627000000     3.400023000000
    C               -0.663556000000     1.733781000000     3.400063000000
    H                1.415200000000     1.780954000000     3.400092000000
    O                1.988816000000    -0.594965000000     3.399991000000
    N               -2.530425000000    -1.368706000000     3.399929000000
    H               -2.356075000000    -2.396574000000     3.399871000000
    H               -3.478136000000    -1.025139000000     3.399928000000
    H               -2.759817000000     1.329225000000     3.400032000000
    H               -0.732189000000     2.820853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4__n0.8___', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.189061000000     1.197258000000     3.400058000000
    C                0.417664000000    -0.193037000000     3.400004000000
    N               -0.667609000000    -1.008972000000     3.399965000000
    C               -1.910962000000    -0.503844000000     3.399971000000
    C               -2.152449000000     0.921627000000     3.400023000000
    C               -1.063556000000     1.733781000000     3.400063000000
    H                1.015200000000     1.780954000000     3.400092000000
    O                1.588816000000    -0.594965000000     3.399991000000
    N               -2.930425000000    -1.368706000000     3.399929000000
    H               -2.756075000000    -2.396574000000     3.399871000000
    H               -3.878136000000    -1.025139000000     3.399928000000
    H               -3.159817000000     1.329225000000     3.400032000000
    H               -1.132189000000     2.820853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4_n1.2____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000    -0.002742000000     3.400058000000
    C                1.217664000000    -1.393037000000     3.400004000000
    N                0.132391000000    -2.208972000000     3.399965000000
    C               -1.110962000000    -1.703844000000     3.399971000000
    C               -1.352449000000    -0.278373000000     3.400023000000
    C               -0.263556000000     0.533781000000     3.400063000000
    H                1.815200000000     0.580954000000     3.400092000000
    O                2.388816000000    -1.794965000000     3.399991000000
    N               -2.130425000000    -2.568706000000     3.399929000000
    H               -1.956075000000    -3.596574000000     3.399871000000
    H               -3.078136000000    -2.225139000000     3.399928000000
    H               -2.359817000000     0.129225000000     3.400032000000
    H               -0.332189000000     1.620853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4_n1.6____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000    -0.402742000000     3.400058000000
    C                1.217664000000    -1.793037000000     3.400004000000
    N                0.132391000000    -2.608972000000     3.399965000000
    C               -1.110962000000    -2.103844000000     3.399971000000
    C               -1.352449000000    -0.678373000000     3.400023000000
    C               -0.263556000000     0.133781000000     3.400063000000
    H                1.815200000000     0.180954000000     3.400092000000
    O                2.388816000000    -2.194965000000     3.399991000000
    N               -2.130425000000    -2.968706000000     3.399929000000
    H               -1.956075000000    -3.996574000000     3.399871000000
    H               -3.078136000000    -2.625139000000     3.399928000000
    H               -2.359817000000    -0.270775000000     3.400032000000
    H               -0.332189000000     1.220853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4_n2.0____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000    -0.802742000000     3.400058000000
    C                1.217664000000    -2.193037000000     3.400004000000
    N                0.132391000000    -3.008972000000     3.399965000000
    C               -1.110962000000    -2.503844000000     3.399971000000
    C               -1.352449000000    -1.078373000000     3.400023000000
    C               -0.263556000000    -0.266219000000     3.400063000000
    H                1.815200000000    -0.219046000000     3.400092000000
    O                2.388816000000    -2.594965000000     3.399991000000
    N               -2.130425000000    -3.368706000000     3.399929000000
    H               -1.956075000000    -4.396574000000     3.399871000000
    H               -3.078136000000    -3.025139000000     3.399928000000
    H               -2.359817000000    -0.670775000000     3.400032000000
    H               -0.332189000000     0.820853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4_n0.4____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     0.797258000000     3.400058000000
    C                1.217664000000    -0.593037000000     3.400004000000
    N                0.132391000000    -1.408972000000     3.399965000000
    C               -1.110962000000    -0.903844000000     3.399971000000
    C               -1.352449000000     0.521627000000     3.400023000000
    C               -0.263556000000     1.333781000000     3.400063000000
    H                1.815200000000     1.380954000000     3.400092000000
    O                2.388816000000    -0.994965000000     3.399991000000
    N               -2.130425000000    -1.768706000000     3.399929000000
    H               -1.956075000000    -2.796574000000     3.399871000000
    H               -3.078136000000    -1.425139000000     3.399928000000
    H               -2.359817000000     0.929225000000     3.400032000000
    H               -0.332189000000     2.420853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.4_n0.8____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     0.397258000000     3.400058000000
    C                1.217664000000    -0.993037000000     3.400004000000
    N                0.132391000000    -1.808972000000     3.399965000000
    C               -1.110962000000    -1.303844000000     3.399971000000
    C               -1.352449000000     0.121627000000     3.400023000000
    C               -0.263556000000     0.933781000000     3.400063000000
    H                1.815200000000     0.980954000000     3.400092000000
    O                2.388816000000    -1.394965000000     3.399991000000
    N               -2.130425000000    -2.168706000000     3.399929000000
    H               -1.956075000000    -3.196574000000     3.399871000000
    H               -3.078136000000    -1.825139000000     3.399928000000
    H               -2.359817000000     0.529225000000     3.400032000000
    H               -0.332189000000     2.020853000000     3.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6____12_', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.967460000000     1.197258000000     3.394419000000
    C                1.191056000000    -0.193037000000     3.346837000000
    N                0.129491000000    -1.008972000000     3.572440000000
    C               -1.086691000000    -0.503844000000     3.830954000000
    C               -1.322890000000     0.921627000000     3.881212000000
    C               -0.257784000000     1.733781000000     3.654858000000
    H                1.775553000000     1.780954000000     3.222689000000
    O                2.336613000000    -0.594965000000     3.103328000000
    N               -2.083885000000    -1.368706000000     4.042871000000
    H               -1.913357000000    -2.396574000000     4.006565000000
    H               -3.010886000000    -1.025139000000     4.239910000000
    H               -2.308243000000     1.329225000000     4.090665000000
    H               -0.324908000000     2.820853000000     3.669170000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6____16_', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.950762000000     1.197258000000     3.327434000000
    C                1.170495000000    -0.193037000000     3.264370000000
    N                0.127253000000    -1.008972000000     3.563474000000
    C               -1.067933000000    -0.503844000000     3.906195000000
    C               -1.300051000000     0.921627000000     3.972808000000
    C               -0.253329000000     1.733781000000     3.672706000000
    H                1.744908000000     1.780954000000     3.099752000000
    O                2.296275000000    -0.594965000000     2.941544000000
    N               -2.047916000000    -1.368706000000     4.187156000000
    H               -1.880336000000    -2.396574000000     4.139043000000
    H               -2.958914000000    -1.025139000000     4.448380000000
    H               -2.268393000000     1.329225000000     4.250484000000
    H               -0.319291000000     2.820853000000     3.691666000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6____20_', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.929433000000     1.197258000000     3.261776000000
    C                1.144231000000    -0.193037000000     3.183538000000
    N                0.124395000000    -1.008972000000     3.554687000000
    C               -1.043973000000    -0.503844000000     3.979944000000
    C               -1.270878000000     0.921627000000     4.062586000000
    C               -0.247640000000     1.733781000000     3.690201000000
    H                1.705762000000     1.780954000000     2.979251000000
    O                2.244750000000    -0.594965000000     2.782968000000
    N               -2.001969000000    -1.368706000000     4.328582000000
    H               -1.838153000000    -2.396574000000     4.268896000000
    H               -2.892526000000    -1.025139000000     4.652717000000
    H               -2.217492000000     1.329225000000     4.407135000000
    H               -0.312119000000     2.820853000000     3.713715000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6____4_', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.986656000000     1.197258000000     3.531064000000
    C                1.214698000000    -0.193037000000     3.515064000000
    N                0.132066000000    -1.008972000000     3.590730000000
    C               -1.108258000000    -0.503844000000     3.677468000000
    C               -1.349153000000     0.921627000000     3.694365000000
    C               -0.262910000000     1.733781000000     3.618448000000
    H                1.810785000000     1.780954000000     3.473470000000
    O                2.382996000000    -0.594965000000     3.433356000000
    N               -2.125240000000    -1.368706000000     3.748540000000
    H               -1.951319000000    -2.396574000000     3.736320000000
    H               -3.070643000000    -1.025139000000     3.814648000000
    H               -2.354066000000     1.329225000000     3.764644000000
    H               -0.331372000000     2.820853000000     3.623278000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6____8_', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.979444000000     1.197258000000     3.462407000000
    C                1.205814000000    -0.193037000000     3.430538000000
    N                0.131098000000    -1.008972000000     3.581540000000
    C               -1.100154000000    -0.503844000000     3.754587000000
    C               -1.339284000000     0.921627000000     3.788247000000
    C               -0.260982000000     1.733781000000     3.636742000000
    H                1.797547000000     1.780954000000     3.347464000000
    O                2.365567000000    -0.594965000000     3.267532000000
    N               -2.109702000000    -1.368706000000     3.896428000000
    H               -1.937057000000    -2.396574000000     3.872105000000
    H               -3.048190000000    -1.025139000000     4.028322000000
    H               -2.336847000000     1.329225000000     3.928455000000
    H               -0.328941000000     2.820853000000     3.646337000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6_____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.197258000000     3.600058000000
    C                1.217664000000    -0.193037000000     3.600004000000
    N                0.132391000000    -1.008972000000     3.599965000000
    C               -1.110962000000    -0.503844000000     3.599971000000
    C               -1.352449000000     0.921627000000     3.600023000000
    C               -0.263556000000     1.733781000000     3.600063000000
    H                1.815200000000     1.780954000000     3.600092000000
    O                2.388816000000    -0.594965000000     3.599991000000
    N               -2.130425000000    -1.368706000000     3.599929000000
    H               -1.956075000000    -2.396574000000     3.599871000000
    H               -3.078136000000    -1.025139000000     3.599928000000
    H               -2.359817000000     1.329225000000     3.600032000000
    H               -0.332189000000     2.820853000000     3.600106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6_____12', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.171083000000     3.848981000000
    C                1.217664000000    -0.188820000000     3.559869000000
    N                0.132391000000    -0.986916000000     3.390189000000
    C               -1.110962000000    -0.492828000000     3.495217000000
    C               -1.352449000000     0.901482000000     3.791640000000
    C               -0.263556000000     1.695881000000     3.960535000000
    H                1.815200000000     1.742017000000     3.970371000000
    O                2.388816000000    -0.581962000000     3.476291000000
    N               -2.130425000000    -1.338782000000     3.315361000000
    H               -1.956075000000    -2.344176000000     3.101598000000
    H               -3.078136000000    -1.002722000000     3.386791000000
    H               -2.359817000000     1.300172000000     3.876393000000
    H               -0.332189000000     2.759189000000     4.186592000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6_____16', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.150862000000     3.930065000000
    C                1.217664000000    -0.185560000000     3.546796000000
    N                0.132391000000    -0.969876000000     3.321856000000
    C               -1.110962000000    -0.484318000000     3.461094000000
    C               -1.352449000000     0.885918000000     3.854057000000
    C               -0.263556000000     1.666600000000     4.077955000000
    H                1.815200000000     1.711938000000     4.090986000000
    O                2.388816000000    -0.571915000000     3.435997000000
    N               -2.130425000000    -1.315665000000     3.222665000000
    H               -1.956075000000    -2.303699000000     2.939291000000
    H               -3.078136000000    -0.985407000000     3.317364000000
    H               -2.359817000000     1.277724000000     3.966415000000
    H               -0.332189000000     2.711549000000     4.377634000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6_____20', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.125035000000     4.009541000000
    C                1.217664000000    -0.181397000000     3.533981000000
    N                0.132391000000    -0.948112000000     3.254878000000
    C               -1.110962000000    -0.473449000000     3.427648000000
    C               -1.352449000000     0.866038000000     3.915237000000
    C               -0.263556000000     1.629200000000     4.193047000000
    H                1.815200000000     1.673518000000     4.209209000000
    O                2.388816000000    -0.559081000000     3.396502000000
    N               -2.130425000000    -1.286139000000     3.131808000000
    H               -1.956075000000    -2.251999000000     2.780202000000
    H               -3.078136000000    -0.963291000000     3.249314000000
    H               -2.359817000000     1.249052000000     4.054652000000
    H               -0.332189000000     2.650698000000     4.564888000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6_____4', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.194337000000     3.683574000000
    C                1.217664000000    -0.192567000000     3.586538000000
    N                0.132391000000    -1.006512000000     3.529583000000
    C               -1.110962000000    -0.502615000000     3.564825000000
    C               -1.352449000000     0.919380000000     3.664312000000
    C               -0.263556000000     1.729553000000     3.721005000000
    H                1.815200000000     1.776609000000     3.724325000000
    O                2.388816000000    -0.593515000000     3.558488000000
    N               -2.130425000000    -1.365367000000     3.504453000000
    H               -1.956075000000    -2.390727000000     3.432695000000
    H               -3.078136000000    -1.022637000000     3.528418000000
    H               -2.359817000000     1.325985000000     3.692754000000
    H               -0.332189000000     2.813974000000     3.796879000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6_____8', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.185598000000     3.766684000000
    C                1.217664000000    -0.191159000000     3.573138000000
    N                0.132391000000    -0.999148000000     3.459544000000
    C               -1.110962000000    -0.498937000000     3.529850000000
    C               -1.352449000000     0.912655000000     3.728288000000
    C               -0.263556000000     1.716899000000     3.841358000000
    H                1.815200000000     1.763609000000     3.847952000000
    O                2.388816000000    -0.589174000000     3.517188000000
    N               -2.130425000000    -1.355376000000     3.409443000000
    H               -1.956075000000    -2.373233000000     3.266334000000
    H               -3.078136000000    -1.015152000000     3.457257000000
    H               -2.359817000000     1.316285000000     3.785024000000
    H               -0.332189000000     2.793386000000     3.992692000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6_____n12', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.171107000000     3.351133000000
    C                1.217664000000    -0.188818000000     3.640139000000
    N                0.132391000000    -0.986931000000     3.809743000000
    C               -1.110962000000    -0.492840000000     3.704727000000
    C               -1.352449000000     0.901492000000     3.408405000000
    C               -0.263556000000     1.695907000000     3.239588000000
    H                1.815200000000     1.742055000000     3.229809000000
    O                2.388816000000    -0.581965000000     3.723691000000
    N               -2.130425000000    -1.338811000000     3.884501000000
    H               -1.956075000000    -2.344230000000     4.098150000000
    H               -3.078136000000    -1.002752000000     3.813068000000
    H               -2.359817000000     1.300185000000     3.323670000000
    H               -0.332189000000     2.759233000000     3.013615000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6_____n16', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.150894000000     3.270047000000
    C                1.217664000000    -0.185558000000     3.653212000000
    N                0.132391000000    -0.969896000000     3.878077000000
    C               -1.110962000000    -0.484334000000     3.738850000000
    C               -1.352449000000     0.885931000000     3.345987000000
    C               -0.263556000000     1.666635000000     3.122166000000
    H                1.815200000000     1.711988000000     3.109191000000
    O                2.388816000000    -0.571920000000     3.763986000000
    N               -2.130425000000    -1.315704000000     3.977198000000
    H               -1.956075000000    -2.303770000000     4.260461000000
    H               -3.078136000000    -0.985447000000     3.882497000000
    H               -2.359817000000     1.277742000000     3.233647000000
    H               -0.332189000000     2.711607000000     2.822569000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6_____n20', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.125074000000     3.190568000000
    C                1.217664000000    -0.181394000000     3.666026000000
    N                0.132391000000    -0.948136000000     3.945056000000
    C               -1.110962000000    -0.473468000000     3.772298000000
    C               -1.352449000000     0.866054000000     3.284807000000
    C               -0.263556000000     1.629243000000     3.007071000000
    H                1.815200000000     1.673581000000     2.990964000000
    O                2.388816000000    -0.559087000000     3.803482000000
    N               -2.130425000000    -1.286187000000     4.068058000000
    H               -1.956075000000    -2.252087000000     4.419555000000
    H               -3.078136000000    -0.963340000000     3.950551000000
    H               -2.359817000000     1.249074000000     3.145408000000
    H               -0.332189000000     2.650771000000     2.635311000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6_____n4', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.194346000000     3.516541000000
    C                1.217664000000    -0.192566000000     3.613470000000
    N                0.132391000000    -1.006517000000     3.670347000000
    C               -1.110962000000    -0.502619000000     3.635117000000
    C               -1.352449000000     0.919384000000     3.535733000000
    C               -0.263556000000     1.729562000000     3.479120000000
    H                1.815200000000     1.776622000000     3.475859000000
    O                2.388816000000    -0.593516000000     3.641494000000
    N               -2.130425000000    -1.365377000000     3.695405000000
    H               -1.956075000000    -2.390745000000     3.767048000000
    H               -3.078136000000    -1.022647000000     3.671438000000
    H               -2.359817000000     1.325989000000     3.507310000000
    H               -0.332189000000     2.813989000000     3.403333000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6_____n8', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.185614000000     3.433431000000
    C                1.217664000000    -0.191158000000     3.626870000000
    N                0.132391000000    -0.999158000000     3.740387000000
    C               -1.110962000000    -0.498945000000     3.670093000000
    C               -1.352449000000     0.912661000000     3.471757000000
    C               -0.263556000000     1.716917000000     3.358767000000
    H                1.815200000000     1.763635000000     3.352230000000
    O                2.388816000000    -0.589176000000     3.682794000000
    N               -2.130425000000    -1.355396000000     3.790417000000
    H               -1.956075000000    -2.373269000000     3.933411000000
    H               -3.078136000000    -1.015172000000     3.742600000000
    H               -2.359817000000     1.316294000000     3.415039000000
    H               -0.332189000000     2.793415000000     3.207518000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6____n12_', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.967436000000     1.197258000000     3.805694000000
    C                1.191054000000    -0.193037000000     3.853170000000
    N                0.129505000000    -1.008972000000     3.627491000000
    C               -1.086679000000    -0.503844000000     3.368990000000
    C               -1.322900000000     0.921627000000     3.318833000000
    C               -0.257810000000     1.733781000000     3.545265000000
    H                1.775514000000     1.780954000000     3.977491000000
    O                2.336617000000    -0.594965000000     4.096654000000
    N               -2.083855000000    -1.368706000000     3.156990000000
    H               -1.913303000000    -2.396574000000     3.193183000000
    H               -3.010856000000    -1.025139000000     2.959949000000
    H               -2.308256000000     1.329225000000     3.109398000000
    H               -0.324952000000     2.820853000000     3.531038000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6____n16_', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.950730000000     1.197258000000     3.872678000000
    C                1.170493000000    -0.193037000000     3.935638000000
    N                0.127272000000    -1.008972000000     3.636458000000
    C               -1.067917000000    -0.503844000000     3.293749000000
    C               -1.300064000000     0.921627000000     3.227237000000
    C               -0.253364000000     1.733781000000     3.527415000000
    H                1.744857000000     1.780954000000     4.100425000000
    O                2.296280000000    -0.594965000000     4.258438000000
    N               -2.047876000000    -1.368706000000     3.012707000000
    H               -1.880264000000    -2.396574000000     3.060709000000
    H               -2.958874000000    -1.025139000000     2.751482000000
    H               -2.268411000000     1.329225000000     2.949577000000
    H               -0.319350000000     2.820853000000     3.508538000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6____n20_', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.929393000000     1.197258000000     3.938333000000
    C                1.144229000000    -0.193037000000     4.016469000000
    N                0.124419000000    -1.008972000000     3.645247000000
    C               -1.043953000000    -0.503844000000     3.220001000000
    C               -1.270894000000     0.921627000000     3.137457000000
    C               -0.247683000000     1.733781000000     3.509918000000
    H                1.705699000000     1.780954000000     4.220921000000
    O                2.244756000000    -0.594965000000     4.417015000000
    N               -2.001920000000    -1.368706000000     2.871285000000
    H               -1.838065000000    -2.396574000000     2.930862000000
    H               -2.892477000000    -1.025139000000     2.547148000000
    H               -2.217514000000     1.329225000000     2.792925000000
    H               -0.312192000000     2.820853000000     3.486484000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6____n4_', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.986648000000     1.197258000000     3.669051000000
    C                1.214698000000    -0.193037000000     3.684944000000
    N                0.132071000000    -1.008972000000     3.609200000000
    C               -1.108254000000    -0.503844000000     3.522474000000
    C               -1.349156000000     0.921627000000     3.505681000000
    C               -0.262918000000     1.733781000000     3.581678000000
    H                1.810772000000     1.780954000000     3.726714000000
    O                2.382998000000    -0.594965000000     3.766626000000
    N               -2.125230000000    -1.368706000000     3.451318000000
    H               -1.951301000000    -2.396574000000     3.463422000000
    H               -3.070633000000    -1.025139000000     3.385208000000
    H               -2.354071000000     1.329225000000     3.435419000000
    H               -0.331387000000     2.820853000000     3.576933000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.6____n8_', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.979427000000     1.197258000000     3.737708000000
    C                1.205813000000    -0.193037000000     3.769470000000
    N                0.131107000000    -1.008972000000     3.618391000000
    C               -1.100146000000    -0.503844000000     3.445355000000
    C               -1.339290000000     0.921627000000     3.411798000000
    C               -0.261000000000     1.733781000000     3.563382000000
    H                1.797522000000     1.780954000000     3.852718000000
    O                2.365569000000    -0.594965000000     3.932450000000
    N               -2.109682000000    -1.368706000000     3.303432000000
    H               -1.937021000000    -2.396574000000     3.327639000000
    H               -3.048170000000    -1.025139000000     3.171535000000
    H               -2.336856000000     1.329225000000     3.271609000000
    H               -0.328971000000     2.820853000000     3.553873000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-3.8_____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.197258000000     3.800058000000
    C                1.217664000000    -0.193037000000     3.800004000000
    N                0.132391000000    -1.008972000000     3.799965000000
    C               -1.110962000000    -0.503844000000     3.799971000000
    C               -1.352449000000     0.921627000000     3.800023000000
    C               -0.263556000000     1.733781000000     3.800063000000
    H                1.815200000000     1.780954000000     3.800092000000
    O                2.388816000000    -0.594965000000     3.799991000000
    N               -2.130425000000    -1.368706000000     3.799929000000
    H               -1.956075000000    -2.396574000000     3.799871000000
    H               -3.078136000000    -1.025139000000     3.799928000000
    H               -2.359817000000     1.329225000000     3.800032000000
    H               -0.332189000000     2.820853000000     3.800106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-4.0_____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.197258000000     4.000058000000
    C                1.217664000000    -0.193037000000     4.000004000000
    N                0.132391000000    -1.008972000000     3.999965000000
    C               -1.110962000000    -0.503844000000     3.999971000000
    C               -1.352449000000     0.921627000000     4.000023000000
    C               -0.263556000000     1.733781000000     4.000063000000
    H                1.815200000000     1.780954000000     4.000092000000
    O                2.388816000000    -0.594965000000     3.999991000000
    N               -2.130425000000    -1.368706000000     3.999929000000
    H               -1.956075000000    -2.396574000000     3.999871000000
    H               -3.078136000000    -1.025139000000     3.999928000000
    H               -2.359817000000     1.329225000000     4.000032000000
    H               -0.332189000000     2.820853000000     4.000106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-4.4_____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.197258000000     4.400058000000
    C                1.217664000000    -0.193037000000     4.400004000000
    N                0.132391000000    -1.008972000000     4.399965000000
    C               -1.110962000000    -0.503844000000     4.399971000000
    C               -1.352449000000     0.921627000000     4.400023000000
    C               -0.263556000000     1.733781000000     4.400063000000
    H                1.815200000000     1.780954000000     4.400092000000
    O                2.388816000000    -0.594965000000     4.399991000000
    N               -2.130425000000    -1.368706000000     4.399929000000
    H               -1.956075000000    -2.396574000000     4.399871000000
    H               -3.078136000000    -1.025139000000     4.399928000000
    H               -2.359817000000     1.329225000000     4.400032000000
    H               -0.332189000000     2.820853000000     4.400106000000

""")

GEOS['%s-%s-%s' % (dbse, 'AC-5.0_____', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               -0.214386000000    -1.993021000000     0.000158000000
    C               -1.483458000000    -1.537649000000     0.000039000000
    N               -1.896987000000    -0.269765000000    -0.000126000000
    C               -0.859092000000     0.583332000000    -0.000130000000
    C                0.502351000000     0.258224000000     0.000025000000
    C                0.818987000000    -1.120415000000     0.000169000000
    N                1.293025000000     1.395362000000    -0.000018000000
    C                0.430884000000     2.385528000000    -0.000191000000
    N               -0.884630000000     1.960299000000    -0.000268000000
    H               -2.247335000000    -2.317004000000     0.000055000000
    N                2.078153000000    -1.590279000000     0.000313000000
    H                2.839337000000    -0.928535000000     0.000354000000
    H                2.257625000000    -2.599226000000     0.000507000000
    H                0.686428000000     3.440469000000    -0.000272000000
    H               -1.715093000000     2.535836000000    -0.000408000000
    --
    0 1
    N                0.989061000000     1.197258000000     5.000058000000
    C                1.217664000000    -0.193037000000     5.000004000000
    N                0.132391000000    -1.008972000000     4.999965000000
    C               -1.110962000000    -0.503844000000     4.999971000000
    C               -1.352449000000     0.921627000000     5.000023000000
    C               -0.263556000000     1.733781000000     5.000063000000
    H                1.815200000000     1.780954000000     5.000092000000
    O                2.388816000000    -0.594965000000     4.999991000000
    N               -2.130425000000    -1.368706000000     4.999929000000
    H               -1.956075000000    -2.396574000000     4.999871000000
    H               -3.078136000000    -1.025139000000     4.999928000000
    H               -2.359817000000     1.329225000000     5.000032000000
    H               -0.332189000000     2.820853000000     5.000106000000

""")

# <<< Derived Geometry Strings >>>
for rxn in HRXN:
    GEOS['%s-%s-monoA-unCP' % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1)
    GEOS['%s-%s-monoB-unCP' % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2)
    GEOS['%s-%s-monoA-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1, 2)
    GEOS['%s-%s-monoB-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2, 1)

#########################################################################

# <<< Supplementary Quantum Chemical Results >>>
DATA = {}

# SAPT0S-SA-jadz
DATA['SAPT ELST ENERGY'] = {}
DATA['SAPT ELST ENERGY']['ACHC-AC-3.0_____'] =   -11.1483
DATA['SAPT ELST ENERGY']['ACHC-AC-3.2_____'] =    -5.5807
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4_____'] =    -2.6372
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6_____'] =    -1.1011
DATA['SAPT ELST ENERGY']['ACHC-AC-3.8_____'] =    -0.3125
DATA['SAPT ELST ENERGY']['ACHC-AC-4.0_____'] =     0.0827
DATA['SAPT ELST ENERGY']['ACHC-AC-4.4_____'] =     0.3568
DATA['SAPT ELST ENERGY']['ACHC-AC-5.0_____'] =     0.3804
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4_n2.0____'] =    -3.4418
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4_n1.6____'] =    -3.2686
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4_n1.2____'] =    -3.0860
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4_n0.8____'] =    -2.9389
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4_n0.4____'] =    -2.8105
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4_0.4____'] =    -2.3464
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4_0.8____'] =    -1.9075
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4_1.2____'] =    -1.3601
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4_1.6____'] =    -0.7859
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4_2.0____'] =    -0.2572
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4__n2.0___'] =     0.5043
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4__n1.6___'] =    -0.0676
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4__n1.2___'] =    -0.6574
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4__n0.8___'] =    -1.2457
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4__n0.4___'] =    -1.8909
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4__0.4___'] =    -3.4229
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4__0.8___'] =    -4.0740
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4__1.2___'] =    -4.4120
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4__1.6___'] =    -4.3818
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4__2.0___'] =    -4.0795
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4___30__'] =    -1.8667
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4___60__'] =    -1.6087
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4___90__'] =    -3.9821
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4___120__'] =    -4.6142
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4___150__'] =    -5.0397
DATA['SAPT ELST ENERGY']['ACHC-AC-3.4___180__'] =    -3.5749
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6____4_'] =    -0.9018
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6____8_'] =    -0.8909
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6____12_'] =    -1.0935
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6____16_'] =    -1.5752
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6____20_'] =    -2.4526
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6____n4_'] =    -1.4972
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6____n8_'] =    -2.1318
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6____n12_'] =    -3.0837
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6____n16_'] =    -4.4792
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6____n20_'] =    -6.5050
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6_____n20'] =    -4.0502
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6_____n16'] =    -2.9317
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6_____n12'] =    -2.1510
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6_____n8'] =    -1.6215
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6_____n4'] =    -1.2846
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6_____4'] =    -1.0458
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6_____8'] =    -1.1042
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6_____12'] =    -1.2707
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6_____16'] =    -1.5473
DATA['SAPT ELST ENERGY']['ACHC-AC-3.6_____20'] =    -1.9429
DATA['SAPT EXCH ENERGY'] = {}
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.0_____'] =    29.8907
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.2_____'] =    15.8159
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4_____'] =     8.3132
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6_____'] =     4.3463
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.8_____'] =     2.2599
DATA['SAPT EXCH ENERGY']['ACHC-AC-4.0_____'] =     1.1666
DATA['SAPT EXCH ENERGY']['ACHC-AC-4.4_____'] =     0.3005
DATA['SAPT EXCH ENERGY']['ACHC-AC-5.0_____'] =     0.0347
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4_n2.0____'] =     6.6598
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4_n1.6____'] =     7.5900
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4_n1.2____'] =     8.2348
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4_n0.8____'] =     8.4123
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4_n0.4____'] =     8.3161
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4_0.4____'] =     8.4787
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4_0.8____'] =     8.4679
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4_1.2____'] =     7.9023
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4_1.6____'] =     6.7918
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4_2.0____'] =     5.5093
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4__n2.0___'] =     4.7678
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4__n1.6___'] =     5.4155
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4__n1.2___'] =     6.0419
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4__n0.8___'] =     6.8025
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4__n0.4___'] =     7.6555
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4__0.4___'] =     8.5175
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4__0.8___'] =     8.3062
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4__1.2___'] =     7.9305
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4__1.6___'] =     7.5405
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4__2.0___'] =     7.0462
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4___30__'] =     8.1036
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4___60__'] =     8.7155
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4___90__'] =     8.8794
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4___120__'] =     9.2743
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4___150__'] =     8.1433
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.4___180__'] =     8.2616
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6____4_'] =     4.1620
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6____8_'] =     4.3834
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6____12_'] =     5.0716
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6____16_'] =     6.3743
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6____20_'] =     8.5532
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6____n4_'] =     4.9486
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6____n8_'] =     6.0547
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6____n12_'] =     7.8363
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6____n16_'] =    10.5781
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6____n20_'] =    14.7215
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6_____n20'] =    11.4364
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6_____n16'] =     8.6205
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6_____n12'] =     6.7170
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6_____n8'] =     5.4700
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6_____n4'] =     4.7139
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6_____4'] =     4.3095
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6_____8'] =     4.5797
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6_____12'] =     5.1606
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6_____16'] =     6.0809
DATA['SAPT EXCH ENERGY']['ACHC-AC-3.6_____20'] =     7.3945
DATA['SAPT IND ENERGY'] = {}
DATA['SAPT IND ENERGY']['ACHC-AC-3.0_____'] =    -3.0264
DATA['SAPT IND ENERGY']['ACHC-AC-3.2_____'] =    -1.9747
DATA['SAPT IND ENERGY']['ACHC-AC-3.4_____'] =    -1.3189
DATA['SAPT IND ENERGY']['ACHC-AC-3.6_____'] =    -0.9201
DATA['SAPT IND ENERGY']['ACHC-AC-3.8_____'] =    -0.6735
DATA['SAPT IND ENERGY']['ACHC-AC-4.0_____'] =    -0.5132
DATA['SAPT IND ENERGY']['ACHC-AC-4.4_____'] =    -0.3215
DATA['SAPT IND ENERGY']['ACHC-AC-5.0_____'] =    -0.1774
DATA['SAPT IND ENERGY']['ACHC-AC-3.4_n2.0____'] =    -1.0633
DATA['SAPT IND ENERGY']['ACHC-AC-3.4_n1.6____'] =    -1.1476
DATA['SAPT IND ENERGY']['ACHC-AC-3.4_n1.2____'] =    -1.2205
DATA['SAPT IND ENERGY']['ACHC-AC-3.4_n0.8____'] =    -1.3018
DATA['SAPT IND ENERGY']['ACHC-AC-3.4_n0.4____'] =    -1.3542
DATA['SAPT IND ENERGY']['ACHC-AC-3.4_0.4____'] =    -1.2025
DATA['SAPT IND ENERGY']['ACHC-AC-3.4_0.8____'] =    -1.0951
DATA['SAPT IND ENERGY']['ACHC-AC-3.4_1.2____'] =    -1.0738
DATA['SAPT IND ENERGY']['ACHC-AC-3.4_1.6____'] =    -1.1093
DATA['SAPT IND ENERGY']['ACHC-AC-3.4_2.0____'] =    -1.1104
DATA['SAPT IND ENERGY']['ACHC-AC-3.4__n2.0___'] =    -1.2374
DATA['SAPT IND ENERGY']['ACHC-AC-3.4__n1.6___'] =    -1.3211
DATA['SAPT IND ENERGY']['ACHC-AC-3.4__n1.2___'] =    -1.3976
DATA['SAPT IND ENERGY']['ACHC-AC-3.4__n0.8___'] =    -1.4264
DATA['SAPT IND ENERGY']['ACHC-AC-3.4__n0.4___'] =    -1.3858
DATA['SAPT IND ENERGY']['ACHC-AC-3.4__0.4___'] =    -1.2915
DATA['SAPT IND ENERGY']['ACHC-AC-3.4__0.8___'] =    -1.3006
DATA['SAPT IND ENERGY']['ACHC-AC-3.4__1.2___'] =    -1.2765
DATA['SAPT IND ENERGY']['ACHC-AC-3.4__1.6___'] =    -1.1757
DATA['SAPT IND ENERGY']['ACHC-AC-3.4__2.0___'] =    -1.0255
DATA['SAPT IND ENERGY']['ACHC-AC-3.4___30__'] =    -1.2941
DATA['SAPT IND ENERGY']['ACHC-AC-3.4___60__'] =    -1.1950
DATA['SAPT IND ENERGY']['ACHC-AC-3.4___90__'] =    -1.4091
DATA['SAPT IND ENERGY']['ACHC-AC-3.4___120__'] =    -1.4116
DATA['SAPT IND ENERGY']['ACHC-AC-3.4___150__'] =    -1.3259
DATA['SAPT IND ENERGY']['ACHC-AC-3.4___180__'] =    -1.2483
DATA['SAPT IND ENERGY']['ACHC-AC-3.6____4_'] =    -0.9913
DATA['SAPT IND ENERGY']['ACHC-AC-3.6____8_'] =    -1.1136
DATA['SAPT IND ENERGY']['ACHC-AC-3.6____12_'] =    -1.2950
DATA['SAPT IND ENERGY']['ACHC-AC-3.6____16_'] =    -1.5494
DATA['SAPT IND ENERGY']['ACHC-AC-3.6____20_'] =    -1.8952
DATA['SAPT IND ENERGY']['ACHC-AC-3.6____n4_'] =    -0.8984
DATA['SAPT IND ENERGY']['ACHC-AC-3.6____n8_'] =    -0.9308
DATA['SAPT IND ENERGY']['ACHC-AC-3.6____n12_'] =    -1.0287
DATA['SAPT IND ENERGY']['ACHC-AC-3.6____n16_'] =    -1.2099
DATA['SAPT IND ENERGY']['ACHC-AC-3.6____n20_'] =    -1.4974
DATA['SAPT IND ENERGY']['ACHC-AC-3.6_____n20'] =    -1.4786
DATA['SAPT IND ENERGY']['ACHC-AC-3.6_____n16'] =    -1.2490
DATA['SAPT IND ENERGY']['ACHC-AC-3.6_____n12'] =    -1.0939
DATA['SAPT IND ENERGY']['ACHC-AC-3.6_____n8'] =    -0.9944
DATA['SAPT IND ENERGY']['ACHC-AC-3.6_____n4'] =    -0.9390
DATA['SAPT IND ENERGY']['ACHC-AC-3.6_____4'] =    -0.9335
DATA['SAPT IND ENERGY']['ACHC-AC-3.6_____8'] =    -0.9769
DATA['SAPT IND ENERGY']['ACHC-AC-3.6_____12'] =    -1.0506
DATA['SAPT IND ENERGY']['ACHC-AC-3.6_____16'] =    -1.1576
DATA['SAPT IND ENERGY']['ACHC-AC-3.6_____20'] =    -1.3047
DATA['SAPT DISP ENERGY'] = {}
DATA['SAPT DISP ENERGY']['ACHC-AC-3.0_____'] =   -18.8219
DATA['SAPT DISP ENERGY']['ACHC-AC-3.2_____'] =   -13.9170
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4_____'] =   -10.3667
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6_____'] =    -7.7868
DATA['SAPT DISP ENERGY']['ACHC-AC-3.8_____'] =    -5.9018
DATA['SAPT DISP ENERGY']['ACHC-AC-4.0_____'] =    -4.5149
DATA['SAPT DISP ENERGY']['ACHC-AC-4.4_____'] =    -2.7170
DATA['SAPT DISP ENERGY']['ACHC-AC-5.0_____'] =    -1.3555
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4_n2.0____'] =    -8.4427
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4_n1.6____'] =    -9.1766
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4_n1.2____'] =    -9.7427
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4_n0.8____'] =   -10.1596
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4_n0.4____'] =   -10.3811
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4_0.4____'] =   -10.1367
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4_0.8____'] =    -9.7193
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4_1.2____'] =    -9.1474
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4_1.6____'] =    -8.4538
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4_2.0____'] =    -7.6519
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4__n2.0___'] =    -7.2773
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4__n1.6___'] =    -8.1449
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4__n1.2___'] =    -8.9401
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4__n0.8___'] =    -9.6107
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4__n0.4___'] =   -10.1016
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4__0.4___'] =   -10.3896
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4__0.8___'] =   -10.1754
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4__1.2___'] =    -9.7367
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4__1.6___'] =    -9.1010
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4__2.0___'] =    -8.3173
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4___30__'] =   -10.3093
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4___60__'] =   -10.3641
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4___90__'] =   -10.5414
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4___120__'] =   -10.4560
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4___150__'] =   -10.1888
DATA['SAPT DISP ENERGY']['ACHC-AC-3.4___180__'] =   -10.1853
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6____4_'] =    -7.6820
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6____8_'] =    -7.7311
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6____12_'] =    -7.9390
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6____16_'] =    -8.3191
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6____20_'] =    -8.8946
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6____n4_'] =    -8.0484
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6____n8_'] =    -8.4772
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6____n12_'] =    -9.0916
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6____n16_'] =    -9.9205
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6____n20_'] =   -11.0063
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6_____n20'] =   -10.2955
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6_____n16'] =    -9.4471
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6_____n12'] =    -8.7982
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6_____n8'] =    -8.3188
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6_____n4'] =    -7.9870
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6_____4'] =    -7.7069
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6_____8'] =    -7.7399
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6_____12'] =    -7.8818
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6_____16'] =    -8.1313
DATA['SAPT DISP ENERGY']['ACHC-AC-3.6_____20'] =    -8.4897
