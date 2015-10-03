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

"""
| Database (Hobza) of interaction energies for dissociation curves of bimolecular complexes.
| Geometries and reference interaction energies from Grafova et al. JCTC 6 2365 (2010).
| Note that the S22by5-N-1.0 members are essentially the same geometries as S22-N (there's trivial round-off error) but the reference interaction energies for S22by5 are of lower quality than those of S22.
| Extended to 0.7 and 0.8 by Daniel G. A. Smith. Reference upgraded uniformly to CCSD(T)/[aTQZ;aDZ] TODO confirm by same.

- **cp**  ``'off'`` || ``'on'``

- **rlxd** ``'off'``

- **subset**

  - ``'small'``
  - ``'large'``
  - ``'equilibrium'``
  - ``'mol1'`` five-point (0.9, 1.0, 1.2, 1.5, 2.0) :math:`\\times R_{eq}` dissociation curve for molecule 1
  - ...
  - ``'mol22'`` five-point (0.9, 1.0, 1.2, 1.5, 2.0) :math:`\\times R_{eq}` dissociation curve for molecule 22

"""
import re
import qcdb

# <<< S22by5 Database Module >>>
dbse = 'S22by7'

# <<< Database Members >>>
AXIS_Rrat = {}

dist = [0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0]
mol1 = ['1-' + str(d) for d in dist]
mol2 = ['2-' + str(d) for d in dist]
mol3 = ['3-' + str(d) for d in dist]
mol4 = ['4-' + str(d) for d in dist]
mol5 = ['5-' + str(d) for d in dist]
mol6 = ['6-' + str(d) for d in dist]
mol7 = ['7-' + str(d) for d in dist]
mol8 = ['8-' + str(d) for d in dist]
mol9 = ['9-' + str(d) for d in dist]
mol10 = ['10-' + str(d) for d in dist]
mol11 = ['11-' + str(d) for d in dist]
mol12 = ['12-' + str(d) for d in dist]
mol13 = ['13-' + str(d) for d in dist]
mol14 = ['14-' + str(d) for d in dist]
mol15 = ['15-' + str(d) for d in dist]
mol16 = ['16-' + str(d) for d in dist]
mol17 = ['17-' + str(d) for d in dist]
mol18 = ['18-' + str(d) for d in dist]
mol19 = ['19-' + str(d) for d in dist]
mol20 = ['20-' + str(d) for d in dist]
mol21 = ['21-' + str(d) for d in dist]
mol22 = ['22-' + str(d) for d in dist]
AXIS_Rrat.update(dict(zip(mol1, dist)))
AXIS_Rrat.update(dict(zip(mol2, dist)))
AXIS_Rrat.update(dict(zip(mol3, dist)))
AXIS_Rrat.update(dict(zip(mol4, dist)))
AXIS_Rrat.update(dict(zip(mol5, dist)))
AXIS_Rrat.update(dict(zip(mol6, dist)))
AXIS_Rrat.update(dict(zip(mol7, dist)))
AXIS_Rrat.update(dict(zip(mol8, dist)))
AXIS_Rrat.update(dict(zip(mol9, dist)))
AXIS_Rrat.update(dict(zip(mol10, dist)))
AXIS_Rrat.update(dict(zip(mol11, dist)))
AXIS_Rrat.update(dict(zip(mol12, dist)))
AXIS_Rrat.update(dict(zip(mol13, dist)))
AXIS_Rrat.update(dict(zip(mol14, dist)))
AXIS_Rrat.update(dict(zip(mol15, dist)))
AXIS_Rrat.update(dict(zip(mol16, dist)))
AXIS_Rrat.update(dict(zip(mol17, dist)))
AXIS_Rrat.update(dict(zip(mol18, dist)))
AXIS_Rrat.update(dict(zip(mol19, dist)))
AXIS_Rrat.update(dict(zip(mol20, dist)))
AXIS_Rrat.update(dict(zip(mol21, dist)))
AXIS_Rrat.update(dict(zip(mol22, dist)))

temp = [mol1, mol2, mol3, mol4,  mol5, mol6, mol7, mol8, mol9, mol10, mol11,
        mol12, mol13, mol14, mol15,  mol16, mol17, mol18, mol19, mol20, mol21, mol22]
HRXN = sum(temp, [])

HRXN_S22BY5 = [rxn for rxn in HRXN if any(rxn.endswith(d) for d in ('0.9', '1.0', '1.2', '1.5', '2.0'))]
HRXN_SM = ['1-0.9', '2-1.0', '8-1.5', '16-2.0']
HRXN_LG = ['15-0.9']
HRXN_EQ = [str(m) + '-1.0' for m in range(1, 23)]
HB = sum([mol1, mol2, mol3, mol4, mol5, mol6, mol7], [])
MX = sum([mol13, mol15, mol16, mol17, mol18, mol19, mol21, mol22], [])
DD = sum([mol8, mol9, mol10, mol11, mol12, mol14, mol20], [])

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV_CP = {}  # order of active reagents per counterpoise-corrected reaction
ACTV_SA = {}  # order of active reagents for non-supramolecular calculations

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

# <<< Reference Values >>>
BIND = {}
# TODO by7
BIND['%s-1-0.7'  % (dbse)] =   None
BIND['%s-1-0.8'  % (dbse)] =   None
BIND['%s-1-0.9'  % (dbse)] =  -2.41
BIND['%s-1-1.0'  % (dbse)] =  -3.14
BIND['%s-1-1.2'  % (dbse)] =  -2.36
BIND['%s-1-1.5'  % (dbse)] =  -1.11
BIND['%s-1-2.0'  % (dbse)] =  -0.36
BIND['%s-2-0.7'  % (dbse)] =   None
BIND['%s-2-0.8'  % (dbse)] =   None
BIND['%s-2-0.9'  % (dbse)] =  -4.32
BIND['%s-2-1.0'  % (dbse)] =  -4.97
BIND['%s-2-1.2'  % (dbse)] =  -4.04
BIND['%s-2-1.5'  % (dbse)] =  -2.29
BIND['%s-2-2.0'  % (dbse)] =  -0.96
BIND['%s-3-0.7'  % (dbse)] =   None
BIND['%s-3-0.8'  % (dbse)] =   None
BIND['%s-3-0.9'  % (dbse)] = -16.34
BIND['%s-3-1.0'  % (dbse)] = -18.59
BIND['%s-3-1.2'  % (dbse)] = -15.62
BIND['%s-3-1.5'  % (dbse)] =  -9.24
BIND['%s-3-2.0'  % (dbse)] =  -3.63
BIND['%s-4-0.7'  % (dbse)] =   None
BIND['%s-4-0.8'  % (dbse)] =   None
BIND['%s-4-0.9'  % (dbse)] = -14.14
BIND['%s-4-1.0'  % (dbse)] = -15.95
BIND['%s-4-1.2'  % (dbse)] = -13.40
BIND['%s-4-1.5'  % (dbse)] =  -8.10
BIND['%s-4-2.0'  % (dbse)] =  -3.51
BIND['%s-5-0.7'  % (dbse)] =   None
BIND['%s-5-0.8'  % (dbse)] =   None
BIND['%s-5-0.9'  % (dbse)] = -18.73
BIND['%s-5-1.0'  % (dbse)] = -20.46
BIND['%s-5-1.2'  % (dbse)] = -17.16
BIND['%s-5-1.5'  % (dbse)] = -10.46
BIND['%s-5-2.0'  % (dbse)] =  -4.58
BIND['%s-6-0.7'  % (dbse)] =   None
BIND['%s-6-0.8'  % (dbse)] =   None
BIND['%s-6-0.9'  % (dbse)] = -15.13
BIND['%s-6-1.0'  % (dbse)] = -16.70
BIND['%s-6-1.2'  % (dbse)] = -13.93
BIND['%s-6-1.5'  % (dbse)] =  -8.18
BIND['%s-6-2.0'  % (dbse)] =  -3.26
BIND['%s-7-0.7'  % (dbse)] =   None
BIND['%s-7-0.8'  % (dbse)] =   None
BIND['%s-7-0.9'  % (dbse)] = -15.02
BIND['%s-7-1.0'  % (dbse)] = -16.37
BIND['%s-7-1.2'  % (dbse)] = -13.30
BIND['%s-7-1.5'  % (dbse)] =  -7.43
BIND['%s-7-2.0'  % (dbse)] =  -2.59
BIND['%s-8-0.7'  % (dbse)] =   None
BIND['%s-8-0.8'  % (dbse)] =   None
BIND['%s-8-0.9'  % (dbse)] =  -0.34
BIND['%s-8-1.0'  % (dbse)] =  -0.53
BIND['%s-8-1.2'  % (dbse)] =  -0.25
BIND['%s-8-1.5'  % (dbse)] =  -0.06
BIND['%s-8-2.0'  % (dbse)] =  -0.01
BIND['%s-9-0.7'  % (dbse)] =   None
BIND['%s-9-0.8'  % (dbse)] =   None
BIND['%s-9-0.9'  % (dbse)] =  -0.68
BIND['%s-9-1.0'  % (dbse)] =  -1.48
BIND['%s-9-1.2'  % (dbse)] =  -0.81
BIND['%s-9-1.5'  % (dbse)] =  -0.20
BIND['%s-9-2.0'  % (dbse)] =  -0.03
BIND['%s-10-0.7' % (dbse)] =   None
BIND['%s-10-0.8' % (dbse)] =   None
BIND['%s-10-0.9' % (dbse)] =  -1.09
BIND['%s-10-1.0' % (dbse)] =  -1.50
BIND['%s-10-1.2' % (dbse)] =  -1.13
BIND['%s-10-1.5' % (dbse)] =  -0.48
BIND['%s-10-2.0' % (dbse)] =  -0.12
BIND['%s-11-0.7' % (dbse)] =   None
BIND['%s-11-0.8' % (dbse)] =   None
BIND['%s-11-0.9' % (dbse)] =  -0.15
BIND['%s-11-1.0' % (dbse)] =  -2.81
BIND['%s-11-1.2' % (dbse)] =  -1.92
BIND['%s-11-1.5' % (dbse)] =  -0.53
BIND['%s-11-2.0' % (dbse)] =  -0.07
BIND['%s-12-0.7' % (dbse)] =   None
BIND['%s-12-0.8' % (dbse)] =   None
BIND['%s-12-0.9' % (dbse)] =  -1.69
BIND['%s-12-1.0' % (dbse)] =  -4.51
BIND['%s-12-1.2' % (dbse)] =  -3.02
BIND['%s-12-1.5' % (dbse)] =  -0.98
BIND['%s-12-2.0' % (dbse)] =  -0.19
BIND['%s-13-0.7' % (dbse)] =   None
BIND['%s-13-0.8' % (dbse)] =   None
BIND['%s-13-0.9' % (dbse)] =  -6.76
BIND['%s-13-1.0' % (dbse)] =  -9.87
BIND['%s-13-1.2' % (dbse)] =  -6.26
BIND['%s-13-1.5' % (dbse)] =  -2.42
BIND['%s-13-2.0' % (dbse)] =  -0.69
BIND['%s-14-0.7' % (dbse)] =   None
BIND['%s-14-0.8' % (dbse)] =   None
BIND['%s-14-0.9' % (dbse)] =  -2.13
BIND['%s-14-1.0' % (dbse)] =  -5.18
BIND['%s-14-1.2' % (dbse)] =  -3.61
BIND['%s-14-1.5' % (dbse)] =  -1.08
BIND['%s-14-2.0' % (dbse)] =  -0.10
BIND['%s-15-0.7' % (dbse)] =   None
BIND['%s-15-0.8' % (dbse)] =   None
BIND['%s-15-0.9' % (dbse)] =  -7.99
BIND['%s-15-1.0' % (dbse)] = -12.22
BIND['%s-15-1.2' % (dbse)] =  -8.23
BIND['%s-15-1.5' % (dbse)] =  -3.25
BIND['%s-15-2.0' % (dbse)] =  -0.92
BIND['%s-16-0.7' % (dbse)] =   None
BIND['%s-16-0.8' % (dbse)] =   None
BIND['%s-16-0.9' % (dbse)] =  -1.17
BIND['%s-16-1.0' % (dbse)] =  -1.49
BIND['%s-16-1.2' % (dbse)] =  -1.08
BIND['%s-16-1.5' % (dbse)] =  -0.49
BIND['%s-16-2.0' % (dbse)] =  -0.15
BIND['%s-17-0.7' % (dbse)] =   None
BIND['%s-17-0.8' % (dbse)] =   None
BIND['%s-17-0.9' % (dbse)] =  -3.01
BIND['%s-17-1.0' % (dbse)] =  -3.27
BIND['%s-17-1.2' % (dbse)] =  -2.47
BIND['%s-17-1.5' % (dbse)] =  -1.30
BIND['%s-17-2.0' % (dbse)] =  -0.49
BIND['%s-18-0.7' % (dbse)] =   None
BIND['%s-18-0.8' % (dbse)] =   None
BIND['%s-18-0.9' % (dbse)] =  -2.04
BIND['%s-18-1.0' % (dbse)] =  -2.35
BIND['%s-18-1.2' % (dbse)] =  -1.75
BIND['%s-18-1.5' % (dbse)] =  -0.85
BIND['%s-18-2.0' % (dbse)] =  -0.28
BIND['%s-19-0.7' % (dbse)] =   None
BIND['%s-19-0.8' % (dbse)] =   None
BIND['%s-19-0.9' % (dbse)] =  -4.02
BIND['%s-19-1.0' % (dbse)] =  -4.52
BIND['%s-19-1.2' % (dbse)] =  -3.68
BIND['%s-19-1.5' % (dbse)] =  -2.09
BIND['%s-19-2.0' % (dbse)] =  -0.85
BIND['%s-20-0.7' % (dbse)] =   None
BIND['%s-20-0.8' % (dbse)] =   None
BIND['%s-20-0.9' % (dbse)] =  -2.20
BIND['%s-20-1.0' % (dbse)] =  -2.80
BIND['%s-20-1.2' % (dbse)] =  -2.25
BIND['%s-20-1.5' % (dbse)] =  -1.12
BIND['%s-20-2.0' % (dbse)] =  -0.35
BIND['%s-21-0.7' % (dbse)] =   None
BIND['%s-21-0.8' % (dbse)] =   None
BIND['%s-21-0.9' % (dbse)] =  -4.99
BIND['%s-21-1.0' % (dbse)] =  -5.74
BIND['%s-21-1.2' % (dbse)] =  -4.88
BIND['%s-21-1.5' % (dbse)] =  -2.80
BIND['%s-21-2.0' % (dbse)] =  -1.10
BIND['%s-22-0.7' % (dbse)] =   None
BIND['%s-22-0.8' % (dbse)] =   None
BIND['%s-22-0.9' % (dbse)] =  -6.42
BIND['%s-22-1.0' % (dbse)] =  -7.05
BIND['%s-22-1.2' % (dbse)] =  -5.79
BIND['%s-22-1.5' % (dbse)] =  -3.41
BIND['%s-22-2.0' % (dbse)] =  -1.38

# <<< Comment Lines >>>
TAGL = {}
rxnpattern = re.compile(r'^(.+)-(.+)$')
for item in mol1:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'HB-1 Ammonia Dimer at %s Req, C2H' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Ammonia Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Ammonia from Ammonia Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Ammonia from Ammonia Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Ammonia from Ammonia Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Ammonia from Ammonia Dimer at %s Req' % (molname.group(2))

for item in mol2:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'HB-2 Water Dimer at %s Req, CS' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Water Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Water from Water Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Water from Water Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Water from Water Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Water from Water Dimer at %s Req' % (molname.group(2))

for item in mol3:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'HB-3 Formic Acid Dimer at %s Req, C2H' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Formic Acid Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Formic Acid from Formic Acid Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Formic Acid from Formic Acid Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Formic Acid from Formic Acid Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Formic Acid from Formic Acid Dimer at %s Req' % (molname.group(2))

for item in mol4:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'HB-4 Formamide Dimer at %s Req, C2H' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Formamide Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Formamide from Formamide Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Formamide from Formamide Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Formamide from Formamide Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Formamide from Formamide Dimer at %s Req' % (molname.group(2))

for item in mol5:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'HB-5 Uracil Dimer HB at %s Req, C2H' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Uracil Dimer HB at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Uracil from Uracil Dimer HB at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Uracil from Uracil Dimer HB at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Uracil from Uracil Dimer HB at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Uracil from Uracil Dimer HB at %s Req' % (molname.group(2))

for item in mol6:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'HB-6 2-Pyridone-2-Aminopyridine Complex at %s Req, C1' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      '2-Pyridone-2-Aminopyridine Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      '2-Pyridone from 2-Pyridone-2-Aminopyridine Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      '2-Aminopyridine from 2-Pyridone-2-Aminopyridine Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      '2-Pyridone from 2-Pyridone-2-Aminopyridine Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      '2-Aminopyridine from 2-Pyridone-2-Aminopyridine Complex at %s Req' % (molname.group(2))

for item in mol7:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'HB-7 Adenine-Thymine Complex WC at %s Req, C1' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Adenine-Thymine Complex WC at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Adenine from Adenine-Thymine Complex WC at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Thymine from Adenine-Thymine Complex WC at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Adenine from Adenine-Thymine Complex WC at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Thymine from Adenine-Thymine Complex WC at %s Req' % (molname.group(2))

for item in mol8:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'DD-1 Methane Dimer at %s Req, D3D' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Methane Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Methane from Methane Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Methane from Methane Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Methane from Methane Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Methane from Methane Dimer at %s Req' % (molname.group(2))

for item in mol9:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'DD-2 Ethene Dimer at %s Req, D2D' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Ethene Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Ethene from Ethene Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Ethene from Ethene Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Ethene from Ethene Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Ethene from Ethene Dimer at %s Req' % (molname.group(2))

for item in mol10:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'DD-3 Benzene-Methane Complex at %s Req, C3' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Benzene-Methane Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Benzene-Methane Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Methane from Benzene-Methane Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Benzene-Methane Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Methane from Benzene-Methane Complex at %s Req' % (molname.group(2))

for item in mol11:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'DD-4 Benzene Dimer Parallel-Disp at %s Req, C2H' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Benzene Dimer PD at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Benzene Dimer PD at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Benzene from Benzene Dimer PD at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Benzene Dimer PD at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Benzene from Benzene Dimer PD at %s Req' % (molname.group(2))

for item in mol12:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'DD-6 Pyrazine Dimer at %s Req, CS' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Pyrazine Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Pyrazine from Pyrazine Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Pyrazine from Pyrazine Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Pyrazine from Pyrazine Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Pyrazine from Pyrazine Dimer at %s Req' % (molname.group(2))

for item in mol13:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'MX-5 Uracil Dimer Stack at %s Req, C2' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Uracil Dimer Stack at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Uracil from Uracil Dimer Stack at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Uracil from Uracil Dimer Stack at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Uracil from Uracil Dimer Stack at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Uracil from Uracil Dimer Stack at %s Req' % (molname.group(2))

for item in mol14:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'DD-7 Indole-Benzene Complex Stack at %s Req, C1' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Indole-Benzene Complex Stack at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Indole-Benzene Complex Stack at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Indole from Indole-Benzene Complex Stack at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Indole-Benzene Complex Stack at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Indole from Indole-Benzene Complex Stack at %s Req' % (molname.group(2))

for item in mol15:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'MX-8 Adenine-Thymine Complex Stack at %s Req, C1' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Adenine-Thymine Complex Stack at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Adenine from Adenine-Thymine Complex Stack at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Thymine from Adenine-Thymine Complex Stack at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Adenine from Adenine-Thymine Complex Stack at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Thymine from Adenine-Thymine Complex Stack at %s Req' % (molname.group(2))

for item in mol16:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'MX-1 Ethene-Ethine Complex at %s Req, C2V' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Ethene-Ethine Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Ethene from Ethene-Ethine Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Ethine from Ethene-Ethine Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Ethene from Ethene-Ethine Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Ethine from Ethene-Ethine Complex at %s Req' % (molname.group(2))

for item in mol17:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'MX-2 Benzene-Water Complex at %s Req, CS' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Benzene-Water Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Benzene-Water Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Water from Benzene-Water Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Benzene-Water Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Water from Benzene-Water Complex at %s Req' % (molname.group(2))

for item in mol18:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'MX-3 Benzene-Ammonia Complex at %s Req, CS' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Benzene-Ammonia Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Benzene-Ammonia Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Ammonia from Benzene-Ammonia Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Benzene-Ammonia Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Ammonia from Benzene-Ammonia Complex at %s Req' % (molname.group(2))

for item in mol19:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'MX-4 Benzene-HCN Complex at %s Req, CS' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Benzene-HCN Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Benzene-HCN Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'HCN from Benzene-HCN Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Benzene-HCN Complex at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'HCN from Benzene-HCN Complex at %s Req' % (molname.group(2))

for item in mol20:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'DD-5 Benzene Dimer T-Shape at %s Req, C2V' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Benzene Dimer T-Shape at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Benzene Dimer T-Shape at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Benzene from Benzene Dimer T-Shape at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Benzene Dimer T-Shape at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Benzene from Benzene Dimer T-Shape at %s Req' % (molname.group(2))

for item in mol21:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'MX-6 Indole-Benzene Complex T-Shape at %s Req, C1' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Indole-Benzene Complex T-Shape at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Benzene from Indole-Benzene Complex T-Shape at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Indole from Indole-Benzene Complex T-Shape at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Benzene from Indole-Benzene Complex T-Shape at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Indole from Indole-Benzene Complex T-Shape at %s Req' % (molname.group(2))

for item in mol22:
    molname = rxnpattern.match(item)
    TAGL['%s-%s'            % (dbse, item)] = 'MX-7 Phenol Dimer at %s Req, C1' % (molname.group(2))
    TAGL['%s-%s-dimer'      % (dbse, item)] =      'Phenol Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-CP'   % (dbse, item)] =      'Phenol from Phenol Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-CP'   % (dbse, item)] =      'Phenol from Phenol Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoA-unCP' % (dbse, item)] =      'Phenol from Phenol Dimer at %s Req' % (molname.group(2))
    TAGL['%s-%s-monoB-unCP' % (dbse, item)] =      'Phenol from Phenol Dimer at %s Req' % (molname.group(2))

TAGL['dbse'] = 'interaction energies for dissociation curves of organic bimolecular complexes'
TAGL['default'] = 'entire database'
TAGL['small'] = 'few computationally quick systems'
TAGL['large'] = 'most computationally expensive systems'
TAGL['equilibrium'] = 'minimum-energy systems on dissociation curves'
TAGL['HB'] = 'hydrogen-bonded systems'
TAGL['MX'] = 'mixed-influence systems'
TAGL['DD'] = 'dispersion-dominated systems'
TAGL['mol1'] = 'dissociation curve for ammonia dimer'
TAGL['mol2'] = 'dissociation curve for water dimer'
TAGL['mol3'] = 'dissociation curve for formic acid dimer'
TAGL['mol4'] = 'dissociation curve for formamide dimer'
TAGL['mol5'] = 'dissociation curve for uracil dimer, HB'
TAGL['mol6'] = 'dissociation curve for 2-Pyridone-2-Aminopyridine'
TAGL['mol7'] = 'dissociation curve for adenine-thymine, WC'
TAGL['mol8'] = 'dissociation curve for methane dimer'
TAGL['mol9'] = 'dissociation curve for ethene dimer'
TAGL['mol10'] = 'dissociation curve for benzene-methane'
TAGL['mol11'] = 'dissociation curve for benzene-dimer, parallel-displaced'
TAGL['mol12'] = 'dissociation curve for pyrazine'
TAGL['mol13'] = 'dissociation curve for uracil dimer, stacked'
TAGL['mol14'] = 'dissociation curve for indole-benzene, stacked'
TAGL['mol15'] = 'dissociation curve for adenine-thymine, stacked'
TAGL['mol16'] = 'dissociation curve for ethene-ethyne'
TAGL['mol17'] = 'dissociation curve for benzene-water'
TAGL['mol18'] = 'dissociation curve for benzene-ammonia'
TAGL['mol19'] = 'dissociation curve for benzene-HCN'
TAGL['mol20'] = 'dissociation curve for benzene dimer, t-shaped'
TAGL['mol21'] = 'dissociation curve for indole-benzene, t-shaped'
TAGL['mol22'] = 'dissociation curve for phenol dimer'

# <<< Geometry Specification Strings >>>
GEOS = {}

GEOS['%s-%s-dimer' % (dbse, '1-0.7')] = qcdb.Molecule("""
0 1
N         -0.53502055 -0.86157001  0.00000000
H         -1.14205870 -0.82574073 -0.80956500
H         -1.14205870 -0.82574073  0.80956500
H          0.00000000  0.00000000  0.00000000
--
0 1
N          1.75281655  0.00000000  0.00000000
H          2.35985469 -0.03582927 -0.80956500
H          1.21779599 -0.86157001  0.00000000
H          2.35985469 -0.03582927  0.80956500
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '1-0.8')] = qcdb.Molecule("""
0 1
N         -0.53502055 -0.86157001  0.00000000
H         -1.14205870 -0.82574073 -0.80956500
H         -1.14205870 -0.82574073  0.80956500
H          0.00000000  0.00000000  0.00000000
--
0 1
N          2.00321891  0.00000000  0.00000000
H          2.61025706 -0.03582927 -0.80956500
H          1.46819836 -0.86157001  0.00000000
H          2.61025706 -0.03582927  0.80956500
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '1-0.9')] = qcdb.Molecule("""
0 1
N   -0.535020551  -0.861570006   0.000000000
H   -1.142058700  -0.825740733  -0.809565000
H   -1.142058700  -0.825740733   0.809565000
H    0.000000000   0.000000000   0.000000000
--
0 1
N    2.253621272   0.000000000   0.000000000
H    2.860659421  -0.035829274  -0.809565000
H    1.718600721  -0.861570006   0.000000000
H    2.860659421  -0.035829274   0.809565000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '1-1.0')] = qcdb.Molecule("""
0 1
N   -1.578718000  -0.046611000   0.000000000
H   -2.158621000   0.136396000  -0.809565000
H   -2.158621000   0.136396000   0.809565000
H   -0.849471000   0.658193000   0.000000000
--
0 1
N    1.578718000   0.046611000   0.000000000
H    2.158621000  -0.136396000  -0.809565000
H    0.849471000  -0.658193000   0.000000000
H    2.158621000  -0.136396000   0.809565000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '1-1.2')] = qcdb.Molecule("""
0 1
N   -0.535020551  -0.861570006   0.000000000
H   -1.142058700  -0.825740733  -0.809565000
H   -1.142058700  -0.825740733   0.809565000
H    0.000000000   0.000000000   0.000000000
--
0 1
N    3.004828362   0.000000000   0.000000000
H    3.611866511  -0.035829274  -0.809565000
H    2.469807811  -0.861570006   0.000000000
H    3.611866511  -0.035829274   0.809565000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '1-1.5')] = qcdb.Molecule("""
0 1
N   -0.535020551  -0.861570006   0.000000000
H   -1.142058700  -0.825740733  -0.809565000
H   -1.142058700  -0.825740733   0.809565000
H    0.000000000   0.000000000   0.000000000
--
0 1
N    3.756035452   0.000000000   0.000000000
H    4.363073601  -0.035829274  -0.809565000
H    3.221014901  -0.861570006   0.000000000
H    4.363073601  -0.035829274   0.809565000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '1-2.0')] = qcdb.Molecule("""
0 1
N   -0.535020551  -0.861570006   0.000000000
H   -1.142058700  -0.825740733  -0.809565000
H   -1.142058700  -0.825740733   0.809565000
H    0.000000000   0.000000000   0.000000000
--
0 1
N    5.008047270   0.000000000   0.000000000
H    5.615085419  -0.035829274  -0.809565000
H    4.473026719  -0.861570006   0.000000000
H    5.615085419  -0.035829274   0.809565000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '2-0.7')] = qcdb.Molecule("""
0 1
O         -0.95633265 -0.12063836  0.00000000
H         -1.30753517  0.76970327  0.00000000
H          0.00000000  0.00000000  0.00000000
--
0 1
O          1.36610958  0.00000000  0.00000000
H          1.67807391 -0.49684729 -0.75856100
H          1.67807391 -0.49684729  0.75856100
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '2-0.8')] = qcdb.Molecule("""
0 1
O         -0.95633265 -0.12063836  0.00000000
H         -1.30753517  0.76970327  0.00000000
H          0.00000000  0.00000000  0.00000000
--
0 1
O          1.56126809  0.00000000  0.00000000
H          1.87323242 -0.49684729 -0.75856100
H          1.87323242 -0.49684729  0.75856100
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '2-0.9')] = qcdb.Molecule("""
0 1
O   -0.956332646  -0.120638358   0.000000000
H   -1.307535174   0.769703274   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
O    1.756426600   0.000000000   0.000000000
H    2.068390928  -0.496847294  -0.758561000
H    2.068390928  -0.496847294   0.758561000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '2-1.0')] = qcdb.Molecule("""
0 1
O   -1.551007000  -0.114520000   0.000000000
H   -1.934259000   0.762503000   0.000000000
H   -0.599677000   0.040712000   0.000000000
--
0 1
O    1.350625000   0.111469000   0.000000000
H    1.680398000  -0.373741000  -0.758561000
H    1.680398000  -0.373741000   0.758561000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '2-1.2')] = qcdb.Molecule("""
0 1
O   -0.956332646  -0.120638358   0.000000000
H   -1.307535174   0.769703274   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
O    2.341902133   0.000000000   0.000000000
H    2.653866461  -0.496847294  -0.758561000
H    2.653866461  -0.496847294   0.758561000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '2-1.5')] = qcdb.Molecule("""
0 1
O   -0.956332646  -0.120638358   0.000000000
H   -1.307535174   0.769703274   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
O    2.927377666   0.000000000   0.000000000
H    3.239341994  -0.496847294  -0.758561000
H    3.239341994  -0.496847294   0.758561000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '2-2.0')] = qcdb.Molecule("""
0 1
O   -0.956332646  -0.120638358   0.000000000
H   -1.307535174   0.769703274   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
O    3.903170222   0.000000000   0.000000000
H    4.215134550  -0.496847294  -0.758561000
H    4.215134550  -0.496847294   0.758561000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '3-0.7')] = qcdb.Molecule("""
0 1
C  -1.43494426e+00  -1.23664395e+00   0.00000000
O  -9.95009531e-01   1.87669300e-03   0.00000000
O  -7.52030700e-01  -2.24846554e+00   0.00000000
H  -2.52766058e+00  -1.27695058e+00   0.00000000
H   0.00000000e+00   0.00000000e+00   0.00000000
--
0 1
C   1.85214060e+00  -1.01182159e+00   0.00000000
O   1.41220587e+00  -2.25034224e+00   0.00000000
O   1.16922704e+00   0.00000000e+00   0.00000000
H   2.94485692e+00  -9.71514961e-01   0.00000000
H   4.17196342e-01  -2.24846554e+00   0.00000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '3-0.8')] = qcdb.Molecule("""
0 1
C  -1.43494426e+00  -1.23664395e+00   0.00000000
O  -9.95009531e-01   1.87669300e-03   0.00000000
O  -7.52030700e-01  -2.24846554e+00   0.00000000
H  -2.52766058e+00  -1.27695058e+00   0.00000000
H   0.00000000e+00   0.00000000e+00   0.00000000
--
0 1
C   2.01917304e+00  -1.01182159e+00   0.00000000
O   1.57923831e+00  -2.25034224e+00   0.00000000
O   1.33625948e+00   0.00000000e+00   0.00000000
H   3.11188936e+00  -9.71514961e-01   0.00000000
H   5.84228776e-01  -2.24846554e+00   0.00000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '3-0.9')] = qcdb.Molecule("""
0 1
C   -1.434944263  -1.236643950   0.000000000
O   -0.995009531   0.001876693   0.000000000
O   -0.752030700  -2.248465543   0.000000000
H   -2.527660580  -1.276950582   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
C    2.186205474  -1.011821594   0.000000000
O    1.746270742  -2.250342236   0.000000000
O    1.503291911   0.000000000   0.000000000
H    3.278921791  -0.971514961   0.000000000
H    0.751261211  -2.248465543   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '3-1.0')] = qcdb.Molecule("""
0 1
C   -1.888896000  -0.179692000   0.000000000
O   -1.493280000   1.073689000   0.000000000
O   -1.170435000  -1.166590000   0.000000000
H   -2.979488000  -0.258829000   0.000000000
H   -0.498833000   1.107195000   0.000000000
--
0 1
C    1.888896000   0.179692000   0.000000000
O    1.493280000  -1.073689000   0.000000000
O    1.170435000   1.166590000   0.000000000
H    2.979488000   0.258829000   0.000000000
H    0.498833000  -1.107195000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '3-1.2')] = qcdb.Molecule("""
0 1
C   -1.434944263  -1.236643950   0.000000000
O   -0.995009531   0.001876693   0.000000000
O   -0.752030700  -2.248465543   0.000000000
H   -2.527660580  -1.276950582   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
C    2.687302778  -1.011821594   0.000000000
O    2.247368046  -2.250342236   0.000000000
O    2.004389215   0.000000000   0.000000000
H    3.780019095  -0.971514961   0.000000000
H    1.252358515  -2.248465543   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '3-1.5')] = qcdb.Molecule("""
0 1
C   -1.434944263  -1.236643950   0.000000000
O   -0.995009531   0.001876693   0.000000000
O   -0.752030700  -2.248465543   0.000000000
H   -2.527660580  -1.276950582   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
C    3.188400082  -1.011821594   0.000000000
O    2.748465350  -2.250342236   0.000000000
O    2.505486519   0.000000000   0.000000000
H    4.281116399  -0.971514961   0.000000000
H    1.753455819  -2.248465543   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '3-2.0')] = qcdb.Molecule("""
0 1
C   -1.434944263  -1.236643950   0.000000000
O   -0.995009531   0.001876693   0.000000000
O   -0.752030700  -2.248465543   0.000000000
H   -2.527660580  -1.276950582   0.000000000
H    0.000000000   0.000000000   0.000000000
--
0 1
C    4.023562255  -1.011821594   0.000000000
O    3.583627523  -2.250342236   0.000000000
O    3.340648692   0.000000000   0.000000000
H    5.116278572  -0.971514961   0.000000000
H    2.588617992  -2.248465543   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '4-0.7')] = qcdb.Molecule("""
0 1
C         -0.60412015 -1.07034623  0.00000000
O          0.00000000  0.00000000  0.00000000
N         -0.03527368 -2.28627761  0.00000000
H         -0.62084753 -3.10091587  0.00000000
H          0.98235653 -2.38710371  0.00000000
H         -1.70418544 -1.09860749  0.00000000
--
0 1
C          2.87487022 -1.31675748  0.00000000
O          2.27075007 -2.38710371  0.00000000
N          2.30602374 -0.10082610  0.00000000
H          2.89159759  0.71381216  0.00000000
H          1.28839354  0.00000000  0.00000000
H          3.97493551 -1.28849622  0.00000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '4-0.8')] = qcdb.Molecule("""
0 1
C         -0.60412015 -1.07034623  0.00000000
O          0.00000000  0.00000000  0.00000000
N         -0.03527368 -2.28627761  0.00000000
H         -0.62084753 -3.10091587  0.00000000
H          0.98235653 -2.38710371  0.00000000
H         -1.70418544 -1.09860749  0.00000000
--
0 1
C          3.05892644 -1.31675748  0.00000000
O          2.45480629 -2.38710371  0.00000000
N          2.49007996 -0.10082610  0.00000000
H          3.07565381  0.71381216  0.00000000
H          1.47244976  0.00000000  0.00000000
H          4.15899173 -1.28849622  0.00000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '4-0.9')] = qcdb.Molecule("""
0 1
C   -0.604120150  -1.070346233   0.000000000
O    0.000000000   0.000000000   0.000000000
N   -0.035273679  -2.286277608   0.000000000
H   -0.620847527  -3.100915874   0.000000000
H    0.982356530  -2.387103713   0.000000000
H   -1.704185444  -1.098607493   0.000000000
--
0 1
C    3.242982655  -1.316757480   0.000000000
O    2.638862505  -2.387103713   0.000000000
N    2.674136184  -0.100826104   0.000000000
H    3.259710032   0.713812161   0.000000000
H    1.656505975   0.000000000   0.000000000
H    4.343047949  -1.288496220   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '4-1.0')] = qcdb.Molecule("""
0 1
C   -2.018649000   0.052883000   0.000000000
O   -1.452200000   1.143634000   0.000000000
N   -1.407770000  -1.142484000   0.000000000
H   -1.964596000  -1.977036000   0.000000000
H   -0.387244000  -1.207782000   0.000000000
H   -3.117061000  -0.013701000   0.000000000
--
0 1
C    2.018649000  -0.052883000   0.000000000
O    1.452200000  -1.143634000   0.000000000
N    1.407770000   1.142484000   0.000000000
H    1.964596000   1.977036000   0.000000000
H    0.387244000   1.207782000   0.000000000
H    3.117061000   0.013701000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '4-1.2')] = qcdb.Molecule("""
0 1
C   -0.604120150  -1.070346233   0.000000000
O    0.000000000   0.000000000   0.000000000
N   -0.035273679  -2.286277608   0.000000000
H   -0.620847527  -3.100915874   0.000000000
H    0.982356530  -2.387103713   0.000000000
H   -1.704185444  -1.098607493   0.000000000
--
0 1
C    3.795151314  -1.316757480   0.000000000
O    3.191031164  -2.387103713   0.000000000
N    3.226304843  -0.100826104   0.000000000
H    3.811878691   0.713812161   0.000000000
H    2.208674634   0.000000000   0.000000000
H    4.895216608  -1.288496220   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '4-1.5')] = qcdb.Molecule("""
0 1
C   -0.604120150  -1.070346233   0.000000000
O    0.000000000   0.000000000   0.000000000
N   -0.035273679  -2.286277608   0.000000000
H   -0.620847527  -3.100915874   0.000000000
H    0.982356530  -2.387103713   0.000000000
H   -1.704185444  -1.098607493   0.000000000
--
0 1
C    4.347319973  -1.316757480   0.000000000
O    3.743199823  -2.387103713   0.000000000
N    3.778473502  -0.100826104   0.000000000
H    4.364047350   0.713812161   0.000000000
H    2.760843293   0.000000000   0.000000000
H    5.447385267  -1.288496220   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '4-2.0')] = qcdb.Molecule("""
0 1
C   -0.604120150  -1.070346233   0.000000000
O    0.000000000   0.000000000   0.000000000
N   -0.035273679  -2.286277608   0.000000000
H   -0.620847527  -3.100915874   0.000000000
H    0.982356530  -2.387103713   0.000000000
H   -1.704185444  -1.098607493   0.000000000
--
0 1
C    5.267601070  -1.316757480   0.000000000
O    4.663480920  -2.387103713   0.000000000
N    4.698754599  -0.100826104   0.000000000
H    5.284328447   0.713812161   0.000000000
H    3.681124390   0.000000000   0.000000000
H    6.367666364  -1.288496220   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '5-0.7')] = qcdb.Molecule("""
0 1
O          0.00000000  0.00000000  0.00000000
C         -0.66424394  1.03687915  0.00000000
N         -0.10866344  2.28638952  0.00000000
C         -0.86469194  3.42752195  0.00000000
C         -2.21423160  3.40390953  0.00000000
C         -2.90986986  2.13180389  0.00000000
N         -2.03492462  1.02930119  0.00000000
O         -4.11552152  1.95873396  0.00000000
H         -2.79384033  4.31079935  0.00000000
H          0.91790819  2.33432991  0.00000000
H         -2.46932580  0.11655133  0.00000000
H         -0.30003763  4.34802404  0.00000000
--
0 1
O          2.16009777  2.33432991  0.00000000
C          2.82434171  1.29745076  0.00000000
N          2.26876121  0.04794039  0.00000000
C          3.02478971 -1.09319205  0.00000000
C          4.37432937 -1.06957963  0.00000000
C          5.06996763  0.20252601  0.00000000
N          4.19502240  1.30502871  0.00000000
O          6.27561930  0.37559595  0.00000000
H          4.95393811 -1.97646944  0.00000000
H          1.24218958  0.00000000  0.00000000
H          4.62942358  2.21777858  0.00000000
H          2.46013541 -2.01369414  0.00000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '5-0.8')] = qcdb.Molecule("""
0 1
O          0.00000000  0.00000000  0.00000000
C         -0.66424394  1.03687915  0.00000000
N         -0.10866344  2.28638952  0.00000000
C         -0.86469194  3.42752195  0.00000000
C         -2.21423160  3.40390953  0.00000000
C         -2.90986986  2.13180389  0.00000000
N         -2.03492462  1.02930119  0.00000000
O         -4.11552152  1.95873396  0.00000000
H         -2.79384033  4.31079935  0.00000000
H          0.91790819  2.33432991  0.00000000
H         -2.46932580  0.11655133  0.00000000
H         -0.30003763  4.34802404  0.00000000
--
0 1
O          2.33755343  2.33432991  0.00000000
C          3.00179737  1.29745076  0.00000000
N          2.44621687  0.04794039  0.00000000
C          3.20224537 -1.09319205  0.00000000
C          4.55178503 -1.06957963  0.00000000
C          5.24742329  0.20252601  0.00000000
N          4.37247805  1.30502871  0.00000000
O          6.45307495  0.37559595  0.00000000
H          5.13139376 -1.97646944  0.00000000
H          1.41964524  0.00000000  0.00000000
H          4.80687923  2.21777858  0.00000000
H          2.63759106 -2.01369414  0.00000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '5-0.9')] = qcdb.Molecule("""
0 1
O    0.000000000   0.000000000   0.000000000
C   -0.664243938   1.036879148   0.000000000
N   -0.108663437   2.286389518   0.000000000
C   -0.864691937   3.427521953   0.000000000
C   -2.214231597   3.403909532   0.000000000
C   -2.909869859   2.131803891   0.000000000
N   -2.034924624   1.029301194   0.000000000
O   -4.115521524   1.958733959   0.000000000
H   -2.793840332   4.310799346   0.000000000
H    0.917908194   2.334329905   0.000000000
H   -2.469325804   0.116551326   0.000000000
H   -0.300037631   4.348024043   0.000000000
--
0 1
O    2.515009084   2.334329905   0.000000000
C    3.179253022   1.297450757   0.000000000
N    2.623672521   0.047940387   0.000000000
C    3.379701020  -1.093192048   0.000000000
C    4.729240680  -1.069579627   0.000000000
C    5.424878943   0.202526014   0.000000000
N    4.549933708   1.305028711   0.000000000
O    6.630530608   0.375595946   0.000000000
H    5.308849416  -1.976469441   0.000000000
H    1.597100890   0.000000000   0.000000000
H    4.984334888   2.217778579   0.000000000
H    2.815046715  -2.013694138   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '5-1.0')] = qcdb.Molecule("""
0 1
O   -1.466332000   1.012169000   0.000000000
C   -0.628146000   1.914268000   0.000000000
N    0.720509000   1.688269000   0.000000000
C    1.636729000   2.705276000   0.000000000
C    1.276904000   4.006176000   0.000000000
C   -0.128601000   4.362155000   0.000000000
N   -0.977723000   3.239643000   0.000000000
O   -0.597223000   5.486407000   0.000000000
H    2.010350000   4.793864000   0.000000000
H    1.023251000   0.706182000   0.000000000
H   -1.970027000   3.432385000   0.000000000
H    2.669062000   2.388342000   0.000000000
--
0 1
O    1.466332000  -1.012169000   0.000000000
C    0.628146000  -1.914268000   0.000000000
N   -0.720509000  -1.688269000   0.000000000
C   -1.636729000  -2.705276000   0.000000000
C   -1.276904000  -4.006176000   0.000000000
C    0.128601000  -4.362155000   0.000000000
N    0.977723000  -3.239643000   0.000000000
O    0.597223000  -5.486407000   0.000000000
H   -2.010350000  -4.793864000   0.000000000
H   -1.023251000  -0.706182000   0.000000000
H    1.970027000  -3.432385000   0.000000000
H   -2.669062000  -2.388342000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '5-1.2')] = qcdb.Molecule("""
0 1
O    0.000000000   0.000000000   0.000000000
C   -0.664243938   1.036879148   0.000000000
N   -0.108663437   2.286389518   0.000000000
C   -0.864691937   3.427521953   0.000000000
C   -2.214231597   3.403909532   0.000000000
C   -2.909869859   2.131803891   0.000000000
N   -2.034924624   1.029301194   0.000000000
O   -4.115521524   1.958733959   0.000000000
H   -2.793840332   4.310799346   0.000000000
H    0.917908194   2.334329905   0.000000000
H   -2.469325804   0.116551326   0.000000000
H   -0.300037631   4.348024043   0.000000000
--
0 1
O    3.047376048   2.334329905   0.000000000
C    3.711619986   1.297450757   0.000000000
N    3.156039485   0.047940387   0.000000000
C    3.912067984  -1.093192048   0.000000000
C    5.261607644  -1.069579627   0.000000000
C    5.957245907   0.202526014   0.000000000
N    5.082300672   1.305028711   0.000000000
O    7.162897572   0.375595946   0.000000000
H    5.841216380  -1.976469441   0.000000000
H    2.129467854   0.000000000   0.000000000
H    5.516701852   2.217778579   0.000000000
H    3.347413679  -2.013694138   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '5-1.5')] = qcdb.Molecule("""
0 1
O    0.000000000   0.000000000   0.000000000
C   -0.664243938   1.036879148   0.000000000
N   -0.108663437   2.286389518   0.000000000
C   -0.864691937   3.427521953   0.000000000
C   -2.214231597   3.403909532   0.000000000
C   -2.909869859   2.131803891   0.000000000
N   -2.034924624   1.029301194   0.000000000
O   -4.115521524   1.958733959   0.000000000
H   -2.793840332   4.310799346   0.000000000
H    0.917908194   2.334329905   0.000000000
H   -2.469325804   0.116551326   0.000000000
H   -0.300037631   4.348024043   0.000000000
--
0 1
O    3.579743012   2.334329905   0.000000000
C    4.243986950   1.297450757   0.000000000
N    3.688406449   0.047940387   0.000000000
C    4.444434948  -1.093192048   0.000000000
C    5.793974608  -1.069579627   0.000000000
C    6.489612871   0.202526014   0.000000000
N    5.614667636   1.305028711   0.000000000
O    7.695264536   0.375595946   0.000000000
H    6.373583344  -1.976469441   0.000000000
H    2.661834818   0.000000000   0.000000000
H    6.049068816   2.217778579   0.000000000
H    3.879780643  -2.013694138   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '5-2.0')] = qcdb.Molecule("""
0 1
O    0.000000000   0.000000000   0.000000000
C   -0.664243938   1.036879148   0.000000000
N   -0.108663437   2.286389518   0.000000000
C   -0.864691937   3.427521953   0.000000000
C   -2.214231597   3.403909532   0.000000000
C   -2.909869859   2.131803891   0.000000000
N   -2.034924624   1.029301194   0.000000000
O   -4.115521524   1.958733959   0.000000000
H   -2.793840332   4.310799346   0.000000000
H    0.917908194   2.334329905   0.000000000
H   -2.469325804   0.116551326   0.000000000
H   -0.300037631   4.348024043   0.000000000
--
0 1
O    4.467021284   2.334329905   0.000000000
C    5.131265222   1.297450757   0.000000000
N    4.575684721   0.047940387   0.000000000
C    5.331713220  -1.093192048   0.000000000
C    6.681252880  -1.069579627   0.000000000
C    7.376891143   0.202526014   0.000000000
N    6.501945908   1.305028711   0.000000000
O    8.582542808   0.375595946   0.000000000
H    7.260861616  -1.976469441   0.000000000
H    3.549113090   0.000000000   0.000000000
H    6.936347088   2.217778579   0.000000000
H    4.767058915  -2.013694138   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '6-0.7')] = qcdb.Molecule("""
0 1
O  -9.69652624e-01  -2.24561116e+00  -3.86822525e-01
N  -1.03778979e+00   4.50875300e-03  -1.13112700e-03
C  -3.75926130e+00   1.40280680e-02  -1.83757600e-02
C  -3.05772706e+00   1.22163116e+00   2.04402100e-01
C  -1.69239288e+00   1.17200070e+00   2.05277859e-01
C  -1.65006801e+00  -1.22251475e+00  -2.17981663e-01
C  -3.08826439e+00  -1.16182823e+00  -2.21825966e-01
H  -4.84130076e+00   1.67084980e-02  -2.68920470e-02
H  -3.56722182e+00   2.15683108e+00   3.69386687e-01
H  -1.06806457e+00   2.03877945e+00   3.67771502e-01
H  -3.61208850e+00  -2.09070100e+00  -3.90563867e-01
H   0.00000000e+00   0.00000000e+00   0.00000000e+00
--
0 1
N   1.30160597e+00   0.00000000e+00   0.00000000e+00
C   1.98020601e+00  -1.14532421e+00   1.92591910e-01
C   3.38857185e+00  -1.16867747e+00   1.96637005e-01
C   4.08768558e+00   5.47708300e-03  -1.72323900e-03
C   3.38329557e+00   1.19444766e+00  -2.02961469e-01
C   2.00100662e+00   1.13032803e+00  -1.92845808e-01
H   3.90738672e+00  -2.10397523e+00   3.56345736e-01
H   5.16911435e+00  -3.10336700e-03  -1.91123500e-03
H   3.88787775e+00   2.13463205e+00  -3.64687797e-01
H   1.41022754e+00   2.02525842e+00  -3.49790900e-01
N   1.24832878e+00  -2.27220155e+00   4.35153550e-01
H   1.72973150e+00  -3.14588817e+00   3.15408858e-01
H   2.72633521e-01  -2.27044207e+00   1.33172072e-01
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '6-0.8')] = qcdb.Molecule("""
0 1
O  -9.69652624e-01  -2.24561116e+00  -3.86822525e-01
N  -1.03778979e+00   4.50875300e-03  -1.13112700e-03
C  -3.75926130e+00   1.40280680e-02  -1.83757600e-02
C  -3.05772706e+00   1.22163116e+00   2.04402100e-01
C  -1.69239288e+00   1.17200070e+00   2.05277859e-01
C  -1.65006801e+00  -1.22251475e+00  -2.17981663e-01
C  -3.08826439e+00  -1.16182823e+00  -2.21825966e-01
H  -4.84130076e+00   1.67084980e-02  -2.68920470e-02
H  -3.56722182e+00   2.15683108e+00   3.69386687e-01
H  -1.06806457e+00   2.03877945e+00   3.67771502e-01
H  -3.61208850e+00  -2.09070100e+00  -3.90563867e-01
H   0.00000000e+00   0.00000000e+00   0.00000000e+00
--
0 1
N   1.48754968e+00   0.00000000e+00   0.00000000e+00
C   2.16614972e+00  -1.14532421e+00   1.92591910e-01
C   3.57451556e+00  -1.16867747e+00   1.96637005e-01
C   4.27362929e+00   5.47708300e-03  -1.72323900e-03
C   3.56923928e+00   1.19444766e+00  -2.02961469e-01
C   2.18695033e+00   1.13032803e+00  -1.92845808e-01
H   4.09333042e+00  -2.10397523e+00   3.56345736e-01
H   5.35505806e+00  -3.10336700e-03  -1.91123500e-03
H   4.07382146e+00   2.13463205e+00  -3.64687797e-01
H   1.59617125e+00   2.02525842e+00  -3.49790900e-01
N   1.43427249e+00  -2.27220155e+00   4.35153550e-01
H   1.91567521e+00  -3.14588817e+00   3.15408858e-01
H   4.58577231e-01  -2.27044207e+00   1.33172072e-01
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '6-0.9')] = qcdb.Molecule("""
0 1
O   -0.969652624  -2.245611164  -0.386822525
N   -1.037789793   0.004508753  -0.001131127
C   -3.759261297   0.014028068  -0.018375760
C   -3.057727058   1.221631156   0.204402100
C   -1.692392879   1.172000703   0.205277859
C   -1.650068007  -1.222514751  -0.217981663
C   -3.088264390  -1.161828225  -0.221825966
H   -4.841300764   0.016708498  -0.026892047
H   -3.567221821   2.156831083   0.369386687
H   -1.068064568   2.038779450   0.367771502
H   -3.612088503  -2.090701001  -0.390563867
H    0.000000000   0.000000000   0.000000000
--
0 1
N    1.673493386   0.000000000   0.000000000
C    2.352093429  -1.145324213   0.192591910
C    3.760459273  -1.168677470   0.196637005
C    4.459573002   0.005477083  -0.001723239
C    3.755182987   1.194447664  -0.202961469
C    2.372894041   1.130328028  -0.192845808
H    4.279274134  -2.103975233   0.356345736
H    5.541001766  -0.003103367  -0.001911235
H    4.259765167   2.134632052  -0.364687797
H    1.782114958   2.025258423  -0.349790900
N    1.620216197  -2.272201547   0.435153550
H    2.101618920  -3.145888174   0.315408858
H    0.644520940  -2.270442069   0.133172072
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '6-1.0')] = qcdb.Molecule("""
0 1
O   -1.397621000  -1.885837000  -0.367306000
N   -1.464255000   0.364183000   0.019230000
C   -4.185740000   0.369667000   0.036096000
C   -3.483260000   1.578311000   0.250075000
C   -2.117950000   1.530705000   0.233838000
C   -2.077383000  -0.863749000  -0.189941000
C   -3.515603000  -0.805195000  -0.175759000
H   -5.267804000   0.370743000   0.041142000
H   -3.992033000   2.512756000   0.421441000
H   -1.492920000   2.398410000   0.388502000
H   -4.040123000  -1.734845000  -0.337927000
H   -0.426527000   0.361213000   0.007354000
--
0 1
N    1.432762000   0.363970000  -0.015951000
C    2.115420000  -0.780345000   0.168110000
C    3.523759000  -0.801610000   0.154503000
C    4.218590000   0.373578000  -0.052593000
C    3.509971000   1.561501000  -0.244976000
C    2.128014000   1.495332000  -0.217537000
H    4.045921000  -1.736136000   0.307688000
H    5.299943000   0.366601000  -0.066335000
H    4.011092000   2.502431000  -0.413005000
H    1.533988000   2.389384000  -0.367057000
N    1.388312000  -1.908304000   0.419815000
H    1.869471000  -2.781277000   0.294038000
H    0.408907000  -1.907994000   0.130086000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '6-1.2')] = qcdb.Molecule("""
0 1
O   -0.969652624  -2.245611164  -0.386822525
N   -1.037789793   0.004508753  -0.001131127
C   -3.759261297   0.014028068  -0.018375760
C   -3.057727058   1.221631156   0.204402100
C   -1.692392879   1.172000703   0.205277859
C   -1.650068007  -1.222514751  -0.217981663
C   -3.088264390  -1.161828225  -0.221825966
H   -4.841300764   0.016708498  -0.026892047
H   -3.567221821   2.156831083   0.369386687
H   -1.068064568   2.038779450   0.367771502
H   -3.612088503  -2.090701001  -0.390563867
H    0.000000000   0.000000000   0.000000000
--
0 1
N    2.231324514   0.000000000   0.000000000
C    2.909924557  -1.145324213   0.192591910
C    4.318290401  -1.168677470   0.196637005
C    5.017404130   0.005477083  -0.001723239
C    4.313014115   1.194447664  -0.202961469
C    2.930725169   1.130328028  -0.192845808
H    4.837105262  -2.103975233   0.356345736
H    6.098832894  -0.003103367  -0.001911235
H    4.817596295   2.134632052  -0.364687797
H    2.339946086   2.025258423  -0.349790900
N    2.178047325  -2.272201547   0.435153550
H    2.659450048  -3.145888174   0.315408858
H    1.202352068  -2.270442069   0.133172072
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '6-1.5')] = qcdb.Molecule("""
0 1
O   -0.969652624  -2.245611164  -0.386822525
N   -1.037789793   0.004508753  -0.001131127
C   -3.759261297   0.014028068  -0.018375760
C   -3.057727058   1.221631156   0.204402100
C   -1.692392879   1.172000703   0.205277859
C   -1.650068007  -1.222514751  -0.217981663
C   -3.088264390  -1.161828225  -0.221825966
H   -4.841300764   0.016708498  -0.026892047
H   -3.567221821   2.156831083   0.369386687
H   -1.068064568   2.038779450   0.367771502
H   -3.612088503  -2.090701001  -0.390563867
H    0.000000000   0.000000000   0.000000000
--
0 1
N    2.789155642   0.000000000   0.000000000
C    3.467755685  -1.145324213   0.192591910
C    4.876121529  -1.168677470   0.196637005
C    5.575235258   0.005477083  -0.001723239
C    4.870845243   1.194447664  -0.202961469
C    3.488556297   1.130328028  -0.192845808
H    5.394936390  -2.103975233   0.356345736
H    6.656664022  -0.003103367  -0.001911235
H    5.375427423   2.134632052  -0.364687797
H    2.897777214   2.025258423  -0.349790900
N    2.735878453  -2.272201547   0.435153550
H    3.217281176  -3.145888174   0.315408858
H    1.760183196  -2.270442069   0.133172072
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '6-2.0')] = qcdb.Molecule("""
0 1
O   -0.969652624  -2.245611164  -0.386822525
N   -1.037789793   0.004508753  -0.001131127
C   -3.759261297   0.014028068  -0.018375760
C   -3.057727058   1.221631156   0.204402100
C   -1.692392879   1.172000703   0.205277859
C   -1.650068007  -1.222514751  -0.217981663
C   -3.088264390  -1.161828225  -0.221825966
H   -4.841300764   0.016708498  -0.026892047
H   -3.567221821   2.156831083   0.369386687
H   -1.068064568   2.038779450   0.367771502
H   -3.612088503  -2.090701001  -0.390563867
H    0.000000000   0.000000000   0.000000000
--
0 1
N    3.718874190   0.000000000   0.000000000
C    4.397474233  -1.145324213   0.192591910
C    5.805840077  -1.168677470   0.196637005
C    6.504953806   0.005477083  -0.001723239
C    5.800563791   1.194447664  -0.202961469
C    4.418274845   1.130328028  -0.192845808
H    6.324654938  -2.103975233   0.356345736
H    7.586382570  -0.003103367  -0.001911235
H    6.305145971   2.134632052  -0.364687797
H    3.827495762   2.025258423  -0.349790900
N    3.665597001  -2.272201547   0.435153550
H    4.146999724  -3.145888174   0.315408858
H    2.689901744  -2.270442069   0.133172072
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '7-0.7')] = qcdb.Molecule("""
0 1
N          0.00000000  0.00000000  0.00000000
C         -0.73868506 -0.15788977  1.11035541
C         -2.13945288 -0.16805356  0.96471256
C         -2.62949719 -0.00866579 -0.33120135
N         -1.91830983  0.15263475 -1.45484404
C         -0.61426222  0.14365987 -1.19354712
N         -3.15298100 -0.31069720  1.88351867
C         -4.24746601 -0.23720033  1.14487498
N         -3.99425073 -0.05660450 -0.18703010
N         -0.13617941 -0.28943385  2.30042802
H          0.05516135  0.26595902 -2.03565509
H         -5.25258545 -0.30895833  1.52540657
H         -4.66840486  0.02624532 -0.92965682
H          0.87687643 -0.32910573  2.35981141
H         -0.70858132 -0.45240707  3.10824060
--
0 1
N          4.31030617  0.15562755 -1.12807516
C          5.00317680 -0.03157353  0.03965251
C          4.38156100 -0.21318055  1.22599931
C          2.92591998 -0.20545954  1.23795900
N          2.31505277 -0.00891377  0.01310903
C          2.92866234  0.17623919 -1.20541710
C          5.10083273 -0.41995094  2.51700092
O          2.25753790 -0.36203166  2.26165430
O          2.33043291  0.34250657 -2.25336777
H          4.79061194  0.28845835 -2.00230090
H          1.27319653  0.00000000  0.00000000
H          6.08042149 -0.02477987 -0.04965000
H          4.83125252  0.35484120  3.23301874
H          4.82014459 -1.37309824  2.96239753
H          6.17860422 -0.40361701  2.36838509
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '7-0.8')] = qcdb.Molecule("""
0 1
N          0.00000000  0.00000000  0.00000000
C         -0.73868506 -0.15788977  1.11035541
C         -2.13945288 -0.16805356  0.96471256
C         -2.62949719 -0.00866579 -0.33120135
N         -1.91830983  0.15263475 -1.45484404
C         -0.61426222  0.14365987 -1.19354712
N         -3.15298100 -0.31069720  1.88351867
C         -4.24746601 -0.23720033  1.14487498
N         -3.99425073 -0.05660450 -0.18703010
N         -0.13617941 -0.28943385  2.30042802
H          0.05516135  0.26595902 -2.03565509
H         -5.25258545 -0.30895833  1.52540657
H         -4.66840486  0.02624532 -0.92965682
H          0.87687643 -0.32910573  2.35981141
H         -0.70858132 -0.45240707  3.10824060
--
0 1
N          4.49219139  0.15562755 -1.12807516
C          5.18506202 -0.03157353  0.03965251
C          4.56344622 -0.21318055  1.22599931
C          3.10780520 -0.20545954  1.23795900
N          2.49693799 -0.00891377  0.01310903
C          3.11054756  0.17623919 -1.20541710
C          5.28271795 -0.41995094  2.51700092
O          2.43942312 -0.36203166  2.26165430
O          2.51231813  0.34250657 -2.25336777
H          4.97249716  0.28845835 -2.00230090
H          1.45508175  0.00000000  0.00000000
H          6.26230671 -0.02477987 -0.04965000
H          5.01313774  0.35484120  3.23301874
H          5.00202981 -1.37309824  2.96239753
H          6.36048944 -0.40361701  2.36838509
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '7-0.9')] = qcdb.Molecule("""
0 1
N    0.000000000   0.000000000   0.000000000
C   -0.738685058  -0.157889771   1.110355410
C   -2.139452884  -0.168053559   0.964712563
C   -2.629497187  -0.008665792  -0.331201352
N   -1.918309833   0.152634753  -1.454844039
C   -0.614262216   0.143659867  -1.193547121
N   -3.152980999  -0.310697201   1.883518666
C   -4.247466012  -0.237200328   1.144874976
N   -3.994250734  -0.056604504  -0.187030096
N   -0.136179412  -0.289433845   2.300428025
H    0.055161346   0.265959015  -2.035655088
H   -5.252585445  -0.308958331   1.525406574
H   -4.668404863   0.026245320  -0.929656824
H    0.876876426  -0.329105732   2.359811410
H   -0.708581316  -0.452407073   3.108240602
--
0 1
N    4.674076612   0.155627547  -1.128075158
C    5.366947235  -0.031573530   0.039652507
C    4.745331442  -0.213180550   1.225999310
C    3.289690418  -0.205459536   1.237959001
N    2.678823212  -0.008913767   0.013109028
C    3.292432779   0.176239188  -1.205417098
C    5.464603172  -0.419950938   2.517000917
O    2.621308338  -0.362031655   2.261654302
O    2.694203350   0.342506569  -2.253367774
H    5.154382378   0.288458351  -2.002300903
H    1.636966971   0.000000000   0.000000000
H    6.444191927  -0.024779868  -0.049650000
H    5.195022957   0.354841198   3.233018736
H    5.183915029  -1.373098243   2.962397530
H    6.542374655  -0.403617008   2.368385087
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '7-1.0')] = qcdb.Molecule("""
0 1
N    0.935015000  -0.027980000  -0.378892000
C    1.673964000  -0.035777000   0.742432000
C    3.074796000  -0.009448000   0.599456000
C    3.564611000   0.019545000  -0.705987000
N    2.853151000   0.025803000  -1.840960000
C    1.549076000   0.001257000  -1.580801000
N    4.088582000  -0.005443000   1.528979000
C    5.182992000   0.025397000   0.787218000
N    4.929487000   0.041240000  -0.556727000
N    1.071618000  -0.076537000   1.939139000
H    0.879444000   0.005026000  -2.431571000
H    6.188259000   0.037554000   1.173882000
H    5.603537000   0.064876000  -1.303681000
H    0.058692000  -0.042376000   2.003918000
H    1.644380000  -0.034739000   2.761916000
--
0 1
N   -3.921173000  -0.000965000  -1.516366000
C   -4.613683000   0.016905000  -0.333652000
C   -3.991739000   0.021935000   0.866334000
C   -2.536137000   0.007465000   0.876672000
N   -1.925648000  -0.011059000  -0.363895000
C   -2.539590000  -0.014947000  -1.596236000
C   -4.710613000   0.041337000   2.173864000
O   -1.867473000   0.011209000   1.912083000
O   -1.941678000  -0.029188000  -2.657378000
H   -4.401717000  -0.003608000  -2.400492000
H   -0.883826000  -0.021617000  -0.378427000
H   -5.690922000   0.026935000  -0.422718000
H   -4.443928000  -0.830257000   2.769566000
H   -4.426706000   0.918618000   2.753026000
H   -5.788397000   0.050553000   2.024728000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '7-1.2')] = qcdb.Molecule("""
0 1
N    0.000000000   0.000000000   0.000000000
C   -0.738685058  -0.157889771   1.110355410
C   -2.139452884  -0.168053559   0.964712563
C   -2.629497187  -0.008665792  -0.331201352
N   -1.918309833   0.152634753  -1.454844039
C   -0.614262216   0.143659867  -1.193547121
N   -3.152980999  -0.310697201   1.883518666
C   -4.247466012  -0.237200328   1.144874976
N   -3.994250734  -0.056604504  -0.187030096
N   -0.136179412  -0.289433845   2.300428025
H    0.055161346   0.265959015  -2.035655088
H   -5.252585445  -0.308958331   1.525406574
H   -4.668404863   0.026245320  -0.929656824
H    0.876876426  -0.329105732   2.359811410
H   -0.708581316  -0.452407073   3.108240602
--
0 1
N    5.219732269   0.155627547  -1.128075158
C    5.912602892  -0.031573530   0.039652507
C    5.290987099  -0.213180550   1.225999310
C    3.835346075  -0.205459536   1.237959001
N    3.224478869  -0.008913767   0.013109028
C    3.838088436   0.176239188  -1.205417098
C    6.010258829  -0.419950938   2.517000917
O    3.166963995  -0.362031655   2.261654302
O    3.239859007   0.342506569  -2.253367774
H    5.700038035   0.288458351  -2.002300903
H    2.182622628   0.000000000   0.000000000
H    6.989847584  -0.024779868  -0.049650000
H    5.740678614   0.354841198   3.233018736
H    5.729570686  -1.373098243   2.962397530
H    7.088030312  -0.403617008   2.368385087
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '7-1.5')] = qcdb.Molecule("""
0 1
N    0.000000000   0.000000000   0.000000000
C   -0.738685058  -0.157889771   1.110355410
C   -2.139452884  -0.168053559   0.964712563
C   -2.629497187  -0.008665792  -0.331201352
N   -1.918309833   0.152634753  -1.454844039
C   -0.614262216   0.143659867  -1.193547121
N   -3.152980999  -0.310697201   1.883518666
C   -4.247466012  -0.237200328   1.144874976
N   -3.994250734  -0.056604504  -0.187030096
N   -0.136179412  -0.289433845   2.300428025
H    0.055161346   0.265959015  -2.035655088
H   -5.252585445  -0.308958331   1.525406574
H   -4.668404863   0.026245320  -0.929656824
H    0.876876426  -0.329105732   2.359811410
H   -0.708581316  -0.452407073   3.108240602
--
0 1
N    5.765387926   0.155627547  -1.128075158
C    6.458258549  -0.031573530   0.039652507
C    5.836642756  -0.213180550   1.225999310
C    4.381001732  -0.205459536   1.237959001
N    3.770134526  -0.008913767   0.013109028
C    4.383744093   0.176239188  -1.205417098
C    6.555914486  -0.419950938   2.517000917
O    3.712619652  -0.362031655   2.261654302
O    3.785514664   0.342506569  -2.253367774
H    6.245693692   0.288458351  -2.002300903
H    2.728278285   0.000000000   0.000000000
H    7.535503241  -0.024779868  -0.049650000
H    6.286334271   0.354841198   3.233018736
H    6.275226343  -1.373098243   2.962397530
H    7.633685969  -0.403617008   2.368385087
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '7-2.0')] = qcdb.Molecule("""
0 1
N    0.000000000   0.000000000   0.000000000
C   -0.738685058  -0.157889771   1.110355410
C   -2.139452884  -0.168053559   0.964712563
C   -2.629497187  -0.008665792  -0.331201352
N   -1.918309833   0.152634753  -1.454844039
C   -0.614262216   0.143659867  -1.193547121
N   -3.152980999  -0.310697201   1.883518666
C   -4.247466012  -0.237200328   1.144874976
N   -3.994250734  -0.056604504  -0.187030096
N   -0.136179412  -0.289433845   2.300428025
H    0.055161346   0.265959015  -2.035655088
H   -5.252585445  -0.308958331   1.525406574
H   -4.668404863   0.026245320  -0.929656824
H    0.876876426  -0.329105732   2.359811410
H   -0.708581316  -0.452407073   3.108240602
--
0 1
N    6.674814021   0.155627547  -1.128075158
C    7.367684644  -0.031573530   0.039652507
C    6.746068851  -0.213180550   1.225999310
C    5.290427827  -0.205459536   1.237959001
N    4.679560621  -0.008913767   0.013109028
C    5.293170188   0.176239188  -1.205417098
C    7.465340581  -0.419950938   2.517000917
O    4.622045747  -0.362031655   2.261654302
O    4.694940759   0.342506569  -2.253367774
H    7.155119787   0.288458351  -2.002300903
H    3.637704380   0.000000000   0.000000000
H    8.444929336  -0.024779868  -0.049650000
H    7.195760366   0.354841198   3.233018736
H    7.184652438  -1.373098243   2.962397530
H    8.543112064  -0.403617008   2.368385087
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '8-0.7')] = qcdb.Molecule("""
0 1
C   0.00000000e+00   0.00000000e+00   0.00000000e+00
H   3.64514644e-01   5.13239461e-01  -8.88512354e-01
H   3.64514644e-01   5.13105641e-01   8.88589641e-01
H   3.64215723e-01  -1.02622643e+00  -7.72780000e-05
H  -1.08912298e+00   3.11014000e-04   2.30000000e-08
--
0 1
C   2.60282541e+00   0.00000000e+00   0.00000000e+00
H   3.69194839e+00  -3.11014000e-04  -2.30000000e-08
H   2.23831076e+00  -5.13105641e-01  -8.88589641e-01
H   2.23831076e+00  -5.13239461e-01   8.88512354e-01
H   2.23860968e+00   1.02622643e+00   7.72780000e-05
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '8-0.8')] = qcdb.Molecule("""
0 1
C   0.00000000e+00   0.00000000e+00   0.00000000e+00
H   3.64514644e-01   5.13239461e-01  -8.88512354e-01
H   3.64514644e-01   5.13105641e-01   8.88589641e-01
H   3.64215723e-01  -1.02622643e+00  -7.72780000e-05
H  -1.08912298e+00   3.11014000e-04   2.30000000e-08
--
0 1
C   2.97465761e+00   0.00000000e+00   0.00000000e+00
H   4.06378059e+00  -3.11014000e-04  -2.30000000e-08
H   2.61014296e+00  -5.13105641e-01  -8.88589641e-01
H   2.61014296e+00  -5.13239461e-01   8.88512354e-01
H   2.61044188e+00   1.02622643e+00   7.72780000e-05
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '8-0.9')] = qcdb.Molecule("""
0 1
C    0.000000000   0.000000000   0.000000000
H    0.364514644   0.513239461  -0.888512354
H    0.364514644   0.513105641   0.888589641
H    0.364215723  -1.026226426  -0.000077278
H   -1.089122980   0.000311014   0.000000023
--
0 1
C    3.346489810   0.000000000   0.000000000
H    4.435612789  -0.000311014  -0.000000023
H    2.981975165  -0.513105641  -0.888589641
H    2.981975165  -0.513239461   0.888512354
H    2.982274086   1.026226426   0.000077278
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '8-1.0')] = qcdb.Molecule("""
0 1
C    0.000000000  -0.000140000   1.859161000
H   -0.888551000   0.513060000   1.494685000
H    0.888551000   0.513060000   1.494685000
H    0.000000000  -1.026339000   1.494868000
H    0.000000000   0.000089000   2.948284000
--
0 1
C    0.000000000   0.000140000  -1.859161000
H    0.000000000  -0.000089000  -2.948284000
H   -0.888551000  -0.513060000  -1.494685000
H    0.888551000  -0.513060000  -1.494685000
H    0.000000000   1.026339000  -1.494868000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '8-1.2')] = qcdb.Molecule("""
0 1
C    0.000000000   0.000000000   0.000000000
H    0.364514644   0.513239461  -0.888512354
H    0.364514644   0.513105641   0.888589641
H    0.364215723  -1.026226426  -0.000077278
H   -1.089122980   0.000311014   0.000000023
--
0 1
C    4.461986413   0.000000000   0.000000000
H    5.551109392  -0.000311014  -0.000000023
H    4.097471768  -0.513105641  -0.888589641
H    4.097471768  -0.513239461   0.888512354
H    4.097770689   1.026226426   0.000077278
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '8-1.5')] = qcdb.Molecule("""
0 1
C    0.000000000   0.000000000   0.000000000
H    0.364514644   0.513239461  -0.888512354
H    0.364514644   0.513105641   0.888589641
H    0.364215723  -1.026226426  -0.000077278
H   -1.089122980   0.000311014   0.000000023
--
0 1
C    5.577483016   0.000000000   0.000000000
H    6.666605995  -0.000311014  -0.000000023
H    5.212968371  -0.513105641  -0.888589641
H    5.212968371  -0.513239461   0.888512354
H    5.213267292   1.026226426   0.000077278
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '8-2.0')] = qcdb.Molecule("""
0 1
C    0.000000000   0.000000000   0.000000000
H    0.364514644   0.513239461  -0.888512354
H    0.364514644   0.513105641   0.888589641
H    0.364215723  -1.026226426  -0.000077278
H   -1.089122980   0.000311014   0.000000023
--
0 1
C    7.436644022   0.000000000   0.000000000
H    8.525767001  -0.000311014  -0.000000023
H    7.072129377  -0.513105641  -0.888589641
H    7.072129377  -0.513239461   0.888512354
H    7.072428298   1.026226426   0.000077278
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '9-0.7')] = qcdb.Molecule("""
0 1
C         0.000000  -0.471925   0.471925
C         0.000000   0.471925  -0.471925
H         0.922986  -0.872422   0.872422
H         0.922986   0.872422  -0.872422
H        -0.924197  -0.870464   0.870464
H        -0.924197   0.870464  -0.870464
--
0 1
C         2.6027554  0.471925   0.471925
C         2.6027554 -0.471925  -0.471925
H         1.6797694  0.872422   0.872422
H         1.6797694 -0.872422  -0.872422
H         3.5269524  0.870464   0.870464
H         3.5269524 -0.870464  -0.870464
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '9-0.8')] = qcdb.Molecule("""
0 1
C         0.000000  -0.471925   0.471925
C         0.000000   0.471925  -0.471925
H         0.922986  -0.872422   0.872422
H         0.922986   0.872422  -0.872422
H        -0.924197  -0.870464   0.870464
H        -0.924197   0.870464  -0.870464
--
0 1
C         2.9745776  0.471925   0.471925
C         2.9745776 -0.471925  -0.471925
H         2.0515916  0.872422   0.872422
H         2.0515916 -0.872422  -0.872422
H         3.8987746  0.870464   0.870464
H         3.8987746 -0.870464  -0.870464
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '9-0.9')] = qcdb.Molecule("""
0 1
C    0.000000000  -0.471925000   0.471925000
C    0.000000000   0.471925000  -0.471925000
H    0.922986000  -0.872422000   0.872422000
H    0.922986000   0.872422000  -0.872422000
H   -0.924197000  -0.870464000   0.870464000
H   -0.924197000   0.870464000  -0.870464000
--
0 1
C    3.346399800   0.471925000   0.471925000
C    3.346399800  -0.471925000  -0.471925000
H    2.423413800   0.872422000   0.872422000
H    2.423413800  -0.872422000  -0.872422000
H    4.270596800   0.870464000   0.870464000
H    4.270596800  -0.870464000  -0.870464000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '9-1.0')] = qcdb.Molecule("""
0 1
C   -0.471925000  -0.471925000  -1.859111000
C    0.471925000   0.471925000  -1.859111000
H   -0.872422000  -0.872422000  -0.936125000
H    0.872422000   0.872422000  -0.936125000
H   -0.870464000  -0.870464000  -2.783308000
H    0.870464000   0.870464000  -2.783308000
--
0 1
C   -0.471925000   0.471925000   1.859111000
C    0.471925000  -0.471925000   1.859111000
H   -0.872422000   0.872422000   0.936125000
H    0.872422000  -0.872422000   0.936125000
H   -0.870464000   0.870464000   2.783308000
H    0.870464000  -0.870464000   2.783308000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '9-1.2')] = qcdb.Molecule("""
0 1
C    0.000000000  -0.471925000   0.471925000
C    0.000000000   0.471925000  -0.471925000
H    0.922986000  -0.872422000   0.872422000
H    0.922986000   0.872422000  -0.872422000
H   -0.924197000  -0.870464000   0.870464000
H   -0.924197000   0.870464000  -0.870464000
--
0 1
C    4.461866400   0.471925000   0.471925000
C    4.461866400  -0.471925000  -0.471925000
H    3.538880400   0.872422000   0.872422000
H    3.538880400  -0.872422000  -0.872422000
H    5.386063400   0.870464000   0.870464000
H    5.386063400  -0.870464000  -0.870464000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '9-1.5')] = qcdb.Molecule("""
0 1
C    0.000000000  -0.471925000   0.471925000
C    0.000000000   0.471925000  -0.471925000
H    0.922986000  -0.872422000   0.872422000
H    0.922986000   0.872422000  -0.872422000
H   -0.924197000  -0.870464000   0.870464000
H   -0.924197000   0.870464000  -0.870464000
--
0 1
C    5.577333000   0.471925000   0.471925000
C    5.577333000  -0.471925000  -0.471925000
H    4.654347000   0.872422000   0.872422000
H    4.654347000  -0.872422000  -0.872422000
H    6.501530000   0.870464000   0.870464000
H    6.501530000  -0.870464000  -0.870464000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '9-2.0')] = qcdb.Molecule("""
0 1
C    0.000000000  -0.471925000   0.471925000
C    0.000000000   0.471925000  -0.471925000
H    0.922986000  -0.872422000   0.872422000
H    0.922986000   0.872422000  -0.872422000
H   -0.924197000  -0.870464000   0.870464000
H   -0.924197000   0.870464000  -0.870464000
--
0 1
C    7.436444000   0.471925000   0.471925000
C    7.436444000  -0.471925000  -0.471925000
H    6.513458000   0.872422000   0.872422000
H    6.513458000  -0.872422000  -0.872422000
H    8.360641000   0.870464000   0.870464000
H    8.360641000  -0.870464000  -0.870464000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '10-0.7')] = qcdb.Molecule("""
0 1
C   1.10020000e-05   3.62910780e-02  -1.39321800e+00
C  -1.10750000e-05  -1.18840188e+00  -7.28035925e-01
C   1.09220000e-05  -1.22470779e+00   6.65180078e-01
C  -1.10020000e-05  -3.62967450e-02   1.39320400e+00
C   1.10750000e-05   1.18841621e+00   7.28037925e-01
C  -1.09220000e-05   1.22469913e+00  -6.65168078e-01
H   1.56700400e-03   6.44480100e-02  -2.47427400e+00
H   1.55086600e-03  -2.11054092e+00  -1.29295887e+00
H   1.56686200e-03  -2.17500776e+00   1.18132314e+00
H   1.55099600e-03  -6.44646770e-02   2.47426100e+00
H   1.56713400e-03   2.11056025e+00   1.29295087e+00
H   1.55113800e-03   2.17500609e+00  -1.18130314e+00
--
0 1
C   2.92692170e+00  -6.90000000e-08   0.00000000e+00
H   3.29067975e+00   8.38173871e-01  -5.86878053e-01
H   3.29067971e+00   8.91639730e-02   1.01931899e+00
H   1.84097270e+00   0.00000000e+00   0.00000000e+00
H   3.29067964e+00  -9.27338119e-01  -4.32440941e-01
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '10-0.8')] = qcdb.Molecule("""
0 1
C   1.10020000e-05   3.62910780e-02  -1.39321800e+00
C  -1.10750000e-05  -1.18840188e+00  -7.28035925e-01
C   1.09220000e-05  -1.22470779e+00   6.65180078e-01
C  -1.10020000e-05  -3.62967450e-02   1.39320400e+00
C   1.10750000e-05   1.18841621e+00   7.28037925e-01
C  -1.09220000e-05   1.22469913e+00  -6.65168078e-01
H   1.56700400e-03   6.44480100e-02  -2.47427400e+00
H   1.55086600e-03  -2.11054092e+00  -1.29295887e+00
H   1.56686200e-03  -2.17500776e+00   1.18132314e+00
H   1.55099600e-03  -6.44646770e-02   2.47426100e+00
H   1.56713400e-03   2.11056025e+00   1.29295087e+00
H   1.55113800e-03   2.17500609e+00  -1.18130314e+00
--
0 1
C   3.18991780e+00  -6.90000000e-08   0.00000000e+00
H   3.55367585e+00   8.38173871e-01  -5.86878053e-01
H   3.55367581e+00   8.91639730e-02   1.01931899e+00
H   2.10396880e+00   0.00000000e+00   0.00000000e+00
H   3.55367574e+00  -9.27338119e-01  -4.32440941e-01
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '10-0.9')] = qcdb.Molecule("""
0 1
C    0.000011002   0.036291078  -1.393218002
C   -0.000011075  -1.188401879  -0.728035925
C    0.000010922  -1.224707791   0.665180078
C   -0.000011002  -0.036296745   1.393204002
C    0.000011075   1.188416213   0.728037925
C   -0.000010922   1.224699125  -0.665168078
H    0.001567004   0.064448010  -2.474274004
H    0.001550866  -2.110540915  -1.292958866
H    0.001566862  -2.175007759   1.181323138
H    0.001550996  -0.064464677   2.474261004
H    0.001567134   2.110560249   1.292950866
H    0.001551138   2.175006092  -1.181303138
--
0 1
C    3.452913900  -0.000000069   0.000000000
H    3.816671953   0.838173871  -0.586878053
H    3.816671906   0.089163973   1.019318994
H    2.366964900   0.000000000   0.000000000
H    3.816671841  -0.927338119  -0.432440941
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '10-1.0')] = qcdb.Molecule("""
0 1
C    1.393218000   0.036291000  -0.633280000
C    0.728036000  -1.188402000  -0.633302000
C   -0.665180000  -1.224708000  -0.633280000
C   -1.393204000  -0.036297000  -0.633302000
C   -0.728038000   1.188416000  -0.633280000
C    0.665168000   1.224699000  -0.633302000
H    2.474274000   0.064448000  -0.631724000
H    1.292959000  -2.110541000  -0.631740000
H   -1.181323000  -2.175008000  -0.631724000
H   -2.474261000  -0.064465000  -0.631740000
H   -1.292951000   2.110560000  -0.631724000
H    1.181303000   2.175006000  -0.631740000
--
0 1
C    0.000000000   0.000000000   3.082619000
H    0.586878000   0.838174000   3.446377000
H   -1.019319000   0.089164000   3.446377000
H    0.000000000   0.000000000   1.996670000
H    0.432441000  -0.927338000   3.446377000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '10-1.2')] = qcdb.Molecule("""
0 1
C    0.000011002   0.036291078  -1.393218002
C   -0.000011075  -1.188401879  -0.728035925
C    0.000010922  -1.224707791   0.665180078
C   -0.000011002  -0.036296745   1.393204002
C    0.000011075   1.188416213   0.728037925
C   -0.000010922   1.224699125  -0.665168078
H    0.001567004   0.064448010  -2.474274004
H    0.001550866  -2.110540915  -1.292958866
H    0.001566862  -2.175007759   1.181323138
H    0.001550996  -0.064464677   2.474261004
H    0.001567134   2.110560249   1.292950866
H    0.001551138   2.175006092  -1.181303138
--
0 1
C    4.241902200  -0.000000069   0.000000000
H    4.605660253   0.838173871  -0.586878053
H    4.605660206   0.089163973   1.019318994
H    3.155953200   0.000000000   0.000000000
H    4.605660141  -0.927338119  -0.432440941
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '10-1.5')] = qcdb.Molecule("""
0 1
C    0.000011002   0.036291078  -1.393218002
C   -0.000011075  -1.188401879  -0.728035925
C    0.000010922  -1.224707791   0.665180078
C   -0.000011002  -0.036296745   1.393204002
C    0.000011075   1.188416213   0.728037925
C   -0.000010922   1.224699125  -0.665168078
H    0.001567004   0.064448010  -2.474274004
H    0.001550866  -2.110540915  -1.292958866
H    0.001566862  -2.175007759   1.181323138
H    0.001550996  -0.064464677   2.474261004
H    0.001567134   2.110560249   1.292950866
H    0.001551138   2.175006092  -1.181303138
--
0 1
C    5.030890500  -0.000000069   0.000000000
H    5.394648553   0.838173871  -0.586878053
H    5.394648506   0.089163973   1.019318994
H    3.944941500   0.000000000   0.000000000
H    5.394648441  -0.927338119  -0.432440941
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '10-2.0')] = qcdb.Molecule("""
0 1
C    0.000011002   0.036291078  -1.393218002
C   -0.000011075  -1.188401879  -0.728035925
C    0.000010922  -1.224707791   0.665180078
C   -0.000011002  -0.036296745   1.393204002
C    0.000011075   1.188416213   0.728037925
C   -0.000010922   1.224699125  -0.665168078
H    0.001567004   0.064448010  -2.474274004
H    0.001550866  -2.110540915  -1.292958866
H    0.001566862  -2.175007759   1.181323138
H    0.001550996  -0.064464677   2.474261004
H    0.001567134   2.110560249   1.292950866
H    0.001551138   2.175006092  -1.181303138
--
0 1
C    6.345871000  -0.000000069   0.000000000
H    6.709629053   0.838173871  -0.586878053
H    6.709629006   0.089163973   1.019318994
H    5.259922000   0.000000000   0.000000000
H    6.709628941  -0.927338119  -0.432440941
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '11-0.7')] = qcdb.Molecule("""
0 1
C          0.62905151 -1.24405848  0.000000
C          0.31407229 -0.62213466  1.206205
C          0.31407229 -0.62213466 -1.206205
C         -0.31481355  0.62169924  1.206954
C         -0.62756899  1.24492931  0.000000
C         -0.31481355  0.62169924 -1.206954
H          0.56393058 -1.10277815 -2.142315
H         -0.55938882  1.10408575 -2.143798
H         -1.11689412  2.20968592  0.000000
H         -0.55938882  1.10408575  2.143798
H          0.56393058 -1.10277815  2.142315
H          1.12972171 -2.20246266  0.000000
--
0 1
C          2.00660462  1.24405848  0.000000
C          2.32158383  0.62213466 -1.206205
C          2.32158383  0.62213466  1.206205
C          2.95046967 -0.62169924 -1.206954
C          3.26322512 -1.24492931  0.000000
C          2.95046967 -0.62169924  1.206954
H          1.50593441  2.20246266  0.000000
H          2.07172555  1.10277815  2.142315
H          3.19504494 -1.10408575  2.143798
H          3.75255025 -2.20968592  0.000000
H          3.19504494 -1.10408575 -2.143798
H          2.07172555  1.10277815 -2.142315
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '11-0.8')] = qcdb.Molecule("""
0 1
C          0.62905151 -1.24405848  0.000000
C          0.31407229 -0.62213466  1.206205
C          0.31407229 -0.62213466 -1.206205
C         -0.31481355  0.62169924  1.206954
C         -0.62756899  1.24492931  0.000000
C         -0.31481355  0.62169924 -1.206954
H          0.56393058 -1.10277815 -2.142315
H         -0.55938882  1.10408575 -2.143798
H         -1.11689412  2.20968592  0.000000
H         -0.55938882  1.10408575  2.143798
H          0.56393058 -1.10277815  2.142315
H          1.12972171 -2.20246266  0.000000
--
0 1
C          2.38312692  1.24405848  0.000000
C          2.69810614  0.62213466 -1.206205
C          2.69810614  0.62213466  1.206205
C          3.32699197 -0.62169924 -1.206954
C          3.63974742 -1.24492931  0.000000
C          3.32699197 -0.62169924  1.206954
H          1.88245672  2.20246266  0.000000
H          2.44824785  1.10277815  2.142315
H          3.57156725 -1.10408575  2.143798
H          4.12907255 -2.20968592  0.000000
H          3.57156725 -1.10408575 -2.143798
H          2.44824785  1.10277815 -2.142315
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '11-0.9')] = qcdb.Molecule("""
0 1
C    0.629051507  -1.244058476   0.000000000
C    0.314072291  -0.622134657   1.206205000
C    0.314072291  -0.622134657  -1.206205000
C   -0.314813547   0.621699240   1.206954000
C   -0.627568995   1.244929310   0.000000000
C   -0.314813547   0.621699240  -1.206954000
H    0.563930576  -1.102778154  -2.142315000
H   -0.559388819   1.104085746  -2.143798000
H   -1.116894124   2.209685917   0.000000000
H   -0.559388819   1.104085746   2.143798000
H    0.563930576  -1.102778154   2.142315000
H    1.129721711  -2.202462660   0.000000000
--
0 1
C    2.759649224   1.244058476   0.000000000
C    3.074628440   0.622134657  -1.206205000
C    3.074628440   0.622134657   1.206205000
C    3.703514278  -0.621699240  -1.206954000
C    4.016269727  -1.244929310   0.000000000
C    3.703514278  -0.621699240   1.206954000
H    2.258979020   2.202462660   0.000000000
H    2.824770156   1.102778154   2.142315000
H    3.948089550  -1.104085746   2.143798000
H    4.505594855  -2.209685917   0.000000000
H    3.948089550  -1.104085746  -2.143798000
H    2.824770156   1.102778154  -2.142315000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '11-1.0')] = qcdb.Molecule("""
0 1
C   -1.047825000  -1.421674000   0.000000000
C   -1.454503000  -0.855446000   1.206205000
C   -1.454503000  -0.855446000  -1.206205000
C   -2.266797000   0.277161000   1.206954000
C   -2.671478000   0.845021000   0.000000000
C   -2.266797000   0.277161000  -1.206954000
H   -1.133853000  -1.292059000  -2.142315000
H   -2.582494000   0.716307000  -2.143798000
H   -3.303042000   1.723270000   0.000000000
H   -2.582494000   0.716307000   2.143798000
H   -1.133853000  -1.292059000   2.142315000
H   -0.406025000  -2.291905000   0.000000000
--
0 1
C    1.047825000   1.421674000   0.000000000
C    1.454503000   0.855446000  -1.206205000
C    1.454503000   0.855446000   1.206205000
C    2.266797000  -0.277161000  -1.206954000
C    2.671478000  -0.845021000   0.000000000
C    2.266797000  -0.277161000   1.206954000
H    0.406025000   2.291905000   0.000000000
H    1.133853000   1.292059000   2.142315000
H    2.582494000  -0.716307000   2.143798000
H    3.303042000  -1.723270000   0.000000000
H    2.582494000  -0.716307000  -2.143798000
H    1.133853000   1.292059000  -2.142315000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '11-1.2')] = qcdb.Molecule("""
0 1
C    0.629051507  -1.244058476   0.000000000
C    0.314072291  -0.622134657   1.206205000
C    0.314072291  -0.622134657  -1.206205000
C   -0.314813547   0.621699240   1.206954000
C   -0.627568995   1.244929310   0.000000000
C   -0.314813547   0.621699240  -1.206954000
H    0.563930576  -1.102778154  -2.142315000
H   -0.559388819   1.104085746  -2.143798000
H   -1.116894124   2.209685917   0.000000000
H   -0.559388819   1.104085746   2.143798000
H    0.563930576  -1.102778154   2.142315000
H    1.129721711  -2.202462660   0.000000000
--
0 1
C    3.889216135   1.244058476   0.000000000
C    4.204195351   0.622134657  -1.206205000
C    4.204195351   0.622134657   1.206205000
C    4.833081189  -0.621699240  -1.206954000
C    5.145836638  -1.244929310   0.000000000
C    4.833081189  -0.621699240   1.206954000
H    3.388545931   2.202462660   0.000000000
H    3.954337067   1.102778154   2.142315000
H    5.077656461  -1.104085746   2.143798000
H    5.635161766  -2.209685917   0.000000000
H    5.077656461  -1.104085746  -2.143798000
H    3.954337067   1.102778154  -2.142315000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '11-1.5')] = qcdb.Molecule("""
0 1
C    0.629051507  -1.244058476   0.000000000
C    0.314072291  -0.622134657   1.206205000
C    0.314072291  -0.622134657  -1.206205000
C   -0.314813547   0.621699240   1.206954000
C   -0.627568995   1.244929310   0.000000000
C   -0.314813547   0.621699240  -1.206954000
H    0.563930576  -1.102778154  -2.142315000
H   -0.559388819   1.104085746  -2.143798000
H   -1.116894124   2.209685917   0.000000000
H   -0.559388819   1.104085746   2.143798000
H    0.563930576  -1.102778154   2.142315000
H    1.129721711  -2.202462660   0.000000000
--
0 1
C    5.018783046   1.244058476   0.000000000
C    5.333762262   0.622134657  -1.206205000
C    5.333762262   0.622134657   1.206205000
C    5.962648100  -0.621699240  -1.206954000
C    6.275403549  -1.244929310   0.000000000
C    5.962648100  -0.621699240   1.206954000
H    4.518112842   2.202462660   0.000000000
H    5.083903978   1.102778154   2.142315000
H    6.207223372  -1.104085746   2.143798000
H    6.764728677  -2.209685917   0.000000000
H    6.207223372  -1.104085746  -2.143798000
H    5.083903978   1.102778154  -2.142315000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '11-2.0')] = qcdb.Molecule("""
0 1
C    0.629051507  -1.244058476   0.000000000
C    0.314072291  -0.622134657   1.206205000
C    0.314072291  -0.622134657  -1.206205000
C   -0.314813547   0.621699240   1.206954000
C   -0.627568995   1.244929310   0.000000000
C   -0.314813547   0.621699240  -1.206954000
H    0.563930576  -1.102778154  -2.142315000
H   -0.559388819   1.104085746  -2.143798000
H   -1.116894124   2.209685917   0.000000000
H   -0.559388819   1.104085746   2.143798000
H    0.563930576  -1.102778154   2.142315000
H    1.129721711  -2.202462660   0.000000000
--
0 1
C    6.901394563   1.244058476   0.000000000
C    7.216373779   0.622134657  -1.206205000
C    7.216373779   0.622134657   1.206205000
C    7.845259617  -0.621699240  -1.206954000
C    8.158015066  -1.244929310   0.000000000
C    7.845259617  -0.621699240   1.206954000
H    6.400724359   2.202462660   0.000000000
H    6.966515495   1.102778154   2.142315000
H    8.089834889  -1.104085746   2.143798000
H    8.647340194  -2.209685917   0.000000000
H    8.089834889  -1.104085746  -2.143798000
H    6.966515495   1.102778154  -2.142315000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '12-0.7')] = qcdb.Molecule("""
0 1
C   3.95653045e-01   1.05943214e+00  -6.96139000e-01
C   3.95653045e-01   1.05943214e+00   6.96139000e-01
N  -3.26335700e-03   2.27377000e-04   1.41448000e+00
C  -3.91847355e-01  -1.05969731e+00   6.96729000e-01
C  -3.91847355e-01  -1.05969731e+00  -6.96729000e-01
N  -3.26335700e-03   2.27377000e-04  -1.41448000e+00
H   7.18983381e-01   1.93337025e+00  -1.24728000e+00
H   7.18983381e-01   1.93337025e+00   1.24728000e+00
H  -7.13152254e-01  -1.93436275e+00   1.24756000e+00
H  -7.13152254e-01  -1.93436275e+00  -1.24756000e+00
--
0 1
C   2.70265202e+00   6.43131999e-01   1.13004500e+00
C   2.16690706e+00  -6.42689433e-01   1.13063100e+00
N   1.89388599e+00  -1.30673885e+00   0.00000000e+00
C   2.16690706e+00  -6.42689433e-01  -1.13063100e+00
C   2.70265202e+00   6.43131999e-01  -1.13004500e+00
N   2.98013696e+00   1.30597985e+00   0.00000000e+00
H   2.91361017e+00   1.15247121e+00   2.06186400e+00
H   1.94717154e+00  -1.14774434e+00   2.06239900e+00
H   1.94717154e+00  -1.14774434e+00  -2.06239900e+00
H   2.91361017e+00   1.15247121e+00  -2.06186400e+00
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '12-0.8')] = qcdb.Molecule("""
0 1
C   3.95653045e-01   1.05943214e+00  -6.96139000e-01
C   3.95653045e-01   1.05943214e+00   6.96139000e-01
N  -3.26335700e-03   2.27377000e-04   1.41448000e+00
C  -3.91847355e-01  -1.05969731e+00   6.96729000e-01
C  -3.91847355e-01  -1.05969731e+00  -6.96729000e-01
N  -3.26335700e-03   2.27377000e-04  -1.41448000e+00
H   7.18983381e-01   1.93337025e+00  -1.24728000e+00
H   7.18983381e-01   1.93337025e+00   1.24728000e+00
H  -7.13152254e-01  -1.93436275e+00   1.24756000e+00
H  -7.13152254e-01  -1.93436275e+00  -1.24756000e+00
--
0 1
C   3.05059511e+00   6.43131999e-01   1.13004500e+00
C   2.51485015e+00  -6.42689433e-01   1.13063100e+00
N   2.24182908e+00  -1.30673885e+00   0.00000000e+00
C   2.51485015e+00  -6.42689433e-01  -1.13063100e+00
C   3.05059511e+00   6.43131999e-01  -1.13004500e+00
N   3.32808005e+00   1.30597985e+00   0.00000000e+00
H   3.26155326e+00   1.15247121e+00   2.06186400e+00
H   2.29511463e+00  -1.14774434e+00   2.06239900e+00
H   2.29511463e+00  -1.14774434e+00  -2.06239900e+00
H   3.26155326e+00   1.15247121e+00  -2.06186400e+00
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '12-0.9')] = qcdb.Molecule("""
0 1
C    0.395653045   1.059432142  -0.696139000
C    0.395653045   1.059432142   0.696139000
N   -0.003263357   0.000227377   1.414480000
C   -0.391847355  -1.059697307   0.696729000
C   -0.391847355  -1.059697307  -0.696729000
N   -0.003263357   0.000227377  -1.414480000
H    0.718983381   1.933370245  -1.247280000
H    0.718983381   1.933370245   1.247280000
H   -0.713152254  -1.934362753   1.247560000
H   -0.713152254  -1.934362753  -1.247560000
--
0 1
C    3.398538200   0.643131999   1.130045000
C    2.862793235  -0.642689433   1.130631000
N    2.589772167  -1.306738847   0.000000000
C    2.862793235  -0.642689433  -1.130631000
C    3.398538200   0.643131999  -1.130045000
N    3.676023139   1.305979850   0.000000000
H    3.609496345   1.152471205   2.061864000
H    2.643057716  -1.147744338   2.062399000
H    2.643057716  -1.147744338  -2.062399000
H    3.609496345   1.152471205  -2.061864000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '12-1.0')] = qcdb.Molecule("""
0 1
C   -1.247189000  -1.171821000  -0.696139000
C   -1.247189000  -1.171821000   0.696139000
N   -0.258951000  -1.723577000   1.414480000
C    0.731533000  -2.265222000   0.696729000
C    0.731533000  -2.265222000  -0.696729000
N   -0.258951000  -1.723577000  -1.414480000
H   -2.063436000  -0.722320000  -1.247280000
H   -2.063436000  -0.722320000   1.247280000
H    1.548800000  -2.712828000   1.247560000
H    1.548800000  -2.712828000  -1.247560000
--
0 1
C   -0.338003000   2.080061000   1.130045000
C    0.854025000   1.359347000   1.130631000
N    1.470179000   0.990760000   0.000000000
C    0.854025000   1.359347000  -1.130631000
C   -0.338003000   2.080061000  -1.130045000
N   -0.952306000   2.452884000   0.000000000
H   -0.810376000   2.364303000   2.061864000
H    1.320858000   1.067061000   2.062399000
H    1.320858000   1.067061000  -2.062399000
H   -0.810376000   2.364303000  -2.061864000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '12-1.2')] = qcdb.Molecule("""
0 1
C    0.395653045   1.059432142  -0.696139000
C    0.395653045   1.059432142   0.696139000
N   -0.003263357   0.000227377   1.414480000
C   -0.391847355  -1.059697307   0.696729000
C   -0.391847355  -1.059697307  -0.696729000
N   -0.003263357   0.000227377  -1.414480000
H    0.718983381   1.933370245  -1.247280000
H    0.718983381   1.933370245   1.247280000
H   -0.713152254  -1.934362753   1.247560000
H   -0.713152254  -1.934362753  -1.247560000
--
0 1
C    4.442367465   0.643131999   1.130045000
C    3.906622500  -0.642689433   1.130631000
N    3.633601432  -1.306738847   0.000000000
C    3.906622500  -0.642689433  -1.130631000
C    4.442367465   0.643131999  -1.130045000
N    4.719852404   1.305979850   0.000000000
H    4.653325610   1.152471205   2.061864000
H    3.686886981  -1.147744338   2.062399000
H    3.686886981  -1.147744338  -2.062399000
H    4.653325610   1.152471205  -2.061864000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '12-1.5')] = qcdb.Molecule("""
0 1
C    0.395653045   1.059432142  -0.696139000
C    0.395653045   1.059432142   0.696139000
N   -0.003263357   0.000227377   1.414480000
C   -0.391847355  -1.059697307   0.696729000
C   -0.391847355  -1.059697307  -0.696729000
N   -0.003263357   0.000227377  -1.414480000
H    0.718983381   1.933370245  -1.247280000
H    0.718983381   1.933370245   1.247280000
H   -0.713152254  -1.934362753   1.247560000
H   -0.713152254  -1.934362753  -1.247560000
--
0 1
C    5.486196730   0.643131999   1.130045000
C    4.950451765  -0.642689433   1.130631000
N    4.677430697  -1.306738847   0.000000000
C    4.950451765  -0.642689433  -1.130631000
C    5.486196730   0.643131999  -1.130045000
N    5.763681669   1.305979850   0.000000000
H    5.697154875   1.152471205   2.061864000
H    4.730716246  -1.147744338   2.062399000
H    4.730716246  -1.147744338  -2.062399000
H    5.697154875   1.152471205  -2.061864000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '12-2.0')] = qcdb.Molecule("""
0 1
C    0.395653045   1.059432142  -0.696139000
C    0.395653045   1.059432142   0.696139000
N   -0.003263357   0.000227377   1.414480000
C   -0.391847355  -1.059697307   0.696729000
C   -0.391847355  -1.059697307  -0.696729000
N   -0.003263357   0.000227377  -1.414480000
H    0.718983381   1.933370245  -1.247280000
H    0.718983381   1.933370245   1.247280000
H   -0.713152254  -1.934362753   1.247560000
H   -0.713152254  -1.934362753  -1.247560000
--
0 1
C    7.225912172   0.643131999   1.130045000
C    6.690167207  -0.642689433   1.130631000
N    6.417146139  -1.306738847   0.000000000
C    6.690167207  -0.642689433  -1.130631000
C    7.225912172   0.643131999  -1.130045000
N    7.503397111   1.305979850   0.000000000
H    7.436870317   1.152471205   2.061864000
H    6.470431688  -1.147744338   2.062399000
H    6.470431688  -1.147744338  -2.062399000
H    7.436870317   1.152471205  -2.061864000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '13-0.7')] = qcdb.Molecule("""
0 1
N         -0.27790501  1.29367954  0.17614197
C         -0.31314340  0.77865720 -1.09019403
H         -0.55662845  1.48297630 -1.87143703
C         -0.05442933 -0.52203414 -1.33828003
H         -0.08333918 -0.92007182 -2.33779603
C          0.31574183 -1.40331977 -0.24638003
O          0.65706663 -2.57165556 -0.35183703
N          0.27289252 -0.78328638  1.00884497
H          0.57557519 -1.34248314  1.79757997
C          0.05767640  0.55148208  1.29293597
O          0.16219780  1.03423971  2.40401497
H         -0.35588204  2.28595021  0.33102197
--
0 1
N          2.63363413 -1.29367954  0.17614197
C          2.66887252 -0.77865720 -1.09019403
H          2.91235758 -1.48297630 -1.87143703
C          2.41015845  0.52203414 -1.33828003
H          2.43906830  0.92007182 -2.33779603
C          2.03998729  1.40331977 -0.24638003
O          1.69866249  2.57165556 -0.35183703
N          2.08283661  0.78328638  1.00884497
H          1.78015393  1.34248314  1.79757997
C          2.29805272 -0.55148208  1.29293597
O          2.19353133 -1.03423971  2.40401497
H          2.71161116 -2.28595021  0.33102197
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '13-0.8')] = qcdb.Molecule("""
0 1
N         -0.27790501  1.29367954  0.17614197
C         -0.31314340  0.77865720 -1.09019403
H         -0.55662845  1.48297630 -1.87143703
C         -0.05442933 -0.52203414 -1.33828003
H         -0.08333918 -0.92007182 -2.33779603
C          0.31574183 -1.40331977 -0.24638003
O          0.65706663 -2.57165556 -0.35183703
N          0.27289252 -0.78328638  1.00884497
H          0.57557519 -1.34248314  1.79757997
C          0.05767640  0.55148208  1.29293597
O          0.16219780  1.03423971  2.40401497
H         -0.35588204  2.28595021  0.33102197
--
0 1
N          2.97016686 -1.29367954  0.17614197
C          3.00540525 -0.77865720 -1.09019403
H          3.24889031 -1.48297630 -1.87143703
C          2.74669118  0.52203414 -1.33828003
H          2.77560103  0.92007182 -2.33779603
C          2.37652002  1.40331977 -0.24638003
O          2.03519522  2.57165556 -0.35183703
N          2.41936934  0.78328638  1.00884497
H          2.11668667  1.34248314  1.79757997
C          2.63458546 -0.55148208  1.29293597
O          2.53006406 -1.03423971  2.40401497
H          3.04814390 -2.28595021  0.33102197
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '13-0.9')] = qcdb.Molecule("""
0 1
N   -0.277905006   1.293679543   0.176141970
C   -0.313143400   0.778657200  -1.090194030
H   -0.556628453   1.482976305  -1.871437030
C   -0.054429325  -0.522034140  -1.338280030
H   -0.083339176  -0.920071815  -2.337796030
C    0.315741834  -1.403319766  -0.246380030
O    0.657066634  -2.571655559  -0.351837030
N    0.272892517  -0.783286382   1.008844970
H    0.575575188  -1.342483138   1.797579970
C    0.057676398   0.551482081   1.292935970
O    0.162197796   1.034239706   2.404014970
H   -0.355882042   2.285950208   0.331021970
--
0 1
N    3.306699593  -1.293679543   0.176141970
C    3.341937987  -0.778657200  -1.090194030
H    3.585423040  -1.482976305  -1.871437030
C    3.083223911   0.522034140  -1.338280030
H    3.112133763   0.920071815  -2.337796030
C    2.713052753   1.403319766  -0.246380030
O    2.371727953   2.571655559  -0.351837030
N    2.755902070   0.783286382   1.008844970
H    2.453219399   1.342483138   1.797579970
C    2.971118189  -0.551482081   1.292935970
O    2.866596791  -1.034239706   2.404014970
H    3.384676629  -2.285950208   0.331021970
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '13-1.0')] = qcdb.Molecule("""
0 1
N    2.011359000  -1.213207000  -0.098067000
C    2.025708000  -0.697180000  -1.364403000
H    2.297521000  -1.391059000  -2.145646000
C    1.714523000   0.591965000  -1.612489000
H    1.727287000   0.990847000  -2.612005000
C    1.308960000   1.457534000  -0.520589000
O    0.920593000   2.611086000  -0.626046000
N    1.376888000   0.839745000   0.734636000
H    1.051804000   1.386223000   1.523371000
C    1.645991000  -0.485211000   1.018727000
O    1.561109000  -0.971806000   2.129806000
H    2.129463000  -2.201505000   0.056813000
--
0 1
N   -2.011359000   1.213207000  -0.098067000
C   -2.025708000   0.697180000  -1.364403000
H   -2.297521000   1.391059000  -2.145646000
C   -1.714523000  -0.591965000  -1.612489000
H   -1.727287000  -0.990847000  -2.612005000
C   -1.308960000  -1.457534000  -0.520589000
O   -0.920593000  -2.611086000  -0.626046000
N   -1.376888000  -0.839745000   0.734636000
H   -1.051804000  -1.386223000   1.523371000
C   -1.645991000   0.485211000   1.018727000
O   -1.561109000   0.971806000   2.129806000
H   -2.129463000   2.201505000   0.056813000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '13-1.2')] = qcdb.Molecule("""
0 1
N   -0.277905006   1.293679543   0.176141970
C   -0.313143400   0.778657200  -1.090194030
H   -0.556628453   1.482976305  -1.871437030
C   -0.054429325  -0.522034140  -1.338280030
H   -0.083339176  -0.920071815  -2.337796030
C    0.315741834  -1.403319766  -0.246380030
O    0.657066634  -2.571655559  -0.351837030
N    0.272892517  -0.783286382   1.008844970
H    0.575575188  -1.342483138   1.797579970
C    0.057676398   0.551482081   1.292935970
O    0.162197796   1.034239706   2.404014970
H   -0.355882042   2.285950208   0.331021970
--
0 1
N    4.316297789  -1.293679543   0.176141970
C    4.351536183  -0.778657200  -1.090194030
H    4.595021236  -1.482976305  -1.871437030
C    4.092822107   0.522034140  -1.338280030
H    4.121731959   0.920071815  -2.337796030
C    3.722650949   1.403319766  -0.246380030
O    3.381326149   2.571655559  -0.351837030
N    3.765500266   0.783286382   1.008844970
H    3.462817595   1.342483138   1.797579970
C    3.980716385  -0.551482081   1.292935970
O    3.876194987  -1.034239706   2.404014970
H    4.394274825  -2.285950208   0.331021970
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '13-1.5')] = qcdb.Molecule("""
0 1
N   -0.277905006   1.293679543   0.176141970
C   -0.313143400   0.778657200  -1.090194030
H   -0.556628453   1.482976305  -1.871437030
C   -0.054429325  -0.522034140  -1.338280030
H   -0.083339176  -0.920071815  -2.337796030
C    0.315741834  -1.403319766  -0.246380030
O    0.657066634  -2.571655559  -0.351837030
N    0.272892517  -0.783286382   1.008844970
H    0.575575188  -1.342483138   1.797579970
C    0.057676398   0.551482081   1.292935970
O    0.162197796   1.034239706   2.404014970
H   -0.355882042   2.285950208   0.331021970
--
0 1
N    5.325895984  -1.293679543   0.176141970
C    5.361134378  -0.778657200  -1.090194030
H    5.604619431  -1.482976305  -1.871437030
C    5.102420302   0.522034140  -1.338280030
H    5.131330154   0.920071815  -2.337796030
C    4.732249144   1.403319766  -0.246380030
O    4.390924344   2.571655559  -0.351837030
N    4.775098461   0.783286382   1.008844970
H    4.472415790   1.342483138   1.797579970
C    4.990314580  -0.551482081   1.292935970
O    4.885793182  -1.034239706   2.404014970
H    5.403873020  -2.285950208   0.331021970
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '13-2.0')] = qcdb.Molecule("""
0 1
N   -0.277905006   1.293679543   0.176141970
C   -0.313143400   0.778657200  -1.090194030
H   -0.556628453   1.482976305  -1.871437030
C   -0.054429325  -0.522034140  -1.338280030
H   -0.083339176  -0.920071815  -2.337796030
C    0.315741834  -1.403319766  -0.246380030
O    0.657066634  -2.571655559  -0.351837030
N    0.272892517  -0.783286382   1.008844970
H    0.575575188  -1.342483138   1.797579970
C    0.057676398   0.551482081   1.292935970
O    0.162197796   1.034239706   2.404014970
H   -0.355882042   2.285950208   0.331021970
--
0 1
N    7.008559644  -1.293679543   0.176141970
C    7.043798038  -0.778657200  -1.090194030
H    7.287283091  -1.482976305  -1.871437030
C    6.785083962   0.522034140  -1.338280030
H    6.813993814   0.920071815  -2.337796030
C    6.414912804   1.403319766  -0.246380030
O    6.073588004   2.571655559  -0.351837030
N    6.457762121   0.783286382   1.008844970
H    6.155079450   1.342483138   1.797579970
C    6.672978240  -0.551482081   1.292935970
O    6.568456842  -1.034239706   2.404014970
H    7.086536680  -2.285950208   0.331021970
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '14-0.7')] = qcdb.Molecule("""
0 1
C            0.00000000  0.00000000  0.00000000
C           -0.04448565 -1.17797863  0.74316010
C           -0.01082464 -2.41120852  0.09533314
C            0.06415077 -2.46693379 -1.29562360
C            0.10095090 -1.28743705 -2.03895997
C            0.06735680 -0.05350021 -1.39137626
H           -0.01379774  0.95688159  0.50334833
H           -0.09134697 -1.13445800  1.82239892
H           -0.03975401 -3.32568027  0.67235867
H            0.08538953 -3.42484902 -1.79837382
H            0.14644278 -1.33017254 -3.11951477
H            0.10085283  0.86245624 -1.96494557
--
0 1
H            2.06737101 -0.57805685  3.49490475
C            2.14311338 -0.57196987  2.41575396
C            2.10265932  0.63365013  1.73434956
H            1.99554084  1.56703853  2.27203610
C            2.20540983  0.62434756  0.33333966
C            2.19524253  1.63366203 -0.67349928
H            2.11161860  2.69803059 -0.53325175
C            2.32582959  0.99280815 -1.88451747
N            2.43153522 -0.36008660 -1.67542289
C            2.34735531 -0.62434756 -0.33333966
C            2.39589311 -1.83984299  0.35175494
H            2.50271193 -2.78021794 -0.17294023
C            2.29112185 -1.79621168  1.73303617
H            2.32275342 -2.71826144  2.29763493
H            2.45348129 -1.05644621 -2.39897877
H            2.36204661  1.39803628 -2.88180774
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '14-0.8')] = qcdb.Molecule("""
0 1
C            0.00000000  0.00000000  0.00000000
C           -0.04448565 -1.17797863  0.74316010
C           -0.01082464 -2.41120852  0.09533314
C            0.06415077 -2.46693379 -1.29562360
C            0.10095090 -1.28743705 -2.03895997
C            0.06735680 -0.05350021 -1.39137626
H           -0.01379774  0.95688159  0.50334833
H           -0.09134697 -1.13445800  1.82239892
H           -0.03975401 -3.32568027  0.67235867
H            0.08538953 -3.42484902 -1.79837382
H            0.14644278 -1.33017254 -3.11951477
H            0.10085283  0.86245624 -1.96494557
--
0 1
H            2.39256852 -0.57805685  3.49490475
C            2.46831089 -0.57196987  2.41575396
C            2.42785683  0.63365013  1.73434956
H            2.32073835  1.56703853  2.27203610
C            2.53060734  0.62434756  0.33333966
C            2.52044004  1.63366203 -0.67349928
H            2.43681611  2.69803059 -0.53325175
C            2.65102710  0.99280815 -1.88451747
N            2.75673273 -0.36008660 -1.67542289
C            2.67255282 -0.62434756 -0.33333966
C            2.72109062 -1.83984299  0.35175494
H            2.82790944 -2.78021794 -0.17294023
C            2.61631936 -1.79621168  1.73303617
H            2.64795093 -2.71826144  2.29763493
H            2.77867880 -1.05644621 -2.39897877
H            2.68724412  1.39803628 -2.88180774
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '14-0.9')] = qcdb.Molecule("""
0 1
C    0.000000000   0.000000000   0.000000000
C   -0.044485647  -1.177978626   0.743160105
C   -0.010824638  -2.411208517   0.095333145
C    0.064150773  -2.466933785  -1.295623602
C    0.100950904  -1.287437054  -2.038959973
C    0.067356799  -0.053500209  -1.391376263
H   -0.013797739   0.956881587   0.503348328
H   -0.091346970  -1.134458005   1.822398921
H   -0.039754009  -3.325680275   0.672358669
H    0.085389531  -3.424849020  -1.798373823
H    0.146442780  -1.330172544  -3.119514770
H    0.100852832   0.862456237  -1.964945566
--
0 1
H    2.717766027  -0.578056849   3.494904751
C    2.793508398  -0.571969873   2.415753956
C    2.753054336   0.633650134   1.734349558
H    2.645935858   1.567038531   2.272036098
C    2.855804852   0.624347564   0.333339655
C    2.845637545   1.633662034  -0.673499279
H    2.762013625   2.698030593  -0.533251753
C    2.976224608   0.992808148  -1.884517470
N    3.081930238  -0.360086596  -1.675422891
C    2.997750328  -0.624347564  -0.333339655
C    3.046288127  -1.839842986   0.351754941
H    3.153106953  -2.780217935  -0.172940228
C    2.941516868  -1.796211682   1.733036170
H    2.973148444  -2.718261443   2.297634930
H    3.103876306  -1.056446212  -2.398978775
H    3.012441631   1.398036276  -2.881807744
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '14-1.0')] = qcdb.Molecule("""
0 1
C   -0.021074000   1.531861000  -1.363935000
C   -1.274679000   0.974103000  -1.607410000
C   -1.378305000  -0.225698000  -2.308415000
C   -0.228943000  -0.866405000  -2.768794000
C    1.024788000  -0.303517000  -2.531241000
C    1.129000000   0.896679000  -1.829983000
H    0.060074000   2.456563000  -0.809396000
H   -2.165100000   1.465452000  -1.240568000
H   -2.350973000  -0.661612000  -2.492670000
H   -0.310342000  -1.795576000  -3.317270000
H    1.916585000  -0.794084000  -2.899394000
H    2.100035000   1.332676000  -1.640042000
--
0 1
H   -2.941765000   0.895383000   2.223905000
C   -2.022067000   0.425854000   1.901355000
C   -0.814942000   1.074045000   2.106698000
H   -0.785153000   2.044381000   2.585609000
C    0.370429000   0.449285000   1.684746000
C    1.750862000   0.803894000   1.719400000
H    2.187011000   1.699828000   2.127590000
C    2.445136000  -0.231074000   1.135331000
N    1.564646000  -1.213781000   0.755538000
C    0.286121000  -0.826949000   1.061875000
C   -0.928467000  -1.485312000   0.860694000
H   -0.972920000  -2.455485000   0.383401000
C   -2.079285000  -0.841767000   1.287644000
H   -3.038997000  -1.320385000   1.146840000
H    1.807574000  -2.036696000   0.233304000
H    3.502879000  -0.348534000   0.969523000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '14-1.2')] = qcdb.Molecule("""
0 1
C    0.000000000   0.000000000   0.000000000
C   -0.044485647  -1.177978626   0.743160105
C   -0.010824638  -2.411208517   0.095333145
C    0.064150773  -2.466933785  -1.295623602
C    0.100950904  -1.287437054  -2.038959973
C    0.067356799  -0.053500209  -1.391376263
H   -0.013797739   0.956881587   0.503348328
H   -0.091346970  -1.134458005   1.822398921
H   -0.039754009  -3.325680275   0.672358669
H    0.085389531  -3.424849020  -1.798373823
H    0.146442780  -1.330172544  -3.119514770
H    0.100852832   0.862456237  -1.964945566
--
0 1
H    3.693358557  -0.578056849   3.494904751
C    3.769100928  -0.571969873   2.415753956
C    3.728646866   0.633650134   1.734349558
H    3.621528388   1.567038531   2.272036098
C    3.831397382   0.624347564   0.333339655
C    3.821230075   1.633662034  -0.673499279
H    3.737606155   2.698030593  -0.533251753
C    3.951817138   0.992808148  -1.884517470
N    4.057522768  -0.360086596  -1.675422891
C    3.973342858  -0.624347564  -0.333339655
C    4.021880657  -1.839842986   0.351754941
H    4.128699483  -2.780217935  -0.172940228
C    3.917109398  -1.796211682   1.733036170
H    3.948740974  -2.718261443   2.297634930
H    4.079468836  -1.056446212  -2.398978775
H    3.988034161   1.398036276  -2.881807744
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '14-1.5')] = qcdb.Molecule("""
0 1
C    0.000000000   0.000000000   0.000000000
C   -0.044485647  -1.177978626   0.743160105
C   -0.010824638  -2.411208517   0.095333145
C    0.064150773  -2.466933785  -1.295623602
C    0.100950904  -1.287437054  -2.038959973
C    0.067356799  -0.053500209  -1.391376263
H   -0.013797739   0.956881587   0.503348328
H   -0.091346970  -1.134458005   1.822398921
H   -0.039754009  -3.325680275   0.672358669
H    0.085389531  -3.424849020  -1.798373823
H    0.146442780  -1.330172544  -3.119514770
H    0.100852832   0.862456237  -1.964945566
--
0 1
H    4.668951087  -0.578056849   3.494904751
C    4.744693458  -0.571969873   2.415753956
C    4.704239396   0.633650134   1.734349558
H    4.597120918   1.567038531   2.272036098
C    4.806989912   0.624347564   0.333339655
C    4.796822605   1.633662034  -0.673499279
H    4.713198685   2.698030593  -0.533251753
C    4.927409668   0.992808148  -1.884517470
N    5.033115298  -0.360086596  -1.675422891
C    4.948935388  -0.624347564  -0.333339655
C    4.997473187  -1.839842986   0.351754941
H    5.104292013  -2.780217935  -0.172940228
C    4.892701928  -1.796211682   1.733036170
H    4.924333504  -2.718261443   2.297634930
H    5.055061366  -1.056446212  -2.398978775
H    4.963626691   1.398036276  -2.881807744
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '14-2.0')] = qcdb.Molecule("""
0 1
C    0.000000000   0.000000000   0.000000000
C   -0.044485647  -1.177978626   0.743160105
C   -0.010824638  -2.411208517   0.095333145
C    0.064150773  -2.466933785  -1.295623602
C    0.100950904  -1.287437054  -2.038959973
C    0.067356799  -0.053500209  -1.391376263
H   -0.013797739   0.956881587   0.503348328
H   -0.091346970  -1.134458005   1.822398921
H   -0.039754009  -3.325680275   0.672358669
H    0.085389531  -3.424849020  -1.798373823
H    0.146442780  -1.330172544  -3.119514770
H    0.100852832   0.862456237  -1.964945566
--
0 1
H    6.294938637  -0.578056849   3.494904751
C    6.370681008  -0.571969873   2.415753956
C    6.330226946   0.633650134   1.734349558
H    6.223108468   1.567038531   2.272036098
C    6.432977462   0.624347564   0.333339655
C    6.422810155   1.633662034  -0.673499279
H    6.339186235   2.698030593  -0.533251753
C    6.553397218   0.992808148  -1.884517470
N    6.659102848  -0.360086596  -1.675422891
C    6.574922938  -0.624347564  -0.333339655
C    6.623460737  -1.839842986   0.351754941
H    6.730279563  -2.780217935  -0.172940228
C    6.518689478  -1.796211682   1.733036170
H    6.550321054  -2.718261443   2.297634930
H    6.681048916  -1.056446212  -2.398978775
H    6.589614241   1.398036276  -2.881807744
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '15-0.7')] = qcdb.Molecule("""
0 1
N          0.06739076  1.21380610 -1.17119251
C         -0.03444069  0.16091603 -2.03517969
H         -0.03790910  0.30769467 -3.10231144
N         -0.12228650 -1.01421449 -1.43165939
C         -0.06127815 -0.69015606 -0.09773853
C         -0.08386647 -1.48000643  1.06512198
N         -0.20755129 -2.83016786  1.00846628
H          0.02023600 -3.31829451  1.85849278
H          0.10082398 -3.26183982  0.15179183
N         -0.01510729 -0.87288624  2.25482044
C          0.09553444  0.46847359  2.28659214
H          0.14844366  0.90243354  3.27705554
N          0.15079163  1.33081754  1.26823241
C          0.06127815  0.69015606  0.09773853
H          0.21312382  2.17853204 -1.42008256
--
0 1
N          2.35826983  1.31891257  0.11516933
C          2.39658658  0.54413479  1.24823546
H          2.52974924  1.08421646  2.17449125
C          2.27593596 -0.80203603  1.21330635
C          2.32838658 -1.66422779  2.42938073
H          1.37260336 -2.16186744  2.58503772
H          3.08922865 -2.43503398  2.31548757
H          2.55194105 -1.07062898  3.31353818
C          2.08145720 -1.44032645 -0.08037966
O          1.92105789 -2.64008185 -0.25503382
N          2.09265212 -0.56083789 -1.16848449
H          1.91696323 -0.97799874 -2.07261756
C          2.17759451  0.81416973 -1.15279815
O          2.09492605  1.51385406 -2.14916326
H          2.39663592  2.32251674  0.17911856
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '15-0.8')] = qcdb.Molecule("""
0 1
N          0.06739076  1.21380610 -1.17119251
C         -0.03444069  0.16091603 -2.03517969
H         -0.03790910  0.30769467 -3.10231144
N         -0.12228650 -1.01421449 -1.43165939
C         -0.06127815 -0.69015606 -0.09773853
C         -0.08386647 -1.48000643  1.06512198
N         -0.20755129 -2.83016786  1.00846628
H          0.02023600 -3.31829451  1.85849278
H          0.10082398 -3.26183982  0.15179183
N         -0.01510729 -0.87288624  2.25482044
C          0.09553444  0.46847359  2.28659214
H          0.14844366  0.90243354  3.27705554
N          0.15079163  1.33081754  1.26823241
C          0.06127815  0.69015606  0.09773853
H          0.21312382  2.17853204 -1.42008256
--
0 1
N          2.67686354  1.31891257  0.11516933
C          2.71518029  0.54413479  1.24823546
H          2.84834294  1.08421646  2.17449125
C          2.59452966 -0.80203603  1.21330635
C          2.64698029 -1.66422779  2.42938073
H          1.69119707 -2.16186744  2.58503772
H          3.40782236 -2.43503398  2.31548757
H          2.87053476 -1.07062898  3.31353818
C          2.40005091 -1.44032645 -0.08037966
O          2.23965160 -2.64008185 -0.25503382
N          2.41124583 -0.56083789 -1.16848449
H          2.23555694 -0.97799874 -2.07261756
C          2.49618822  0.81416973 -1.15279815
O          2.41351976  1.51385406 -2.14916326
H          2.71522963  2.32251674  0.17911856
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '15-0.9')] = qcdb.Molecule("""
0 1
N    0.067390759   1.213806097  -1.171192513
C   -0.034440687   0.160916029  -2.035179690
H   -0.037909102   0.307694674  -3.102311444
N   -0.122286497  -1.014214485  -1.431659388
C   -0.061278153  -0.690156063  -0.097738525
C   -0.083866474  -1.480006435   1.065121981
N   -0.207551291  -2.830167865   1.008466281
H    0.020236002  -3.318294510   1.858492777
H    0.100823981  -3.261839820   0.151791829
N   -0.015107287  -0.872886238   2.254820437
C    0.095534438   0.468473589   2.286592142
H    0.148443656   0.902433537   3.277055537
N    0.150791629   1.330817541   1.268232413
C    0.061278153   0.690156063   0.097738525
H    0.213123816   2.178532043  -1.420082564
--
0 1
N    2.995457244   1.318912569   0.115169333
C    3.033773997   0.544134785   1.248235461
H    3.166936649   1.084216460   2.174491246
C    2.913123372  -0.802036026   1.213306349
C    2.965573998  -1.664227788   2.429380731
H    2.009790775  -2.161867438   2.585037720
H    3.726416066  -2.435033978   2.315487569
H    3.189128467  -1.070628980   3.313538183
C    2.718644614  -1.440326451  -0.080379664
O    2.558245305  -2.640081851  -0.255033817
N    2.729839539  -0.560837886  -1.168484485
H    2.554150647  -0.977998743  -2.072617562
C    2.814781928   0.814169728  -1.152798148
O    2.732113465   1.513854058  -2.149163262
H    3.033823338   2.322516737   0.179118562
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '15-1.0')] = qcdb.Molecule("""
0 1
N    0.279301000   2.406839000  -0.605752000
C   -1.084857000   2.445746000  -0.551161000
H   -1.659440000   3.023029000  -1.256090000
N   -1.597712000   1.717988000   0.428754000
C   -0.489725000   1.171436000   1.030191000
C   -0.346137000   0.291471000   2.117234000
N   -1.418709000  -0.167777000   2.810144000
H   -1.238875000  -0.959480000   3.404758000
H   -2.291873000  -0.178822000   2.307362000
N    0.885763000  -0.070076000   2.491949000
C    1.935235000   0.407288000   1.796802000
H    2.906033000   0.078841000   2.145818000
N    1.940978000   1.224202000   0.740220000
C    0.695219000   1.577986000   0.406398000
H    0.861007000   2.829804000  -1.310450000
--
0 1
N    1.275461000  -0.647899000  -1.977910000
C    1.413053000  -1.553685000  -0.955067000
H    2.425877000  -1.867078000  -0.746878000
C    0.357598000  -2.023950000  -0.253057000
C    0.482129000  -3.017949000   0.852122000
H    0.175770000  -2.575607000   1.798628000
H   -0.160169000  -3.877041000   0.663950000
H    1.511244000  -3.357277000   0.951366000
C   -0.968471000  -1.529811000  -0.593979000
O   -2.002928000  -1.839696000  -0.019945000
N   -0.995692000  -0.638387000  -1.672042000
H   -1.901406000  -0.250172000  -1.898576000
C    0.068470000  -0.119176000  -2.376376000
O   -0.039788000   0.722701000  -3.253108000
H    2.085329000  -0.276018000  -2.445458000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '15-1.2')] = qcdb.Molecule("""
0 1
N    0.067390759   1.213806097  -1.171192513
C   -0.034440687   0.160916029  -2.035179690
H   -0.037909102   0.307694674  -3.102311444
N   -0.122286497  -1.014214485  -1.431659388
C   -0.061278153  -0.690156063  -0.097738525
C   -0.083866474  -1.480006435   1.065121981
N   -0.207551291  -2.830167865   1.008466281
H    0.020236002  -3.318294510   1.858492777
H    0.100823981  -3.261839820   0.151791829
N   -0.015107287  -0.872886238   2.254820437
C    0.095534438   0.468473589   2.286592142
H    0.148443656   0.902433537   3.277055537
N    0.150791629   1.330817541   1.268232413
C    0.061278153   0.690156063   0.097738525
H    0.213123816   2.178532043  -1.420082564
--
0 1
N    3.951238365   1.318912569   0.115169333
C    3.989555118   0.544134785   1.248235461
H    4.122717770   1.084216460   2.174491246
C    3.868904493  -0.802036026   1.213306349
C    3.921355119  -1.664227788   2.429380731
H    2.965571896  -2.161867438   2.585037720
H    4.682197187  -2.435033978   2.315487569
H    4.144909588  -1.070628980   3.313538183
C    3.674425735  -1.440326451  -0.080379664
O    3.514026426  -2.640081851  -0.255033817
N    3.685620660  -0.560837886  -1.168484485
H    3.509931768  -0.977998743  -2.072617562
C    3.770563049   0.814169728  -1.152798148
O    3.687894586   1.513854058  -2.149163262
H    3.989604459   2.322516737   0.179118562
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '15-1.5')] = qcdb.Molecule("""
0 1
N    0.067390759   1.213806097  -1.171192513
C   -0.034440687   0.160916029  -2.035179690
H   -0.037909102   0.307694674  -3.102311444
N   -0.122286497  -1.014214485  -1.431659388
C   -0.061278153  -0.690156063  -0.097738525
C   -0.083866474  -1.480006435   1.065121981
N   -0.207551291  -2.830167865   1.008466281
H    0.020236002  -3.318294510   1.858492777
H    0.100823981  -3.261839820   0.151791829
N   -0.015107287  -0.872886238   2.254820437
C    0.095534438   0.468473589   2.286592142
H    0.148443656   0.902433537   3.277055537
N    0.150791629   1.330817541   1.268232413
C    0.061278153   0.690156063   0.097738525
H    0.213123816   2.178532043  -1.420082564
--
0 1
N    4.907019487   1.318912569   0.115169333
C    4.945336240   0.544134785   1.248235461
H    5.078498892   1.084216460   2.174491246
C    4.824685615  -0.802036026   1.213306349
C    4.877136241  -1.664227788   2.429380731
H    3.921353018  -2.161867438   2.585037720
H    5.637978309  -2.435033978   2.315487569
H    5.100690710  -1.070628980   3.313538183
C    4.630206857  -1.440326451  -0.080379664
O    4.469807548  -2.640081851  -0.255033817
N    4.641401782  -0.560837886  -1.168484485
H    4.465712890  -0.977998743  -2.072617562
C    4.726344171   0.814169728  -1.152798148
O    4.643675708   1.513854058  -2.149163262
H    4.945385581   2.322516737   0.179118562
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '15-2.0')] = qcdb.Molecule("""
0 1
N    0.067390759   1.213806097  -1.171192513
C   -0.034440687   0.160916029  -2.035179690
H   -0.037909102   0.307694674  -3.102311444
N   -0.122286497  -1.014214485  -1.431659388
C   -0.061278153  -0.690156063  -0.097738525
C   -0.083866474  -1.480006435   1.065121981
N   -0.207551291  -2.830167865   1.008466281
H    0.020236002  -3.318294510   1.858492777
H    0.100823981  -3.261839820   0.151791829
N   -0.015107287  -0.872886238   2.254820437
C    0.095534438   0.468473589   2.286592142
H    0.148443656   0.902433537   3.277055537
N    0.150791629   1.330817541   1.268232413
C    0.061278153   0.690156063   0.097738525
H    0.213123816   2.178532043  -1.420082564
--
0 1
N    6.499988023   1.318912569   0.115169333
C    6.538304776   0.544134785   1.248235461
H    6.671467428   1.084216460   2.174491246
C    6.417654151  -0.802036026   1.213306349
C    6.470104777  -1.664227788   2.429380731
H    5.514321554  -2.161867438   2.585037720
H    7.230946845  -2.435033978   2.315487569
H    6.693659246  -1.070628980   3.313538183
C    6.223175393  -1.440326451  -0.080379664
O    6.062776084  -2.640081851  -0.255033817
N    6.234370318  -0.560837886  -1.168484485
H    6.058681426  -0.977998743  -2.072617562
C    6.319312707   0.814169728  -1.152798148
O    6.236644244   1.513854058  -2.149163262
H    6.538354117   2.322516737   0.179118562
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '16-0.7')] = qcdb.Molecule("""
0 1
C   0.00000000e+00  -6.67578000e-01   0.00000000e+00
C   0.00000000e+00   6.67578000e-01   0.00000000e+00
H  -1.52600000e-03  -1.23225300e+00  -9.23621000e-01
H  -1.52600000e-03  -1.23225300e+00   9.23621000e-01
H  -1.52600000e-03   1.23225300e+00   9.23621000e-01
H  -1.52600000e-03   1.23225300e+00  -9.23621000e-01
--
0 1
C   4.19955870e+00   0.00000000e+00   0.00000000e+00
C   2.99229570e+00   0.00000000e+00   0.00000000e+00
H   1.92640770e+00   0.00000000e+00   0.00000000e+00
H   5.26298470e+00   0.00000000e+00   0.00000000e+00
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '16-0.8')] = qcdb.Molecule("""
0 1
C   0.00000000e+00  -6.67578000e-01   0.00000000e+00
C   0.00000000e+00   6.67578000e-01   0.00000000e+00
H  -1.52600000e-03  -1.23225300e+00  -9.23621000e-01
H  -1.52600000e-03  -1.23225300e+00   9.23621000e-01
H  -1.52600000e-03   1.23225300e+00   9.23621000e-01
H  -1.52600000e-03   1.23225300e+00  -9.23621000e-01
--
0 1
C   4.47475980e+00   0.00000000e+00   0.00000000e+00
C   3.26749680e+00   0.00000000e+00   0.00000000e+00
H   2.20160880e+00   0.00000000e+00   0.00000000e+00
H   5.53818580e+00   0.00000000e+00   0.00000000e+00
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '16-0.9')] = qcdb.Molecule("""
0 1
C    0.000000000  -0.667578000   0.000000000
C    0.000000000   0.667578000   0.000000000
H   -0.001526000  -1.232253000  -0.923621000
H   -0.001526000  -1.232253000   0.923621000
H   -0.001526000   1.232253000   0.923621000
H   -0.001526000   1.232253000  -0.923621000
--
0 1
C    4.749960900   0.000000000   0.000000000
C    3.542697900   0.000000000   0.000000000
H    2.476809900   0.000000000   0.000000000
H    5.813386900   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '16-1.0')] = qcdb.Molecule("""
0 1
C    0.000000000  -0.667578000  -2.124659000
C    0.000000000   0.667578000  -2.124659000
H    0.923621000  -1.232253000  -2.126185000
H   -0.923621000  -1.232253000  -2.126185000
H   -0.923621000   1.232253000  -2.126185000
H    0.923621000   1.232253000  -2.126185000
--
0 1
C    0.000000000   0.000000000   2.900503000
C    0.000000000   0.000000000   1.693240000
H    0.000000000   0.000000000   0.627352000
H    0.000000000   0.000000000   3.963929000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '16-1.2')] = qcdb.Molecule("""
0 1
C    0.000000000  -0.667578000   0.000000000
C    0.000000000   0.667578000   0.000000000
H   -0.001526000  -1.232253000  -0.923621000
H   -0.001526000  -1.232253000   0.923621000
H   -0.001526000   1.232253000   0.923621000
H   -0.001526000   1.232253000  -0.923621000
--
0 1
C    5.575564200   0.000000000   0.000000000
C    4.368301200   0.000000000   0.000000000
H    3.302413200   0.000000000   0.000000000
H    6.638990200   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '16-1.5')] = qcdb.Molecule("""
0 1
C    0.000000000  -0.667578000   0.000000000
C    0.000000000   0.667578000   0.000000000
H   -0.001526000  -1.232253000  -0.923621000
H   -0.001526000  -1.232253000   0.923621000
H   -0.001526000   1.232253000   0.923621000
H   -0.001526000   1.232253000  -0.923621000
--
0 1
C    6.401167500   0.000000000   0.000000000
C    5.193904500   0.000000000   0.000000000
H    4.128016500   0.000000000   0.000000000
H    7.464593500   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '16-2.0')] = qcdb.Molecule("""
0 1
C    0.000000000  -0.667578000   0.000000000
C    0.000000000   0.667578000   0.000000000
H   -0.001526000  -1.232253000  -0.923621000
H   -0.001526000  -1.232253000   0.923621000
H   -0.001526000   1.232253000   0.923621000
H   -0.001526000   1.232253000  -0.923621000
--
0 1
C    7.777173000   0.000000000   0.000000000
C    6.569910000   0.000000000   0.000000000
H    5.504022000   0.000000000   0.000000000
H    8.840599000   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '17-0.7')] = qcdb.Molecule("""
0 1
C          0.06873616  1.39238384 -1.20754300
C          0.00000000  0.00000000 -1.20790400
C         -0.03480730 -0.69643588  0.00000000
C          0.00000000  0.00000000  1.20790400
C          0.06873616  1.39238384  1.20754300
C          0.10258114  2.08831334  0.00000000
H          0.09647711  1.93199935 -2.14414800
H         -0.02281541 -0.54039795 -2.14405500
H         -0.08669494 -1.77649774  0.00000000
H         -0.02281541 -0.54039795  2.14405500
H          0.09647711  1.93199935  2.14414800
H          0.15343075  3.16857919  0.00000000
--
0 1
O          2.68148004  0.12436973  0.00000000
H          2.77175628  1.07911799  0.00000000
H          1.72753553  0.00000000  0.00000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '17-0.8')] = qcdb.Molecule("""
0 1
C          0.06873616  1.39238384 -1.20754300
C          0.00000000  0.00000000 -1.20790400
C         -0.03480730 -0.69643588  0.00000000
C          0.00000000  0.00000000  1.20790400
C          0.06873616  1.39238384  1.20754300
C          0.10258114  2.08831334  0.00000000
H          0.09647711  1.93199935 -2.14414800
H         -0.02281541 -0.54039795 -2.14405500
H         -0.08669494 -1.77649774  0.00000000
H         -0.02281541 -0.54039795  2.14405500
H          0.09647711  1.93199935  2.14414800
H          0.15343075  3.16857919  0.00000000
--
0 1
O          2.92827083  0.12436973  0.00000000
H          3.01854707  1.07911799  0.00000000
H          1.97432633  0.00000000  0.00000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '17-0.9')] = qcdb.Molecule("""
0 1
C    0.068736158   1.392383840  -1.207543000
C    0.000000000   0.000000000  -1.207904000
C   -0.034807303  -0.696435878   0.000000000
C    0.000000000   0.000000000   1.207904000
C    0.068736158   1.392383840   1.207543000
C    0.102581137   2.088313342   0.000000000
H    0.096477114   1.931999350  -2.144148000
H   -0.022815407  -0.540397951  -2.144055000
H   -0.086694943  -1.776497744   0.000000000
H   -0.022815407  -0.540397951   2.144055000
H    0.096477114   1.931999350   2.144148000
H    0.153430751   3.168579194   0.000000000
--
0 1
O    3.175061618   0.124369730   0.000000000
H    3.265337861   1.079117991   0.000000000
H    2.221117117   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '17-1.0')] = qcdb.Molecule("""
0 1
C    0.780612000  -0.609888000  -1.207543000
C    0.478404000   0.751041000  -1.207904000
C    0.327659000   1.431857000   0.000000000
C    0.478404000   0.751041000   1.207904000
C    0.780612000  -0.609888000   1.207543000
C    0.932151000  -1.289961000   0.000000000
H    0.896669000  -1.137605000  -2.144148000
H    0.357390000   1.278209000  -2.144055000
H    0.091859000   2.487141000   0.000000000
H    0.357390000   1.278209000   2.144055000
H    0.896669000  -1.137605000   2.144148000
H    1.169006000  -2.345167000   0.000000000
--
0 1
O   -2.788527000  -0.274485000   0.000000000
H   -2.622911000  -1.219083000   0.000000000
H   -1.901510000   0.097911000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '17-1.2')] = qcdb.Molecule("""
0 1
C    0.068736158   1.392383840  -1.207543000
C    0.000000000   0.000000000  -1.207904000
C   -0.034807303  -0.696435878   0.000000000
C    0.000000000   0.000000000   1.207904000
C    0.068736158   1.392383840   1.207543000
C    0.102581137   2.088313342   0.000000000
H    0.096477114   1.931999350  -2.144148000
H   -0.022815407  -0.540397951  -2.144055000
H   -0.086694943  -1.776497744   0.000000000
H   -0.022815407  -0.540397951   2.144055000
H    0.096477114   1.931999350   2.144148000
H    0.153430751   3.168579194   0.000000000
--
0 1
O    3.915433991   0.124369730   0.000000000
H    4.005710234   1.079117991   0.000000000
H    2.961489490   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '17-1.5')] = qcdb.Molecule("""
0 1
C    0.068736158   1.392383840  -1.207543000
C    0.000000000   0.000000000  -1.207904000
C   -0.034807303  -0.696435878   0.000000000
C    0.000000000   0.000000000   1.207904000
C    0.068736158   1.392383840   1.207543000
C    0.102581137   2.088313342   0.000000000
H    0.096477114   1.931999350  -2.144148000
H   -0.022815407  -0.540397951  -2.144055000
H   -0.086694943  -1.776497744   0.000000000
H   -0.022815407  -0.540397951   2.144055000
H    0.096477114   1.931999350   2.144148000
H    0.153430751   3.168579194   0.000000000
--
0 1
O    4.655806363   0.124369730   0.000000000
H    4.746082606   1.079117991   0.000000000
H    3.701861862   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '17-2.0')] = qcdb.Molecule("""
0 1
C    0.068736158   1.392383840  -1.207543000
C    0.000000000   0.000000000  -1.207904000
C   -0.034807303  -0.696435878   0.000000000
C    0.000000000   0.000000000   1.207904000
C    0.068736158   1.392383840   1.207543000
C    0.102581137   2.088313342   0.000000000
H    0.096477114   1.931999350  -2.144148000
H   -0.022815407  -0.540397951  -2.144055000
H   -0.086694943  -1.776497744   0.000000000
H   -0.022815407  -0.540397951   2.144055000
H    0.096477114   1.931999350   2.144148000
H    0.153430751   3.168579194   0.000000000
--
0 1
O    5.889760317   0.124369730   0.000000000
H    5.980036560   1.079117991   0.000000000
H    4.935815816   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '18-0.7')] = qcdb.Molecule("""
0 1
C          0.00000000  0.00000000 -1.20710800
C         -0.09472391 -0.69068717  0.00000000
C          0.00000000  0.00000000  1.20710800
C          0.18929305  1.38119484  1.20707300
C          0.28420947  2.07177137  0.00000000
C          0.18929305  1.38119484 -1.20707300
H         -0.07088443 -0.53645471 -2.14328900
H         -0.23533516 -1.76264080  0.00000000
H         -0.07088443 -0.53645471  2.14328900
H          0.26243423  1.91683009  2.14369500
H          0.43037381  3.14325787  0.00000000
H          0.26243423  1.91683009 -2.14369500
--
0 1
N          2.80591307 -0.17515845  0.00000000
H          3.16920386  0.31696099 -0.80607300
H          3.16920386  0.31696099  0.80607300
H          1.80781864  0.00000000  0.00000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '18-0.8')] = qcdb.Molecule("""
0 1
C          0.00000000  0.00000000 -1.20710800
C         -0.09472391 -0.69068717  0.00000000
C          0.00000000  0.00000000  1.20710800
C          0.18929305  1.38119484  1.20707300
C          0.28420947  2.07177137  0.00000000
C          0.18929305  1.38119484 -1.20707300
H         -0.07088443 -0.53645471 -2.14328900
H         -0.23533516 -1.7626408   0.00000000
H         -0.07088443 -0.53645471  2.14328900
H          0.26243423  1.91683009  2.14369500
H          0.43037381  3.14325787  0.00000000
H          0.26243423  1.91683009 -2.14369500
--
0 1
N          3.06417287 -0.17515845  0.00000000
H          3.42746366  0.31696099 -0.80607300
H          3.42746366  0.31696099  0.80607300
H          2.06607844  0.00000000  0.00000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '18-0.9')] = qcdb.Molecule("""
0 1
C    0.000000000   0.000000000  -1.207108000
C   -0.094723910  -0.690687169   0.000000000
C    0.000000000   0.000000000   1.207108000
C    0.189293052   1.381194838   1.207073000
C    0.284209467   2.071771374   0.000000000
C    0.189293052   1.381194838  -1.207073000
H   -0.070884435  -0.536454706  -2.143289000
H   -0.235335157  -1.762640796   0.000000000
H   -0.070884435  -0.536454706   2.143289000
H    0.262434233   1.916830087   2.143695000
H    0.430373810   3.143257869   0.000000000
H    0.262434233   1.916830087  -2.143695000
--
0 1
N    3.322432676  -0.175158455   0.000000000
H    3.685723470   0.316960994  -0.806073000
H    3.685723470   0.316960994   0.806073000
H    2.324338249   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '18-1.0')] = qcdb.Molecule("""
0 1
C   -0.739281000   0.515879000  -1.207108000
C   -1.426144000   0.396545000   0.000000000
C   -0.739281000   0.515879000   1.207108000
C    0.634227000   0.754640000   1.207073000
C    1.321043000   0.873757000   0.000000000
C    0.634227000   0.754640000  -1.207073000
H   -1.271950000   0.420632000  -2.143289000
H   -2.490220000   0.205238000   0.000000000
H   -1.271950000   0.420632000   2.143289000
H    1.166800000   0.847488000   2.143695000
H    2.386359000   1.059631000   0.000000000
H    1.166800000   0.847488000  -2.143695000
--
0 1
N    0.180393000  -2.949123000   0.000000000
H    0.759549000  -3.145948000  -0.806073000
H    0.759549000  -3.145948000   0.806073000
H    0.044417000  -1.944940000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '18-1.2')] = qcdb.Molecule("""
0 1
C    0.000000000   0.000000000  -1.207108000
C   -0.094723910  -0.690687169   0.000000000
C    0.000000000   0.000000000   1.207108000
C    0.189293052   1.381194838   1.207073000
C    0.284209467   2.071771374   0.000000000
C    0.189293052   1.381194838  -1.207073000
H   -0.070884435  -0.536454706  -2.143289000
H   -0.235335157  -1.762640796   0.000000000
H   -0.070884435  -0.536454706   2.143289000
H    0.262434233   1.916830087   2.143695000
H    0.430373810   3.143257869   0.000000000
H    0.262434233   1.916830087  -2.143695000
--
0 1
N    4.097212092  -0.175158455   0.000000000
H    4.460502886   0.316960994  -0.806073000
H    4.460502886   0.316960994   0.806073000
H    3.099117665   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '18-1.5')] = qcdb.Molecule("""
0 1
C    0.000000000   0.000000000  -1.207108000
C   -0.094723910  -0.690687169   0.000000000
C    0.000000000   0.000000000   1.207108000
C    0.189293052   1.381194838   1.207073000
C    0.284209467   2.071771374   0.000000000
C    0.189293052   1.381194838  -1.207073000
H   -0.070884435  -0.536454706  -2.143289000
H   -0.235335157  -1.762640796   0.000000000
H   -0.070884435  -0.536454706   2.143289000
H    0.262434233   1.916830087   2.143695000
H    0.430373810   3.143257869   0.000000000
H    0.262434233   1.916830087  -2.143695000
--
0 1
N    4.871991508  -0.175158455   0.000000000
H    5.235282302   0.316960994  -0.806073000
H    5.235282302   0.316960994   0.806073000
H    3.873897081   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '18-2.0')] = qcdb.Molecule("""
0 1
C    0.000000000   0.000000000  -1.207108000
C   -0.094723910  -0.690687169   0.000000000
C    0.000000000   0.000000000   1.207108000
C    0.189293052   1.381194838   1.207073000
C    0.284209467   2.071771374   0.000000000
C    0.189293052   1.381194838  -1.207073000
H   -0.070884435  -0.536454706  -2.143289000
H   -0.235335157  -1.762640796   0.000000000
H   -0.070884435  -0.536454706   2.143289000
H    0.262434233   1.916830087   2.143695000
H    0.430373810   3.143257869   0.000000000
H    0.262434233   1.916830087  -2.143695000
--
0 1
N    6.163290535  -0.175158455   0.000000000
H    6.526581329   0.316960994  -0.806073000
H    6.526581329   0.316960994   0.806073000
H    5.165196108   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '19-0.7')] = qcdb.Molecule("""
0 1
C         -0.02310095  0.69697859  1.20770200
C         -0.04616033  1.39380803  0.00000000
C         -0.02310095  0.69697859 -1.20770200
C          0.02308582 -0.69689511 -1.20786500
C          0.04619059 -1.39397501  0.00000000
C          0.02308582 -0.69689511  1.20786500
H         -0.03862462  1.23736918  2.14405100
H         -0.07914868  2.47449307  0.00000000
H         -0.03862462  1.23736918 -2.14405100
H          0.04283969 -1.23714251 -2.14425600
H          0.08340142 -2.47459358  0.00000000
H          0.04283969 -1.23714251  2.14425600
--
0 1
N          3.84277676  0.30453686  0.00000000
C          2.68628601  0.14576395  0.00000000
H          1.62840272  0.00000000  0.00000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '19-0.8')] = qcdb.Molecule("""
0 1
C         -0.02310095  0.69697859  1.20770200
C         -0.04616033  1.39380803  0.00000000
C         -0.02310095  0.69697859 -1.20770200
C          0.02308582 -0.69689511 -1.20786500
C          0.04619059 -1.39397501  0.00000000
C          0.02308582 -0.69689511  1.20786500
H         -0.03862462  1.23736918  2.14405100
H         -0.07914868  2.47449307  0.00000000
H         -0.03862462  1.23736918 -2.14405100
H          0.04283969 -1.23714251 -2.14425600
H          0.08340142 -2.47459358  0.00000000
H          0.04283969 -1.23714251  2.14425600
--
0 1
N          4.07540572  0.30453686  0.00000000
C          2.91891497  0.14576395  0.00000000
H          1.86103168  0.00000000  0.00000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '19-0.9')] = qcdb.Molecule("""
0 1
C   -0.023100946   0.696978594   1.207702000
C   -0.046160335   1.393808033   0.000000000
C   -0.023100946   0.696978594  -1.207702000
C    0.023085816  -0.696895106  -1.207865000
C    0.046190594  -1.393975010   0.000000000
C    0.023085816  -0.696895106   1.207865000
H   -0.038624622   1.237369182   2.144051000
H   -0.079148681   2.474493071   0.000000000
H   -0.038624622   1.237369182  -2.144051000
H    0.042839694  -1.237142510  -2.144256000
H    0.083401415  -2.474593580   0.000000000
H    0.042839694  -1.237142510   2.144256000
--
0 1
N    4.308034683   0.304536859   0.000000000
C    3.151543935   0.145763954   0.000000000
H    2.093660645   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '19-1.0')] = qcdb.Molecule("""
0 1
C   -0.709774000  -0.990423000   1.207702000
C   -1.406534000  -0.965353000   0.000000000
C   -0.709774000  -0.990423000  -1.207702000
C    0.683965000  -1.040510000  -1.207865000
C    1.380978000  -1.065552000   0.000000000
C    0.683965000  -1.040510000   1.207865000
H   -1.249948000  -0.968628000   2.144051000
H   -2.486920000  -0.923706000   0.000000000
H   -1.249948000  -0.968628000  -2.144051000
H    1.224288000  -1.058075000  -2.144256000
H    2.461589000  -1.102982000   0.000000000
H    1.224288000  -1.058075000   2.144256000
--
0 1
N   -0.003412000   3.535393000   0.000000000
C    0.075196000   2.370704000   0.000000000
H    0.147629000   1.305285000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '19-1.2')] = qcdb.Molecule("""
0 1
C   -0.023100946   0.696978594   1.207702000
C   -0.046160335   1.393808033   0.000000000
C   -0.023100946   0.696978594  -1.207702000
C    0.023085816  -0.696895106  -1.207865000
C    0.046190594  -1.393975010   0.000000000
C    0.023085816  -0.696895106   1.207865000
H   -0.038624622   1.237369182   2.144051000
H   -0.079148681   2.474493071   0.000000000
H   -0.038624622   1.237369182  -2.144051000
H    0.042839694  -1.237142510  -2.144256000
H    0.083401415  -2.474593580   0.000000000
H    0.042839694  -1.237142510   2.144256000
--
0 1
N    5.005921565   0.304536859   0.000000000
C    3.849430817   0.145763954   0.000000000
H    2.791547527   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '19-1.5')] = qcdb.Molecule("""
0 1
C   -0.023100946   0.696978594   1.207702000
C   -0.046160335   1.393808033   0.000000000
C   -0.023100946   0.696978594  -1.207702000
C    0.023085816  -0.696895106  -1.207865000
C    0.046190594  -1.393975010   0.000000000
C    0.023085816  -0.696895106   1.207865000
H   -0.038624622   1.237369182   2.144051000
H   -0.079148681   2.474493071   0.000000000
H   -0.038624622   1.237369182  -2.144051000
H    0.042839694  -1.237142510  -2.144256000
H    0.083401415  -2.474593580   0.000000000
H    0.042839694  -1.237142510   2.144256000
--
0 1
N    5.703808447   0.304536859   0.000000000
C    4.547317699   0.145763954   0.000000000
H    3.489434409   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '19-2.0')] = qcdb.Molecule("""
0 1
C   -0.023100946   0.696978594   1.207702000
C   -0.046160335   1.393808033   0.000000000
C   -0.023100946   0.696978594  -1.207702000
C    0.023085816  -0.696895106  -1.207865000
C    0.046190594  -1.393975010   0.000000000
C    0.023085816  -0.696895106   1.207865000
H   -0.038624622   1.237369182   2.144051000
H   -0.079148681   2.474493071   0.000000000
H   -0.038624622   1.237369182  -2.144051000
H    0.042839694  -1.237142510  -2.144256000
H    0.083401415  -2.474593580   0.000000000
H    0.042839694  -1.237142510   2.144256000
--
0 1
N    6.866953250   0.304536859   0.000000000
C    5.710462502   0.145763954   0.000000000
H    4.652579212   0.000000000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '20-0.7')] = qcdb.Molecule("""
0 1
C         -1.080615    0.000000    0.000000
C         -1.779254   -1.206008    0.000000
C         -3.173171   -1.207177    0.000000
C         -3.870155    0.000000    0.000000
C         -3.173171    1.207177    0.000000
C         -1.779254    1.206008    0.000000
H          0.000000    0.000000    0.000000
H         -1.236002   -2.141639    0.000000
H         -3.714575   -2.143566    0.000000
H         -4.951730    0.000000    0.000000
H         -3.714575    2.143566    0.000000
H         -1.236002    2.141639    0.000000
--
0 1
C          1.7027052   0.000000   -1.394063
C          1.7031812   1.207238   -0.697047
C          1.7031812   1.207238    0.697047
C          1.7027052   0.000000    1.394063
C          1.7031812  -1.207238    0.697047
C          1.7031812  -1.207238   -0.697047
H          1.6988752   0.000000   -2.475399
H          1.7022292   2.143565   -1.238232
H          1.7022292   2.143565    1.238232
H          1.6988752   0.000000    2.475399
H          1.7022292  -2.143565    1.238232
H          1.7022292  -2.143565   -1.238232
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '20-0.8')] = qcdb.Molecule("""
0 1
C         -1.080615    0.000000    0.000000
C         -1.779254   -1.206008    0.000000
C         -3.173171   -1.207177    0.000000
C         -3.870155    0.000000    0.000000
C         -3.173171    1.207177    0.000000
C         -1.779254    1.206008    0.000000
H          0.000000    0.000000    0.000000
H         -1.236002   -2.141639    0.000000
H         -3.714575   -2.143566    0.000000
H         -4.951730    0.000000    0.000000
H         -3.714575    2.143566    0.000000
H         -1.236002    2.141639    0.000000
--
0 1
C          1.94599413  0.000000   -1.394063
C          1.94647013  1.207238   -0.697047
C          1.94647013  1.207238    0.697047
C          1.94599413  0.000000    1.394063
C          1.94647013 -1.207238    0.697047
C          1.94647013 -1.207238   -0.697047
H          1.94216413  0.000000   -2.475399
H          1.94551813  2.143565   -1.238232
H          1.94551813  2.143565    1.238232
H          1.94216413  0.000000    2.475399
H          1.94551813 -2.143565    1.238232
H          1.94551813 -2.143565   -1.238232
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '20-0.9')] = qcdb.Molecule("""
0 1
C   -1.080615000   0.000000000   0.000000000
C   -1.779254000  -1.206008000   0.000000000
C   -3.173171000  -1.207177000   0.000000000
C   -3.870155000   0.000000000   0.000000000
C   -3.173171000   1.207177000   0.000000000
C   -1.779254000   1.206008000   0.000000000
H    0.000000000   0.000000000   0.000000000
H   -1.236002000  -2.141639000   0.000000000
H   -3.714575000  -2.143566000   0.000000000
H   -4.951730000   0.000000000   0.000000000
H   -3.714575000   2.143566000   0.000000000
H   -1.236002000   2.141639000   0.000000000
--
0 1
C    2.189283067   0.000000000  -1.394063000
C    2.189759067   1.207238000  -0.697047000
C    2.189759067   1.207238000   0.697047000
C    2.189283067   0.000000000   1.394063000
C    2.189759067  -1.207238000   0.697047000
C    2.189759067  -1.207238000  -0.697047000
H    2.185453067   0.000000000  -2.475399000
H    2.188807067   2.143565000  -1.238232000
H    2.188807067   2.143565000   1.238232000
H    2.185453067   0.000000000   2.475399000
H    2.188807067  -2.143565000   1.238232000
H    2.188807067  -2.143565000  -1.238232000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '20-1.0')] = qcdb.Molecule("""
0 1
C    0.000000000   0.000000000   1.059035000
C    0.000000000  -1.206008000   1.757674000
C    0.000000000  -1.207177000   3.151591000
C    0.000000000   0.000000000   3.848575000
C    0.000000000   1.207177000   3.151591000
C    0.000000000   1.206008000   1.757674000
H    0.000000000   0.000000000  -0.021580000
H    0.000000000  -2.141639000   1.214422000
H    0.000000000  -2.143566000   3.692995000
H    0.000000000   0.000000000   4.930150000
H    0.000000000   2.143566000   3.692995000
H    0.000000000   2.141639000   1.214422000
--
0 1
C   -1.394063000   0.000000000  -2.454152000
C   -0.697047000   1.207238000  -2.454628000
C    0.697047000   1.207238000  -2.454628000
C    1.394063000   0.000000000  -2.454152000
C    0.697047000  -1.207238000  -2.454628000
C   -0.697047000  -1.207238000  -2.454628000
H   -2.475399000   0.000000000  -2.450322000
H   -1.238232000   2.143565000  -2.453676000
H    1.238232000   2.143565000  -2.453676000
H    2.475399000   0.000000000  -2.450322000
H    1.238232000  -2.143565000  -2.453676000
H   -1.238232000  -2.143565000  -2.453676000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '20-1.2')] = qcdb.Molecule("""
0 1
C   -1.080615000   0.000000000   0.000000000
C   -1.779254000  -1.206008000   0.000000000
C   -3.173171000  -1.207177000   0.000000000
C   -3.870155000   0.000000000   0.000000000
C   -3.173171000   1.207177000   0.000000000
C   -1.779254000   1.206008000   0.000000000
H    0.000000000   0.000000000   0.000000000
H   -1.236002000  -2.141639000   0.000000000
H   -3.714575000  -2.143566000   0.000000000
H   -4.951730000   0.000000000   0.000000000
H   -3.714575000   2.143566000   0.000000000
H   -1.236002000   2.141639000   0.000000000
--
0 1
C    2.919149867   0.000000000  -1.394063000
C    2.919625867   1.207238000  -0.697047000
C    2.919625867   1.207238000   0.697047000
C    2.919149867   0.000000000   1.394063000
C    2.919625867  -1.207238000   0.697047000
C    2.919625867  -1.207238000  -0.697047000
H    2.915319867   0.000000000  -2.475399000
H    2.918673867   2.143565000  -1.238232000
H    2.918673867   2.143565000   1.238232000
H    2.915319867   0.000000000   2.475399000
H    2.918673867  -2.143565000   1.238232000
H    2.918673867  -2.143565000  -1.238232000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '20-1.5')] = qcdb.Molecule("""
0 1
C   -1.080615000   0.000000000   0.000000000
C   -1.779254000  -1.206008000   0.000000000
C   -3.173171000  -1.207177000   0.000000000
C   -3.870155000   0.000000000   0.000000000
C   -3.173171000   1.207177000   0.000000000
C   -1.779254000   1.206008000   0.000000000
H    0.000000000   0.000000000   0.000000000
H   -1.236002000  -2.141639000   0.000000000
H   -3.714575000  -2.143566000   0.000000000
H   -4.951730000   0.000000000   0.000000000
H   -3.714575000   2.143566000   0.000000000
H   -1.236002000   2.141639000   0.000000000
--
0 1
C    3.649016667   0.000000000  -1.394063000
C    3.649492667   1.207238000  -0.697047000
C    3.649492667   1.207238000   0.697047000
C    3.649016667   0.000000000   1.394063000
C    3.649492667  -1.207238000   0.697047000
C    3.649492667  -1.207238000  -0.697047000
H    3.645186667   0.000000000  -2.475399000
H    3.648540667   2.143565000  -1.238232000
H    3.648540667   2.143565000   1.238232000
H    3.645186667   0.000000000   2.475399000
H    3.648540667  -2.143565000   1.238232000
H    3.648540667  -2.143565000  -1.238232000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '20-2.0')] = qcdb.Molecule("""
0 1
C   -1.080615000   0.000000000   0.000000000
C   -1.779254000  -1.206008000   0.000000000
C   -3.173171000  -1.207177000   0.000000000
C   -3.870155000   0.000000000   0.000000000
C   -3.173171000   1.207177000   0.000000000
C   -1.779254000   1.206008000   0.000000000
H    0.000000000   0.000000000   0.000000000
H   -1.236002000  -2.141639000   0.000000000
H   -3.714575000  -2.143566000   0.000000000
H   -4.951730000   0.000000000   0.000000000
H   -3.714575000   2.143566000   0.000000000
H   -1.236002000   2.141639000   0.000000000
--
0 1
C    4.865461333   0.000000000  -1.394063000
C    4.865937333   1.207238000  -0.697047000
C    4.865937333   1.207238000   0.697047000
C    4.865461333   0.000000000   1.394063000
C    4.865937333  -1.207238000   0.697047000
C    4.865937333  -1.207238000  -0.697047000
H    4.861631333   0.000000000  -2.475399000
H    4.864985333   2.143565000  -1.238232000
H    4.864985333   2.143565000   1.238232000
H    4.861631333   0.000000000   2.475399000
H    4.864985333  -2.143565000   1.238232000
H    4.864985333  -2.143565000  -1.238232000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '21-0.7')] = qcdb.Molecule("""
0 1
C  -5.26520770e-02  -1.39322578e+00   0.00000000
C  -2.55433470e-02  -6.96940104e-01  -1.20829200
C   2.63482540e-02   6.96724226e-01  -1.20836500
C   5.10422630e-02   1.39365754e+00   0.00000000
C   2.63482540e-02   6.96724226e-01   1.20836500
C  -2.55433470e-02  -6.96940104e-01   1.20829200
H  -9.74306610e-02  -2.47365597e+00   0.00000000
H  -4.05097560e-02  -1.23736007e+00  -2.14459000
H   5.09555750e-02   1.23653129e+00  -2.14483800
H   8.96576450e-02   2.47441242e+00   0.00000000
H   5.09555750e-02   1.23653129e+00   2.14483800
H  -4.05097560e-02  -1.23736007e+00   2.14459000
--
0 1
H   1.56162022e+00   0.00000000e+00   0.00000000
N   2.56893762e+00   5.05638800e-03   0.00000000
C   3.35059181e+00   1.13260494e+00   0.00000000
C   4.67947653e+00   7.72354616e-01   0.00000000
C   4.72087002e+00  -6.53193161e-01   0.00000000
C   3.37102538e+00  -1.10492088e+00   0.00000000
C   3.03636571e+00  -2.46209497e+00   0.00000000
C   4.07855802e+00  -3.37617889e+00   0.00000000
C   5.42288146e+00  -2.95164129e+00   0.00000000
C   5.75322134e+00  -1.60670557e+00   0.00000000
H   2.89689758e+00   2.10959476e+00   0.00000000
H   5.51486634e+00   1.45148992e+00   0.00000000
H   2.00397677e+00  -2.78573081e+00   0.00000000
H   3.85684057e+00  -4.43482278e+00   0.00000000
H   6.20894638e+00  -3.69457014e+00   0.00000000
H   6.78954712e+00  -1.29459388e+00   0.00000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '21-0.8')] = qcdb.Molecule("""
0 1
C  -5.26520770e-02  -1.39322578e+00   0.00000000
C  -2.55433470e-02  -6.96940104e-01  -1.20829200
C   2.63482540e-02   6.96724226e-01  -1.20836500
C   5.10422630e-02   1.39365754e+00   0.00000000
C   2.63482540e-02   6.96724226e-01   1.20836500
C  -2.55433470e-02  -6.96940104e-01   1.20829200
H  -9.74306610e-02  -2.47365597e+00   0.00000000
H  -4.05097560e-02  -1.23736007e+00  -2.14459000
H   5.09555750e-02   1.23653129e+00  -2.14483800
H   8.96576450e-02   2.47441242e+00   0.00000000
H   5.09555750e-02   1.23653129e+00   2.14483800
H  -4.05097560e-02  -1.23736007e+00   2.14459000
--
0 1
H   1.78470882e+00   0.00000000e+00   0.00000000
N   2.79202623e+00   5.05638800e-03   0.00000000
C   3.57368041e+00   1.13260494e+00   0.00000000
C   4.90256514e+00   7.72354616e-01   0.00000000
C   4.94395862e+00  -6.53193161e-01   0.00000000
C   3.59411399e+00  -1.10492088e+00   0.00000000
C   3.25945432e+00  -2.46209497e+00   0.00000000
C   4.30164662e+00  -3.37617889e+00   0.00000000
C   5.64597006e+00  -2.95164129e+00   0.00000000
C   5.97630994e+00  -1.60670557e+00   0.00000000
H   3.11998618e+00   2.10959476e+00   0.00000000
H   5.73795494e+00   1.45148992e+00   0.00000000
H   2.22706538e+00  -2.78573081e+00   0.00000000
H   4.07992918e+00  -4.43482278e+00   0.00000000
H   6.43203498e+00  -3.69457014e+00   0.00000000
H   7.01263572e+00  -1.29459388e+00   0.00000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '21-0.9')] = qcdb.Molecule("""
0 1
C   -0.052652077  -1.393225783   0.000000000
C   -0.025543347  -0.696940104  -1.208292000
C    0.026348254   0.696724226  -1.208365000
C    0.051042263   1.393657541   0.000000000
C    0.026348254   0.696724226   1.208365000
C   -0.025543347  -0.696940104   1.208292000
H   -0.097430661  -2.473655966   0.000000000
H   -0.040509756  -1.237360068  -2.144590000
H    0.050955575   1.236531293  -2.144838000
H    0.089657645   2.474412421   0.000000000
H    0.050955575   1.236531293   2.144838000
H   -0.040509756  -1.237360068   2.144590000
--
0 1
H    2.007797424   0.000000000   0.000000000
N    3.015114828   0.005056388   0.000000000
C    3.796769012   1.132604937   0.000000000
C    5.125653739   0.772354616   0.000000000
C    5.167047225  -0.653193161   0.000000000
C    3.817202589  -1.104920876   0.000000000
C    3.482542920  -2.462094972   0.000000000
C    4.524735226  -3.376178892   0.000000000
C    5.869058665  -2.951641292   0.000000000
C    6.199398544  -1.606705567   0.000000000
H    3.343074787   2.109594763   0.000000000
H    5.961043541   1.451489921   0.000000000
H    2.450153978  -2.785730808   0.000000000
H    4.303017780  -4.434822780   0.000000000
H    6.655123584  -3.694570139   0.000000000
H    7.235724321  -1.294593877   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '21-1.0')] = qcdb.Molecule("""
0 1
C    2.511900000   1.625015000   0.000000000
C    2.713009000   0.957854000  -1.208292000
C    3.117782000  -0.376744000  -1.208365000
C    3.321385000  -1.043731000   0.000000000
C    3.117782000  -0.376744000   1.208365000
C    2.713009000   0.957854000   1.208292000
H    2.202404000   2.661136000   0.000000000
H    2.551176000   1.473691000  -2.144590000
H    3.270300000  -0.895141000  -2.144838000
H    3.636814000  -2.078152000   0.000000000
H    3.270300000  -0.895141000   2.144838000
H    2.551176000   1.473691000   2.144590000
--
0 1
H    0.806524000  -0.435887000   0.000000000
N   -0.144241000  -0.768693000   0.000000000
C   -0.516112000  -2.089322000   0.000000000
C   -1.889876000  -2.181449000   0.000000000
C   -2.393232000  -0.847083000   0.000000000
C   -1.264065000   0.019589000   0.000000000
C   -1.389600000   1.411767000   0.000000000
C   -2.672650000   1.936645000   0.000000000
C   -3.805451000   1.097479000   0.000000000
C   -3.679817000  -0.281721000   0.000000000
H    0.231002000  -2.865317000   0.000000000
H   -2.458576000  -3.095605000   0.000000000
H   -0.518873000   2.053952000   0.000000000
H   -2.807757000   3.009786000   0.000000000
H   -4.790599000   1.543937000   0.000000000
H   -4.558019000  -0.914292000   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '21-1.2')] = qcdb.Molecule("""
0 1
C   -0.052652077  -1.393225783   0.000000000
C   -0.025543347  -0.696940104  -1.208292000
C    0.026348254   0.696724226  -1.208365000
C    0.051042263   1.393657541   0.000000000
C    0.026348254   0.696724226   1.208365000
C   -0.025543347  -0.696940104   1.208292000
H   -0.097430661  -2.473655966   0.000000000
H   -0.040509756  -1.237360068  -2.144590000
H    0.050955575   1.236531293  -2.144838000
H    0.089657645   2.474412421   0.000000000
H    0.050955575   1.236531293   2.144838000
H   -0.040509756  -1.237360068   2.144590000
--
0 1
H    2.677063232   0.000000000   0.000000000
N    3.684380636   0.005056388   0.000000000
C    4.466034820   1.132604937   0.000000000
C    5.794919547   0.772354616   0.000000000
C    5.836313033  -0.653193161   0.000000000
C    4.486468397  -1.104920876   0.000000000
C    4.151808728  -2.462094972   0.000000000
C    5.194001034  -3.376178892   0.000000000
C    6.538324473  -2.951641292   0.000000000
C    6.868664352  -1.606705567   0.000000000
H    4.012340595   2.109594763   0.000000000
H    6.630309349   1.451489921   0.000000000
H    3.119419786  -2.785730808   0.000000000
H    4.972283588  -4.434822780   0.000000000
H    7.324389392  -3.694570139   0.000000000
H    7.904990129  -1.294593877   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '21-1.5')] = qcdb.Molecule("""
0 1
C   -0.052652077  -1.393225783   0.000000000
C   -0.025543347  -0.696940104  -1.208292000
C    0.026348254   0.696724226  -1.208365000
C    0.051042263   1.393657541   0.000000000
C    0.026348254   0.696724226   1.208365000
C   -0.025543347  -0.696940104   1.208292000
H   -0.097430661  -2.473655966   0.000000000
H   -0.040509756  -1.237360068  -2.144590000
H    0.050955575   1.236531293  -2.144838000
H    0.089657645   2.474412421   0.000000000
H    0.050955575   1.236531293   2.144838000
H   -0.040509756  -1.237360068   2.144590000
--
0 1
H    3.346329040   0.000000000   0.000000000
N    4.353646444   0.005056388   0.000000000
C    5.135300628   1.132604937   0.000000000
C    6.464185355   0.772354616   0.000000000
C    6.505578841  -0.653193161   0.000000000
C    5.155734205  -1.104920876   0.000000000
C    4.821074536  -2.462094972   0.000000000
C    5.863266842  -3.376178892   0.000000000
C    7.207590281  -2.951641292   0.000000000
C    7.537930160  -1.606705567   0.000000000
H    4.681606403   2.109594763   0.000000000
H    7.299575157   1.451489921   0.000000000
H    3.788685594  -2.785730808   0.000000000
H    5.641549396  -4.434822780   0.000000000
H    7.993655200  -3.694570139   0.000000000
H    8.574255937  -1.294593877   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '21-2.0')] = qcdb.Molecule("""
0 1
C   -0.052652077  -1.393225783   0.000000000
C   -0.025543347  -0.696940104  -1.208292000
C    0.026348254   0.696724226  -1.208365000
C    0.051042263   1.393657541   0.000000000
C    0.026348254   0.696724226   1.208365000
C   -0.025543347  -0.696940104   1.208292000
H   -0.097430661  -2.473655966   0.000000000
H   -0.040509756  -1.237360068  -2.144590000
H    0.050955575   1.236531293  -2.144838000
H    0.089657645   2.474412421   0.000000000
H    0.050955575   1.236531293   2.144838000
H   -0.040509756  -1.237360068   2.144590000
--
0 1
H    4.461772054   0.000000000   0.000000000
N    5.469089458   0.005056388   0.000000000
C    6.250743642   1.132604937   0.000000000
C    7.579628369   0.772354616   0.000000000
C    7.621021855  -0.653193161   0.000000000
C    6.271177219  -1.104920876   0.000000000
C    5.936517550  -2.462094972   0.000000000
C    6.978709856  -3.376178892   0.000000000
C    8.323033295  -2.951641292   0.000000000
C    8.653373174  -1.606705567   0.000000000
H    5.797049417   2.109594763   0.000000000
H    8.415018171   1.451489921   0.000000000
H    4.904128608  -2.785730808   0.000000000
H    6.756992410  -4.434822780   0.000000000
H    9.109098214  -3.694570139   0.000000000
H    9.689698951  -1.294593877   0.000000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '22-0.7')] = qcdb.Molecule("""
0 1
C         -1.44596736 -1.22106586  0.26580875
O         -0.94522991 -0.04731809 -0.20946756
H          0.00000000  0.00000000  0.00000000
C         -0.68314270 -2.12778520  1.00510901
C         -1.25779840 -3.31409098  1.45654066
C         -2.59062773 -3.60542792  1.17905167
C         -3.34850062 -2.69511685  0.44328611
C         -2.78254941 -1.50970190 -0.01328725
H          0.35278643 -1.90546397  1.22478105
H         -0.65634919 -4.00957603  2.02623132
H         -3.03299319 -4.52638433  1.53108506
H         -4.38551290 -2.90731744  0.22101793
H         -3.35788896 -0.79601701 -0.58623496
--
0 1
O          1.35604706  0.00000000  0.00000000
C          1.95453947 -1.14289879 -0.48373245
H          1.95539652  0.41760444  0.62804116
C          1.25804307 -1.86762267 -1.44721153
C          1.81729768 -3.03591279 -1.95456799
C          3.06185406 -3.47935031 -1.50964741
C          3.74916754 -2.74469642 -0.54741031
C          3.19686752 -1.57495261 -0.02943675
H          0.29401278 -1.51302849 -1.78446706
H          1.27428716 -3.60008236 -2.69989621
H          3.49051400 -4.38751129 -1.90820423
H          4.71518108 -3.07749715 -0.19400516
H          3.72884791 -1.00425164  0.72233320
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '22-0.8')] = qcdb.Molecule("""
0 1
C         -1.44596736 -1.22106586  0.26580875
O         -0.94522991 -0.04731809 -0.20946756
H          0.00000000  0.00000000  0.00000000
C         -0.68314270 -2.12778520  1.00510901
C         -1.25779840 -3.31409098  1.45654066
C         -2.59062773 -3.60542792  1.17905167
C         -3.34850062 -2.69511685  0.44328611
C         -2.78254941 -1.50970190 -0.01328725
H          0.35278643 -1.90546397  1.22478105
H         -0.65634919 -4.00957603  2.02623132
H         -3.03299319 -4.52638433  1.53108506
H         -4.38551290 -2.90731744  0.22101793
H         -3.35788896 -0.79601701 -0.58623496
--
0 1
O          1.54976807  0.00000000  0.00000000
C          2.14826048 -1.14289879 -0.48373245
H          2.14911752  0.41760444  0.62804116
C          1.45176408 -1.86762267 -1.44721153
C          2.01101869 -3.03591279 -1.95456799
C          3.25557507 -3.47935031 -1.50964741
C          3.94288855 -2.74469642 -0.54741031
C          3.39058853 -1.57495261 -0.02943675
H          0.48773379 -1.51302849 -1.78446706
H          1.46800817 -3.60008236 -2.69989621
H          3.68423500 -4.38751129 -1.90820423
H          4.90890209 -3.07749715 -0.19400516
H          3.92256892 -1.00425164  0.72233320
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '22-0.9')] = qcdb.Molecule("""
0 1
C   -1.445967355  -1.221065858   0.265808750
O   -0.945229913  -0.047318091  -0.209467563
H    0.000000000   0.000000000   0.000000000
C   -0.683142700  -2.127785201   1.005109011
C   -1.257798399  -3.314090975   1.456540663
C   -2.590627730  -3.605427919   1.179051667
C   -3.348500619  -2.695116849   0.443286115
C   -2.782549405  -1.509701903  -0.013287247
H    0.352786431  -1.905463972   1.224781047
H   -0.656349187  -4.009576034   2.026231320
H   -3.032993188  -4.526384329   1.531085059
H   -4.385512900  -2.907317436   0.221017935
H   -3.357888956  -0.796017014  -0.586234960
--
0 1
O    1.743489077   0.000000000   0.000000000
C    2.341981491  -1.142898789  -0.483732445
H    2.342838533   0.417604441   0.628041164
C    1.645485086  -1.867622674  -1.447211527
C    2.204739700  -3.035912794  -1.954567993
C    3.449296078  -3.479350313  -1.509647408
C    4.136609561  -2.744696418  -0.547410307
C    3.584309534  -1.574952605  -0.029436748
H    0.681454799  -1.513028491  -1.784467064
H    1.661729182  -3.600082357  -2.699896207
H    3.877956013  -4.387511286  -1.908204233
H    5.102623102  -3.077497147  -0.194005162
H    4.116289930  -1.004251641   0.722333197
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '22-1.0')] = qcdb.Molecule("""
0 1
C   -2.007106000   0.763846000  -0.108351000
O   -1.388504000   1.929852000  -0.443121000
H   -0.523812000   1.964652000  -0.006461000
C   -1.463081000  -0.151912000   0.794993000
C   -2.147579000  -1.329509000   1.088368000
C   -3.374321000  -1.603143000   0.489586000
C   -3.914373000  -0.683855000  -0.409103000
C   -3.237050000   0.492961000  -0.709613000
H   -0.510651000   0.056657000   1.264256000
H   -1.715113000  -2.032145000   1.787842000
H   -3.902466000  -2.517387000   0.719795000
H   -4.867073000  -0.882294000  -0.881132000
H   -3.643166000   1.213434000  -1.405759000
--
0 1
O    1.353117000   1.938272000   0.472313000
C    2.036975000   0.786504000   0.149549000
H    1.784285000   2.348749000   1.229711000
C    1.590403000   0.069686000  -0.957415000
C    2.241737000  -1.106977000  -1.312811000
C    3.331567000  -1.566560000  -0.574864000
C    3.769684000  -0.839690000   0.528644000
C    3.122484000   0.338350000   0.896049000
H    0.744551000   0.436798000  -1.521858000
H    1.892146000  -1.664973000  -2.170184000
H    3.833023000  -2.481154000  -0.856667000
H    4.613763000  -1.185010000   1.109263000
H    3.459885000   0.903038000   1.756949000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '22-1.2')] = qcdb.Molecule("""
0 1
C   -1.445967355  -1.221065858   0.265808750
O   -0.945229913  -0.047318091  -0.209467563
H    0.000000000   0.000000000   0.000000000
C   -0.683142700  -2.127785201   1.005109011
C   -1.257798399  -3.314090975   1.456540663
C   -2.590627730  -3.605427919   1.179051667
C   -3.348500619  -2.695116849   0.443286115
C   -2.782549405  -1.509701903  -0.013287247
H    0.352786431  -1.905463972   1.224781047
H   -0.656349187  -4.009576034   2.026231320
H   -3.032993188  -4.526384329   1.531085059
H   -4.385512900  -2.907317436   0.221017935
H   -3.357888956  -0.796017014  -0.586234960
--
0 1
O    2.324652103   0.000000000   0.000000000
C    2.923144517  -1.142898789  -0.483732445
H    2.924001559   0.417604441   0.628041164
C    2.226648112  -1.867622674  -1.447211527
C    2.785902726  -3.035912794  -1.954567993
C    4.030459104  -3.479350313  -1.509647408
C    4.717772587  -2.744696418  -0.547410307
C    4.165472560  -1.574952605  -0.029436748
H    1.262617825  -1.513028491  -1.784467064
H    2.242892208  -3.600082357  -2.699896207
H    4.459119039  -4.387511286  -1.908204233
H    5.683786128  -3.077497147  -0.194005162
H    4.697452956  -1.004251641   0.722333197
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '22-1.5')] = qcdb.Molecule("""
0 1
C   -1.445967355  -1.221065858   0.265808750
O   -0.945229913  -0.047318091  -0.209467563
H    0.000000000   0.000000000   0.000000000
C   -0.683142700  -2.127785201   1.005109011
C   -1.257798399  -3.314090975   1.456540663
C   -2.590627730  -3.605427919   1.179051667
C   -3.348500619  -2.695116849   0.443286115
C   -2.782549405  -1.509701903  -0.013287247
H    0.352786431  -1.905463972   1.224781047
H   -0.656349187  -4.009576034   2.026231320
H   -3.032993188  -4.526384329   1.531085059
H   -4.385512900  -2.907317436   0.221017935
H   -3.357888956  -0.796017014  -0.586234960
--
0 1
O    2.905815129   0.000000000   0.000000000
C    3.504307543  -1.142898789  -0.483732445
H    3.505164585   0.417604441   0.628041164
C    2.807811138  -1.867622674  -1.447211527
C    3.367065752  -3.035912794  -1.954567993
C    4.611622130  -3.479350313  -1.509647408
C    5.298935613  -2.744696418  -0.547410307
C    4.746635586  -1.574952605  -0.029436748
H    1.843780851  -1.513028491  -1.784467064
H    2.824055234  -3.600082357  -2.699896207
H    5.040282065  -4.387511286  -1.908204233
H    6.264949154  -3.077497147  -0.194005162
H    5.278615982  -1.004251641   0.722333197
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '22-2.0')] = qcdb.Molecule("""
0 1
C   -1.445967355  -1.221065858   0.265808750
O   -0.945229913  -0.047318091  -0.209467563
H    0.000000000   0.000000000   0.000000000
C   -0.683142700  -2.127785201   1.005109011
C   -1.257798399  -3.314090975   1.456540663
C   -2.590627730  -3.605427919   1.179051667
C   -3.348500619  -2.695116849   0.443286115
C   -2.782549405  -1.509701903  -0.013287247
H    0.352786431  -1.905463972   1.224781047
H   -0.656349187  -4.009576034   2.026231320
H   -3.032993188  -4.526384329   1.531085059
H   -4.385512900  -2.907317436   0.221017935
H   -3.357888956  -0.796017014  -0.586234960
--
0 1
O    3.874420172   0.000000000   0.000000000
C    4.472912586  -1.142898789  -0.483732445
H    4.473769628   0.417604441   0.628041164
C    3.776416181  -1.867622674  -1.447211527
C    4.335670795  -3.035912794  -1.954567993
C    5.580227173  -3.479350313  -1.509647408
C    6.267540656  -2.744696418  -0.547410307
C    5.715240629  -1.574952605  -0.029436748
H    2.812385894  -1.513028491  -1.784467064
H    3.792660277  -3.600082357  -2.699896207
H    6.008887108  -4.387511286  -1.908204233
H    7.233554197  -3.077497147  -0.194005162
H    6.247221025  -1.004251641   0.722333197
units angstrom
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

DATA['NUCLEAR REPULSION ENERGY'] = {}
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-0.9-dimer'             ] =      41.68443604
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-0.9-monoA-unCP'        ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-0.9-monoB-unCP'        ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-1.0-dimer'             ] =      40.31423984
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-1.0-monoA-unCP'        ] =      11.94743172
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-1.0-monoB-unCP'        ] =      11.94743172
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-1.2-dimer'             ] =      38.12133822
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-1.2-monoA-unCP'        ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-1.2-monoB-unCP'        ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-1.5-dimer'             ] =      35.74458853
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-1.5-monoA-unCP'        ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-1.5-monoB-unCP'        ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-2.0-dimer'             ] =      33.16044361
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-2.0-monoA-unCP'        ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-2.0-monoB-unCP'        ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-0.9-dimer'             ] =      38.01573644
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-0.9-monoA-unCP'        ] =       9.16383015
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-0.9-monoB-unCP'        ] =       9.17803891
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-1.0-dimer'             ] =      36.66284785
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-1.0-monoA-unCP'        ] =       9.16383015
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-1.0-monoB-unCP'        ] =       9.17803890
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-1.2-dimer'             ] =      34.45708029
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-1.2-monoA-unCP'        ] =       9.16383015
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-1.2-monoB-unCP'        ] =       9.17803891
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-1.5-dimer'             ] =      32.00182680
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-1.5-monoA-unCP'        ] =       9.16383015
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-1.5-monoB-unCP'        ] =       9.17803891
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-2.0-dimer'             ] =      29.24367241
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-2.0-monoA-unCP'        ] =       9.16383015
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-2.0-monoB-unCP'        ] =       9.17803891
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-0.9-dimer'             ] =     241.06935387
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-0.9-monoA-unCP'        ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-0.9-monoB-unCP'        ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-1.0-dimer'             ] =     235.94662032
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-1.0-monoA-unCP'        ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-1.0-monoB-unCP'        ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-1.2-dimer'             ] =     227.13906126
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-1.2-monoA-unCP'        ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-1.2-monoB-unCP'        ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-1.5-dimer'             ] =     216.59223679
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-1.5-monoA-unCP'        ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-1.5-monoB-unCP'        ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-2.0-dimer'             ] =     203.69065134
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-2.0-monoA-unCP'        ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-2.0-monoB-unCP'        ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-0.9-dimer'             ] =     235.63918473
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-0.9-monoA-unCP'        ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-0.9-monoB-unCP'        ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-1.0-dimer'             ] =     230.79485521
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-1.0-monoA-unCP'        ] =      71.07286375
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-1.0-monoB-unCP'        ] =      71.07286375
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-1.2-dimer'             ] =     222.48256856
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-1.2-monoA-unCP'        ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-1.2-monoB-unCP'        ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-1.5-dimer'             ] =     212.56291415
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-1.5-monoA-unCP'        ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-1.5-monoB-unCP'        ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-2.0-dimer'             ] =     200.48924225
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-2.0-monoA-unCP'        ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-2.0-monoB-unCP'        ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-0.9-dimer'             ] =    1043.41428619
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-0.9-monoA-unCP'        ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-0.9-monoB-unCP'        ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-1.0-dimer'             ] =    1032.28187517
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-1.0-monoA-unCP'        ] =     357.22675307
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-1.0-monoB-unCP'        ] =     357.22675307
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-1.2-dimer'             ] =    1012.32214892
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-1.2-monoA-unCP'        ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-1.2-monoB-unCP'        ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-1.5-dimer'             ] =     986.94222381
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-1.5-monoA-unCP'        ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-1.5-monoB-unCP'        ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-2.0-dimer'             ] =     953.38226556
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-2.0-monoA-unCP'        ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-2.0-monoB-unCP'        ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-0.9-dimer'             ] =     822.43713935
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-0.9-monoA-unCP'        ] =     275.70186301
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-0.9-monoB-unCP'        ] =     275.67198279
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-1.0-dimer'             ] =     812.28851500
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-1.0-monoA-unCP'        ] =     275.70186300
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-1.0-monoB-unCP'        ] =     275.67198277
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-1.2-dimer'             ] =     794.18088651
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-1.2-monoA-unCP'        ] =     275.70186301
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-1.2-monoB-unCP'        ] =     275.67198279
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-1.5-dimer'             ] =     771.37080294
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-1.5-monoA-unCP'        ] =     275.70186301
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-1.5-monoB-unCP'        ] =     275.67198279
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-2.0-dimer'             ] =     741.67097558
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-2.0-monoA-unCP'        ] =     275.70186301
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-2.0-monoB-unCP'        ] =     275.67198279
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-0.9-dimer'             ] =    1379.46455665
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-0.9-monoA-unCP'        ] =     503.39628584
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-0.9-monoB-unCP'        ] =     440.30157444
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-1.0-dimer'             ] =    1365.23225970
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-1.0-monoA-unCP'        ] =     503.39628585
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-1.0-monoB-unCP'        ] =     440.30157446
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-1.2-dimer'             ] =    1339.53568799
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-1.2-monoA-unCP'        ] =     503.39628584
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-1.2-monoB-unCP'        ] =     440.30157444
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-1.5-dimer'             ] =    1306.57985526
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-1.5-monoA-unCP'        ] =     503.39628584
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-1.5-monoB-unCP'        ] =     440.30157444
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-2.0-dimer'             ] =    1262.60816943
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-2.0-monoA-unCP'        ] =     503.39628584
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-2.0-monoB-unCP'        ] =     440.30157444
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-0.9-dimer'             ] =      42.51958713
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-0.9-monoA-unCP'        ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-0.9-monoB-unCP'        ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-1.0-dimer'             ] =      41.00026380
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-1.0-monoA-unCP'        ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-1.0-monoB-unCP'        ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-1.2-dimer'             ] =      38.69237776
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-1.2-monoA-unCP'        ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-1.2-monoB-unCP'        ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-1.5-dimer'             ] =      36.35739726
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-1.5-monoA-unCP'        ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-1.5-monoB-unCP'        ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-2.0-dimer'             ] =      34.00360791
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-2.0-monoA-unCP'        ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-2.0-monoB-unCP'        ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-0.9-dimer'             ] =     105.78748629
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-0.9-monoA-unCP'        ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-0.9-monoB-unCP'        ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-1.0-dimer'             ] =     102.16530928
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-1.0-monoA-unCP'        ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-1.0-monoB-unCP'        ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-1.2-dimer'             ] =      96.54599785
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-1.2-monoA-unCP'        ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-1.2-monoB-unCP'        ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-1.5-dimer'             ] =      90.75195168
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-1.5-monoA-unCP'        ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-1.5-monoB-unCP'        ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-2.0-dimer'             ] =      84.83436234
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-2.0-monoA-unCP'        ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-2.0-monoB-unCP'        ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-0.9-dimer'            ] =     276.01204527
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-0.9-monoA-unCP'       ] =     203.70797334
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-0.9-monoB-unCP'       ] =      13.48552804
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-1.0-dimer'            ] =     272.46180693
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-1.0-monoA-unCP'       ] =     203.70797334
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-1.0-monoB-unCP'       ] =      13.48552804
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-1.2-dimer'            ] =     266.43391366
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-1.2-monoA-unCP'       ] =     203.70797334
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-1.2-monoB-unCP'       ] =      13.48552804
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-1.5-dimer'            ] =     259.41406121
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-1.5-monoA-unCP'       ] =     203.70797334
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-1.5-monoB-unCP'       ] =      13.48552804
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-2.0-dimer'            ] =     251.20587713
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-2.0-monoA-unCP'       ] =     203.70797334
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-2.0-monoB-unCP'       ] =      13.48552804
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-0.9-dimer'            ] =     648.07922043
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-0.9-monoA-unCP'       ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-0.9-monoB-unCP'       ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-1.0-dimer'            ] =     628.97202476
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-1.0-monoA-unCP'       ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-1.0-monoB-unCP'       ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-1.2-dimer'            ] =     597.97029184
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-1.2-monoA-unCP'       ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-1.2-monoB-unCP'       ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-1.5-dimer'            ] =     564.14537638
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-1.5-monoA-unCP'       ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-1.5-monoB-unCP'       ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-2.0-dimer'            ] =     527.66680634
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-2.0-monoA-unCP'       ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-2.0-monoB-unCP'       ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-0.9-dimer'            ] =     674.23986713
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-0.9-monoA-unCP'       ] =     208.63967419
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-0.9-monoB-unCP'       ] =     208.62628028
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-1.0-dimer'            ] =     654.13200064
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-1.0-monoA-unCP'       ] =     208.63967421
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-1.0-monoB-unCP'       ] =     208.62628027
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-1.2-dimer'            ] =     621.43592234
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-1.2-monoA-unCP'       ] =     208.63967419
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-1.2-monoB-unCP'       ] =     208.62628028
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-1.5-dimer'            ] =     585.61705547
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-1.5-monoA-unCP'       ] =     208.63967419
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-1.5-monoB-unCP'       ] =     208.62628028
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-2.0-dimer'            ] =     546.77330917
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-2.0-monoA-unCP'       ] =     208.63967419
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-2.0-monoB-unCP'       ] =     208.62628028
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-0.9-dimer'            ] =    1195.17642590
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-0.9-monoA-unCP'       ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-0.9-monoB-unCP'       ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-1.0-dimer'            ] =    1161.47071638
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-1.0-monoA-unCP'       ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-1.0-monoB-unCP'       ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-1.2-dimer'            ] =    1105.10369422
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-1.2-monoA-unCP'       ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-1.2-monoB-unCP'       ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-1.5-dimer'            ] =    1041.01622426
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-1.5-monoA-unCP'       ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-1.5-monoB-unCP'       ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-2.0-dimer'            ] =     968.73081054
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-2.0-monoA-unCP'       ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-2.0-monoB-unCP'       ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-0.9-dimer'            ] =     958.32945282
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-0.9-monoA-unCP'       ] =     203.66956609
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-0.9-monoB-unCP'       ] =     401.14359309
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-1.0-dimer'            ] =     935.53014764
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-1.0-monoA-unCP'       ] =     203.66956608
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-1.0-monoB-unCP'       ] =     401.14359309
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-1.2-dimer'            ] =     896.86323089
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-1.2-monoA-unCP'       ] =     203.66956609
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-1.2-monoB-unCP'       ] =     401.14359309
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-1.5-dimer'            ] =     852.00454788
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-1.5-monoA-unCP'       ] =     203.66956609
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-1.5-monoB-unCP'       ] =     401.14359309
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-2.0-dimer'            ] =     800.10296417
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-2.0-monoA-unCP'       ] =     203.66956609
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-2.0-monoB-unCP'       ] =     401.14359309
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-0.9-dimer'            ] =    1583.70519732
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-0.9-monoA-unCP'       ] =     503.36563835
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-0.9-monoB-unCP'       ] =     440.14698892
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-1.0-dimer'            ] =    1542.14304855
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-1.0-monoA-unCP'       ] =     503.36563836
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-1.0-monoB-unCP'       ] =     440.14698895
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-1.2-dimer'            ] =    1471.85303122
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-1.2-monoA-unCP'       ] =     503.36563835
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-1.2-monoB-unCP'       ] =     440.14698892
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-1.5-dimer'            ] =    1390.54653938
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-1.5-monoA-unCP'       ] =     503.36563835
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-1.5-monoB-unCP'       ] =     440.14698892
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-2.0-dimer'            ] =    1296.67402837
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-2.0-monoA-unCP'       ] =     503.36563835
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-2.0-monoB-unCP'       ] =     440.14698892
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-0.9-dimer'            ] =      87.03519106
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-0.9-monoA-unCP'       ] =      33.35807208
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-0.9-monoB-unCP'       ] =      24.69794610
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-1.0-dimer'            ] =      85.18906420
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-1.0-monoA-unCP'       ] =      33.35807208
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-1.0-monoB-unCP'       ] =      24.69794610
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-1.2-dimer'            ] =      82.12691022
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-1.2-monoA-unCP'       ] =      33.35807208
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-1.2-monoB-unCP'       ] =      24.69794610
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-1.5-dimer'            ] =      78.64799841
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-1.5-monoA-unCP'       ] =      33.35807208
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-1.5-monoB-unCP'       ] =      24.69794610
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-2.0-dimer'            ] =      74.65765021
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-2.0-monoA-unCP'       ] =      33.35807208
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-2.0-monoB-unCP'       ] =      24.69794610
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-0.9-dimer'            ] =     277.22718975
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-0.9-monoA-unCP'       ] =     203.63369790
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-0.9-monoB-unCP'       ] =       9.16734035
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-1.0-dimer'            ] =     273.32940378
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-1.0-monoA-unCP'       ] =     203.63369789
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-1.0-monoB-unCP'       ] =       9.16734036
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-1.2-dimer'            ] =     266.69472124
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-1.2-monoA-unCP'       ] =     203.63369790
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-1.2-monoB-unCP'       ] =       9.16734035
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-1.5-dimer'            ] =     258.94739117
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-1.5-monoA-unCP'       ] =     203.63369790
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-1.5-monoB-unCP'       ] =       9.16734035
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-2.0-dimer'            ] =     249.87743360
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-2.0-monoA-unCP'       ] =     203.63369790
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-2.0-monoB-unCP'       ] =       9.16734035
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-0.9-dimer'            ] =     276.96602136
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-0.9-monoA-unCP'       ] =     203.67277418
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-0.9-monoB-unCP'       ] =      11.96105518
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-1.0-dimer'            ] =     273.27963627
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-1.0-monoA-unCP'       ] =     203.67277417
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-1.0-monoB-unCP'       ] =      11.96105518
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-1.2-dimer'            ] =     266.99529350
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-1.2-monoA-unCP'       ] =     203.67277418
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-1.2-monoB-unCP'       ] =      11.96105518
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-1.5-dimer'            ] =     259.64238957
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-1.5-monoA-unCP'       ] =     203.67277418
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-1.5-monoB-unCP'       ] =      11.96105518
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-2.0-dimer'            ] =     251.01630493
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-2.0-monoA-unCP'       ] =     203.67277418
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-2.0-monoB-unCP'       ] =      11.96105518
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-0.9-dimer'            ] =     307.56164739
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-0.9-monoA-unCP'       ] =     203.59513507
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-0.9-monoB-unCP'       ] =      23.66987363
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-1.0-dimer'            ] =     303.28139253
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-1.0-monoA-unCP'       ] =     203.59513507
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-1.0-monoB-unCP'       ] =      23.66987364
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-1.2-dimer'            ] =     295.87805947
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-1.2-monoA-unCP'       ] =     203.59513507
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-1.2-monoB-unCP'       ] =      23.66987363
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-1.5-dimer'            ] =     287.02868743
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-1.5-monoA-unCP'       ] =     203.59513507
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-1.5-monoB-unCP'       ] =      23.66987363
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-2.0-dimer'            ] =     276.34590865
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-2.0-monoA-unCP'       ] =     203.59513507
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-2.0-monoB-unCP'       ] =      23.66987363
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-0.9-dimer'            ] =     601.46920410
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-0.9-monoA-unCP'       ] =     203.68142723
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-0.9-monoB-unCP'       ] =     203.66408617
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-1.0-dimer'            ] =     592.41663921
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-1.0-monoA-unCP'       ] =     203.68142723
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-1.0-monoB-unCP'       ] =     203.66408617
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-1.2-dimer'            ] =     576.54530095
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-1.2-monoA-unCP'       ] =     203.68142723
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-1.2-monoB-unCP'       ] =     203.66408617
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-1.5-dimer'            ] =     557.14862254
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-1.5-monoA-unCP'       ] =     203.68142723
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-1.5-monoB-unCP'       ] =     203.66408617
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-2.0-dimer'            ] =     532.99122415
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-2.0-monoA-unCP'       ] =     203.68142723
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-2.0-monoB-unCP'       ] =     203.66408617
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-0.9-dimer'            ] =     888.79508333
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-0.9-monoA-unCP'       ] =     203.56579265
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-0.9-monoB-unCP'       ] =     401.05660150
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-1.0-dimer'            ] =     876.91918503
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-1.0-monoA-unCP'       ] =     203.56579265
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-1.0-monoB-unCP'       ] =     401.05660150
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-1.2-dimer'            ] =     855.79228809
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-1.2-monoA-unCP'       ] =     203.56579265
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-1.2-monoB-unCP'       ] =     401.05660150
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-1.5-dimer'            ] =     829.42534245
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-1.5-monoA-unCP'       ] =     203.56579265
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-1.5-monoB-unCP'       ] =     401.05660150
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-2.0-dimer'            ] =     795.71041545
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-2.0-monoA-unCP'       ] =     203.56579265
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-2.0-monoB-unCP'       ] =     401.05660150
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-0.9-dimer'            ] =     814.74763476
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-0.9-monoA-unCP'       ] =     271.43868469
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-0.9-monoB-unCP'       ] =     271.34619735
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-1.0-dimer'            ] =     805.11772632
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-1.0-monoA-unCP'       ] =     271.43868470
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-1.0-monoB-unCP'       ] =     271.34619734
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-1.2-dimer'            ] =     787.64120113
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-1.2-monoA-unCP'       ] =     271.43868469
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-1.2-monoB-unCP'       ] =     271.34619735
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-1.5-dimer'            ] =     765.14830882
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-1.5-monoA-unCP'       ] =     271.43868469
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-1.5-monoB-unCP'       ] =     271.34619735
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-2.0-dimer'            ] =     735.23075037
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-2.0-monoA-unCP'       ] =     271.43868469
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-2.0-monoB-unCP'       ] =     271.34619735
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-0.9-monoA-CP'          ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-0.9-monoB-CP'          ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-1.0-monoA-CP'          ] =      11.94743172
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-1.0-monoB-CP'          ] =      11.94743172
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-1.2-monoA-CP'          ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-1.2-monoB-CP'          ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-1.5-monoA-CP'          ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-1.5-monoB-CP'          ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-2.0-monoA-CP'          ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-1-2.0-monoB-CP'          ] =      11.94743173
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-0.9-monoA-CP'          ] =       9.16383015
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-0.9-monoB-CP'          ] =       9.17803891
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-1.0-monoA-CP'          ] =       9.16383015
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-1.0-monoB-CP'          ] =       9.17803890
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-1.2-monoA-CP'          ] =       9.16383015
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-1.2-monoB-CP'          ] =       9.17803891
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-1.5-monoA-CP'          ] =       9.16383015
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-1.5-monoB-CP'          ] =       9.17803891
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-2.0-monoA-CP'          ] =       9.16383015
DATA['NUCLEAR REPULSION ENERGY']['S22by5-2-2.0-monoB-CP'          ] =       9.17803891
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-0.9-monoA-CP'          ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-0.9-monoB-CP'          ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-1.0-monoA-CP'          ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-1.0-monoB-CP'          ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-1.2-monoA-CP'          ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-1.2-monoB-CP'          ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-1.5-monoA-CP'          ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-1.5-monoB-CP'          ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-2.0-monoA-CP'          ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-3-2.0-monoB-CP'          ] =      70.11578330
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-0.9-monoA-CP'          ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-0.9-monoB-CP'          ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-1.0-monoA-CP'          ] =      71.07286375
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-1.0-monoB-CP'          ] =      71.07286375
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-1.2-monoA-CP'          ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-1.2-monoB-CP'          ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-1.5-monoA-CP'          ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-1.5-monoB-CP'          ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-2.0-monoA-CP'          ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-4-2.0-monoB-CP'          ] =      71.07286374
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-0.9-monoA-CP'          ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-0.9-monoB-CP'          ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-1.0-monoA-CP'          ] =     357.22675307
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-1.0-monoB-CP'          ] =     357.22675307
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-1.2-monoA-CP'          ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-1.2-monoB-CP'          ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-1.5-monoA-CP'          ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-1.5-monoB-CP'          ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-2.0-monoA-CP'          ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-5-2.0-monoB-CP'          ] =     357.22675306
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-0.9-monoA-CP'          ] =     275.70186301
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-0.9-monoB-CP'          ] =     275.67198279
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-1.0-monoA-CP'          ] =     275.70186300
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-1.0-monoB-CP'          ] =     275.67198277
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-1.2-monoA-CP'          ] =     275.70186301
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-1.2-monoB-CP'          ] =     275.67198279
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-1.5-monoA-CP'          ] =     275.70186301
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-1.5-monoB-CP'          ] =     275.67198279
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-2.0-monoA-CP'          ] =     275.70186301
DATA['NUCLEAR REPULSION ENERGY']['S22by5-6-2.0-monoB-CP'          ] =     275.67198279
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-0.9-monoA-CP'          ] =     503.39628584
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-0.9-monoB-CP'          ] =     440.30157444
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-1.0-monoA-CP'          ] =     503.39628585
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-1.0-monoB-CP'          ] =     440.30157446
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-1.2-monoA-CP'          ] =     503.39628584
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-1.2-monoB-CP'          ] =     440.30157444
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-1.5-monoA-CP'          ] =     503.39628584
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-1.5-monoB-CP'          ] =     440.30157444
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-2.0-monoA-CP'          ] =     503.39628584
DATA['NUCLEAR REPULSION ENERGY']['S22by5-7-2.0-monoB-CP'          ] =     440.30157444
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-0.9-monoA-CP'          ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-0.9-monoB-CP'          ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-1.0-monoA-CP'          ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-1.0-monoB-CP'          ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-1.2-monoA-CP'          ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-1.2-monoB-CP'          ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-1.5-monoA-CP'          ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-1.5-monoB-CP'          ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-2.0-monoA-CP'          ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-8-2.0-monoB-CP'          ] =      13.44804227
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-0.9-monoA-CP'          ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-0.9-monoB-CP'          ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-1.0-monoA-CP'          ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-1.0-monoB-CP'          ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-1.2-monoA-CP'          ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-1.2-monoB-CP'          ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-1.5-monoA-CP'          ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-1.5-monoB-CP'          ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-2.0-monoA-CP'          ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-9-2.0-monoB-CP'          ] =      33.36026958
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-0.9-monoA-CP'         ] =     203.70797334
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-0.9-monoB-CP'         ] =      13.48552804
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-1.0-monoA-CP'         ] =     203.70797334
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-1.0-monoB-CP'         ] =      13.48552804
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-1.2-monoA-CP'         ] =     203.70797334
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-1.2-monoB-CP'         ] =      13.48552804
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-1.5-monoA-CP'         ] =     203.70797334
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-1.5-monoB-CP'         ] =      13.48552804
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-2.0-monoA-CP'         ] =     203.70797334
DATA['NUCLEAR REPULSION ENERGY']['S22by5-10-2.0-monoB-CP'         ] =      13.48552804
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-0.9-monoA-CP'         ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-0.9-monoB-CP'         ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-1.0-monoA-CP'         ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-1.0-monoB-CP'         ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-1.2-monoA-CP'         ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-1.2-monoB-CP'         ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-1.5-monoA-CP'         ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-1.5-monoB-CP'         ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-2.0-monoA-CP'         ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-11-2.0-monoB-CP'         ] =     203.71090864
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-0.9-monoA-CP'         ] =     208.63967419
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-0.9-monoB-CP'         ] =     208.62628028
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-1.0-monoA-CP'         ] =     208.63967421
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-1.0-monoB-CP'         ] =     208.62628027
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-1.2-monoA-CP'         ] =     208.63967419
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-1.2-monoB-CP'         ] =     208.62628028
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-1.5-monoA-CP'         ] =     208.63967419
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-1.5-monoB-CP'         ] =     208.62628028
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-2.0-monoA-CP'         ] =     208.63967419
DATA['NUCLEAR REPULSION ENERGY']['S22by5-12-2.0-monoB-CP'         ] =     208.62628028
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-0.9-monoA-CP'         ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-0.9-monoB-CP'         ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-1.0-monoA-CP'         ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-1.0-monoB-CP'         ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-1.2-monoA-CP'         ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-1.2-monoB-CP'         ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-1.5-monoA-CP'         ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-1.5-monoB-CP'         ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-2.0-monoA-CP'         ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-13-2.0-monoB-CP'         ] =     357.16045924
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-0.9-monoA-CP'         ] =     203.66956609
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-0.9-monoB-CP'         ] =     401.14359309
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-1.0-monoA-CP'         ] =     203.66956608
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-1.0-monoB-CP'         ] =     401.14359309
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-1.2-monoA-CP'         ] =     203.66956609
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-1.2-monoB-CP'         ] =     401.14359309
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-1.5-monoA-CP'         ] =     203.66956609
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-1.5-monoB-CP'         ] =     401.14359309
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-2.0-monoA-CP'         ] =     203.66956609
DATA['NUCLEAR REPULSION ENERGY']['S22by5-14-2.0-monoB-CP'         ] =     401.14359309
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-0.9-monoA-CP'         ] =     503.36563835
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-0.9-monoB-CP'         ] =     440.14698892
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-1.0-monoA-CP'         ] =     503.36563836
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-1.0-monoB-CP'         ] =     440.14698895
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-1.2-monoA-CP'         ] =     503.36563835
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-1.2-monoB-CP'         ] =     440.14698892
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-1.5-monoA-CP'         ] =     503.36563835
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-1.5-monoB-CP'         ] =     440.14698892
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-2.0-monoA-CP'         ] =     503.36563835
DATA['NUCLEAR REPULSION ENERGY']['S22by5-15-2.0-monoB-CP'         ] =     440.14698892
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-0.9-monoA-CP'         ] =      33.35807208
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-0.9-monoB-CP'         ] =      24.69794610
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-1.0-monoA-CP'         ] =      33.35807208
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-1.0-monoB-CP'         ] =      24.69794610
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-1.2-monoA-CP'         ] =      33.35807208
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-1.2-monoB-CP'         ] =      24.69794610
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-1.5-monoA-CP'         ] =      33.35807208
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-1.5-monoB-CP'         ] =      24.69794610
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-2.0-monoA-CP'         ] =      33.35807208
DATA['NUCLEAR REPULSION ENERGY']['S22by5-16-2.0-monoB-CP'         ] =      24.69794610
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-0.9-monoA-CP'         ] =     203.63369790
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-0.9-monoB-CP'         ] =       9.16734035
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-1.0-monoA-CP'         ] =     203.63369789
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-1.0-monoB-CP'         ] =       9.16734036
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-1.2-monoA-CP'         ] =     203.63369790
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-1.2-monoB-CP'         ] =       9.16734035
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-1.5-monoA-CP'         ] =     203.63369790
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-1.5-monoB-CP'         ] =       9.16734035
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-2.0-monoA-CP'         ] =     203.63369790
DATA['NUCLEAR REPULSION ENERGY']['S22by5-17-2.0-monoB-CP'         ] =       9.16734035
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-0.9-monoA-CP'         ] =     203.67277418
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-0.9-monoB-CP'         ] =      11.96105518
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-1.0-monoA-CP'         ] =     203.67277417
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-1.0-monoB-CP'         ] =      11.96105518
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-1.2-monoA-CP'         ] =     203.67277418
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-1.2-monoB-CP'         ] =      11.96105518
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-1.5-monoA-CP'         ] =     203.67277418
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-1.5-monoB-CP'         ] =      11.96105518
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-2.0-monoA-CP'         ] =     203.67277418
DATA['NUCLEAR REPULSION ENERGY']['S22by5-18-2.0-monoB-CP'         ] =      11.96105518
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-0.9-monoA-CP'         ] =     203.59513507
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-0.9-monoB-CP'         ] =      23.66987363
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-1.0-monoA-CP'         ] =     203.59513507
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-1.0-monoB-CP'         ] =      23.66987364
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-1.2-monoA-CP'         ] =     203.59513507
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-1.2-monoB-CP'         ] =      23.66987363
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-1.5-monoA-CP'         ] =     203.59513507
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-1.5-monoB-CP'         ] =      23.66987363
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-2.0-monoA-CP'         ] =     203.59513507
DATA['NUCLEAR REPULSION ENERGY']['S22by5-19-2.0-monoB-CP'         ] =      23.66987363
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-0.9-monoA-CP'         ] =     203.68142723
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-0.9-monoB-CP'         ] =     203.66408617
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-1.0-monoA-CP'         ] =     203.68142723
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-1.0-monoB-CP'         ] =     203.66408617
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-1.2-monoA-CP'         ] =     203.68142723
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-1.2-monoB-CP'         ] =     203.66408617
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-1.5-monoA-CP'         ] =     203.68142723
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-1.5-monoB-CP'         ] =     203.66408617
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-2.0-monoA-CP'         ] =     203.68142723
DATA['NUCLEAR REPULSION ENERGY']['S22by5-20-2.0-monoB-CP'         ] =     203.66408617
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-0.9-monoA-CP'         ] =     203.56579265
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-0.9-monoB-CP'         ] =     401.05660150
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-1.0-monoA-CP'         ] =     203.56579265
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-1.0-monoB-CP'         ] =     401.05660150
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-1.2-monoA-CP'         ] =     203.56579265
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-1.2-monoB-CP'         ] =     401.05660150
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-1.5-monoA-CP'         ] =     203.56579265
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-1.5-monoB-CP'         ] =     401.05660150
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-2.0-monoA-CP'         ] =     203.56579265
DATA['NUCLEAR REPULSION ENERGY']['S22by5-21-2.0-monoB-CP'         ] =     401.05660150
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-0.9-monoA-CP'         ] =     271.43868469
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-0.9-monoB-CP'         ] =     271.34619735
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-1.0-monoA-CP'         ] =     271.43868470
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-1.0-monoB-CP'         ] =     271.34619734
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-1.2-monoA-CP'         ] =     271.43868469
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-1.2-monoB-CP'         ] =     271.34619735
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-1.5-monoA-CP'         ] =     271.43868469
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-1.5-monoB-CP'         ] =     271.34619735
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-2.0-monoA-CP'         ] =     271.43868469
DATA['NUCLEAR REPULSION ENERGY']['S22by5-22-2.0-monoB-CP'         ] =     271.34619735
DATA['SAPT ELST ENERGY'] = {}
DATA['SAPT ELST ENERGY']['S22by7-1-0.7'] =   -27.4508
DATA['SAPT ELST ENERGY']['S22by7-1-0.8'] =   -14.7915
DATA['SAPT ELST ENERGY']['S22by7-1-0.9'] =    -8.4599
DATA['SAPT ELST ENERGY']['S22by7-1-1.0'] =    -5.1992
DATA['SAPT ELST ENERGY']['S22by7-1-1.2'] =    -2.3971
DATA['SAPT ELST ENERGY']['S22by7-1-1.5'] =    -1.0331
DATA['SAPT ELST ENERGY']['S22by7-1-2.0'] =    -0.3600
DATA['SAPT ELST ENERGY']['S22by7-2-0.7'] =   -29.8544
DATA['SAPT ELST ENERGY']['S22by7-2-0.8'] =   -19.0881
DATA['SAPT ELST ENERGY']['S22by7-2-0.9'] =   -12.6676
DATA['SAPT ELST ENERGY']['S22by7-2-1.0'] =    -8.8322
DATA['SAPT ELST ENERGY']['S22by7-2-1.2'] =    -4.8719
DATA['SAPT ELST ENERGY']['S22by7-2-1.5'] =    -2.4653
DATA['SAPT ELST ENERGY']['S22by7-2-2.0'] =    -1.0873
DATA['SAPT ELST ENERGY']['S22by7-3-0.7'] =   -91.5955
DATA['SAPT ELST ENERGY']['S22by7-3-0.8'] =   -65.3747
DATA['SAPT ELST ENERGY']['S22by7-3-0.9'] =   -46.3889
DATA['SAPT ELST ENERGY']['S22by7-3-1.0'] =   -33.4534
DATA['SAPT ELST ENERGY']['S22by7-3-1.2'] =   -18.8110
DATA['SAPT ELST ENERGY']['S22by7-3-1.5'] =    -9.2299
DATA['SAPT ELST ENERGY']['S22by7-3-2.0'] =    -3.6064
DATA['SAPT ELST ENERGY']['S22by7-4-0.7'] =  -352.8075
DATA['SAPT ELST ENERGY']['S22by7-4-0.8'] =  -340.6278
DATA['SAPT ELST ENERGY']['S22by7-4-0.9'] =   -36.9683
DATA['SAPT ELST ENERGY']['S22by7-4-1.0'] =   -26.2900
DATA['SAPT ELST ENERGY']['S22by7-4-1.2'] =   -14.9422
DATA['SAPT ELST ENERGY']['S22by7-4-1.5'] =    -7.7996
DATA['SAPT ELST ENERGY']['S22by7-4-2.0'] =    -3.5297
DATA['SAPT ELST ENERGY']['S22by7-5-0.7'] =   -89.1571
DATA['SAPT ELST ENERGY']['S22by7-5-0.8'] =   -61.6784
DATA['SAPT ELST ENERGY']['S22by7-5-0.9'] =   -43.4890
DATA['SAPT ELST ENERGY']['S22by7-5-1.0'] =   -31.7248
DATA['SAPT ELST ENERGY']['S22by7-5-1.2'] =   -18.7656
DATA['SAPT ELST ENERGY']['S22by7-5-1.5'] =   -10.1627
DATA['SAPT ELST ENERGY']['S22by7-5-2.0'] =    -4.6908
DATA['SAPT ELST ENERGY']['S22by7-6-0.7'] =   -84.3469
DATA['SAPT ELST ENERGY']['S22by7-6-0.8'] =   -57.1674
DATA['SAPT ELST ENERGY']['S22by7-6-0.9'] =   -39.1790
DATA['SAPT ELST ENERGY']['S22by7-6-1.0'] =   -27.5327
DATA['SAPT ELST ENERGY']['S22by7-6-1.2'] =   -14.9048
DATA['SAPT ELST ENERGY']['S22by7-6-1.5'] =    -7.1257
DATA['SAPT ELST ENERGY']['S22by7-6-2.0'] =    -2.8109
DATA['SAPT ELST ENERGY']['S22by7-7-0.7'] =   -83.8073
DATA['SAPT ELST ENERGY']['S22by7-7-0.8'] =   -57.0062
DATA['SAPT ELST ENERGY']['S22by7-7-0.9'] =   -39.1540
DATA['SAPT ELST ENERGY']['S22by7-7-1.0'] =   -27.4781
DATA['SAPT ELST ENERGY']['S22by7-7-1.2'] =   -14.6006
DATA['SAPT ELST ENERGY']['S22by7-7-1.5'] =    -6.5308
DATA['SAPT ELST ENERGY']['S22by7-7-2.0'] =    -2.1167
DATA['SAPT ELST ENERGY']['S22by7-8-0.7'] =    -7.4376
DATA['SAPT ELST ENERGY']['S22by7-8-0.8'] =    -2.0215
DATA['SAPT ELST ENERGY']['S22by7-8-0.9'] =    -0.5443
DATA['SAPT ELST ENERGY']['S22by7-8-1.0'] =    -0.1399
DATA['SAPT ELST ENERGY']['S22by7-8-1.2'] =    -0.0035
DATA['SAPT ELST ENERGY']['S22by7-8-1.5'] =     0.0014
DATA['SAPT ELST ENERGY']['S22by7-8-2.0'] =     0.0002
DATA['SAPT ELST ENERGY']['S22by7-9-0.7'] =   -21.9239
DATA['SAPT ELST ENERGY']['S22by7-9-0.8'] =    -8.0217
DATA['SAPT ELST ENERGY']['S22by7-9-0.9'] =    -2.9597
DATA['SAPT ELST ENERGY']['S22by7-9-1.0'] =    -1.1398
DATA['SAPT ELST ENERGY']['S22by7-9-1.2'] =    -0.2123
DATA['SAPT ELST ENERGY']['S22by7-9-1.5'] =    -0.0316
DATA['SAPT ELST ENERGY']['S22by7-9-2.0'] =    -0.0035
DATA['SAPT ELST ENERGY']['S22by7-10-0.7'] =    -8.4453
DATA['SAPT ELST ENERGY']['S22by7-10-0.8'] =    -4.1143
DATA['SAPT ELST ENERGY']['S22by7-10-0.9'] =    -2.0088
DATA['SAPT ELST ENERGY']['S22by7-10-1.0'] =    -0.9989
DATA['SAPT ELST ENERGY']['S22by7-10-1.2'] =    -0.2812
DATA['SAPT ELST ENERGY']['S22by7-10-1.5'] =    -0.0800
DATA['SAPT ELST ENERGY']['S22by7-10-2.0'] =    -0.0257
DATA['SAPT ELST ENERGY']['S22by7-11-0.7'] =   -80.7532
DATA['SAPT ELST ENERGY']['S22by7-11-0.8'] =   -28.9536
DATA['SAPT ELST ENERGY']['S22by7-11-0.9'] =    -9.8574
DATA['SAPT ELST ENERGY']['S22by7-11-1.0'] =    -2.9403
DATA['SAPT ELST ENERGY']['S22by7-11-1.2'] =     0.1788
DATA['SAPT ELST ENERGY']['S22by7-11-1.5'] =     0.2620
DATA['SAPT ELST ENERGY']['S22by7-11-2.0'] =     0.0664
DATA['SAPT ELST ENERGY']['S22by7-12-0.7'] =   -88.4220
DATA['SAPT ELST ENERGY']['S22by7-12-0.8'] =   -32.9343
DATA['SAPT ELST ENERGY']['S22by7-12-0.9'] =   -12.4365
DATA['SAPT ELST ENERGY']['S22by7-12-1.0'] =    -4.8584
DATA['SAPT ELST ENERGY']['S22by7-12-1.2'] =    -0.9494
DATA['SAPT ELST ENERGY']['S22by7-12-1.5'] =    -0.2143
DATA['SAPT ELST ENERGY']['S22by7-12-2.0'] =    -0.0541
DATA['SAPT ELST ENERGY']['S22by7-13-0.7'] =  -121.8987
DATA['SAPT ELST ENERGY']['S22by7-13-0.8'] =   -46.3659
DATA['SAPT ELST ENERGY']['S22by7-13-0.9'] =   -19.1217
DATA['SAPT ELST ENERGY']['S22by7-13-1.0'] =    -9.1074
DATA['SAPT ELST ENERGY']['S22by7-13-1.2'] =    -3.4398
DATA['SAPT ELST ENERGY']['S22by7-13-1.5'] =    -1.5690
DATA['SAPT ELST ENERGY']['S22by7-13-2.0'] =    -0.6242
DATA['SAPT ELST ENERGY']['S22by7-14-0.7'] =   -92.9501
DATA['SAPT ELST ENERGY']['S22by7-14-0.8'] =   -36.9642
DATA['SAPT ELST ENERGY']['S22by7-14-0.9'] =   -14.0687
DATA['SAPT ELST ENERGY']['S22by7-14-1.0'] =    -4.8321
DATA['SAPT ELST ENERGY']['S22by7-14-1.2'] =     0.1278
DATA['SAPT ELST ENERGY']['S22by7-14-1.5'] =     0.6068
DATA['SAPT ELST ENERGY']['S22by7-14-2.0'] =     0.2903
DATA['SAPT ELST ENERGY']['S22by7-15-0.7'] =  -146.8569
DATA['SAPT ELST ENERGY']['S22by7-15-0.8'] =   -59.2073
DATA['SAPT ELST ENERGY']['S22by7-15-0.9'] =   -24.9484
DATA['SAPT ELST ENERGY']['S22by7-15-1.0'] =   -11.5015
DATA['SAPT ELST ENERGY']['S22by7-15-1.2'] =    -3.7393
DATA['SAPT ELST ENERGY']['S22by7-15-1.5'] =    -1.5563
DATA['SAPT ELST ENERGY']['S22by7-15-2.0'] =    -0.6233
DATA['SAPT ELST ENERGY']['S22by7-16-0.7'] =   -10.5236
DATA['SAPT ELST ENERGY']['S22by7-16-0.8'] =    -5.7851
DATA['SAPT ELST ENERGY']['S22by7-16-0.9'] =    -3.3850
DATA['SAPT ELST ENERGY']['S22by7-16-1.0'] =    -2.1200
DATA['SAPT ELST ENERGY']['S22by7-16-1.2'] =    -0.9890
DATA['SAPT ELST ENERGY']['S22by7-16-1.5'] =    -0.4134
DATA['SAPT ELST ENERGY']['S22by7-16-2.0'] =    -0.1347
DATA['SAPT ELST ENERGY']['S22by7-17-0.7'] =   -10.7239
DATA['SAPT ELST ENERGY']['S22by7-17-0.8'] =    -6.5632
DATA['SAPT ELST ENERGY']['S22by7-17-0.9'] =    -4.3823
DATA['SAPT ELST ENERGY']['S22by7-17-1.0'] =    -3.1953
DATA['SAPT ELST ENERGY']['S22by7-17-1.2'] =    -1.9958
DATA['SAPT ELST ENERGY']['S22by7-17-1.5'] =    -1.1453
DATA['SAPT ELST ENERGY']['S22by7-17-2.0'] =    -0.5174
DATA['SAPT ELST ENERGY']['S22by7-18-0.7'] =    -9.4801
DATA['SAPT ELST ENERGY']['S22by7-18-0.8'] =    -5.1597
DATA['SAPT ELST ENERGY']['S22by7-18-0.9'] =    -3.0489
DATA['SAPT ELST ENERGY']['S22by7-18-1.0'] =    -1.9880
DATA['SAPT ELST ENERGY']['S22by7-18-1.2'] =    -1.0648
DATA['SAPT ELST ENERGY']['S22by7-18-1.5'] =    -0.5504
DATA['SAPT ELST ENERGY']['S22by7-18-2.0'] =    -0.2231
DATA['SAPT ELST ENERGY']['S22by7-19-0.7'] =   -12.2280
DATA['SAPT ELST ENERGY']['S22by7-19-0.8'] =    -8.1149
DATA['SAPT ELST ENERGY']['S22by7-19-0.9'] =    -5.8015
DATA['SAPT ELST ENERGY']['S22by7-19-1.0'] =    -4.4374
DATA['SAPT ELST ENERGY']['S22by7-19-1.2'] =    -2.9398
DATA['SAPT ELST ENERGY']['S22by7-19-1.5'] =    -1.7996
DATA['SAPT ELST ENERGY']['S22by7-19-2.0'] =    -0.8624
DATA['SAPT ELST ENERGY']['S22by7-20-0.7'] =   -12.5300
DATA['SAPT ELST ENERGY']['S22by7-20-0.8'] =    -6.7388
DATA['SAPT ELST ENERGY']['S22by7-20-0.9'] =    -3.7099
DATA['SAPT ELST ENERGY']['S22by7-20-1.0'] =    -2.1426
DATA['SAPT ELST ENERGY']['S22by7-20-1.2'] =    -0.8815
DATA['SAPT ELST ENERGY']['S22by7-20-1.5'] =    -0.3838
DATA['SAPT ELST ENERGY']['S22by7-20-2.0'] =    -0.1561
DATA['SAPT ELST ENERGY']['S22by7-21-0.7'] =  -699.8211
DATA['SAPT ELST ENERGY']['S22by7-21-0.8'] =  -695.9204
DATA['SAPT ELST ENERGY']['S22by7-21-0.9'] =    -6.3493
DATA['SAPT ELST ENERGY']['S22by7-21-1.0'] =    -4.5110
DATA['SAPT ELST ENERGY']['S22by7-21-1.2'] =    -2.7504
DATA['SAPT ELST ENERGY']['S22by7-21-1.5'] =    -1.6350
DATA['SAPT ELST ENERGY']['S22by7-21-2.0'] =    -0.7857
DATA['SAPT ELST ENERGY']['S22by7-22-0.7'] =   -31.0889
DATA['SAPT ELST ENERGY']['S22by7-22-0.8'] =   -19.9434
DATA['SAPT ELST ENERGY']['S22by7-22-0.9'] =   -13.1881
DATA['SAPT ELST ENERGY']['S22by7-22-1.0'] =    -9.0960
DATA['SAPT ELST ENERGY']['S22by7-22-1.2'] =    -4.8581
DATA['SAPT ELST ENERGY']['S22by7-22-1.5'] =    -2.3498
DATA['SAPT ELST ENERGY']['S22by7-22-2.0'] =    -0.9867
DATA['SAPT EXCH ENERGY'] = {}
DATA['SAPT EXCH ENERGY']['S22by7-1-0.7'] =    54.6575
DATA['SAPT EXCH ENERGY']['S22by7-1-0.8'] =    23.9218
DATA['SAPT EXCH ENERGY']['S22by7-1-0.9'] =    10.2933
DATA['SAPT EXCH ENERGY']['S22by7-1-1.0'] =     4.3850
DATA['SAPT EXCH ENERGY']['S22by7-1-1.2'] =     0.7859
DATA['SAPT EXCH ENERGY']['S22by7-1-1.5'] =     0.0572
DATA['SAPT EXCH ENERGY']['S22by7-1-2.0'] =     0.0006
DATA['SAPT EXCH ENERGY']['S22by7-2-0.7'] =    56.6680
DATA['SAPT EXCH ENERGY']['S22by7-2-0.8'] =    28.6313
DATA['SAPT EXCH ENERGY']['S22by7-2-0.9'] =    14.3476
DATA['SAPT EXCH ENERGY']['S22by7-2-1.0'] =     7.1453
DATA['SAPT EXCH ENERGY']['S22by7-2-1.2'] =     1.7408
DATA['SAPT EXCH ENERGY']['S22by7-2-1.5'] =     0.2033
DATA['SAPT EXCH ENERGY']['S22by7-2-2.0'] =     0.0052
DATA['SAPT EXCH ENERGY']['S22by7-3-0.7'] =   200.6937
DATA['SAPT EXCH ENERGY']['S22by7-3-0.8'] =   114.1081
DATA['SAPT EXCH ENERGY']['S22by7-3-0.9'] =    64.2381
DATA['SAPT EXCH ENERGY']['S22by7-3-1.0'] =    35.9954
DATA['SAPT EXCH ENERGY']['S22by7-3-1.2'] =    11.2342
DATA['SAPT EXCH ENERGY']['S22by7-3-1.5'] =     1.9365
DATA['SAPT EXCH ENERGY']['S22by7-3-2.0'] =     0.1036
DATA['SAPT EXCH ENERGY']['S22by7-4-0.7'] =  1013.2877
DATA['SAPT EXCH ENERGY']['S22by7-4-0.8'] =   980.7164
DATA['SAPT EXCH ENERGY']['S22by7-4-0.9'] =    46.1334
DATA['SAPT EXCH ENERGY']['S22by7-4-1.0'] =    24.5192
DATA['SAPT EXCH ENERGY']['S22by7-4-1.2'] =     6.8755
DATA['SAPT EXCH ENERGY']['S22by7-4-1.5'] =     1.0253
DATA['SAPT EXCH ENERGY']['S22by7-4-2.0'] =     0.0429
DATA['SAPT EXCH ENERGY']['S22by7-5-0.7'] =   177.1874
DATA['SAPT EXCH ENERGY']['S22by7-5-0.8'] =    96.8761
DATA['SAPT EXCH ENERGY']['S22by7-5-0.9'] =    52.4165
DATA['SAPT EXCH ENERGY']['S22by7-5-1.0'] =    28.1856
DATA['SAPT EXCH ENERGY']['S22by7-5-1.2'] =     8.0463
DATA['SAPT EXCH ENERGY']['S22by7-5-1.5'] =     1.2109
DATA['SAPT EXCH ENERGY']['S22by7-5-2.0'] =     0.0499
DATA['SAPT EXCH ENERGY']['S22by7-6-0.7'] =   166.8985
DATA['SAPT EXCH ENERGY']['S22by7-6-0.8'] =    93.5326
DATA['SAPT EXCH ENERGY']['S22by7-6-0.9'] =    51.8081
DATA['SAPT EXCH ENERGY']['S22by7-6-1.0'] =    28.4753
DATA['SAPT EXCH ENERGY']['S22by7-6-1.2'] =     8.4683
DATA['SAPT EXCH ENERGY']['S22by7-6-1.5'] =     1.3486
DATA['SAPT EXCH ENERGY']['S22by7-6-2.0'] =     0.0604
DATA['SAPT EXCH ENERGY']['S22by7-7-0.7'] =   165.3393
DATA['SAPT EXCH ENERGY']['S22by7-7-0.8'] =    92.2677
DATA['SAPT EXCH ENERGY']['S22by7-7-0.9'] =    51.0380
DATA['SAPT EXCH ENERGY']['S22by7-7-1.0'] =    28.0713
DATA['SAPT EXCH ENERGY']['S22by7-7-1.2'] =     8.3980
DATA['SAPT EXCH ENERGY']['S22by7-7-1.5'] =     1.3555
DATA['SAPT EXCH ENERGY']['S22by7-7-2.0'] =     0.0626
DATA['SAPT EXCH ENERGY']['S22by7-8-0.7'] =    25.0846
DATA['SAPT EXCH ENERGY']['S22by7-8-0.8'] =     7.2243
DATA['SAPT EXCH ENERGY']['S22by7-8-0.9'] =     2.0145
DATA['SAPT EXCH ENERGY']['S22by7-8-1.0'] =     0.5470
DATA['SAPT EXCH ENERGY']['S22by7-8-1.2'] =     0.0377
DATA['SAPT EXCH ENERGY']['S22by7-8-1.5'] =     0.0007
DATA['SAPT EXCH ENERGY']['S22by7-8-2.0'] =    -0.0000
DATA['SAPT EXCH ENERGY']['S22by7-9-0.7'] =    63.4942
DATA['SAPT EXCH ENERGY']['S22by7-9-0.8'] =    22.2637
DATA['SAPT EXCH ENERGY']['S22by7-9-0.9'] =     7.3171
DATA['SAPT EXCH ENERGY']['S22by7-9-1.0'] =     2.2805
DATA['SAPT EXCH ENERGY']['S22by7-9-1.2'] =     0.1965
DATA['SAPT EXCH ENERGY']['S22by7-9-1.5'] =     0.0042
DATA['SAPT EXCH ENERGY']['S22by7-9-2.0'] =     0.0000
DATA['SAPT EXCH ENERGY']['S22by7-10-0.7'] =    24.8009
DATA['SAPT EXCH ENERGY']['S22by7-10-0.8'] =    11.8299
DATA['SAPT EXCH ENERGY']['S22by7-10-0.9'] =     5.4609
DATA['SAPT EXCH ENERGY']['S22by7-10-1.0'] =     2.4581
DATA['SAPT EXCH ENERGY']['S22by7-10-1.2'] =     0.4732
DATA['SAPT EXCH ENERGY']['S22by7-10-1.5'] =     0.0351
DATA['SAPT EXCH ENERGY']['S22by7-10-2.0'] =     0.0003
DATA['SAPT EXCH ENERGY']['S22by7-11-0.7'] =   189.8517
DATA['SAPT EXCH ENERGY']['S22by7-11-0.8'] =    70.9068
DATA['SAPT EXCH ENERGY']['S22by7-11-0.9'] =    25.8116
DATA['SAPT EXCH ENERGY']['S22by7-11-1.0'] =     9.1767
DATA['SAPT EXCH ENERGY']['S22by7-11-1.2'] =     1.1215
DATA['SAPT EXCH ENERGY']['S22by7-11-1.5'] =     0.0426
DATA['SAPT EXCH ENERGY']['S22by7-11-2.0'] =     0.0002
DATA['SAPT EXCH ENERGY']['S22by7-12-0.7'] =   211.9560
DATA['SAPT EXCH ENERGY']['S22by7-12-0.8'] =    78.5144
DATA['SAPT EXCH ENERGY']['S22by7-12-0.9'] =    28.3076
DATA['SAPT EXCH ENERGY']['S22by7-12-1.0'] =     9.9625
DATA['SAPT EXCH ENERGY']['S22by7-12-1.2'] =     1.1750
DATA['SAPT EXCH ENERGY']['S22by7-12-1.5'] =     0.0400
DATA['SAPT EXCH ENERGY']['S22by7-12-2.0'] =     0.0001
DATA['SAPT EXCH ENERGY']['S22by7-13-0.7'] =   283.0975
DATA['SAPT EXCH ENERGY']['S22by7-13-0.8'] =    98.7720
DATA['SAPT EXCH ENERGY']['S22by7-13-0.9'] =    33.2213
DATA['SAPT EXCH ENERGY']['S22by7-13-1.0'] =    10.8815
DATA['SAPT EXCH ENERGY']['S22by7-13-1.2'] =     1.1037
DATA['SAPT EXCH ENERGY']['S22by7-13-1.5'] =     0.0280
DATA['SAPT EXCH ENERGY']['S22by7-13-2.0'] =    -0.0001
DATA['SAPT EXCH ENERGY']['S22by7-14-0.7'] =   214.0774
DATA['SAPT EXCH ENERGY']['S22by7-14-0.8'] =    86.6796
DATA['SAPT EXCH ENERGY']['S22by7-14-0.9'] =    34.2899
DATA['SAPT EXCH ENERGY']['S22by7-14-1.0'] =    13.3016
DATA['SAPT EXCH ENERGY']['S22by7-14-1.2'] =     1.9405
DATA['SAPT EXCH ENERGY']['S22by7-14-1.5'] =     0.0975
DATA['SAPT EXCH ENERGY']['S22by7-14-2.0'] =     0.0003
DATA['SAPT EXCH ENERGY']['S22by7-15-0.7'] =   338.5822
DATA['SAPT EXCH ENERGY']['S22by7-15-0.8'] =   129.0756
DATA['SAPT EXCH ENERGY']['S22by7-15-0.9'] =    47.6112
DATA['SAPT EXCH ENERGY']['S22by7-15-1.0'] =    17.0933
DATA['SAPT EXCH ENERGY']['S22by7-15-1.2'] =     2.0764
DATA['SAPT EXCH ENERGY']['S22by7-15-1.5'] =     0.0706
DATA['SAPT EXCH ENERGY']['S22by7-15-2.0'] =     0.0001
DATA['SAPT EXCH ENERGY']['S22by7-16-0.7'] =    25.9486
DATA['SAPT EXCH ENERGY']['S22by7-16-0.8'] =    11.7870
DATA['SAPT EXCH ENERGY']['S22by7-16-0.9'] =     5.2675
DATA['SAPT EXCH ENERGY']['S22by7-16-1.0'] =     2.3246
DATA['SAPT EXCH ENERGY']['S22by7-16-1.2'] =     0.4434
DATA['SAPT EXCH ENERGY']['S22by7-16-1.5'] =     0.0352
DATA['SAPT EXCH ENERGY']['S22by7-16-2.0'] =     0.0005
DATA['SAPT EXCH ENERGY']['S22by7-17-0.7'] =    28.8925
DATA['SAPT EXCH ENERGY']['S22by7-17-0.8'] =    13.9581
DATA['SAPT EXCH ENERGY']['S22by7-17-0.9'] =     6.5986
DATA['SAPT EXCH ENERGY']['S22by7-17-1.0'] =     3.0621
DATA['SAPT EXCH ENERGY']['S22by7-17-1.2'] =     0.6268
DATA['SAPT EXCH ENERGY']['S22by7-17-1.5'] =     0.0510
DATA['SAPT EXCH ENERGY']['S22by7-17-2.0'] =     0.0008
DATA['SAPT EXCH ENERGY']['S22by7-18-0.7'] =    26.3818
DATA['SAPT EXCH ENERGY']['S22by7-18-0.8'] =    12.4175
DATA['SAPT EXCH ENERGY']['S22by7-18-0.9'] =     5.7075
DATA['SAPT EXCH ENERGY']['S22by7-18-1.0'] =     2.5749
DATA['SAPT EXCH ENERGY']['S22by7-18-1.2'] =     0.5031
DATA['SAPT EXCH ENERGY']['S22by7-18-1.5'] =     0.0388
DATA['SAPT EXCH ENERGY']['S22by7-18-2.0'] =     0.0005
DATA['SAPT EXCH ENERGY']['S22by7-19-0.7'] =    33.4802
DATA['SAPT EXCH ENERGY']['S22by7-19-0.8'] =    17.4681
DATA['SAPT EXCH ENERGY']['S22by7-19-0.9'] =     8.8476
DATA['SAPT EXCH ENERGY']['S22by7-19-1.0'] =     4.3692
DATA['SAPT EXCH ENERGY']['S22by7-19-1.2'] =     1.0104
DATA['SAPT EXCH ENERGY']['S22by7-19-1.5'] =     0.1034
DATA['SAPT EXCH ENERGY']['S22by7-19-2.0'] =     0.0018
DATA['SAPT EXCH ENERGY']['S22by7-20-0.7'] =    35.2484
DATA['SAPT EXCH ENERGY']['S22by7-20-0.8'] =    17.9713
DATA['SAPT EXCH ENERGY']['S22by7-20-0.9'] =     8.8941
DATA['SAPT EXCH ENERGY']['S22by7-20-1.0'] =     4.2914
DATA['SAPT EXCH ENERGY']['S22by7-20-1.2'] =     0.9465
DATA['SAPT EXCH ENERGY']['S22by7-20-1.5'] =     0.0888
DATA['SAPT EXCH ENERGY']['S22by7-20-2.0'] =     0.0014
DATA['SAPT EXCH ENERGY']['S22by7-21-0.7'] =  1501.3669
DATA['SAPT EXCH ENERGY']['S22by7-21-0.8'] =  1490.5265
DATA['SAPT EXCH ENERGY']['S22by7-21-0.9'] =    12.8764
DATA['SAPT EXCH ENERGY']['S22by7-21-1.0'] =     6.6375
DATA['SAPT EXCH ENERGY']['S22by7-21-1.2'] =     1.6900
DATA['SAPT EXCH ENERGY']['S22by7-21-1.5'] =     0.1999
DATA['SAPT EXCH ENERGY']['S22by7-21-2.0'] =     0.0048
DATA['SAPT EXCH ENERGY']['S22by7-22-0.7'] =    65.4982
DATA['SAPT EXCH ENERGY']['S22by7-22-0.8'] =    34.7280
DATA['SAPT EXCH ENERGY']['S22by7-22-0.9'] =    18.4734
DATA['SAPT EXCH ENERGY']['S22by7-22-1.0'] =     9.8875
DATA['SAPT EXCH ENERGY']['S22by7-22-1.2'] =     2.8938
DATA['SAPT EXCH ENERGY']['S22by7-22-1.5'] =     0.4846
DATA['SAPT EXCH ENERGY']['S22by7-22-2.0'] =     0.0266
DATA['SAPT IND ENERGY'] = {}
DATA['SAPT IND ENERGY']['S22by7-1-0.7'] =    -7.1135
DATA['SAPT IND ENERGY']['S22by7-1-0.8'] =    -3.5206
DATA['SAPT IND ENERGY']['S22by7-1-0.9'] =    -1.6683
DATA['SAPT IND ENERGY']['S22by7-1-1.0'] =    -0.8002
DATA['SAPT IND ENERGY']['S22by7-1-1.2'] =    -0.2094
DATA['SAPT IND ENERGY']['S22by7-1-1.5'] =    -0.0419
DATA['SAPT IND ENERGY']['S22by7-1-2.0'] =    -0.0065
DATA['SAPT IND ENERGY']['S22by7-2-0.7'] =   -12.7266
DATA['SAPT IND ENERGY']['S22by7-2-0.8'] =    -7.1119
DATA['SAPT IND ENERGY']['S22by7-2-0.9'] =    -3.8533
DATA['SAPT IND ENERGY']['S22by7-2-1.0'] =    -2.0897
DATA['SAPT IND ENERGY']['S22by7-2-1.2'] =    -0.6535
DATA['SAPT IND ENERGY']['S22by7-2-1.5'] =    -0.1415
DATA['SAPT IND ENERGY']['S22by7-2-2.0'] =    -0.0214
DATA['SAPT IND ENERGY']['S22by7-3-0.7'] =   -67.0024
DATA['SAPT IND ENERGY']['S22by7-3-0.8'] =   -43.1788
DATA['SAPT IND ENERGY']['S22by7-3-0.9'] =   -26.5766
DATA['SAPT IND ENERGY']['S22by7-3-1.0'] =   -16.1137
DATA['SAPT IND ENERGY']['S22by7-3-1.2'] =    -6.0022
DATA['SAPT IND ENERGY']['S22by7-3-1.5'] =    -1.5422
DATA['SAPT IND ENERGY']['S22by7-3-2.0'] =    -0.2361
DATA['SAPT IND ENERGY']['S22by7-4-0.7'] =  -753.2927
DATA['SAPT IND ENERGY']['S22by7-4-0.8'] =  -813.7292
DATA['SAPT IND ENERGY']['S22by7-4-0.9'] =   -16.3758
DATA['SAPT IND ENERGY']['S22by7-4-1.0'] =    -9.6651
DATA['SAPT IND ENERGY']['S22by7-4-1.2'] =    -3.5675
DATA['SAPT IND ENERGY']['S22by7-4-1.5'] =    -0.9729
DATA['SAPT IND ENERGY']['S22by7-4-2.0'] =    -0.1831
DATA['SAPT IND ENERGY']['S22by7-5-0.7'] =   -54.0349
DATA['SAPT IND ENERGY']['S22by7-5-0.8'] =   -33.5997
DATA['SAPT IND ENERGY']['S22by7-5-0.9'] =   -20.2525
DATA['SAPT IND ENERGY']['S22by7-5-1.0'] =   -12.2118
DATA['SAPT IND ENERGY']['S22by7-5-1.2'] =    -4.6969
DATA['SAPT IND ENERGY']['S22by7-5-1.5'] =    -1.3697
DATA['SAPT IND ENERGY']['S22by7-5-2.0'] =    -0.2904
DATA['SAPT IND ENERGY']['S22by7-6-0.7'] =   -45.6700
DATA['SAPT IND ENERGY']['S22by7-6-0.8'] =   -29.2539
DATA['SAPT IND ENERGY']['S22by7-6-0.9'] =   -17.9777
DATA['SAPT IND ENERGY']['S22by7-6-1.0'] =   -10.9344
DATA['SAPT IND ENERGY']['S22by7-6-1.2'] =    -4.1520
DATA['SAPT IND ENERGY']['S22by7-6-1.5'] =    -1.1270
DATA['SAPT IND ENERGY']['S22by7-6-2.0'] =    -0.2050
DATA['SAPT IND ENERGY']['S22by7-7-0.7'] =   -45.2036
DATA['SAPT IND ENERGY']['S22by7-7-0.8'] =   -28.3875
DATA['SAPT IND ENERGY']['S22by7-7-0.9'] =   -17.1590
DATA['SAPT IND ENERGY']['S22by7-7-1.0'] =   -10.2389
DATA['SAPT IND ENERGY']['S22by7-7-1.2'] =    -3.6743
DATA['SAPT IND ENERGY']['S22by7-7-1.5'] =    -0.8641
DATA['SAPT IND ENERGY']['S22by7-7-2.0'] =    -0.1088
DATA['SAPT IND ENERGY']['S22by7-8-0.7'] =    -0.9922
DATA['SAPT IND ENERGY']['S22by7-8-0.8'] =    -0.3599
DATA['SAPT IND ENERGY']['S22by7-8-0.9'] =    -0.1195
DATA['SAPT IND ENERGY']['S22by7-8-1.0'] =    -0.0398
DATA['SAPT IND ENERGY']['S22by7-8-1.2'] =    -0.0040
DATA['SAPT IND ENERGY']['S22by7-8-1.5'] =    -0.0008
DATA['SAPT IND ENERGY']['S22by7-8-2.0'] =    -0.0002
DATA['SAPT IND ENERGY']['S22by7-9-0.7'] =    -3.4873
DATA['SAPT IND ENERGY']['S22by7-9-0.8'] =    -1.6701
DATA['SAPT IND ENERGY']['S22by7-9-0.9'] =    -0.6348
DATA['SAPT IND ENERGY']['S22by7-9-1.0'] =    -0.2192
DATA['SAPT IND ENERGY']['S22by7-9-1.2'] =    -0.0229
DATA['SAPT IND ENERGY']['S22by7-9-1.5'] =    -0.0019
DATA['SAPT IND ENERGY']['S22by7-9-2.0'] =     0.0008
DATA['SAPT IND ENERGY']['S22by7-10-0.7'] =    -1.5068
DATA['SAPT IND ENERGY']['S22by7-10-0.8'] =    -0.9565
DATA['SAPT IND ENERGY']['S22by7-10-0.9'] =    -0.5349
DATA['SAPT IND ENERGY']['S22by7-10-1.0'] =    -0.2921
DATA['SAPT IND ENERGY']['S22by7-10-1.2'] =    -0.0930
DATA['SAPT IND ENERGY']['S22by7-10-1.5'] =    -0.0191
DATA['SAPT IND ENERGY']['S22by7-10-2.0'] =    -0.0018
DATA['SAPT IND ENERGY']['S22by7-11-0.7'] =     3.1321
DATA['SAPT IND ENERGY']['S22by7-11-0.8'] =    -3.1214
DATA['SAPT IND ENERGY']['S22by7-11-0.9'] =    -1.9723
DATA['SAPT IND ENERGY']['S22by7-11-1.0'] =    -0.9162
DATA['SAPT IND ENERGY']['S22by7-11-1.2'] =    -0.2117
DATA['SAPT IND ENERGY']['S22by7-11-1.5'] =    -0.0475
DATA['SAPT IND ENERGY']['S22by7-11-2.0'] =    -0.0039
DATA['SAPT IND ENERGY']['S22by7-12-0.7'] =    -0.8102
DATA['SAPT IND ENERGY']['S22by7-12-0.8'] =    -4.8258
DATA['SAPT IND ENERGY']['S22by7-12-0.9'] =    -2.4739
DATA['SAPT IND ENERGY']['S22by7-12-1.0'] =    -0.9892
DATA['SAPT IND ENERGY']['S22by7-12-1.2'] =    -0.1424
DATA['SAPT IND ENERGY']['S22by7-12-1.5'] =    -0.0134
DATA['SAPT IND ENERGY']['S22by7-12-2.0'] =    -0.0007
DATA['SAPT IND ENERGY']['S22by7-13-0.7'] =    -9.5727
DATA['SAPT IND ENERGY']['S22by7-13-0.8'] =    -7.6168
DATA['SAPT IND ENERGY']['S22by7-13-0.9'] =    -3.5284
DATA['SAPT IND ENERGY']['S22by7-13-1.0'] =    -1.5619
DATA['SAPT IND ENERGY']['S22by7-13-1.2'] =    -0.4171
DATA['SAPT IND ENERGY']['S22by7-13-1.5'] =    -0.1145
DATA['SAPT IND ENERGY']['S22by7-13-2.0'] =    -0.0272
DATA['SAPT IND ENERGY']['S22by7-14-0.7'] =     2.7917
DATA['SAPT IND ENERGY']['S22by7-14-0.8'] =    -3.9170
DATA['SAPT IND ENERGY']['S22by7-14-0.9'] =    -2.7807
DATA['SAPT IND ENERGY']['S22by7-14-1.0'] =    -1.4610
DATA['SAPT IND ENERGY']['S22by7-14-1.2'] =    -0.4060
DATA['SAPT IND ENERGY']['S22by7-14-1.5'] =    -0.1041
DATA['SAPT IND ENERGY']['S22by7-14-2.0'] =    -0.0184
DATA['SAPT IND ENERGY']['S22by7-15-0.7'] =   -10.1385
DATA['SAPT IND ENERGY']['S22by7-15-0.8'] =   -10.1703
DATA['SAPT IND ENERGY']['S22by7-15-0.9'] =    -5.1448
DATA['SAPT IND ENERGY']['S22by7-15-1.0'] =    -2.2901
DATA['SAPT IND ENERGY']['S22by7-15-1.2'] =    -0.5087
DATA['SAPT IND ENERGY']['S22by7-15-1.5'] =    -0.1117
DATA['SAPT IND ENERGY']['S22by7-15-2.0'] =    -0.0265
DATA['SAPT IND ENERGY']['S22by7-16-0.7'] =    -4.4074
DATA['SAPT IND ENERGY']['S22by7-16-0.8'] =    -2.2479
DATA['SAPT IND ENERGY']['S22by7-16-0.9'] =    -1.1119
DATA['SAPT IND ENERGY']['S22by7-16-1.0'] =    -0.5519
DATA['SAPT IND ENERGY']['S22by7-16-1.2'] =    -0.1453
DATA['SAPT IND ENERGY']['S22by7-16-1.5'] =    -0.0252
DATA['SAPT IND ENERGY']['S22by7-16-2.0'] =    -0.0022
DATA['SAPT IND ENERGY']['S22by7-17-0.7'] =    -4.9147
DATA['SAPT IND ENERGY']['S22by7-17-0.8'] =    -2.7811
DATA['SAPT IND ENERGY']['S22by7-17-0.9'] =    -1.5521
DATA['SAPT IND ENERGY']['S22by7-17-1.0'] =    -0.8866
DATA['SAPT IND ENERGY']['S22by7-17-1.2'] =    -0.3228
DATA['SAPT IND ENERGY']['S22by7-17-1.5'] =    -0.0897
DATA['SAPT IND ENERGY']['S22by7-17-2.0'] =    -0.0180
DATA['SAPT IND ENERGY']['S22by7-18-0.7'] =    -2.8378
DATA['SAPT IND ENERGY']['S22by7-18-0.8'] =    -1.5889
DATA['SAPT IND ENERGY']['S22by7-18-0.9'] =    -0.8549
DATA['SAPT IND ENERGY']['S22by7-18-1.0'] =    -0.4682
DATA['SAPT IND ENERGY']['S22by7-18-1.2'] =    -0.1591
DATA['SAPT IND ENERGY']['S22by7-18-1.5'] =    -0.0418
DATA['SAPT IND ENERGY']['S22by7-18-2.0'] =    -0.0088
DATA['SAPT IND ENERGY']['S22by7-19-0.7'] =    -6.9076
DATA['SAPT IND ENERGY']['S22by7-19-0.8'] =    -4.4948
DATA['SAPT IND ENERGY']['S22by7-19-0.9'] =    -2.8055
DATA['SAPT IND ENERGY']['S22by7-19-1.0'] =    -1.7444
DATA['SAPT IND ENERGY']['S22by7-19-1.2'] =    -0.7157
DATA['SAPT IND ENERGY']['S22by7-19-1.5'] =    -0.2333
DATA['SAPT IND ENERGY']['S22by7-19-2.0'] =    -0.0577
DATA['SAPT IND ENERGY']['S22by7-20-0.7'] =    -2.7243
DATA['SAPT IND ENERGY']['S22by7-20-0.8'] =    -1.8359
DATA['SAPT IND ENERGY']['S22by7-20-0.9'] =    -1.0936
DATA['SAPT IND ENERGY']['S22by7-20-1.0'] =    -0.6284
DATA['SAPT IND ENERGY']['S22by7-20-1.2'] =    -0.2155
DATA['SAPT IND ENERGY']['S22by7-20-1.5'] =    -0.0533
DATA['SAPT IND ENERGY']['S22by7-20-2.0'] =    -0.0118
DATA['SAPT IND ENERGY']['S22by7-21-0.7'] =    44.4600
DATA['SAPT IND ENERGY']['S22by7-21-0.8'] =    51.1890
DATA['SAPT IND ENERGY']['S22by7-21-0.9'] =    -2.6288
DATA['SAPT IND ENERGY']['S22by7-21-1.0'] =    -1.6770
DATA['SAPT IND ENERGY']['S22by7-21-1.2'] =    -0.7139
DATA['SAPT IND ENERGY']['S22by7-21-1.5'] =    -0.2327
DATA['SAPT IND ENERGY']['S22by7-21-2.0'] =    -0.0567
DATA['SAPT IND ENERGY']['S22by7-22-0.7'] =   -14.9227
DATA['SAPT IND ENERGY']['S22by7-22-0.8'] =    -8.5504
DATA['SAPT IND ENERGY']['S22by7-22-0.9'] =    -4.8291
DATA['SAPT IND ENERGY']['S22by7-22-1.0'] =    -2.7692
DATA['SAPT IND ENERGY']['S22by7-22-1.2'] =    -1.0007
DATA['SAPT IND ENERGY']['S22by7-22-1.5'] =    -0.2824
DATA['SAPT IND ENERGY']['S22by7-22-2.0'] =    -0.0626
DATA['SAPT DISP ENERGY'] = {}
DATA['SAPT DISP ENERGY']['S22by7-1-0.7'] =    -6.4602
DATA['SAPT DISP ENERGY']['S22by7-1-0.8'] =    -3.5943
DATA['SAPT DISP ENERGY']['S22by7-1-0.9'] =    -2.0397
DATA['SAPT DISP ENERGY']['S22by7-1-1.0'] =    -1.1882
DATA['SAPT DISP ENERGY']['S22by7-1-1.2'] =    -0.4362
DATA['SAPT DISP ENERGY']['S22by7-1-1.5'] =    -0.1209
DATA['SAPT DISP ENERGY']['S22by7-1-2.0'] =    -0.0238
DATA['SAPT DISP ENERGY']['S22by7-2-0.7'] =    -5.2583
DATA['SAPT DISP ENERGY']['S22by7-2-0.8'] =    -3.2507
DATA['SAPT DISP ENERGY']['S22by7-2-0.9'] =    -2.0149
DATA['SAPT DISP ENERGY']['S22by7-2-1.0'] =    -1.2650
DATA['SAPT DISP ENERGY']['S22by7-2-1.2'] =    -0.5298
DATA['SAPT DISP ENERGY']['S22by7-2-1.5'] =    -0.1647
DATA['SAPT DISP ENERGY']['S22by7-2-2.0'] =    -0.0339
DATA['SAPT DISP ENERGY']['S22by7-3-0.7'] =   -18.9855
DATA['SAPT DISP ENERGY']['S22by7-3-0.8'] =   -12.9622
DATA['SAPT DISP ENERGY']['S22by7-3-0.9'] =    -8.8676
DATA['SAPT DISP ENERGY']['S22by7-3-1.0'] =    -6.1040
DATA['SAPT DISP ENERGY']['S22by7-3-1.2'] =    -3.0050
DATA['SAPT DISP ENERGY']['S22by7-3-1.5'] =    -1.1629
DATA['SAPT DISP ENERGY']['S22by7-3-2.0'] =    -0.3060
DATA['SAPT DISP ENERGY']['S22by7-4-0.7'] =  -281.8221
DATA['SAPT DISP ENERGY']['S22by7-4-0.8'] =  -313.9454
DATA['SAPT DISP ENERGY']['S22by7-4-0.9'] =    -7.2690
DATA['SAPT DISP ENERGY']['S22by7-4-1.0'] =    -4.8909
DATA['SAPT DISP ENERGY']['S22by7-4-1.2'] =    -2.3552
DATA['SAPT DISP ENERGY']['S22by7-4-1.5'] =    -0.8986
DATA['SAPT DISP ENERGY']['S22by7-4-2.0'] =    -0.2367
DATA['SAPT DISP ENERGY']['S22by7-5-0.7'] =   -19.6378
DATA['SAPT DISP ENERGY']['S22by7-5-0.8'] =   -13.2683
DATA['SAPT DISP ENERGY']['S22by7-5-0.9'] =    -9.0714
DATA['SAPT DISP ENERGY']['S22by7-5-1.0'] =    -6.3054
DATA['SAPT DISP ENERGY']['S22by7-5-1.2'] =    -3.2431
DATA['SAPT DISP ENERGY']['S22by7-5-1.5'] =    -1.3697
DATA['SAPT DISP ENERGY']['S22by7-5-2.0'] =    -0.4277
DATA['SAPT DISP ENERGY']['S22by7-6-0.7'] =   -20.8538
DATA['SAPT DISP ENERGY']['S22by7-6-0.8'] =   -14.4083
DATA['SAPT DISP ENERGY']['S22by7-6-0.9'] =   -10.0590
DATA['SAPT DISP ENERGY']['S22by7-6-1.0'] =    -7.1240
DATA['SAPT DISP ENERGY']['S22by7-6-1.2'] =    -3.7639
DATA['SAPT DISP ENERGY']['S22by7-6-1.5'] =    -1.6185
DATA['SAPT DISP ENERGY']['S22by7-6-2.0'] =    -0.5040
DATA['SAPT DISP ENERGY']['S22by7-7-0.7'] =   -21.1227
DATA['SAPT DISP ENERGY']['S22by7-7-0.8'] =   -14.5745
DATA['SAPT DISP ENERGY']['S22by7-7-0.9'] =   -10.1938
DATA['SAPT DISP ENERGY']['S22by7-7-1.0'] =    -7.2460
DATA['SAPT DISP ENERGY']['S22by7-7-1.2'] =    -3.8659
DATA['SAPT DISP ENERGY']['S22by7-7-1.5'] =    -1.6938
DATA['SAPT DISP ENERGY']['S22by7-7-2.0'] =    -0.5436
DATA['SAPT DISP ENERGY']['S22by7-8-0.7'] =    -4.7430
DATA['SAPT DISP ENERGY']['S22by7-8-0.8'] =    -2.1629
DATA['SAPT DISP ENERGY']['S22by7-8-0.9'] =    -1.0480
DATA['SAPT DISP ENERGY']['S22by7-8-1.0'] =    -0.5362
DATA['SAPT DISP ENERGY']['S22by7-8-1.2'] =    -0.1663
DATA['SAPT DISP ENERGY']['S22by7-8-1.5'] =    -0.0412
DATA['SAPT DISP ENERGY']['S22by7-8-2.0'] =    -0.0071
DATA['SAPT DISP ENERGY']['S22by7-9-0.7'] =   -12.2603
DATA['SAPT DISP ENERGY']['S22by7-9-0.8'] =    -6.1967
DATA['SAPT DISP ENERGY']['S22by7-9-0.9'] =    -3.1900
DATA['SAPT DISP ENERGY']['S22by7-9-1.0'] =    -1.6871
DATA['SAPT DISP ENERGY']['S22by7-9-1.2'] =    -0.5234
DATA['SAPT DISP ENERGY']['S22by7-9-1.5'] =    -0.1246
DATA['SAPT DISP ENERGY']['S22by7-9-2.0'] =    -0.0207
DATA['SAPT DISP ENERGY']['S22by7-10-0.7'] =    -7.7858
DATA['SAPT DISP ENERGY']['S22by7-10-0.8'] =    -5.1093
DATA['SAPT DISP ENERGY']['S22by7-10-0.9'] =    -3.3622
DATA['SAPT DISP ENERGY']['S22by7-10-1.0'] =    -2.2308
DATA['SAPT DISP ENERGY']['S22by7-10-1.2'] =    -1.0220
DATA['SAPT DISP ENERGY']['S22by7-10-1.5'] =    -0.3600
DATA['SAPT DISP ENERGY']['S22by7-10-2.0'] =    -0.0855
DATA['SAPT DISP ENERGY']['S22by7-11-0.7'] =   -41.0665
DATA['SAPT DISP ENERGY']['S22by7-11-0.8'] =   -23.7973
DATA['SAPT DISP ENERGY']['S22by7-11-0.9'] =   -13.7895
DATA['SAPT DISP ENERGY']['S22by7-11-1.0'] =    -8.0975
DATA['SAPT DISP ENERGY']['S22by7-11-1.2'] =    -2.9657
DATA['SAPT DISP ENERGY']['S22by7-11-1.5'] =    -0.7919
DATA['SAPT DISP ENERGY']['S22by7-11-2.0'] =    -0.1376
DATA['SAPT DISP ENERGY']['S22by7-12-0.7'] =   -42.7538
DATA['SAPT DISP ENERGY']['S22by7-12-0.8'] =   -24.7652
DATA['SAPT DISP ENERGY']['S22by7-12-0.9'] =   -14.3893
DATA['SAPT DISP ENERGY']['S22by7-12-1.0'] =    -8.4782
DATA['SAPT DISP ENERGY']['S22by7-12-1.2'] =    -3.1367
DATA['SAPT DISP ENERGY']['S22by7-12-1.5'] =    -0.8586
DATA['SAPT DISP ENERGY']['S22by7-12-2.0'] =    -0.1547
DATA['SAPT DISP ENERGY']['S22by7-13-0.7'] =   -48.0247
DATA['SAPT DISP ENERGY']['S22by7-13-0.8'] =   -27.4589
DATA['SAPT DISP ENERGY']['S22by7-13-0.9'] =   -15.8291
DATA['SAPT DISP ENERGY']['S22by7-13-1.0'] =    -9.3275
DATA['SAPT DISP ENERGY']['S22by7-13-1.2'] =    -3.5187
DATA['SAPT DISP ENERGY']['S22by7-13-1.5'] =    -1.0106
DATA['SAPT DISP ENERGY']['S22by7-13-2.0'] =    -0.1964
DATA['SAPT DISP ENERGY']['S22by7-14-0.7'] =   -48.8384
DATA['SAPT DISP ENERGY']['S22by7-14-0.8'] =   -30.5780
DATA['SAPT DISP ENERGY']['S22by7-14-0.9'] =   -19.1726
DATA['SAPT DISP ENERGY']['S22by7-14-1.0'] =   -12.1545
DATA['SAPT DISP ENERGY']['S22by7-14-1.2'] =    -5.1218
DATA['SAPT DISP ENERGY']['S22by7-14-1.5'] =    -1.6123
DATA['SAPT DISP ENERGY']['S22by7-14-2.0'] =    -0.3325
DATA['SAPT DISP ENERGY']['S22by7-15-0.7'] =   -68.1930
DATA['SAPT DISP ENERGY']['S22by7-15-0.8'] =   -41.2848
DATA['SAPT DISP ENERGY']['S22by7-15-0.9'] =   -25.0048
DATA['SAPT DISP ENERGY']['S22by7-15-1.0'] =   -15.3543
DATA['SAPT DISP ENERGY']['S22by7-15-1.2'] =    -6.1536
DATA['SAPT DISP ENERGY']['S22by7-15-1.5'] =    -1.8671
DATA['SAPT DISP ENERGY']['S22by7-15-2.0'] =    -0.3835
DATA['SAPT DISP ENERGY']['S22by7-16-0.7'] =    -4.5434
DATA['SAPT DISP ENERGY']['S22by7-16-0.8'] =    -2.6986
DATA['SAPT DISP ENERGY']['S22by7-16-0.9'] =    -1.6573
DATA['SAPT DISP ENERGY']['S22by7-16-1.0'] =    -1.0476
DATA['SAPT DISP ENERGY']['S22by7-16-1.2'] =    -0.4444
DATA['SAPT DISP ENERGY']['S22by7-16-1.5'] =    -0.1424
DATA['SAPT DISP ENERGY']['S22by7-16-2.0'] =    -0.0325
DATA['SAPT DISP ENERGY']['S22by7-17-0.7'] =    -7.2155
DATA['SAPT DISP ENERGY']['S22by7-17-0.8'] =    -4.7128
DATA['SAPT DISP ENERGY']['S22by7-17-0.9'] =    -3.0892
DATA['SAPT DISP ENERGY']['S22by7-17-1.0'] =    -2.0443
DATA['SAPT DISP ENERGY']['S22by7-17-1.2'] =    -0.9326
DATA['SAPT DISP ENERGY']['S22by7-17-1.5'] =    -0.3246
DATA['SAPT DISP ENERGY']['S22by7-17-2.0'] =    -0.0756
DATA['SAPT DISP ENERGY']['S22by7-18-0.7'] =    -7.5472
DATA['SAPT DISP ENERGY']['S22by7-18-0.8'] =    -4.9200
DATA['SAPT DISP ENERGY']['S22by7-18-0.9'] =    -3.2261
DATA['SAPT DISP ENERGY']['S22by7-18-1.0'] =    -2.1376
DATA['SAPT DISP ENERGY']['S22by7-18-1.2'] =    -0.9793
DATA['SAPT DISP ENERGY']['S22by7-18-1.5'] =    -0.3450
DATA['SAPT DISP ENERGY']['S22by7-18-2.0'] =    -0.0818
DATA['SAPT DISP ENERGY']['S22by7-19-0.7'] =    -9.6485
DATA['SAPT DISP ENERGY']['S22by7-19-0.8'] =    -6.6019
DATA['SAPT DISP ENERGY']['S22by7-19-0.9'] =    -4.5198
DATA['SAPT DISP ENERGY']['S22by7-19-1.0'] =    -3.1114
DATA['SAPT DISP ENERGY']['S22by7-19-1.2'] =    -1.5181
DATA['SAPT DISP ENERGY']['S22by7-19-1.5'] =    -0.5735
DATA['SAPT DISP ENERGY']['S22by7-19-2.0'] =    -0.1507
DATA['SAPT DISP ENERGY']['S22by7-20-0.7'] =   -12.7067
DATA['SAPT DISP ENERGY']['S22by7-20-0.8'] =    -8.8065
DATA['SAPT DISP ENERGY']['S22by7-20-0.9'] =    -6.1170
DATA['SAPT DISP ENERGY']['S22by7-20-1.0'] =    -4.2769
DATA['SAPT DISP ENERGY']['S22by7-20-1.2'] =    -2.1575
DATA['SAPT DISP ENERGY']['S22by7-20-1.5'] =    -0.8582
DATA['SAPT DISP ENERGY']['S22by7-20-2.0'] =    -0.2419
DATA['SAPT DISP ENERGY']['S22by7-21-0.7'] =  -230.5171
DATA['SAPT DISP ENERGY']['S22by7-21-0.8'] =  -222.3505
DATA['SAPT DISP ENERGY']['S22by7-21-0.9'] =    -8.1660
DATA['SAPT DISP ENERGY']['S22by7-21-1.0'] =    -5.8902
DATA['SAPT DISP ENERGY']['S22by7-21-1.2'] =    -3.1564
DATA['SAPT DISP ENERGY']['S22by7-21-1.5'] =    -1.3517
DATA['SAPT DISP ENERGY']['S22by7-21-2.0'] =    -0.4121
DATA['SAPT DISP ENERGY']['S22by7-22-0.7'] =   -12.8597
DATA['SAPT DISP ENERGY']['S22by7-22-0.8'] =    -9.3135
DATA['SAPT DISP ENERGY']['S22by7-22-0.9'] =    -6.8389
DATA['SAPT DISP ENERGY']['S22by7-22-1.0'] =    -5.1004
DATA['SAPT DISP ENERGY']['S22by7-22-1.2'] =    -2.9592
DATA['SAPT DISP ENERGY']['S22by7-22-1.5'] =    -1.4183
DATA['SAPT DISP ENERGY']['S22by7-22-2.0'] =    -0.4948
