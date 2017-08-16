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
| Database of <description of members and reference energy type>.
| Geometries from <Reference>.
| Reference interaction energies from <Reference>.


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

# <<< UBQ Database Module >>>
dbse = 'UBQ'

# <<< Database Members >>>
HRXN = ['ala28-asp32', 'ala46-lys48', 'arg42-val70', 'asn25-lys29', 'asp21-leu56',
        'gln2-glu64-small', 'gln31-gly35', 'gln40-arg72', 'gln41-arg72', 'gln62-ser65',
        'glu16-lys29', 'glu18-asp21-C', 'glu18-asp21-N', 'glu24-ala28', 'glu24-asp52',
        'glu51-arg54', 'glu51-tyr59', 'ile13-leu15', 'ile13-lys33', 'ile23-arg54',
        'ile23-leu50', 'ile23-leu56', 'ile23-lys27', 'ile3-leu15-big', 'ile3-leu15-part1',
        'ile3-leu15-part2', 'ile3-leu15', 'ile3-leu67', 'ile3-val17', 'ile30-glu34',
        'ile30-leu43', 'ile36-gln41', 'ile36-leu69', 'ile36-leu71', 'ile44-gly47',
        'ile44-his68-big', 'ile44-his68-small', 'ile44-his68', 'ile61-leu67', 'leu15-ile30',
        'leu15-val26', 'leu43-leu50', 'leu43-leu69', 'leu50-ile61', 'leu50-tyr59',
        'leu56-ile61np', 'leu56-ile61p', 'leu56-tyr59', 'leu69-leu71', 'leu71-leu73',
        'lys11-glu34', 'lys27-asp52', 'lys27-gln31', 'lys27-leu43-nonh3', 'lys27-leu43',
        'lys29-lys33', 'lys6-thr66', 'met1-ile3', 'met1-lys63', 'met1-val17-bi',
        'met1-val17-mono', 'met1-val17', 'phe4-leu67', 'phe4-ser65', 'phe4-thr12',
        'phe4-thr14', 'phe45-ala46', 'phe45-ile61', 'phe45-leu67', 'phe45-lys48',
        'pro19-ser57', 'pro37-gln40', 'ser57-asn60', 'thr14-lys33', 'thr22-asn25',
        'thr22-thr55-big', 'thr22-thr55-small', 'thr22-val26', 'thr55-asp58-backbone', 'thr55-asp58-sidechain-small',
        'thr55-asp58-sidechain', 'thr7-gly10', 'thr7-ile13', 'thr7-lys11', 'thr7-thr9-v2',
        'thr7-thr9', 'tyr59-ile61', 'val17-leu56', 'val26-ile30', 'val26-leu43',
        'val26-leu56', 'val5-ile13', 'val5-leu15', 'val5-leu69']
#HRXN_SM = []
#HRXN_LG = []

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
BIND = {}
# Original publication
# Current revision
BIND_UBQREF = {}
BIND_UBQREF['%s-%s'            % (dbse, 'ala28-asp32'                )] =   -7.775
BIND_UBQREF['%s-%s'            % (dbse, 'ala46-lys48'                )] =  -25.767
BIND_UBQREF['%s-%s'            % (dbse, 'arg42-val70'                )] =  -12.312
BIND_UBQREF['%s-%s'            % (dbse, 'asn25-lys29'                )] =   -7.189
BIND_UBQREF['%s-%s'            % (dbse, 'asp21-leu56'                )] =   -7.573
BIND_UBQREF['%s-%s'            % (dbse, 'gln2-glu64-small'           )] =   -7.047
BIND_UBQREF['%s-%s'            % (dbse, 'gln31-gly35'                )] =   -6.530
BIND_UBQREF['%s-%s'            % (dbse, 'gln40-arg72'                )] =   -6.843
BIND_UBQREF['%s-%s'            % (dbse, 'gln41-arg72'                )] =  -12.327
BIND_UBQREF['%s-%s'            % (dbse, 'gln62-ser65'                )] =  -11.211
BIND_UBQREF['%s-%s'            % (dbse, 'glu16-lys29'                )] =  -22.047
BIND_UBQREF['%s-%s'            % (dbse, 'glu18-asp21-C'              )] =   -6.026
BIND_UBQREF['%s-%s'            % (dbse, 'glu18-asp21-N'              )] =  -15.748
BIND_UBQREF['%s-%s'            % (dbse, 'glu24-ala28'                )] =   -6.496
BIND_UBQREF['%s-%s'            % (dbse, 'glu24-asp52'                )] =   -7.775
BIND_UBQREF['%s-%s'            % (dbse, 'glu51-arg54'                )] =  -51.003
BIND_UBQREF['%s-%s'            % (dbse, 'glu51-tyr59'                )] =   -4.517
BIND_UBQREF['%s-%s'            % (dbse, 'ile13-leu15'                )] =   -1.037
BIND_UBQREF['%s-%s'            % (dbse, 'ile13-lys33'                )] =   -1.993
BIND_UBQREF['%s-%s'            % (dbse, 'ile23-arg54'                )] =   -7.280
BIND_UBQREF['%s-%s'            % (dbse, 'ile23-leu50'                )] =   -1.584
BIND_UBQREF['%s-%s'            % (dbse, 'ile23-leu56'                )] =   -0.380
BIND_UBQREF['%s-%s'            % (dbse, 'ile23-lys27'                )] =   -2.624
BIND_UBQREF['%s-%s'            % (dbse, 'ile3-leu15-big'             )] =  -10.054
BIND_UBQREF['%s-%s'            % (dbse, 'ile3-leu15-part1'           )] =   -6.974
BIND_UBQREF['%s-%s'            % (dbse, 'ile3-leu15-part2'           )] =   -6.993
BIND_UBQREF['%s-%s'            % (dbse, 'ile3-leu15'                 )] =   -0.948
BIND_UBQREF['%s-%s'            % (dbse, 'ile3-leu67'                 )] =   -1.142
BIND_UBQREF['%s-%s'            % (dbse, 'ile3-val17'                 )] =   -1.471
BIND_UBQREF['%s-%s'            % (dbse, 'ile30-glu34'                )] =   -7.522
BIND_UBQREF['%s-%s'            % (dbse, 'ile30-leu43'                )] =   -0.801
BIND_UBQREF['%s-%s'            % (dbse, 'ile36-gln41'                )] =   -6.429
BIND_UBQREF['%s-%s'            % (dbse, 'ile36-leu69'                )] =   -0.667
BIND_UBQREF['%s-%s'            % (dbse, 'ile36-leu71'                )] =   -1.352
BIND_UBQREF['%s-%s'            % (dbse, 'ile44-gly47'                )] =   -1.525
BIND_UBQREF['%s-%s'            % (dbse, 'ile44-his68-big'            )] =   -0.229
BIND_UBQREF['%s-%s'            % (dbse, 'ile44-his68-small'          )] =    0.213
BIND_UBQREF['%s-%s'            % (dbse, 'ile44-his68'                )] =  -14.224
BIND_UBQREF['%s-%s'            % (dbse, 'ile61-leu67'                )] =   -0.997
BIND_UBQREF['%s-%s'            % (dbse, 'leu15-ile30'                )] =   -0.864
BIND_UBQREF['%s-%s'            % (dbse, 'leu15-val26'                )] =   -0.422
BIND_UBQREF['%s-%s'            % (dbse, 'leu43-leu50'                )] =   -7.446
BIND_UBQREF['%s-%s'            % (dbse, 'leu43-leu69'                )] =   -1.101
BIND_UBQREF['%s-%s'            % (dbse, 'leu50-ile61'                )] =   -0.648
BIND_UBQREF['%s-%s'            % (dbse, 'leu50-tyr59'                )] =   -2.350
BIND_UBQREF['%s-%s'            % (dbse, 'leu56-ile61np'              )] =   -1.277
BIND_UBQREF['%s-%s'            % (dbse, 'leu56-ile61p'               )] =   -5.159
BIND_UBQREF['%s-%s'            % (dbse, 'leu56-tyr59'                )] =   -4.417
BIND_UBQREF['%s-%s'            % (dbse, 'leu69-leu71'                )] =   -1.160
BIND_UBQREF['%s-%s'            % (dbse, 'leu71-leu73'                )] =   -0.865
BIND_UBQREF['%s-%s'            % (dbse, 'lys11-glu34'                )] = -101.500
BIND_UBQREF['%s-%s'            % (dbse, 'lys27-asp52'                )] = -110.569
BIND_UBQREF['%s-%s'            % (dbse, 'lys27-gln31'                )] =   -6.660
BIND_UBQREF['%s-%s'            % (dbse, 'lys27-leu43-nonh3'          )] =   -1.356
BIND_UBQREF['%s-%s'            % (dbse, 'lys27-leu43'                )] =   -2.139
BIND_UBQREF['%s-%s'            % (dbse, 'lys29-lys33'                )] =   -6.710
BIND_UBQREF['%s-%s'            % (dbse, 'lys6-thr66'                 )] =   -0.685
BIND_UBQREF['%s-%s'            % (dbse, 'met1-ile3'                  )] =   -0.957
BIND_UBQREF['%s-%s'            % (dbse, 'met1-lys63'                 )] =   -2.995
BIND_UBQREF['%s-%s'            % (dbse, 'met1-val17-bi'              )] =  -21.054
BIND_UBQREF['%s-%s'            % (dbse, 'met1-val17-mono'            )] =   -7.047
BIND_UBQREF['%s-%s'            % (dbse, 'met1-val17'                 )] =   -0.754
BIND_UBQREF['%s-%s'            % (dbse, 'phe4-leu67'                 )] =   -7.379
BIND_UBQREF['%s-%s'            % (dbse, 'phe4-ser65'                 )] =   -7.558
BIND_UBQREF['%s-%s'            % (dbse, 'phe4-thr12'                 )] =   -0.815
BIND_UBQREF['%s-%s'            % (dbse, 'phe4-thr14'                 )] =   -0.519
BIND_UBQREF['%s-%s'            % (dbse, 'phe45-ala46'                )] =   -0.865
BIND_UBQREF['%s-%s'            % (dbse, 'phe45-ile61'                )] =   -2.201
BIND_UBQREF['%s-%s'            % (dbse, 'phe45-leu67'                )] =   -1.034
BIND_UBQREF['%s-%s'            % (dbse, 'phe45-lys48'                )] =  -11.022
BIND_UBQREF['%s-%s'            % (dbse, 'pro19-ser57'                )] =   -7.002
BIND_UBQREF['%s-%s'            % (dbse, 'pro37-gln40'                )] =   -6.569
BIND_UBQREF['%s-%s'            % (dbse, 'ser57-asn60'                )] =   -6.898
BIND_UBQREF['%s-%s'            % (dbse, 'thr14-lys33'                )] =  -16.355
BIND_UBQREF['%s-%s'            % (dbse, 'thr22-asn25'                )] =   -5.169
BIND_UBQREF['%s-%s'            % (dbse, 'thr22-thr55-big'            )] =  -22.482
BIND_UBQREF['%s-%s'            % (dbse, 'thr22-thr55-small'          )] =   -1.031
BIND_UBQREF['%s-%s'            % (dbse, 'thr22-val26'                )] =   -6.544
BIND_UBQREF['%s-%s'            % (dbse, 'thr55-asp58-backbone'       )] =   -4.796
BIND_UBQREF['%s-%s'            % (dbse, 'thr55-asp58-sidechain-small')] =  -19.806
BIND_UBQREF['%s-%s'            % (dbse, 'thr55-asp58-sidechain'      )] =  -15.470
BIND_UBQREF['%s-%s'            % (dbse, 'thr7-gly10'                 )] =   -5.230
BIND_UBQREF['%s-%s'            % (dbse, 'thr7-ile13'                 )] =   -0.486
BIND_UBQREF['%s-%s'            % (dbse, 'thr7-lys11'                 )] =   -7.229
BIND_UBQREF['%s-%s'            % (dbse, 'thr7-thr9-v2'               )] =   -1.281
BIND_UBQREF['%s-%s'            % (dbse, 'thr7-thr9'                  )] =   -3.895
BIND_UBQREF['%s-%s'            % (dbse, 'tyr59-ile61'                )] =   -0.938
BIND_UBQREF['%s-%s'            % (dbse, 'val17-leu56'                )] =   -1.033
BIND_UBQREF['%s-%s'            % (dbse, 'val26-ile30'                )] =   -6.728
BIND_UBQREF['%s-%s'            % (dbse, 'val26-leu43'                )] =   -0.635
BIND_UBQREF['%s-%s'            % (dbse, 'val26-leu56'                )] =   -0.575
BIND_UBQREF['%s-%s'            % (dbse, 'val5-ile13'                 )] =  -14.475
BIND_UBQREF['%s-%s'            % (dbse, 'val5-leu15'                 )] =   -0.128
BIND_UBQREF['%s-%s'            % (dbse, 'val5-leu69'                 )] =   -0.685
# Set default
BIND = BIND_UBQREF
# Reference information
BINDINFO_UBQREF = {}
for rxn in HRXN:
    if rxn in ['glu51-tyr59', 'ile44-his68-big', 'leu50-tyr59', 'phe4-thr14', 'phe45-ala46', 'phe45-ile61', 'phe45-leu67', 'tyr59-ile61']:
        BINDINFO_UBQREF['%s-%s' % (dbse, rxn)] = {'citation': '1ubq', 'method': 'CCSDT', 'mode': 'CP', 'basis': 'atqzhadtz'}
    else:
        BINDINFO_UBQREF['%s-%s' % (dbse, rxn)] = {'citation': '1ubq', 'method': 'MP2', 'mode': 'CP', 'basis': 'atqz'}

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, 'ala28-asp32'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ala28-asp32'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ala28-asp32'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ala28-asp32'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ala28-asp32'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ala28-asp32'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ala46-lys48'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ala46-lys48'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ala46-lys48'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ala46-lys48'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ala46-lys48'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ala46-lys48'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'arg42-val70'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'arg42-val70'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'arg42-val70'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'arg42-val70'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'arg42-val70'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'arg42-val70'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'asn25-lys29'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'asn25-lys29'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'asn25-lys29'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'asn25-lys29'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'asn25-lys29'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'asn25-lys29'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'asp21-leu56'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'asp21-leu56'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'asp21-leu56'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'asp21-leu56'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'asp21-leu56'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'asp21-leu56'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'gln2-glu64-small'      )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'gln2-glu64-small'      )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'gln2-glu64-small'      )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'gln2-glu64-small'      )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'gln2-glu64-small'      )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'gln2-glu64-small'      )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'gln31-gly35'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'gln31-gly35'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'gln31-gly35'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'gln31-gly35'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'gln31-gly35'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'gln31-gly35'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'gln40-arg72'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'gln40-arg72'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'gln40-arg72'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'gln40-arg72'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'gln40-arg72'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'gln40-arg72'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'gln41-arg72'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'gln41-arg72'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'gln41-arg72'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'gln41-arg72'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'gln41-arg72'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'gln41-arg72'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'gln62-ser65'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'gln62-ser65'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'gln62-ser65'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'gln62-ser65'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'gln62-ser65'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'gln62-ser65'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'glu16-lys29'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'glu16-lys29'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'glu16-lys29'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'glu16-lys29'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'glu16-lys29'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'glu16-lys29'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'glu18-asp21-C'         )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'glu18-asp21-C'         )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'glu18-asp21-C'         )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'glu18-asp21-C'         )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'glu18-asp21-C'         )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'glu18-asp21-C'         )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'glu18-asp21-N'         )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'glu18-asp21-N'         )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'glu18-asp21-N'         )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'glu18-asp21-N'         )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'glu18-asp21-N'         )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'glu18-asp21-N'         )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'glu24-ala28'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'glu24-ala28'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'glu24-ala28'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'glu24-ala28'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'glu24-ala28'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'glu24-ala28'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'glu24-asp52'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'glu24-asp52'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'glu24-asp52'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'glu24-asp52'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'glu24-asp52'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'glu24-asp52'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'glu51-arg54'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'glu51-arg54'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'glu51-arg54'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'glu51-arg54'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'glu51-arg54'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'glu51-arg54'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'glu51-tyr59'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'glu51-tyr59'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'glu51-tyr59'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'glu51-tyr59'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'glu51-tyr59'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'glu51-tyr59'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile13-leu15'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile13-leu15'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile13-leu15'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile13-leu15'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile13-leu15'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile13-leu15'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile13-lys33'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile13-lys33'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile13-lys33'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile13-lys33'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile13-lys33'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile13-lys33'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile23-arg54'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile23-arg54'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile23-arg54'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile23-arg54'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile23-arg54'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile23-arg54'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile23-leu50'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile23-leu50'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile23-leu50'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile23-leu50'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile23-leu50'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile23-leu50'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile23-leu56'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile23-leu56'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile23-leu56'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile23-leu56'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile23-leu56'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile23-leu56'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile23-lys27'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile23-lys27'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile23-lys27'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile23-lys27'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile23-lys27'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile23-lys27'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile3-leu15-big'        )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile3-leu15-big'        )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile3-leu15-big'        )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile3-leu15-big'        )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile3-leu15-big'        )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile3-leu15-big'        )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile3-leu15-part1'      )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile3-leu15-part1'      )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile3-leu15-part1'      )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile3-leu15-part1'      )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile3-leu15-part1'      )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile3-leu15-part1'      )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile3-leu15-part2'      )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile3-leu15-part2'      )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile3-leu15-part2'      )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile3-leu15-part2'      )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile3-leu15-part2'      )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile3-leu15-part2'      )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile3-leu15'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile3-leu15'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile3-leu15'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile3-leu15'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile3-leu15'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile3-leu15'            )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile3-leu67'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile3-leu67'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile3-leu67'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile3-leu67'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile3-leu67'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile3-leu67'            )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile3-val17'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile3-val17'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile3-val17'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile3-val17'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile3-val17'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile3-val17'            )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile30-glu34'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile30-glu34'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile30-glu34'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile30-glu34'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile30-glu34'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile30-glu34'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile30-leu43'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile30-leu43'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile30-leu43'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile30-leu43'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile30-leu43'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile30-leu43'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile36-gln41'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile36-gln41'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile36-gln41'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile36-gln41'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile36-gln41'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile36-gln41'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile36-leu69'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile36-leu69'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile36-leu69'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile36-leu69'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile36-leu69'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile36-leu69'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile36-leu71'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile36-leu71'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile36-leu71'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile36-leu71'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile36-leu71'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile36-leu71'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile44-gly47'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile44-gly47'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile44-gly47'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile44-gly47'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile44-gly47'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile44-gly47'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile44-his68-big'       )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile44-his68-big'       )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile44-his68-big'       )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile44-his68-big'       )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile44-his68-big'       )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile44-his68-big'       )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile44-his68-small'     )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile44-his68-small'     )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile44-his68-small'     )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile44-his68-small'     )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile44-his68-small'     )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile44-his68-small'     )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile44-his68'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile44-his68'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile44-his68'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile44-his68'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile44-his68'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile44-his68'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ile61-leu67'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ile61-leu67'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ile61-leu67'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ile61-leu67'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ile61-leu67'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ile61-leu67'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'leu15-ile30'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'leu15-ile30'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'leu15-ile30'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'leu15-ile30'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'leu15-ile30'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'leu15-ile30'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'leu15-val26'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'leu15-val26'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'leu15-val26'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'leu15-val26'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'leu15-val26'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'leu15-val26'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'leu43-leu50'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'leu43-leu50'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'leu43-leu50'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'leu43-leu50'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'leu43-leu50'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'leu43-leu50'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'leu43-leu69'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'leu43-leu69'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'leu43-leu69'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'leu43-leu69'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'leu43-leu69'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'leu43-leu69'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'leu50-ile61'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'leu50-ile61'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'leu50-ile61'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'leu50-ile61'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'leu50-ile61'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'leu50-ile61'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'leu50-tyr59'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'leu50-tyr59'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'leu50-tyr59'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'leu50-tyr59'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'leu50-tyr59'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'leu50-tyr59'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'leu56-ile61np'         )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'leu56-ile61np'         )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'leu56-ile61np'         )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'leu56-ile61np'         )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'leu56-ile61np'         )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'leu56-ile61np'         )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'leu56-ile61p'          )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'leu56-ile61p'          )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'leu56-ile61p'          )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'leu56-ile61p'          )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'leu56-ile61p'          )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'leu56-ile61p'          )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'leu56-tyr59'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'leu56-tyr59'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'leu56-tyr59'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'leu56-tyr59'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'leu56-tyr59'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'leu56-tyr59'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'leu69-leu71'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'leu69-leu71'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'leu69-leu71'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'leu69-leu71'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'leu69-leu71'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'leu69-leu71'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'leu71-leu73'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'leu71-leu73'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'leu71-leu73'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'leu71-leu73'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'leu71-leu73'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'leu71-leu73'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'lys11-glu34'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'lys11-glu34'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'lys11-glu34'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'lys11-glu34'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'lys11-glu34'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'lys11-glu34'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'lys27-asp52'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'lys27-asp52'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'lys27-asp52'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'lys27-asp52'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'lys27-asp52'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'lys27-asp52'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'lys27-gln31'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'lys27-gln31'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'lys27-gln31'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'lys27-gln31'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'lys27-gln31'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'lys27-gln31'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'lys27-leu43-nonh3'     )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'lys27-leu43-nonh3'     )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'lys27-leu43-nonh3'     )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'lys27-leu43-nonh3'     )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'lys27-leu43-nonh3'     )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'lys27-leu43-nonh3'     )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'lys27-leu43'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'lys27-leu43'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'lys27-leu43'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'lys27-leu43'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'lys27-leu43'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'lys27-leu43'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'lys29-lys33'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'lys29-lys33'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'lys29-lys33'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'lys29-lys33'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'lys29-lys33'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'lys29-lys33'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'lys6-thr66'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'lys6-thr66'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'lys6-thr66'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'lys6-thr66'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'lys6-thr66'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'lys6-thr66'            )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'met1-ile3'             )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'met1-ile3'             )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'met1-ile3'             )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'met1-ile3'             )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'met1-ile3'             )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'met1-ile3'             )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'met1-lys63'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'met1-lys63'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'met1-lys63'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'met1-lys63'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'met1-lys63'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'met1-lys63'            )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'met1-val17-bi'         )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'met1-val17-bi'         )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'met1-val17-bi'         )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'met1-val17-bi'         )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'met1-val17-bi'         )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'met1-val17-bi'         )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'met1-val17-mono'       )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'met1-val17-mono'       )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'met1-val17-mono'       )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'met1-val17-mono'       )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'met1-val17-mono'       )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'met1-val17-mono'       )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'met1-val17'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'met1-val17'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'met1-val17'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'met1-val17'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'met1-val17'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'met1-val17'            )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'phe4-leu67'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'phe4-leu67'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'phe4-leu67'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'phe4-leu67'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'phe4-leu67'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'phe4-leu67'            )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'phe4-ser65'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'phe4-ser65'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'phe4-ser65'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'phe4-ser65'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'phe4-ser65'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'phe4-ser65'            )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'phe4-thr12'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'phe4-thr12'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'phe4-thr12'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'phe4-thr12'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'phe4-thr12'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'phe4-thr12'            )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'phe4-thr14'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'phe4-thr14'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'phe4-thr14'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'phe4-thr14'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'phe4-thr14'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'phe4-thr14'            )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'phe45-ala46'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'phe45-ala46'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'phe45-ala46'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'phe45-ala46'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'phe45-ala46'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'phe45-ala46'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'phe45-ile61'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'phe45-ile61'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'phe45-ile61'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'phe45-ile61'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'phe45-ile61'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'phe45-ile61'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'phe45-leu67'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'phe45-leu67'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'phe45-leu67'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'phe45-leu67'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'phe45-leu67'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'phe45-leu67'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'phe45-lys48'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'phe45-lys48'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'phe45-lys48'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'phe45-lys48'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'phe45-lys48'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'phe45-lys48'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'pro19-ser57'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'pro19-ser57'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'pro19-ser57'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'pro19-ser57'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'pro19-ser57'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'pro19-ser57'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'pro37-gln40'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'pro37-gln40'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'pro37-gln40'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'pro37-gln40'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'pro37-gln40'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'pro37-gln40'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'ser57-asn60'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'ser57-asn60'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'ser57-asn60'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'ser57-asn60'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'ser57-asn60'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'ser57-asn60'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'thr14-lys33'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'thr14-lys33'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'thr14-lys33'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'thr14-lys33'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'thr14-lys33'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'thr14-lys33'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'thr22-asn25'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'thr22-asn25'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'thr22-asn25'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'thr22-asn25'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'thr22-asn25'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'thr22-asn25'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'thr22-thr55-big'       )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'thr22-thr55-big'       )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'thr22-thr55-big'       )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'thr22-thr55-big'       )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'thr22-thr55-big'       )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'thr22-thr55-big'       )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'thr22-thr55-small'     )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'thr22-thr55-small'     )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'thr22-thr55-small'     )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'thr22-thr55-small'     )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'thr22-thr55-small'     )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'thr22-thr55-small'     )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'thr22-val26'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'thr22-val26'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'thr22-val26'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'thr22-val26'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'thr22-val26'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'thr22-val26'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'thr55-asp58-backbone'  )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'thr55-asp58-backbone'  )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'thr55-asp58-backbone'  )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'thr55-asp58-backbone'  )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'thr55-asp58-backbone'  )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'thr55-asp58-backbone'  )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'thr55-asp58-sidechain-small' )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'thr55-asp58-sidechain-small' )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'thr55-asp58-sidechain-small' )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'thr55-asp58-sidechain-small' )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'thr55-asp58-sidechain-small' )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'thr55-asp58-sidechain-small' )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'thr55-asp58-sidechain' )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'thr55-asp58-sidechain' )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'thr55-asp58-sidechain' )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'thr55-asp58-sidechain' )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'thr55-asp58-sidechain' )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'thr55-asp58-sidechain' )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'thr7-gly10'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'thr7-gly10'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'thr7-gly10'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'thr7-gly10'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'thr7-gly10'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'thr7-gly10'            )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'thr7-ile13'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'thr7-ile13'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'thr7-ile13'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'thr7-ile13'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'thr7-ile13'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'thr7-ile13'            )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'thr7-lys11'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'thr7-lys11'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'thr7-lys11'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'thr7-lys11'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'thr7-lys11'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'thr7-lys11'            )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'thr7-thr9-v2'          )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'thr7-thr9-v2'          )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'thr7-thr9-v2'          )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'thr7-thr9-v2'          )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'thr7-thr9-v2'          )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'thr7-thr9-v2'          )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'thr7-thr9'             )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'thr7-thr9'             )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'thr7-thr9'             )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'thr7-thr9'             )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'thr7-thr9'             )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'thr7-thr9'             )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'tyr59-ile61'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'tyr59-ile61'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'tyr59-ile61'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'tyr59-ile61'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'tyr59-ile61'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'tyr59-ile61'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'val17-leu56'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'val17-leu56'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'val17-leu56'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'val17-leu56'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'val17-leu56'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'val17-leu56'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'val26-ile30'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'val26-ile30'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'val26-ile30'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'val26-ile30'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'val26-ile30'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'val26-ile30'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'val26-leu43'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'val26-leu43'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'val26-leu43'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'val26-leu43'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'val26-leu43'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'val26-leu43'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'val26-leu56'           )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'val26-leu56'           )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'val26-leu56'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'val26-leu56'           )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'val26-leu56'           )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'val26-leu56'           )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'val5-ile13'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'val5-ile13'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'val5-ile13'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'val5-ile13'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'val5-ile13'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'val5-ile13'            )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'val5-leu15'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'val5-leu15'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'val5-leu15'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'val5-leu15'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'val5-leu15'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'val5-leu15'            )] = """Monomer B from  """
TAGL['%s-%s'            % (dbse, 'val5-leu69'            )] = """ """
TAGL['%s-%s-dimer'      % (dbse, 'val5-leu69'            )] = """Dimer from  """
TAGL['%s-%s-monoA-CP'   % (dbse, 'val5-leu69'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-CP'   % (dbse, 'val5-leu69'            )] = """Monomer B from  """
TAGL['%s-%s-monoA-unCP' % (dbse, 'val5-leu69'            )] = """Monomer A from  """
TAGL['%s-%s-monoB-unCP' % (dbse, 'val5-leu69'            )] = """Monomer B from  """

TAGL['dbse'] = 'interaction energies for bimolecular complexes from native 1UBQ protein fold'
TAGL['HB'] = 'hydrogen-bonded systems'
TAGL['MX'] = 'mixed-influence systems'
TAGL['DD'] = 'dispersion-dominated systems'
#TAGL['large'] = 'most computationally expensive systems'
TAGL['default'] = 'entire database'
#TAGL['small'] = 'few computationally quick systems'

# <<< Geometry Specification Strings >>>
GEOS = {}

GEOS['%s-%s-%s' % (dbse, 'ala28-asp32', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               38.794000000000    25.761000000000    13.880000000000
    H               39.533000000000    26.227000000000    14.535000000000
    C               38.728000000000    26.591000000000    12.611000000000
    O               39.704000000000    27.346000000000    12.277000000000
    N               37.633000000000    26.543000000000    11.867000000000
    H               36.939000000000    25.831000000000    12.081000000000
    H               37.512390000000    27.174350000000    10.974330000000
    H               39.144810000000    24.742860000000    13.655650000000
    H               37.823590000000    25.747510000000    14.397800000000
    --
    0 1
    C               41.092000000000    30.808000000000    12.851000000000
    O               41.828000000000    31.808000000000    12.681000000000
    N               41.001000000000    29.878000000000    11.931000000000
    H               40.442000000000    29.048000000000    12.087000000000
    C               41.718000000000    30.022000000000    10.643000000000
    H               42.786000000000    30.042000000000    10.870000000000
    H               41.486150000000    30.978480000000    10.151680000000
    H               41.487830000000    29.128670000000    10.043850000000
    H               40.503330000000    30.593420000000    13.755110000000

""")

GEOS['%s-%s-%s' % (dbse, 'ala46-lys48', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               19.172000000000    29.808000000000    21.243000000000
    H               19.577000000000    30.558000000000    20.693000000000
    C               18.453000000000    28.941000000000    20.591000000000
    O               17.860000000000    27.994000000000    21.128000000000
    C               18.443000000000    29.143000000000    19.083000000000
    H               17.837000000000    30.025000000000    18.873000000000
    N               19.810000000000    29.378000000000    18.578000000000
    H               19.931000000000    30.108000000000    17.886000000000
    C               20.835000000000    28.629000000000    18.904000000000
    O               20.821000000000    27.734000000000    19.749000000000
    H               18.017260000000    28.285640000000    18.541080000000
    H               21.756670000000    28.938130000000    18.389260000000
    H               19.346280000000    29.874030000000    22.327100000000
    --
    1 1
    C               19.634000000000    23.885000000000    21.531000000000
    H               20.606000000000    23.515000000000    21.199000000000
    H               19.157000000000    23.068000000000    22.076000000000
    C               18.791000000000    24.221000000000    20.313000000000
    H               19.267000000000    25.002000000000    19.710000000000
    H               18.697000000000    23.306000000000    19.723000000000
    N               17.440000000000    24.655000000000    20.827000000000
    H               17.514000000000    25.570000000000    21.262000000000
    H               16.813000000000    24.738000000000    20.037000000000
    H               17.066000000000    23.982000000000    21.479000000000
    H               19.727020000000    24.751270000000    22.202510000000

""")

GEOS['%s-%s-%s' % (dbse, 'arg42-val70', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               33.589000000000    30.189000000000    22.181000000000
    O               33.580000000000    29.009000000000    22.499000000000
    N               32.478000000000    30.917000000000    22.269000000000
    H               32.500000000000    31.907000000000    22.035000000000
    C               31.200000000000    30.329000000000    22.780000000000
    H               31.300000000000    29.251000000000    22.923000000000
    C               30.210000000000    30.509000000000    21.650000000000
    O               29.978000000000    31.726000000000    21.269000000000
    N               29.694000000000    29.436000000000    21.054000000000
    H               29.926000000000    28.506000000000    21.391000000000
    H               29.003640000000    29.537480000000    20.203640000000
    H               34.423830000000    30.687430000000    21.666590000000
    H               30.942690000000    30.767820000000    23.755310000000
    --
    0 1
    C               30.052000000000    35.042000000000    20.004000000000
    O               30.105000000000    36.305000000000    19.788000000000
    N               30.124000000000    34.533000000000    21.191000000000
    H               30.079000000000    33.523000000000    21.295000000000
    C               30.479000000000    35.369000000000    22.374000000000
    H               30.505000000000    36.428000000000    22.117000000000
    C               31.901000000000    34.910000000000    22.728000000000
    O               32.190000000000    33.696000000000    22.635000000000
    N               32.763000000000    35.831000000000    23.090000000000
    H               32.466000000000    36.787000000000    23.210000000000
    H               29.867860000000    34.383940000000    19.141990000000
    H               29.750620000000    35.233020000000    23.187000000000
    H               33.789860000000    35.564250000000    23.380520000000

""")

GEOS['%s-%s-%s' % (dbse, 'asn25-lys29', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               35.590000000000    21.945000000000    13.302000000000
    H               36.681000000000    21.885000000000    13.291000000000
    C               35.238000000000    23.382000000000    12.920000000000
    O               36.066000000000    24.109000000000    12.333000000000
    N               34.007000000000    23.745000000000    13.250000000000
    H               33.397000000000    23.059000000000    13.687000000000
    H               33.649460000000    24.764840000000    13.044830000000
    H               35.252140000000    21.704520000000    14.320830000000
    H               35.212480000000    21.235900000000    12.550560000000
    --
    0 1
    C               38.728000000000    26.591000000000    12.611000000000
    O               39.704000000000    27.346000000000    12.277000000000
    N               37.633000000000    26.543000000000    11.867000000000
    H               36.939000000000    25.831000000000    12.081000000000
    C               37.471000000000    27.391000000000    10.668000000000
    H               38.331000000000    27.243000000000    10.009000000000
    H               38.775830000000    25.989460000000    13.530710000000
    H               36.547490000000    27.150370000000    10.120980000000
    H               37.449570000000    28.456040000000    10.942300000000

""")

GEOS['%s-%s-%s' % (dbse, 'asp21-leu56', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               29.258000000000    17.318000000000     9.984000000000
    O               29.930000000000    16.477000000000    10.606000000000
    N               29.599000000000    18.599000000000     9.828000000000
    H               29.037000000000    19.248000000000     9.293000000000
    C               30.796000000000    19.083000000000    10.566000000000
    H               31.634000000000    18.398000000000    10.426000000000
    C               30.491000000000    19.162000000000    12.040000000000
    O               29.367000000000    19.523000000000    12.441000000000
    N               31.510000000000    18.936000000000    12.852000000000
    H               32.418000000000    18.706000000000    12.477000000000
    H               31.424680000000    19.033500000000    13.944340000000
    H               28.381280000000    16.966290000000     9.420391000000
    H               31.048400000000    20.089800000000    10.201810000000
    --
    0 1
    C               26.482000000000    19.280000000000    14.432000000000
    O               25.609000000000    19.388000000000    15.287000000000
    N               26.585000000000    20.063000000000    13.378000000000
    H               27.373000000000    19.916000000000    12.751000000000
    C               25.594000000000    21.109000000000    13.072000000000
    H               25.497000000000    21.778000000000    13.926000000000
    C               24.241000000000    20.436000000000    12.857000000000
    O               23.264000000000    20.951000000000    13.329000000000
    N               24.240000000000    19.233000000000    12.246000000000
    H               25.095000000000    18.802000000000    11.905000000000
    H               25.943230000000    21.664210000000    12.188950000000
    H               27.258440000000    18.506400000000    14.525140000000
    H               23.264740000000    18.751300000000    12.082220000000

""")

GEOS['%s-%s-%s' % (dbse, 'gln2-glu64-small', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               26.100000000000    29.253000000000     5.202000000000
    O               24.865000000000    29.024000000000     5.330000000000
    N               26.849000000000    29.656000000000     6.217000000000
    H               27.846000000000    29.802000000000     6.083000000000
    C               26.850000000000    29.021000000000     3.898000000000
    H               27.909000000000    28.907000000000     4.126000000000
    H               26.471480000000    28.101520000000     3.427602000000
    H               26.764580000000    29.843840000000     3.172994000000
    H               26.391250000000    29.955700000000     7.171275000000
    --
    0 1
    N               22.099000000000    29.163000000000     5.605000000000
    H               23.043000000000    28.873000000000     5.376000000000
    C               21.127000000000    28.240000000000     5.574000000000
    O               19.958000000000    28.465000000000     5.842000000000
    C               21.907000000000    30.563000000000     5.881000000000
    H               22.892000000000    31.025000000000     5.806000000000
    H               21.281980000000    31.029650000000     5.105378000000
    H               21.583670000000    30.848940000000     6.892780000000
    H               21.508060000000    27.236550000000     5.333403000000

""")

GEOS['%s-%s-%s' % (dbse, 'gln31-gly35', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               40.269000000000    30.508000000000    14.115000000000
    H               40.236000000000    31.406000000000    14.731000000000
    C               41.092000000000    30.808000000000    12.851000000000
    O               41.828000000000    31.808000000000    12.681000000000
    N               41.001000000000    29.878000000000    11.931000000000
    H               40.442000000000    29.048000000000    12.087000000000
    H               41.533500000000    29.984950000000    10.974440000000
    H               40.793930000000    29.707250000000    14.656530000000
    H               39.250560000000    30.215550000000    13.819610000000
    --
    0 1
    C               40.675000000000    35.527000000000    13.200000000000
    O               40.814000000000    36.528000000000    13.911000000000
    N               41.317000000000    34.393000000000    13.432000000000
    H               41.204000000000    33.614000000000    12.789000000000
    C               42.345000000000    34.269000000000    14.431000000000
    H               42.983000000000    33.429000000000    14.158000000000
    H               42.969000000000    35.163000000000    14.392000000000
    H               42.050310000000    34.125380000000    15.481010000000
    H               39.945760000000    35.541600000000    12.376590000000

""")

GEOS['%s-%s-%s' % (dbse, 'gln40-arg72', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               37.738000000000    31.637000000000    23.712000000000
    H               37.651000000000    31.725000000000    24.796000000000
    C               36.334000000000    31.742000000000    23.087000000000
    O               35.574000000000    32.618000000000    23.483000000000
    N               36.000000000000    30.860000000000    22.172000000000
    H               36.664000000000    30.123000000000    21.952000000000
    H               35.037800000000    30.871440000000    21.639050000000
    H               38.243140000000    30.699400000000    23.436800000000
    H               38.300570000000    32.503640000000    23.334580000000
    --
    0 1
    C               34.239000000000    35.353000000000    24.979000000000
    O               33.707000000000    36.197000000000    25.728000000000
    N               34.930000000000    34.384000000000    25.451000000000
    H               35.349000000000    33.731000000000    24.799000000000
    C               35.161000000000    34.174000000000    26.896000000000
    H               34.663000000000    34.917000000000    27.522000000000
    H               36.248630000000    34.261880000000    27.035020000000
    H               34.842310000000    33.159070000000    27.175930000000
    H               34.170330000000    35.439940000000    23.884590000000

""")

GEOS['%s-%s-%s' % (dbse, 'gln41-arg72', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               34.738000000000    30.875000000000    21.473000000000
    H               34.456000000000    31.916000000000    21.318000000000
    C               33.589000000000    30.189000000000    22.181000000000
    O               33.580000000000    29.009000000000    22.499000000000
    N               32.478000000000    30.917000000000    22.269000000000
    H               32.500000000000    31.907000000000    22.035000000000
    H               31.538740000000    30.484850000000    22.644560000000
    H               34.835870000000    30.422530000000    20.475150000000
    H               35.700200000000    30.863560000000    22.005950000000
    --
    1 1
    C               35.612000000000    30.577000000000    28.044000000000
    H               36.605000000000    30.127000000000    28.125000000000
    H               35.004000000000    30.265000000000    28.894000000000
    N               35.040000000000    30.252000000000    26.730000000000
    H               35.512000000000    30.595000000000    25.905000000000
    C               34.338000000000    29.103000000000    26.650000000000
    N               34.110000000000    28.437000000000    27.768000000000
    H               34.564000000000    28.732000000000    28.609000000000
    H               33.689000000000    27.527000000000    27.717000000000
    N               34.014000000000    28.657000000000    25.457000000000
    H               34.369000000000    29.108000000000    24.627000000000
    H               33.716000000000    27.705000000000    25.340000000000
    H               35.715480000000    31.668760000000    28.129740000000

""")

GEOS['%s-%s-%s' % (dbse, 'gln62-ser65', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               21.110000000000    24.111000000000    10.173000000000
    O               21.841000000000    23.198000000000     9.686000000000
    N               20.291000000000    24.875000000000     9.507000000000
    H               19.702000000000    25.510000000000    10.033000000000
    C               20.081000000000    24.773000000000     8.033000000000
    H               20.519000000000    23.844000000000     7.665000000000
    C               20.822000000000    25.914000000000     7.332000000000
    O               21.323000000000    26.830000000000     8.008000000000
    N               20.924000000000    25.862000000000     6.006000000000
    H               20.502000000000    25.078000000000     5.521000000000
    H               19.004050000000    24.746110000000     7.810633000000
    H               21.480590000000    26.610970000000     5.423553000000
    H               21.163330000000    24.220540000000    11.266230000000
    --
    0 1
    C               21.466000000000    30.953000000000     7.261000000000
    O               21.066000000000    32.112000000000     7.533000000000
    N               21.674000000000    30.034000000000     8.191000000000
    H               21.986000000000    29.123000000000     7.887000000000
    C               21.419000000000    30.253000000000     9.620000000000
    H               20.418000000000    30.649000000000     9.793000000000
    C               21.637000000000    28.923000000000    10.353000000000
    H               22.580000000000    28.468000000000    10.043000000000
    H               21.669000000000    29.098000000000    11.428000000000
    O               20.544000000000    28.047000000000    10.059000000000
    H               20.792000000000    27.597000000000     9.215000000000
    H               22.190350000000    30.946150000000     9.986835000000
    H               21.789330000000    30.667060000000     6.249220000000

""")

GEOS['%s-%s-%s' % (dbse, 'glu16-lys29', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               31.022000000000    29.288000000000     5.665000000000
    O               29.809000000000    29.395000000000     5.545000000000
    N               31.834000000000    28.412000000000     5.125000000000
    H               32.837000000000    28.458000000000     5.260000000000
    C               31.220000000000    27.341000000000     4.275000000000
    H               30.149000000000    27.499000000000     4.153000000000
    C               31.827000000000    27.262000000000     2.894000000000
    H               32.916000000000    27.263000000000     2.970000000000
    H               31.513000000000    26.319000000000     2.443000000000
    O               32.576000000000    25.802000000000     5.461000000000
    C               31.440000000000    26.079000000000     5.080000000000
    N               30.310000000000    25.458000000000     5.384000000000
    H               29.430000000000    25.804000000000     5.018000000000
    H               30.293400000000    24.542960000000     5.994275000000
    H               31.497660000000    28.076820000000     2.232488000000
    H               31.491820000000    29.995950000000     6.363626000000
    --
    1 1
    C               34.758000000000    25.280000000000     8.900000000000
    H               34.198000000000    24.848000000000     9.730000000000
    H               34.230000000000    26.175000000000     8.568000000000
    C               34.793000000000    24.264000000000     7.767000000000
    H               35.625000000000    23.572000000000     7.930000000000
    H               33.866000000000    23.679000000000     7.782000000000
    N               34.914000000000    24.944000000000     6.441000000000
    H               35.767000000000    25.474000000000     6.361000000000
    H               34.880000000000    24.229000000000     5.720000000000
    H               34.097000000000    25.520000000000     6.254000000000
    H               35.765290000000    25.525510000000     9.267536000000

""")

GEOS['%s-%s-%s' % (dbse, 'glu18-asp21-C', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               29.279000000000    23.227000000000     5.641000000000
    O               28.478000000000    23.522000000000     4.725000000000
    N               29.380000000000    22.057000000000     6.232000000000
    H               30.114000000000    21.905000000000     6.927000000000
    C               28.468000000000    20.940000000000     5.980000000000
    H               27.743000000000    21.232000000000     5.229000000000
    C               29.213000000000    19.697000000000     5.506000000000
    H               30.051000000000    19.495000000000     6.176000000000
    H               28.530000000000    18.849000000000     5.568000000000
    C               27.819000000000    20.609000000000     7.316000000000
    O               28.449000000000    20.674000000000     8.360000000000
    N               26.559000000000    20.220000000000     7.288000000000
    H               25.932730000000    20.136990000000     6.387503000000
    H               26.010510000000    19.923220000000     8.194135000000
    H               29.581800000000    19.738530000000     4.470499000000
    H               30.001620000000    23.956060000000     6.036327000000
    --
    0 1
    C               29.258000000000    17.318000000000     9.984000000000
    O               29.930000000000    16.477000000000    10.606000000000
    N               29.599000000000    18.599000000000     9.828000000000
    H               29.037000000000    19.248000000000     9.293000000000
    C               30.796000000000    19.083000000000    10.566000000000
    H               31.634000000000    18.398000000000    10.426000000000
    C               31.155000000000    20.515000000000    10.048000000000
    H               30.266000000000    21.137000000000     9.956000000000
    H               31.808000000000    20.995000000000    10.779000000000
    C               30.491000000000    19.162000000000    12.040000000000
    O               29.367000000000    19.523000000000    12.441000000000
    N               31.510000000000    18.936000000000    12.852000000000
    H               32.418000000000    18.706000000000    12.477000000000
    H               31.715970000000    20.457300000000     9.103552000000
    H               28.381280000000    16.966290000000     9.420391000000
    H               31.424680000000    19.033500000000    13.944340000000

""")

GEOS['%s-%s-%s' % (dbse, 'glu18-asp21-N', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               29.279000000000    23.227000000000     5.641000000000
    O               28.478000000000    23.522000000000     4.725000000000
    N               29.380000000000    22.057000000000     6.232000000000
    H               30.114000000000    21.905000000000     6.927000000000
    C               28.468000000000    20.940000000000     5.980000000000
    H               27.743000000000    21.232000000000     5.229000000000
    C               29.213000000000    19.697000000000     5.506000000000
    H               30.051000000000    19.495000000000     6.176000000000
    H               28.530000000000    18.849000000000     5.568000000000
    C               27.819000000000    20.609000000000     7.316000000000
    O               28.449000000000    20.674000000000     8.360000000000
    N               26.559000000000    20.220000000000     7.288000000000
    H               26.010510000000    19.923220000000     8.194135000000
    H               25.932730000000    20.136990000000     6.387503000000
    H               30.001620000000    23.956060000000     6.036327000000
    H               29.581800000000    19.738530000000     4.470499000000
    --
    -1 1
    C               31.155000000000    20.515000000000    10.048000000000
    H               30.266000000000    21.137000000000     9.956000000000
    H               31.808000000000    20.995000000000    10.779000000000
    C               31.923000000000    20.436000000000     8.755000000000
    O               32.493000000000    19.374000000000     8.456000000000
    O               31.838000000000    21.402000000000     7.968000000000
    H               30.902600000000    19.508200000000    10.412190000000

""")

GEOS['%s-%s-%s' % (dbse, 'glu24-ala28', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               35.031000000000    21.722000000000    17.069000000000
    H               35.253000000000    22.459000000000    17.842000000000
    C               35.615000000000    22.190000000000    15.759000000000
    O               36.532000000000    23.046000000000    15.724000000000
    N               35.139000000000    21.624000000000    14.662000000000
    H               34.436000000000    20.901000000000    14.762000000000
    H               35.476860000000    21.864480000000    13.643170000000
    H               33.943920000000    21.578330000000    16.981770000000
    H               35.488310000000    20.759200000000    17.340800000000
    --
    0 1
    C               36.975000000000    26.826000000000    15.107000000000
    O               37.579000000000    27.926000000000    15.159000000000
    N               37.499000000000    25.743000000000    14.571000000000
    H               37.012000000000    24.858000000000    14.664000000000
    C               38.794000000000    25.761000000000    13.880000000000
    H               39.533000000000    26.227000000000    14.535000000000
    H               39.144810000000    24.742860000000    13.655650000000
    H               38.746170000000    26.362540000000    12.960290000000
    H               35.976870000000    26.745660000000    15.562280000000

""")

GEOS['%s-%s-%s' % (dbse, 'glu24-asp52', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               32.776000000000    22.519000000000    16.577000000000
    O               33.233000000000    23.659000000000    16.384000000000
    N               33.548000000000    21.526000000000    16.950000000000
    H               33.142000000000    20.631000000000    17.211000000000
    C               35.031000000000    21.722000000000    17.069000000000
    H               35.253000000000    22.459000000000    17.842000000000
    H               31.706190000000    22.290370000000    16.461970000000
    H               35.488310000000    20.759200000000    17.340800000000
    H               35.456800000000    22.063220000000    16.113870000000
    --
    0 1
    C               32.128000000000    19.364000000000    19.750000000000
    O               32.546000000000    19.317000000000    18.558000000000
    N               31.697000000000    18.311000000000    20.406000000000
    H               31.456000000000    18.438000000000    21.382000000000
    C               32.262000000000    20.670000000000    20.514000000000
    H               32.330000000000    21.430000000000    19.743000000000
    H               31.600760000000    17.304590000000    19.972550000000
    H               31.486120000000    20.986690000000    21.226540000000
    H               33.233880000000    20.702490000000    21.028190000000

""")

GEOS['%s-%s-%s' % (dbse, 'glu51-arg54', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    -1 1
    C               29.247000000000    19.456000000000    23.705000000000
    H               29.233000000000    18.605000000000    23.022000000000
    H               30.189000000000    19.974000000000    23.884000000000
    C               28.722000000000    19.047000000000    25.066000000000
    O               29.139000000000    18.132000000000    25.746000000000
    O               27.777000000000    19.842000000000    25.367000000000
    H               28.607540000000    20.227050000000    23.250480000000
    --
    1 1
    N               26.975000000000    15.521000000000    20.942000000000
    H               27.242000000000    15.890000000000    21.837000000000
    C               27.603000000000    14.423000000000    20.655000000000
    N               27.479000000000    13.733000000000    19.537000000000
    H               26.860000000000    14.107000000000    18.824000000000
    H               28.115000000000    12.997000000000    19.302000000000
    N               28.519000000000    13.967000000000    21.550000000000
    H               28.711000000000    14.480000000000    22.391000000000
    H               29.094000000000    13.183000000000    21.299000000000
    H               26.240790000000    16.166950000000    20.438330000000

""")

GEOS['%s-%s-%s' % (dbse, 'glu51-tyr59', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               27.994000000000    23.781000000000    21.643000000000
    O               28.904000000000    24.444000000000    22.098000000000
    N               27.942000000000    22.448000000000    21.648000000000
    H               27.152000000000    21.956000000000    21.249000000000
    C               29.015000000000    21.657000000000    22.288000000000
    H               29.591000000000    22.266000000000    22.985000000000
    H               27.148660000000    24.316570000000    21.186310000000
    H               29.693110000000    21.253940000000    21.521380000000
    H               28.530990000000    20.846440000000    22.852560000000
    --
    0 1
    C               22.945000000000    21.951000000000    17.785000000000
    C               24.272000000000    21.544000000000    17.644000000000
    H               24.693000000000    21.388000000000    16.659000000000
    C               25.052000000000    21.285000000000    18.776000000000
    H               26.081000000000    20.990000000000    18.667000000000
    C               24.517000000000    21.470000000000    20.030000000000
    O               25.248000000000    21.302000000000    21.191000000000
    H               24.728000000000    21.394000000000    21.993000000000
    C               23.204000000000    21.907000000000    20.192000000000
    H               22.798000000000    22.050000000000    21.180000000000
    C               22.437000000000    22.157000000000    19.065000000000
    H               21.415000000000    22.483000000000    19.173000000000
    H               22.318640000000    22.171680000000    16.908090000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile13-leu15', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               32.995000000000    35.883000000000     9.934000000000
    H               33.813000000000    36.384000000000     9.412000000000
    C               33.109000000000    36.381000000000    11.435000000000
    H               32.449000000000    35.800000000000    12.081000000000
    H               34.130000000000    36.321000000000    11.799000000000
    H               32.849000000000    37.436000000000    11.509000000000
    C               33.306000000000    34.381000000000     9.840000000000
    H               32.459000000000    33.813000000000    10.228000000000
    H               33.449000000000    34.098000000000     8.802000000000
    C               34.535000000000    34.028000000000    10.720000000000
    H               34.329000000000    34.129000000000    11.785000000000
    H               34.809000000000    32.988000000000    10.545000000000
    H               35.346000000000    34.662000000000    10.400000000000
    H               32.083020000000    36.173400000000     9.391819000000
    --
    0 1
    C               31.562000000000    29.686000000000     8.045000000000
    H               30.790000000000    30.265000000000     8.559000000000
    H               31.073000000000    28.716000000000     7.948000000000
    C               32.631000000000    29.444000000000     9.060000000000
    H               33.028000000000    28.443000000000     8.889000000000
    C               33.814000000000    30.390000000000     9.030000000000
    H               33.465000000000    31.415000000000     9.143000000000
    H               34.503000000000    30.145000000000     9.839000000000
    H               34.339000000000    30.288000000000     8.080000000000
    C               31.945000000000    29.449000000000    10.436000000000
    H               31.086000000000    28.781000000000    10.437000000000
    H               32.653000000000    29.096000000000    11.183000000000
    H               31.642000000000    30.461000000000    10.683000000000
    H               31.644750000000    30.109820000000     7.033303000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile13-lys33', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               32.995000000000    35.883000000000     9.934000000000
    H               33.813000000000    36.384000000000     9.412000000000
    C               33.109000000000    36.381000000000    11.435000000000
    H               32.449000000000    35.800000000000    12.081000000000
    H               34.130000000000    36.321000000000    11.799000000000
    H               32.849000000000    37.436000000000    11.509000000000
    C               33.306000000000    34.381000000000     9.840000000000
    H               32.459000000000    33.813000000000    10.228000000000
    H               33.449000000000    34.098000000000     8.802000000000
    C               34.535000000000    34.028000000000    10.720000000000
    H               34.329000000000    34.129000000000    11.785000000000
    H               34.809000000000    32.988000000000    10.545000000000
    H               35.346000000000    34.662000000000    10.400000000000
    H               32.083020000000    36.173400000000     9.391819000000
    --
    1 1
    C               37.220000000000    32.822000000000     8.827000000000
    H               36.843000000000    31.804000000000     8.932000000000
    H               37.193000000000    33.302000000000     9.802000000000
    C               36.351000000000    33.613000000000     7.838000000000
    H               35.342000000000    33.632000000000     8.206000000000
    H               36.727000000000    34.638000000000     7.774000000000
    C               36.322000000000    32.944000000000     6.477000000000
    H               37.337000000000    32.652000000000     6.186000000000
    H               35.709000000000    32.038000000000     6.517000000000
    N               35.768000000000    33.945000000000     5.489000000000
    H               36.271000000000    34.819000000000     5.550000000000
    H               35.864000000000    33.572000000000     4.552000000000
    H               34.779000000000    34.081000000000     5.667000000000
    H               38.253730000000    32.806440000000     8.451299000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile23-arg54', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               31.593000000000    20.553000000000    14.655000000000
    O               32.159000000000    21.311000000000    13.861000000000
    N               31.113000000000    20.863000000000    15.860000000000
    H               30.636000000000    20.148000000000    16.405000000000
    C               31.288000000000    22.201000000000    16.417000000000
    H               30.905000000000    22.936000000000    15.712000000000
    H               31.454290000000    19.493820000000    14.392520000000
    H               32.357810000000    22.429630000000    16.532030000000
    H               30.744270000000    22.271090000000    17.370650000000
    --
    0 1
    C               28.108000000000    17.439000000000    18.276000000000
    H               27.847000000000    16.388000000000    18.150000000000
    C               28.375000000000    17.999000000000    16.887000000000
    O               29.326000000000    18.786000000000    16.690000000000
    N               27.510000000000    17.689000000000    15.954000000000
    H               26.729000000000    17.079000000000    16.191000000000
    H               27.557550000000    18.062720000000    14.920520000000
    H               29.025470000000    17.554430000000    18.871760000000
    H               27.258410000000    17.979520000000    18.718760000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile23-leu50', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               30.520000000000    22.300000000000    17.764000000000
    H               30.863000000000    21.515000000000    18.413000000000
    C               30.832000000000    23.699000000000    18.358000000000
    H               30.685000000000    24.476000000000    17.610000000000
    H               30.172000000000    23.874000000000    19.207000000000
    H               31.859000000000    23.726000000000    18.727000000000
    C               29.006000000000    22.043000000000    17.442000000000
    H               28.831000000000    21.016000000000    17.160000000000
    H               28.443000000000    22.207000000000    18.360000000000
    H               28.795000000000    22.680000000000    15.384000000000
    H               27.326000000000    22.805000000000    16.351000000000
    H               28.618000000000    23.995000000000    16.580000000000
    C               28.407000000000    22.948000000000    16.366000000000
    H               31.063730000000    22.229910000000    16.810350000000
    --
    0 1
    C               27.043000000000    24.992000000000    19.571000000000
    H               27.973000000000    25.560000000000    19.529000000000
    H               27.174000000000    24.103000000000    18.954000000000
    C               25.931000000000    25.844000000000    18.959000000000
    H               25.888000000000    26.813000000000    19.450000000000
    C               26.203000000000    26.083000000000    17.471000000000
    H               26.174000000000    25.135000000000    16.931000000000
    H               25.450000000000    26.753000000000    17.058000000000
    H               27.188000000000    26.530000000000    17.344000000000
    C               24.577000000000    25.190000000000    19.079000000000
    H               24.281000000000    25.115000000000    20.124000000000
    H               23.829000000000    25.764000000000    18.560000000000
    H               24.614000000000    24.192000000000    18.645000000000
    H               26.887140000000    24.653700000000    20.606020000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile23-leu56', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               30.520000000000    22.300000000000    17.764000000000
    H               30.863000000000    21.515000000000    18.413000000000
    C               30.832000000000    23.699000000000    18.358000000000
    H               30.685000000000    24.476000000000    17.610000000000
    H               30.172000000000    23.874000000000    19.207000000000
    H               31.859000000000    23.726000000000    18.727000000000
    C               29.006000000000    22.043000000000    17.442000000000
    H               28.831000000000    21.016000000000    17.160000000000
    H               28.443000000000    22.207000000000    18.360000000000
    H               28.795000000000    22.680000000000    15.384000000000
    H               27.326000000000    22.805000000000    16.351000000000
    H               28.618000000000    23.995000000000    16.580000000000
    C               28.407000000000    22.948000000000    16.366000000000
    H               31.063730000000    22.229910000000    16.810350000000
    --
    0 1
    C               26.084000000000    21.888000000000    11.833000000000
    H               26.170000000000    21.168000000000    11.029000000000
    H               25.313000000000    22.606000000000    11.554000000000
    C               27.426000000000    22.616000000000    11.902000000000
    H               28.230000000000    21.920000000000    12.129000000000
    C               27.718000000000    23.341000000000    10.578000000000
    H               26.935000000000    24.062000000000    10.350000000000
    H               28.683000000000    23.842000000000    10.643000000000
    H               27.773000000000    22.608000000000     9.776000000000
    C               27.380000000000    23.721000000000    12.955000000000
    H               27.051000000000    23.309000000000    13.905000000000
    H               28.369000000000    24.127000000000    13.111000000000
    H               26.698000000000    24.514000000000    12.653000000000
    H               25.734770000000    21.332790000000    12.716050000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile23-lys27', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               32.776000000000    22.519000000000    16.577000000000
    O               33.233000000000    23.659000000000    16.384000000000
    N               33.548000000000    21.526000000000    16.950000000000
    H               33.142000000000    20.631000000000    17.211000000000
    C               31.288000000000    22.201000000000    16.417000000000
    H               30.905000000000    22.936000000000    15.712000000000
    H               34.635080000000    21.669670000000    17.037230000000
    H               31.156140000000    21.192800000000    15.997290000000
    H               30.744270000000    22.271090000000    17.370650000000
    --
    0 1
    C               34.441000000000    26.099000000000    13.684000000000
    O               34.883000000000    27.090000000000    13.093000000000
    N               34.734000000000    25.822000000000    14.949000000000
    H               34.360000000000    24.976000000000    15.377000000000
    C               35.596000000000    26.715000000000    15.736000000000
    H               35.175000000000    27.721000000000    15.746000000000
    H               33.786230000000    25.376450000000    13.174890000000
    H               35.681600000000    26.346700000000    16.768970000000
    H               36.594130000000    26.795340000000    15.280720000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile3-leu15-big', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               26.100000000000    29.253000000000     5.202000000000
    O               24.865000000000    29.024000000000     5.330000000000
    N               26.849000000000    29.656000000000     6.217000000000
    H               27.846000000000    29.802000000000     6.083000000000
    C               26.235000000000    30.058000000000     7.497000000000
    H               25.173000000000    30.262000000000     7.348000000000
    C               26.882000000000    31.428000000000     7.862000000000
    O               27.906000000000    31.711000000000     7.264000000000
    N               26.214000000000    32.097000000000     8.771000000000
    H               25.335000000000    31.759000000000     9.150000000000
    H               26.313280000000    29.334060000000     8.321488000000
    H               26.642020000000    29.085330000000     4.259608000000
    H               26.619990000000    33.071220000000     9.080946000000
    --
    0 1
    C               31.409000000000    32.680000000000     6.446000000000
    O               32.619000000000    32.812000000000     6.125000000000
    N               30.884000000000    31.485000000000     6.666000000000
    H               29.884000000000    31.395000000000     6.821000000000
    C               31.677000000000    30.275000000000     6.639000000000
    H               32.719000000000    30.440000000000     6.361000000000
    C               31.022000000000    29.288000000000     5.665000000000
    O               29.809000000000    29.395000000000     5.545000000000
    N               31.834000000000    28.412000000000     5.125000000000
    H               32.837000000000    28.458000000000     5.260000000000
    H               30.749160000000    33.558800000000     6.494174000000
    H               31.594250000000    29.851180000000     7.650697000000
    H               31.383390000000    27.626000000000     4.501187000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile3-leu15-part1', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               26.235000000000    30.058000000000     7.497000000000
    H               25.173000000000    30.262000000000     7.348000000000
    C               26.882000000000    31.428000000000     7.862000000000
    O               27.906000000000    31.711000000000     7.264000000000
    N               26.214000000000    32.097000000000     8.771000000000
    H               25.335000000000    31.759000000000     9.150000000000
    H               26.692750000000    29.758300000000     6.542726000000
    H               26.313280000000    29.334060000000     8.321488000000
    H               26.619990000000    33.071220000000     9.080946000000
    --
    0 1
    C               31.409000000000    32.680000000000     6.446000000000
    O               32.619000000000    32.812000000000     6.125000000000
    N               30.884000000000    31.485000000000     6.666000000000
    H               29.884000000000    31.395000000000     6.821000000000
    C               31.677000000000    30.275000000000     6.639000000000
    H               32.719000000000    30.440000000000     6.361000000000
    H               30.749160000000    33.558800000000     6.494174000000
    H               31.207180000000    29.567050000000     5.940373000000
    H               31.594250000000    29.851180000000     7.650697000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile3-leu15-part2', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               26.100000000000    29.253000000000     5.202000000000
    O               24.865000000000    29.024000000000     5.330000000000
    N               26.849000000000    29.656000000000     6.217000000000
    H               27.846000000000    29.802000000000     6.083000000000
    C               26.235000000000    30.058000000000     7.497000000000
    H               25.173000000000    30.262000000000     7.348000000000
    H               26.691680000000    31.024990000000     7.754630000000
    H               26.313280000000    29.334060000000     8.321488000000
    H               26.642020000000    29.085330000000     4.259608000000
    --
    0 1
    H               32.719000000000    30.440000000000     6.361000000000
    C               31.677000000000    30.275000000000     6.639000000000
    C               31.022000000000    29.288000000000     5.665000000000
    O               29.809000000000    29.395000000000     5.545000000000
    N               31.834000000000    28.412000000000     5.125000000000
    H               32.837000000000    28.458000000000     5.260000000000
    H               31.383390000000    27.626000000000     4.501187000000
    H               31.074150000000    31.194860000000     6.659526000000
    H               31.594250000000    29.851180000000     7.650697000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile3-leu15', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               26.344000000000    29.050000000000     8.645000000000
    H               25.879000000000    29.504000000000     9.522000000000
    C               25.491000000000    27.771000000000     8.287000000000
    H               25.945000000000    27.245000000000     7.447000000000
    H               25.448000000000    27.099000000000     9.144000000000
    H               24.475000000000    28.063000000000     8.019000000000
    C               27.810000000000    28.748000000000     8.999000000000
    H               28.241000000000    28.093000000000     8.242000000000
    H               28.394000000000    29.667000000000     9.015000000000
    C               27.967000000000    28.087000000000    10.417000000000
    H               27.385000000000    27.167000000000    10.473000000000
    H               29.013000000000    27.846000000000    10.595000000000
    H               27.625000000000    28.782000000000    11.183000000000
    H               26.265720000000    29.773940000000     7.820512000000
    --
    0 1
    C               31.562000000000    29.686000000000     8.045000000000
    H               30.790000000000    30.265000000000     8.559000000000
    H               31.073000000000    28.716000000000     7.948000000000
    C               32.631000000000    29.444000000000     9.060000000000
    H               33.028000000000    28.443000000000     8.889000000000
    C               33.814000000000    30.390000000000     9.030000000000
    H               33.465000000000    31.415000000000     9.143000000000
    H               34.503000000000    30.145000000000     9.839000000000
    H               34.339000000000    30.288000000000     8.080000000000
    C               31.945000000000    29.449000000000    10.436000000000
    H               31.086000000000    28.781000000000    10.437000000000
    H               32.653000000000    29.096000000000    11.183000000000
    H               31.642000000000    30.461000000000    10.683000000000
    H               31.644750000000    30.109820000000     7.033303000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile3-leu67', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               26.344000000000    29.050000000000     8.645000000000
    H               25.879000000000    29.504000000000     9.522000000000
    C               25.491000000000    27.771000000000     8.287000000000
    H               25.945000000000    27.245000000000     7.447000000000
    H               25.448000000000    27.099000000000     9.144000000000
    H               24.475000000000    28.063000000000     8.019000000000
    C               27.810000000000    28.748000000000     8.999000000000
    H               28.241000000000    28.093000000000     8.242000000000
    H               28.394000000000    29.667000000000     9.015000000000
    C               27.967000000000    28.087000000000    10.417000000000
    H               27.385000000000    27.167000000000    10.473000000000
    H               29.013000000000    27.846000000000    10.595000000000
    H               27.625000000000    28.782000000000    11.183000000000
    H               26.265720000000    29.773940000000     7.820512000000
    --
    0 1
    C               26.310000000000    30.594000000000    14.967000000000
    H               27.244000000000    31.135000000000    14.803000000000
    H               26.368000000000    30.156000000000    15.965000000000
    C               26.290000000000    29.480000000000    13.960000000000
    H               26.459000000000    29.906000000000    12.971000000000
    C               27.393000000000    28.442000000000    14.229000000000
    H               27.257000000000    27.992000000000    15.212000000000
    H               27.352000000000    27.657000000000    13.474000000000
    H               28.371000000000    28.921000000000    14.181000000000
    C               24.942000000000    28.807000000000    13.952000000000
    H               24.133000000000    29.498000000000    13.737000000000
    H               24.937000000000    28.054000000000    13.165000000000
    H               24.769000000000    28.313000000000    14.905000000000
    H               25.481890000000    31.317980000000    14.976270000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile3-val17', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               26.344000000000    29.050000000000     8.645000000000
    H               25.879000000000    29.504000000000     9.522000000000
    C               25.491000000000    27.771000000000     8.287000000000
    H               25.945000000000    27.245000000000     7.447000000000
    H               25.448000000000    27.099000000000     9.144000000000
    H               24.475000000000    28.063000000000     8.019000000000
    C               27.810000000000    28.748000000000     8.999000000000
    H               28.241000000000    28.093000000000     8.242000000000
    H               28.394000000000    29.667000000000     9.015000000000
    C               27.967000000000    28.087000000000    10.417000000000
    H               27.385000000000    27.167000000000    10.473000000000
    H               29.013000000000    27.846000000000    10.595000000000
    H               27.625000000000    28.782000000000    11.183000000000
    H               26.265720000000    29.773940000000     7.820512000000
    --
    0 1
    C               29.903000000000    24.590000000000     7.665000000000
    H               29.900000000000    23.663000000000     8.237000000000
    C               30.862000000000    25.496000000000     8.389000000000
    H               30.895000000000    26.467000000000     7.902000000000
    H               30.547000000000    25.617000000000     9.425000000000
    H               31.853000000000    25.045000000000     8.370000000000
    C               28.476000000000    25.135000000000     7.705000000000
    H               27.762000000000    24.367000000000     7.412000000000
    H               28.235000000000    25.443000000000     8.721000000000
    H               28.371000000000    26.000000000000     7.051000000000
    H               30.174450000000    24.346750000000     6.627144000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile30-glu34', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               36.731000000000    30.570000000000    12.645000000000
    H               36.443000000000    31.212000000000    11.812000000000
    C               38.148000000000    30.981000000000    13.069000000000
    O               38.544000000000    32.150000000000    12.856000000000
    N               38.883000000000    30.110000000000    13.713000000000
    H               38.533000000000    29.186000000000    13.937000000000
    H               39.901440000000    30.402450000000    14.008390000000
    H               36.790720000000    29.524970000000    12.306860000000
    H               36.010140000000    30.715160000000    13.463100000000
    --
    0 1
    C               39.837000000000    34.271000000000     9.995000000000
    O               40.164000000000    35.323000000000     9.345000000000
    N               39.655000000000    34.335000000000    11.285000000000
    H               39.422000000000    33.465000000000    11.753000000000
    C               39.676000000000    35.547000000000    12.072000000000
    H               39.933000000000    36.404000000000    11.449000000000
    H               39.815550000000    33.326570000000     9.431448000000
    H               40.405240000000    35.532400000000    12.895400000000
    H               38.688610000000    35.737210000000    12.517960000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile30-leu43', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               35.708000000000    30.776000000000    13.806000000000
    H               35.901000000000    29.999000000000    14.547000000000
    C               35.874000000000    32.138000000000    14.512000000000
    H               35.708000000000    32.951000000000    13.805000000000
    H               35.191000000000    32.233000000000    15.353000000000
    H               36.876000000000    32.231000000000    14.929000000000
    C               34.228000000000    30.630000000000    13.319000000000
    H               33.937000000000    31.481000000000    12.702000000000
    H               34.143000000000    29.718000000000    12.732000000000
    C               33.284000000000    30.504000000000    14.552000000000
    H               33.209000000000    31.442000000000    15.098000000000
    H               32.287000000000    30.236000000000    14.204000000000
    H               33.632000000000    29.710000000000    15.211000000000
    H               36.428860000000    30.630840000000    12.987900000000
    --
    0 1
    C               29.151000000000    28.655000000000    18.755000000000
    H               29.215000000000    27.648000000000    19.172000000000
    H               28.322000000000    28.637000000000    18.045000000000
    C               30.416000000000    28.912000000000    17.980000000000
    H               31.242000000000    29.080000000000    18.671000000000
    C               30.738000000000    27.693000000000    17.122000000000
    H               29.921000000000    27.482000000000    16.432000000000
    H               31.647000000000    27.884000000000    16.550000000000
    H               30.912000000000    26.825000000000    17.756000000000
    C               30.205000000000    30.168000000000    17.129000000000
    H               30.017000000000    31.034000000000    17.761000000000
    H               31.110000000000    30.356000000000    16.551000000000
    H               29.373000000000    30.026000000000    16.440000000000
    H               28.870000000000    29.318130000000    19.586440000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile36-gln41', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               40.226000000000    33.716000000000    17.509000000000
    H               40.748000000000    34.446000000000    18.116000000000
    C               40.449000000000    32.278000000000    17.945000000000
    O               39.936000000000    31.336000000000    17.315000000000
    N               41.189000000000    32.085000000000    19.031000000000
    H               40.537020000000    33.865530000000    16.464530000000
    H               39.161530000000    33.986800000000    17.568720000000
    H               41.655400000000    32.897700000000    19.607190000000
    H               41.392090000000    31.088980000000    19.451360000000
    --
    0 1
    C               36.012000000000    30.860000000000    19.221000000000
    H               35.839000000000    31.925000000000    19.081000000000
    H               36.960000000000    30.720000000000    19.726000000000
    C               36.083000000000    30.194000000000    17.875000000000
    O               35.048000000000    29.702000000000    17.393000000000
    N               37.228000000000    30.126000000000    17.233000000000
    H               37.261000000000    29.548000000000    16.400000000000
    H               38.088000000000    30.532000000000    17.583000000000
    H               35.204150000000    30.416960000000    19.821910000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile36-leu69', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               38.693000000000    34.106000000000    17.595000000000
    H               38.151000000000    33.411000000000    16.950000000000
    C               38.146000000000    33.932000000000    19.027000000000
    H               38.621000000000    34.647000000000    19.698000000000
    H               37.073000000000    34.080000000000    19.020000000000
    H               38.337000000000    32.923000000000    19.390000000000
    C               38.471000000000    35.546000000000    17.045000000000
    H               38.776000000000    36.283000000000    17.789000000000
    H               39.030000000000    35.735000000000    16.131000000000
    C               36.958000000000    35.746000000000    16.680000000000
    H               36.344000000000    35.711000000000    17.578000000000
    H               36.831000000000    36.722000000000    16.212000000000
    H               36.635000000000    34.976000000000    15.978000000000
    H               39.757470000000    33.835190000000    17.535280000000
    --
    0 1
    C               30.925000000000    34.304000000000    17.753000000000
    H               30.790000000000    33.539000000000    16.990000000000
    H               30.827000000000    35.284000000000    17.283000000000
    C               32.345000000000    34.183000000000    18.358000000000
    H               32.492000000000    34.899000000000    19.164000000000
    C               32.555000000000    32.783000000000    18.870000000000
    H               32.403000000000    32.055000000000    18.071000000000
    H               33.579000000000    32.699000000000    19.225000000000
    H               31.891000000000    32.566000000000    19.702000000000
    C               33.361000000000    34.491000000000    17.245000000000
    H               33.262000000000    35.526000000000    16.923000000000
    H               34.372000000000    34.336000000000    17.619000000000
    H               33.201000000000    33.829000000000    16.393000000000
    H               30.134520000000    34.192180000000    18.509730000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile36-leu71', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               38.693000000000    34.106000000000    17.595000000000
    H               38.151000000000    33.411000000000    16.950000000000
    C               38.146000000000    33.932000000000    19.027000000000
    H               38.621000000000    34.647000000000    19.698000000000
    H               37.073000000000    34.080000000000    19.020000000000
    H               38.337000000000    32.923000000000    19.390000000000
    C               38.471000000000    35.546000000000    17.045000000000
    H               38.776000000000    36.283000000000    17.789000000000
    H               39.030000000000    35.735000000000    16.131000000000
    C               36.958000000000    35.746000000000    16.680000000000
    H               36.344000000000    35.711000000000    17.578000000000
    H               36.831000000000    36.722000000000    16.212000000000
    H               36.635000000000    34.976000000000    15.978000000000
    H               39.757470000000    33.835190000000    17.535280000000
    --
    0 1
    C               35.114000000000    36.564000000000    22.907000000000
    H               34.545000000000    37.433000000000    22.572000000000
    H               35.807000000000    36.898000000000    23.679000000000
    C               35.926000000000    35.979000000000    21.737000000000
    H               36.734000000000    35.364000000000    22.133000000000
    C               35.003000000000    35.084000000000    20.920000000000
    H               34.004000000000    35.502000000000    20.834000000000
    H               35.381000000000    35.032000000000    19.904000000000
    H               34.994000000000    34.061000000000    21.285000000000
    C               36.533000000000    37.087000000000    20.917000000000
    H               37.167000000000    37.709000000000    21.550000000000
    H               37.146000000000    36.666000000000    20.120000000000
    H               35.752000000000    37.705000000000    20.474000000000
    H               34.434530000000    35.798280000000    23.309490000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile44-gly47', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               24.639000000000    31.426000000000    21.286000000000
    H               24.663000000000    32.301000000000    20.671000000000
    C               23.181000000000    31.309000000000    21.824000000000
    H               23.022000000000    30.347000000000    22.311000000000
    H               22.976000000000    32.113000000000    22.530000000000
    H               22.472000000000    31.400000000000    21.001000000000
    C               25.646000000000    31.670000000000    22.421000000000
    H               26.626000000000    31.926000000000    22.016000000000
    H               25.304000000000    32.521000000000    23.011000000000
    C               25.778000000000    30.436000000000    23.356000000000
    H               26.403000000000    29.690000000000    22.898000000000
    H               26.269000000000    30.746000000000    24.274000000000
    H               24.808000000000    30.016000000000    23.614000000000
    H               24.912880000000    30.555120000000    20.672360000000
    --
    0 1
    C               18.453000000000    28.941000000000    20.591000000000
    O               17.860000000000    27.994000000000    21.128000000000
    N               19.172000000000    29.808000000000    21.243000000000
    H               19.577000000000    30.558000000000    20.693000000000
    C               19.399000000000    29.894000000000    22.655000000000
    H               20.004000000000    30.779000000000    22.848000000000
    H               18.438000000000    30.033000000000    23.151000000000
    C               20.083000000000    28.729000000000    23.321000000000
    O               19.991000000000    28.584000000000    24.561000000000
    N               20.801000000000    27.931000000000    22.578000000000
    H               20.819000000000    28.069000000000    21.573000000000
    H               18.445770000000    29.087040000000    19.500760000000
    H               21.361950000000    27.080960000000    22.993660000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile44-his68-big', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               24.639000000000    31.426000000000    21.286000000000
    H               24.663000000000    32.301000000000    20.671000000000
    C               23.181000000000    31.309000000000    21.824000000000
    H               23.022000000000    30.347000000000    22.311000000000
    H               22.976000000000    32.113000000000    22.530000000000
    H               22.472000000000    31.400000000000    21.001000000000
    C               25.646000000000    31.670000000000    22.421000000000
    H               26.626000000000    31.926000000000    22.016000000000
    H               25.304000000000    32.521000000000    23.011000000000
    C               25.778000000000    30.436000000000    23.356000000000
    H               26.403000000000    29.690000000000    22.898000000000
    H               26.269000000000    30.746000000000    24.274000000000
    H               24.808000000000    30.016000000000    23.614000000000
    H               24.912880000000    30.555120000000    20.672360000000
    --
    0 1
    C               25.214000000000    34.565000000000    18.780000000000
    H               24.940000000000    33.724000000000    19.379000000000
    H               25.680000000000    35.333000000000    19.400000000000
    C               23.978000000000    35.121000000000    18.126000000000
    N               23.853000000000    36.432000000000    17.781000000000
    C               22.674000000000    36.627000000000    17.200000000000
    H               22.295000000000    37.565000000000    16.823000000000
    N               22.045000000000    35.455000000000    17.173000000000
    H               21.128000000000    35.293000000000    16.781000000000
    C               22.824000000000    34.514000000000    17.782000000000
    H               22.603000000000    33.463000000000    17.912000000000
    H               25.899200000000    34.254000000000    17.977640000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile44-his68-small', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               24.639000000000    31.426000000000    21.286000000000
    H               24.663000000000    32.301000000000    20.671000000000
    C               23.181000000000    31.309000000000    21.824000000000
    H               23.022000000000    30.347000000000    22.311000000000
    H               22.976000000000    32.113000000000    22.530000000000
    H               22.472000000000    31.400000000000    21.001000000000
    C               25.646000000000    31.670000000000    22.421000000000
    H               26.626000000000    31.926000000000    22.016000000000
    H               25.304000000000    32.521000000000    23.011000000000
    C               25.778000000000    30.436000000000    23.356000000000
    H               26.403000000000    29.690000000000    22.898000000000
    H               26.269000000000    30.746000000000    24.274000000000
    H               24.808000000000    30.016000000000    23.614000000000
    H               24.912880000000    30.555120000000    20.672360000000
    --
    0 1
    C               25.214000000000    34.565000000000    18.780000000000
    H               24.940000000000    33.724000000000    19.379000000000
    H               25.680000000000    35.333000000000    19.400000000000
    H               24.310520000000    34.971420000000    18.301940000000
    H               25.899200000000    34.254000000000    17.977640000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile44-his68', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               27.331000000000    29.317000000000    20.364000000000
    O               27.101000000000    28.346000000000    21.097000000000
    N               26.436000000000    30.232000000000    20.004000000000
    H               26.749000000000    31.024000000000    19.452000000000
    C               25.034000000000    30.170000000000    20.401000000000
    H               24.865000000000    29.269000000000    20.985000000000
    C               24.101000000000    30.149000000000    19.196000000000
    O               24.196000000000    30.948000000000    18.287000000000
    N               23.141000000000    29.187000000000    19.241000000000
    H               23.123000000000    28.529000000000    20.017000000000
    H               22.382230000000    29.093560000000    18.450090000000
    H               24.760120000000    31.040880000000    21.014640000000
    H               28.363770000000    29.501760000000    20.033460000000
    --
    0 1
    N               28.525000000000    34.447000000000    18.189000000000
    H               28.437000000000    35.312000000000    17.662000000000
    O               27.507000000000    32.587000000000    18.958000000000
    C               27.475000000000    33.651000000000    18.304000000000
    C               26.179000000000    34.127000000000    17.650000000000
    H               26.376000000000    34.966000000000    16.979000000000
    N               25.621000000000    32.945000000000    16.950000000000
    H               25.216000000000    32.182000000000    17.480000000000
    C               25.698000000000    32.876000000000    15.669000000000
    O               26.158000000000    33.730000000000    14.894000000000
    H               29.486960000000    34.219330000000    18.671490000000
    H               25.306670000000    31.972870000000    15.177870000000
    H               25.493800000000    34.438000000000    18.452360000000

""")

GEOS['%s-%s-%s' % (dbse, 'ile61-leu67', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               22.650000000000    24.516000000000    12.172000000000
    H               23.215000000000    23.598000000000    12.003000000000
    C               23.376000000000    25.645000000000    11.409000000000
    H               22.954000000000    26.611000000000    11.687000000000
    H               24.444000000000    25.631000000000    11.623000000000
    H               23.277000000000    25.511000000000    10.333000000000
    C               22.662000000000    24.819000000000    13.699000000000
    H               22.117000000000    25.741000000000    13.891000000000
    H               22.184000000000    24.007000000000    14.242000000000
    C               24.123000000000    24.981000000000    14.195000000000
    H               24.588000000000    25.879000000000    13.799000000000
    H               24.114000000000    25.061000000000    15.282000000000
    H               24.702000000000    24.100000000000    13.921000000000
    H               21.618800000000    24.338040000000    11.832960000000
    --
    0 1
    C               26.310000000000    30.594000000000    14.967000000000
    H               27.244000000000    31.135000000000    14.803000000000
    H               26.368000000000    30.156000000000    15.965000000000
    C               26.290000000000    29.480000000000    13.960000000000
    H               26.459000000000    29.906000000000    12.971000000000
    C               27.393000000000    28.442000000000    14.229000000000
    H               27.257000000000    27.992000000000    15.212000000000
    H               27.352000000000    27.657000000000    13.474000000000
    H               28.371000000000    28.921000000000    14.181000000000
    C               24.942000000000    28.807000000000    13.952000000000
    H               24.133000000000    29.498000000000    13.737000000000
    H               24.937000000000    28.054000000000    13.165000000000
    H               24.769000000000    28.313000000000    14.905000000000
    H               25.481890000000    31.317980000000    14.976270000000

""")

GEOS['%s-%s-%s' % (dbse, 'leu15-ile30', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               31.562000000000    29.686000000000     8.045000000000
    H               30.790000000000    30.265000000000     8.559000000000
    H               31.073000000000    28.716000000000     7.948000000000
    C               32.631000000000    29.444000000000     9.060000000000
    H               33.028000000000    28.443000000000     8.889000000000
    C               33.814000000000    30.390000000000     9.030000000000
    H               33.465000000000    31.415000000000     9.143000000000
    H               34.503000000000    30.145000000000     9.839000000000
    H               34.339000000000    30.288000000000     8.080000000000
    C               31.945000000000    29.449000000000    10.436000000000
    H               31.086000000000    28.781000000000    10.437000000000
    H               32.653000000000    29.096000000000    11.183000000000
    H               31.642000000000    30.461000000000    10.683000000000
    H               31.644750000000    30.109820000000     7.033303000000
    --
    0 1
    C               35.708000000000    30.776000000000    13.806000000000
    H               35.901000000000    29.999000000000    14.547000000000
    C               35.874000000000    32.138000000000    14.512000000000
    H               35.708000000000    32.951000000000    13.805000000000
    H               35.191000000000    32.233000000000    15.353000000000
    H               36.876000000000    32.231000000000    14.929000000000
    C               34.228000000000    30.630000000000    13.319000000000
    H               33.937000000000    31.481000000000    12.702000000000
    H               34.143000000000    29.718000000000    12.732000000000
    C               33.284000000000    30.504000000000    14.552000000000
    H               33.209000000000    31.442000000000    15.098000000000
    H               32.287000000000    30.236000000000    14.204000000000
    H               33.632000000000    29.710000000000    15.211000000000
    H               36.428860000000    30.630840000000    12.987900000000

""")

GEOS['%s-%s-%s' % (dbse, 'leu15-val26', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               31.562000000000    29.686000000000     8.045000000000
    H               30.790000000000    30.265000000000     8.559000000000
    H               31.073000000000    28.716000000000     7.948000000000
    C               32.631000000000    29.444000000000     9.060000000000
    H               33.028000000000    28.443000000000     8.889000000000
    C               33.814000000000    30.390000000000     9.030000000000
    H               33.465000000000    31.415000000000     9.143000000000
    H               34.503000000000    30.145000000000     9.839000000000
    H               34.339000000000    30.288000000000     8.080000000000
    C               31.945000000000    29.449000000000    10.436000000000
    H               31.086000000000    28.781000000000    10.437000000000
    H               32.653000000000    29.096000000000    11.183000000000
    H               31.642000000000    30.461000000000    10.683000000000
    H               31.644750000000    30.109820000000     7.033303000000
    --
    0 1
    C               32.060000000000    25.257000000000    13.364000000000
    H               31.920000000000    24.919000000000    14.391000000000
    C               31.684000000000    26.749000000000    13.342000000000
    H               31.902000000000    27.168000000000    12.361000000000
    H               30.619000000000    26.858000000000    13.547000000000
    H               32.222000000000    27.303000000000    14.110000000000
    C               31.152000000000    24.421000000000    12.477000000000
    H               31.066000000000    23.409000000000    12.871000000000
    H               30.173000000000    24.880000000000    12.423000000000
    H               31.551000000000    24.370000000000    11.464000000000
    H               33.118250000000    25.142050000000    13.086690000000

""")

GEOS['%s-%s-%s' % (dbse, 'leu43-leu50', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               28.762000000000    29.573000000000    19.906000000000
    H               28.801000000000    30.593000000000    19.529000000000
    C               27.331000000000    29.317000000000    20.364000000000
    O               27.101000000000    28.346000000000    21.097000000000
    N               26.436000000000    30.232000000000    20.004000000000
    H               26.749000000000    31.024000000000    19.452000000000
    H               25.378570000000    30.185240000000    20.303430000000
    H               29.452360000000    29.471520000000    20.756360000000
    H               29.043000000000    28.909870000000    19.074560000000
    --
    0 1
    C               26.826000000000    24.521000000000    21.012000000000
    H               25.983000000000    23.830000000000    21.007000000000
    N               26.465000000000    25.689000000000    21.833000000000
    H               26.766000000000    26.605000000000    21.524000000000
    C               25.743000000000    25.586000000000    22.922000000000
    O               25.325000000000    24.489000000000    23.378000000000
    H               25.459000000000    26.512960000000    23.441700000000
    H               26.981860000000    24.859300000000    19.976980000000
    H               27.671340000000    23.985430000000    21.468690000000

""")

GEOS['%s-%s-%s' % (dbse, 'leu43-leu69', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               29.151000000000    28.655000000000    18.755000000000
    H               29.215000000000    27.648000000000    19.172000000000
    H               28.322000000000    28.637000000000    18.045000000000
    C               30.416000000000    28.912000000000    17.980000000000
    H               31.242000000000    29.080000000000    18.671000000000
    C               30.738000000000    27.693000000000    17.122000000000
    H               29.921000000000    27.482000000000    16.432000000000
    H               31.647000000000    27.884000000000    16.550000000000
    H               30.912000000000    26.825000000000    17.756000000000
    C               30.205000000000    30.168000000000    17.129000000000
    H               30.017000000000    31.034000000000    17.761000000000
    H               31.110000000000    30.356000000000    16.551000000000
    H               29.373000000000    30.026000000000    16.440000000000
    H               28.870000000000    29.318130000000    19.586440000000
    --
    0 1
    C               30.925000000000    34.304000000000    17.753000000000
    H               30.790000000000    33.539000000000    16.990000000000
    H               30.827000000000    35.284000000000    17.283000000000
    C               32.345000000000    34.183000000000    18.358000000000
    H               32.492000000000    34.899000000000    19.164000000000
    C               32.555000000000    32.783000000000    18.870000000000
    H               32.403000000000    32.055000000000    18.071000000000
    H               33.579000000000    32.699000000000    19.225000000000
    H               31.891000000000    32.566000000000    19.702000000000
    C               33.361000000000    34.491000000000    17.245000000000
    H               33.262000000000    35.526000000000    16.923000000000
    H               34.372000000000    34.336000000000    17.619000000000
    H               33.201000000000    33.829000000000    16.393000000000
    H               30.134520000000    34.192180000000    18.509730000000

""")

GEOS['%s-%s-%s' % (dbse, 'leu50-ile61', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               27.043000000000    24.992000000000    19.571000000000
    H               27.973000000000    25.560000000000    19.529000000000
    H               27.174000000000    24.103000000000    18.954000000000
    C               25.931000000000    25.844000000000    18.959000000000
    H               25.888000000000    26.813000000000    19.450000000000
    C               26.203000000000    26.083000000000    17.471000000000
    H               26.174000000000    25.135000000000    16.931000000000
    H               25.450000000000    26.753000000000    17.058000000000
    H               27.188000000000    26.530000000000    17.344000000000
    C               24.577000000000    25.190000000000    19.079000000000
    H               24.281000000000    25.115000000000    20.124000000000
    H               23.829000000000    25.764000000000    18.560000000000
    H               24.614000000000    24.192000000000    18.645000000000
    H               26.887140000000    24.653700000000    20.606020000000
    --
    0 1
    C               22.650000000000    24.516000000000    12.172000000000
    H               23.215000000000    23.598000000000    12.003000000000
    C               23.376000000000    25.645000000000    11.409000000000
    H               22.954000000000    26.611000000000    11.687000000000
    H               24.444000000000    25.631000000000    11.623000000000
    H               23.277000000000    25.511000000000    10.333000000000
    C               22.662000000000    24.819000000000    13.699000000000
    H               22.117000000000    25.741000000000    13.891000000000
    H               22.184000000000    24.007000000000    14.242000000000
    C               24.123000000000    24.981000000000    14.195000000000
    H               24.588000000000    25.879000000000    13.799000000000
    H               24.114000000000    25.061000000000    15.282000000000
    H               24.702000000000    24.100000000000    13.921000000000
    H               21.618800000000    24.338040000000    11.832960000000

""")

GEOS['%s-%s-%s' % (dbse, 'leu50-tyr59', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               27.043000000000    24.992000000000    19.571000000000
    H               27.973000000000    25.560000000000    19.529000000000
    H               27.174000000000    24.103000000000    18.954000000000
    C               25.931000000000    25.844000000000    18.959000000000
    H               25.888000000000    26.813000000000    19.450000000000
    C               26.203000000000    26.083000000000    17.471000000000
    H               26.174000000000    25.135000000000    16.931000000000
    H               25.450000000000    26.753000000000    17.058000000000
    H               27.188000000000    26.530000000000    17.344000000000
    C               24.577000000000    25.190000000000    19.079000000000
    H               24.281000000000    25.115000000000    20.124000000000
    H               23.829000000000    25.764000000000    18.560000000000
    H               24.614000000000    24.192000000000    18.645000000000
    H               26.887140000000    24.653700000000    20.606020000000
    --
    0 1
    C               22.945000000000    21.951000000000    17.785000000000
    C               24.272000000000    21.544000000000    17.644000000000
    H               24.693000000000    21.388000000000    16.659000000000
    C               25.052000000000    21.285000000000    18.776000000000
    H               26.081000000000    20.990000000000    18.667000000000
    C               24.517000000000    21.470000000000    20.030000000000
    O               25.248000000000    21.302000000000    21.191000000000
    H               24.728000000000    21.394000000000    21.993000000000
    C               23.204000000000    21.907000000000    20.192000000000
    H               22.798000000000    22.050000000000    21.180000000000
    C               22.437000000000    22.157000000000    19.065000000000
    H               21.415000000000    22.483000000000    19.173000000000
    H               22.318640000000    22.171680000000    16.908090000000

""")

GEOS['%s-%s-%s' % (dbse, 'leu56-ile61np', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               26.084000000000    21.888000000000    11.833000000000
    H               26.170000000000    21.168000000000    11.029000000000
    H               25.313000000000    22.606000000000    11.554000000000
    C               27.426000000000    22.616000000000    11.902000000000
    H               28.230000000000    21.920000000000    12.129000000000
    C               27.718000000000    23.341000000000    10.578000000000
    H               26.935000000000    24.062000000000    10.350000000000
    H               28.683000000000    23.842000000000    10.643000000000
    H               27.773000000000    22.608000000000     9.776000000000
    C               27.380000000000    23.721000000000    12.955000000000
    H               27.051000000000    23.309000000000    13.905000000000
    H               28.369000000000    24.127000000000    13.111000000000
    H               26.698000000000    24.514000000000    12.653000000000
    H               25.734770000000    21.332790000000    12.716050000000
    --
    0 1
    C               22.650000000000    24.516000000000    12.172000000000
    H               23.215000000000    23.598000000000    12.003000000000
    C               23.376000000000    25.645000000000    11.409000000000
    H               22.954000000000    26.611000000000    11.687000000000
    H               24.444000000000    25.631000000000    11.623000000000
    H               23.277000000000    25.511000000000    10.333000000000
    C               22.662000000000    24.819000000000    13.699000000000
    H               22.117000000000    25.741000000000    13.891000000000
    H               22.184000000000    24.007000000000    14.242000000000
    C               24.123000000000    24.981000000000    14.195000000000
    H               24.588000000000    25.879000000000    13.799000000000
    H               24.114000000000    25.061000000000    15.282000000000
    H               24.702000000000    24.100000000000    13.921000000000
    H               21.618800000000    24.338040000000    11.832960000000

""")

GEOS['%s-%s-%s' % (dbse, 'leu56-ile61p', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               24.240000000000    19.233000000000    12.246000000000
    H               25.095000000000    18.802000000000    11.905000000000
    C               24.241000000000    20.436000000000    12.857000000000
    O               23.264000000000    20.951000000000    13.329000000000
    C               25.594000000000    21.109000000000    13.072000000000
    H               25.497000000000    21.778000000000    13.926000000000
    H               26.334040000000    20.327890000000    13.300510000000
    H               25.943230000000    21.664210000000    12.188950000000
    H               23.264740000000    18.751300000000    12.082220000000
    --
    0 1
    C               19.442000000000    22.745000000000    12.510000000000
    O               18.571000000000    23.610000000000    12.289000000000
    N               20.717000000000    22.964000000000    12.260000000000
    H               21.395000000000    22.241000000000    12.466000000000
    C               21.184000000000    24.263000000000    11.690000000000
    H               20.547000000000    25.089000000000    12.012000000000
    H               19.169840000000    21.739360000000    12.863020000000
    H               21.130670000000    24.153460000000    10.596770000000
    H               22.215200000000    24.440960000000    12.029040000000

""")

GEOS['%s-%s-%s' % (dbse, 'leu56-tyr59', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               25.594000000000    21.109000000000    13.072000000000
    H               25.497000000000    21.778000000000    13.926000000000
    C               24.241000000000    20.436000000000    12.857000000000
    O               23.264000000000    20.951000000000    13.329000000000
    N               24.240000000000    19.233000000000    12.246000000000
    H               25.095000000000    18.802000000000    11.905000000000
    H               23.264740000000    18.751300000000    12.082220000000
    H               25.943230000000    21.664210000000    12.188950000000
    H               26.334040000000    20.327890000000    13.300510000000
    --
    0 1
    C               21.460000000000    18.737000000000    16.163000000000
    O               20.497000000000    18.506000000000    16.900000000000
    N               21.846000000000    19.954000000000    15.905000000000
    H               22.698000000000    20.061000000000    15.367000000000
    C               21.079000000000    21.149000000000    16.251000000000
    H               20.475000000000    20.949000000000    17.136000000000
    H               22.147940000000    17.947810000000    15.825490000000
    H               20.397430000000    21.469780000000    15.449400000000
    H               21.802100000000    21.943260000000    16.488200000000

""")

GEOS['%s-%s-%s' % (dbse, 'leu69-leu71', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               30.925000000000    34.304000000000    17.753000000000
    H               30.790000000000    33.539000000000    16.990000000000
    H               30.827000000000    35.284000000000    17.283000000000
    C               32.345000000000    34.183000000000    18.358000000000
    H               32.492000000000    34.899000000000    19.164000000000
    C               32.555000000000    32.783000000000    18.870000000000
    H               32.403000000000    32.055000000000    18.071000000000
    H               33.579000000000    32.699000000000    19.225000000000
    H               31.891000000000    32.566000000000    19.702000000000
    C               33.361000000000    34.491000000000    17.245000000000
    H               33.262000000000    35.526000000000    16.923000000000
    H               34.372000000000    34.336000000000    17.619000000000
    H               33.201000000000    33.829000000000    16.393000000000
    H               30.134520000000    34.192180000000    18.509730000000
    --
    0 1
    C               35.114000000000    36.564000000000    22.907000000000
    H               34.545000000000    37.433000000000    22.572000000000
    H               35.807000000000    36.898000000000    23.679000000000
    C               35.926000000000    35.979000000000    21.737000000000
    H               36.734000000000    35.364000000000    22.133000000000
    C               35.003000000000    35.084000000000    20.920000000000
    H               34.004000000000    35.502000000000    20.834000000000
    H               35.381000000000    35.032000000000    19.904000000000
    H               34.994000000000    34.061000000000    21.285000000000
    C               36.533000000000    37.087000000000    20.917000000000
    H               37.167000000000    37.709000000000    21.550000000000
    H               37.146000000000    36.666000000000    20.120000000000
    H               35.752000000000    37.705000000000    20.474000000000
    H               34.434530000000    35.798280000000    23.309490000000

""")

GEOS['%s-%s-%s' % (dbse, 'leu71-leu73', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               35.114000000000    36.564000000000    22.907000000000
    H               34.545000000000    37.433000000000    22.572000000000
    H               35.807000000000    36.898000000000    23.679000000000
    C               35.926000000000    35.979000000000    21.737000000000
    H               36.734000000000    35.364000000000    22.133000000000
    C               35.003000000000    35.084000000000    20.920000000000
    H               34.004000000000    35.502000000000    20.834000000000
    H               35.381000000000    35.032000000000    19.904000000000
    H               34.994000000000    34.061000000000    21.285000000000
    C               36.533000000000    37.087000000000    20.917000000000
    H               37.167000000000    37.709000000000    21.550000000000
    H               37.146000000000    36.666000000000    20.120000000000
    H               35.752000000000    37.705000000000    20.474000000000
    H               34.434530000000    35.798280000000    23.309490000000
    --
    0 1
    C               39.080000000000    36.941000000000    27.406000000000
    H               38.272000000000    37.593000000000    27.742000000000
    H               39.928000000000    37.175000000000    28.054000000000
    C               39.502000000000    37.340000000000    26.002000000000
    H               40.549000000000    37.076000000000    25.856000000000
    C               38.684000000000    36.647000000000    24.923000000000
    H               37.628000000000    36.851000000000    25.089000000000
    H               38.971000000000    37.044000000000    23.951000000000
    H               38.853000000000    35.575000000000    24.951000000000
    C               39.337000000000    38.854000000000    25.862000000000
    H               39.924000000000    39.357000000000    26.633000000000
    H               39.708000000000    39.175000000000    24.889000000000
    H               38.290000000000    39.140000000000    25.962000000000
    H               38.782170000000    35.900780000000    27.604070000000

""")

GEOS['%s-%s-%s' % (dbse, 'lys11-glu34', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    1 1
    C               34.762000000000    40.799000000000    11.470000000000
    H               35.174000000000    41.532000000000    12.167000000000
    H               34.864000000000    39.816000000000    11.936000000000
    C               35.614000000000    40.847000000000    10.240000000000
    H               35.727000000000    41.894000000000     9.947000000000
    H               36.609000000000    40.468000000000    10.512000000000
    N               35.100000000000    40.073000000000     9.101000000000
    H               34.148000000000    40.315000000000     8.835000000000
    H               35.703000000000    40.203000000000     8.302000000000
    H               35.178000000000    39.085000000000     9.334000000000
    H               33.695790000000    41.005480000000    11.295180000000
    --
    -1 1
    C               38.290000000000    35.814000000000    12.698000000000
    H               38.041000000000    34.959000000000    13.329000000000
    H               38.325000000000    36.685000000000    13.353000000000
    C               37.156000000000    35.985000000000    11.688000000000
    H               37.250000000000    35.207000000000    10.937000000000
    H               36.216000000000    35.866000000000    12.226000000000
    C               37.192000000000    37.361000000000    11.033000000000
    O               37.519000000000    38.360000000000    11.645000000000
    O               36.861000000000    37.320000000000     9.822000000000
    H               39.277390000000    35.623790000000    12.252040000000

""")

GEOS['%s-%s-%s' % (dbse, 'lys27-asp52', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    1 1
    C               34.509000000000    26.077000000000    19.360000000000
    H               34.787000000000    25.022000000000    19.430000000000
    H               35.306000000000    26.680000000000    19.800000000000
    C               33.206000000000    26.311000000000    20.122000000000
    H               32.920000000000    27.363000000000    20.056000000000
    H               32.425000000000    25.682000000000    19.686000000000
    N               33.455000000000    25.910000000000    21.546000000000
    H               34.179000000000    26.474000000000    21.962000000000
    H               32.614000000000    25.985000000000    22.091000000000
    H               33.744000000000    24.923000000000    21.545000000000
    H               34.388610000000    26.343890000000    18.299680000000
    --
    -1 1
    C               33.638000000000    20.716000000000    21.242000000000
    H               33.550000000000    20.270000000000    22.234000000000
    H               34.364000000000    20.133000000000    20.673000000000
    C               34.174000000000    22.129000000000    21.354000000000
    O               35.252000000000    22.322000000000    21.958000000000
    O               33.544000000000    23.086000000000    20.883000000000
    H               32.666120000000    20.683510000000    20.727810000000

""")

GEOS['%s-%s-%s' % (dbse, 'lys27-gln31', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               35.596000000000    26.715000000000    15.736000000000
    H               35.175000000000    27.721000000000    15.746000000000
    C               36.975000000000    26.826000000000    15.107000000000
    O               37.579000000000    27.926000000000    15.159000000000
    N               37.499000000000    25.743000000000    14.571000000000
    H               37.012000000000    24.858000000000    14.664000000000
    H               38.469410000000    25.756490000000    14.053200000000
    H               34.950810000000    26.046610000000    15.146950000000
    H               35.681600000000    26.346700000000    16.768970000000
    --
    0 1
    C               38.148000000000    30.981000000000    13.069000000000
    O               38.544000000000    32.150000000000    12.856000000000
    N               38.883000000000    30.110000000000    13.713000000000
    H               38.533000000000    29.186000000000    13.937000000000
    C               40.269000000000    30.508000000000    14.115000000000
    H               40.236000000000    31.406000000000    14.731000000000
    H               40.793930000000    29.707250000000    14.656530000000
    H               40.857670000000    30.722580000000    13.210890000000
    H               37.132640000000    30.686490000000    12.765180000000

""")

GEOS['%s-%s-%s' % (dbse, 'lys27-leu43-nonh3', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               35.715000000000    26.203000000000    17.172000000000
    H               35.982000000000    25.146000000000    17.194000000000
    H               36.483000000000    26.785000000000    17.680000000000
    C               34.343000000000    26.445000000000    17.898000000000
    H               34.052000000000    27.490000000000    17.825000000000
    H               33.556000000000    25.841000000000    17.455000000000
    C               34.509000000000    26.077000000000    19.360000000000
    H               34.787000000000    25.022000000000    19.430000000000
    H               35.306000000000    26.680000000000    19.800000000000
    C               33.206000000000    26.311000000000    20.122000000000
    H               32.920000000000    27.363000000000    20.056000000000
    H               32.425000000000    25.682000000000    19.686000000000
    H               35.629400000000    26.571300000000    16.139030000000
    H               33.388580000000    26.016970000000    21.166130000000
    --
    0 1
    C               29.151000000000    28.655000000000    18.755000000000
    H               29.215000000000    27.648000000000    19.172000000000
    H               28.322000000000    28.637000000000    18.045000000000
    C               30.416000000000    28.912000000000    17.980000000000
    H               31.242000000000    29.080000000000    18.671000000000
    C               30.738000000000    27.693000000000    17.122000000000
    H               29.921000000000    27.482000000000    16.432000000000
    H               31.647000000000    27.884000000000    16.550000000000
    H               30.912000000000    26.825000000000    17.756000000000
    C               30.205000000000    30.168000000000    17.129000000000
    H               30.017000000000    31.034000000000    17.761000000000
    H               31.110000000000    30.356000000000    16.551000000000
    H               29.373000000000    30.026000000000    16.440000000000
    H               28.870000000000    29.318130000000    19.586440000000

""")

GEOS['%s-%s-%s' % (dbse, 'lys27-leu43', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    1 1
    C               35.715000000000    26.203000000000    17.172000000000
    H               35.982000000000    25.146000000000    17.194000000000
    H               36.483000000000    26.785000000000    17.680000000000
    C               34.343000000000    26.445000000000    17.898000000000
    H               34.052000000000    27.490000000000    17.825000000000
    H               33.556000000000    25.841000000000    17.455000000000
    C               34.509000000000    26.077000000000    19.360000000000
    H               34.787000000000    25.022000000000    19.430000000000
    H               35.306000000000    26.680000000000    19.800000000000
    C               33.206000000000    26.311000000000    20.122000000000
    H               32.920000000000    27.363000000000    20.056000000000
    H               32.425000000000    25.682000000000    19.686000000000
    N               33.455000000000    25.910000000000    21.546000000000
    H               34.179000000000    26.474000000000    21.962000000000
    H               32.614000000000    25.985000000000    22.091000000000
    H               33.744000000000    24.923000000000    21.545000000000
    H               35.629400000000    26.571300000000    16.139030000000
    --
    0 1
    C               29.151000000000    28.655000000000    18.755000000000
    H               29.215000000000    27.648000000000    19.172000000000
    H               28.322000000000    28.637000000000    18.045000000000
    C               30.416000000000    28.912000000000    17.980000000000
    H               31.242000000000    29.080000000000    18.671000000000
    C               30.738000000000    27.693000000000    17.122000000000
    H               29.921000000000    27.482000000000    16.432000000000
    H               31.647000000000    27.884000000000    16.550000000000
    H               30.912000000000    26.825000000000    17.756000000000
    C               30.205000000000    30.168000000000    17.129000000000
    H               30.017000000000    31.034000000000    17.761000000000
    H               31.110000000000    30.356000000000    16.551000000000
    H               29.373000000000    30.026000000000    16.440000000000
    H               28.870000000000    29.318130000000    19.586440000000

""")

GEOS['%s-%s-%s' % (dbse, 'lys29-lys33', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               37.471000000000    27.391000000000    10.668000000000
    H               38.331000000000    27.243000000000    10.009000000000
    C               37.441000000000    28.882000000000    11.052000000000
    O               38.020000000000    29.772000000000    10.382000000000
    N               36.811000000000    29.170000000000    12.192000000000
    H               36.331000000000    28.429000000000    12.692000000000
    H               36.751280000000    30.215030000000    12.530140000000
    H               37.591610000000    26.759650000000    11.560670000000
    H               36.547490000000    27.150370000000    10.120980000000
    --
    0 1
    C               41.399000000000    31.338000000000     9.967000000000
    O               42.260000000000    32.036000000000     9.381000000000
    N               40.117000000000    31.750000000000     9.988000000000
    H               39.395000000000    31.145000000000    10.363000000000
    C               39.808000000000    32.994000000000     9.233000000000
    H               40.603000000000    33.161000000000     8.505000000000
    H               39.829450000000    33.938430000000     9.796552000000
    H               38.941580000000    32.853830000000     8.569927000000
    H               41.630850000000    30.381520000000    10.458320000000

""")

GEOS['%s-%s-%s' % (dbse, 'lys6-thr66', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               26.219000000000    37.684000000000    14.307000000000
    H               25.721000000000    37.423000000000    13.371000000000
    H               25.789000000000    37.056000000000    15.091000000000
    C               25.884000000000    39.139000000000    14.615000000000
    H               26.301000000000    39.420000000000    15.583000000000
    H               26.304000000000    39.777000000000    13.837000000000
    C               24.348000000000    39.296000000000    14.642000000000
    H               23.937000000000    38.867000000000    13.726000000000
    H               23.971000000000    38.722000000000    15.490000000000
    H               23.996220000000    40.335320000000    14.719930000000
    H               27.279810000000    37.418080000000    14.188810000000
    --
    0 1
    C               22.699000000000    34.267000000000    11.985000000000
    H               21.747000000000    34.305000000000    12.516000000000
    C               23.727000000000    35.131000000000    12.722000000000
    H               24.714000000000    34.993000000000    12.277000000000
    H               23.447000000000    36.180000000000    12.638000000000
    H               23.766000000000    34.863000000000    13.778000000000
    O               22.495000000000    34.690000000000    10.589000000000
    H               22.365000000000    35.639000000000    10.549000000000
    H               23.053280000000    33.227640000000    11.920080000000

""")

GEOS['%s-%s-%s' % (dbse, 'met1-ile3', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               25.353000000000    24.860000000000     5.134000000000
    H               26.271000000000    24.325000000000     5.365000000000
    H               25.400000000000    25.873000000000     5.533000000000
    S               23.930000000000    23.959000000000     5.904000000000
    C               24.447000000000    23.984000000000     7.620000000000
    H               24.441000000000    25.008000000000     7.993000000000
    H               23.766000000000    23.377000000000     8.211000000000
    H               25.451000000000    23.568000000000     7.707000000000
    H               25.176800000000    24.874620000000     4.048302000000
    --
    0 1
    C               26.344000000000    29.050000000000     8.645000000000
    H               25.879000000000    29.504000000000     9.522000000000
    C               25.491000000000    27.771000000000     8.287000000000
    H               25.945000000000    27.245000000000     7.447000000000
    H               25.448000000000    27.099000000000     9.144000000000
    H               24.475000000000    28.063000000000     8.019000000000
    C               27.810000000000    28.748000000000     8.999000000000
    H               28.241000000000    28.093000000000     8.242000000000
    H               28.394000000000    29.667000000000     9.015000000000
    C               27.967000000000    28.087000000000    10.417000000000
    H               27.385000000000    27.167000000000    10.473000000000
    H               29.013000000000    27.846000000000    10.595000000000
    H               27.625000000000    28.782000000000    11.183000000000
    H               26.265720000000    29.773940000000     7.820512000000

""")

GEOS['%s-%s-%s' % (dbse, 'met1-lys63', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               25.112000000000    24.880000000000     3.649000000000
    H               24.278000000000    25.546000000000     3.464000000000
    H               24.837000000000    23.888000000000     3.296000000000
    C               25.353000000000    24.860000000000     5.134000000000
    H               26.271000000000    24.325000000000     5.365000000000
    H               25.400000000000    25.873000000000     5.533000000000
    S               23.930000000000    23.959000000000     5.904000000000
    C               24.447000000000    23.984000000000     7.620000000000
    H               24.441000000000    25.008000000000     7.993000000000
    H               23.766000000000    23.377000000000     8.211000000000
    H               25.451000000000    23.568000000000     7.707000000000
    H               25.955080000000    25.269390000000     3.059431000000
    --
    0 1
    O               21.323000000000    26.830000000000     8.008000000000
    C               20.822000000000    25.914000000000     7.332000000000
    N               20.924000000000    25.862000000000     6.006000000000
    H               20.502000000000    25.078000000000     5.521000000000
    C               21.656000000000    26.847000000000     5.240000000000
    H               22.699000000000    26.813000000000     5.562000000000
    C               21.631000000000    26.642000000000     3.731000000000
    H               22.082000000000    27.514000000000     3.259000000000
    H               22.206000000000    25.755000000000     3.485000000000
    C               21.127000000000    28.240000000000     5.574000000000
    O               19.958000000000    28.465000000000     5.842000000000
    N               22.099000000000    29.163000000000     5.605000000000
    H               23.043000000000    28.873000000000     5.376000000000
    H               20.289420000000    25.093930000000     7.835830000000
    H               21.952310000000    30.232590000000     5.815862000000
    H               20.617010000000    26.485730000000     3.334254000000

""")

GEOS['%s-%s-%s' % (dbse, 'met1-val17-bi', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    1 1
    N               27.340000000000    24.430000000000     2.614000000000
    H               28.235000000000    24.895000000000     2.553000000000
    H               27.193000000000    23.902000000000     1.767000000000
    H               27.397000000000    23.781000000000     3.387000000000
    C               26.266000000000    25.413000000000     2.842000000000
    H               25.897000000000    25.738000000000     1.870000000000
    C               26.913000000000    26.639000000000     3.531000000000
    O               27.886000000000    26.463000000000     4.263000000000
    N               26.335000000000    27.770000000000     3.258000000000
    H               25.505000000000    27.818000000000     2.673000000000
    H               25.422920000000    25.023610000000     3.431569000000
    H               26.713520000000    28.689480000000     3.728398000000
    --
    0 1
    C               31.440000000000    26.079000000000     5.080000000000
    O               32.576000000000    25.802000000000     5.461000000000
    N               30.310000000000    25.458000000000     5.384000000000
    H               29.430000000000    25.804000000000     5.018000000000
    C               30.288000000000    24.245000000000     6.193000000000
    H               31.271000000000    23.769000000000     6.195000000000
    C               29.279000000000    23.227000000000     5.641000000000
    O               28.478000000000    23.522000000000     4.725000000000
    N               29.380000000000    22.057000000000     6.232000000000
    H               30.114000000000    21.905000000000     6.927000000000
    H               30.016550000000    24.488250000000     7.230856000000
    H               28.694700000000    21.217650000000     6.042640000000
    H               31.280050000000    26.996530000000     4.494726000000

""")

GEOS['%s-%s-%s' % (dbse, 'met1-val17-mono', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               26.913000000000    26.639000000000     3.531000000000
    O               27.886000000000    26.463000000000     4.263000000000
    N               26.335000000000    27.770000000000     3.258000000000
    H               25.505000000000    27.818000000000     2.673000000000
    H               26.453260000000    25.767830000000     3.041411000000
    H               26.713520000000    28.689480000000     3.728398000000
    --
    0 1
    C               31.440000000000    26.079000000000     5.080000000000
    O               32.576000000000    25.802000000000     5.461000000000
    N               30.310000000000    25.458000000000     5.384000000000
    H               29.430000000000    25.804000000000     5.018000000000
    H               31.280050000000    26.996530000000     4.494726000000
    H               30.293400000000    24.542960000000     5.994275000000

""")

GEOS['%s-%s-%s' % (dbse, 'met1-val17', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               25.353000000000    24.860000000000     5.134000000000
    H               26.271000000000    24.325000000000     5.365000000000
    H               25.400000000000    25.873000000000     5.533000000000
    S               23.930000000000    23.959000000000     5.904000000000
    C               24.447000000000    23.984000000000     7.620000000000
    H               24.441000000000    25.008000000000     7.993000000000
    H               23.766000000000    23.377000000000     8.211000000000
    H               25.451000000000    23.568000000000     7.707000000000
    H               25.176800000000    24.874620000000     4.048302000000
    --
    0 1
    C               29.903000000000    24.590000000000     7.665000000000
    H               29.900000000000    23.663000000000     8.237000000000
    C               30.862000000000    25.496000000000     8.389000000000
    H               30.895000000000    26.467000000000     7.902000000000
    H               30.547000000000    25.617000000000     9.425000000000
    H               31.853000000000    25.045000000000     8.370000000000
    C               28.476000000000    25.135000000000     7.705000000000
    H               27.762000000000    24.367000000000     7.412000000000
    H               28.235000000000    25.443000000000     8.721000000000
    H               28.371000000000    26.000000000000     7.051000000000
    H               30.174450000000    24.346750000000     6.627144000000

""")

GEOS['%s-%s-%s' % (dbse, 'phe4-leu67', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               26.772000000000    33.436000000000     9.197000000000
    H               27.659000000000    33.695000000000     8.619000000000
    C               27.151000000000    33.362000000000    10.650000000000
    O               26.350000000000    32.778000000000    11.395000000000
    N               28.260000000000    33.943000000000    11.096000000000
    H               28.901000000000    34.360000000000    10.429000000000
    H               26.366010000000    32.461780000000     8.887053000000
    H               25.999310000000    34.197930000000     9.016922000000
    H               28.521930000000    33.959700000000    12.164230000000
    --
    0 1
    C               23.509000000000    32.224000000000    13.290000000000
    O               22.544000000000    31.942000000000    14.034000000000
    N               24.790000000000    32.021000000000    13.618000000000
    H               25.521000000000    32.273000000000    12.959000000000
    C               25.149000000000    31.609000000000    14.980000000000
    H               24.299000000000    31.206000000000    15.534000000000
    H               23.295190000000    32.611300000000    12.282880000000
    H               25.977110000000    30.885020000000    14.970730000000
    H               25.540330000000    32.512130000000    15.471130000000

""")

GEOS['%s-%s-%s' % (dbse, 'phe4-ser65', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               26.882000000000    31.428000000000     7.862000000000
    O               27.906000000000    31.711000000000     7.264000000000
    N               26.214000000000    32.097000000000     8.771000000000
    H               25.335000000000    31.759000000000     9.150000000000
    C               26.772000000000    33.436000000000     9.197000000000
    H               27.659000000000    33.695000000000     8.619000000000
    H               25.999310000000    34.197930000000     9.016922000000
    H               27.049300000000    33.381860000000    10.260100000000
    H               26.425330000000    30.461010000000     7.604371000000
    --
    0 1
    C               21.419000000000    30.253000000000     9.620000000000
    H               20.418000000000    30.649000000000     9.793000000000
    C               22.504000000000    31.228000000000    10.136000000000
    O               23.579000000000    31.321000000000     9.554000000000
    N               22.241000000000    31.873000000000    11.241000000000
    H               21.334000000000    31.762000000000    11.658000000000
    H               22.968480000000    32.539050000000    11.727980000000
    H               21.610080000000    30.088900000000     8.549224000000
    H               21.575300000000    29.299400000000    10.145560000000

""")

GEOS['%s-%s-%s' % (dbse, 'phe4-thr12', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               26.882000000000    31.428000000000     7.862000000000
    O               27.906000000000    31.711000000000     7.264000000000
    N               26.214000000000    32.097000000000     8.771000000000
    H               25.335000000000    31.759000000000     9.150000000000
    C               26.772000000000    33.436000000000     9.197000000000
    H               27.659000000000    33.695000000000     8.619000000000
    C               25.695000000000    34.498000000000     8.946000000000
    H               24.818000000000    34.250000000000     9.540000000000
    H               26.056000000000    35.461000000000     9.292000000000
    C               27.151000000000    33.362000000000    10.650000000000
    O               26.350000000000    32.778000000000    11.395000000000
    N               28.260000000000    33.943000000000    11.096000000000
    H               28.901000000000    34.360000000000    10.429000000000
    H               26.425330000000    30.461010000000     7.604371000000
    H               28.521930000000    33.959700000000    12.164230000000
    H               25.397970000000    34.579010000000     7.889965000000
    --
    0 1
    C               28.113000000000    39.049000000000    10.015000000000
    H               28.132000000000    39.632000000000     9.094000000000
    C               27.588000000000    37.635000000000     9.715000000000
    H               27.589000000000    37.035000000000    10.624000000000
    H               26.569000000000    37.716000000000     9.340000000000
    H               28.195000000000    37.151000000000     8.950000000000
    O               27.280000000000    39.722000000000    10.996000000000
    H               27.694000000000    40.576000000000    11.171000000000
    H               29.117270000000    39.028620000000    10.463370000000

""")

GEOS['%s-%s-%s' % (dbse, 'phe4-thr14', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               25.695000000000    34.498000000000     8.946000000000
    H               24.818000000000    34.250000000000     9.540000000000
    H               26.056000000000    35.461000000000     9.292000000000
    C               25.288000000000    34.609000000000     7.499000000000
    C               24.147000000000    33.966000000000     7.038000000000
    H               23.517000000000    33.410000000000     7.721000000000
    C               23.812000000000    34.031000000000     5.677000000000
    H               22.916000000000    33.556000000000     5.304000000000
    C               24.620000000000    34.778000000000     4.853000000000
    H               24.348000000000    34.851000000000     3.809000000000
    C               25.810000000000    35.392000000000     5.267000000000
    H               26.421000000000    35.939000000000     4.563000000000
    C               26.136000000000    35.346000000000     6.640000000000
    H               26.994000000000    35.858000000000     7.039000000000
    H               26.467690000000    33.736080000000     9.126078000000
    --
    0 1
    C               30.091000000000    34.393000000000     5.078000000000
    H               29.482000000000    33.653000000000     4.557000000000
    C               29.420000000000    35.756000000000     5.119000000000
    H               30.164000000000    36.545000000000     5.227000000000
    H               28.874000000000    35.909000000000     4.188000000000
    H               28.740000000000    35.800000000000     5.956000000000
    O               31.440000000000    34.513000000000     4.487000000000
    H               31.819000000000    33.626000000000     4.434000000000
    H               30.379780000000    34.037950000000     6.078273000000

""")

GEOS['%s-%s-%s' % (dbse, 'phe45-ala46', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               21.447000000000    27.869000000000    16.026000000000
    C               21.325000000000    28.813000000000    15.005000000000
    H               21.969000000000    29.678000000000    14.997000000000
    C               20.369000000000    28.648000000000    14.001000000000
    H               20.257000000000    29.367000000000    13.204000000000
    C               19.593000000000    27.465000000000    14.021000000000
    H               18.874000000000    27.293000000000    13.231000000000
    C               19.677000000000    26.539000000000    15.051000000000
    H               19.057000000000    25.651000000000    15.039000000000
    C               20.638000000000    26.735000000000    16.053000000000
    H               20.746000000000    25.999000000000    16.838000000000
    H               22.205680000000    28.005230000000    16.810760000000
    --
    0 1
    C               17.864000000000    27.977000000000    18.346000000000
    H               18.428000000000    27.068000000000    18.559000000000
    H               16.825000000000    27.830000000000    18.643000000000
    H               17.893000000000    28.164000000000    17.271000000000
    H               18.289740000000    28.834360000000    18.887920000000

""")

GEOS['%s-%s-%s' % (dbse, 'phe45-ile61', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               22.494000000000    28.057000000000    17.109000000000
    H               23.431000000000    28.368000000000    16.652000000000
    H               22.654000000000    27.097000000000    17.579000000000
    C               21.447000000000    27.869000000000    16.026000000000
    C               21.325000000000    28.813000000000    15.005000000000
    H               21.969000000000    29.678000000000    14.997000000000
    C               20.369000000000    28.648000000000    14.001000000000
    H               20.257000000000    29.367000000000    13.204000000000
    C               19.593000000000    27.465000000000    14.021000000000
    H               18.874000000000    27.293000000000    13.231000000000
    C               19.677000000000    26.539000000000    15.051000000000
    H               19.057000000000    25.651000000000    15.039000000000
    C               20.638000000000    26.735000000000    16.053000000000
    H               20.746000000000    25.999000000000    16.838000000000
    H               22.227020000000    28.786110000000    17.888170000000
    --
    0 1
    C               22.650000000000    24.516000000000    12.172000000000
    H               23.215000000000    23.598000000000    12.003000000000
    C               23.376000000000    25.645000000000    11.409000000000
    H               22.954000000000    26.611000000000    11.687000000000
    H               24.444000000000    25.631000000000    11.623000000000
    H               23.277000000000    25.511000000000    10.333000000000
    C               22.662000000000    24.819000000000    13.699000000000
    H               22.117000000000    25.741000000000    13.891000000000
    H               22.184000000000    24.007000000000    14.242000000000
    C               24.123000000000    24.981000000000    14.195000000000
    H               24.588000000000    25.879000000000    13.799000000000
    H               24.114000000000    25.061000000000    15.282000000000
    H               24.702000000000    24.100000000000    13.921000000000
    H               21.618800000000    24.338040000000    11.832960000000

""")

GEOS['%s-%s-%s' % (dbse, 'phe45-leu67', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               22.494000000000    28.057000000000    17.109000000000
    H               23.431000000000    28.368000000000    16.652000000000
    H               22.654000000000    27.097000000000    17.579000000000
    C               21.447000000000    27.869000000000    16.026000000000
    C               21.325000000000    28.813000000000    15.005000000000
    H               21.969000000000    29.678000000000    14.997000000000
    C               20.369000000000    28.648000000000    14.001000000000
    H               20.257000000000    29.367000000000    13.204000000000
    C               19.593000000000    27.465000000000    14.021000000000
    H               18.874000000000    27.293000000000    13.231000000000
    C               19.677000000000    26.539000000000    15.051000000000
    H               19.057000000000    25.651000000000    15.039000000000
    C               20.638000000000    26.735000000000    16.053000000000
    H               20.746000000000    25.999000000000    16.838000000000
    H               22.227020000000    28.786110000000    17.888170000000
    --
    0 1
    C               26.310000000000    30.594000000000    14.967000000000
    H               27.244000000000    31.135000000000    14.803000000000
    H               26.368000000000    30.156000000000    15.965000000000
    C               26.290000000000    29.480000000000    13.960000000000
    H               26.459000000000    29.906000000000    12.971000000000
    C               27.393000000000    28.442000000000    14.229000000000
    H               27.257000000000    27.992000000000    15.212000000000
    H               27.352000000000    27.657000000000    13.474000000000
    H               28.371000000000    28.921000000000    14.181000000000
    C               24.942000000000    28.807000000000    13.952000000000
    H               24.133000000000    29.498000000000    13.737000000000
    H               24.937000000000    28.054000000000    13.165000000000
    H               24.769000000000    28.313000000000    14.905000000000
    H               25.481890000000    31.317980000000    14.976270000000

""")

GEOS['%s-%s-%s' % (dbse, 'phe45-lys48', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               19.810000000000    29.378000000000    18.578000000000
    H               19.931000000000    30.108000000000    17.886000000000
    O               20.821000000000    27.734000000000    19.749000000000
    C               20.835000000000    28.629000000000    18.904000000000
    C               22.126000000000    29.062000000000    18.183000000000
    H               21.966000000000    30.031000000000    17.708000000000
    N               23.141000000000    29.187000000000    19.241000000000
    H               23.123000000000    28.529000000000    20.017000000000
    C               24.101000000000    30.149000000000    19.196000000000
    O               24.196000000000    30.948000000000    18.287000000000
    H               24.774370000000    30.164160000000    20.065680000000
    H               18.791320000000    29.202880000000    18.954320000000
    H               22.392980000000    28.332890000000    17.403830000000
    --
    0 1
    N               23.880000000000    26.727000000000    23.851000000000
    H               23.552000000000    26.164000000000    24.620000000000
    O               23.383000000000    27.627000000000    21.870000000000
    C               23.046000000000    27.087000000000    22.913000000000
    C               21.550000000000    26.796000000000    23.133000000000
    H               21.355000000000    26.660000000000    24.198000000000
    N               20.801000000000    27.931000000000    22.578000000000
    H               20.819000000000    28.069000000000    21.573000000000
    C               20.083000000000    28.729000000000    23.321000000000
    O               19.991000000000    28.584000000000    24.561000000000
    H               19.583460000000    29.579820000000    22.834610000000
    H               21.325440000000    25.864930000000    22.592000000000
    H               24.963970000000    26.833990000000    23.697520000000

""")

GEOS['%s-%s-%s' % (dbse, 'pro19-ser57', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               27.819000000000    20.609000000000     7.316000000000
    O               28.449000000000    20.674000000000     8.360000000000
    N               26.559000000000    20.220000000000     7.288000000000
    C               25.829000000000    19.825000000000     8.494000000000
    H               25.680000000000    20.700000000000     9.116000000000
    C               26.541000000000    18.732000000000     9.251000000000
    O               26.333000000000    18.536000000000    10.457000000000
    N               27.361000000000    17.959000000000     8.559000000000
    H               27.523000000000    18.111000000000     7.575000000000
    H               27.878790000000    17.119180000000     9.045408000000
    H               28.288140000000    20.848270000000     6.350256000000
    H               25.932730000000    20.136990000000     6.387503000000
    H               24.860590000000    19.473950000000     8.108061000000
    --
    0 1
    C               24.241000000000    20.436000000000    12.857000000000
    O               23.264000000000    20.951000000000    13.329000000000
    N               24.240000000000    19.233000000000    12.246000000000
    H               25.095000000000    18.802000000000    11.905000000000
    C               22.924000000000    18.583000000000    12.025000000000
    H               22.278000000000    19.280000000000    11.494000000000
    C               22.229000000000    18.244000000000    13.325000000000
    O               20.963000000000    18.253000000000    13.395000000000
    N               22.997000000000    17.978000000000    14.366000000000
    H               24.005000000000    17.928000000000    14.255000000000
    C               23.059000000000    17.326000000000    11.154000000000
    H               22.077000000000    16.875000000000    11.002000000000
    H               23.472000000000    17.606000000000    10.185000000000
    O               23.914000000000    16.395000000000    11.755000000000
    H               23.376000000000    15.700000000000    12.202000000000
    H               25.216070000000    20.921010000000    13.011940000000
    H               22.568740000000    17.726520000000    15.347510000000

""")

GEOS['%s-%s-%s' % (dbse, 'pro37-gln40', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               41.461000000000    30.751000000000    19.594000000000
    H               42.130000000000    30.205000000000    18.936000000000
    C               40.168000000000    30.026000000000    19.918000000000
    O               39.264000000000    30.662000000000    20.521000000000
    N               40.059000000000    28.758000000000    19.607000000000
    H               40.782270000000    28.128850000000    19.067510000000
    H               39.130860000000    28.206500000000    19.817740000000
    H               41.257910000000    31.747020000000    19.173640000000
    H               41.978800000000    31.026830000000    20.524480000000
    --
    0 1
    C               37.738000000000    31.637000000000    23.712000000000
    H               37.651000000000    31.725000000000    24.796000000000
    N               38.419000000000    30.373000000000    23.341000000000
    H               38.981000000000    30.319000000000    22.499000000000
    C               38.365000000000    29.335000000000    24.159000000000
    O               37.684000000000    29.390000000000    25.221000000000
    H               38.869030000000    28.416480000000    23.823940000000
    H               36.735410000000    31.711980000000    23.265690000000
    H               38.300570000000    32.503640000000    23.334580000000

""")

GEOS['%s-%s-%s' % (dbse, 'ser57-asn60', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               22.924000000000    18.583000000000    12.025000000000
    H               22.278000000000    19.280000000000    11.494000000000
    C               22.229000000000    18.244000000000    13.325000000000
    O               20.963000000000    18.253000000000    13.395000000000
    N               22.997000000000    17.978000000000    14.366000000000
    H               24.005000000000    17.928000000000    14.255000000000
    H               22.568740000000    17.726520000000    15.347510000000
    H               23.020730000000    17.682350000000    11.400920000000
    H               23.899260000000    19.064700000000    12.188780000000
    --
    0 1
    C               20.142000000000    21.590000000000    15.149000000000
    O               19.499000000000    22.645000000000    15.321000000000
    N               19.993000000000    20.884000000000    14.049000000000
    H               20.459000000000    19.985000000000    13.959000000000
    C               19.065000000000    21.352000000000    12.999000000000
    H               19.219000000000    20.677000000000    12.159000000000
    H               18.016100000000    21.302360000000    13.326650000000
    H               20.823570000000    21.269220000000    15.950590000000
    H               19.337160000000    22.357640000000    12.645980000000

""")

GEOS['%s-%s-%s' % (dbse, 'thr14-lys33', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               30.955000000000    35.211000000000     8.459000000000
    O               30.025000000000    34.618000000000     9.040000000000
    N               31.244000000000    34.986000000000     7.197000000000
    H               31.951000000000    35.496000000000     6.692000000000
    C               30.505000000000    33.884000000000     6.512000000000
    H               29.587000000000    33.607000000000     7.028000000000
    C               30.091000000000    34.393000000000     5.078000000000
    H               29.482000000000    33.653000000000     4.557000000000
    C               29.420000000000    35.756000000000     5.119000000000
    H               30.164000000000    36.545000000000     5.227000000000
    H               28.874000000000    35.909000000000     4.188000000000
    H               28.740000000000    35.800000000000     5.956000000000
    O               31.440000000000    34.513000000000     4.487000000000
    H               31.819000000000    33.626000000000     4.434000000000
    C               31.409000000000    32.680000000000     6.446000000000
    O               32.619000000000    32.812000000000     6.125000000000
    N               30.884000000000    31.485000000000     6.666000000000
    H               29.884000000000    31.395000000000     6.821000000000
    H               31.486850000000    30.565140000000     6.645474000000
    H               31.514580000000    35.999540000000     8.983472000000
    --
    1 1
    C               36.351000000000    33.613000000000     7.838000000000
    H               35.342000000000    33.632000000000     8.206000000000
    H               36.727000000000    34.638000000000     7.774000000000
    C               36.322000000000    32.944000000000     6.477000000000
    H               37.337000000000    32.652000000000     6.186000000000
    H               35.709000000000    32.038000000000     6.517000000000
    N               35.768000000000    33.945000000000     5.489000000000
    H               36.271000000000    34.819000000000     5.550000000000
    H               35.864000000000    33.572000000000     4.552000000000
    H               34.779000000000    34.081000000000     5.667000000000
    H               36.973380000000    33.046490000000     8.546318000000

""")

GEOS['%s-%s-%s' % (dbse, 'thr22-asn25', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               32.492000000000    18.193000000000    14.995000000000
    H               32.417000000000    18.313000000000    16.076000000000
    C               32.352000000000    16.700000000000    14.630000000000
    H               32.515000000000    16.546000000000    13.564000000000
    H               33.089000000000    16.126000000000    15.191000000000
    H               31.355000000000    16.351000000000    14.901000000000
    O               33.778000000000    18.739000000000    14.516000000000
    H               34.474000000000    18.249000000000    14.975000000000
    H               31.724450000000    18.804090000000    14.497570000000
    --
    0 1
    C               35.615000000000    22.190000000000    15.759000000000
    O               36.532000000000    23.046000000000    15.724000000000
    N               35.139000000000    21.624000000000    14.662000000000
    H               34.436000000000    20.901000000000    14.762000000000
    C               35.590000000000    21.945000000000    13.302000000000
    H               36.681000000000    21.885000000000    13.291000000000
    H               35.212480000000    21.235900000000    12.550560000000
    H               35.336600000000    22.979490000000    13.027000000000
    H               35.189200000000    21.848780000000    16.714120000000

""")

GEOS['%s-%s-%s' % (dbse, 'thr22-thr55-big', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               30.491000000000    19.162000000000    12.040000000000
    O               29.367000000000    19.523000000000    12.441000000000
    N               31.510000000000    18.936000000000    12.852000000000
    H               32.418000000000    18.706000000000    12.477000000000
    C               31.398000000000    19.064000000000    14.286000000000
    H               30.436000000000    18.719000000000    14.622000000000
    C               32.492000000000    18.193000000000    14.995000000000
    H               32.417000000000    18.313000000000    16.076000000000
    C               32.352000000000    16.700000000000    14.630000000000
    H               32.515000000000    16.546000000000    13.564000000000
    H               33.089000000000    16.126000000000    15.191000000000
    H               31.355000000000    16.351000000000    14.901000000000
    O               33.778000000000    18.739000000000    14.516000000000
    H               34.474000000000    18.249000000000    14.975000000000
    C               31.593000000000    20.553000000000    14.655000000000
    O               32.159000000000    21.311000000000    13.861000000000
    N               31.113000000000    20.863000000000    15.860000000000
    H               30.636000000000    20.148000000000    16.405000000000
    H               31.244870000000    21.871200000000    16.279710000000
    H               30.713580000000    19.104350000000    10.964300000000
    --
    0 1
    C               28.375000000000    17.999000000000    16.887000000000
    O               29.326000000000    18.786000000000    16.690000000000
    N               27.510000000000    17.689000000000    15.954000000000
    H               26.729000000000    17.079000000000    16.191000000000
    C               27.574000000000    18.192000000000    14.563000000000
    H               28.529000000000    18.639000000000    14.366000000000
    C               27.299000000000    17.055000000000    13.533000000000
    H               27.287000000000    17.474000000000    12.529000000000
    C               28.236000000000    15.864000000000    13.558000000000
    H               28.192000000000    15.374000000000    14.531000000000
    H               27.945000000000    15.155000000000    12.783000000000
    H               29.254000000000    16.201000000000    13.365000000000
    O               25.925000000000    16.611000000000    13.913000000000
    H               25.537000000000    16.148000000000    13.159000000000
    C               26.482000000000    19.280000000000    14.432000000000
    O               25.609000000000    19.388000000000    15.287000000000
    N               26.585000000000    20.063000000000    13.378000000000
    H               27.373000000000    19.916000000000    12.751000000000
    H               28.181940000000    17.594070000000    17.891370000000
    H               25.844960000000    20.844110000000    13.149490000000

""")

GEOS['%s-%s-%s' % (dbse, 'thr22-thr55-small', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    N               31.510000000000    18.936000000000    12.852000000000
    H               32.418000000000    18.706000000000    12.477000000000
    C               31.398000000000    19.064000000000    14.286000000000
    H               30.436000000000    18.719000000000    14.622000000000
    C               32.492000000000    18.193000000000    14.995000000000
    H               32.417000000000    18.313000000000    16.076000000000
    C               32.352000000000    16.700000000000    14.630000000000
    H               32.515000000000    16.546000000000    13.564000000000
    H               33.089000000000    16.126000000000    15.191000000000
    H               31.355000000000    16.351000000000    14.901000000000
    O               33.778000000000    18.739000000000    14.516000000000
    H               34.474000000000    18.249000000000    14.975000000000
    C               31.593000000000    20.553000000000    14.655000000000
    O               32.159000000000    21.311000000000    13.861000000000
    H               31.197080000000    20.808700000000    15.648920000000
    H               30.662380000000    19.123990000000    12.176570000000
    --
    0 1
    N               27.510000000000    17.689000000000    15.954000000000
    H               26.729000000000    17.079000000000    16.191000000000
    C               27.574000000000    18.192000000000    14.563000000000
    H               28.529000000000    18.639000000000    14.366000000000
    C               27.299000000000    17.055000000000    13.533000000000
    H               27.287000000000    17.474000000000    12.529000000000
    C               28.236000000000    15.864000000000    13.558000000000
    H               28.192000000000    15.374000000000    14.531000000000
    H               27.945000000000    15.155000000000    12.783000000000
    H               29.254000000000    16.201000000000    13.365000000000
    O               25.925000000000    16.611000000000    13.913000000000
    H               25.537000000000    16.148000000000    13.159000000000
    C               26.482000000000    19.280000000000    14.432000000000
    O               25.609000000000    19.388000000000    15.287000000000
    H               28.236610000000    17.949400000000    16.737730000000
    H               26.568030000000    19.933960000000    13.551700000000

""")

GEOS['%s-%s-%s' % (dbse, 'thr22-val26', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               31.398000000000    19.064000000000    14.286000000000
    H               30.436000000000    18.719000000000    14.622000000000
    C               31.593000000000    20.553000000000    14.655000000000
    O               32.159000000000    21.311000000000    13.861000000000
    N               31.113000000000    20.863000000000    15.860000000000
    H               30.636000000000    20.148000000000    16.405000000000
    H               31.244870000000    21.871200000000    16.279710000000
    H               32.165550000000    18.452910000000    14.783430000000
    H               31.483320000000    18.966500000000    13.193660000000
    --
    0 1
    C               35.238000000000    23.382000000000    12.920000000000
    O               36.066000000000    24.109000000000    12.333000000000
    N               34.007000000000    23.745000000000    13.250000000000
    H               33.397000000000    23.059000000000    13.687000000000
    C               33.533000000000    25.097000000000    12.978000000000
    H               33.628000000000    25.293000000000    11.911000000000
    H               32.474750000000    25.211950000000    13.255310000000
    H               34.187770000000    25.819560000000    13.487110000000
    H               35.491400000000    22.347510000000    13.195000000000

""")

GEOS['%s-%s-%s' % (dbse, 'thr55-asp58-backbone', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               27.299000000000    17.055000000000    13.533000000000
    H               27.287000000000    17.474000000000    12.529000000000
    C               28.236000000000    15.864000000000    13.558000000000
    H               28.192000000000    15.374000000000    14.531000000000
    H               27.945000000000    15.155000000000    12.783000000000
    H               29.254000000000    16.201000000000    13.365000000000
    O               25.925000000000    16.611000000000    13.913000000000
    H               25.537000000000    16.148000000000    13.159000000000
    H               27.493080000000    17.857440000000    14.259930000000
    --
    0 1
    C               22.229000000000    18.244000000000    13.325000000000
    O               20.963000000000    18.253000000000    13.395000000000
    N               22.997000000000    17.978000000000    14.366000000000
    H               24.005000000000    17.928000000000    14.255000000000
    C               22.418000000000    17.638000000000    15.693000000000
    H               21.801000000000    16.746000000000    15.572000000000
    H               21.730060000000    18.427190000000    16.030510000000
    H               23.177750000000    17.414370000000    16.456390000000
    H               22.734420000000    18.490530000000    12.379610000000

""")

GEOS['%s-%s-%s' % (dbse, 'thr55-asp58-sidechain-small', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               28.375000000000    17.999000000000    16.887000000000
    O               29.326000000000    18.786000000000    16.690000000000
    N               27.510000000000    17.689000000000    15.954000000000
    H               26.729000000000    17.079000000000    16.191000000000
    C               27.574000000000    18.192000000000    14.563000000000
    H               28.529000000000    18.639000000000    14.366000000000
    H               26.797560000000    18.965600000000    14.469860000000
    H               27.379920000000    17.389560000000    13.836070000000
    H               28.181940000000    17.594070000000    17.891370000000
    --
    -1 1
    C               23.461000000000    17.331000000000    16.741000000000
    H               24.175000000000    18.147000000000    16.796000000000
    H               22.956000000000    17.309000000000    17.708000000000
    C               24.184000000000    16.016000000000    16.619000000000
    O               25.303000000000    15.894000000000    17.152000000000
    O               23.572000000000    15.107000000000    15.975000000000
    H               22.701250000000    17.554630000000    15.977610000000

""")

GEOS['%s-%s-%s' % (dbse, 'thr55-asp58-sidechain', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               28.375000000000    17.999000000000    16.887000000000
    O               29.326000000000    18.786000000000    16.690000000000
    N               27.510000000000    17.689000000000    15.954000000000
    H               26.729000000000    17.079000000000    16.191000000000
    C               27.574000000000    18.192000000000    14.563000000000
    H               28.529000000000    18.639000000000    14.366000000000
    C               26.482000000000    19.280000000000    14.432000000000
    O               25.609000000000    19.388000000000    15.287000000000
    N               26.585000000000    20.063000000000    13.378000000000
    H               27.373000000000    19.916000000000    12.751000000000
    H               25.844960000000    20.844110000000    13.149490000000
    H               27.379920000000    17.389560000000    13.836070000000
    H               28.181940000000    17.594070000000    17.891370000000
    --
    -1 1
    C               23.461000000000    17.331000000000    16.741000000000
    H               24.175000000000    18.147000000000    16.796000000000
    H               22.956000000000    17.309000000000    17.708000000000
    C               24.184000000000    16.016000000000    16.619000000000
    O               25.303000000000    15.894000000000    17.152000000000
    O               23.572000000000    15.107000000000    15.975000000000
    H               22.701250000000    17.554630000000    15.977610000000

""")

GEOS['%s-%s-%s' % (dbse, 'thr7-gly10', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               30.225000000000    38.643000000000    16.662000000000
    H               30.179000000000    37.771000000000    17.313000000000
    C               29.664000000000    39.839000000000    17.434000000000
    O               28.850000000000    40.565000000000    16.859000000000
    N               30.132000000000    40.069000000000    18.642000000000
    H               30.855000000000    39.482000000000    19.022000000000
    H               29.627320000000    38.483670000000    15.752390000000
    H               31.282870000000    38.807360000000    16.409200000000
    H               29.741810000000    40.894710000000    19.255150000000
    --
    0 1
    C               30.755000000000    44.351000000000    16.277000000000
    O               31.207000000000    45.268000000000    15.566000000000
    N               29.721000000000    43.673000000000    15.885000000000
    H               29.439000000000    42.883000000000    16.446000000000
    C               28.978000000000    43.960000000000    14.678000000000
    H               27.906000000000    43.924000000000    14.672000000000
    H               29.055000000000    45.044000000000    14.581000000000
    H               29.437240000000    43.627670000000    13.735310000000
    H               31.245010000000    44.049060000000    17.214400000000

""")

GEOS['%s-%s-%s' % (dbse, 'thr7-ile13', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               31.744000000000    38.879000000000    16.299000000000
    H               32.375000000000    38.789000000000    17.181000000000
    C               32.260000000000    37.969000000000    15.171000000000
    H               31.795000000000    38.224000000000    14.220000000000
    H               33.339000000000    38.089000000000    15.075000000000
    H               32.039000000000    36.929000000000    15.407000000000
    O               31.737000000000    40.257000000000    15.824000000000
    H               32.648000000000    40.543000000000    15.685000000000
    H               30.686130000000    38.714650000000    16.551800000000
    --
    0 1
    C               32.995000000000    35.883000000000     9.934000000000
    H               33.813000000000    36.384000000000     9.412000000000
    C               33.109000000000    36.381000000000    11.435000000000
    H               32.449000000000    35.800000000000    12.081000000000
    H               34.130000000000    36.321000000000    11.799000000000
    H               32.849000000000    37.436000000000    11.509000000000
    C               33.306000000000    34.381000000000     9.840000000000
    H               32.459000000000    33.813000000000    10.228000000000
    H               33.449000000000    34.098000000000     8.802000000000
    C               34.535000000000    34.028000000000    10.720000000000
    H               34.329000000000    34.129000000000    11.785000000000
    H               34.809000000000    32.988000000000    10.545000000000
    H               35.346000000000    34.662000000000    10.400000000000
    H               32.083020000000    36.173400000000     9.391819000000

""")

GEOS['%s-%s-%s' % (dbse, 'thr7-lys11', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               28.469000000000    37.475000000000    15.420000000000
    O               28.213000000000    36.753000000000    16.411000000000
    N               29.426000000000    38.430000000000    15.446000000000
    H               29.612000000000    39.012000000000    14.635000000000
    C               30.225000000000    38.643000000000    16.662000000000
    H               30.179000000000    37.771000000000    17.313000000000
    H               27.899930000000    37.357970000000    14.485940000000
    H               31.282870000000    38.807360000000    16.409200000000
    H               29.821690000000    39.502830000000    17.217010000000
    --
    0 1
    C               31.191000000000    42.012000000000    12.331000000000
    H               31.052000000000    42.639000000000    11.449000000000
    C               30.459000000000    40.666000000000    12.130000000000
    O               30.253000000000    39.991000000000    13.133000000000
    N               30.163000000000    40.338000000000    10.886000000000
    H               30.353000000000    40.972000000000    10.114000000000
    H               29.700030000000    39.355410000000    10.712290000000
    H               30.716850000000    42.473320000000    13.209840000000
    H               32.262710000000    41.798530000000    12.456910000000

""")

GEOS['%s-%s-%s' % (dbse, 'thr7-thr9-v2', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               31.744000000000    38.879000000000    16.299000000000
    H               32.375000000000    38.789000000000    17.181000000000
    C               32.260000000000    37.969000000000    15.171000000000
    H               31.795000000000    38.224000000000    14.220000000000
    H               33.339000000000    38.089000000000    15.075000000000
    H               32.039000000000    36.929000000000    15.407000000000
    O               31.737000000000    40.257000000000    15.824000000000
    H               32.648000000000    40.543000000000    15.685000000000
    H               30.686130000000    38.714650000000    16.551800000000
    --
    0 1
    C               32.979000000000    43.918000000000    17.445000000000
    H               33.355000000000    44.924000000000    17.255000000000
    C               33.657000000000    43.319000000000    18.672000000000
    H               33.468000000000    42.248000000000    18.740000000000
    H               34.733000000000    43.473000000000    18.588000000000
    H               33.299000000000    43.817000000000    19.573000000000
    O               33.174000000000    43.067000000000    16.265000000000
    H               33.740000000000    43.561000000000    15.655000000000
    H               31.881750000000    43.933500000000    17.521110000000

""")

GEOS['%s-%s-%s' % (dbse, 'thr7-thr9', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               31.744000000000    38.879000000000    16.299000000000
    H               32.375000000000    38.789000000000    17.181000000000
    C               32.260000000000    37.969000000000    15.171000000000
    H               31.795000000000    38.224000000000    14.220000000000
    H               33.339000000000    38.089000000000    15.075000000000
    H               32.039000000000    36.929000000000    15.407000000000
    O               31.737000000000    40.257000000000    15.824000000000
    H               32.648000000000    40.543000000000    15.685000000000
    H               30.686130000000    38.714650000000    16.551800000000
    --
    0 1
    C               30.075000000000    42.538000000000    18.984000000000
    O               29.586000000000    43.570000000000    19.483000000000
    N               30.991000000000    42.571000000000    17.998000000000
    H               31.378000000000    41.738000000000    17.572000000000
    C               31.422000000000    43.940000000000    17.553000000000
    H               31.183000000000    44.701000000000    18.293000000000
    H               29.735290000000    41.552260000000    19.334600000000
    H               32.519260000000    43.924500000000    17.476890000000
    H               30.931990000000    44.241940000000    16.615600000000

""")

GEOS['%s-%s-%s' % (dbse, 'tyr59-ile61', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               22.085000000000    22.254000000000    16.581000000000
    H               22.734000000000    22.409000000000    15.719000000000
    H               21.556000000000    23.189000000000    16.770000000000
    C               22.945000000000    21.951000000000    17.785000000000
    C               24.272000000000    21.544000000000    17.644000000000
    H               24.693000000000    21.388000000000    16.659000000000
    C               25.052000000000    21.285000000000    18.776000000000
    H               26.081000000000    20.990000000000    18.667000000000
    C               24.517000000000    21.470000000000    20.030000000000
    O               25.248000000000    21.302000000000    21.191000000000
    H               24.728000000000    21.394000000000    21.993000000000
    C               23.204000000000    21.907000000000    20.192000000000
    H               22.798000000000    22.050000000000    21.180000000000
    C               22.437000000000    22.157000000000    19.065000000000
    H               21.415000000000    22.483000000000    19.173000000000
    H               21.361890000000    21.459740000000    16.343800000000
    --
    0 1
    C               22.650000000000    24.516000000000    12.172000000000
    H               23.215000000000    23.598000000000    12.003000000000
    C               23.376000000000    25.645000000000    11.409000000000
    H               22.954000000000    26.611000000000    11.687000000000
    H               24.444000000000    25.631000000000    11.623000000000
    H               23.277000000000    25.511000000000    10.333000000000
    C               22.662000000000    24.819000000000    13.699000000000
    H               22.117000000000    25.741000000000    13.891000000000
    H               22.184000000000    24.007000000000    14.242000000000
    C               24.123000000000    24.981000000000    14.195000000000
    H               24.588000000000    25.879000000000    13.799000000000
    H               24.114000000000    25.061000000000    15.282000000000
    H               24.702000000000    24.100000000000    13.921000000000
    H               21.618800000000    24.338040000000    11.832960000000

""")

GEOS['%s-%s-%s' % (dbse, 'val17-leu56', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               29.903000000000    24.590000000000     7.665000000000
    H               29.900000000000    23.663000000000     8.237000000000
    C               30.862000000000    25.496000000000     8.389000000000
    H               30.895000000000    26.467000000000     7.902000000000
    H               30.547000000000    25.617000000000     9.425000000000
    H               31.853000000000    25.045000000000     8.370000000000
    C               28.476000000000    25.135000000000     7.705000000000
    H               27.762000000000    24.367000000000     7.412000000000
    H               28.235000000000    25.443000000000     8.721000000000
    H               28.371000000000    26.000000000000     7.051000000000
    H               30.174450000000    24.346750000000     6.627144000000
    --
    0 1
    C               26.084000000000    21.888000000000    11.833000000000
    H               26.170000000000    21.168000000000    11.029000000000
    H               25.313000000000    22.606000000000    11.554000000000
    C               27.426000000000    22.616000000000    11.902000000000
    H               28.230000000000    21.920000000000    12.129000000000
    C               27.718000000000    23.341000000000    10.578000000000
    H               26.935000000000    24.062000000000    10.350000000000
    H               28.683000000000    23.842000000000    10.643000000000
    H               27.773000000000    22.608000000000     9.776000000000
    C               27.380000000000    23.721000000000    12.955000000000
    H               27.051000000000    23.309000000000    13.905000000000
    H               28.369000000000    24.127000000000    13.111000000000
    H               26.698000000000    24.514000000000    12.653000000000
    H               25.734770000000    21.332790000000    12.716050000000

""")

GEOS['%s-%s-%s' % (dbse, 'val26-ile30', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               33.533000000000    25.097000000000    12.978000000000
    H               33.628000000000    25.293000000000    11.911000000000
    C               34.441000000000    26.099000000000    13.684000000000
    O               34.883000000000    27.090000000000    13.093000000000
    N               34.734000000000    25.822000000000    14.949000000000
    H               34.360000000000    24.976000000000    15.377000000000
    H               35.379190000000    26.490390000000    15.538050000000
    H               33.890540000000    24.077160000000    13.183170000000
    H               32.474750000000    25.211950000000    13.255310000000
    --
    0 1
    C               37.441000000000    28.882000000000    11.052000000000
    O               38.020000000000    29.772000000000    10.382000000000
    N               36.811000000000    29.170000000000    12.192000000000
    H               36.331000000000    28.429000000000    12.692000000000
    C               36.731000000000    30.570000000000    12.645000000000
    H               36.443000000000    31.212000000000    11.812000000000
    H               37.462430000000    27.816960000000    10.777710000000
    H               37.746360000000    30.864510000000    12.948820000000
    H               36.010140000000    30.715160000000    13.463100000000

""")

GEOS['%s-%s-%s' % (dbse, 'val26-leu43', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               32.060000000000    25.257000000000    13.364000000000
    H               31.920000000000    24.919000000000    14.391000000000
    C               31.684000000000    26.749000000000    13.342000000000
    H               31.902000000000    27.168000000000    12.361000000000
    H               30.619000000000    26.858000000000    13.547000000000
    H               32.222000000000    27.303000000000    14.110000000000
    C               31.152000000000    24.421000000000    12.477000000000
    H               31.066000000000    23.409000000000    12.871000000000
    H               30.173000000000    24.880000000000    12.423000000000
    H               31.551000000000    24.370000000000    11.464000000000
    H               33.118250000000    25.142050000000    13.086690000000
    --
    0 1
    C               29.151000000000    28.655000000000    18.755000000000
    H               29.215000000000    27.648000000000    19.172000000000
    H               28.322000000000    28.637000000000    18.045000000000
    C               30.416000000000    28.912000000000    17.980000000000
    H               31.242000000000    29.080000000000    18.671000000000
    C               30.738000000000    27.693000000000    17.122000000000
    H               29.921000000000    27.482000000000    16.432000000000
    H               31.647000000000    27.884000000000    16.550000000000
    H               30.912000000000    26.825000000000    17.756000000000
    C               30.205000000000    30.168000000000    17.129000000000
    H               30.017000000000    31.034000000000    17.761000000000
    H               31.110000000000    30.356000000000    16.551000000000
    H               29.373000000000    30.026000000000    16.440000000000
    H               28.870000000000    29.318130000000    19.586440000000

""")

GEOS['%s-%s-%s' % (dbse, 'val26-leu56', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               32.060000000000    25.257000000000    13.364000000000
    H               31.920000000000    24.919000000000    14.391000000000
    C               31.684000000000    26.749000000000    13.342000000000
    H               31.902000000000    27.168000000000    12.361000000000
    H               30.619000000000    26.858000000000    13.547000000000
    H               32.222000000000    27.303000000000    14.110000000000
    C               31.152000000000    24.421000000000    12.477000000000
    H               31.066000000000    23.409000000000    12.871000000000
    H               30.173000000000    24.880000000000    12.423000000000
    H               31.551000000000    24.370000000000    11.464000000000
    H               33.118250000000    25.142050000000    13.086690000000
    --
    0 1
    C               26.084000000000    21.888000000000    11.833000000000
    H               26.170000000000    21.168000000000    11.029000000000
    H               25.313000000000    22.606000000000    11.554000000000
    C               27.426000000000    22.616000000000    11.902000000000
    H               28.230000000000    21.920000000000    12.129000000000
    C               27.718000000000    23.341000000000    10.578000000000
    H               26.935000000000    24.062000000000    10.350000000000
    H               28.683000000000    23.842000000000    10.643000000000
    H               27.773000000000    22.608000000000     9.776000000000
    C               27.380000000000    23.721000000000    12.955000000000
    H               27.051000000000    23.309000000000    13.905000000000
    H               28.369000000000    24.127000000000    13.111000000000
    H               26.698000000000    24.514000000000    12.653000000000
    H               25.734770000000    21.332790000000    12.716050000000

""")

GEOS['%s-%s-%s' % (dbse, 'val5-ile13', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               27.151000000000    33.362000000000    10.650000000000
    O               26.350000000000    32.778000000000    11.395000000000
    N               28.260000000000    33.943000000000    11.096000000000
    H               28.901000000000    34.360000000000    10.429000000000
    C               28.605000000000    33.965000000000    12.503000000000
    H               27.847000000000    33.457000000000    13.097000000000
    C               28.638000000000    35.461000000000    12.900000000000
    O               29.522000000000    36.103000000000    12.320000000000
    N               27.751000000000    35.867000000000    13.740000000000
    H               27.107000000000    35.195000000000    14.151000000000
    H               26.873700000000    33.416140000000     9.586903000000
    H               29.577220000000    33.501080000000    12.725650000000
    H               27.707120000000    36.925880000000    14.034700000000
    --
    0 1
    C               30.494000000000    38.261000000000     9.729000000000
    O               30.849000000000    38.850000000000     8.706000000000
    N               30.795000000000    37.015000000000    10.095000000000
    H               30.443000000000    36.586000000000    10.942000000000
    C               31.720000000000    36.289000000000     9.176000000000
    H               32.071000000000    36.955000000000     8.389000000000
    C               30.955000000000    35.211000000000     8.459000000000
    O               30.025000000000    34.618000000000     9.040000000000
    N               31.244000000000    34.986000000000     7.197000000000
    H               31.951000000000    35.496000000000     6.692000000000
    H               30.699610000000    34.174210000000     6.692391000000
    H               32.631980000000    35.998600000000     9.718181000000
    H               29.808860000000    38.807240000000    10.393990000000

""")

GEOS['%s-%s-%s' % (dbse, 'val5-leu15', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               29.963000000000    33.317000000000    12.814000000000
    H               30.765000000000    33.836000000000    12.289000000000
    C               30.211000000000    33.394000000000    14.304000000000
    H               29.387000000000    32.936000000000    14.852000000000
    H               31.131000000000    32.860000000000    14.542000000000
    H               30.335000000000    34.429000000000    14.622000000000
    C               29.957000000000    31.838000000000    12.352000000000
    H               29.843000000000    31.787000000000    11.269000000000
    H               30.898000000000    31.368000000000    12.633000000000
    H               29.135000000000    31.300000000000    12.823000000000
    H               28.990780000000    33.780920000000    12.591350000000
    --
    0 1
    C               31.562000000000    29.686000000000     8.045000000000
    H               30.790000000000    30.265000000000     8.559000000000
    H               31.073000000000    28.716000000000     7.948000000000
    C               32.631000000000    29.444000000000     9.060000000000
    H               33.028000000000    28.443000000000     8.889000000000
    C               33.814000000000    30.390000000000     9.030000000000
    H               33.465000000000    31.415000000000     9.143000000000
    H               34.503000000000    30.145000000000     9.839000000000
    H               34.339000000000    30.288000000000     8.080000000000
    C               31.945000000000    29.449000000000    10.436000000000
    H               31.086000000000    28.781000000000    10.437000000000
    H               32.653000000000    29.096000000000    11.183000000000
    H               31.642000000000    30.461000000000    10.683000000000
    H               31.644750000000    30.109820000000     7.033303000000

""")

GEOS['%s-%s-%s' % (dbse, 'val5-leu69', 'dimer')] = qcdb.Molecule("""
    units Angstrom
    no_com
    no_reorient
    0 1
    C               29.963000000000    33.317000000000    12.814000000000
    H               30.765000000000    33.836000000000    12.289000000000
    C               30.211000000000    33.394000000000    14.304000000000
    H               29.387000000000    32.936000000000    14.852000000000
    H               31.131000000000    32.860000000000    14.542000000000
    H               30.335000000000    34.429000000000    14.622000000000
    C               29.957000000000    31.838000000000    12.352000000000
    H               29.843000000000    31.787000000000    11.269000000000
    H               30.898000000000    31.368000000000    12.633000000000
    H               29.135000000000    31.300000000000    12.823000000000
    H               28.990780000000    33.780920000000    12.591350000000
    --
    0 1
    C               30.925000000000    34.304000000000    17.753000000000
    H               30.790000000000    33.539000000000    16.990000000000
    H               30.827000000000    35.284000000000    17.283000000000
    C               32.345000000000    34.183000000000    18.358000000000
    H               32.492000000000    34.899000000000    19.164000000000
    C               32.555000000000    32.783000000000    18.870000000000
    H               32.403000000000    32.055000000000    18.071000000000
    H               33.579000000000    32.699000000000    19.225000000000
    H               31.891000000000    32.566000000000    19.702000000000
    C               33.361000000000    34.491000000000    17.245000000000
    H               33.262000000000    35.526000000000    16.923000000000
    H               34.372000000000    34.336000000000    17.619000000000
    H               33.201000000000    33.829000000000    16.393000000000
    H               30.134520000000    34.192180000000    18.509730000000

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

DATA['SAPT ELST ENERGY'] = {}
DATA['SAPT ELST ENERGY']['UBQ-ala28-asp32'] =   -11.7248
DATA['SAPT ELST ENERGY']['UBQ-ala46-lys48'] =   -24.2119
DATA['SAPT ELST ENERGY']['UBQ-arg42-val70'] =   -21.2677
DATA['SAPT ELST ENERGY']['UBQ-asn25-lys29'] =   -10.2528
DATA['SAPT ELST ENERGY']['UBQ-asp21-leu56'] =    -9.2462
DATA['SAPT ELST ENERGY']['UBQ-gln2-glu64-small'] =   -11.9284
DATA['SAPT ELST ENERGY']['UBQ-gln31-gly35'] =   -12.6573
DATA['SAPT ELST ENERGY']['UBQ-gln40-arg72'] =   -15.1020
DATA['SAPT ELST ENERGY']['UBQ-gln41-arg72'] =   -11.3788
DATA['SAPT ELST ENERGY']['UBQ-gln62-ser65'] =   -28.0651
DATA['SAPT ELST ENERGY']['UBQ-glu16-lys29'] =   -23.6167
DATA['SAPT ELST ENERGY']['UBQ-glu18-asp21-C'] =   -11.3551
DATA['SAPT ELST ENERGY']['UBQ-glu18-asp21-N'] =   -13.3323
DATA['SAPT ELST ENERGY']['UBQ-glu24-ala28'] =    -7.6454
DATA['SAPT ELST ENERGY']['UBQ-glu24-asp52'] =   -11.4119
DATA['SAPT ELST ENERGY']['UBQ-glu51-arg54'] =   -49.0845
DATA['SAPT ELST ENERGY']['UBQ-glu51-tyr59'] =    -7.1708
DATA['SAPT ELST ENERGY']['UBQ-ile13-leu15'] =    -0.3948
DATA['SAPT ELST ENERGY']['UBQ-ile13-lys33'] =    -1.3968
DATA['SAPT ELST ENERGY']['UBQ-ile23-arg54'] =   -10.9911
DATA['SAPT ELST ENERGY']['UBQ-ile23-leu50'] =    -0.7202
DATA['SAPT ELST ENERGY']['UBQ-ile23-leu56'] =    -0.6154
DATA['SAPT ELST ENERGY']['UBQ-ile23-lys27'] =    -9.1942
DATA['SAPT ELST ENERGY']['UBQ-ile3-leu15-big'] =   -13.2646
DATA['SAPT ELST ENERGY']['UBQ-ile3-leu15-part1'] =    -8.7290
DATA['SAPT ELST ENERGY']['UBQ-ile3-leu15-part2'] =    -8.5937
DATA['SAPT ELST ENERGY']['UBQ-ile3-leu15'] =    -0.4229
DATA['SAPT ELST ENERGY']['UBQ-ile3-leu67'] =    -0.4147
DATA['SAPT ELST ENERGY']['UBQ-ile3-val17'] =    -0.6497
DATA['SAPT ELST ENERGY']['UBQ-ile30-glu34'] =   -10.6798
DATA['SAPT ELST ENERGY']['UBQ-ile30-leu43'] =    -0.3023
DATA['SAPT ELST ENERGY']['UBQ-ile36-gln41'] =    -8.4861
DATA['SAPT ELST ENERGY']['UBQ-ile36-leu69'] =    -0.2999
DATA['SAPT ELST ENERGY']['UBQ-ile36-leu71'] =    -0.8897
DATA['SAPT ELST ENERGY']['UBQ-ile44-gly47'] =    -0.3886
DATA['SAPT ELST ENERGY']['UBQ-ile44-his68-big'] =    -0.5677
DATA['SAPT ELST ENERGY']['UBQ-ile44-his68-small'] =    -0.5516
DATA['SAPT ELST ENERGY']['UBQ-ile44-his68'] =   -25.5115
DATA['SAPT ELST ENERGY']['UBQ-ile61-leu67'] =    -0.5310
DATA['SAPT ELST ENERGY']['UBQ-leu15-ile30'] =    -0.4560
DATA['SAPT ELST ENERGY']['UBQ-leu15-val26'] =    -0.3602
DATA['SAPT ELST ENERGY']['UBQ-leu43-leu50'] =   -12.1545
DATA['SAPT ELST ENERGY']['UBQ-leu43-leu69'] =    -0.4344
DATA['SAPT ELST ENERGY']['UBQ-leu50-ile61'] =    -0.2081
DATA['SAPT ELST ENERGY']['UBQ-leu50-tyr59'] =    -1.2697
DATA['SAPT ELST ENERGY']['UBQ-leu56-ile61np'] =    -0.9028
DATA['SAPT ELST ENERGY']['UBQ-leu56-ile61p'] =    -5.1228
DATA['SAPT ELST ENERGY']['UBQ-leu56-tyr59'] =    -4.7242
DATA['SAPT ELST ENERGY']['UBQ-leu69-leu71'] =    -0.5102
DATA['SAPT ELST ENERGY']['UBQ-leu71-leu73'] =    -0.3149
DATA['SAPT ELST ENERGY']['UBQ-lys11-glu34'] =   -95.6878
DATA['SAPT ELST ENERGY']['UBQ-lys27-asp52'] =  -105.3789
DATA['SAPT ELST ENERGY']['UBQ-lys27-gln31'] =    -9.2986
DATA['SAPT ELST ENERGY']['UBQ-lys27-leu43-nonh3'] =    -0.3533
DATA['SAPT ELST ENERGY']['UBQ-lys27-leu43'] =     0.1407
DATA['SAPT ELST ENERGY']['UBQ-lys29-lys33'] =    -9.7287
DATA['SAPT ELST ENERGY']['UBQ-lys6-thr66'] =    -0.2149
DATA['SAPT ELST ENERGY']['UBQ-met1-ile3'] =    -0.4614
DATA['SAPT ELST ENERGY']['UBQ-met1-lys63'] =    -3.2075
DATA['SAPT ELST ENERGY']['UBQ-met1-val17-bi'] =   -31.5827
DATA['SAPT ELST ENERGY']['UBQ-met1-val17-mono'] =   -11.7346
DATA['SAPT ELST ENERGY']['UBQ-met1-val17'] =    -0.3404
DATA['SAPT ELST ENERGY']['UBQ-phe4-leu67'] =   -12.3241
DATA['SAPT ELST ENERGY']['UBQ-phe4-ser65'] =   -11.5522
DATA['SAPT ELST ENERGY']['UBQ-phe4-thr12'] =    -0.5898
DATA['SAPT ELST ENERGY']['UBQ-phe4-thr14'] =    -1.0980
DATA['SAPT ELST ENERGY']['UBQ-phe45-ala46'] =    -0.3407
DATA['SAPT ELST ENERGY']['UBQ-phe45-ile61'] =    -1.3694
DATA['SAPT ELST ENERGY']['UBQ-phe45-leu67'] =    -0.5102
DATA['SAPT ELST ENERGY']['UBQ-phe45-lys48'] =   -16.8229
DATA['SAPT ELST ENERGY']['UBQ-pro19-ser57'] =   -11.0448
DATA['SAPT ELST ENERGY']['UBQ-pro37-gln40'] =    -8.6466
DATA['SAPT ELST ENERGY']['UBQ-ser57-asn60'] =   -10.7975
DATA['SAPT ELST ENERGY']['UBQ-thr14-lys33'] =   -13.5642
DATA['SAPT ELST ENERGY']['UBQ-thr22-asn25'] =    -5.8638
DATA['SAPT ELST ENERGY']['UBQ-thr22-thr55-big'] =   -28.3087
DATA['SAPT ELST ENERGY']['UBQ-thr22-thr55-small'] =     0.0132
DATA['SAPT ELST ENERGY']['UBQ-thr22-val26'] =    -7.7691
DATA['SAPT ELST ENERGY']['UBQ-thr55-asp58-backbone'] =    -5.1253
DATA['SAPT ELST ENERGY']['UBQ-thr55-asp58-sidechain-small'] =   -19.2450
DATA['SAPT ELST ENERGY']['UBQ-thr55-asp58-sidechain'] =   -14.8812
DATA['SAPT ELST ENERGY']['UBQ-thr7-gly10'] =    -5.0323
DATA['SAPT ELST ENERGY']['UBQ-thr7-ile13'] =    -0.0568
DATA['SAPT ELST ENERGY']['UBQ-thr7-lys11'] =   -10.2814
DATA['SAPT ELST ENERGY']['UBQ-thr7-thr9-v2'] =    -0.8037
DATA['SAPT ELST ENERGY']['UBQ-thr7-thr9'] =    -4.1104
DATA['SAPT ELST ENERGY']['UBQ-tyr59-ile61'] =    -0.4914
DATA['SAPT ELST ENERGY']['UBQ-val17-leu56'] =    -0.7235
DATA['SAPT ELST ENERGY']['UBQ-val26-ile30'] =    -9.0470
DATA['SAPT ELST ENERGY']['UBQ-val26-leu43'] =    -0.1397
DATA['SAPT ELST ENERGY']['UBQ-val26-leu56'] =    -0.8839
DATA['SAPT ELST ENERGY']['UBQ-val5-ile13'] =   -26.3890
DATA['SAPT ELST ENERGY']['UBQ-val5-leu15'] =    -0.8344
DATA['SAPT ELST ENERGY']['UBQ-val5-leu69'] =    -0.6049
DATA['SAPT EXCH ENERGY'] = {}
DATA['SAPT EXCH ENERGY']['UBQ-ala28-asp32'] =     9.0182
DATA['SAPT EXCH ENERGY']['UBQ-ala46-lys48'] =     1.1511
DATA['SAPT EXCH ENERGY']['UBQ-arg42-val70'] =    19.0227
DATA['SAPT EXCH ENERGY']['UBQ-asn25-lys29'] =     7.1097
DATA['SAPT EXCH ENERGY']['UBQ-asp21-leu56'] =     4.8948
DATA['SAPT EXCH ENERGY']['UBQ-gln2-glu64-small'] =    10.4545
DATA['SAPT EXCH ENERGY']['UBQ-gln31-gly35'] =    12.9391
DATA['SAPT EXCH ENERGY']['UBQ-gln40-arg72'] =    16.1738
DATA['SAPT EXCH ENERGY']['UBQ-gln41-arg72'] =     3.6047
DATA['SAPT EXCH ENERGY']['UBQ-gln62-ser65'] =    33.3887
DATA['SAPT EXCH ENERGY']['UBQ-glu16-lys29'] =    13.3615
DATA['SAPT EXCH ENERGY']['UBQ-glu18-asp21-C'] =    14.6063
DATA['SAPT EXCH ENERGY']['UBQ-glu18-asp21-N'] =     8.5660
DATA['SAPT EXCH ENERGY']['UBQ-glu24-ala28'] =     3.8277
DATA['SAPT EXCH ENERGY']['UBQ-glu24-asp52'] =     8.6559
DATA['SAPT EXCH ENERGY']['UBQ-glu51-arg54'] =     0.0103
DATA['SAPT EXCH ENERGY']['UBQ-glu51-tyr59'] =     6.9113
DATA['SAPT EXCH ENERGY']['UBQ-ile13-leu15'] =     1.1927
DATA['SAPT EXCH ENERGY']['UBQ-ile13-lys33'] =     5.8960
DATA['SAPT EXCH ENERGY']['UBQ-ile23-arg54'] =     7.9228
DATA['SAPT EXCH ENERGY']['UBQ-ile23-leu50'] =     2.4150
DATA['SAPT EXCH ENERGY']['UBQ-ile23-leu56'] =     2.2041
DATA['SAPT EXCH ENERGY']['UBQ-ile23-lys27'] =     5.4804
DATA['SAPT EXCH ENERGY']['UBQ-ile3-leu15-big'] =     8.3880
DATA['SAPT EXCH ENERGY']['UBQ-ile3-leu15-part1'] =     4.5347
DATA['SAPT EXCH ENERGY']['UBQ-ile3-leu15-part2'] =     4.2336
DATA['SAPT EXCH ENERGY']['UBQ-ile3-leu15'] =     1.8195
DATA['SAPT EXCH ENERGY']['UBQ-ile3-leu67'] =     1.3970
DATA['SAPT EXCH ENERGY']['UBQ-ile3-val17'] =     2.2465
DATA['SAPT EXCH ENERGY']['UBQ-ile30-glu34'] =     7.4299
DATA['SAPT EXCH ENERGY']['UBQ-ile30-leu43'] =     0.9351
DATA['SAPT EXCH ENERGY']['UBQ-ile36-gln41'] =     5.7076
DATA['SAPT EXCH ENERGY']['UBQ-ile36-leu69'] =     1.2421
DATA['SAPT EXCH ENERGY']['UBQ-ile36-leu71'] =     3.1877
DATA['SAPT EXCH ENERGY']['UBQ-ile44-gly47'] =     0.5502
DATA['SAPT EXCH ENERGY']['UBQ-ile44-his68-big'] =     2.5237
DATA['SAPT EXCH ENERGY']['UBQ-ile44-his68-small'] =     2.5058
DATA['SAPT EXCH ENERGY']['UBQ-ile44-his68'] =    23.2787
DATA['SAPT EXCH ENERGY']['UBQ-ile61-leu67'] =     1.7875
DATA['SAPT EXCH ENERGY']['UBQ-leu15-ile30'] =     1.7255
DATA['SAPT EXCH ENERGY']['UBQ-leu15-val26'] =     1.2921
DATA['SAPT EXCH ENERGY']['UBQ-leu43-leu50'] =    10.1907
DATA['SAPT EXCH ENERGY']['UBQ-leu43-leu69'] =     1.3689
DATA['SAPT EXCH ENERGY']['UBQ-leu50-ile61'] =     0.7652
DATA['SAPT EXCH ENERGY']['UBQ-leu50-tyr59'] =     3.6014
DATA['SAPT EXCH ENERGY']['UBQ-leu56-ile61np'] =     3.3407
DATA['SAPT EXCH ENERGY']['UBQ-leu56-ile61p'] =     1.2178
DATA['SAPT EXCH ENERGY']['UBQ-leu56-tyr59'] =     3.1583
DATA['SAPT EXCH ENERGY']['UBQ-leu69-leu71'] =     1.7974
DATA['SAPT EXCH ENERGY']['UBQ-leu71-leu73'] =     1.1239
DATA['SAPT EXCH ENERGY']['UBQ-lys11-glu34'] =     4.4314
DATA['SAPT EXCH ENERGY']['UBQ-lys27-asp52'] =    10.1333
DATA['SAPT EXCH ENERGY']['UBQ-lys27-gln31'] =     6.7191
DATA['SAPT EXCH ENERGY']['UBQ-lys27-leu43-nonh3'] =     1.3172
DATA['SAPT EXCH ENERGY']['UBQ-lys27-leu43'] =     1.1785
DATA['SAPT EXCH ENERGY']['UBQ-lys29-lys33'] =     7.2810
DATA['SAPT EXCH ENERGY']['UBQ-lys6-thr66'] =     0.8261
DATA['SAPT EXCH ENERGY']['UBQ-met1-ile3'] =     1.4611
DATA['SAPT EXCH ENERGY']['UBQ-met1-lys63'] =     6.3857
DATA['SAPT EXCH ENERGY']['UBQ-met1-val17-bi'] =    25.6990
DATA['SAPT EXCH ENERGY']['UBQ-met1-val17-mono'] =     9.0562
DATA['SAPT EXCH ENERGY']['UBQ-met1-val17'] =     1.2724
DATA['SAPT EXCH ENERGY']['UBQ-phe4-leu67'] =    10.2946
DATA['SAPT EXCH ENERGY']['UBQ-phe4-ser65'] =     9.0278
DATA['SAPT EXCH ENERGY']['UBQ-phe4-thr12'] =     2.3162
DATA['SAPT EXCH ENERGY']['UBQ-phe4-thr14'] =     3.7487
DATA['SAPT EXCH ENERGY']['UBQ-phe45-ala46'] =     1.1145
DATA['SAPT EXCH ENERGY']['UBQ-phe45-ile61'] =     3.4056
DATA['SAPT EXCH ENERGY']['UBQ-phe45-leu67'] =     2.3317
DATA['SAPT EXCH ENERGY']['UBQ-phe45-lys48'] =    14.5638
DATA['SAPT EXCH ENERGY']['UBQ-pro19-ser57'] =    10.2061
DATA['SAPT EXCH ENERGY']['UBQ-pro37-gln40'] =     6.4593
DATA['SAPT EXCH ENERGY']['UBQ-ser57-asn60'] =    10.5916
DATA['SAPT EXCH ENERGY']['UBQ-thr14-lys33'] =     0.8302
DATA['SAPT EXCH ENERGY']['UBQ-thr22-asn25'] =     3.0100
DATA['SAPT EXCH ENERGY']['UBQ-thr22-thr55-big'] =    23.3514
DATA['SAPT EXCH ENERGY']['UBQ-thr22-thr55-small'] =     3.2781
DATA['SAPT EXCH ENERGY']['UBQ-thr22-val26'] =     3.6528
DATA['SAPT EXCH ENERGY']['UBQ-thr55-asp58-backbone'] =     2.3760
DATA['SAPT EXCH ENERGY']['UBQ-thr55-asp58-sidechain-small'] =     6.3731
DATA['SAPT EXCH ENERGY']['UBQ-thr55-asp58-sidechain'] =     8.6530
DATA['SAPT EXCH ENERGY']['UBQ-thr7-gly10'] =     1.5770
DATA['SAPT EXCH ENERGY']['UBQ-thr7-ile13'] =     0.2074
DATA['SAPT EXCH ENERGY']['UBQ-thr7-lys11'] =     7.3040
DATA['SAPT EXCH ENERGY']['UBQ-thr7-thr9-v2'] =     1.0384
DATA['SAPT EXCH ENERGY']['UBQ-thr7-thr9'] =     2.6207
DATA['SAPT EXCH ENERGY']['UBQ-tyr59-ile61'] =     1.7098
DATA['SAPT EXCH ENERGY']['UBQ-val17-leu56'] =     2.4373
DATA['SAPT EXCH ENERGY']['UBQ-val26-ile30'] =     5.7454
DATA['SAPT EXCH ENERGY']['UBQ-val26-leu43'] =     0.6367
DATA['SAPT EXCH ENERGY']['UBQ-val26-leu56'] =     2.9796
DATA['SAPT EXCH ENERGY']['UBQ-val5-ile13'] =    25.3482
DATA['SAPT EXCH ENERGY']['UBQ-val5-leu15'] =     2.8750
DATA['SAPT EXCH ENERGY']['UBQ-val5-leu69'] =     2.0609
DATA['SAPT IND ENERGY'] = {}
DATA['SAPT IND ENERGY']['UBQ-ala28-asp32'] =    -3.6612
DATA['SAPT IND ENERGY']['UBQ-ala46-lys48'] =    -3.9156
DATA['SAPT IND ENERGY']['UBQ-arg42-val70'] =    -6.9990
DATA['SAPT IND ENERGY']['UBQ-asn25-lys29'] =    -2.9061
DATA['SAPT IND ENERGY']['UBQ-asp21-leu56'] =    -2.3513
DATA['SAPT IND ENERGY']['UBQ-gln2-glu64-small'] =    -3.7270
DATA['SAPT IND ENERGY']['UBQ-gln31-gly35'] =    -3.5317
DATA['SAPT IND ENERGY']['UBQ-gln40-arg72'] =    -5.3124
DATA['SAPT IND ENERGY']['UBQ-gln41-arg72'] =    -2.8596
DATA['SAPT IND ENERGY']['UBQ-gln62-ser65'] =   -11.3106
DATA['SAPT IND ENERGY']['UBQ-glu16-lys29'] =   -10.2970
DATA['SAPT IND ENERGY']['UBQ-glu18-asp21-C'] =    -4.5482
DATA['SAPT IND ENERGY']['UBQ-glu18-asp21-N'] =    -7.2726
DATA['SAPT IND ENERGY']['UBQ-glu24-ala28'] =    -1.8014
DATA['SAPT IND ENERGY']['UBQ-glu24-asp52'] =    -2.9730
DATA['SAPT IND ENERGY']['UBQ-glu51-arg54'] =    -0.9686
DATA['SAPT IND ENERGY']['UBQ-glu51-tyr59'] =    -1.8882
DATA['SAPT IND ENERGY']['UBQ-ile13-leu15'] =    -0.1009
DATA['SAPT IND ENERGY']['UBQ-ile13-lys33'] =    -1.8558
DATA['SAPT IND ENERGY']['UBQ-ile23-arg54'] =    -3.0772
DATA['SAPT IND ENERGY']['UBQ-ile23-leu50'] =    -0.1968
DATA['SAPT IND ENERGY']['UBQ-ile23-leu56'] =    -0.1614
DATA['SAPT IND ENERGY']['UBQ-ile23-lys27'] =    -2.4411
DATA['SAPT IND ENERGY']['UBQ-ile3-leu15-big'] =    -3.1126
DATA['SAPT IND ENERGY']['UBQ-ile3-leu15-part1'] =    -2.0878
DATA['SAPT IND ENERGY']['UBQ-ile3-leu15-part2'] =    -1.9874
DATA['SAPT IND ENERGY']['UBQ-ile3-leu15'] =    -0.1763
DATA['SAPT IND ENERGY']['UBQ-ile3-leu67'] =    -0.1290
DATA['SAPT IND ENERGY']['UBQ-ile3-val17'] =    -0.1806
DATA['SAPT IND ENERGY']['UBQ-ile30-glu34'] =    -3.0607
DATA['SAPT IND ENERGY']['UBQ-ile30-leu43'] =    -0.0800
DATA['SAPT IND ENERGY']['UBQ-ile36-gln41'] =    -2.2817
DATA['SAPT IND ENERGY']['UBQ-ile36-leu69'] =    -0.0971
DATA['SAPT IND ENERGY']['UBQ-ile36-leu71'] =    -0.3109
DATA['SAPT IND ENERGY']['UBQ-ile44-gly47'] =    -0.2171
DATA['SAPT IND ENERGY']['UBQ-ile44-his68-big'] =    -0.3595
DATA['SAPT IND ENERGY']['UBQ-ile44-his68-small'] =    -0.3489
DATA['SAPT IND ENERGY']['UBQ-ile44-his68'] =    -8.6337
DATA['SAPT IND ENERGY']['UBQ-ile61-leu67'] =    -0.1601
DATA['SAPT IND ENERGY']['UBQ-leu15-ile30'] =    -0.1583
DATA['SAPT IND ENERGY']['UBQ-leu15-val26'] =    -0.1080
DATA['SAPT IND ENERGY']['UBQ-leu43-leu50'] =    -3.8834
DATA['SAPT IND ENERGY']['UBQ-leu43-leu69'] =    -0.1048
DATA['SAPT IND ENERGY']['UBQ-leu50-ile61'] =    -0.0565
DATA['SAPT IND ENERGY']['UBQ-leu50-tyr59'] =    -0.3818
DATA['SAPT IND ENERGY']['UBQ-leu56-ile61np'] =    -0.2737
DATA['SAPT IND ENERGY']['UBQ-leu56-ile61p'] =    -0.8393
DATA['SAPT IND ENERGY']['UBQ-leu56-tyr59'] =    -1.1212
DATA['SAPT IND ENERGY']['UBQ-leu69-leu71'] =    -0.1581
DATA['SAPT IND ENERGY']['UBQ-leu71-leu73'] =    -0.1172
DATA['SAPT IND ENERGY']['UBQ-lys11-glu34'] =    -7.6774
DATA['SAPT IND ENERGY']['UBQ-lys27-asp52'] =   -12.4700
DATA['SAPT IND ENERGY']['UBQ-lys27-gln31'] =    -2.6350
DATA['SAPT IND ENERGY']['UBQ-lys27-leu43-nonh3'] =    -0.0974
DATA['SAPT IND ENERGY']['UBQ-lys27-leu43'] =    -1.1590
DATA['SAPT IND ENERGY']['UBQ-lys29-lys33'] =    -2.8604
DATA['SAPT IND ENERGY']['UBQ-lys6-thr66'] =    -0.0662
DATA['SAPT IND ENERGY']['UBQ-met1-ile3'] =    -0.1633
DATA['SAPT IND ENERGY']['UBQ-met1-lys63'] =    -0.8274
DATA['SAPT IND ENERGY']['UBQ-met1-val17-bi'] =   -11.7954
DATA['SAPT IND ENERGY']['UBQ-met1-val17-mono'] =    -3.4538
DATA['SAPT IND ENERGY']['UBQ-met1-val17'] =    -0.1257
DATA['SAPT IND ENERGY']['UBQ-phe4-leu67'] =    -3.8516
DATA['SAPT IND ENERGY']['UBQ-phe4-ser65'] =    -3.6120
DATA['SAPT IND ENERGY']['UBQ-phe4-thr12'] =    -0.2676
DATA['SAPT IND ENERGY']['UBQ-phe4-thr14'] =    -0.3769
DATA['SAPT IND ENERGY']['UBQ-phe45-ala46'] =    -0.1133
DATA['SAPT IND ENERGY']['UBQ-phe45-ile61'] =    -0.3818
DATA['SAPT IND ENERGY']['UBQ-phe45-leu67'] =    -0.2510
DATA['SAPT IND ENERGY']['UBQ-phe45-lys48'] =    -5.2366
DATA['SAPT IND ENERGY']['UBQ-pro19-ser57'] =    -3.3413
DATA['SAPT IND ENERGY']['UBQ-pro37-gln40'] =    -2.4968
DATA['SAPT IND ENERGY']['UBQ-ser57-asn60'] =    -3.6606
DATA['SAPT IND ENERGY']['UBQ-thr14-lys33'] =    -3.8349
DATA['SAPT IND ENERGY']['UBQ-thr22-asn25'] =    -1.0897
DATA['SAPT IND ENERGY']['UBQ-thr22-thr55-big'] =    -9.1113
DATA['SAPT IND ENERGY']['UBQ-thr22-thr55-small'] =    -0.5978
DATA['SAPT IND ENERGY']['UBQ-thr22-val26'] =    -1.7651
DATA['SAPT IND ENERGY']['UBQ-thr55-asp58-backbone'] =    -0.9030
DATA['SAPT IND ENERGY']['UBQ-thr55-asp58-sidechain-small'] =    -5.1921
DATA['SAPT IND ENERGY']['UBQ-thr55-asp58-sidechain'] =    -6.2431
DATA['SAPT IND ENERGY']['UBQ-thr7-gly10'] =    -0.9265
DATA['SAPT IND ENERGY']['UBQ-thr7-ile13'] =    -0.0168
DATA['SAPT IND ENERGY']['UBQ-thr7-lys11'] =    -2.9891
DATA['SAPT IND ENERGY']['UBQ-thr7-thr9-v2'] =    -0.3505
DATA['SAPT IND ENERGY']['UBQ-thr7-thr9'] =    -0.9363
DATA['SAPT IND ENERGY']['UBQ-tyr59-ile61'] =    -0.1559
DATA['SAPT IND ENERGY']['UBQ-val17-leu56'] =    -0.1770
DATA['SAPT IND ENERGY']['UBQ-val26-ile30'] =    -2.3882
DATA['SAPT IND ENERGY']['UBQ-val26-leu43'] =    -0.0502
DATA['SAPT IND ENERGY']['UBQ-val26-leu56'] =    -0.2945
DATA['SAPT IND ENERGY']['UBQ-val5-ile13'] =    -9.3600
DATA['SAPT IND ENERGY']['UBQ-val5-leu15'] =    -0.2229
DATA['SAPT IND ENERGY']['UBQ-val5-leu69'] =    -0.1576
DATA['SAPT DISP ENERGY'] = {}
DATA['SAPT DISP ENERGY']['UBQ-ala28-asp32'] =    -2.6231
DATA['SAPT DISP ENERGY']['UBQ-ala46-lys48'] =    -1.3880
DATA['SAPT DISP ENERGY']['UBQ-arg42-val70'] =    -5.9718
DATA['SAPT DISP ENERGY']['UBQ-asn25-lys29'] =    -2.2387
DATA['SAPT DISP ENERGY']['UBQ-asp21-leu56'] =    -2.2932
DATA['SAPT DISP ENERGY']['UBQ-gln2-glu64-small'] =    -2.9918
DATA['SAPT DISP ENERGY']['UBQ-gln31-gly35'] =    -3.6543
DATA['SAPT DISP ENERGY']['UBQ-gln40-arg72'] =    -3.5238
DATA['SAPT DISP ENERGY']['UBQ-gln41-arg72'] =    -1.6371
DATA['SAPT DISP ENERGY']['UBQ-gln62-ser65'] =    -7.0350
DATA['SAPT DISP ENERGY']['UBQ-glu16-lys29'] =    -3.3623
DATA['SAPT DISP ENERGY']['UBQ-glu18-asp21-C'] =    -4.9466
DATA['SAPT DISP ENERGY']['UBQ-glu18-asp21-N'] =    -3.9782
DATA['SAPT DISP ENERGY']['UBQ-glu24-ala28'] =    -1.6326
DATA['SAPT DISP ENERGY']['UBQ-glu24-asp52'] =    -2.7900
DATA['SAPT DISP ENERGY']['UBQ-glu51-arg54'] =    -0.1657
DATA['SAPT DISP ENERGY']['UBQ-glu51-tyr59'] =    -2.4052
DATA['SAPT DISP ENERGY']['UBQ-ile13-leu15'] =    -1.3303
DATA['SAPT DISP ENERGY']['UBQ-ile13-lys33'] =    -3.4067
DATA['SAPT DISP ENERGY']['UBQ-ile23-arg54'] =    -2.4429
DATA['SAPT DISP ENERGY']['UBQ-ile23-leu50'] =    -2.3971
DATA['SAPT DISP ENERGY']['UBQ-ile23-leu56'] =    -1.3890
DATA['SAPT DISP ENERGY']['UBQ-ile23-lys27'] =    -1.9782
DATA['SAPT DISP ENERGY']['UBQ-ile3-leu15-big'] =    -3.9753
DATA['SAPT DISP ENERGY']['UBQ-ile3-leu15-part1'] =    -1.8234
DATA['SAPT DISP ENERGY']['UBQ-ile3-leu15-part2'] =    -1.7616
DATA['SAPT DISP ENERGY']['UBQ-ile3-leu15'] =    -1.7068
DATA['SAPT DISP ENERGY']['UBQ-ile3-leu67'] =    -1.5652
DATA['SAPT DISP ENERGY']['UBQ-ile3-val17'] =    -2.2371
DATA['SAPT DISP ENERGY']['UBQ-ile30-glu34'] =    -2.3423
DATA['SAPT DISP ENERGY']['UBQ-ile30-leu43'] =    -1.0335
DATA['SAPT DISP ENERGY']['UBQ-ile36-gln41'] =    -1.9844
DATA['SAPT DISP ENERGY']['UBQ-ile36-leu69'] =    -1.1791
DATA['SAPT DISP ENERGY']['UBQ-ile36-leu71'] =    -2.6527
DATA['SAPT DISP ENERGY']['UBQ-ile44-gly47'] =    -1.1752
DATA['SAPT DISP ENERGY']['UBQ-ile44-his68-big'] =    -1.4508
DATA['SAPT DISP ENERGY']['UBQ-ile44-his68-small'] =    -1.1032
DATA['SAPT DISP ENERGY']['UBQ-ile44-his68'] =    -6.2437
DATA['SAPT DISP ENERGY']['UBQ-ile61-leu67'] =    -1.6368
DATA['SAPT DISP ENERGY']['UBQ-leu15-ile30'] =    -1.5679
DATA['SAPT DISP ENERGY']['UBQ-leu15-val26'] =    -0.9500
DATA['SAPT DISP ENERGY']['UBQ-leu43-leu50'] =    -2.7452
DATA['SAPT DISP ENERGY']['UBQ-leu43-leu69'] =    -1.4963
DATA['SAPT DISP ENERGY']['UBQ-leu50-ile61'] =    -0.8700
DATA['SAPT DISP ENERGY']['UBQ-leu50-tyr59'] =    -3.9275
DATA['SAPT DISP ENERGY']['UBQ-leu56-ile61np'] =    -2.7135
DATA['SAPT DISP ENERGY']['UBQ-leu56-ile61p'] =    -0.9829
DATA['SAPT DISP ENERGY']['UBQ-leu56-tyr59'] =    -1.6519
DATA['SAPT DISP ENERGY']['UBQ-leu69-leu71'] =    -1.8072
DATA['SAPT DISP ENERGY']['UBQ-leu71-leu73'] =    -1.2259
DATA['SAPT DISP ENERGY']['UBQ-lys11-glu34'] =    -2.2248
DATA['SAPT DISP ENERGY']['UBQ-lys27-asp52'] =    -2.9174
DATA['SAPT DISP ENERGY']['UBQ-lys27-gln31'] =    -2.1507
DATA['SAPT DISP ENERGY']['UBQ-lys27-leu43-nonh3'] =    -1.7636
DATA['SAPT DISP ENERGY']['UBQ-lys27-leu43'] =    -1.7160
DATA['SAPT DISP ENERGY']['UBQ-lys29-lys33'] =    -2.2883
DATA['SAPT DISP ENERGY']['UBQ-lys6-thr66'] =    -0.9303
DATA['SAPT DISP ENERGY']['UBQ-met1-ile3'] =    -1.3473
DATA['SAPT DISP ENERGY']['UBQ-met1-lys63'] =    -4.4536
DATA['SAPT DISP ENERGY']['UBQ-met1-val17-bi'] =    -6.5429
DATA['SAPT DISP ENERGY']['UBQ-met1-val17-mono'] =    -2.2646
DATA['SAPT DISP ENERGY']['UBQ-met1-val17'] =    -1.1762
DATA['SAPT DISP ENERGY']['UBQ-phe4-leu67'] =    -2.7411
DATA['SAPT DISP ENERGY']['UBQ-phe4-ser65'] =    -2.5172
DATA['SAPT DISP ENERGY']['UBQ-phe4-thr12'] =    -1.7032
DATA['SAPT DISP ENERGY']['UBQ-phe4-thr14'] =    -2.2947
DATA['SAPT DISP ENERGY']['UBQ-phe45-ala46'] =    -1.2309
DATA['SAPT DISP ENERGY']['UBQ-phe45-ile61'] =    -3.6168
DATA['SAPT DISP ENERGY']['UBQ-phe45-leu67'] =    -2.1299
DATA['SAPT DISP ENERGY']['UBQ-phe45-lys48'] =    -5.1716
DATA['SAPT DISP ENERGY']['UBQ-pro19-ser57'] =    -3.6994
DATA['SAPT DISP ENERGY']['UBQ-pro37-gln40'] =    -2.2360
DATA['SAPT DISP ENERGY']['UBQ-ser57-asn60'] =    -3.1094
DATA['SAPT DISP ENERGY']['UBQ-thr14-lys33'] =    -1.3346
DATA['SAPT DISP ENERGY']['UBQ-thr22-asn25'] =    -1.2953
DATA['SAPT DISP ENERGY']['UBQ-thr22-thr55-big'] =   -10.3148
DATA['SAPT DISP ENERGY']['UBQ-thr22-thr55-small'] =    -2.8377
DATA['SAPT DISP ENERGY']['UBQ-thr22-val26'] =    -1.5610
DATA['SAPT DISP ENERGY']['UBQ-thr55-asp58-backbone'] =    -1.1305
DATA['SAPT DISP ENERGY']['UBQ-thr55-asp58-sidechain-small'] =    -2.2623
DATA['SAPT DISP ENERGY']['UBQ-thr55-asp58-sidechain'] =    -3.4331
DATA['SAPT DISP ENERGY']['UBQ-thr7-gly10'] =    -1.1545
DATA['SAPT DISP ENERGY']['UBQ-thr7-ile13'] =    -0.4556
DATA['SAPT DISP ENERGY']['UBQ-thr7-lys11'] =    -2.2898
DATA['SAPT DISP ENERGY']['UBQ-thr7-thr9-v2'] =    -0.7523
DATA['SAPT DISP ENERGY']['UBQ-thr7-thr9'] =    -1.1990
DATA['SAPT DISP ENERGY']['UBQ-tyr59-ile61'] =    -1.5476
DATA['SAPT DISP ENERGY']['UBQ-val17-leu56'] =    -1.9615
DATA['SAPT DISP ENERGY']['UBQ-val26-ile30'] =    -1.9879
DATA['SAPT DISP ENERGY']['UBQ-val26-leu43'] =    -0.8272
DATA['SAPT DISP ENERGY']['UBQ-val26-leu56'] =    -1.8626
DATA['SAPT DISP ENERGY']['UBQ-val5-ile13'] =    -6.6829
DATA['SAPT DISP ENERGY']['UBQ-val5-leu15'] =    -1.5059
DATA['SAPT DISP ENERGY']['UBQ-val5-leu69'] =    -1.5274
