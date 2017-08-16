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
| Database (Sherrill) of interaction energies for dissociation curves of dispersion-bound bimolecular complexes.
| Geometries and Reference interaction energies from the following articles:
|   Benzene Dimers from Sherrill et al. JPCA 113 10146 (2009).
|   Benzene-Hydrogen Sulfide from Sherrill et al. JPCA 113 10146 (2009).
|   Benzene-Methane from Sherrill et al. JPCA 113 10146 (2009).
|   Methane Dimer from Takatani et al. PCCP 9 6106 (2007).
|   Pyridine Dimers from Hohenstein et al. JPCA 113 878 (2009).
|   Collection into NBC10 from Burns et al. JCP 134 084107 (2011).
|   Revised reference interaction energies (NBC10A) from Marshall et al. JCP 135 194102 (2011).
|   Revised (pure T-zeta CC) reference interaction energies (NBC10B) from Smith et al. JPCL 7 2197 (2016).

- **cp**  ``'off'`` || ``'on'``

- **rlxd** ``'off'``

- **benchmark**

  - ``'NBC100'`` Burns et al. JCP 134 084107 (2011).
  - ``'NBC10A'`` Marshall et al. JCP 135 194102 (2011).
  - |dl| ``'NBC10B'`` |dr| Smith et al. JPCL 7 2197 (2016).

- **subset**

  - ``'small'``
  - ``'large'``
  - ``'equilibrium'``
  - ``'BzBz_S'`` dissociation curve for benzene dimer, sandwich
  - ``'BzBz_T'`` dissociation curve for benzene dimer, t-shaped
  - ``'BzBz_PD34'`` dissociation curve for benzene dimer, parallel displaced by 3.4A
  - ``'BzH2S'`` dissociation curve for benzene-H2S
  - ``'BzMe'`` dissociation curve for benzene-methane
  - ``'MeMe'`` dissociation curve for methane dimer
  - ``'PyPy_S2'`` dissociation curve for pyridine dimer, sandwich
  - ``'PyPy_T3'`` dissociation curve for pyridine dimer, t-shaped
  - ``'BzBz_PD32'`` dissociation curve for benzene dimer, parallel displaced by 3.2A
  - ``'BzBz_PD36'`` dissociation curve for benzene dimer, parallel displaced by 3.6A
  - ``'5min'`` five points on each dissociation curve incl. and surrounding equilibrium
  - ``'MX'`` mixed-influence systems
  - ``'DD'`` dispersion-dominated systems
  - ``'MXDDPP'`` pi-pi-type mixed and dispersion systems
  - ``'MXDDNP'`` non-pi-pi-type mixed and dispersion systems
  - ``'WALL'``

"""
import re
import qcdb

# <<< NBC10 Database Module >>>
dbse = 'NBC1'

# <<< Database Members >>>
AXIS_Rang = {}
AXIS_Rrat = {}

dist = [3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.5, 5.0, 5.5, 6.0, 6.5, 10.0]
BzBz_S = ['BzBz_S-' + str(d) for d in dist]
AXIS_Rang.update(dict(zip(BzBz_S, dist)))
AXIS_Rrat.update(dict(zip(BzBz_S, [d / 3.9 for d in dist])))

dist = [4.3, 4.35, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 6.0, 6.5, 7.0, 7.5, 8.0]
BzBz_T = ['BzBz_T-' + str(d) for d in dist]
AXIS_Rang.update(dict(zip(BzBz_T, dist)))
AXIS_Rrat.update(dict(zip(BzBz_T, [d / 5.0 for d in dist])))

dist = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
BzBz_PD34 = ['BzBz_PD34-' + str(d) for d in dist]
AXIS_Rang.update(dict(zip(BzBz_PD34, dist)))
AXIS_Rrat.update(dict(zip(BzBz_PD34, [d / 1.8 for d in dist])))

dist = [3.15, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.5, 4.75, 5.0, 5.25, 5.5, 6.0, 6.5, 7.0, 7.5]
BzH2S = ['BzH2S-' + str(d) for d in dist]
AXIS_Rang.update(dict(zip(BzH2S, dist)))
AXIS_Rrat.update(dict(zip(BzH2S, [d / 3.8 for d in dist])))

dist = [3.15, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 6.0]
BzMe = ['BzMe-' + str(d) for d in dist]
AXIS_Rang.update(dict(zip(BzMe, dist)))
AXIS_Rrat.update(dict(zip(BzMe, [d / 3.8 for d in dist])))

dist = [3.1, 3.15, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.6, 4.8, 5.0, 5.4, 5.8]
MeMe = ['MeMe-' + str(d) for d in dist]
AXIS_Rang.update(dict(zip(MeMe, dist)))
AXIS_Rrat.update(dict(zip(MeMe, [d / 3.6 for d in dist])))

dist = [3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.7, 5.0, 5.5, 6.0, 6.5, 7.0]
PyPy_S2 = ['PyPy_S2-' + str(d) for d in dist]
AXIS_Rang.update(dict(zip(PyPy_S2, dist)))
AXIS_Rrat.update(dict(zip(PyPy_S2, [d / 3.7 for d in dist])))

dist = [4.1, 4.27, 4.3, 4.36, 4.44, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.7, 6.0, 6.5, 7.0, 8.0, 9.0]
PyPy_T3 = ['PyPy_T3-' + str(d) for d in dist]
AXIS_Rang.update(dict(zip(PyPy_T3, dist)))
AXIS_Rrat.update(dict(zip(PyPy_T3, [d / 4.9 for d in dist])))

dist = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
BzBz_PD32 = ['BzBz_PD32-' + str(d) for d in dist]
AXIS_Rang.update(dict(zip(BzBz_PD32, dist)))
AXIS_Rrat.update(dict(zip(BzBz_PD32, [d / 1.9 for d in dist])))

dist = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
BzBz_PD36 = ['BzBz_PD36-' + str(d) for d in dist]
AXIS_Rang.update(dict(zip(BzBz_PD36, dist)))
AXIS_Rrat.update(dict(zip(BzBz_PD36, [d / 1.7 for d in dist])))

temp = [BzBz_S, BzBz_T, BzBz_PD34, BzH2S, BzMe, MeMe, PyPy_S2, PyPy_T3, BzBz_PD32, BzBz_PD36]
HRXN = sum(temp, [])

HRXN_5MIN = ['BzBz_S-3.7', 'BzBz_S-3.8', 'BzBz_S-3.9', 'BzBz_S-4.0', 'BzBz_S-4.1',
             'BzBz_T-4.8', 'BzBz_T-4.9', 'BzBz_T-5.0', 'BzBz_T-5.1', 'BzBz_T-5.2',
             'BzBz_PD34-1.6', 'BzBz_PD34-1.7', 'BzBz_PD34-1.8', 'BzBz_PD34-1.9', 'BzBz_PD34-2.0',
             'BzH2S-3.6', 'BzH2S-3.7', 'BzH2S-3.8', 'BzH2S-3.9', 'BzH2S-4.0',
             'BzMe-3.6', 'BzMe-3.7', 'BzMe-3.8', 'BzMe-3.9', 'BzMe-4.0',
             'MeMe-3.4', 'MeMe-3.5', 'MeMe-3.6', 'MeMe-3.7', 'MeMe-3.8',
             'PyPy_S2-3.5', 'PyPy_S2-3.6', 'PyPy_S2-3.7', 'PyPy_S2-3.8', 'PyPy_S2-3.9',
             'PyPy_T3-4.7', 'PyPy_T3-4.8', 'PyPy_T3-4.9', 'PyPy_T3-5.0', 'PyPy_T3-5.1',
             'BzBz_PD32-1.7', 'BzBz_PD32-1.8', 'BzBz_PD32-1.9', 'BzBz_PD32-2.0', 'BzBz_PD32-2.2',
             'BzBz_PD36-1.5', 'BzBz_PD36-1.6', 'BzBz_PD36-1.7', 'BzBz_PD36-1.8', 'BzBz_PD36-1.9']
HRXN_SM = ['BzMe-6.0', 'MeMe-5.0']
HRXN_LG = ['BzBz_T-5.0']
HRXN_EQ = ['BzBz_S-3.9', 'BzBz_T-5.0', 'BzBz_PD34-1.8', 'BzH2S-3.8', 'BzMe-3.8',
           'MeMe-3.6', 'PyPy_S2-3.7', 'PyPy_T3-4.9', 'BzBz_PD32-1.9', 'BzBz_PD36-1.7']
#HRXN_WALL = ['BzBz_T-4.3', 'BzBz_T-4.35', 'BzH2S-3.15', 'BzH2S-3.3', 'BzMe-3.15', 'MeMe-3.1',
#             'MeMe-3.15', 'PyPy_S2-3.2', 'PyPy_T3-4.27', 'PyPy_T3-4.36', 'PyPy_T3-4.44']
MX = sum([BzH2S, PyPy_T3, BzBz_PD32], [])
DD = sum([BzBz_S, BzBz_T, BzBz_PD34, BzMe, MeMe, PyPy_S2, BzBz_PD36], [])
MXDDPP = sum([BzBz_S, BzBz_PD34, PyPy_S2, BzBz_PD32, BzBz_PD36], [])
MXDDNP = sum([BzBz_T, BzH2S, BzMe, MeMe, PyPy_T3], [])

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV_CP = {}  # order of active reagents per counterpoise-corrected reaction
ACTV_SA = {}  # order of active reagents for non-supramolecular calculations
for rxn in HRXN:

    if (rxn in BzBz_S) or (rxn in BzBz_PD34) or (rxn in BzBz_PD32) or (rxn in BzBz_PD36):
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -2,
                                          '%s-Bz-mono-unCP'  % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Bz-mono-unCP'  % (dbse) ]

    elif rxn in BzBz_T:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                          '%s-Bz-mono-unCP'  % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn),
                                          '%s-%s-monoB-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Bz-mono-unCP'  % (dbse) ]

    elif rxn in BzH2S:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                          '%s-Bz-mono-unCP'  % (dbse)      : -1,
                                          '%s-H2S-mono-unCP' % (dbse)      : -1 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn),
                                          '%s-%s-monoB-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Bz-mono-unCP'  % (dbse),
                                          '%s-H2S-mono-unCP' % (dbse) ]

    elif rxn in BzMe:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                          '%s-Bz2-mono-unCP' % (dbse)      : -1,
                                          '%s-Me-mono-unCP'  % (dbse)      : -1 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn),
                                          '%s-%s-monoB-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Bz2-mono-unCP' % (dbse),
                                          '%s-Me-mono-unCP'  % (dbse) ]

    elif rxn in MeMe:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -2,
                                          '%s-Me-mono-unCP'  % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Me-mono-unCP'  % (dbse) ]

    elif rxn in PyPy_S2:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -2,
                                          '%s-Py-mono-unCP'  % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Py-mono-unCP'  % (dbse) ]

    elif rxn in PyPy_T3:
        RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                          '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                          '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                          '%s-Py-mono-unCP'  % (dbse)      : -2 }

        ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

        ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-%s-monoA-CP'   % (dbse, rxn),
                                          '%s-%s-monoB-CP'   % (dbse, rxn) ]

        ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                          '%s-Py-mono-unCP'  % (dbse) ]

# <<< Reference Values >>>
BIND = {}
# Original publication
#BIND_NBC100 = {}
# Current published revision
BIND_NBC10A = {}
BIND_NBC10A['%s-BzBz_S-3.2'  % (dbse)] =  3.462
BIND_NBC10A['%s-BzBz_S-3.3'  % (dbse)] =  1.484
BIND_NBC10A['%s-BzBz_S-3.4'  % (dbse)] =  0.147
BIND_NBC10A['%s-BzBz_S-3.5'  % (dbse)] = -0.724
BIND_NBC10A['%s-BzBz_S-3.6'  % (dbse)] = -1.259
BIND_NBC10A['%s-BzBz_S-3.7'  % (dbse)] = -1.558
BIND_NBC10A['%s-BzBz_S-3.8'  % (dbse)] = -1.693
BIND_NBC10A['%s-BzBz_S-3.9'  % (dbse)] = -1.717  # BzBz_S minimum
BIND_NBC10A['%s-BzBz_S-4.0'  % (dbse)] = -1.669
BIND_NBC10A['%s-BzBz_S-4.1'  % (dbse)] = -1.577
BIND_NBC10A['%s-BzBz_S-4.2'  % (dbse)] = -1.459
BIND_NBC10A['%s-BzBz_S-4.5'  % (dbse)] = -1.066
BIND_NBC10A['%s-BzBz_S-5.0'  % (dbse)] = -0.546
BIND_NBC10A['%s-BzBz_S-5.5'  % (dbse)] = -0.251
BIND_NBC10A['%s-BzBz_S-6.0'  % (dbse)] = -0.101
BIND_NBC10A['%s-BzBz_S-6.5'  % (dbse)] = -0.029
BIND_NBC10A['%s-BzBz_S-10.0' % (dbse)] =  0.018

BIND_NBC10A['%s-BzBz_T-4.3' % (dbse)] =  2.633
BIND_NBC10A['%s-BzBz_T-4.35' % (dbse)] = 1.535
BIND_NBC10A['%s-BzBz_T-4.4' % (dbse)] =  0.617
BIND_NBC10A['%s-BzBz_T-4.5' % (dbse)] = -0.769
BIND_NBC10A['%s-BzBz_T-4.6' % (dbse)] = -1.682
BIND_NBC10A['%s-BzBz_T-4.7' % (dbse)] = -2.246
BIND_NBC10A['%s-BzBz_T-4.8' % (dbse)] = -2.559
BIND_NBC10A['%s-BzBz_T-4.9' % (dbse)] = -2.693
BIND_NBC10A['%s-BzBz_T-5.0' % (dbse)] = -2.703  # BzBz_T minimum
BIND_NBC10A['%s-BzBz_T-5.1' % (dbse)] = -2.630
BIND_NBC10A['%s-BzBz_T-5.2' % (dbse)] = -2.506
BIND_NBC10A['%s-BzBz_T-5.3' % (dbse)] = -2.351
BIND_NBC10A['%s-BzBz_T-5.4' % (dbse)] = -2.181
BIND_NBC10A['%s-BzBz_T-5.5' % (dbse)] = -2.006
BIND_NBC10A['%s-BzBz_T-5.6' % (dbse)] = -1.834
BIND_NBC10A['%s-BzBz_T-6.0' % (dbse)] = -1.242
BIND_NBC10A['%s-BzBz_T-6.5' % (dbse)] = -0.752
BIND_NBC10A['%s-BzBz_T-7.0' % (dbse)] = -0.468
BIND_NBC10A['%s-BzBz_T-7.5' % (dbse)] = -0.302
BIND_NBC10A['%s-BzBz_T-8.0' % (dbse)] = -0.203

BIND_NBC10A['%s-BzBz_PD34-0.2' % (dbse)] =  0.029
BIND_NBC10A['%s-BzBz_PD34-0.4' % (dbse)] = -0.298
BIND_NBC10A['%s-BzBz_PD34-0.6' % (dbse)] = -0.768
BIND_NBC10A['%s-BzBz_PD34-0.8' % (dbse)] = -1.298
BIND_NBC10A['%s-BzBz_PD34-1.0' % (dbse)] = -1.802
BIND_NBC10A['%s-BzBz_PD34-1.2' % (dbse)] = -2.213
BIND_NBC10A['%s-BzBz_PD34-1.4' % (dbse)] = -2.497
BIND_NBC10A['%s-BzBz_PD34-1.5' % (dbse)] = -2.586
BIND_NBC10A['%s-BzBz_PD34-1.6' % (dbse)] = -2.643
BIND_NBC10A['%s-BzBz_PD34-1.7' % (dbse)] = -2.668
BIND_NBC10A['%s-BzBz_PD34-1.8' % (dbse)] = -2.670  # BzBz_PD34 minimum
BIND_NBC10A['%s-BzBz_PD34-1.9' % (dbse)] = -2.649
BIND_NBC10A['%s-BzBz_PD34-2.0' % (dbse)] = -2.611
BIND_NBC10A['%s-BzBz_PD34-2.2' % (dbse)] = -2.501
BIND_NBC10A['%s-BzBz_PD34-2.4' % (dbse)] = -2.377
BIND_NBC10A['%s-BzBz_PD34-2.6' % (dbse)] = -2.260
BIND_NBC10A['%s-BzBz_PD34-2.8' % (dbse)] = -2.163
BIND_NBC10A['%s-BzBz_PD34-3.0' % (dbse)] = -2.093

BIND_NBC10A['%s-BzH2S-3.15' % (dbse)] =  2.379
BIND_NBC10A['%s-BzH2S-3.2'  % (dbse)] =  1.236
BIND_NBC10A['%s-BzH2S-3.3'  % (dbse)] = -0.471
BIND_NBC10A['%s-BzH2S-3.4'  % (dbse)] = -1.584
BIND_NBC10A['%s-BzH2S-3.5'  % (dbse)] = -2.269
BIND_NBC10A['%s-BzH2S-3.6'  % (dbse)] = -2.649
BIND_NBC10A['%s-BzH2S-3.7'  % (dbse)] = -2.818
BIND_NBC10A['%s-BzH2S-3.8'  % (dbse)] = -2.843  # BzH2S minimum
BIND_NBC10A['%s-BzH2S-3.9'  % (dbse)] = -2.773
BIND_NBC10A['%s-BzH2S-4.0'  % (dbse)] = -2.645
BIND_NBC10A['%s-BzH2S-4.1'  % (dbse)] = -2.483
BIND_NBC10A['%s-BzH2S-4.2'  % (dbse)] = -2.305
BIND_NBC10A['%s-BzH2S-4.5'  % (dbse)] = -1.771
BIND_NBC10A['%s-BzH2S-4.75' % (dbse)] = -1.393
BIND_NBC10A['%s-BzH2S-5.0'  % (dbse)] = -1.092
BIND_NBC10A['%s-BzH2S-5.25' % (dbse)] = -0.859
BIND_NBC10A['%s-BzH2S-5.5'  % (dbse)] = -0.682
BIND_NBC10A['%s-BzH2S-6.0'  % (dbse)] = -0.444
BIND_NBC10A['%s-BzH2S-6.5'  % (dbse)] = -0.301
BIND_NBC10A['%s-BzH2S-7.0'  % (dbse)] = -0.212
BIND_NBC10A['%s-BzH2S-7.5'  % (dbse)] = -0.154

BIND_NBC10A['%s-BzMe-3.15' % (dbse)] = 1.283
BIND_NBC10A['%s-BzMe-3.2' % (dbse)] =  0.686
BIND_NBC10A['%s-BzMe-3.3' % (dbse)] = -0.213
BIND_NBC10A['%s-BzMe-3.4' % (dbse)] = -0.805
BIND_NBC10A['%s-BzMe-3.5' % (dbse)] = -1.173
BIND_NBC10A['%s-BzMe-3.6' % (dbse)] = -1.378
BIND_NBC10A['%s-BzMe-3.7' % (dbse)] = -1.470
BIND_NBC10A['%s-BzMe-3.8' % (dbse)] = -1.484  # BzMe minimum
BIND_NBC10A['%s-BzMe-3.9' % (dbse)] = -1.445
BIND_NBC10A['%s-BzMe-4.0' % (dbse)] = -1.374
BIND_NBC10A['%s-BzMe-4.1' % (dbse)] = -1.284
BIND_NBC10A['%s-BzMe-4.2' % (dbse)] = -1.185
BIND_NBC10A['%s-BzMe-4.4' % (dbse)] = -0.984
BIND_NBC10A['%s-BzMe-4.6' % (dbse)] = -0.800
BIND_NBC10A['%s-BzMe-4.8' % (dbse)] = -0.643
BIND_NBC10A['%s-BzMe-5.0' % (dbse)] = -0.515
BIND_NBC10A['%s-BzMe-5.2' % (dbse)] = -0.413
BIND_NBC10A['%s-BzMe-5.4' % (dbse)] = -0.332
BIND_NBC10A['%s-BzMe-5.6' % (dbse)] = -0.268
BIND_NBC10A['%s-BzMe-6.0' % (dbse)] = -0.177

BIND_NBC10A['%s-MeMe-3.1' % (dbse)] =  0.570
BIND_NBC10A['%s-MeMe-3.15' % (dbse)] = 0.291
BIND_NBC10A['%s-MeMe-3.2' % (dbse)] =  0.069
BIND_NBC10A['%s-MeMe-3.3' % (dbse)] = -0.239
BIND_NBC10A['%s-MeMe-3.4' % (dbse)] = -0.417
BIND_NBC10A['%s-MeMe-3.5' % (dbse)] = -0.508
BIND_NBC10A['%s-MeMe-3.6' % (dbse)] = -0.541  # MeMe minimum
BIND_NBC10A['%s-MeMe-3.7' % (dbse)] = -0.539
BIND_NBC10A['%s-MeMe-3.8' % (dbse)] = -0.515
BIND_NBC10A['%s-MeMe-3.9' % (dbse)] = -0.480
BIND_NBC10A['%s-MeMe-4.0' % (dbse)] = -0.439
BIND_NBC10A['%s-MeMe-4.1' % (dbse)] = -0.396
BIND_NBC10A['%s-MeMe-4.2' % (dbse)] = -0.354
BIND_NBC10A['%s-MeMe-4.3' % (dbse)] = -0.315
BIND_NBC10A['%s-MeMe-4.4' % (dbse)] = -0.279
BIND_NBC10A['%s-MeMe-4.6' % (dbse)] = -0.217
BIND_NBC10A['%s-MeMe-4.8' % (dbse)] = -0.168
BIND_NBC10A['%s-MeMe-5.0' % (dbse)] = -0.130
BIND_NBC10A['%s-MeMe-5.4' % (dbse)] = -0.080
BIND_NBC10A['%s-MeMe-5.8' % (dbse)] = -0.050

BIND_NBC10A['%s-PyPy_S2-3.1' % (dbse)] =  2.387
BIND_NBC10A['%s-PyPy_S2-3.2' % (dbse)] =  0.245
BIND_NBC10A['%s-PyPy_S2-3.3' % (dbse)] = -1.165
BIND_NBC10A['%s-PyPy_S2-3.4' % (dbse)] = -2.050
BIND_NBC10A['%s-PyPy_S2-3.5' % (dbse)] = -2.562
BIND_NBC10A['%s-PyPy_S2-3.6' % (dbse)] = -2.815
BIND_NBC10A['%s-PyPy_S2-3.7' % (dbse)] = -2.890  # PyPy_S2 minimum
BIND_NBC10A['%s-PyPy_S2-3.8' % (dbse)] = -2.849
BIND_NBC10A['%s-PyPy_S2-3.9' % (dbse)] = -2.733
BIND_NBC10A['%s-PyPy_S2-4.0' % (dbse)] = -2.573
BIND_NBC10A['%s-PyPy_S2-4.1' % (dbse)] = -2.391
BIND_NBC10A['%s-PyPy_S2-4.2' % (dbse)] = -2.201
BIND_NBC10A['%s-PyPy_S2-4.3' % (dbse)] = -2.012
BIND_NBC10A['%s-PyPy_S2-4.4' % (dbse)] = -1.830
BIND_NBC10A['%s-PyPy_S2-4.5' % (dbse)] = -1.660
BIND_NBC10A['%s-PyPy_S2-4.7' % (dbse)] = -1.357
BIND_NBC10A['%s-PyPy_S2-5.0' % (dbse)] = -1.002
BIND_NBC10A['%s-PyPy_S2-5.5' % (dbse)] = -0.619
BIND_NBC10A['%s-PyPy_S2-6.0' % (dbse)] = -0.402
BIND_NBC10A['%s-PyPy_S2-6.5' % (dbse)] = -0.276
BIND_NBC10A['%s-PyPy_S2-7.0' % (dbse)] = -0.200

BIND_NBC10A['%s-PyPy_T3-4.1' % (dbse)] =  9.341
BIND_NBC10A['%s-PyPy_T3-4.27' % (dbse)] = 2.778
BIND_NBC10A['%s-PyPy_T3-4.3' % (dbse)] =  1.991
BIND_NBC10A['%s-PyPy_T3-4.36' % (dbse)] = 0.671
BIND_NBC10A['%s-PyPy_T3-4.44' % (dbse)] = -0.650
BIND_NBC10A['%s-PyPy_T3-4.5' % (dbse)] = -1.377
BIND_NBC10A['%s-PyPy_T3-4.6' % (dbse)] = -2.203
BIND_NBC10A['%s-PyPy_T3-4.7' % (dbse)] = -2.673
BIND_NBC10A['%s-PyPy_T3-4.8' % (dbse)] = -2.896
BIND_NBC10A['%s-PyPy_T3-4.9' % (dbse)] = -2.954  # PyPy_T3 minimum
BIND_NBC10A['%s-PyPy_T3-5.0' % (dbse)] = -2.903
BIND_NBC10A['%s-PyPy_T3-5.1' % (dbse)] = -2.783
BIND_NBC10A['%s-PyPy_T3-5.2' % (dbse)] = -2.625
BIND_NBC10A['%s-PyPy_T3-5.3' % (dbse)] = -2.447
BIND_NBC10A['%s-PyPy_T3-5.4' % (dbse)] = -2.262
BIND_NBC10A['%s-PyPy_T3-5.5' % (dbse)] = -2.080
BIND_NBC10A['%s-PyPy_T3-5.7' % (dbse)] = -1.741
BIND_NBC10A['%s-PyPy_T3-6.0' % (dbse)] = -1.323
BIND_NBC10A['%s-PyPy_T3-6.5' % (dbse)] = -0.852
BIND_NBC10A['%s-PyPy_T3-7.0' % (dbse)] = -0.573
BIND_NBC10A['%s-PyPy_T3-8.0' % (dbse)] = -0.296
BIND_NBC10A['%s-PyPy_T3-9.0' % (dbse)] = -0.174

BIND_NBC10A['%s-BzBz_PD32-0.2' % (dbse)] =  3.241
BIND_NBC10A['%s-BzBz_PD32-0.4' % (dbse)] =  2.619
BIND_NBC10A['%s-BzBz_PD32-0.6' % (dbse)] =  1.726
BIND_NBC10A['%s-BzBz_PD32-0.8' % (dbse)] =  0.726
BIND_NBC10A['%s-BzBz_PD32-1.0' % (dbse)] = -0.222
BIND_NBC10A['%s-BzBz_PD32-1.2' % (dbse)] = -1.002
BIND_NBC10A['%s-BzBz_PD32-1.4' % (dbse)] = -1.553
BIND_NBC10A['%s-BzBz_PD32-1.5' % (dbse)] = -1.738
BIND_NBC10A['%s-BzBz_PD32-1.6' % (dbse)] = -1.868
BIND_NBC10A['%s-BzBz_PD32-1.7' % (dbse)] = -1.949
BIND_NBC10A['%s-BzBz_PD32-1.8' % (dbse)] = -1.988
BIND_NBC10A['%s-BzBz_PD32-1.9' % (dbse)] = -1.992  # BzBz_PD32 minimum
BIND_NBC10A['%s-BzBz_PD32-2.0' % (dbse)] = -1.971
BIND_NBC10A['%s-BzBz_PD32-2.2' % (dbse)] = -1.891
BIND_NBC10A['%s-BzBz_PD32-2.4' % (dbse)] = -1.795
BIND_NBC10A['%s-BzBz_PD32-2.6' % (dbse)] = -1.727
BIND_NBC10A['%s-BzBz_PD32-2.8' % (dbse)] = -1.702
BIND_NBC10A['%s-BzBz_PD32-3.0' % (dbse)] = -1.725

BIND_NBC10A['%s-BzBz_PD36-0.2' % (dbse)] = -1.321
BIND_NBC10A['%s-BzBz_PD36-0.4' % (dbse)] = -1.490
BIND_NBC10A['%s-BzBz_PD36-0.6' % (dbse)] = -1.735
BIND_NBC10A['%s-BzBz_PD36-0.8' % (dbse)] = -2.011
BIND_NBC10A['%s-BzBz_PD36-1.0' % (dbse)] = -2.273
BIND_NBC10A['%s-BzBz_PD36-1.2' % (dbse)] = -2.482
BIND_NBC10A['%s-BzBz_PD36-1.4' % (dbse)] = -2.619
BIND_NBC10A['%s-BzBz_PD36-1.5' % (dbse)] = -2.657
BIND_NBC10A['%s-BzBz_PD36-1.6' % (dbse)] = -2.674
BIND_NBC10A['%s-BzBz_PD36-1.7' % (dbse)] = -2.675  # BzBz_PD36 minimum
BIND_NBC10A['%s-BzBz_PD36-1.8' % (dbse)] = -2.662
BIND_NBC10A['%s-BzBz_PD36-1.9' % (dbse)] = -2.633
BIND_NBC10A['%s-BzBz_PD36-2.0' % (dbse)] = -2.593
BIND_NBC10A['%s-BzBz_PD36-2.2' % (dbse)] = -2.489
BIND_NBC10A['%s-BzBz_PD36-2.4' % (dbse)] = -2.371
BIND_NBC10A['%s-BzBz_PD36-2.6' % (dbse)] = -2.253
BIND_NBC10A['%s-BzBz_PD36-2.8' % (dbse)] = -2.143
BIND_NBC10A['%s-BzBz_PD36-3.0' % (dbse)] = -2.046
# Current revision (pure triple-zeta)
BIND_NBC10B = {}
BIND_NBC10B['%s-BzBz_S-3.2'    % (dbse)] =  3.459
BIND_NBC10B['%s-BzBz_S-3.3'    % (dbse)] =  1.484
BIND_NBC10B['%s-BzBz_S-3.4'    % (dbse)] =  0.149
BIND_NBC10B['%s-BzBz_S-3.5'    % (dbse)] = -0.721
BIND_NBC10B['%s-BzBz_S-3.6'    % (dbse)] = -1.256
BIND_NBC10B['%s-BzBz_S-3.7'    % (dbse)] = -1.556
BIND_NBC10B['%s-BzBz_S-3.8'    % (dbse)] = -1.693
BIND_NBC10B['%s-BzBz_S-3.9'    % (dbse)] = -1.719  # BzBz_S minimum; CCSD(T)/[aQ5Z; D:haTZ]
BIND_NBC10B['%s-BzBz_S-4.0'    % (dbse)] = -1.672
BIND_NBC10B['%s-BzBz_S-4.1'    % (dbse)] = -1.582
BIND_NBC10B['%s-BzBz_S-4.2'    % (dbse)] = -1.464
BIND_NBC10B['%s-BzBz_S-4.5'    % (dbse)] = -1.072
BIND_NBC10B['%s-BzBz_S-5.0'    % (dbse)] = -0.550
BIND_NBC10B['%s-BzBz_S-5.5'    % (dbse)] = -0.252
BIND_NBC10B['%s-BzBz_S-6.0'    % (dbse)] = -0.101
BIND_NBC10B['%s-BzBz_S-6.5'    % (dbse)] = -0.029
BIND_NBC10B['%s-BzBz_S-10.0'   % (dbse)] =  0.018

BIND_NBC10B['%s-BzBz_T-4.3'    % (dbse)] =  2.631
BIND_NBC10B['%s-BzBz_T-4.35'   % (dbse)] =  1.533
BIND_NBC10B['%s-BzBz_T-4.4'    % (dbse)] =  0.616
BIND_NBC10B['%s-BzBz_T-4.5'    % (dbse)] = -0.770
BIND_NBC10B['%s-BzBz_T-4.6'    % (dbse)] = -1.682
BIND_NBC10B['%s-BzBz_T-4.7'    % (dbse)] = -2.246
BIND_NBC10B['%s-BzBz_T-4.8'    % (dbse)] = -2.559
BIND_NBC10B['%s-BzBz_T-4.9'    % (dbse)] = -2.693
BIND_NBC10B['%s-BzBz_T-5.0'    % (dbse)] = -2.702  # BzBz_T minimum; CCSD(T)/[aQ5Z; D:haTZ]
BIND_NBC10B['%s-BzBz_T-5.1'    % (dbse)] = -2.629
BIND_NBC10B['%s-BzBz_T-5.2'    % (dbse)] = -2.504
BIND_NBC10B['%s-BzBz_T-5.3'    % (dbse)] = -2.349
BIND_NBC10B['%s-BzBz_T-5.4'    % (dbse)] = -2.180
BIND_NBC10B['%s-BzBz_T-5.5'    % (dbse)] = -2.006
BIND_NBC10B['%s-BzBz_T-5.6'    % (dbse)] = -1.834
BIND_NBC10B['%s-BzBz_T-6.0'    % (dbse)] = -1.243
BIND_NBC10B['%s-BzBz_T-6.5'    % (dbse)] = -0.753
BIND_NBC10B['%s-BzBz_T-7.0'    % (dbse)] = -0.468
BIND_NBC10B['%s-BzBz_T-7.5'    % (dbse)] = -0.302
BIND_NBC10B['%s-BzBz_T-8.0'    % (dbse)] = -0.203

BIND_NBC10B['%s-BzBz_PD34-0.2' % (dbse)] =  0.032
BIND_NBC10B['%s-BzBz_PD34-0.4' % (dbse)] = -0.294
BIND_NBC10B['%s-BzBz_PD34-0.6' % (dbse)] = -0.766
BIND_NBC10B['%s-BzBz_PD34-0.8' % (dbse)] = -1.296
BIND_NBC10B['%s-BzBz_PD34-1.0' % (dbse)] = -1.800
BIND_NBC10B['%s-BzBz_PD34-1.2' % (dbse)] = -2.210
BIND_NBC10B['%s-BzBz_PD34-1.4' % (dbse)] = -2.494
BIND_NBC10B['%s-BzBz_PD34-1.5' % (dbse)] = -2.582
BIND_NBC10B['%s-BzBz_PD34-1.6' % (dbse)] = -2.639
BIND_NBC10B['%s-BzBz_PD34-1.7' % (dbse)] = -2.664
BIND_NBC10B['%s-BzBz_PD34-1.8' % (dbse)] = -2.666  # BzBz_PD34 minimum; CCSD(T)/[aQ5Z; D:haTZ]
BIND_NBC10B['%s-BzBz_PD34-1.9' % (dbse)] = -2.646
BIND_NBC10B['%s-BzBz_PD34-2.0' % (dbse)] = -2.608
BIND_NBC10B['%s-BzBz_PD34-2.2' % (dbse)] = -2.500
BIND_NBC10B['%s-BzBz_PD34-2.4' % (dbse)] = -2.375
BIND_NBC10B['%s-BzBz_PD34-2.6' % (dbse)] = -2.259
BIND_NBC10B['%s-BzBz_PD34-2.8' % (dbse)] = -2.162
BIND_NBC10B['%s-BzBz_PD34-3.0' % (dbse)] = -2.092

BIND_NBC10B['%s-BzH2S-3.15'    % (dbse)] =  2.386
BIND_NBC10B['%s-BzH2S-3.2'     % (dbse)] =  1.246
BIND_NBC10B['%s-BzH2S-3.3'     % (dbse)] = -0.458
BIND_NBC10B['%s-BzH2S-3.4'     % (dbse)] = -1.570
BIND_NBC10B['%s-BzH2S-3.5'     % (dbse)] = -2.254
BIND_NBC10B['%s-BzH2S-3.6'     % (dbse)] = -2.635
BIND_NBC10B['%s-BzH2S-3.7'     % (dbse)] = -2.805
BIND_NBC10B['%s-BzH2S-3.8'     % (dbse)] = -2.831  # BzH2S minimum; CCSD(T)/[aQ5Z; D:aTZ]
BIND_NBC10B['%s-BzH2S-3.9'     % (dbse)] = -2.763
BIND_NBC10B['%s-BzH2S-4.0'     % (dbse)] = -2.636
BIND_NBC10B['%s-BzH2S-4.1'     % (dbse)] = -2.476
BIND_NBC10B['%s-BzH2S-4.2'     % (dbse)] = -2.299
BIND_NBC10B['%s-BzH2S-4.5'     % (dbse)] = -1.768
BIND_NBC10B['%s-BzH2S-4.75'    % (dbse)] = -1.390
BIND_NBC10B['%s-BzH2S-5.0'     % (dbse)] = -1.089
BIND_NBC10B['%s-BzH2S-5.25'    % (dbse)] = -0.857
BIND_NBC10B['%s-BzH2S-5.5'     % (dbse)] = -0.680
BIND_NBC10B['%s-BzH2S-6.0'     % (dbse)] = -0.442
BIND_NBC10B['%s-BzH2S-6.5'     % (dbse)] = -0.300
BIND_NBC10B['%s-BzH2S-7.0'     % (dbse)] = -0.212
BIND_NBC10B['%s-BzH2S-7.5'     % (dbse)] = -0.154

BIND_NBC10B['%s-BzMe-3.15'     % (dbse)] =  1.281
BIND_NBC10B['%s-BzMe-3.2'      % (dbse)] =  0.685
BIND_NBC10B['%s-BzMe-3.3'      % (dbse)] = -0.213
BIND_NBC10B['%s-BzMe-3.4'      % (dbse)] = -0.804
BIND_NBC10B['%s-BzMe-3.5'      % (dbse)] = -1.171
BIND_NBC10B['%s-BzMe-3.6'      % (dbse)] = -1.376
BIND_NBC10B['%s-BzMe-3.7'      % (dbse)] = -1.468
BIND_NBC10B['%s-BzMe-3.8'      % (dbse)] = -1.481  # BzMe minimum; CCSD(T)/[aQ5Z; D:aTZ]
BIND_NBC10B['%s-BzMe-3.9'      % (dbse)] = -1.443
BIND_NBC10B['%s-BzMe-4.0'      % (dbse)] = -1.373
BIND_NBC10B['%s-BzMe-4.1'      % (dbse)] = -1.283
BIND_NBC10B['%s-BzMe-4.2'      % (dbse)] = -1.185
BIND_NBC10B['%s-BzMe-4.4'      % (dbse)] = -0.985
BIND_NBC10B['%s-BzMe-4.6'      % (dbse)] = -0.801
BIND_NBC10B['%s-BzMe-4.8'      % (dbse)] = -0.645
BIND_NBC10B['%s-BzMe-5.0'      % (dbse)] = -0.517
BIND_NBC10B['%s-BzMe-5.2'      % (dbse)] = -0.414
BIND_NBC10B['%s-BzMe-5.4'      % (dbse)] = -0.332
BIND_NBC10B['%s-BzMe-5.6'      % (dbse)] = -0.268
BIND_NBC10B['%s-BzMe-6.0'      % (dbse)] = -0.177

BIND_NBC10B['%s-MeMe-3.1'      % (dbse)] =  0.570
BIND_NBC10B['%s-MeMe-3.15'     % (dbse)] =  0.291
BIND_NBC10B['%s-MeMe-3.2'      % (dbse)] =  0.069
BIND_NBC10B['%s-MeMe-3.3'      % (dbse)] = -0.239
BIND_NBC10B['%s-MeMe-3.4'      % (dbse)] = -0.417
BIND_NBC10B['%s-MeMe-3.5'      % (dbse)] = -0.508
BIND_NBC10B['%s-MeMe-3.6'      % (dbse)] = -0.541  # MeMe minimum; CCSD(T)/aTQZ; unchgd NBC10A
BIND_NBC10B['%s-MeMe-3.7'      % (dbse)] = -0.539
BIND_NBC10B['%s-MeMe-3.8'      % (dbse)] = -0.515
BIND_NBC10B['%s-MeMe-3.9'      % (dbse)] = -0.480
BIND_NBC10B['%s-MeMe-4.0'      % (dbse)] = -0.439
BIND_NBC10B['%s-MeMe-4.1'      % (dbse)] = -0.396
BIND_NBC10B['%s-MeMe-4.2'      % (dbse)] = -0.354
BIND_NBC10B['%s-MeMe-4.3'      % (dbse)] = -0.315
BIND_NBC10B['%s-MeMe-4.4'      % (dbse)] = -0.279
BIND_NBC10B['%s-MeMe-4.6'      % (dbse)] = -0.217
BIND_NBC10B['%s-MeMe-4.8'      % (dbse)] = -0.168
BIND_NBC10B['%s-MeMe-5.0'      % (dbse)] = -0.130
BIND_NBC10B['%s-MeMe-5.4'      % (dbse)] = -0.080
BIND_NBC10B['%s-MeMe-5.8'      % (dbse)] = -0.050

BIND_NBC10B['%s-PyPy_S2-3.1'   % (dbse)] =  2.368
BIND_NBC10B['%s-PyPy_S2-3.2'   % (dbse)] =  0.233
BIND_NBC10B['%s-PyPy_S2-3.3'   % (dbse)] = -1.170
BIND_NBC10B['%s-PyPy_S2-3.4'   % (dbse)] = -2.051
BIND_NBC10B['%s-PyPy_S2-3.5'   % (dbse)] = -2.562
BIND_NBC10B['%s-PyPy_S2-3.6'   % (dbse)] = -2.814
BIND_NBC10B['%s-PyPy_S2-3.7'   % (dbse)] = -2.890  # PyPy_S2 minimum; CCSD(T)/[aQ5Z; D:haTZ]
BIND_NBC10B['%s-PyPy_S2-3.8'   % (dbse)] = -2.849
BIND_NBC10B['%s-PyPy_S2-3.9'   % (dbse)] = -2.734
BIND_NBC10B['%s-PyPy_S2-4.0'   % (dbse)] = -2.575
BIND_NBC10B['%s-PyPy_S2-4.1'   % (dbse)] = -2.394
BIND_NBC10B['%s-PyPy_S2-4.2'   % (dbse)] = -2.204
BIND_NBC10B['%s-PyPy_S2-4.3'   % (dbse)] = -2.016
BIND_NBC10B['%s-PyPy_S2-4.4'   % (dbse)] = -1.835
BIND_NBC10B['%s-PyPy_S2-4.5'   % (dbse)] = -1.664
BIND_NBC10B['%s-PyPy_S2-4.7'   % (dbse)] = -1.361
BIND_NBC10B['%s-PyPy_S2-5.0'   % (dbse)] = -1.004
BIND_NBC10B['%s-PyPy_S2-5.5'   % (dbse)] = -0.620
BIND_NBC10B['%s-PyPy_S2-6.0'   % (dbse)] = -0.403
BIND_NBC10B['%s-PyPy_S2-6.5'   % (dbse)] = -0.277
BIND_NBC10B['%s-PyPy_S2-7.0'   % (dbse)] = -0.200

BIND_NBC10B['%s-PyPy_T3-4.1'   % (dbse)] =  9.342
BIND_NBC10B['%s-PyPy_T3-4.27'  % (dbse)] =  2.792
BIND_NBC10B['%s-PyPy_T3-4.3'   % (dbse)] =  2.006
BIND_NBC10B['%s-PyPy_T3-4.36'  % (dbse)] =  0.689
BIND_NBC10B['%s-PyPy_T3-4.44'  % (dbse)] = -0.629
BIND_NBC10B['%s-PyPy_T3-4.5'   % (dbse)] = -1.355
BIND_NBC10B['%s-PyPy_T3-4.6'   % (dbse)] = -2.180
BIND_NBC10B['%s-PyPy_T3-4.7'   % (dbse)] = -2.650
BIND_NBC10B['%s-PyPy_T3-4.8'   % (dbse)] = -2.875
BIND_NBC10B['%s-PyPy_T3-4.9'   % (dbse)] = -2.934  # PyPy_T3 minimum; CCSD(T)/[aQ5Z; D:haTZ]
BIND_NBC10B['%s-PyPy_T3-5.0'   % (dbse)] = -2.885
BIND_NBC10B['%s-PyPy_T3-5.1'   % (dbse)] = -2.768
BIND_NBC10B['%s-PyPy_T3-5.2'   % (dbse)] = -2.612
BIND_NBC10B['%s-PyPy_T3-5.3'   % (dbse)] = -2.436
BIND_NBC10B['%s-PyPy_T3-5.4'   % (dbse)] = -2.254
BIND_NBC10B['%s-PyPy_T3-5.5'   % (dbse)] = -2.073
BIND_NBC10B['%s-PyPy_T3-5.7'   % (dbse)] = -1.737
BIND_NBC10B['%s-PyPy_T3-6.0'   % (dbse)] = -1.322
BIND_NBC10B['%s-PyPy_T3-6.5'   % (dbse)] = -0.853
BIND_NBC10B['%s-PyPy_T3-7.0'   % (dbse)] = -0.574
BIND_NBC10B['%s-PyPy_T3-8.0'   % (dbse)] = -0.297
BIND_NBC10B['%s-PyPy_T3-9.0'   % (dbse)] = -0.175

BIND_NBC10B['%s-BzBz_PD32-0.2' % (dbse)] =  3.239
BIND_NBC10B['%s-BzBz_PD32-0.4' % (dbse)] =  2.617
BIND_NBC10B['%s-BzBz_PD32-0.6' % (dbse)] =  1.724
BIND_NBC10B['%s-BzBz_PD32-0.8' % (dbse)] =  0.725
BIND_NBC10B['%s-BzBz_PD32-1.0' % (dbse)] = -0.223
BIND_NBC10B['%s-BzBz_PD32-1.2' % (dbse)] = -1.001
BIND_NBC10B['%s-BzBz_PD32-1.4' % (dbse)] = -1.551
BIND_NBC10B['%s-BzBz_PD32-1.5' % (dbse)] = -1.735
BIND_NBC10B['%s-BzBz_PD32-1.6' % (dbse)] = -1.865
BIND_NBC10B['%s-BzBz_PD32-1.7' % (dbse)] = -1.947
BIND_NBC10B['%s-BzBz_PD32-1.8' % (dbse)] = -1.986
BIND_NBC10B['%s-BzBz_PD32-1.9' % (dbse)] = -1.991  # BzBz_PD32 minimum; CCSD(T)/[aQ5Z; D:haTZ]
BIND_NBC10B['%s-BzBz_PD32-2.0' % (dbse)] = -1.971
BIND_NBC10B['%s-BzBz_PD32-2.2' % (dbse)] = -1.891
BIND_NBC10B['%s-BzBz_PD32-2.4' % (dbse)] = -1.795
BIND_NBC10B['%s-BzBz_PD32-2.6' % (dbse)] = -1.727
BIND_NBC10B['%s-BzBz_PD32-2.8' % (dbse)] = -1.702
BIND_NBC10B['%s-BzBz_PD32-3.0' % (dbse)] = -1.724

BIND_NBC10B['%s-BzBz_PD36-0.2' % (dbse)] = -1.318
BIND_NBC10B['%s-BzBz_PD36-0.4' % (dbse)] = -1.486
BIND_NBC10B['%s-BzBz_PD36-0.6' % (dbse)] = -1.732
BIND_NBC10B['%s-BzBz_PD36-0.8' % (dbse)] = -2.008
BIND_NBC10B['%s-BzBz_PD36-1.0' % (dbse)] = -2.271
BIND_NBC10B['%s-BzBz_PD36-1.2' % (dbse)] = -2.479
BIND_NBC10B['%s-BzBz_PD36-1.4' % (dbse)] = -2.616
BIND_NBC10B['%s-BzBz_PD36-1.5' % (dbse)] = -2.654
BIND_NBC10B['%s-BzBz_PD36-1.6' % (dbse)] = -2.671
BIND_NBC10B['%s-BzBz_PD36-1.7' % (dbse)] = -2.672  # BzBz_PD36 minimum; CCSD(T)/[aQ5Z; D:haTZ]
BIND_NBC10B['%s-BzBz_PD36-1.8' % (dbse)] = -2.659
BIND_NBC10B['%s-BzBz_PD36-1.9' % (dbse)] = -2.631
BIND_NBC10B['%s-BzBz_PD36-2.0' % (dbse)] = -2.591
BIND_NBC10B['%s-BzBz_PD36-2.2' % (dbse)] = -2.488
BIND_NBC10B['%s-BzBz_PD36-2.4' % (dbse)] = -2.370
BIND_NBC10B['%s-BzBz_PD36-2.6' % (dbse)] = -2.253
BIND_NBC10B['%s-BzBz_PD36-2.8' % (dbse)] = -2.142
BIND_NBC10B['%s-BzBz_PD36-3.0' % (dbse)] = -2.046
# Set default
BIND = BIND_NBC10B
# Reference information
BINDINFO_NBC10A = {}
for rxn in HRXN:
    if (rxn in BzBz_S) or (rxn in BzBz_T) or (rxn in BzBz_PD34) or \
       (rxn in PyPy_S2) or (rxn in BzBz_PD32) or (rxn in BzBz_PD36):
        BINDINFO_NBC10A['%s-%s' % (dbse, rxn)] = {'citation': 's22b', 'method': 'CCSDT', 'mode': 'CP', 'basis': 'atqzhatz'}
    elif (rxn in BzH2S) or (rxn in BzMe):
        BINDINFO_NBC10A['%s-%s' % (dbse, rxn)] = {'citation': 's22b', 'method': 'CCSDT', 'mode': 'CP', 'basis': 'atqzatz'}
    elif rxn in MeMe:
        BINDINFO_NBC10A['%s-%s' % (dbse, rxn)] = {'citation': 's22b', 'method': 'CCSDT', 'mode': 'CP', 'basis': 'atqz'}
    elif rxn in PyPy_T3:
        BINDINFO_NBC10A['%s-%s' % (dbse, rxn)] = {'citation': 's22b', 'method': 'CCSDT', 'mode': 'CP', 'basis': 'atqzadz'}
BINDINFO_NBC10B = {}
for rxn in HRXN:
    if (rxn in BzBz_S) or (rxn in BzBz_T) or (rxn in BzBz_PD34) or \
       (rxn in PyPy_S2) or (rxn in PyPy_T3) or (rxn in BzBz_PD32) or \
       (rxn in BzBz_PD36):
        BINDINFO_NBC10B['%s-%s' % (dbse, rxn)] = {'citation': 'dfit', 'method': 'CCSDT', 'mode': 'CP', 'basis': 'aq5zhatz'}
    elif (rxn in BzH2S) or (rxn in BzMe):
        BINDINFO_NBC10B['%s-%s' % (dbse, rxn)] = {'citation': 'dfit', 'method': 'CCSDT', 'mode': 'CP', 'basis': 'aq5zatz'}
    elif rxn in MeMe:
        BINDINFO_NBC10B['%s-%s' % (dbse, rxn)] = {'citation': 's22b', 'method': 'CCSDT', 'mode': 'CP', 'basis': 'atqz'}

# <<< Comment Lines >>>
TAGL = {}
rxnpattern = re.compile(r'^(.+)-(.+)$')

for item in BzBz_S:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Sandwich Benzene Dimer at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Sandwich Benzene Dimer at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Benzene from Sandwich Benzene Dimer at %s A' % (distance.group(2))

for item in BzBz_T:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'T-shaped Benzene Dimer at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'T-shaped Benzene Dimer at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Benzene from T-shaped Benzene Dimer at %s A' % (distance.group(2))
    TAGL['%s-%s-monoB-CP' % (dbse, item)] = 'Benzene from T-shaped Benzene Dimer at %s A' % (distance.group(2))

for item in BzBz_PD34:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Parallel Displaced Benzene Dimer Interplane 3.4 at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Parallel Displaced Benzene Dimer Interplane 3.4 at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Benzene from Parallel Displaced Benzene Dimer Interplane 3.4 at %s A' % (distance.group(2))

for item in BzH2S:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Benzene-H2S at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Benzene-H2S at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Benzene from Benzene-Methane at %s A' % (distance.group(2))
    TAGL['%s-%s-monoB-CP' % (dbse, item)] = 'Hydrogen Sulfide from Benzene-Methane at %s A' % (distance.group(2))

for item in BzMe:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Benzene-Methane at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Benzene-Methane at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Benzene from Benzene-Methane at %s A' % (distance.group(2))
    TAGL['%s-%s-monoB-CP' % (dbse, item)] = 'Methane from Benzene-Methane at %s A' % (distance.group(2))

for item in MeMe:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Methane Dimer at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Methane Dimer at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Methane from Methane Dimer at %s A' % (distance.group(2))

for item in PyPy_S2:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Pyridine Dimer S2 Configuration at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Pyridine Dimer S2 Configuration at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Pyridine from Pyridine Dimer S2 Configuration at %s A' % (distance.group(2))

for item in PyPy_T3:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Pyridine Dimer T3 Configuration at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Pyridine Dimer T3 Configuration at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Pyridine from Pyridine Dimer T3 Configuration at %s A' % (distance.group(2))
    TAGL['%s-%s-monoB-CP' % (dbse, item)] = 'Pyridine from Pyridine Dimer T3 Configuration at %s A' % (distance.group(2))

for item in BzBz_PD32:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Parallel Displaced Benzene Dimer Interplane 3.2 at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Parallel Displaced Benzene Dimer Interplane 3.2 at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Benzene from Parallel Displaced Benzene Dimer Interplane 3.2 at %s A' % (distance.group(2))

for item in BzBz_PD36:
    distance = rxnpattern.match(item)
    TAGL['%s-%s'          % (dbse, item)] = 'Parallel Displaced Benzene Dimer Interplane 3.6 at %s A' % (distance.group(2))
    TAGL['%s-%s-dimer'    % (dbse, item)] = 'Parallel Displaced Benzene Dimer Interplane 3.6 at %s A' % (distance.group(2))
    TAGL['%s-%s-monoA-CP' % (dbse, item)] = 'Benzene from Parallel Displaced Benzene Dimer Interplane 3.6 at %s A' % (distance.group(2))

TAGL['%s-Bz-mono-unCP'  % (dbse)] = 'Benzene'
TAGL['%s-H2S-mono-unCP' % (dbse)] = 'Hydrogen Sulfide'
TAGL['%s-Bz2-mono-unCP' % (dbse)] = 'Benzene (alt. geometry)'
TAGL['%s-Me-mono-unCP'  % (dbse)] = 'Methane'
TAGL['%s-Py-mono-unCP'  % (dbse)] = 'Pyridine'

TAGL['dbse'] = 'interaction energies for dissociation curves of dispersion-bound bimolecular complexes'
TAGL['BzBz_S'] = 'dissociation curve for benzene dimer, sandwich'
TAGL['BzBz_T'] = 'dissociation curve for benzene dimer, t-shaped'
TAGL['BzBz_PD34'] = 'dissociation curve for benzene dimer, parallel displaced by 3.4A'
TAGL['BzH2S'] = 'dissociation curve for benzene-H2S'
TAGL['BzMe'] = 'dissociation curve for benzene-methane'
TAGL['MeMe'] = 'dissociation curve for methane dimer'
TAGL['PyPy_S2'] = 'dissociation curve for pyridine dimer, sandwich'
TAGL['PyPy_T3'] = 'dissociation curve for pyridine dimer, t-shaped'
TAGL['BzBz_PD32'] = 'dissociation curve for benzene dimer, parallel displaced by 3.2A'
TAGL['BzBz_PD36'] = 'dissociation curve for benzene dimer, parallel displaced by 3.6A'
TAGL['5min'] = 'five points on each dissociation curve incl. and surrounding equilibrium'
TAGL['MX'] = 'mixed-influence systems'
TAGL['DD'] = 'dispersion-dominated systems'
TAGL['MXDDPP'] = 'pi-pi-type mixed and dispersion systems'
TAGL['MXDDNP'] = 'non-pi-pi-type mixed and dispersion systems'
TAGL['large'] = 'most computationally expensive systems'
TAGL['default'] = 'entire database'
TAGL['equilibrium'] = 'minimum-energy systems on dissociation curves'
TAGL['small'] = 'few computationally quick systems'

#<<< Geometry Specification Strings >>>
GEOS = {}

for rxn in BzBz_S:
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
X
X  1  RXX
X  2  RXX   1  90.0
C  3  RCC   2  90.0  1  0.0
C  3  RCC   2  90.0  1  60.0
C  3  RCC   2  90.0  1  120.0
C  3  RCC   2  90.0  1  180.0
C  3  RCC   2  90.0  1  240.0
C  3  RCC   2  90.0  1  300.0
H  3  RCH   2  90.0  1  0.0
H  3  RCH   2  90.0  1  60.0
H  3  RCH   2  90.0  1  120.0
H  3  RCH   2  90.0  1  180.0
H  3  RCH   2  90.0  1  240.0
H  3  RCH   2  90.0  1  300.0
--
0 1
X  3  RXX   2  90.0  1  0.0
X  3  R     16 90.0  2  180.0
X  3  DRXX  16 90.0  2  180.0
X  18 RXX   17 90.0  16 0.0
C  17 RCC   18 90.0  19 0.0
C  17 RCC   18 90.0  19 60.0
C  17 RCC   18 90.0  19 120.0
C  17 RCC   18 90.0  19 180.0
C  17 RCC   18 90.0  19 240.0
C  17 RCC   18 90.0  19 300.0
H  17 RCH   18 90.0  19 0.0
H  17 RCH   18 90.0  19 60.0
H  17 RCH   18 90.0  19 120.0
H  17 RCH   18 90.0  19 180.0
H  17 RCH   18 90.0  19 240.0
H  17 RCH   18 90.0  19 300.0

RXX    = 1.0
DRXX   = 12.0
RCC    = 1.3915
RCH    = 2.4715
R      = %(Rval)s
units angstrom
""" % vars())

for rxn in BzBz_T:
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
X
X  1  RXX
X  2  RXX  1  90.0
C  3  RCC  2  90.0  1  0.0
C  3  RCC  2  90.0  1  60.0
C  3  RCC  2  90.0  1  120.0
C  3  RCC  2  90.0  1  180.0
C  3  RCC  2  90.0  1  240.0
C  3  RCC  2  90.0  1  300.0
H  3  RCH  2  90.0  1  0.0
H  3  RCH  2  90.0  1  60.0
H  3  RCH  2  90.0  1  120.0
H  3  RCH  2  90.0  1  180.0
H  3  RCH  2  90.0  1  240.0
H  3  RCH  2  90.0  1  300.0
--
0 1
X  3  RXX  2  90.0  1  0.0
X  3  R    16 90.0  1  0.0
X  17 RXX  3  90.0  16 180.0
X  18 RXX  17 90.0  3  0.0
C  17 RCC  18 90.0  19 0.0
C  17 RCC  18 90.0  19 60.0
C  17 RCC  18 90.0  19 120.0
C  17 RCC  18 90.0  19 180.0
C  17 RCC  18 90.0  19 240.0
C  17 RCC  18 90.0  19 300.0
H  17 RCH  18 90.0  19 0.0
H  17 RCH  18 90.0  19 60.0
H  17 RCH  18 90.0  19 120.0
H  17 RCH  18 90.0  19 180.0
H  17 RCH  18 90.0  19 240.0
H  17 RCH  18 90.0  19 300.0

RXX    = 1.0
RCC    = 1.3915
RCH    = 2.4715
R      = %(Rval)s
units angstrom
""" % vars())

for rxn in sum([BzBz_PD32, BzBz_PD34, BzBz_PD36], []):
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))
    if rxn in BzBz_PD32:
        R2val = 3.2
    elif rxn in BzBz_PD34:
        R2val = 3.4
    elif rxn in BzBz_PD36:
        R2val = 3.6

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
X
X  1  RXX
X  2  R2   1  90.0
C  3  RCC  2  90.0  1  0.0
C  3  RCC  2  90.0  1  60.0
C  3  RCC  2  90.0  1  120.0
C  3  RCC  2  90.0  1  180.0
C  3  RCC  2  90.0  1  240.0
C  3  RCC  2  90.0  1  300.0
H  3  RCH  2  90.0  1  0.0
H  3  RCH  2  90.0  1  60.0
H  3  RCH  2  90.0  1  120.0
H  3  RCH  2  90.0  1  180.0
H  3  RCH  2  90.0  1  240.0
H  3  RCH  2  90.0  1  300.0
--
0 1
X  3  RXX  2  90.0  1  0.0
X  2  R    3  90.0  16 90.0
X  17 RXX  2  90.0  1  90.0
X  18 RXX  17 90.0  2  90.0
C  17 RCC  18 90.0  19 0.0
C  17 RCC  18 90.0  19 60.0
C  17 RCC  18 90.0  19 120.0
C  17 RCC  18 90.0  19 180.0
C  17 RCC  18 90.0  19 240.0
C  17 RCC  18 90.0  19 300.0
H  17 RCH  18 90.0  19 0.0
H  17 RCH  18 90.0  19 60.0
H  17 RCH  18 90.0  19 120.0
H  17 RCH  18 90.0  19 180.0
H  17 RCH  18 90.0  19 240.0
H  17 RCH  18 90.0  19 300.0

RXX    = 1.0
RCC    = 1.3915
RCH    = 2.4715
R      = %(Rval)s
R2     = %(R2val)s
units angstrom
""" % vars())

for rxn in BzH2S:
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
X
X  1  1.0
C  2  BZCX  1 90.0
C  2  BZCX  1 90.0   3 60.0
C  2  BZCX  1 90.0   4 60.0
C  2  BZCX  1 90.0   5 60.0
C  2  BZCX  1 90.0   6 60.0
C  2  BZCX  1 90.0   7 60.0
X  3  1.0   2 90.0   1 0.0
H  3  BZHC  9 90.0   2 180.0
H  4  BZHC  3 120.0  2 180.0
H  5  BZHC  4 120.0  2 180.0
H  6  BZHC  5 120.0  2 180.0
H  7  BZHC  6 120.0  2 180.0
H  8  BZHC  7 120.0  2 180.0
--
0 1
S  2  R     3 90.0   4 90.0
H  16 HS    2 HHSH   9 180.0
H  16 HS    2 HHSH   9 0.0

BZCX   = 1.3915
BZHC   = 1.0800
HS     = 1.3356
HHSH   = 46.06
R      = %(Rval)s
units angstrom
""" % vars())

for rxn in BzMe:
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
X
X  1  1.0
C  2  CQ   1  90.0
C  3  CQ   2  60.0   1  90.0
C  4  CQ   2  60.0   1  90.0
C  5  CQ   2  60.0   1  90.0
C  6  CQ   2  60.0   1  90.0
C  7  CQ   2  60.0   1  90.0
X  3  1.0  2  90.0   1  0.0
H  3  CH1  9  90.0   2  180.0
H  4  CH1  3  120.0  2  180.0
H  5  CH1  4  120.0  2  180.0
H  6  CH1  5  120.0  2  180.0
H  7  CH1  6  120.0  2  180.0
H  8  CH1  7  120.0  2  180.0
--
0 1
C  2  R    3  90.0   9  0.0
H  16 CH2  2  0.0    3  0.0
H  16 CH2  2  HCH    3  0.0
H  16 CH2  17 HCH    18 120.0
H  16 CH2  17 HCH    18 240.0

CQ     = 1.405731
CH1    = 1.095210
CH2    = 1.099503
HCH    = 109.471209
R      = %(Rval)s
units angstrom
""" % vars())

for rxn in MeMe:
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
C
H  1 CH2
H  1 CH2  2 HCH
H  1 CH2  2 HCH    3  120.0
H  1 CH2  2 HCH    3  240.0
--
0 1
C  1 R    2 180.0  4  120.0
H  6 CH2  2 180.0  4  120.0
H  6 CH2  7 HCH    3  180.0
H  6 CH2  7 HCH    4  180.0
H  6 CH2  7 HCH    5  180.0

CH2    = 1.099503
HCH    = 109.471209
R      = %(Rval)s
units angstrom
""" % vars())

for rxn in PyPy_S2:
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
X
X  1  R
N  1  1.3980380  2  90.0
C  1  1.3371053  2  90.0    3  -58.504950
C  1  1.3822904  2  90.0    4  -61.640500
C  1  1.4067471  2  90.0    5  -59.854550
C  1  1.3822904  2  90.0    6  -59.854550
C  1  1.3371053  2  90.0    7  -61.640500
H  4  1.08650    3  116.01  8  180.0
H  5  1.08260    4  120.12  3  180.0
H  6  1.08180    3  180.00  4  0.0
H  7  1.08260    8  120.12  3  180.0
H  8  1.08650    3  116.01  4  180.0
--
0 1
N  2  1.3980380  1  90.0    3  theta
C  2  1.3371053  1  90.0    14 -58.504950
C  2  1.3822904  1  90.0    15 -61.640500
C  2  1.4067471  1  90.0    16 -59.854550
C  2  1.3822904  1  90.0    17 -59.854550
C  2  1.3371053  1  90.0    18 -61.640500
H  15 1.08650    14 116.01  19 180.0
H  16 1.08260    15 120.12  14 180.0
H  17 1.08180    14 180.00  15 0.0
H  18 1.08260    19 120.12  14 180.0
H  19 1.08650    14 116.01  15 180.0

theta  = 180.0
R      = %(Rval)s
units angstrom
""" % vars())

for rxn in PyPy_T3:
    molname = rxnpattern.match(rxn)
    Rval = float(molname.group(2))

    GEOS['%s-%s-%s' % (dbse, rxn, 'dimer')] = qcdb.Molecule("""
0 1
X
X  1  R
N  1  1.3980380  2   90.0
C  1  1.3371053  2   90.0   3  -58.504950
C  1  1.3822904  2   90.0   4  -61.640500
C  1  1.4067471  2   90.0   5  -59.854550
C  1  1.3822904  2   90.0   6  -59.854550
C  1  1.3371053  2   90.0   7  -61.640500
H  4  1.08650    3  116.01  8  180.0
H  5  1.08260    4  120.12  3  180.0
H  6  1.08180    3  180.00  4    0.0
H  7  1.08260    8  120.12  3  180.0
H  8  1.08650    3  116.01  4  180.0
--
0 1
X  2  2.0000000  1   90.0   3  theta
N  2  1.3980380  14  90.0   1  updown
C  2  1.3371053  14  90.0   15 -58.504950
C  2  1.3822904  14  90.0   16 -61.640500
C  2  1.4067471  14  90.0   17 -59.854550
C  2  1.3822904  14  90.0   18 -59.854550
C  2  1.3371053  14  90.0   19 -61.640500
H  16 1.08650    15 116.01  20 180.0
H  17 1.08260    16 120.12  15 180.0
H  18 1.08180    15 180.00  16 0.0
H  19 1.08260    20 120.12  15 180.0
H  20 1.08650    15 116.01  16 180.0

theta  = 90.0
updown = 270.0
R      = %(Rval)s
units angstrom
""" % vars())

GEOS['%s-%s-%s' % (dbse, 'Bz', 'mono-unCP')] = qcdb.Molecule("""
0 1
X
X  1  RXX
X  2  RXX 1  90.0
C  3  RCC 2  90.0 1  0.0
C  3  RCC 2  90.0 1  60.0
C  3  RCC 2  90.0 1  120.0
C  3  RCC 2  90.0 1  180.0
C  3  RCC 2  90.0 1  240.0
C  3  RCC 2  90.0 1  300.0
H  3  RCH 2  90.0 1  0.0
H  3  RCH 2  90.0 1  60.0
H  3  RCH 2  90.0 1  120.0
H  3  RCH 2  90.0 1  180.0
H  3  RCH 2  90.0 1  240.0
H  3  RCH 2  90.0 1  300.0

RXX = 1.0
RCC = 1.3915
RCH = 2.4715
units angstrom
""" % vars())

GEOS['%s-%s-%s' % (dbse, 'H2S', 'mono-unCP')] = qcdb.Molecule("""
0 1
S
H  1 HS
H  1 HS  2 HSH

HS  = 1.3356
HSH = 92.12
units angstrom
""" % vars())

GEOS['%s-%s-%s' % (dbse, 'Bz2', 'mono-unCP')] = qcdb.Molecule("""
0 1
X
X  1  1.0
C  2  CQ   1  90.0
C  3  CQ   2  60.0  1  90.0
C  4  CQ   2  60.0  1  90.0
C  5  CQ   2  60.0  1  90.0
C  6  CQ   2  60.0  1  90.0
C  7  CQ   2  60.0  1  90.0
X  3  1.0  2  90.0  1  0.0
H  3  CH1  9  90.0  2  180.0
H  4  CH1  3  120.0 2  180.0
H  5  CH1  4  120.0 2  180.0
H  6  CH1  5  120.0 2  180.0
H  7  CH1  6  120.0 2  180.0
H  8  CH1  7  120.0 2  180.0

CQ  = 1.405731
CH1 = 1.095210
units angstrom
""" % vars())

GEOS['%s-%s-%s' % (dbse, 'Me', 'mono-unCP')] = qcdb.Molecule("""
0 1
C
H  1  CH2
H  1  CH2  2  HCH
H  1  CH2  2  HCH   3  120.0
H  1  CH2  2  HCH   3  240.0

CH2 = 1.099503
HCH = 109.471209
units angstrom
""" % vars())

GEOS['%s-%s-%s' % (dbse, 'Py', 'mono-unCP')] = qcdb.Molecule("""
0 1
X
X  1  RXX
N  1  1.3980380  2  90.0
C  1  1.3371053  2  90.0    3  -58.504950
C  1  1.3822904  2  90.0    4  -61.640500
C  1  1.4067471  2  90.0    5  -59.854550
C  1  1.3822904  2  90.0    6  -59.854550
C  1  1.3371053  2  90.0    7  -61.640500
H  4  1.08650    3  116.01  8  180.0
H  5  1.08260    4  120.12  3  180.0
H  6  1.08180    3  180.00  4  0.0
H  7  1.08260    8  120.12  3  180.0
H  8  1.08650    3  116.01  4  180.0

RXX = 1.0
units angstrom
""" % vars())

# <<< Derived Geometry Strings >>>
for rxn in HRXN:
    GEOS['%s-%s-monoA-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1, 2)
    if rxn in sum([BzBz_T, BzH2S, BzMe, PyPy_T3], []):
        GEOS['%s-%s-monoB-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2, 1)

#########################################################################
